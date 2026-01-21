
from __future__ import annotations

import mmap
import re
import sys
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd


@dataclass(frozen=True)
class _ParsedHeader:
    """A struct to hold the parsed data from a contig header."""
    name: str
    length: int
    is_circular: bool
    mylo_depth: str
    is_duplicated: bool


class AbundanceDB:
    """
    A memory-efficient, database-like interface for a single abundance TSV file,
    using memory mapping (mmap) to avoid loading the entire file into memory.
    This class is a context manager to ensure proper resource cleanup.
    """
    def __init__(self, depth_file_path: Path):
        self.depth_file_path = depth_file_path
        self._file = None
        self._mm: mmap.mmap = field(init=False)

        self._file = open(self.depth_file_path, "rb")
        self._mm = mmap.mmap(self._file.fileno(), 0, access=mmap.ACCESS_READ)
        
        # Read and process header once.
        self._mm.seek(0)
        header_line = self._mm.readline().decode().strip()
        self.header_fields = header_line.split('\t')
        
        original_var_cols = {h for h in self.header_fields if h.endswith("-var")}
        cleaned_var_names = {h.removesuffix("-var") for h in original_var_cols}
        
        self.original_keys_to_extract = []
        self.target_data_cols = []
        
        for hf in self.header_fields:
            if hf not in original_var_cols and hf in cleaned_var_names:
                self.original_keys_to_extract.append(hf)
                self.target_data_cols.append(hf.removesuffix("_merged.bam"))

        if not self.target_data_cols:
            raise ValueError(f"No valid sample columns found in header of '{self.depth_file_path}'.")

    def get_abund_for_contigs(self, contig_names_to_find: set[str]) -> dict[str, dict[str, float]]:
        """
        Finds and parses abundance data for a specific set of contigs in a single
        pass over the mmap'd file, using a substring search for each contig name.
        """
        results = {}
        contigs_left = contig_names_to_find.copy()

        self._mm.seek(0)
        _ = self._mm.readline()  # Skip header

        while contigs_left and (line_bytes := self._mm.readline()):
            found_name = None
            # Check if any of the required contig names exist as a substring in the current line.
            for name in contigs_left:
                if name.encode() in line_bytes:
                    found_name = name
                    break # Found a match, stop checking names for this line
            
            if found_name:
                try:
                    line_str = line_bytes.decode('utf-8').strip()
                    line_values = line_str.split('\t')
                    full_row_dict = dict(zip(self.header_fields, line_values))

                    result_dict = {
                        target_key: float(full_row_dict[original_key])
                        for target_key, original_key in zip(self.target_data_cols, self.original_keys_to_extract)
                        if original_key in full_row_dict
                    }

                    results[found_name] = result_dict
                    contigs_left.remove(found_name)
                except (UnicodeDecodeError, IndexError, ValueError):
                    # If line is malformed, skip it and continue.
                    continue

        return results

    def __enter__(self): return self
    def __exit__(self, exc_type, exc_val, exc_tb): self.close()

    def close(self):
        if self._mm and not self._mm.closed:
            self._mm.close()
        if self._file and not self._file.closed:
            self._file.close()


@dataclass
class ContigID:
    """
    Represents a single contig, parsing its metadata from the FASTA header
    and holding its abundance information.
    """
    _name: str  # The full header line from the FASTA file.
    _abund_info: dict[str, float]  # Abundance data for this contig across samples.
    _header: _ParsedHeader = field(init=False, repr=False)  # Parsed data from the header.

    HEADER_PATTERN = re.compile(
        r"^>?(?P<name>\w+)_len-(?P<length>\d+)_circular-(?P<circular>yes|possibly|no)_depth-(?P<depth>[\d.-]+)_duplicated-(?P<duplicated>yes|possibly|no)$"
    )

    def __post_init__(self):
        header_string = self._name.split(" ")[0]
        match = self.HEADER_PATTERN.match(header_string)

        if not match:
            raise ValueError(f"Unexpected contig header format: '{header_string}'")

        groups = match.groupdict()
        self._header = _ParsedHeader(
            name=groups['name'],
            length=int(groups['length']),
            is_circular=groups['circular'] == 'yes',
            mylo_depth=groups['depth'],
            is_duplicated=groups['duplicated'] == 'yes'
        )

    @staticmethod
    def parse_name_from_header(header_string: str) -> str | None:
        """Parses just the contig name from a full header string."""
        match = ContigID.HEADER_PATTERN.match(header_string.split(" ")[0])
        if match:
            return match.group('name')
        return None

    # Properties are now simple, fast accessors to the parsed data.
    @property
    def name(self) -> str: return self._header.name
    @property
    def length(self) -> int: return self._header.length
    @property
    def mylo_depth(self) -> str: return self._header.mylo_depth
    def is_circular(self) -> bool: return self._header.is_circular
    def is_duplicated(self) -> bool: return self._header.is_duplicated

    def __repr__(self) -> str: return ( f"{self.__class__.__name__}(name={self.name}, "
        f"length={self.length}, "
        f"is_circular={self.is_circular()}, "
        f"depth={self.mylo_depth}, "
        f"duplicated={self.is_duplicated()})" )

    def get_abund_info(self, rename=True) -> dict[str, float]:
        """Returns the pre-fetched abundance information for this contig."""
        return self._abund_info

    @property
    def depth_from_all_samples(self) -> float:
        # The dict may be empty if the contig was not found in the abundance file
        if not self._abund_info:
            return 0.0
        return sum(self._abund_info.values())

@dataclass
class Mag:
    """
    Represents a Metagenome-Assembled Genome (MAG), including its summary statistics,
    filesystem path, and a list of its constituent contigs.
    """
    name: str  # The unique name of the MAG (e.g., 'C1A3_A_metabat.872').
    _data: dict  # Raw data dictionary from the summary TSV file.
    _fp: Path  # Filesystem path to the MAG's FASTA file.
    contigids: list[ContigID]  # List of ContigID objects belonging to this MAG.

    @property
    def fp(self): return str(self._fp)

    @property
    def sample(self): return self.name.split("_")[0]
    @property
    def long_sample(self): return f"C00{self.sample[1]}_{self.sample[2:4]}"
    @property
    def binner(self): return self.name.split("_")[2].split(".")[0]
    @property
    def assembler(self): return { "A": "myloasm", "M": "medaka", "P": "proovframe" }[self.name.split("_")[1]]
    @property
    def completeness(self): return float(self._data[ 'Completeness' ])
    @property
    def contamination(self): return float(self._data[ 'Contamination' ])
    @property
    def contig_n50(self): return self._data[ 'Contig_N50' ]
    @property
    def gc_content(self): return float(self._data[ 'GC_Content' ])
    @property
    def max_contig_length(self): return self._data[ 'Max_Contig_Length' ]
    @property
    def genome_size(self): return self._data[ 'genome_size' ]
    @property
    def average_coverage(self): return self._data[ 'average_coverage' ]
    @property
    def total_contigs(self): return self._data[ 'Total_Contigs' ]
    @property
    def trna_counts(self): return self._data[ 'tRNA counts' ]
    @property
    def classification(self): return self._data[ 'classification' ]
    @property
    def red_value(self): return self._data[ 'red_value' ]
    def has_16s_rrna(self): return self._data[ '16S_rRNA' ]
    def has_23s_rrna(self): return self._data[ '23S_rRNA' ]
    def has_5s_rrna(self): return self._data[ '5S_rRNA' ]

    def __repr__(self): return ( f"{self.__class__.__name__}"
                                f"(name={self.name}, "
                                f"completeness={self.completeness}, "
                                f"contamination={self.contamination}, "
                                f"assembler={self.assembler}, "
                                f"sample={self.sample}, "
                                f"binner={self.binner}, "
                                f"total_contigs={self.total_contigs})" )

    def get_depth_df(self) -> pd.DataFrame:
        """
        Generates a pandas DataFrame of contig depths for all samples.

        Returns:
            pd.DataFrame: A DataFrame where rows are contig names and columns
                          are sample names, with values being the depth.
        """
        depths = [ contigid.get_abund_info()
                for contigid in self.contigids ]
        indices = [contigid.name for contigid in self.contigids]
        return pd.DataFrame(
            depths, index=pd.Index(indices)
        )

    def get_depth_per_contig(self, contigid: ContigID, samplename=None):
        if samplename:
            return contigid.get_abund_info()[samplename]
        return contigid.depth_from_all_samples

    @property
    def average_coverage_total(self) -> float:
        # values = [ (contigid.name, contigid.length, self.get_depth_per_contig(contigid) )
        #           for contigid in self.contigids ]
        total_length = 0
        total_bases = 0
        for c in self.contigids:
            total_bases += c.depth_from_all_samples * c.length
            total_length += c.length
        # for v in values:
        #     total_bases += v[1] * v[2]
        return total_bases / total_length

    def average_coverage_per_sample(self, samplename: Optional[str] = None ) -> float:
        if samplename == None:
            samplename = self.long_sample
        values = [ (contigid.name, contigid.length, self.get_depth_per_contig(contigid, samplename) )
                  for contigid in self.contigids ]
        total_bases = 0
        for v in values:
            total_bases += v[1] * v[2]
        return total_bases / sum([v[1] for v in values])

    def is_high_quality(self): return self.completeness >= 90 and self.contamination <= 5


@dataclass
class SummaryRepository:
    """
    Loads and provides access to the main summary data file, which contains
    statistics for all MAGs.
    """
    _summary_path: Path  # Path to the summary TSV file.
    summary: pd.DataFrame = field(init=False)

    def __post_init__(self):
        self.summary = pd.read_csv(self._summary_path, sep="\t", index_col=0)
        # Add MAG attributes as columns for querying
        self.summary['binner'] = self.summary.index.map(lambda x: x.split("_")[2].split(".")[0])
        self.summary['assembler'] = self.summary.index.map(lambda x: { "A": "myloasm", "M": "medaka", "P": "proovframe" }.get(x.split("_")[1]))


    def get_mag_data(self, mag_name: str) -> dict:
        return self.summary.loc[mag_name, :].to_dict()

    def get_summary_index(self) -> list[str]: return self.summary.index.tolist()

    def get_mags_by_query(self, query_string: str) -> list[str]:
        """
        Filters MAGs based on a pandas query string.

        Args:
            query_string: A string that is a valid pandas query.

        Returns:
            A list of MAG names that match the query.
        """
        try:
            return self.summary.query(query_string).index.tolist()
        except Exception as e:
            print(f"Error executing query: {e}", file=sys.stderr)
            return []

    def get_hq_mag_names(self) -> list[str]:
        """A convenience shortcut to get high-quality MAGs."""
        hq_query = "Completeness >= 90 and Contamination <= 5"
        return self.get_mags_by_query(hq_query)

@dataclass
class FilesystemLocator:
    """
    A helper class to locate MAG and abundance files on the filesystem.
    """
    mag_dir: Path  # Directory containing MAG FASTA files.
    abund_dir: Path  # Directory containing abundance TSV files.

    def find_mag_file(self, mag_name: str) -> Path:
        check_path = self.mag_dir / f"{mag_name}.fasta"
        if check_path.exists():
            return check_path
        raise FileNotFoundError(f"MAG file not found for {mag_name} in {self.mag_dir}")


class SessionManager:
    """
    The main orchestrator class.

    It ties together the summary data, filesystem locations, and abundance
    information to construct fully-realized Mag objects on demand.
    """
    def __init__(self, summary_path: Path, mag_dir: Path, abund_dir: Path) -> None:
        self.summary_repo = SummaryRepository(summary_path)
        self.locator = FilesystemLocator(mag_dir, abund_dir)

    def get_mag(self, mag_name: str) -> Mag:
        # Get initial data from repositories
        mag_data = self.summary_repo.get_mag_data(mag_name)
        mag_fasta_path = self.locator.find_mag_file(mag_name)

        # Process contigs by reading headers and fetching abundance data
        header_strings = self._read_fasta_headers(mag_fasta_path)
        contig_names = {
            name for h in header_strings if (name := ContigID.parse_name_from_header(h))
        }

        abundance_data = self._get_abundance_data(mag_name, contig_names)
        contig_ids = self._create_contig_ids(header_strings, abundance_data)

        # Create and return the final Mag object
        return Mag(mag_name, mag_data, mag_fasta_path, contig_ids)

    def get_mags_by_query(self,
                          query_string: str,
                          restrict_to: Iterable[str] | None = None) -> list[str]:
        """Filters MAGs using a query string on the summary data."""
        mags = set(self.summary_repo.get_mags_by_query(query_string))
        if restrict_to is not None:
            mags &= set(restrict_to)
        return list(mags)

    def get_hq_mag_names(self) -> list[str]:
        """A convenience shortcut to get high-quality MAGs."""
        return self.summary_repo.get_hq_mag_names()

    def _read_fasta_headers(self, fasta_path: Path) -> list[str]:
        """Reads all header lines from a FASTA file."""
        with open(fasta_path, "r") as f:
            return [line.strip() for line in f if line.startswith(">")]

    def _get_abundance_data(self, mag_name: str, contig_names: set[str]) -> dict[str, dict[str, float]]:
        """Finds and retrieves abundance data for a set of contigs."""
        abund_dir = self.locator.abund_dir
        long_sample = self._get_long_sample_name_from_mag_name(mag_name)
        assembler = self._get_assembler_from_mag_name(mag_name)
        depth_file_path = abund_dir / f"{assembler}__{long_sample}.tsv"

        if not depth_file_path.exists() or not contig_names:
            return {}

        with AbundanceDB(depth_file_path) as abundance_db:
            return abundance_db.get_abund_for_contigs(contig_names)

    def _create_contig_ids(self, header_strings: list[str], all_abund_data: dict) -> list[ContigID]:
        """Creates a list of ContigID objects from headers and abundance data."""
        contigids = []
        for header in header_strings:
            name = ContigID.parse_name_from_header(header)
            # If a contig was in the FASTA but not the abundance file, it gets an empty dict.
            abund_data = all_abund_data.get(name, {}) if name else {}
            contigids.append(ContigID(header, abund_data))
        return contigids

    def get_summary_index(self):
        return self.summary_repo.get_summary_index()

    def _get_long_sample_name_from_mag_name(self, mag_name: str) -> str:
        sample = mag_name.split("_")[0]
        return f"C00{sample[1]}_{sample[2:4]}"

    def _get_assembler_from_mag_name(self, mag_name: str) -> str:
        assembler_char = mag_name.split("_")[1]
        return { "A": "myloasm", "M": "medaka", "P": "proovframe" }.get(assembler_char, "unknown") # Default to 'unknown' if char not found

