import sys
from collections.abc import Iterable
from pathlib import Path

import pandas as pd

from .models import AbundanceDB, ContigID, Mag


class SummaryRepository:
    """
    Loads and provides access to the main summary data file, which contains
    statistics for all MAGs.
    """
    def __init__(self, _summary_path: Path):
        self.summary = pd.read_csv(_summary_path, sep="\t", index_col=0)
        # Add MAG attributes as columns for querying
        self.summary['binner'] = self.summary.index.map(lambda x: x.split("_")[2].split(".")[0])
        self.summary['assembler'] = self.summary.index.map(self._get_assembler_from_mag_name)

    @staticmethod
    def _get_assembler_from_mag_name(mag_name: str) -> str:
        assembler_char = mag_name.split("_")[1]
        return { "A": "myloasm", "M": "medaka", "P": "proovframe" }.get(assembler_char, "unknown") # Default to 'unknown' if char not found

    def get_mag_data(self, mag_name: str) -> dict:
        return self.summary.loc[mag_name, :].to_dict()

    def get_summary_index(self) -> list[str]:
        return self.summary.index.tolist()    # pyright: ignore

    def get_mags_by_query(self, query_string: str) -> list[str]:
        """
        Filters MAGs based on a pandas query string.

        Args:
            query_string: A string that is a valid pandas query.

        Returns:
            A list of MAG names that match the query.
        """
        try:
            return self.summary.query(query_string).index.tolist()    # pyright: ignore
        except Exception as e:
            print(f"Error executing query: {e}", file=sys.stderr)
            return []

    def get_hq_mag_names(self) -> list[str]:
        """A convenience shortcut to get high-quality MAGs."""
        hq_query = "Completeness >= 90 and Contamination <= 5"
        return self.get_mags_by_query(hq_query)

class FilesystemLocator:
    """
    A helper class to locate MAG and abundance files on the filesystem.
    """
    def __init__(self, mag_dir: Path, abund_dir: Path) -> None:
        self.mag_dir = mag_dir        # Directory containing MAG FASTA files.
        self.abund_dir = abund_dir    # Directory containing abundance TSV files.

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

