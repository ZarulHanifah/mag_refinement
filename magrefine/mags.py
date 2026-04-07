from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Generator, Optional, Sequence

import pandas as pd

from magrefine.contigids import BaseContigID, MyloContigID


class BaseMag(ABC):
    """
    Abstract Base Class for a Metagenome-Assembled Genome (MAG).
    Provides common trajectory, metadata parsing, and relationship logic.
    """
    name: str
    _fp: Path
    contigids: Sequence[BaseContigID]
    completeness: float
    contamination: float
    total_contigs: int

    @property
    def fp(self) -> str: 
        return str(self._fp)

    @property
    def sample(self) -> str: 
        return self.name.split("_")[0]

    @property
    def long_sample(self) -> str: 
        return f"C00{self.sample[1]}_{self.sample[2:4]}"

    @property
    def binner(self) -> str: 
        return self.name.split("_")[2].split(".")[0]

    @property
    def assembler(self) -> str: 
        return { "A": "myloasm", "M": "medaka", "P": "proovframe" }[self.name.split("_")[1]]

    # Subclasses must implement completeness, contamination, and total_contigs

    def is_high_quality(self) -> bool:
        return self.completeness >= 90 and self.contamination <= 5

    def compare_to(self, other: BaseMag) -> dict[str, float]:
        """
        Compare this MAG's quality stats to another MAG (e.g. comparing Refined to Initial).
        Returns the straight difference.
        Meaning:
         - positive completeness_diff means this MAG is MORE complete.
         - negative contamination_diff means this MAG is LESS contaminated (better).
         - negative total_contigs_diff means this MAG is MORE contiguous (better).
        """
        return {
            "completeness_diff": round(self.completeness - other.completeness, 2),
            "contamination_diff": round(self.contamination - other.contamination, 2),
            "total_contigs_diff": round(self.total_contigs - other.total_contigs, 2)
        }

    @property
    def children(self) -> list[BaseMag]:
        if not hasattr(self, "_children"):
            self._children = []
        return self._children

    def add_child(self, child: BaseMag):
        self.children.append(child)

    def tree_report(self) -> str:
        """
        Generates a top-down tree showing this MAG and all its descendant derivations.
        """
        lines = []
        straight_down = "│"

        def recurse(mag: BaseMag, prefix: str, is_last: bool):
            comp, contam, contigs = mag.completeness, mag.contamination, mag.total_contigs

            if not prefix:
                lines.append(f"└── {mag.name} (Comp: {comp:.2f}% | Contam: {contam:.2f}% | Contigs: {contigs})")
                child_prefix = "    "
            else:
                branch = "└──" if is_last else "├──"
                lines.append(f"{prefix}{branch} {mag.name} (Comp: {comp:.2f}% | Contam: {contam:.2f}% | Contigs: {contigs})")
                child_prefix = prefix + ("    " if is_last else "│   ")

            for idx, child in enumerate(mag.children):
                child_is_last = (idx == len(mag.children) - 1)

                diffs = child.compare_to(mag)
                comp_diff = diffs["completeness_diff"]
                contam_diff = diffs["contamination_diff"]
                contigs_diff = diffs["total_contigs_diff"]

                comp_color = "green" if comp_diff > 0 else ("red" if comp_diff < 0 else "white")
                contam_color = "red" if contam_diff > 0 else ("green" if contam_diff < 0 else "white")
                contigs_color = "red" if contigs_diff > 0 else ("green" if contigs_diff < 0 else "white")

                comp_str = f"[{comp_color}]{comp_diff:+.2f}%[/{comp_color}]"
                contam_str = f"[{contam_color}]{contam_diff:+.2f}%[/{contam_color}]"
                contigs_str = f"[{contigs_color}]{contigs_diff:+d}[/{contigs_color}]"

                lines.append(f"{child_prefix}{straight_down}")
                lines.append(f"{child_prefix}↓ (comp: {comp_str}, contam: {contam_str}, contigs: {contigs_str})")
                recurse(child, child_prefix, child_is_last)

        recurse(self, "", True)
        return "\n".join(lines)

    def lineage_report(self) -> str:
        """
        Generates a tree-like string showing the history of improvements 
        from the original root MAG to this current MAG.
        """
        chain = []
        current = self
        while current is not None:
            chain.append(current)
            current = getattr(current, "parent", None)

        chain.reverse()

        lines = []
        for i, mag in enumerate(chain):
            indent = "    " * i

            comp = mag.completeness
            contam = mag.contamination
            contigs = mag.total_contigs
            mag_str = f"{indent}└── {mag.name} (Comp: {comp:.2f}% | Contam: {contam:.2f}% | Contigs: {contigs})"

            if i > 0:
                prev_mag = chain[i - 1]
                diffs = mag.compare_to(prev_mag)
                comp_diff = diffs["completeness_diff"]
                contam_diff = diffs["contamination_diff"]
                contigs_diff = diffs["total_contigs_diff"]

                comp_color = "green" if comp_diff > 0 else ("red" if comp_diff < 0 else "white")
                contam_color = "red" if contam_diff > 0 else ("green" if contam_diff < 0 else "white")
                contigs_color = "red" if contigs_diff > 0 else ("green" if contigs_diff < 0 else "white")

                comp_str = f"[{comp_color}]{comp_diff:+.2f}%[/{comp_color}]"
                contam_str = f"[{contam_color}]{contam_diff:+.2f}%[/{contam_color}]"
                contigs_str = f"[{contigs_color}]{contigs_diff:+d}[/{contigs_color}]"

                arrow_indent = "    " * i
                delta_str = f"{arrow_indent}↓ (comp: {comp_str}, contam: {contam_str}, contigs: {contigs_str})"
                lines.append(delta_str)

            lines.append(mag_str)

        return "\n".join(lines)

    def get_depth_df(self) -> pd.DataFrame:
        """
        Generates a pandas DataFrame of contig depths for all samples.
        """
        depths = [ contigid.get_abund_info() for contigid in self.contigids ]
        indices = [contigid.name for contigid in self.contigids]
        return pd.DataFrame(depths, index=pd.Index(indices))

    def get_depth_per_contig(self, contigid: BaseContigID, samplename=None):
        if samplename:
            return contigid.get_abund_info()[samplename]
        return contigid.depth_from_all_samples

    @property
    def average_coverage_total(self) -> float:
        total_length = 0
        total_bases = 0
        for c in self.contigids:
            length = c.length or 0  # Fallback gracefully
            total_bases += c.depth_from_all_samples * length
            total_length += length
        if total_length == 0: return 0.0
        return total_bases / total_length

    def average_coverage_per_sample(self, samplename: Optional[str] = None ) -> float:
        if samplename == None:
            samplename = self.long_sample
        values = [ (contigid.name, contigid.length or 0, self.get_depth_per_contig(contigid, samplename) )
                  for contigid in self.contigids ]
        total_bases = 0
        for v in values:
            total_bases += v[1] * v[2]

        total_len = sum([v[1] for v in values])
        if total_len == 0: return 0.0
        return total_bases / total_len


@dataclass
class Mag(BaseMag):
    """
    Original MAG representation, capturing summary metrics explicitly.
    """
    name: str
    _fp: Path
    contigids: Sequence[BaseContigID]
    completeness: float
    contamination: float
    total_contigs: int

    contig_n50: int | float = 0
    gc_content: float = 0.0
    max_contig_length: int = 0
    genome_size: int = 0
    average_coverage: float = 0.0
    trna_counts: int = 0
    classification: str = ""
    red_value: float = 0.0
    _has_16s: bool = False
    _has_23s: bool = False
    _has_5s: bool = False

    @classmethod
    def from_summary_data(cls, name: str, data: dict, fp: Path, contigids: Sequence[BaseContigID]) -> "Mag":
        """
        Alternate constructor to instantiate a Mag directly from a CheckM summary dictionary.
        """
        return cls(
            name=name,
            _fp=fp,
            contigids=contigids,
            completeness=float(data['Completeness']),
            contamination=float(data['Contamination']),
            total_contigs=int(data['Total_Contigs']),
            contig_n50=data['Contig_N50'],
            gc_content=float(data['GC_Content']),
            max_contig_length=data['Max_Contig_Length'],
            genome_size=data['genome_size'],
            average_coverage=float(data.get('average_coverage', 0.0)),
            trna_counts=data.get('tRNA counts', 0),
            classification=data.get('classification', ''),
            red_value=data.get('red_value', 0.0),
            _has_16s=bool(data.get('16S_rRNA', False)),
            _has_23s=bool(data.get('23S_rRNA', False)),
            _has_5s=bool(data.get('5S_rRNA', False))
        )

    def has_16s_rrna(self): return self._has_16s
    def has_23s_rrna(self): return self._has_23s
    def has_5s_rrna(self):  return self._has_5s

    def __repr__(self): 
        return (
            f"{self.__class__.__name__}"
            f"(name={self.name}, "
            f"completeness={self.completeness}, "
            f"contamination={self.contamination}, "
            f"assembler={self.assembler}, "
            f"sample={self.sample}, "
            f"binner={self.binner}, "
            f"total_contigs={self.total_contigs})"
        )


@dataclass
class RefinedMag(BaseMag):
    """
    A MAG generated via refinement operations (e.g. scaffolding) whose metrics are explicitly given.
    """
    name: str
    _fp: Path
    contigids: Sequence[BaseContigID]
    completeness: float
    contamination: float
    parent: Optional[BaseMag] = None
    total_contigs: int = field(init=False)

    def __post_init__(self):
        self.total_contigs = len(self.contigids)
        if self.parent is not None:
            self.parent.add_child(self)

    def improvement_report(self) -> dict[str, float] | None:
        """
        Returns stats relative to the parent MAG, if a parent was assigned.
        """
        if self.parent is None:
            return None
        return self.compare_to(self.parent)

    @staticmethod
    def read_checkm_df(mag_name, checkmqual: Path):
        df = pd.read_csv(checkmqual, sep="\t", index_col=0)
        df = df.loc[[mag_name], :]
        df_dict =  df.to_dict(orient="index")
        return df_dict[mag_name]

    @staticmethod
    def extract_header_from_fasta_file(fasta_file: Path) -> Generator[str, Any, Any]:
        with open(fasta_file, "r") as f:
            for line in f:
                if line.startswith(">"):
                    yield line.removeprefix(">")

    @classmethod
    def get_mylocontigid_from_fasta_file(cls, fasta_file: Path):
        return [
            MyloContigID(contig_header)
            for contig_header in cls.extract_header_from_fasta_file(fasta_file)
        ]

    @classmethod
    def from_checkmqual(cls,
                        mag_name: str,
                        _fp: Path,
                        checkmqual: Path,
                        parent: Optional[BaseMag] = None) -> RefinedMag:
        checkm_dict = RefinedMag.read_checkm_df(mag_name, checkmqual)
        # how checkm_dict looks like:
        # {'Completeness': 77.4, 'Contamination': 7.54, 'Completeness_Model_Used': 'Gradient Boost (General Model)', 'T ranslation_Table_Used': 11, 'Coding_Density': 0.815, 'Contig_N50': 31122, 'Average_Gene_Length': 218.8590924599113, 'Genome_Size': 23 47852, 'GC_Content': 0.47, 'Total_Coding_Sequences': 2931, 'Total_Contigs': 145, 'Max_Contig_Length': 154524, 'Additional_Notes': nan }

        contigids = RefinedMag.get_mylocontigid_from_fasta_file(_fp)
        return cls(
            mag_name,
            _fp,
            contigids,
            checkm_dict["Completeness"],
            checkm_dict["Contamination"],
            parent
        )

    def __repr__(self): 
        return (
            f"{self.__class__.__name__}"
            f"(name={self.name}, "
            f"completeness={self.completeness}, "
            f"contamination={self.contamination}, "
            f"assembler={self.assembler}, "
            f"sample={self.sample}, "
            f"binner={self.binner}, "
            f"total_contigs={self.total_contigs})"
        )


if __name__ == "__main__":
    fp = Path("tests/input_folder/dereplicated_genomes/C1E5_M_metabat.1297.fasta")
    checkmqual = Path("tests/input_folder/checkm2/C1E5_M_metabat.1297.fasta/quality_report.tsv")
    root = RefinedMag.from_checkmqual("C1E5_M_metabat.1297", fp, checkmqual, None)
    childF1 = RefinedMag.from_checkmqual("C1E5_M_metabat.1297", fp, checkmqual, root)
    
