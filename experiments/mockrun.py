from pathlib import Path

from rich import print

from magrefine.mags import RefinedMag

if __name__ == "__main__":
    fp = Path("tests/input_folder/dereplicated_genomes/C1E5_M_metabat.1297.fasta")
    checkmqual = Path("tests/input_folder/checkm2/C1E5_M_metabat.1297/quality_report.tsv")
    checkmqualdir = Path("tests/input_folder/checkm2/C1E5_M_metabat.1297")

    root = RefinedMag.from_checkm2qual("C1E5_M_metabat.1297", fp, checkmqual, None)
    root.name = "root"
    root.completeness = 50.0
    root.contamination = 10.0

    childF1 = RefinedMag.from_checkm2qual("C1E5_M_metabat.1297", fp, checkmqual, root)
    childF1.name = "childF1"

    childF2 = RefinedMag.from_checkm2qual("C1E5_M_metabat.1297", fp, checkmqualdir, root)
    childF2.name = "childF2"
    childF2.contamination = 15.0

    print("IMPROVEMENT REPORT")
    print(childF1.improvement_report())
    print()

    print("LINEAGE REPORT")
    print(childF1.lineage_report())
    print()

    print("TREE REPORT")
    print(root.tree_report())
    print()

