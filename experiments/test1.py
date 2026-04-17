from pathlib import Path

import os

from magrefine.mags import Mag, RefinedMag
from magrefine.sessionmanager import SessionManager

from rich import print

summary_path = Path("./input_folder/summary.tsv")
mag_dir_path = Path("./input_folder/dereplicated_genomes/")
abund_dir_path = Path("./input_folder/mtbt_gen_depth/")

sesh = SessionManager(
    summary_path=summary_path,
    mag_dir=Path(mag_dir_path),
    abund_dir=Path(abund_dir_path)
)

# mag_name = "C1E5_M_metabat.1297"

for mag_name in os.listdir("results/flye_fq/"):
    root = sesh.get_mag(mag_name)

    assemblers = [
        "flye_fq",
        "hifiasm_fq",
        "myloasm_fq",
        "longstitch"
    ]

    for assem in assemblers:
        short_assem = assem.split("_")[0]
        assem_fp = Path(f"results/{assem}/{mag_name}/{mag_name}.{short_assem}.fasta")
        assem_checkmqual = Path(f"results/checkm2/{assem}/{mag_name}/quality_report.tsv")
        childF1 = RefinedMag.from_checkm2qual(mag_name, assem_fp, assem_checkmqual, root)
        childF1.name = assem

    print(mag_name, root.classification)
    print(root.tree_report())
