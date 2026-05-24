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
mag_names = os.listdir("results/flye/")

for mag_name in mag_names:
    root = sesh.get_mag(mag_name)

    assemblers = [
        "flye",
        "hifiasm",
        "myloasm",
        "longstitch",
        # "medaka",
    ]

    for assem in assemblers:
        childF1 = RefinedMag.from_checkm2qual(
          mag_name,
          Path(f"results/{assem}/{mag_name}/{mag_name}.{assem}.fasta"),
          Path(f"results/checkm2/{assem}/{mag_name}/quality_report.tsv"),
          root
        )
        childF1.name = assem

    # racon
    assem = "racon"
    medaka_rd1 = RefinedMag.from_checkm2qual(
      mag_name,
      Path(f"results/{assem}/{mag_name}.fasta"),
      Path(f"results/checkm2/{assem}/{mag_name}/quality_report.tsv"),
      root,
    )
    medaka_rd1.name = "racon"

    # medaka rd1
    assem = "medaka"
    medaka_rd1 = RefinedMag.from_checkm2qual(
      mag_name,
      Path(f"results/{assem}/{mag_name}.fasta"),
      # Path(f"results/{assem}/{mag_name}/{mag_name}.{assem}.fasta"),
      Path(f"results/checkm2/{assem}/{mag_name}/quality_report.tsv"),
      root,
    )
    medaka_rd1.name = "medaka_rd1"

    # medaka rd1 --> proovframe
    assem = "proovframe"
    proovframe = RefinedMag.from_checkm2qual(
      mag_name,
      Path(f"results/proovframe/assem/{mag_name}.fasta"),
      Path(f"results/checkm2/proovframe/{mag_name}/quality_report.tsv"),
      medaka_rd1,
    )
    proovframe.name = "proovframe"

    # medaka rd1 --> proovframe
    assem = "medaka_rd2"
    medaka_rd2 = RefinedMag.from_checkm2qual(
      mag_name,
      Path(f"results/{assem}/{mag_name}.fasta"),
      Path(f"results/checkm2/{assem}/{mag_name}/quality_report.tsv"),
      medaka_rd1,
    )
    medaka_rd2.name = "medaka_rd2"

    # print(mag_name, root.classification)
    print(mag_name)
    print(root.tree_report())
