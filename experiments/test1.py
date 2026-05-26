from pathlib import Path

import os

from magrefine.mags import Mag, RefinedMag
from magrefine.graph import RefinementGraph
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
    graph = RefinementGraph()
    root = sesh.get_mag(mag_name)
    graph.add_mag(root)

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
          Path(f"results/checkm2/{assem}/{mag_name}/quality_report.tsv")
        )
        childF1.name = assem
        graph.add_derivation(parent=root, child=childF1)

    # racon
    assem = "racon"
    racon_polish = RefinedMag.from_checkm2qual(
      mag_name,
      Path(f"results/{assem}/{mag_name}.fasta"),
      Path(f"results/checkm2/{assem}/{mag_name}/quality_report.tsv")
    )
    racon_polish.name = "racon"
    graph.add_derivation(parent=root, child=racon_polish)

    # medaka rd1
    assem = "medaka"
    medaka_rd1 = RefinedMag.from_checkm2qual(
      mag_name,
      Path(f"results/{assem}/{mag_name}.fasta"),
      Path(f"results/checkm2/{assem}/{mag_name}/quality_report.tsv")
    )
    medaka_rd1.name = "medaka_rd1"
    graph.add_derivation(parent=root, child=medaka_rd1)

    # medaka rd1 --> proovframe
    assem = "proovframe"
    proovframe = RefinedMag.from_checkm2qual(
      mag_name,
      Path(f"results/proovframe/assem/{mag_name}.fasta"),
      Path(f"results/checkm2/proovframe/{mag_name}/quality_report.tsv")
    )
    proovframe.name = "proovframe"
    graph.add_derivation(parent=medaka_rd1, child=proovframe)

    # medaka rd1 --> medaka rd2
    assem = "medaka_rd2"
    medaka_rd2 = RefinedMag.from_checkm2qual(
      mag_name,
      Path(f"results/{assem}/{mag_name}.fasta"),
      Path(f"results/checkm2/{assem}/{mag_name}/quality_report.tsv")
    )
    medaka_rd2.name = "medaka_rd2"
    graph.add_derivation(parent=medaka_rd1, child=medaka_rd2)

    # print(mag_name, root.classification)
    print(mag_name)
    print(graph.tree_report(root))
