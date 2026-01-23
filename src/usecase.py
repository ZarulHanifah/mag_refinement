import os
import sys

print(os.getcwd())
sys.path.append('./src/')

from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path
from pprint import pprint
from typing import Generator, Optional
import subprocess
from string import ascii_lowercase

from tqdm import tqdm

from models import ContigID, Mag, SessionManager


def is_mag_single_circular(mag: Mag) -> bool:
    return len(mag.contigids) == 1 and mag.contigids[0].is_circular()

def is_mag_multiple_contigs_but_got_circular_contig(mag: Mag) -> bool:
    return len(mag.contigids) > 1 and  any([ c.is_circular() for c in mag.contigids ])

def get_mags_single_circular(mags: list[Mag]) -> Generator[Mag, None, None]:
    for mag in mags:
        if is_mag_single_circular(mag):
            yield mag

def get_mags_single_linear(mags: list[Mag]) -> Generator[Mag, None, None]:
    for mag in mags:
        if not is_mag_single_circular(mag):
            yield mag

def get_mags_multiple_contigs_but_with_circular_contig(mags: list[Mag]) -> Generator[Mag, None, None]:
    for mag in mags:
        if is_mag_multiple_contigs_but_got_circular_contig(mag):
            yield mag

def get_circular_contigs_from_mag(mag: Mag) -> Generator[ContigID, None, None]:
    for contig in mag.contigids:
        if contig.is_circular():
            yield contig

def split_ctg_fasta(mag: Mag):
    seqtk_path = "/home/ahbui/Zarul/Software/miniconda3/envs/bio1/bin/seqtk"
    
    print(mag.name)

    idx = 0
    for ctg in get_circular_contigs_from_mag(mag):
        print(ctg.name)
        source_fasta = mag.fp
        target_name = f"{mag.name}{ascii_lowercase[idx]}"
        ctg_path = f"/tmp/{target_name}.id"
        target_path = f"temp_path/{target_name}.fasta"

        with open(ctg_path, "w") as f:
            f.write(ctg.name)

        with open(target_path, "w") as f_out:
            subprocess.run([seqtk_path, "subseq", source_fasta, ctg_path], stdout=f_out)
        idx += 1

def check_mag_depths(mag: Mag, depth_cutoff: int) -> tuple[Mag, bool, bool]:
    is_above_ind = mag.average_coverage_per_sample(mag.long_sample) >= depth_cutoff
    is_above_total = mag.average_coverage_total >= depth_cutoff
    return mag, is_above_ind, is_above_total

def get_mag(mag_name: str, sesh: SessionManager) -> Optional[Mag]:
    """Worker function to get a MAG and check if it's high quality."""
    return sesh.get_mag(mag_name)

def get_hq_mags_and_above_cutoff(sesh: SessionManager,
                                 depth_cutoff: int = 30,
                                 max_workers: int = 4):
    mag_names = sesh.get_hq_mag_names() 

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        get_mag_func = partial(get_mag, sesh=sesh)
        hq_mag_list = list(
            tqdm(
                executor.map(get_mag_func, mag_names),
                total=len(mag_names),
                desc="Loading HQ MAGs",
                unit="MAG",
            )
        )

    check_func = partial(check_mag_depths, depth_cutoff=depth_cutoff)

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = list(
            tqdm(executor.map(check_func, hq_mag_list), total=len(hq_mag_list),
                 desc="Checking MAG quality", unit="MAGs")
        )

    above_ind_cutoff = []
    above_total_cutoff = []
    below_cutoff = []

    for mag, is_above_ind, is_above_total in results:
        if is_above_ind:
            above_ind_cutoff.append(mag)
        if is_above_total:
            above_total_cutoff.append(mag)
        if not is_above_ind and not is_above_total:
            below_cutoff.append(mag)

    return above_ind_cutoff, above_total_cutoff, below_cutoff

if __name__ == "__main__":
    summary_repo = Path('../input_folder/summary.tsv')
    mag_dir = Path('../input_folder/dereplicated_genomes/')
    abund_dir = Path('../input_folder/mtbt_gen_depth/')

    sesh = SessionManager(summary_repo, mag_dir, abund_dir)
    above_ind_cutoff, above_total_cutoff, below_cutoff = get_hq_mags_and_above_cutoff(sesh, 30, 6)

    print("MAGs with multiple contigs but with circular contig:")
    for mag in get_mags_multiple_contigs_but_with_circular_contig(above_ind_cutoff):
        print(f"MAG: {mag.name}, {mag.classification}")
        for contig in get_circular_contigs_from_mag(mag):
            print(f"  - Circular Contig: {contig.name}, Length: {contig.length}")
        split_ctg_fasta(mag)
        print()
