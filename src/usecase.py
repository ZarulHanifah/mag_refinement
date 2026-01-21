from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path
from pprint import pprint
from typing import Optional

from tqdm import tqdm

from models import Mag, SessionManager


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

def main():
    summary_repo = Path('../input_folder/summary.tsv')
    mag_dir = Path('../input_folder//dereplicated_genomes/')
    abund_dir = Path('../input_folder/mtbt_gen_depth/')

    sesh = SessionManager(summary_repo, mag_dir, abund_dir)
    above_ind_cutoff, above_total_cutoff, below_cutoff = get_hq_mags_and_above_cutoff(sesh, 30, 6)

    for i in [above_ind_cutoff, above_total_cutoff, below_cutoff]:
        i = [m.name for m in i]
        varname = f"{i=}".split("=")[0]
        print(varname.center(30, "="))
        pprint(i)
        print()

if __name__ == "__main__":
    main()
