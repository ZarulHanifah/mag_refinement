# execute_workflow.py
import argparse
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path
from typing import List, Optional

import snakemake
import yaml
from tqdm import tqdm

# This allows the script to find your 'models.py' file in the 'src' directory.
sys.path.append(str(Path(__file__).parent / "magrefine"))
from magrefine.models import Mag
from magrefine.sessionmanager import SessionManager


def get_mag_worker(mag_name: str, session: SessionManager) -> Mag:
    """Helper function to instantiate a single Mag object for parallel processing."""
    return session.get_mag(mag_name)


def apply_custom_filters(mag: Mag) -> bool:
    """
    Apply custom filtering logic to a single Mag object.
    <<< THIS IS THE SECTION YOU WILL MODIFY OFTEN FOR YOUR FILTERING CRITERIA. >>>
    Returns True if the MAG passes the filters, False otherwise.
    """
    # FIXIT:  CHANGE HERE!
    # Example 1: Filter by substring in classification
    if "Asgard" in mag.classification:
        return True

    # Example 2: Filter based on a computed property and other criteria
    # if mag.assembler == "myloasm" and mag.average_coverage_total > 15:
    #     return True

    # Example 3: Filter on rRNA presence (uncomment and adjust as needed)
    # if mag.has_16s_rrna() and mag.completeness > 90:
    #     return True

    # If no criteria above are met, the MAG does not pass the filters
    return False


def select_mags(session: SessionManager, 
                completeness: int = 50,
                contamination: int = 10,
                max_workers: int = 4) -> list[str]:
    """
    This function contains all the logic for selecting which MAGs to process.
    
    <<< THIS IS THE SECTION YOU WILL MODIFY OFTEN. >>>

    It returns a final list of MAG names to be processed by the workflow.
    """
    print("Starting MAG selection...")
    

    # --- Step 1: Pre-filter using a fast query (Optional, but Recommended) ---
    # FIXIT:  CHANGE HERE!
    pre_filter_query = (
            f"Completeness >= {completeness}"
            f" and Contamination <= {contamination}"
            f" and classification.str.contains('Asgard')"
    )

    candidate_mag_names = session.get_mags_by_query(pre_filter_query)
    print(f"Found {len(candidate_mag_names)} candidates after pre-filtering with query: '{pre_filter_query}'")

    # --- Step 2: Instantiate full Mag objects in parallel (Fast and Scalable) ---
    print(f"Loading {len(candidate_mag_names)} full Mag objects in parallel...")
    mag_objects = []
    # Use max_workers=None to default to the number of processors on your machine
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Create a partial function to pass the session object to each worker
        worker_func = partial(get_mag_worker, session=session)
        
        mag_objects = list(
            tqdm(
                executor.map(worker_func, candidate_mag_names),
                total=len(candidate_mag_names),
                desc="Loading Mag objects"
            )
        )

    # --- Step 3: Post-filter on the in-memory Mag objects (Flexible and Powerful) ---
    final_mag_names = []
    print("Applying flexible post-filters on in-memory objects...")

    for mag in tqdm(mag_objects, desc="Applying custom filters"):
        if apply_custom_filters(mag):
            final_mag_names.append(mag.name)

    print(f"Found {len(final_mag_names)} MAGs that passed the final flexible filters.")
    return sorted(list(set(final_mag_names))) # Return unique, sorted list


def run_workflow(mags_to_process: list[str],
                 chunk_size: int,
                 base_config: dict,
                 njobs: int,
                 dry_run: bool,
                 extra_snakemake_args: List[str]):
    if not mags_to_process:
        print("The filtering criteria resulted in an empty list of MAGs. Nothing to run.")
        return

    chunks = [
        mags_to_process[i:i + chunk_size]
        for i in range(0, len(mags_to_process), chunk_size)
    ]
    print(f"\nProcessing {len(mags_to_process)} MAGs in {len(chunks)} chunks of up to {chunk_size} each.")

    for i, chunk in enumerate(chunks, 1):
        print("-" * 80)
        print(f"Starting chunk {i}/{len(chunks)} with {len(chunk)} MAGs...")

        # Write temp config for this chunk
        snake_config = dict(base_config)
        snake_config["mags_to_process"] = chunk

        snake_config_path = Path("config_snakefile.yaml")
        with open(snake_config_path, "w") as f:
            yaml.safe_dump(snake_config, f)

        # Build the base snakemake command
        base_cmd = [
            "snakemake",
            "--snakefile", "Snakefile",
            "--configfile", str(snake_config_path),
        ]
        
        # 1. Unlock the directory
        unlock_cmd = base_cmd + ["--unlock"]
        print(f"Unlocking directory... {' '.join(unlock_cmd)}")
        subprocess.run(unlock_cmd, check=True)

        # 2. Build the main execution command
        run_cmd = base_cmd + [
            "-k",
            "--rerun-trigger", "mtime",
            "--jobs", str(njobs),
            "--use-conda",
            "--cluster-config", base_config["cluster_config"],
            "--cluster",
            "sbatch --cpus-per-task={threads} "
            "--output={cluster.output} --error={cluster.error} "
            "--job-name={cluster.jobname} {cluster.etc}",
            "--rerun-incomplete",
        ] + extra_snakemake_args

        if dry_run:
            run_cmd.append("--dryrun")

        print(f"Executing command: {' '.join(run_cmd)}")
        result = subprocess.run(run_cmd)
        if result.returncode != 0:
            print(f"Workflow failed on chunk {i}. Aborting.", file=sys.stderr)
            sys.exit(1)

        print(f"Finished chunk {i}/{len(chunks)} successfully.")

    print("-" * 80)
    print("All chunks processed successfully!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Flexible Snakemake workflow orchestrator for MAGs.",
        epilog="Any unrecognized arguments will be passed directly to Snakemake."
    )
    # General arguments
    parser.add_argument("--configfile", default="config/config.yaml", help="Path to the base YAML configuration file.")
    parser.add_argument("--list-mags", action="store_true", help="List the selected MAGs without running the workflow.")
    
    # Filtering arguments
    parser.add_argument("--completeness", type=int, default=50, help="Minimum completeness for pre-filtering.")
    parser.add_argument("--contamination", type=int, default=10, help="Maximum contamination for pre-filtering.")

    # Execution arguments
    parser.add_argument("-j", "--jobs", type=int, help="Number of parallel jobs for Snakemake (overrides config).")
    parser.add_argument("--chunk-size", type=int, default=10, help="Number of MAGs to process per Snakemake run.")
    parser.add_argument("-n", "--dry-run", action="store_true", help="Perform a dry run (passes --dryrun to Snakemake).")

    # Path override arguments
    parser.add_argument("--summary-path", type=str, help="Path to summary TSV file (overrides config).")
    parser.add_argument("--mag-dir", type=str, help="Path to MAGs directory (overrides config).")
    parser.add_argument("--abund-dir", type=str, help="Path to abundance directory (overrides config).")

    args, snakemake_args = parser.parse_known_args()

    # --- Initialize ---
    base_config = yaml.safe_load(Path(args.configfile).read_text())

    # Determine paths: CLI args override config file
    summary_path_str = args.summary_path if args.summary_path else base_config["summary_path"]
    mag_dir_str = args.mag_dir if args.mag_dir else base_config["mag_dir"]
    abund_dir_str = args.abund_dir if args.abund_dir else base_config["abund_dir"]

    summary_path = Path(summary_path_str)
    if not summary_path.is_file():
        print(f"Error: Summary file not found at '{summary_path}'", file=sys.stderr)
        sys.exit(1)

    session = SessionManager(
        summary_path=summary_path, 
        mag_dir=Path(mag_dir_str), 
        abund_dir=Path(abund_dir_str)
    )

    # --- Configuration for the Workflow Run ---
    njobs = args.jobs if args.jobs is not None else base_config.get('njobs', 20)
    
    # 1. Select MAGs using your custom logic
    selected_mags = select_mags(session, completeness=args.completeness, contamination=args.contamination)

    # 2. List MAGs or run the workflow
    if args.list_mags:
        for mag_name in selected_mags:
            print(mag_name)
    else:
        run_workflow(
            mags_to_process=selected_mags,
            chunk_size=args.chunk_size,
            base_config=base_config,
            njobs=njobs,
            dry_run=args.dry_run,
            extra_snakemake_args=snakemake_args
        )
