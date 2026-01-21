# create_manifest_from_directory.py
import argparse
import sys
from pathlib import Path
import pandas as pd
import csv

def get_abundance_path(mag_name: str, abund_dir: Path) -> str:
    """
    Replicates the logic from SessionManager to determine the correct
    abundance file path based on the original MAG naming convention.
    
    Example MAG name: C1A3_A_metabat.872
    - C1A3 -> Sample
    - A -> Assembler (A: myloasm, M: medaka, P: proovframe)
    """
    try:
        parts = mag_name.split('_')
        sample_short = parts[0]
        assembler_char = parts[1]

        # Reconstruct long sample name (e.g., C1A3 -> C001_A3)
        long_sample = f"C00{sample_short[1]}_{sample_short[2:4]}"
        
        assembler_map = {"A": "myloasm", "M": "medaka", "P": "proovframe"}
        assembler = assembler_map.get(assembler_char)

        if not assembler:
            print(f"Warning: Unknown assembler character '{assembler_char}' for MAG '{mag_name}'. Skipping abundance path.", file=sys.stderr)
            return "NA"

        # Construct the expected path
        expected_path = abund_dir / f"{assembler}__{long_sample}.tsv"
        
        # Return the path as a string if it exists
        if expected_path.exists():
            return str(expected_path)
        else:
            print(f"Warning: Could not find abundance file for MAG '{mag_name}' at expected path: {expected_path}", file=sys.stderr)
            return "NA"

    except IndexError:
        print(f"Warning: Could not parse MAG name '{mag_name}' for abundance info. Skipping abundance path.", file=sys.stderr)
        return "NA"

def create_manifest(mag_dir: Path, abund_dir: Path, summary_path: Path, output_path: Path):
    """
    Generates a manifest file from existing MAG directories and a summary file.
    """
    print(f"Reading summary file from: {summary_path}")
    if not summary_path.is_file():
        print(f"Error: Summary file not found at '{summary_path}'", file=sys.stderr)
        sys.exit(1)
        
    summary_df = pd.read_csv(summary_path, sep="\t", index_col=0)
    mag_names = summary_df.index.tolist()
    print(f"Found {len(mag_names)} MAGs in summary file.")

    manifest_data = []

    for name in mag_names:
        fasta_path = mag_dir / f"{name}.fasta"
        if not fasta_path.is_file():
            print(f"Warning: FASTA file for '{name}' not found at '{fasta_path}'. Skipping.", file=sys.stderr)
            continue
            
        abundance_path = get_abundance_path(name, abund_dir)

        manifest_data.append({
            "mag_id": name,          # Default mag_id is the original name
            "parent_id": "NA",       # No parent for initial MAGs
            "fasta_path": str(fasta_path),
            "abundance_path": abundance_path,
            "summary_name": name,    # Links back to the row in the summary file
        })

    if not manifest_data:
        print("No valid MAGs found to create a manifest.", file=sys.stderr)
        sys.exit(1)

    print(f"Writing manifest for {len(manifest_data)} MAGs to: {output_path}")
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=manifest_data[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(manifest_data)

    print("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create a MAG manifest file from a directory of MAGs and a summary file."
    )
    parser.add_argument(
        "--summary-file", 
        type=Path, 
        required=True,
        help="Path to the summary.tsv file from CheckM2."
    )
    parser.add_argument(
        "--mag-dir", 
        type=Path, 
        required=True,
        help="Path to the directory containing MAG .fasta files."
    )
    parser.add_argument(
        "--abund-dir", 
        type=Path, 
        required=True,
        help="Path to the directory containing abundance .tsv files."
    )
    parser.add_argument(
        "-o", "--output", 
        type=Path, 
        default="mags_manifest.tsv",
        help="Path for the output manifest TSV file."
    )
    
    args = parser.parse_args()

    create_manifest(
        mag_dir=args.mag_dir,
        abund_dir=args.abund_dir,
        summary_path=args.summary_file,
        output_path=args.output
    )
