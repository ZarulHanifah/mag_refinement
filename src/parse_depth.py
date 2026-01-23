#!/usr/bin/env python

import sys
import argparse
import pandas as pd


def main():
    # Argument parser setup
    parser = argparse.ArgumentParser(
        description="Filter and transpose a table based on specified MAG names."
    )
    parser.add_argument(
        "-i", "--input-df-path", required=True, help="Path to the input depth table (tab-separated)."
    )
    parser.add_argument(
        "-n", "--names", nargs="+", required=True, help="List of MAG names to filter."
    )
    parser.add_argument(
        "-o", "--output", default=None, help="Path to save the output table (default: stdout)."
    )
    args = parser.parse_args()

    df_path = args.input_df_path
    names = args.names
    output = args.output

    # Input file validation
    try:
        df = pd.read_csv(df_path, sep="\t", index_col=0)
    except FileNotFoundError:
        sys.stderr.write(f"Error: Input file '{df_path}' not found.\n")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        sys.stderr.write(f"Error: Input file '{df_path}' is empty or invalid.\n")
        sys.exit(1)

    # Filter the DataFrame
    missing_names = [name for name in names if name not in df.index]
    if missing_names:
        sys.stderr.write(f"Warning: The following names were not found in the table: {missing_names}\n")

    filtered_df = df.loc[df.index.intersection(names), :]

    if filtered_df.empty:
        sys.stderr.write("Error: None of the provided names were found in the input table.\n")
        sys.exit(1)

    # Transpose the DataFrame
    transposed_df = filtered_df.T

    # Output the result
    if output:
        transposed_df.to_csv(output, sep="\t")
    else:
        transposed_df.to_csv(sys.stdout, sep="\t")


if __name__ == "__main__":
    main()

