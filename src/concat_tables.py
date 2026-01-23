#!/usr/bin/env python

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Concatenate two tables with the same headers")
    parser.add_argument("-i", "--input", nargs="+", required=True, help="Paths to the input table CSV files")
    parser.add_argument("-o", "--output", required=True, help="Path to the output concatenated CSV file")
    args = parser.parse_args()

    input_files = args.input
    output_file = args.output

    try:
        dfs = [pd.read_csv(file, sep = "\t", index_col = 0) for file in input_files]
    except Exception as e:
        print("Error reading input files:", e)
        return

    concatenated_df = pd.concat(dfs)
    
    try:
        concatenated_df.to_csv(output_file, sep = "\t")
        print("Concatenated table saved to", output_file)
    except Exception as e:
        print("Error saving concatenated table:", e)

if __name__ == "__main__":
    main()
