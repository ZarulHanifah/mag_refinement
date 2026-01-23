#!/usr/bin/env python

import sys
import os
import argparse
import pandas as pd

"""
python merge_metabat2_tables.py -t input_table1.tsv input_table2.tsv ... [-o merged_output.tsv] [--outer]
"""

def read_table(df_path, sample_name):
    df = pd.read_csv(df_path, sep="\t", index_col=0)
    df.columns = df.columns.tolist()[:2] + [sample_name] + [f"{sample_name}-var"]
    return df

def main():
    parser = argparse.ArgumentParser(description="Merge multiple tables")
    parser.add_argument("-t", "--tables", required=True, nargs='+', help="Paths to the table TSV files")
    parser.add_argument("-o", "--output", help="Path to the output merged CSV file (prints to stdout if not specified)")
    parser.add_argument("--outer", action='store_true', help="Use outer join instead of inner join")
    args = parser.parse_args()

    try:
        dfs = []
        for df_path in sorted(args.tables):
            sample_name = os.path.dirname(df_path)
            sample_name = os.path.basename(sample_name)
            df = read_table(df_path, sample_name)
            dfs.append(df)
        # dfs = [pd.read_csv(table, sep="\t", index_col=0) for table in sorted(args.tables)]
    except Exception as e:
        print("Error reading input files:", e)
        return

    # Set join type
    join_type = 'outer' if args.outer else 'inner'

    # calculate totalAvgDepth
    avg_dfs = [df.loc[:, ["totalAvgDepth"]] for df in dfs]

    final_avg_df = avg_dfs[0]
    for idx, avg_df in enumerate(avg_dfs[1:]):
        avg_df.columns = [idx]
        final_avg_df = final_avg_df.merge(avg_df, how=join_type, left_index=True, right_index=True)
    final_avg_df = pd.DataFrame(final_avg_df.sum(axis="columns").round(4))
    final_avg_df.columns=["totalAvgDepth"]

    # Merge all tables
    merged_df = dfs[0]
    original_columns_order = merged_df.columns.tolist()
    merged_df.drop(labels=["totalAvgDepth"], axis="columns", inplace=True)
    merged_df = merged_df.merge(final_avg_df, how=join_type, left_index=True, right_index=True)
    merged_df = merged_df.loc[:, original_columns_order]

    for df in dfs[1:]:
        df2 = df.copy()
        df2.drop(labels=["contigLen", "totalAvgDepth"], axis="columns", inplace=True)
        merged_df = merged_df.merge(df2, how=join_type, left_index=True, right_index=True)

    if args.output:
        try:
            merged_df.to_csv(args.output, sep="\t")
            print(f"Merged table saved to {args.output}")
        except Exception as e:
            print("Error saving merged table:", e)
    else:
        print(merged_df.to_csv(sep='\t'))

if __name__ == "__main__":
    main()
