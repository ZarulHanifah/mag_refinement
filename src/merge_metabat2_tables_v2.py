#!/usr/bin/env python

import sys
import os
import argparse
import pandas as pd
import numpy as np
import logging

def read_table(df_path, sample_name):
    df = pd.read_csv(df_path, sep="\t", index_col=0)
    if df.empty:
        raise ValueError(f"Input file {df_path} is empty")
    df.columns = df.columns.tolist()[:2] + [sample_name] + [f"{sample_name}-var"]
    return df

def detect_outliers_mag_pattern(df, threshold=0.5):
    coverage_columns = [col for col in df.columns if not col.endswith('-var') and col not in ['contigLen', 'totalAvgDepth']]
    coverage_data = df[coverage_columns]

    # Calculate median coverage pattern for the MAG
    mag_pattern = coverage_data.median()

    # Calculate relative deviation for each contig
    relative_deviation = (coverage_data - mag_pattern).abs() / (mag_pattern + 1e-10)  # Add small value to avoid division by zero

    # Calculate mean deviation across all samples for each contig
    mean_deviation = relative_deviation.mean(axis=1)

    # Identify outliers
    outliers = mean_deviation > threshold

    return df[~outliers], df[outliers]

def main():
    parser = argparse.ArgumentParser(description="Merge multiple tables")
    parser.add_argument("-i", "--input-tables", required=True, nargs='+', help="Paths to the input table TSV files")
    parser.add_argument("-o", "--output", help="Path to the output merged CSV file (prints to stdout if not specified)")
    parser.add_argument("--inner", action='store_true', help="Use inner join instead of outer join")
    parser.add_argument("--keep-outliers", action='store_true', help="Keep contigs with outlier coverage values")
    parser.add_argument("-t", "--threshold", type=float, default=3.5, 
                        help="Threshold for outlier removal (default: 3.5 for MAD, 1.5 for IQR). Higher values are more lenient.")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(message)s')
    logger = logging.getLogger()

    try:
        dfs = []
        for df_path in sorted(args.input_tables):
            sample_name = os.path.basename(os.path.dirname(df_path))
            df = read_table(df_path, sample_name)
            dfs.append(df)
            # logger.info(f"Read {df_path}: {len(df)} contigs")
    except Exception as e:
        logger.error(f"Error reading input files: {e}")
        return

    join_type = 'inner' if args.inner else 'outer'
    logger.info(f"Using {join_type} join")

    # Calculate totalAvgDepth
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
        # logger.info(f"After merge: {len(merged_df)} contigs")

    if merged_df.empty:
        logger.error("Error: The resulting merged table has no contigs. Please check your input data and parameters.")
        return

    filtered_df, outliers = detect_outliers_mag_pattern(merged_df, threshold=args.threshold)
    if not args.keep_outliers:
        merged_df = filtered_df

    logger.info(f"weird {len(outliers)} contigs with outlier coverage values:")
    logger.info(outliers.index.tolist())


    logger.info(f"Final number of contigs: {len(merged_df)}")

    if args.output:
        try:
            merged_df.to_csv(args.output, sep="\t")
            logger.info(f"Merged table saved to {args.output}")
        except Exception as e:
            logger.error(f"Error saving merged table: {e}")
    else:
        logger.info(merged_df.to_csv(sep='\t'))

if __name__ == "__main__":
    main()

