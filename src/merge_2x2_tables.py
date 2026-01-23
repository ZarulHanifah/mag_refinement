#!/usr/bin/env python

import sys
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Merge multiple 2x2 tables with the same MAG index but different sample names.")
    parser.add_argument("-t", "--tables", required=True, nargs='+', help="Paths to the table TSV files.")
    parser.add_argument("-o", "--output", help="Path to the output merged CSV file (prints to stdout if not specified).")
    parser.add_argument("--outer", action='store_true', help="Use outer join instead of inner join.")
    args = parser.parse_args()

    # Initialize the merge operation
    merged_table = None

    for table_path in args.tables:
        # Load each 2x2 table
        table = pd.read_csv(table_path, sep='\t', index_col=0)
        
        # Merge the current table into the accumulated DataFrame
        if merged_table is None:
            merged_table = table
        else:
            join_type = 'outer' if args.outer else 'inner'
            merged_table = merged_table.join(table, how=join_type)
    
    # Output the merged table
    if args.output:
        merged_table.to_csv(args.output, sep='\t')
    else:
        merged_table.to_csv(sys.stdout, sep='\t')

if __name__ == "__main__":
    main()
