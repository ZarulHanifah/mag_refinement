#!/usr/bin/env python

import sys
import pandas as pd

def main():
    # Read the TSV input from stdin
    df = pd.read_csv(sys.stdin, sep="\t", index_col=0)
    
    # Transpose the DataFrame
    df_transposed = df.T
    
    # Write the transposed DataFrame to stdout as TSV
    df_transposed.to_csv(sys.stdout, sep="\t")

if __name__ == "__main__":
    main()

