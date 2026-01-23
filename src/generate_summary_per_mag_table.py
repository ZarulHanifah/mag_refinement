import pandas as pd
import argparse
import sys

def generate_table(input_file, output_file, sample, mag):
    # Read the input TSV file
    df = pd.read_csv(input_file, sep='\t')
    
    # Ensure the column for SAMPLE.bam is correctly identified
    sample_col = df.columns[3]  # Fourth column
    
    # Calculate the required values for the bottom-right cell
    df['k'] = df['contigLen'] * df[sample_col]
    sum_k = df['k'].sum()
    sum_contig_len = df['contigLen'].sum()
    bottom_right = sum_k / sum_contig_len if sum_contig_len > 0 else 0

    # Create the output 2x2 table
    output_table = pd.DataFrame(
        data=[[bottom_right]],
        index=[mag],
        columns=[sample]
    )
    
    # Write the output table to a file
    output_table.to_csv(output_file, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a 2x2 summary table from a TSV file.")
    parser.add_argument('-i', '--input', help="Input TSV file.", required=True)
    parser.add_argument('-o', '--output', help="Output TSV file.", required=True)
    parser.add_argument('-s', '--sample', help="Sample name for the top right cell.", required=True)
    parser.add_argument('-m', '--mag', help="MAG name for the bottom left cell.", required=True)
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    args = parser.parse_args()
    
    generate_table(args.input, args.output, args.sample, args.mag)

