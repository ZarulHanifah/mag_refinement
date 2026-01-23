import argparse
import pandas as pd
import numpy as np
import re

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process paths for cluster and checkm files.")
    parser.add_argument("--cluster", required=True, help="Path to the cluster file")
    parser.add_argument("--checkm", required=True, help="Path to the checkm file")
    parser.add_argument("-o", "--outdf", required=True, help="Path to the output summary df")
    parser.add_argument("-l", "--logfile", required=True, help="Path to the log of selecting representative of each cluster")
    return parser.parse_args()

def main():
    args = parse_arguments()
    cluster_path = args.cluster
    checkm_path = args.checkm
    outdf = args.outdf
    logfile = args.logfile

    # Now you can use cluster_path and checkm_path in your code
    cluster_df = pd.read_csv(cluster_path, index_col=0)
    cluster_df.index = [re.sub(".fasta", "", i) for i in cluster_df.index.tolist()]

    checkm_df = pd.read_csv(checkm_path, sep="\t", index_col=0)

    def calculate_score(sub_checkm_df):
        return sub_checkm_df["Completeness"] - (5 * sub_checkm_df["Contamination"]) + (0.5 * np.log(sub_checkm_df["Contig_N50"]))

    final_df = []
    final_df_index = []
    final_df_columns = ["Completeness", "Contamination", "Contig_N50", "genome_size", "GC_Content"]

    final_df = pd.DataFrame(columns=final_df_columns)

    with open(logfile, "w") as f:
        for i in cluster_df["secondary_cluster"].unique():
            sub_cluster_df = cluster_df.loc[cluster_df["secondary_cluster"] == i, :]
            
            sub_checkm_df = checkm_df.copy().loc[sub_cluster_df.index.tolist(), :]
            sub_checkm_df["score"] = calculate_score(sub_checkm_df)
            sub_checkm_df.sort_values("score", ascending=False, inplace=True)
            
            f.write(f"Cluster {i}\n")
            f.write(sub_checkm_df.to_string())
            f.write("\n\n")
            first_choice = sub_checkm_df.index.tolist()[0]
            final_df_index.append(first_choice)
            c = sub_checkm_df.loc[[first_choice], final_df_columns]

            final_df = pd.concat([final_df, c])
        
    # final_df = pd.DataFrame(final_df, index=final_df_index, columns=final_df_columns)
    final_df.to_csv(outdf, sep="\t")

if __name__ == "__main__":
    main()

