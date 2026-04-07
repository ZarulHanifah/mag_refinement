#!/usr/bin/env python

import ast
from pprint import pprint
import pandas as pd

name = "C1E5_M_metabat.1297"
path_gene_stats = f"./results/checkm1/original/{name}/storage/bin_stats_ext.tsv"

# with open(path_gene_stats) as h:
#     data = h.read()
#
# data = data.split("\t")[1]
# data = re.sub("\'", "\"", data)
# data = json.loads(data)
#
#
# pprint(data)

df = pd.read_csv(path_gene_stats, header=None, sep="\t", index_col=0)
val = df.loc[name, 1]
val = ast.literal_eval(val)

rename_val_keys = {
    "Completeness": "Completeness",
    "Contamination": "Contamination",
    "Translation table": "Translation_Table_Used",
    "Coding density": "Coding_Density",
    "N50 (contigs)": "Contig_N50",
    "Genome size": "Genome_Size",
    "GC": "GC_Content",
    "# predicted genes": "Total_Coding_Sequences",
    "# contigs": "Total_Contigs",
    "Longest contig": "Max_Contig_Length",
}

val = { rename_val_keys[k]: val[k] for k in rename_val_keys.keys() }
val = pd.DataFrame(val, index=[name])

print(val)
