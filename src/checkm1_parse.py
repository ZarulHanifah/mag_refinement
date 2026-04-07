#!/usr/bin/env python

import re
import json
from pprint import pprint

path_gene_stats = './results/checkm1/original/C1E5_M_metabat.1297/storage/marker_gene_stats.tsv'

with open(path_gene_stats) as h:
    data = h.read()

data = data.split("\t")[1]
data = re.sub("\'", "\"", data)
data = json.loads(data)


pprint(data)
