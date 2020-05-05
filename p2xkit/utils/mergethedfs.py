#!/usr/bin/env python3

import pandas as pd
import glob
import sys
import numpy as np

filenames = glob.glob('*_qpcr.tsv')+glob.glob('gisaid/*_qpcr.tsv')
print(filenames, file=sys.stderr)
dfs = []
for fname in filenames:
    #open only the dfs that have data
    try:
        dfs.append(pd.read_csv(fname, sep="\t"))
    except pd.errors.EmptyDataError:
        pass #   print(fname)

df = pd.concat(dfs)
#get only the rows with a value in probe_match_mismatch
# df = dfs_concat.loc[dfs_concat['probe_match_mismatch'].notnull()]
# print(df.columns)
#add the taxon info
taxa = pd.read_csv(sys.argv[1], sep="\t", header=None)
# print(taxa)

df['taxon'] = df['template_name'].apply(lambda x: f"''" if len(taxa.loc[taxa[0]==x][2].values) == 0 else f"{taxa.loc[taxa[0]==x][2].values[0]}")
print(df.to_csv(sep="\t"))

