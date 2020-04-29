#!/usr/bin/env python3

import pandas as pd
import sys
import glob
import os

fnames = [fname for fname in open('probehits.txt', 'r').readlines()]
print(fnames)

# df = pd.read_csv(sys.argv[1], sep="\t", header=0)
# df.columns = headers
# df.probe_match_mismatch.astype(str)
# for i in df.probe_match_mismatch:
#     if isinstance(i, float):
#         print(i)
# # sys.exit()
# df['probe_match_mismatch'] = df[['probe_match_mismatch']].apply(lambda x: f"{x.values[0].replace('=', '0').replace('X', '1')}", axis=1)
# # print(df.to_csv(sep="\t"))
# # sys.exit()
# headerpre = 'pos'
# df = pd.concat([df, df['probe_match_mismatch'].str.split('', expand=True).rename(columns = lambda x: f"{headerpre}"+str(x+1))], axis=1)
# #drop empty columns
# for col in df.columns:
#     if not list(filter(None, set(df[col]))):
#         df.drop(col, inplace=True, axis=1)

# # sys.exit()
# dfmelted = df.melt(id_vars=['taxon'], value_vars=df.columns.values[4:])
# print(dfmelted.to_csv(sep="\t"))