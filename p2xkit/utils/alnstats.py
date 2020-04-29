#!/usr/bin/env python3

databin = """primer	matchmismatch	accession	taxon
RdRP_SARSr_DE	=========================	MT123293	human
RdRP_SARSr_DE	===========X========X====	NC_045512	SARS2
RdRP_SARSr_DE	===========X========X====	MT123291	SARS2
RdRP_SARSr_DE	===========X=============	MT123292	SARS1"""

import pandas as pd
from io import StringIO
df = pd.read_csv(StringIO(databin), sep="\t", header=0)
df['matchmismatch'] = df[['matchmismatch']].apply(lambda x: f"{x.values[0].replace('=', '0').replace('X', '1')}", axis=1)
headerpre = 'pos'
df = pd.concat([df, df['matchmismatch'].str.split('', expand=True).rename(columns = lambda x: f"{headerpre}"+str(x+1))], axis=1)
for col in df.columns:
    if not list(filter(None, set(df[col]))):
        df.drop(col, inplace=True, axis=1)
dfmelted = df.melt(id_vars=['taxon'], value_vars=df.columns.values[4:])
print(dfmelted.to_csv(sep="\t"))