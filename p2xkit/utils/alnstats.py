#!/usr/bin/env python3

databin = """RdRP_SARSr_DE	=========================	MT123293	human
RdRP_SARSr_DE	===========X========X====	NC_045512	SARS2
RdRP_SARSr_DE	===========X========X====	MT123291	SARS2
RdRP_SARSr_DE	===========X=============	MT123292	SARS1"""

import pandas as pd
from io import StringIO
df = pd.read_csv(StringIO(databin), sep="\t", header=None)
df[1] = df[[1]].apply(lambda x: f"{x.values[0].replace('=', '0').replace('X', '1')}", axis=1)
print(df)
print(df[1].str.split('', expand=True).rename(columns = lambda x: "string"+str(x+1)))
df = pd.concat([df, df[1].str.split('', expand=True).rename(columns = lambda x: "string"+str(x+1))], axis=1)
print(df.to_csv(sep="\t"))

