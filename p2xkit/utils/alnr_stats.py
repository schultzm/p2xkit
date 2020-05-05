#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import glob
import os

# translation=pd.read_csv('trans.tsv', sep="\t")
# # print(translation)
# fnames = [fname.rstrip() for fname in open('probehits.txt', 'r').readlines()]
# # print(fnames)
# dflist = [pd.read_csv(fname, sep="\t", header=0) for fname in fnames]
# df = pd.concat(dflist, axis=0)
# print(df)
#these are the offtarget hits: ['AY278487', 'AY278488', 'AY278489', 'AY278490', 'AY278741', 'AY279354', 'AY304486', 'AY304488', 'AY338174', 'AY338175', 'AY348314', 'AY394850', 'AY572034', 'AY572035', 'AY572038', 'AY686863', 'AY686864', 'AY772062', 'AY864806', 'CS460762', 'DQ071615', 'DQ182595', 'DQ412042', 'DQ412043', 'DQ898174', 'FJ211859', 'FJ588686', 'FJ882926', 'FJ882927', 'FJ882928', 'FJ882929', 'FJ882930', 'FJ882931', 'FJ882932', 'FJ882933', 'FJ882934', 'FJ882935', 'FJ882936', 'FJ882937', 'FJ882938', 'FJ882939', 'FJ882940', 'FJ882941', 'FJ882944', 'FJ882946', 'FJ882947', 'FJ882949', 'FJ882950', 'FJ882954', 'FJ882955', 'FJ882956', 'FJ882960', 'LC522972', 'LC522973', 'LC522974', 'LC522975', 'LC528232', 'LC528233', 'LR757995', 'LR757996', 'LR757997', 'LR757998', 'MN908947', 'MN938384', 'MN975262', 'MN985325', 'MN988668', 'MN988669', 'MN988713', 'MN994467', 'MN994468', 'MN996527', 'MN996528', 'MN996529', 'MN996530', 'MN996531', 'MN997409', 'MT007544', 'MT019529', 'MT019530', 'MT019531', 'MT019532', 'MT019533', 'MT020781', 'MT020880', 'MT020881', 'MT027062', 'MT027063', 'MT027064', 'MT039873', 'MT039887', 'MT039888', 'MT039890', 'MT044257', 'MT044258', 'MT049951', 'MT066175', 'MT066176', 'MT072688', 'MT093571', 'MT093631', 'MT106052', 'MT106053', 'MT106054', 'MT118835', 'MT123290', 'MT123291', 'MT123292', 'MT123293', 'NC_045512']
df = pd.read_csv(sys.argv[1], sep="\t", header=0)
# df.columns = headers
# df.probe_match_mismatch.astype(str)
# for i in df.probe_match_mismatch:
#     if isinstance(i, float):
#         print(i)
# # sys.exit()
quoter ="'"
# df['probe_match_mismatch'] = df[['probe_match_mismatch']].apply(lambda x: "NA" if isinstance(x.values[0], float) else f"{x.values[0].replace(quoter, '').replace('=', '0').replace('X', '1').strip()}", axis=1)
# print(df['probe_match_mismatch'])
# sys.exit()
# # print(df.to_csv(sep="\t"))
# # sys.exit()
headerpre = 'pos'
# print(df['probe_match_mismatch'])
# sys.exit()
df = df.replace(np.nan, '', regex=True).replace("'", '', regex=True)
# df['probe_match_mismatch'] = df['probe_match_mismatch'].apply(lambda x: x.replace("'", ""))#.replace("'", "", inplace=True)
df = pd.concat([df, df['probe_match_mismatch'].str.split('', expand=True).rename(columns = lambda x: f"{headerpre}"+str(x))], axis=1)
# for i in pd.unique(df['template_name']): #df.template_name.unique)
#     print(i)
# # print(df.to_csv(sep="\t"))
#drop empty columns
for col in df.columns:
    if not list(filter(None, set(df[col]))):
        df.drop(col, inplace=True, axis=1)
    # elif list(set(df[col])) == {"'"}:
    #     df.drop(col, inplace=True, axis=1)

# # sys.exit()
# print(df.columns.values[4:])
# print([i if 'pos' in i for i in df.columns.values])
print(df.to_csv(sep="\t"))
# dfmelted = pd.melt(df, id_vars=[('template_name', 'taxon')], value_vars=[('pos1', 'pos2')])
# print(dfmelted.to_csv(sep="\t"))