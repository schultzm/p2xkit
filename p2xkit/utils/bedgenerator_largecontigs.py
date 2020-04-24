#!/usr/bin/env python3

strt=0
pcrupper=35000
maxsearchwindow=5000000
stop = maxsearchwindow
maxcontig=40000000
accession=None

import sys
from Bio import SeqIO
with open(sys.argv[1], 'r') as input_handle:
    records = [seq for seq in list(SeqIO.parse(input_handle, 'fasta'))]
    maxcontig = max([len(seq.seq) for seq in records])
    accessions = '\t'.join([seq.id for seq in records])



tab="\t"
while stop < maxcontig:
    print(f"{accessions}{tab}{strt}{tab}{stop}")
    strt = stop - pcrupper
    stop =strt + maxsearchwindow
else:
    print(f"{accessions}{tab}{strt}{tab}{maxcontig}")