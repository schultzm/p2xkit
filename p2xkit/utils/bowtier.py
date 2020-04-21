
import os
import sys
import shlex
from subprocess import Popen, PIPE
import pysam
from collections import defaultdict
from pathlib import Path, PurePath
from Bio import SeqIO

class Bowtier:
    def __init__(self, templates):
        # self.probes = probes
        self.templates = templates #PurePath

    def indexit(self):
        # Template strand bowtie2 index files generation
        # print()
        template_bowtie2_idx_fnames = list(self.templates.parent.glob(f"{self.templates.stem}.*.bt2"))
        if len(template_bowtie2_idx_fnames) != 6:
            cmd = f"bowtie2-build -f {self.templates} {PurePath(self.templates.parent, self.templates.stem)}"
            os.system(cmd)
        self.bowtieindex = list(self.templates.parent.glob(f"{self.templates.stem}.*.bt2"))

    def bowtieit(self, amplimer_table, probes):
        with open(probes, 'r') as input_handle:
            probes = list(SeqIO.parse(input_handle, 'fasta'))
            probes_dict = defaultdict(list)
            for probe in probes:
                probes_dict[probe.id].append(probe)
            # {probe.id: probe for probe in list(SeqIO.parse(input_handle, 'fasta'))}
            print(probes_dict)
        
        # bowtie2_map_cmd = shlex.split(f"bowtie2 -x {PurePath(template_fasta.parent, template_fasta.stem)} -U {probe.seq} -c -a --end-to-end --very-sensitive")
        pass