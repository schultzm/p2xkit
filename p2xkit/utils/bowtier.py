
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
        # this could all be done by creating a mini file of the amplimer and mapping to that. todo. 
        # instead of mapping to whole template

    def indexit(self):
        # Template strand bowtie2 index files generation
        template_bowtie2_idx_fnames = list(self.templates.parent.glob(f"{self.templates.stem}.*.bt2"))
        if len(template_bowtie2_idx_fnames) != 6:
            cmd = f"bowtie2-build -f {self.templates} {PurePath(self.templates.parent, self.templates.stem)}"
            os.system(cmd)
        self.bowtieindex = list(self.templates.parent.glob(f"{self.templates.stem}.*.bt2"))

    def bowtieit(self, amplimer_table, probes):
        print(amplimer_table.to_csv(sep="\t"))
        probes_dict = defaultdict(list)
        with open(probes, 'r') as input_handle:
            probes = list(SeqIO.parse(input_handle, 'fasta'))
            for probe in probes:
                probes_dict[probe.id].append(probe)
        for primerpair_name, probes_list in probes_dict.items():
            run_cmd = [f"bowtie2 -x {PurePath(self.templates.parent, self.templates.stem)} -U {tqmanprobe.seq} -c -a --end-to-end --very-sensitive" for tqmanprobe in probes_list]
            print(run_cmd)
        print()
        # print(run_cmd)
