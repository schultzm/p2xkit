
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
        # print(amplimer_table.to_csv(sep="\t"))
        map_results_dfs = [] 
        probes_dict = defaultdict(list)
        with open(probes, 'r') as input_handle:
            probes = list(SeqIO.parse(input_handle, 'fasta'))
            for probe in probes:
                print(probe)
                probes_dict[probe.id].append(probe)
        for primerpair_name, probes_list in probes_dict.items():
            # -L <int>           length of seed substrings; must be >3, <32 (22)
            # -c                 <m1>, <m2>, <r> are sequences themselves, not files
            # -U unpairedreads
            # -a/--all           report all alignments; very slow, MAPQ not meaningful 
            # -N <int>           max # mismatches in seed alignment; can be 0 or 1 (0)
            # --np <int>         penalty for non-A/C/G/Ts in read/ref (1)
            # -R <int>           for reads w/ repetitive seeds, try <int> sets of seeds (2)
            map_cmds = [(f"bowtie2 -x {PurePath(self.templates.parent, self.templates.stem)} -U {tqmanprobe.seq} -c --all --end-to-end --very-sensitive -L 3 -N 1 --np 0 -R 10", tqmanprobe.description) for tqmanprobe in probes_list]
            for map_cmd in map_cmds:
                # map_cmd[0] is the bowtie2 cmd
                # map_cmd[1] is the probe name
                proc1 = Popen(shlex.split(map_cmd[0]), stdout=PIPE, stderr=PIPE)
                samfile = proc1.stdout.fileno()
                with pysam.AlignmentFile(samfile, "r") as sam:
                    for rec in sam.fetch():
                        if not rec.is_unmapped:

                            # print(rec.to_dict())
                            # if rec.flag == 16:
                            #         seq_str = probe.seq.reverse_complement()
        #                         id = f'probe_{idx}_rcomp_{str(primerpair)}'
        #                     else:
        #                         seq_str = probe.seq
        #                         id = f'probe_{idx}_{str(primerpair)}'

                    print()
 
