import os
import sys
import shlex
from subprocess import Popen, PIPE
# import pysam
from collections import defaultdict
from pathlib import Path, PurePath

class Mafft_aln:
    def __init__(self, template):
        # self.probes = probe
        self.template = template
        print(self.template)
        pass

    def blastit(self):
        pass
        # Template strand bowtie2 index files generation
        # template_bowtie2_idx_fnames = list(Path(PurePath(self.template).parent).glob(f"{self.template.stem}.*.bt2"))
        # print(template_bowtie2_idx_fnames)
        # if len(template_bowtie2_idx_fnames) != 6:
        #     os.system(f"bowtie2-build -f {template_fasta} {PurePath(template_fasta.parent, template_fasta.stem)}")

    def alignit(self):
        pass