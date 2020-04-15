import os
import sys
import shlex
from subprocess import Popen, PIPE
import pysam
from collections import defaultdict

class Bowtier:
    def __init__(self, probe, template):
        self.probes = probe
        self.template = template
        pass

    def indexit(self):
        # Template strand bowtie2 index files generation
        template_bowtie2_idx_fnames = list(Path(PurePath(template_fasta).parent).glob(f"{template_fasta.stem}.*.bt2"))
        if len(template_bowtie2_idx_fnames) != 6:
            os.system(f"bowtie2-build -f {template_fasta} {PurePath(template_fasta.parent, template_fasta.stem)}")

    def bowtieit(self):
        pass