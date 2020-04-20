
import os
import sys
import shlex
from subprocess import Popen, PIPE
import pysam
from collections import defaultdict
from pathlib import Path, PurePath

class Bowtier:
    def __init__(self, probetable, templates):
        self.probes = probetable
        self.templates = templates #PurePath

    def indexit(self):
        # Template strand bowtie2 index files generation
        # print()
        template_bowtie2_idx_fnames = list(self.templates.parent.glob(f"{self.templates.stem}.*.bt2"))
        if len(template_bowtie2_idx_fnames) != 6:
            cmd = f"bowtie2-build -f {self.templates} {PurePath(self.templates.parent, self.templates.stem)}"
            os.system(cmd)

    def bowtieit(self):
        pass