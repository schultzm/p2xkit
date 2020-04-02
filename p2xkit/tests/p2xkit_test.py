"""
Unit Tests.
"""

import unittest
# from pathlib import Path, PurePath
from .. import (__parent_dir__,
                __test_primers__,
                __test_probes__,
                __test_template__)
import pkg_resources
from Bio.Emboss import PrimerSearch as psearch
from Bio.Emboss.Applications import PrimerSearchCommandline
from io import StringIO
import sys

class BedTestCasePass(unittest.TestCase):
    maxDiff = None
    def setUp(self):
        self.primers = pkg_resources.resource_filename(__parent_dir__, __test_primers__)
        self.template = pkg_resources.resource_filename(__parent_dir__, __test_template__)
    def psearcher(self):
        print("primers table")
        psearchcl = PrimerSearchCommandline()
        psearchcl.seqall = f"{self.template} -snucleotide1"
        psearchcl.infile = self.primers
        psearchcl.mismatchpercent = 20
        psearchcl.outfile = "stdout"
        print(psearchcl, file = sys.stderr)
        stdout, stderr = psearchcl()
        primersearch_results = psearch.read(StringIO(stdout))
        print(primersearch_results)
