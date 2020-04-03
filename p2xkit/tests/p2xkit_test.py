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
from p2xkit.utils.psearcher import Psearcher, Hit


class BedTestCasePass(unittest.TestCase):
    maxDiff = None
    def setUp(self):
        self.template = pkg_resources.resource_filename(__parent_dir__, __test_template__)
        self.primers = pkg_resources.resource_filename(__parent_dir__, __test_primers__)
        self.mismatch = 20

    def psearcher(self):
        result = Psearcher(self.template, self.primers, self.mismatch)
        result.psearchit()
        pcr_results = result.pcr_hits
        self.assertEqual(pcr_results.amplifiers['N_Sarbeco_DE'][0].hit_info,
        """MN908947.3  
	Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome
	CACATTGGCACCCGCAATC hits forward strand at 28706 with 0 mismatches
	GAGGAACGAGAAGAGGCTTG hits reverse strand at [1071] with 0 mismatches""")
        for key, value in pcr_results.amplifiers.items():
            print(key, value)

    def iupaccheck(self):
        inseq = Hit("CA[GAR]ATGTTAAA[GCS]ACACTATTAGCATA")
        inseq.collapsed_iupac()
        self.assertEqual("CARATGTTAAASACACTATTAGCATA", inseq.collapsed)
