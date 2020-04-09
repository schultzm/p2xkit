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
from p2xkit.utils.psearcher import Psearcher#, Hit_parser


class BedTestCasePass(unittest.TestCase):
    maxDiff = None
    def setUp(self):
        self.template = pkg_resources.resource_filename(__parent_dir__, __test_template__)
        self.primers  = pkg_resources.resource_filename(__parent_dir__, __test_primers__)
        self.probes   = pkg_resources.resource_filename(__parent_dir__, __test_probes__)
        self.mismatch = 20

    def psearcher(self):
        reaction = Psearcher(self.template, self.primers, self.mismatch)
        self.pcr_results = reaction.psearchit()
        self.assertEqual(self.pcr_results.amplifiers['N_Sarbeco_DE'][0].hit_info,
        """MN908947.3  
	Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome
	CACATTGGCACCCGCAATC hits forward strand at 28706 with 0 mismatches
	GAGGAACGAGAAGAGGCTTG hits reverse strand at [1071] with 0 mismatches""")
        # for key, value in pcr_result.amplifiers.items():
        #     for val in value:
        #         pass
                # print(val.hit_info)

    def iupaccheck(self):
        #this is weird.
        # need to create a psearcher object,
        # but it's perhaps too much for a simple
        # string parser.  Maybe break out the iupac def
        expanded = Psearcher(self.template,
                          self.primers,
                          self.mismatch)
        collapsed = expanded.collapsed_iupac("CA[GAR]ATGTTAAA[GCS]ACACTATTAGCATA")
        self.assertEqual("CARATGTTAAASACACTATTAGCATA", collapsed)

    def amplifiers_parsed(self):
        reaction = Psearcher(self.template,
                            self.primers,
                            self.mismatch)
        pcr_results = reaction.psearchit()
        for template_name, value in amplifiers.items():
            print(type(value[0]))

        # print(amplifiers)