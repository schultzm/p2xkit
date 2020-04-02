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
from p2xkit.utils.psearcher import Psearcher


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
        # print(result.pcr_hits)

        # # print(result())
        for key, value in pcr_results.amplifiers.items():
            for i in value:
                print(key, "\n", i.hit_info, i.length)
        # self.assertEqual(Psearcher, 'x')