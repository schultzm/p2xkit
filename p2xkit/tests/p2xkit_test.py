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
        self.primers = pkg_resources.resource_filename(__parent_dir__, __test_primers__)
        self.template = pkg_resources.resource_filename(__parent_dir__, __test_template__)
    def psearcher(self):
        assertEqual(Psearcher, 'x')