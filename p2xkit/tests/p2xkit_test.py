"""
Unit Tests.
"""

import unittest
from pathlib import Path
# from .. import __parent_dir__, __test_tree__
import pkg_resources

class TreeTestCasePass(unittest.TestCase):
    maxDiff = None
    def setUp(self):
        print("")
        self.primers = ()
