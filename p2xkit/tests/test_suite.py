import unittest
from ..tests.p2xkit_test import BedTestCasePass


def suite():
    """
    This is the test suite.
    """
    suite = unittest.TestSuite()
    suite.addTest(BedTestCasePass("primers_table"))
    return suite