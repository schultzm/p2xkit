import unittest
from ..tests.budgitree_test import TreeTestCasePass


def suite():
    """
    This is the test suite.
    """
    suite = unittest.TestSuite()
    # suite.addTest(TreeTestCasePass('ete_print_tree'))
    return suite