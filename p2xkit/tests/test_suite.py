import unittest
from ..tests.p2xkit_test import BedTestCasePass


def suite():
    """
    This is the test suite.
    """
    suite = unittest.TestSuite()
    suite.addTest(BedTestCasePass("psearcher"))
    suite.addTest(BedTestCasePass("iupaccheck"))
    suite.addTest(BedTestCasePass("amplifiers_parsed"))
<<<<<<< HEAD
    suite.addTest(BedTestCasePass("bowtie_index"))
    suite.addTest(BedTestCasePass("bowtie_map"))
=======
    suite.addTest(BedTestCasePass("blast_orientate"))
    suite.addTest(BedTestCasePass("mafft_aln"))
>>>>>>> 75e15d93f1f505cfbd5d3ae082b668142e304b76
    return suite