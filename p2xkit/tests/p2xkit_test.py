











"""
Unit Tests.
"""

import unittest
from pathlib import Path, PurePath, PosixPath
from .. import (__parent_dir__,
                __test_primers__,
                __test_probes__,
                __test_templates__)
import pkg_resources
from ..utils.psearcher import Psearcher, _collapsed_iupac, _iupac_zipper #, Hit_parser
from ..utils.bowtier import Bowtier


class BedTestCasePass(unittest.TestCase):
    maxDiff = None
    def setUp(self):
        self.templates = Path(PurePath(pkg_resources.resource_filename(__parent_dir__,
                                                        __test_templates__)))
        self.primers  = Path(PurePath(pkg_resources.resource_filename(__parent_dir__,
                                                        __test_primers__)))
        self.probes   = Path(PurePath(pkg_resources.resource_filename(__parent_dir__,
                                                        __test_probes__)))
        self.mismatch = 20

    def psearcher(self):
        reaction = Psearcher(self.templates, self.primers, self.mismatch,
                             False, 0, 50000, 60000)
        reaction.psearchit()
        self.assertEqual(reaction.pcr_results.amplifiers['N_Sarbeco_DE'][0]. \
                         hit_info,
        """MN908947.3  
	Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome
	CACATTGGCACCCGCAATC hits forward strand at 28706 with 0 mismatches
	GAGGAACGAGAAGAGGCTTG hits reverse strand at [1071] with 0 mismatches""")

    def iupaccheck(self):
        self.assertEqual("CARATGTTAAASACACTATTAGCATA",
                         _collapsed_iupac("CA[GAR]ATGTTAAA[GCS]ACACTATTAGCATA"))

    def iupaczipper(self):
        self.assertEqual(_iupac_zipper('CCAGGTGGWACRTCATCMGGTGATGC', 'CCAGGTGGAACCTCATCAGGAGATGC'),
                         "'===========X========X=====")

    def amplifiers_parsed(self):
        reaction = Psearcher(self.templates,
                            self.primers,
                            self.mismatch,
                            False,
                            0,
                            50000,
                            60000)
        reaction.psearchit() # get PrimerSearch.OutputRecords
        amplimer_table = reaction.amplimer_table() # get the full table
        self.assertEqual(amplimer_table.iloc[4].loc['rev_match_mismatch'], "'========================")

    def bowtie_map(self):
        reaction = Psearcher(self.templates,
                             self.primers,
                             self.mismatch,
                             False,
                             0,
                             50000,
                             60000)
        reaction.psearchit() # get PrimerSearch.OutputRecords
        amplimer_table = reaction.amplimer_table()
        qpcr_setup = Bowtier(amplimer_table, self.templates, self.probes, False)
        mapped = qpcr_setup.bowtieit()
        self.assertEqual(mapped.iloc[4,3], 'N_CN P')
        self.assertEqual(mapped.iloc[4,6], 'N_CN')