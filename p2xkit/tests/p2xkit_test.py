











"""
Unit Tests.
"""

import unittest
from pathlib import Path, PurePath, PosixPath
from .. import (__parent_dir__,
                __test_primers__,
                __test_probes__,
                __test_template__, #maybe redundant with templates in place
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
        # self.templates = pkg_resources.resource_filename(__parent_dir__,
        #                                                  __test_templates__)

    def psearcher(self):
        reaction = Psearcher(self.templates, self.primers, self.mismatch)
        reaction.psearchit()
        self.assertEqual(reaction.pcr_results.amplifiers['N_Sarbeco_DE'][0]. \
                         hit_info,
        """MN908947.3  
	Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome
	CACATTGGCACCCGCAATC hits forward strand at 28706 with 0 mismatches
	GAGGAACGAGAAGAGGCTTG hits reverse strand at [1071] with 0 mismatches""")
        # for key, value in pcr_result.amplifiers.items():
        #     for val in value:
        #         pass
                # print(val.hit_info)

    def iupaccheck(self):
        self.assertEqual("CARATGTTAAASACACTATTAGCATA",
                         _collapsed_iupac("CA[GAR]ATGTTAAA[GCS]ACACTATTAGCATA"))

    def amplifiers_parsed(self):
        reaction = Psearcher(self.templates,
                            self.primers,
                            self.mismatch)
        reaction.psearchit() # get PrimerSearch.OutputRecords
        amplimer_table = reaction.amplimer_table() # get the full table
        # print(reaction.amplimer_tab.to_csv(sep="\t"))
        self.assertEqual(amplimer_table.iloc[4].loc['rev_match_mismatch'], "'===================")

    def bowtie_index(self):
        reaction = Psearcher(self.templates,
                    self.primers,
                    self.mismatch)
        reaction.psearchit() # get PrimerSearch.OutputRecords
        amplimer_table = reaction.amplimer_table()
        indexed = Bowtier(self.templates)
        indexed.indexit()
        # indexes = list(self.templates.parent.glob(f"{self.templates.stem}.*.bt2"))
        # print(indexes)
        self.assertEqual(indexed.bowtieindex[0].suffix, '.bt2')
        # for i in indexes:
        #     i.unlink() #remove all the index files
    def bowtie_map(self):
        reaction = Psearcher(self.templates,
                    self.primers,
                    self.mismatch)
        reaction.psearchit() # get PrimerSearch.OutputRecords
        amplimer_table = reaction.amplimer_table()

        # print(amplimer_table.to_csv(sep="\t"))
        indexed = Bowtier(self.templates)
        indexed.indexit()
        mapped = indexed.bowtieit(amplimer_table, self.probes)
        # mapped.bowtieit(amplimer_table, self.probes)
        
        self.assertEqual(mapped, 'x')
