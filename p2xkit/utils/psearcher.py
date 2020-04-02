from Bio.Emboss import PrimerSearch as psearch
from Bio.Emboss.Applications import PrimerSearchCommandline
from io import StringIO
import sys
import Bio

class Psearcher:
    def __init__(self, template, primers, mismatch):
        self.template = template
        self.primers = primers
        self.mismatch = mismatch

    def psearchit(self):
        """
        run primersearch
        return a dictionary of amplifier objects
            keys are primer-pair names
            values are amplifier
                each amplifier is a list of PCR hits with attributes:
                    hit_info
                    length
        """
        psearchcl = PrimerSearchCommandline()
        psearchcl.seqall = f"{self.template} -snucleotide1"
        psearchcl.infile = self.primers
        psearchcl.mismatchpercent = self.mismatch
        psearchcl.outfile = "stdout"
        # print(psearchcl, file = sys.stderr)
        stdout, stderr = psearchcl()
        self.pcr_hits = psearch.read(StringIO(stdout))
        #  = primersearch_results
