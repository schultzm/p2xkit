from Bio.Emboss import PrimerSearch as psearch
from Bio.Emboss.Applications import PrimerSearchCommandline
from io import StringIO
import sys

class Psearcher:
    def __init__(self, template, primers):
        print(psearch)
        print("primers table")
        psearchcl = PrimerSearchCommandline()
        psearchcl.seqall = f"{self.template} -snucleotide1"
        psearchcl.infile = self.primers
        psearchcl.mismatchpercent = 20
        psearchcl.outfile = "stdout"
        print(psearchcl, file = sys.stderr)
        stdout, stderr = psearchcl()
        primersearch_results = psearch.read(StringIO(stdout))
        print(primersearch_results)
