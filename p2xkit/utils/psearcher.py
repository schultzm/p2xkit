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

class Hit:
    def __init__(self, seqstring):
        self.expanded = seqstring
    
    def collapsed_iupac(self):
        '''
        Given an 'expanded' seqstring, 'CA[GAR]ATGTTAAA[GCS]ACACTATTAGCATA',
        return collapsed seqstring 'CARATGTTAAASACACTATTAGCATA' (with IUPAC codes).
        '''
        IUPAC_codes = '''RYSWKMBDHVNryswkmbdhvn'''
        collapsed = ''
        regions = self.expanded.replace('[', ']').split(']')
        for region in regions:
            IUPAC_letter = None
            for nuc_letter in region:
                if nuc_letter in IUPAC_codes:
                    IUPAC_letter = region[-1]
            if IUPAC_letter:
                collapsed += IUPAC_letter
            else:
                collapsed += region
        self.collapsed = collapsed

