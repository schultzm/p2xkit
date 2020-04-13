from Bio.Emboss import PrimerSearch as psearch
from Bio.Emboss.Applications import PrimerSearchCommandline
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Align import MultipleSeqAlignment as MSA


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
        return a dictionary of {key(primer pair name [str]):
                                value ([list] of amplifiers)}
        """
        psearchcl = PrimerSearchCommandline()
        psearchcl.seqall = f"{self.template} -snucleotide1"
        psearchcl.infile = self.primers
        psearchcl.mismatchpercent = self.mismatch
        psearchcl.outfile = "stdout"
        stdout, stderr = psearchcl()
        # print(stdout)
        self.pcr_results = psearch.read(StringIO(stdout))

    def _collapsed_iupac(self, primerstring):
        '''
        Given an 'expanded' seqstring, 'CA[GAR]ATGTTAAA[GCS]ACACTATTAGCATA',
        return collapsed seqstring 'CARATGTTAAASACACTATTAGCATA' (with IUPAC codes).
        '''
        IUPAC_codes = '''RYSWKMBDHVNryswkmbdhvn'''
        collapsed = ''
        regions = primerstring.replace('[', ']').split(']')
        for region in regions:
            IUPAC_letter = None
            for nuc_letter in region:
                if nuc_letter in IUPAC_codes:
                    IUPAC_letter = region[-1]
            if IUPAC_letter:
                collapsed += IUPAC_letter
            else:
                collapsed += region
        return collapsed

    def amplimer_table(self):
        '''
        An amplimer is a PCR product between a primer-pair.
        The amplifier is a dict, storing these data:
            keys are primer-pair names
            values are lists of amplifiers
                each amplifier in amplifiers has attributes:
                    hit_info (the amplimer)
                    length (of amplimer)
        '''
        results_dfs_list = []
        for primerpair_name, amplifier in self.pcr_results.amplifiers.items():
            for index, amplimer in enumerate(amplifier): # need the indices
                print(primerpair_name, f"Amplimer {index}:", amplimer.hit_info)

        #amplicons = {primerpair: template_list for primerpair, template_list in amplifiers.items()}
        # return(amp_dict)
        # for key, template_list in amp_dict.items:
            
        # for primerpair, template_list in amplifiers.items():
        #     print(primerpair, template_list)
        #     # For each 
        #     # enumerate in this loop to track (index) the amplimer (hit or "amplimer"")
            # for index, amplimer in enumerate(template_list):
        #         # Reformat stdout from primersearch to a list for parsing
        # self.results_dfs_list = pd.concat(results_dfs_list)