from Bio.Emboss import PrimerSearch as psearch
from Bio.Emboss.Applications import PrimerSearchCommandline
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Align import MultipleSeqAlignment as MSA
import pandas as pd

from io import StringIO
import sys
import Bio

class Psearcher:
    def __init__(self, template, primers, mismatch):
        self.template = template
        self.template_seqs = {seq.id: seq for seq in list(SeqIO.parse(open(self.template, 'r'), 'fasta'))}
        self.primers = primers
        self.mismatch = mismatch
        # print(self.template_seqs)

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
        Fix this...
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
            if not amplifier: # i.e., if list is empty
                df = pd.DataFrame({}, index=[primerpair_name])
                results_dfs_list.append(df)
            else:
                for index, amplimer in enumerate(amplifier): # need the indices
                    sub_df = {}
                    hit = amplimer.hit_info.replace('\t', '').rstrip(' ').split('\n')
                    hit = [i.strip() for i in hit]
                    sub_df['amplimer_n'] = index+1
                    sub_df['template_name'] = hit[0]
                    # print('HIT0', hit[0])
                    fwd = hit[-2].split(' ')
                    rev = hit[-1].split(' ')
                    sub_df['fwd_oligo'] = self._collapsed_iupac(fwd[0]).upper()#, alphabet = IUPAC.ambiguous_dna)
                    sub_df['rev_oligo'] = self._collapsed_iupac(rev[0]).upper()#, alphabet = IUPAC.ambiguous_dna)
                    sub_df['fwd_oligo_tmplt_start'] = int(fwd[-4])
                    sub_df['fwd_oligo_tmplt_end']   = sub_df['fwd_oligo_tmplt_start'] + len(sub_df['fwd_oligo'])
                    sub_df['rev_oligo_tmplt_end']   = len(self.template_seqs[hit[0]].seq) - int(rev[-4].replace('[', '').replace(']', ''))
                    sub_df['rev_oligo_tmplt_start'] = sub_df['rev_oligo_tmplt_end'] - len(sub_df['rev_oligo'])
                    sub_df['amplicon'] = str(self.template_seqs[sub_df['template_name']].seq[sub_df['fwd_oligo_tmplt_end']:sub_df['rev_oligo_tmplt_start']].upper())
                    sub_df['fwd_oligo_match'] = str(self.template_seqs[sub_df['template_name']].seq[sub_df['fwd_oligo_tmplt_start']:sub_df['fwd_oligo_tmplt_end']])
                    sub_df['rev_oligo_match'] = str(self.template_seqs[sub_df['template_name']].seq[sub_df['rev_oligo_tmplt_start']:sub_df['rev_oligo_tmplt_end']])
                    sub_df['product']  = str(sub_df['fwd_oligo']+sub_df['amplicon']+Seq(sub_df['rev_oligo'], alphabet = IUPAC.ambiguous_dna).reverse_complement())
                    # print(sub_df['amplicon'])
                    df = pd.DataFrame(sub_df, index=[primerpair_name])
                    results_dfs_list.append(df)
        self.amplimer_tab = pd.concat(results_dfs_list)
