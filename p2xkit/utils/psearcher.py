from Bio.Emboss import PrimerSearch as psearch
from Bio.Emboss.Applications import PrimerSearchCommandline
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Align import MultipleSeqAlignment as MSA
import pandas as pd
from Bio.Data.IUPACData import ambiguous_dna_values
from io import StringIO
import sys
import Bio
# from Bio.Application import _Option #allows me to append options to the primersearchcl

# Moved next two functions to module scope
def _collapsed_iupac(primerstring):
    '''
    Given an 'expanded' seqstring, 'CA[GAR]ATGTTAAA[GCS]ACACTATTAGCATA',
    return collapsed seqstring 'CARATGTTAAASACACTATTAGCATA' (with IUPAC codes).
    '''
    IUPAC_codes = '''RYSWKMBDHVNryswkmbdhvn''' #allows upper and lowercase iupac
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

def _iupac_zipper(seq1, seq2):
    # print('seq1', seq1)
    # print('seq2', seq2)
    zipped = list(zip([i for i in seq1.upper()], [i for i in seq2.upper()]))
    binrep = ''
    for i in zipped:
        a = [j for j in ambiguous_dna_values[i[0]]]
        b = [j for j in ambiguous_dna_values[i[1]]]
        intersection = set(a).intersection(set(b))
        if intersection:
            binrep += f"="
        else:
            binrep += f"X"
    return f"'{binrep}"


class Psearcher:
    def __init__(self, template, primers, mismatch, reverse_complement, begin, end):
        self.template = template
        self.template_seqs = {seq.id: seq for seq in list(SeqIO.parse(open(self.template, 'r'), 'fasta'))}
        self.reverse_complement = reverse_complement
        if self.reverse_complement:
            self.template_seqs = {seq.id: seq.reverse_complement() for seq in list(SeqIO.parse(open(self.template, 'r'), 'fasta'))}
        self.primers = primers
        self.mismatch = mismatch
        self.begin = begin
        self.end = end

    def psearchit(self):
        """
        run primersearch
        return a dictionary of {key(primer pair name [str]):
                                value ([list] of amplifiers)}
        """
        # -sbegin1            integer    Start of each sequence to be used
        # -send1
        sreverse = '-sreverse1=N'
        if self.reverse_complement:
            sreverse = '-sreverse1=Y'

        psearchcl = PrimerSearchCommandline()
        psearchcl.seqall = f"{self.template} -snucleotide1 -sbegin1={self.begin} -send1={self.end} {sreverse}" 
        psearchcl.infile = self.primers
        psearchcl.mismatchpercent = self.mismatch
        psearchcl.outfile = "stdout"
        print(psearchcl, file=sys.stderr)
        stdout, stderr = psearchcl()
        self.pcr_results = psearch.read(StringIO(stdout))


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
                df = pd.DataFrame({'primer_pair': primerpair_name, 'amplicon_insert': None}, index=[primerpair_name])
                results_dfs_list.append(df)
            else:
                for index, amplimer in enumerate(amplifier): # need the indices
                    sub_df = {}
                    hit = amplimer.hit_info.replace('\t', '').rstrip(' ').split('\n')
                    hit = [i.strip() for i in hit]
                    #set primer_pari column
                    sub_df['primer_pair'] = primerpair_name
                    sub_df['amplimer_n'] = f"Amplimer_{index+1}"
                    sub_df['template_name'] = hit[0]
                    fwd = hit[-2].split(' ')
                    rev = hit[-1].split(' ')
                    sub_df['fwd_oligo'] = _collapsed_iupac(fwd[0]).upper()#, alphabet = IUPAC.ambiguous_dna)
                    sub_df['rev_oligo'] = _collapsed_iupac(rev[0]).upper()#, alphabet = IUPAC.ambiguous_dna)
                    sub_df['fwd_mismatches'] = int(fwd[-2]) 
                    sub_df['rev_mismatches'] = int(rev[-2])
                    sub_df['fwd_oligo_tmplt_start'] = int(fwd[-4]) - 1
                    sub_df['fwd_oligo_tmplt_end']   = sub_df['fwd_oligo_tmplt_start'] + len(sub_df['fwd_oligo'])
                    sub_df['rev_oligo_tmplt_end']   = len(self.template_seqs[hit[0]].seq) - int(rev[-4].replace('[', '').replace(']', '')) + 1
                    sub_df['rev_oligo_tmplt_start'] = sub_df['rev_oligo_tmplt_end'] - len(sub_df['rev_oligo'])
                    sub_df['amplicon_insert'] = str(self.template_seqs[sub_df['template_name']].seq[sub_df['fwd_oligo_tmplt_end']:sub_df['rev_oligo_tmplt_start']].upper())
                    sub_df['amplicon_full'] = str(self.template_seqs[sub_df['template_name']].seq[sub_df['fwd_oligo_tmplt_start']:sub_df['rev_oligo_tmplt_end']].upper())
                    sub_df['amplicon_length'] = len(sub_df['amplicon_full'])
                    sub_df['fwd_oligo_match'] = str(self.template_seqs[sub_df['template_name']].seq[sub_df['fwd_oligo_tmplt_start']:sub_df['fwd_oligo_tmplt_end']]).upper()
                    sub_df['rev_oligo_match'] = str(self.template_seqs[sub_df['template_name']].seq[sub_df['rev_oligo_tmplt_start']:sub_df['rev_oligo_tmplt_end']]).upper()
                    sub_df['product'] = str(sub_df['fwd_oligo']+sub_df['amplicon_insert']+Seq(sub_df['rev_oligo'], alphabet=IUPAC.ambiguous_dna).reverse_complement())
                    sub_df['fwd_match_mismatch'] = _iupac_zipper(sub_df['fwd_oligo'], sub_df['fwd_oligo_match'])
                    sub_df['rev_match_mismatch'] = _iupac_zipper(sub_df['rev_oligo'], str(Seq(sub_df['rev_oligo_match'], alphabet=IUPAC.ambiguous_dna).reverse_complement()))
                    df = pd.DataFrame(sub_df, index=[primerpair_name])
                    results_dfs_list.append(df)
        results = pd.concat(results_dfs_list, ignore_index=True)
        return results

# TODO convert float columns to integers