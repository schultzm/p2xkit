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
        print(stdout)
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
            # sub_df = {}
            if not amplifier: # i.e., if list is empty
                df = pd.DataFrame({}, index=[primerpair_name])
                results_dfs_list.append(df)
            else:
                # pass
            # if len(amplifier) > 0:
            #     # sub_dfs_list = []
                # sub_df['primer_pair'] = primerpair_name
                # print(len(amplifier))
                for index, amplimer in enumerate(amplifier): # need the indices
                    sub_df = {}
                    # print('index', index)
                    # print('amplimer', amplimer.hit_info)
                    hit = amplimer.hit_info.replace('\t', '').rstrip(' ').split('\n')
                    hit = [i.strip() for i in hit]
                    # print(hit)
                    # sub_df['primer_pair'] = primerpair_name
                    sub_df['amplimer_n'] = index+1
                    sub_df['template_name'] = hit[0]
                    fwd = hit[-2].split(' ')
                    rev = hit[-1].split(' ')
                    # print(len(self.template_seqs[sub_df['template_name']].seq))

                    # print(str(self.template_seqs[sub_df['template_name']].seq))
                    # sub_df['template_seq'] = str(self.template_seqs[sub_df['template_name']].seq)
                    sub_df['forward_oligo'] = self._collapsed_iupac(fwd[0])#, alphabet = IUPAC.ambiguous_dna)
                    sub_df['reverse_oligo'] = self._collapsed_iupac(rev[0])#, alphabet = IUPAC.ambiguous_dna)
                    sub_df['forward_oligo_template_start'] = int(fwd[-4])
                    sub_df['reverse_oligo_template_end']   = len(self.template_seqs[sub_df['template_name']].seq) - int(rev[-4].replace('[', '').replace(']', ''))
                    sub_df['reverse_oligo_template_start'] = sub_df['reverse_oligo_template_end'] - len(sub_df['reverse_oligo'])
                                            # 'strand': amplimer.hit_info[1][2],
                                            # 'position': int(amplimer.hit_info[1][5]) - 1,
                                            # 'mismatches': amplimer.hit_info[1][7]}
                    # print(sub_df)
                    df = pd.DataFrame(sub_df, index=[primerpair_name])
                    results_dfs_list.append(df)
                    # print(sub_df)
                    # print()
        self.amplimer_tab = pd.concat(results_dfs_list)
                # print(primerpair_name, f"Amplimer {index}:", amplimer.hit_info)

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