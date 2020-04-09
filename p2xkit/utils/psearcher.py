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
        return a dictionary of {key(primer pair name [str]):
                                value ([list] of amplifiers)}
        """
        psearchcl = PrimerSearchCommandline()
        psearchcl.seqall = f"{self.template} -snucleotide1"
        psearchcl.infile = self.primers
        psearchcl.mismatchpercent = self.mismatch
        psearchcl.outfile = "stdout"
        stdout, stderr = psearchcl()
        return psearch.read(StringIO(stdout))

    def collapsed_iupac(self, primerstring):
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

    def amplimer(self, amplifier):
        '''
        An amplimer is a PCR product between a primer-pair.
        The amplifier is a dict, storing these data:
            keys are primer-pair names
            values are lists of amplifiers
                each amplifier in amplifiers has attributes:
                    hit_info (the amplimer)
                    length (of amplimer)
        '''
        amplicons = {primerpair: template_list for primerpair, template_list in amplifiers.items()}
        # return(amp_dict)
        for key, template_list in amp_dict.items:
            
        # for primerpair, template_list in amplifiers.items():
        #     print(primerpair, template_list)
        #     # For each 
        #     # enumerate in this loop to track (index) the amplimer (hit or "amplimer"")
            for index, amplimer in enumerate(template_list):
        #         # Reformat stdout from primersearch to a list for parsing
                amplimer.hit_info = amplimer.hit_info.replace('\t', '').strip().split('\n')
                amplimer.hit_info = [list(filter(None, i.split(' '))) for i in amplimer.hit_info]

                amplicon = {'hit_no': index}
                amplicon['template_name'] = amplimer.hit_info[0][0]
                amplicon['full_template_seq'] = seq[seq_idx[amplicon['template_name']]]
                amplicon['forward_oligo'] = {'primer_seq': Seq(collapse_iupac(amplimer.hit_info[1][0]), alphabet = IUPAC.ambiguous_dna),
                                        'strand': amplimer.hit_info[1][2],
                                        'position': int(amplimer.hit_info[1][5]) - 1,
                                        'mismatches': int(amplimer.hit_info[1][7])}
                amplicon['reverse_oligo'] = {'primer_seq': Seq(collapse_iupac(amplimer.hit_info[2][0]), alphabet = IUPAC.ambiguous_dna).reverse_complement(),
                                        'strand': amplimer.hit_info[2][2],
                                        'position': len(amplicon['full_template_seq'].seq) - int(amplimer.hit_info[2][5].replace('[', '').replace(']', '')) + 1,
                                        'mismatches': int(amplimer.hit_info[2][7])}
                amplicon['target_template_seq'] = amplicon['full_template_seq'][amplicon['forward_oligo']['position']:amplicon['reverse_oligo']['position']]
                amplicon['target_template_seq'].description = f"| {primerpair}~~~Amplimer {index}~~~{amplicon['forward_oligo']['position']}:{amplicon['reverse_oligo']['position']}"
                padded_primer_forward = SeqRecord(Seq(f"{amplicon['forward_oligo']['primer_seq']}{'-'*(len(amplicon['target_template_seq'])-len(amplicon['forward_oligo']['primer_seq']))}", alphabet = IUPAC.ambiguous_dna),#-len(amplimer['reverse_oligo']['primer_seq']))}{amplimer['reverse_oligo']['primer_seq']}"
                                        id = f'forward_oligo_{str(primerpair)}',
                                        description = str(primerpair))
                padded_primer_reverse = SeqRecord(Seq(f"{'-'*(len(amplicon['target_template_seq'])-len(amplicon['reverse_oligo']['primer_seq']))}{amplicon['reverse_oligo']['primer_seq']}", alphabet = IUPAC.ambiguous_dna),#-len(amplimer['reverse_oligo']['primer_seq']))}{amplimer['reverse_oligo']['primer_seq']}"
                                        id = f'reverse_oligo_{str(primerpair)}',
                                        description = str(primerpair))