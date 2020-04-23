
import os
import sys
import shlex
from subprocess import Popen, PIPE
import pysam
from collections import defaultdict
from pathlib import Path, PurePath
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import pandas as pd
from ..utils.psearcher import _iupac_zipper

# Moved to module scope
def indexit(infile):
    # Template bowtie2 index files generation
    template_bowtie2_idx_fnames = list(infile.parent.glob(f"{infile.stem}.*.bt2"))
    if len(template_bowtie2_idx_fnames) != 6:
        cmd = f"bowtie2-build --threads 4 -q -f {infile} {PurePath(infile.parent, infile.stem)}"
        os.system(cmd)
    # return list(infile.parent.glob(f"{infile.stem}.*.bt2"))


class Bowtier:
    def __init__(self, amplimer_table, probes):
        self.probes = probes
        self.amplimer_table = amplimer_table

    def bowtieit(self):
        print(self.amplimer_table)
        probes_dict = defaultdict(list)
        with open(self.probes, 'r') as input_handle:
            probes = list(SeqIO.parse(input_handle, 'fasta'))
            for probe in probes:
                probes_dict[probe.id].append(probe) # Need to specify that probe names must be same as primer_pair with a space and then a probe identifier (e.g., 'RdRP_SARSr_DE P2')
        for rown in self.amplimer_table.index.values:
            if pd.notnull(self.amplimer_table.loc[rown, 'amplicon_insert']):
                subseq = SeqRecord(Seq(self.amplimer_table.loc[rown, 'amplicon_insert'],
                                       alphabet=IUPAC.ambiguous_dna),
                                id=self.amplimer_table.loc[rown, 'primer_pair'],
                                name=self.amplimer_table.loc[rown, 'template_name'],
                                description=f"Amplimer_{self.amplimer_table.loc[rown, 'amplimer_n']:.0f}")
                outhandle = Path(f"{subseq.id}_{subseq.description}.fasta")
                SeqIO.write(subseq, outhandle, 'fasta') #more thought needs to be put into this
                indexit(outhandle) # need to clean this up at end

                map_results_dfs = [] 
                for probe in probes_dict[subseq.id]:
                    probe_id = probe.description.split(' ')[-1] #this is risky if space is not used to delimit probe ID
                    print(probe_id)
                    # Bowtie2 summary of options used
                    # -L <int>           length of seed substrings; must be >3, <32 (22)
                    # -c                 <m1>, <m2>, <r> are sequences themselves, not files
                    # -U unpairedreads
                    # -a/--all           report all alignments; very slow, MAPQ not meaningful # don't use this
                    # -N <int>           max # mismatches in seed alignment; can be 0 or 1 (0)
                    # --np <int>         penalty for non-A/C/G/Ts in read/ref (1)
                    # -R <int>           for reads w/ repetitive seeds, try <int> sets of seeds (2)
                    # -f                 query input files are (multi-)FASTA .fa/.mfa
                    # --sam-no-qname-trunc Suppress standard behavior of truncating readname at first whitespace 
                    #                      at the expense of generating non-standard SAM.
                    map_cmd = f"bowtie2 -x {PurePath(outhandle.parent, outhandle.stem)} -U {probe.seq} -c --sam-no-qname-trunc --end-to-end  -L 7 -D 20"#q --np 0 -R 10"#, tqmanprobe.description) for tqmanprobe in probes_list]
                    proc1 = Popen(shlex.split(map_cmd), stdout=PIPE, stderr=PIPE)
                    samfile = proc1.stdout.fileno()
                    with pysam.AlignmentFile(samfile, "r") as sam:
                        for rec in sam.fetch():
                            if not rec.is_unmapped:
                                probe_range = {'probe_template_start': rec.get_aligned_pairs()[0][1],
                                               'probe_template_end'  : rec.get_aligned_pairs()[-1][1]}
                                keys_to_keep = ['name', 'flag', 'ref_name', 'ref_pos'] # SAM is 1-based, so ref_pos will be probe_template_start+1
                                records = rec.to_dict()
                                records_subset = {key: value for key, value in records.items() if key in keys_to_keep}
                                records_subset = {**probe_range, **records_subset}
                                records_subset['probe_length'] = len(probe.seq)
                                records_subset['probe_length_aligned'] = records_subset['probe_template_end']-records_subset['probe_template_start']+1 #+1 as e.g., template matches at 21,22,23,24,25 are 5 matches but 25-21=4
                                records_subset['probe_globally_aligned'] = ''.join(['True' if records_subset['probe_length_aligned']==records_subset['probe_length'] else 'False'])
                                records_subset['probe_seq'] = str(probe.seq)
                                records_subset['probe_orientation'] = 'FORWARD'
                                records_subset['probe_match'] = subseq[records_subset['probe_template_start']:records_subset['probe_template_end']+1].seq
                                records_subset['primer_pair'] = probe.id
                                if rec.flag == 16: # REVERSED
                                    records_subset['probe_orientation'] = 'REVERSE'
                                    records_subset['probe_match'] = str(records_subset['probe_match'].reverse_complement()).upper()
                                else:
                                    records_subset['probe_match'] = str(records_subset['probe_match']).upper()
                                records_subset['probe_match_mismatch'] = _iupac_zipper(records_subset['probe_seq'], records_subset['probe_match'])
                                sub_df = pd.DataFrame(records_subset, index=[records_subset['name']])
                                print('SUBDF', sub_df.to_csv(sep="\t"))
                                print('!!!', self.amplimer_table.loc[self.amplimer_table['primer_pair'] == records_subset['primer_pair']].to_csv(sep="\t"))
                                map_results_dfs.append(sub_df)
