
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
    return list(infile.parent.glob(f"{infile.stem}.*.bt2"))


class Bowtier:
    def __init__(self, amplimer_table, template, probes, reverse_complement):
        self.template = template
        self.probes = probes
        self.amplimer_table = amplimer_table
        self.reverse_complement = reverse_complement
        self.template_seqs = {seq.id: seq for seq in list(SeqIO.parse(open(self.template, 'r'), 'fasta'))}
        if self.reverse_complement:
            self.template_seqs = {seq.id: seq.reverse_complement() for seq in list(SeqIO.parse(open(self.template, 'r'), 'fasta'))}


    def bowtieit(self):
        self.amplimer_table.set_index('primer_pair')
        probes_dict = defaultdict(list)
        total_df = []
        with open(self.probes, 'r') as input_handle:
            probes = list(SeqIO.parse(input_handle, 'fasta'))
            for probe in probes:
                probes_dict[probe.id].append(probe) # Need to specify in readme that probe names must be same as primer_pair with a space and then a probe identifier (e.g., 'RdRP_SARSr_DE P2')
        for rown in self.amplimer_table.index.values: #iterate through all the amplimers and map probes
            amplimer_n = f"{self.amplimer_table.loc[rown, 'amplimer_n']}"
            amplicon_insert = self.template_seqs[self.amplimer_table.loc[rown, 'template_name']]. \
                              seq[self.amplimer_table.loc[rown, 'fwd_oligo_tmplt_end']: \
                                  self.amplimer_table.loc[rown, 'rev_oligo_tmplt_start']].upper()
            subseq = SeqRecord(amplicon_insert,
                              id=self.amplimer_table.loc[rown, 'primer_pair'],
                              description=f"{amplimer_n} from {self.amplimer_table.loc[rown, 'template_name']}")
            outhandle = Path(f"{self.amplimer_table.loc[rown, 'template_name']}_primerpair{subseq.id}_{amplimer_n}.fasta")
            SeqIO.write(subseq, outhandle, 'fasta') # writes out the amplicon insert to file
            indexed = indexit(outhandle) # indexes the amplicon insert **need to clean this up at end
            for probe in probes_dict[subseq.id]: #Iterate through probes for each primer pair
                map_cmd = f"bowtie2 -x {PurePath(outhandle.parent, outhandle.stem)} -U {probe.seq} -c --sam-no-qname-trunc --end-to-end  -L 7 -D 20"
                proc1 = Popen(shlex.split(map_cmd), stdout=PIPE, stderr=PIPE)
                samfile = proc1.stdout.fileno()
                with pysam.AlignmentFile(samfile, "r") as sam:
                    map_results_dfs = [] # this will store the dfs of mappings
                    for rec in sam.fetch(): # rec is a row in the SAM output
                        if not rec.is_unmapped:
                            records_subset = {'probe_template_start': rec.get_aligned_pairs()[0][1], # these indices are 0-based
                                            'probe_template_end'  : rec.get_aligned_pairs()[-1][1]} # 0-based
                            records_subset['probe_SAM_flag'] = rec.flag
                            records_subset['probe_name'] = probe.description
                            records_subset['probe_id'] = probe.description.split(' ')[-1] #e.g., "P", "P1" or "P2"; this is risky if space is not used to delimit probe ID
                            records_subset['template_name'] = subseq.description.split(' ')[-1]
                            records_subset['primer_pair'] = probe.id
                            records_subset['amplimer_n'] = amplimer_n
                            records_subset['probe_length'] = len(probe.seq)
                            records_subset['probe_seq'] = str(probe.seq)
                            records_subset['probe_length_aligned'] = records_subset['probe_template_end']-records_subset['probe_template_start']+1 #+1 as e.g., number of template matches at 21,22,23,24,25 are 5 but 25-21=4
                            records_subset['probe_globally_aligned'] = ''.join(['True' if records_subset['probe_length_aligned']==records_subset['probe_length'] else 'False'])
                            records_subset['probe_orientation'] = 'FORWARD' # gets converted to reverse if samflag is 16, below
                            records_subset['probe_match'] = subseq[records_subset['probe_template_start']:records_subset['probe_template_end']+1].seq #+1 as this is a slice index; this a Seq object
                            if rec.is_reverse: # REVERSED, do some revcomp and flagging
                                records_subset['probe_orientation'] = 'REVERSE'
                                records_subset['probe_match'] = str(records_subset['probe_match'].reverse_complement()).upper()
                            else: # convert match to string
                                records_subset['probe_match'] = str(records_subset['probe_match']).upper()
                            records_subset['probe_match_mismatch'] = _iupac_zipper(records_subset['probe_seq'], records_subset['probe_match'])

                            sub_df = pd.DataFrame(records_subset, index=[probe.id]) #probe.id is primer_pair
                            to_join = self.amplimer_table.loc[(self.amplimer_table['primer_pair'] == records_subset['primer_pair']) & \
                                                                (self.amplimer_table['amplimer_n'] == records_subset['amplimer_n']) & \
                                                                (self.amplimer_table['template_name'] == records_subset['template_name'])]
                            to_join.set_index('primer_pair', inplace=True)
                            to_join = to_join[[column for column in to_join.columns if column not in sub_df.columns]]
                            output_df = pd.concat([sub_df, to_join], axis=1, join='inner')
                            map_results_dfs.append(output_df)
                        else:
                            sub_df = pd.DataFrame({}, index=[probe.id]) #probe.id is primer_pair
                            to_join = self.amplimer_table.loc[(self.amplimer_table['primer_pair'] == probe.id) & \
                                                                (self.amplimer_table['amplimer_n'] == amplimer_n) & \
                                                                (self.amplimer_table['template_name'] == subseq.description.split(' ')[-1])]
                            to_join.set_index('primer_pair', inplace=True)
                            to_join = to_join[[column for column in to_join.columns if column not in sub_df.columns]]
                            output_df = pd.concat([sub_df, to_join], axis=1, join='inner')
                            map_results_dfs.append(output_df)
                    df = pd.concat(map_results_dfs)
                    total_df.append(df)
            for i in indexed:
                i.unlink() #remove all the index files
            outhandle.unlink() # remove the subseq.fa
        if total_df:
            probes_mapped_table = pd.concat(total_df)#.to_csv(sep="\t"))
            return probes_mapped_table
        else:
            return pd.DataFrame({}, index=['NO qPCR HITS FOUND.  Run programname ispcr cmds instead to check for stage 1 PCR hits']) # TODO: update this