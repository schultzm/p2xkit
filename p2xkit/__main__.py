#!/usr/bin/env python3

# from p2xkit.utils.oligo_pair import Oligo_pair
from Bio.Emboss import  PrimerSearch as psearch
from Bio.Emboss.Applications import PrimerSearchCommandline
from io import StringIO
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from Bio.Alphabet import IUPAC
# from Bio.Align import MultipleSeqAlignment as MSA
# import os
# import sys
# import shlex
# from subprocess import Popen, PIPE
# import pysam
# from collections import defaultdict
# print(sys.argv)
# psearchcl = PrimerSearchCommandline()
# psearchcl.seqall = f"{sys.argv[1]} -snucleotide1"
# psearchcl.infile = sys.argv[2]
# psearchcl.mismatchpercent = 20
# psearchcl.outfile = "stdout"
# print(psearchcl, file = sys.stderr)
# stdout, stderr = psearchcl()
# primersearch_results = psearch.read(StringIO(stdout))

# # Template seqs with key as seq.id
# template_seqrecords = None
# with open(sys.argv[1]) as input_handle:
#     seqs = {seq.id: seq for seq in list(SeqIO.parse(input_handle, 'fasta'))}

# # Probe seqs with key as primerpair_name.  
# # Each probe MUST have same name as corresponding primerpair
# probes = None
# with open(sys.argv[3]) as input_handle:
#     probes = {probe.id: probe for probe in list(SeqIO.parse(input_handle, 'fasta'))}


# print(amplicon_dict_builder(primersearch_results.amplifiers, template_seqrecords))

# Process emboss primersearch 'results' dict:
    # holds template_list under primerpair keys
        # Create a dict of amplimer details

        # Write the template strand to fasta if id doesn't already exist
        # template_fasta = Path(PurePath(sys.argv[1]))
        # # if not template_fasta.is_file():
        # #     with open(template_fasta, "w") as outhandle:
        # #         SeqIO.write(amplimer['full_template_seq'], outhandle, 'fasta')
        # # Template strand bowtie2 index files generation
        # template_bowtie2_idx_fnames = list(Path(PurePath(template_fasta).parent).glob(f"{template_fasta.stem}.*.bt2"))
        # if len(template_bowtie2_idx_fnames) != 6:
        #     os.system(f"bowtie2-build -f {template_fasta} {PurePath(template_fasta.parent, template_fasta.stem)}")

        # else:
        #     # primerpair is the primer_pair name
        #     # idx is the probe number in an amplimer
        #     # probe_number is lookup index of probe in probes SeqRecord list
        #     for idx, probe_number in enumerate(probes_idx[primerpair]):
        #         probe = probes[probe_number]
        #         bowtie2_map_cmd = shlex.split(f"bowtie2 -x {PurePath(template_fasta.parent, template_fasta.stem)} -U {probe.seq} -c -a --end-to-end --very-sensitive")
        #         # need to move the bowtie step outside the loop, create a dict of dicts and access subdicts by probe name then template name
        #         proc1 = Popen(bowtie2_map_cmd, stdout=PIPE, stderr=PIPE)
        #         samfile = proc1.stdout.fileno()
        #         probe_mapped = None
        #         with pysam.AlignmentFile(samfile, "r") as sam:
        #             # bt2_sams = 
        #             # print(sam)
        #             # print(dir(sam))
        #             for rec in sam.fetch():
        #                 # print(dir(rec))
        #                 # print(type(rec))
        #                 # print(rec.to_dict())
        #                 # sys.exit()
        #                 # bt2_hit = 
        #                 if not rec.is_unmapped:
        #                     probe_posn = rec.aligned_pairs[0][-1] # the first ref pos in the aligned pair
        #                     left_pad   = '-' * (probe_posn - amplicon['forward_oligo']['position'])
        #                     right_pad  = '-' * (len(amplicon['target_template_seq']) - (len(left_pad)+len(probe)))
        #                     seq_str = None
        #                     id = None
        #                     if rec.flag == 16:
        #                         seq_str = probe.seq.reverse_complement()
        #                         id = f'probe_{idx}_rcomp_{str(primerpair)}'
        #                     else:
        #                         seq_str = probe.seq
        #                         id = f'probe_{idx}_{str(primerpair)}'
        #                     probe_mapped = SeqRecord(Seq(f"{left_pad}{seq_str}{right_pad}",
        #                                                 alphabet = IUPAC.ambiguous_dna),
        #                                             id = id,
        #                                             description = primerpair)
        #                     # print(amplimer['target_template_seq'].seq)
        #                     # print(str(padded_primer_forward.seq),'\n',
        #                     #         str(padded_primer_reverse.seq),'\n',
        #                     #         str(probe_mapped.seq),'\n')
        #                     aln = MSA([amplicon['target_template_seq'],
        #                             padded_primer_forward,
        #                             padded_primer_reverse,
        #                             probe_mapped])
        #                     outdir = Path(PurePath(Path.cwd(), sys.argv[4], amplicon['target_template_seq'].id, primerpair, f'probe_{idx}'))
        #                     # print(outdir)
        #                     outdir.mkdir(parents=True, exist_ok=True)
        #                     with open(PurePath(outdir, f"amplimer_{amplicon['hit_no']}.fa"), "w") as outhandle:
        #                         AlignIO.write(aln, outhandle, 'fasta')
        #                         # print(aln.format('fasta'))

def main():
    """Perform the main routine."""
    import argparse
    from pathlib import Path, PurePath
    import p2xkit

    parser = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            description="""
                        """)
    subparser1_args = argparse.ArgumentParser(add_help=False)
    subparser1_args.add_argument("template", help = "Template fasta")
    subparser1_args.add_argument("primers", help = "primersearch formatted tab file")
    subparser1_args.add_argument("mismatch", help = "Percent mismatch tolerance")
    subparser1_args.add_argument("probes", help = "Fasta formatted qPCR probes file",)
    subparser_modules = parser.add_subparsers(
        title="Sub-commands help", help="", metavar="", dest="subparser_name")
    subparser_modules.add_parser(
        "bed", help="""Map the primers and probes to ref output
                       bed file.""",
        description="Create a bed file of primer and probe map.",
        parents=[subparser1_args],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparser_modules.add_parser(
        "test", help="""Run p2xkit unittests""",
        description="Run p2xkit unittests",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    args = parser.parse_args()

    if not args.subparser_name:
        parser.print_help()
    elif args.subparser_name == "bed":
        pass
        # from p2xkit import __version__
    elif args.subparser_name == "version":
        print(p2xkit.__version__)
    elif args.subparser_name == "test":
        import unittest
        from .tests.test_suite import suite
        runner = unittest.TextTestRunner(verbosity=2)
        runner.run(suite())


if __name__ == "__main__":
    main()