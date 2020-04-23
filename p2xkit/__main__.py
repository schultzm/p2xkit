#!/usr/bin/env python3


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
    subparser1_args.add_argument("template", help = "Template fasta.")
    subparser1_args.add_argument("primers", help = "primersearch formatted tab-delimited primer file")
    subparser1_args.add_argument('-m', "--mismatch", help = """Percent mismatch
                                 identity tolerance""",
                                 type=int, default=20, required=False)
    subparser2_args = argparse.ArgumentParser(add_help=False)
    subparser2_args.add_argument("probes", help = "Fasta formatted qPCR probes file.",)
    subparser_modules = parser.add_subparsers(
        title="Sub-commands help", help="", metavar="", dest="subparser_name")
    subparser_modules.add_parser(
        "ispcr", help="""Perform in-silico PCR.""",
        description="Perform in-silico PCR.",
        parents=[subparser1_args],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparser_modules.add_parser(
        "qpcr", help="""Perform in-silico qPCR.""",
        description="Perform in-silico qPCR.",
        parents=[subparser1_args, subparser2_args],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparser_modules.add_parser(
        "test", help="""Run p2xkit unittests""",
        description="Run p2xkit unittests",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparser_modules.add_parser(
        "version", help="""Print version to stdout""",
        description="Print version to stdout.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    args = parser.parse_args()

    if not args.subparser_name:
        parser.print_help()
    elif args.subparser_name == "ispcr":
        from .utils.psearcher import Psearcher
        reaction = Psearcher(args.template, args.primers, args.mismatch)
        reaction.psearchit()
        print(reaction.amplimer_table().to_csv(sep="\t"))
    elif args.subparser_name == "qpcr":
        from .utils.psearcher import Psearcher
        from .utils.bowtier import Bowtier
        reaction = Psearcher(args.template, args.primers, args.mismatch)
        reaction.psearchit()
        amplimer_table = reaction.amplimer_table()
        qpcr_setup = Bowtier(amplimer_table, args.probes)
        print(qpcr_setup.bowtieit().to_csv(sep="\t"))
    elif args.subparser_name == "version":
        print(p2xkit.__version__)
    elif args.subparser_name == "test":
        import unittest
        from .tests.test_suite import suite
        runner = unittest.TextTestRunner(verbosity=2)
        runner.run(suite())


if __name__ == "__main__":
    main()