#!/usr/bin/env python3
import sys
import argparse
import analysis as ra

def main():
    Commands()


class Commands(object):
    def __init__(self):
        parser = argparse.ArgumentParser(
            prog="rnaindel",
            usage="""rnaindel <subcommand> [<args>]

subcommands are:
    analysis              Predict somatic indels from tumor RNA-Seq data
    featureCalculation               Calculate and report features for training""",
        )

        parser.add_argument(
            "subcommand",
            help="analysis, featureCalculation",
        )

        args = parser.parse_args(sys.argv[1:2])

        if not hasattr(self, args.subcommand):
            sys.exit("Error: invalid subcommand")

        getattr(self, args.subcommand)()

    def analysis(self):
        ra.analyze("analysis")

    def featureCalculation(self):
        ra.analyze("featureCalculation")



if __name__ == "__main__":
    main()
