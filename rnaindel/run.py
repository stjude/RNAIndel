#!/usr/bin/env python3
import sys
import argparse
import analysis as ra
from version import __version__


def main():
    Commands()


class Commands(object):
    def __init__(self):
        parser = argparse.ArgumentParser(
            prog="rnaindel",
            usage="""rnaindel <subcommand> [<args>]

subcommands are:
    PredictSomaticIndels            Predict somatic indels from tumor RNA-Seq data
    CalculateFeatures               Calculate and report features for training""",
        )

        parser.add_argument(
            "subcommand",
            help="PredictSomaticIndels, CalculateFeatures",
        )

        args = parser.parse_args(sys.argv[1:2])

        if not hasattr(self, args.subcommand):
            sys.exit("Error: invalid subcommand")

        getattr(self, args.subcommand)()
         
    def PredictSomaticIndels(self):
        ra.analyze("PredictSomaticIndels", version=__version__)

    def CalculateFeatures(self):
        ra.analyze("CalculateFeatures")



if __name__ == "__main__":
    main()
