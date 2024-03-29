#!/usr/bin/env python3
import sys
import argparse
import rnaindel.analysis as ra
import rnaindel.training as tr
import rnaindel.occurrence as oc

from .version import __version__


def main():
    Commands()


class Commands(object):
    def __init__(self):
        parser = argparse.ArgumentParser(
            prog="rnaindel",
            usage="""rnaindel <subcommand> [<args>]

subcommands are:
    SetUp                     Initialize prediction models
    PredictIndels             Predict somatic/germline/artifact indels from tumor RNA-Seq data
    CalculateFeatures         Calculate and report features for training
    Train                     Perform model training
    CountOccurrence           Find overrepresented somatic predictions within cohort to filter false positives""",
        )

        parser.add_argument(
            "subcommand",
            help="SetUp, PredictIndels, CalculateFeatures, Train, CountOccurrence",
        )

        parser.add_argument(
            "--version",
            action="version",
            version="%(prog)s {version}".format(version=__version__),
        )

        args = parser.parse_args(sys.argv[1:2])

        if not hasattr(self, args.subcommand):
            sys.exit("Error: invalid subcommand")

        getattr(self, args.subcommand)()

    def SetUp(self):
        tr.setup()
    
    def PredictIndels(self):
        ra.analyze("PredictIndels", version=__version__)

    def CalculateFeatures(self):
        ra.analyze("CalculateFeatures")

    def Train(self):
        tr.train()

    def CountOccurrence(self):
        oc.count()


if __name__ == "__main__":
    sys.exit(main())
