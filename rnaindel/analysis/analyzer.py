import os
import sys
import tempfile
import argparse
from functools import partial

from defaultcaller.defaultcaller import callindel
from .utils import validate_int_inputs


from .preprocessor import preprocess
from .classifier import classify
from .postprocessor import postprocess
from .vcf_writer import write_vcf


def analyze(subcommand, version=None):

    args = get_args(subcommand)

    bam = args.bam
    fasta = args.reference
    data_dir = args.data_dir.rstrip("/")
    mapq = args.uniq_mapq
    
    n_processes = 1 if args.region else args.process_num

    with tempfile.TemporaryDirectory() as tmp_dir:
        callindel(bam, fasta, tmp_dir, args.heap_memory, args.region, n_processes)

        df = preprocess(tmp_dir, fasta, bam, data_dir, mapq, n_processes)

    df = classify(df, "{}/models".format(data_dir), n_processes)

    df = postprocess(df, data_dir)
   
    write_vcf(df, version, args)

    #df.to_csv(args.output_vcf, sep="\t", index=False)


def get_args(subcommand):
    prog = "rnaindel " + subcommand
    parser = argparse.ArgumentParser(prog=prog)

    parser.add_argument(
        "-i",
        "--bam",
        metavar="FILE",
        required=True,
        help="input tumor RNA-Seq BAM file (must be STAR-mapped).",
    )

    parser.add_argument(
        "-r",
        "--reference",
        metavar="FILE",
        required=True,
        help="reference genome FASTA file.",
    )

    parser.add_argument(
        "-d",
        "--data-dir",
        metavar="DIR",
        required=True,
        help="data directory contains databases and models",
    )

    parser.add_argument(
        "-m",
        "--heap-memory",
        metavar="STR",
        default="6000m",
        help="maximum heap space (default: 6000m)",
    )

    parser.add_argument(
        "-p",
        "--process-num",
        metavar="INT",
        default=1,
        type=validate_int_inputs,
        help="number of processes (defaul: 1)",
    )

    parser.add_argument(
        "-q",
        "--uniq-mapq",
        metavar="INT",
        default=255,
        type=partial(validate_int_inputs, preset="mapq"),
        help="STAR mapping quality MAPQ for unique mappers (default: 255)",
    )

    if subcommand == "PredictSomaticIndels":
        parser.add_argument(
            "-o", "--output-vcf", metavar="FILE", required=True, help="output VCF file"
        )

        parser.add_argument(
            "--region",
            metavar="STR",
            default=None,
            help="specify region for target analysis: chrN:start-stop (default: None)",
        )
    else:
        parser.add_argument(
            "-o",
            "--output-tab",
            metavar="FILE",
            required=True,
            help="output tab-delimited file",
        )

    return parser.parse_args(sys.argv[2:])
