import os
import sys
import tempfile
import argparse
from functools import partial

from rnaindel.defaultcaller.defaultcaller import callindel

from .preprocessor import preprocess
from .classifier import classify
from .postprocessor import postprocess
from .vcf_writer import write_vcf
from .utils import validate_int_inputs, validate_file_input, validate_dir_input


def analyze(subcommand, version=None):

    args = get_args(subcommand)

    bam = args.bam
    fasta = args.reference
    data_dir = args.data_dir.rstrip("/")
    mapq = args.uniq_mapq
    region = args.region
    external_vcf = args.vcf_file

    n_processes = 1 if region else args.process_num

    with tempfile.TemporaryDirectory() as tmp_dir:
        callindel(bam, fasta, tmp_dir, args.heap_memory, region, n_processes)

        df = preprocess(
            tmp_dir, fasta, bam, data_dir, mapq, n_processes, region, external_vcf, args.pass_only
        )
        if len(df) == 0:
            write_vcf(df, version, args)
            sys.exit(0)

    df = classify(df, "{}/models".format(data_dir), n_processes)

    df = postprocess(df, data_dir, args.perform_outlier_analysis, args.pon)

    write_vcf(df, version, args)


def get_args(subcommand):
    prog = "rnaindel " + subcommand
    parser = argparse.ArgumentParser(prog=prog)

    parser.add_argument(
        "-i",
        "--bam",
        metavar="FILE",
        required=True,
        type=validate_file_input,
        help="input tumor RNA-Seq BAM file (must be STAR-mapped).",
    )

    if subcommand == "PredictIndels":
        parser.add_argument(
            "-o", "--output-vcf", metavar="FILE", required=True, help="output VCF file"
        )

    parser.add_argument(
        "-r",
        "--reference",
        metavar="FILE",
        required=True,
        type=validate_file_input,
        help="reference genome FASTA file.",
    )

    parser.add_argument(
        "-d",
        "--data-dir",
        metavar="DIR",
        required=True,
        type=validate_dir_input,
        help="data directory contains databases and models",
    )

    parser.add_argument(
        "-v",
        "--vcf-file",
        metavar="FILE",
        default=None,
        type=validate_file_input,
        help="VCF file from external caller. Supply as vcf.gz + index.",
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

    parser.add_argument(
        "-m",
        "--heap-memory",
        metavar="STR",
        default="6000m",
        help="maximum heap space (default: 6000m)",
    )

    if subcommand == "PredictIndels":

        parser.add_argument(
            "--pon",
            metavar="FILE",
            default=None,
            type=validate_file_input,
            help="User defined panel of normals to refine somatic predictions. Supply as vcf.gz + index",
        )

        parser.add_argument(
            "--region",
            metavar="STR",
            default=None,
            help="specify region for target analysis: chrN:start-stop (default: None)",
        )

        parser.add_argument(
            "--include-all-external-calls",
            dest="pass_only",
            action="store_false",
            help="use all calls in external VCF (at default, use calls with PASS in FILTER)",
        )

        parser.set_defaults(pass_only=True)

        parser.add_argument(
            "--skip-homopolyer-outlier-analysis",
            dest="perform_outlier_analysis",
            action="store_false",
            help="skip oulier analysis to rescue somatic homopolyer indels (experimental feature)",
        )

        parser.set_defaults(perform_outlier_analysis=True)

    else:
        parser.add_argument(
            "-o",
            "--output-tab",
            metavar="FILE",
            required=True,
            help="output tab-delimited file",
        )

    return parser.parse_args(sys.argv[2:])
