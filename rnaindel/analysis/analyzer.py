import re
import os
import sys
import pysam
import sklearn
import tempfile
import argparse
import pandas as pd
from functools import partial
from collections import Counter

from rnaindel.defaultcaller.defaultcaller import callindel
from rnaindel.defaultcaller.softclip_realigner import realn_softclips

from .preprocessor import preprocess
from .classifier import classify
from .postprocessor import postprocess
from .cross_checker import cross_check
from .vcf_writer import write_vcf
from .utils import (
    validate_int_inputs,
    validate_file_input,
    validate_dir_input,
    validate_str_input,
)


def analyze(subcommand, version=None):
    args = get_args(subcommand)

    bam = args.bam
    fasta = args.reference
    data_dir = args.data_dir.rstrip("/")
    mapq = args.uniq_mapq
    region = args.region
    external_vcf = args.vcf_file
    tdna = args.tumor_dna
    ndna = args.normal_dna

    setup_check(data_dir)

    n_processes = 1 if region else args.process_num

    check_cosmic(data_dir)

    read_len = read_len_check(bam)
    canonicals = canonical_chromosomes(bam)
    with tempfile.TemporaryDirectory() as tmp_dir:
        callindel(
            bam, fasta, tmp_dir, args.heap_memory, canonicals, region, n_processes
        )
        if not args.deactivate_sensitive_mode:
            realn_softclips(
                bam,
                fasta,
                tmp_dir,
                data_dir,
                canonicals,
                region,
                n_processes,
                args.safety_mode,
                read_len,
            )

        df = preprocess(
            tmp_dir,
            fasta,
            bam,
            data_dir,
            mapq,
            n_processes,
            canonicals,
            region,
            external_vcf,
            args.pass_only,
            args.safety_mode,
        )
        if len(df) == 0:
            write_vcf(df, version, args, tdna, ndna)
            sys.exit(0)

    df = classify(df, "{}/models".format(data_dir), n_processes)

    df = postprocess(df, data_dir, args.perform_outlier_analysis, args.pon)

    if tdna or ndna:
        df = cross_check(df, fasta, tdna, ndna)

    write_vcf(df, version, args, tdna, ndna)


def setup_check(data_dir):
    setup_log = "{}/models/setuplog.txt".format(data_dir)
    if os.path.isfile(setup_log):
        with open(setup_log) as f:
            curr_v = str(sklearn.__version__)
            setup_v = str(f.read().rstrip())
            if curr_v != setup_v:
                print(
                    "Warning: Consider running SetUp command, rnaindel SetUp -d /path/to/data_dir",
                    file=sys.stderr,
                )
    else:
        print(
            "Error: Run SetUp command, rnaindel SetUp -d /path/to/data_dir",
            file=sys.stderr,
        )
        sys.exit(1)


def read_len_check(bam):
    i = 0
    lens = []
    reads = pysam.AlignmentFile(bam).fetch()
    for read in reads:
        if i > 99:
            break
        lens.append(read.query_alignment_length)
        i += 1
    if lens:
        return Counter(lens).most_common(1)[0][0]
    return 75


def canonical_chromosomes(bam):
    chromosomes = pysam.AlignmentFile(bam).references
    canonicals = []
    for chrom in chromosomes:
        chrom = chrom.replace("chr", "")
        if chrom in ["X", "Y"]:
            canonicals.append(chrom)
        elif re.match(r"\d+", chrom):
            canonicals.append(str(chrom))
    return canonicals


def check_cosmic(data_dir):
    path = os.path.join(data_dir, "cosmic")

    # non-human support
    if not os.path.isdir(path):
        return

    vcf_gz_ok = False
    gz_tbi_ok = False
    for _ in os.listdir(path):
        if not vcf_gz_ok and _.endswith(".vcf.gz"):
            vcf_gz_ok = True

        if not gz_tbi_ok and _.endswith(".gz.tbi"):
            gz_tbi_ok = True

    if vcf_gz_ok and gz_tbi_ok:
        pass
    else:
        print(
            "ERROR: check COSMIC VCF file and index under {}/cosmic".format(data_dir),
            file=sys.stderr,
        )
        sys.exit(1)


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
        help="number of processes (default: 1)",
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
        type=validate_str_input,
        help="maximum heap space (default: 6000m)",
    )

    if subcommand == "PredictIndels":
        parser.add_argument(
            "-t",
            "--tumor-dna",
            metavar="FILE",
            default=None,
            type=validate_file_input,
            help="Tumor DNA-Seq BAM file for cross-platform check (default: None)",
        )

        parser.add_argument(
            "-n",
            "--normal-dna",
            metavar="FILE",
            default=None,
            type=validate_file_input,
            help="Normal DNA-Seq BAM file for cross-platform check (default: None)",
        )

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

        parser.add_argument(
            "--safety-mode",
            dest="safety_mode",
            action="store_true",
            help="deactivate parallelism at realignment step (experimental feature)",
        )

        parser.set_defaults(safety_mode=False)

        parser.add_argument(
            "--deactivate-sensitive-mode",
            dest="deactivate_sensitive_mode",
            action="store_true",
            help="deactivate additional realignments around softclips (experimental feature)",
        )

        parser.set_defaults(deactivate_sensitive_mode=False)
    else:
        parser.add_argument(
            "-o",
            "--output-tab",
            metavar="FILE",
            required=True,
            help="output tab-delimited file",
        )

    return parser.parse_args(sys.argv[2:])
