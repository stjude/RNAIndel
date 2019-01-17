#!/usr/bin/env python3 

import os
import sys
import pathlib
import logging
import warnings
import argparse
import tempfile
import pandas as pd
from functools import partial
from .version import __version__

import rnaindel.bambino_lib as bl
import rnaindel.rnaindel_lib as rl

warnings.filterwarnings("ignore", category=UserWarning)

def main():
    args = get_args()
    create_logger(args.log_dir)
    data_dir = args.data_dir.rstrip("/")
    refgene = "{}/refgene/refCodingExon.bed.gz".format(data_dir)
    dbsnp = "{}/dbsnp/dbsnp.indel.vcf.gz".format(data_dir)
    clinvar = "{}/clinvar/clinvar.indel.vcf.gz".format(data_dir)
    model_dir = "{}/models".format(data_dir)
    
    # Preprocessing 
    # Variant calling will be performed if no external VCF is supplied
    if not args.input_vcf:
        # indel calling
        bambino_output = os.path.join(tempfile.mkdtemp(), "bambino.txt")
        bl.bambino(args.bam, args.fasta, bambino_output, args.heap_memory)
        
        # preprocess indels from the built-in caller
        df, chr_prefixed = rl.indel_preprocessor(
            bambino_output, args.bam, refgene, args.fasta
        )
        df = rl.indel_rescuer(
            df, args.fasta, args.bam, chr_prefixed, num_of_processes=args.process_num
        )
        
        # delete the temp file
        os.remove(bambino_output)
    else:
        # preprocess indels from external VCF
        df, chr_prefixed = rl.indel_vcf_preprocessor(
            args.input_vcf, args.bam, refgene, args.fasta
        )
        df = rl.indel_rescuer(
            df,
            args.fasta,
            args.bam,
            chr_prefixed,
            num_of_processes=args.process_num,
            left_aligned=True,
            external_vcf=True,
        )

    # Analysis 1: indel annotation
    df = rl.indel_annotator(df, refgene, args.fasta, chr_prefixed)
    # Analysis 2: feature calculation using
    df, df_filtered_premerge = rl.indel_sequence_processor(
        df, args.fasta, args.bam, args.uniq_mapq, chr_prefixed
    )
    df = rl.indel_protein_processor(df, refgene)
    # Analysis 3: merging equivalent indels
    df, df_filtered_postmerge = rl.indel_equivalence_solver(
        df, args.fasta, refgene, chr_prefixed
    )
    # Analysis 4: dbSNP annotation
    df = rl.indel_snp_annotator(df, args.fasta, dbsnp, clinvar, chr_prefixed)
    # Analysis 5: prediction
    df = rl.indel_classifier(df, model_dir, num_of_processes=args.process_num)

    # Analysis 6: concatenating invalid(filtered) entries
    df_filtered = pd.concat(
        [df_filtered_premerge, df_filtered_postmerge],
        axis=0,
        ignore_index=True,
        sort=True,
    )

    # Analysis 7(Optional): custom refinement of somatic prediction
    if args.non_somatic_panel:
        df = rl.indel_reclassifier(df, args.fasta, chr_prefixed, args.non_somatic_panel)

    # PostProcessing & VCF formatting
    df, df_filtered = rl.indel_postprocessor(
        df, df_filtered, refgene, args.fasta, chr_prefixed
    )
    rl.indel_vcf_writer(
        df,
        df_filtered,
        args.bam,
        args.fasta,
        chr_prefixed,
        args.output_vcf,
        __version__,
    )

    print("rnaindel completed successfully.", file=sys.stderr)


def get_args():
    parser = argparse.ArgumentParser(prog="rnaindel")
    parser.add_argument(
        "-b",
        "--bam",
        metavar="FILE",
        required=True,
        type=partial(check_file, file_name="BAM file (.bam)"),
        help="input tumor RNA-Seq bam file (must be STAR-mapped).",
    )
    parser.add_argument(
        "-f",
        "--fasta",
        metavar="FILE",
        required=True,
        type=partial(check_file, file_name="FASTA file"),
        help="reference genome FASTA file.",
    )
    parser.add_argument(
        "-d",
        "--data-dir",
        metavar="DIR",
        required=True,
        help="data directory contains refgene, dbsnp and clinvar databases and models",
        type=check_folder_existence,
    )
    parser.add_argument(
        "-o", "--output-vcf", metavar="FILE", required=True, help="output vcf file"
    )
    # input VCF from other callers (optional)
    parser.add_argument(
        "-c",
        "--input-vcf",
        metavar="FILE",
        type=partial(check_file, file_name="VCF (.vcf) file"),
        help="input vcf file from other callers",
    )
    parser.add_argument(
        "-q",
        "--uniq-mapq",
        metavar="INT",
        default=255,
        type=check_mapq,
        help="STAR mapping quality MAPQ for unique mappers (default: 255)",
    )
    parser.add_argument(
        "-p",
        "--process-num",
        metavar="INT",
        default=1,
        type=check_pos_int,
        help="number of processes (default: 1)",
    )
    parser.add_argument(
        "-n",
        "--non-somatic-panel",
        metavar="FILE",
        type=partial(check_file, file_name="Panel of non-somatic (.vcf)"),
        help="user-defined panel of non-somatic indels in VCF format",
    )
    parser.add_argument(
        "-m",
        "--heap-memory",
        metavar="STR",
        default="6000m",
        help="maximum heap space (defalt: 6000m)",
    )
    parser.add_argument(
        "-l",
        "--log-dir",
        metavar="DIR",
        type=check_folder_existence,
        help="directory for storing log files",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=__version__),
    )
    args = parser.parse_args()
    return args


def create_logger(log_dir):
    logger = logging.getLogger("")
    logger.setLevel(logging.INFO)

    if log_dir:
        fh = logging.FileHandler(os.path.join(log_dir, "rna_indel.log"), delay=True)
        fh.setLevel(logging.INFO)
        fh_formatter = logging.Formatter(
            "%(asctime)s %(module)-12s %(levelname)-8s %(message)s"
        )
        fh.setFormatter(fh_formatter)
        logger.addHandler(fh)

    sh = logging.StreamHandler()
    sh.setLevel(logging.WARNING)
    logger.addHandler(sh)
    return logger


def check_pos_int(val):
    val = int(val)
    if val <= 0:
        sys.exit("Error: The number of processes must be a positive integer.")
    return val


def check_mapq(val):
    val = int(val)
    if not 0 <= val <= 255:
        sys.exit("Error: the MAPQ value must be between 0 and 255.")
    return val


def check_folder_existence(folder):
    p = pathlib.Path(folder)
    if not p.exists():
        sys.exit("Error: {} directory Not Found.".format(folder))
    return folder


def check_file(file_path, file_name):
    if not os.path.isfile(file_path):
        sys.exit("Error: {} Not Found.".format(file_name))
    return file_path


if __name__ == "__main__":
    main()
