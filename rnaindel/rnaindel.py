#!/usr/bin/env python3

import os
import sys
import pysam
import logging
import pathlib
import argparse
import warnings
import tempfile
import pandas as pd
from functools import partial
from .version import __version__

import rnaindel.bambino_lib as bl
import rnaindel.rnaindel_lib as rl
import rnaindel.training_lib as tl

warnings.filterwarnings("ignore", category=UserWarning)

def main():
    Commands()


class Commands(object):
    def __init__(self):
        parser = argparse.ArgumentParser(
            prog="rnaindel",
            usage="""rnaindel <command> [<args>]

commands are:
    analysis    Predict somatic indels from tumor RNA-Seq data
    feature     Calculate and report features for training
    training    Train models""",
        )

        parser.add_argument("command", help="analysis, feature, or training")
        parser.add_argument(
            "--version",
            action="version",
            version="%(prog)s {version}".format(version=__version__),
        )

        args = parser.parse_args(sys.argv[1:2])

        if not hasattr(self, args.command):
            sys.exit("Error: invalid command")

        getattr(self, args.command)()

    def analysis(self):
        run("analysis")

    def feature(self):
        run("feature")

    def training(self):
        run("training")


def run(command):
    args = get_args(command)
    data_dir = args.data_dir.rstrip("/")
    log_dir = args.log_dir.rstrip("/")
    model_dir = "{}/models".format(data_dir)

    if command == "training":
        df = tl.input_validator(args.training_data)

        # downsampling
        artifact_ratio, ds_f_beta, ds_precision = tl.downsampler(
            df, args.k_fold, args.indel_class, args.ds_beta, args.process_num
        )

        # feature_selection
        selected_features, fs_f_beta, fs_precision = tl.selector(
            df,
            args.k_fold,
            args.indel_class,
            artifact_ratio,
            args.fs_beta,
            args.process_num,
        )

        # parameter tuning
        feature_lst = selected_features.split(";")
        max_features, pt_f_beta, pt_precision = tl.tuner(
            df,
            args.k_fold,
            args.indel_class,
            artifact_ratio,
            feature_lst,
            args.pt_beta,
            args.process_num,
        )

        # update models
        tl.updater(
            df, args.indel_class, artifact_ratio, feature_lst, max_features, model_dir
        )

        # make report
        tl.reporter(
            args.indel_class,
            args.ds_beta,
            ds_f_beta,
            ds_precision,
            artifact_ratio,
            args.fs_beta,
            fs_f_beta,
            fs_precision,
            selected_features,
            args.pt_beta,
            pt_f_beta,
            pt_precision,
            max_features,
            args.log_dir,
        )

        msg = (
            "single-nucleotide indels"
            if args.indel_class == "s"
            else "multi-nucleotide indels"
        )

        print(
            "rnaindel training for " + msg + " completed successfully.", file=sys.stdout
        )

    else:
        create_logger(log_dir)

        genome = pysam.FastaFile(args.fasta)
        alignments = pysam.AlignmentFile(args.bam)
        refgene = "{}/refgene/refCodingExon.bed.gz".format(data_dir)
        exons = pysam.TabixFile(refgene)
        protein = "{}/protein/proteinConservedDomains.txt".format(data_dir)
        dbsnp = pysam.TabixFile("{}/dbsnp/dbsnp.indel.vcf.gz".format(data_dir))
        clinvar = pysam.TabixFile("{}/clinvar/clinvar.indel.vcf.gz".format(data_dir))

        # preprocessing
        # variant calling will be performed if no external VCF is supplied
        if not args.input_vcf:
            # indel calling
            bambino_output = os.path.join(tempfile.mkdtemp(), "bambino.txt")
            bl.bambino(args.bam, args.fasta, bambino_output, args.heap_memory)

            # preprocess indels from the built-in caller
            df, chr_prefixed = rl.indel_preprocessor(
                bambino_output, genome, alignments, exons
            )

            df = rl.indel_rescuer(
                df, args.fasta, args.bam, chr_prefixed, args.process_num
            )

            # delete the temp file
            os.remove(bambino_output)
        else:
            # preprocess indels from external VCF
            df, chr_prefixed = rl.indel_vcf_preprocessor(
                args.input_vcf, genome, alignments, exons
            )

            df = rl.indel_rescuer(
                df,
                args.fasta,
                args.bam,
                chr_prefixed,
                args.process_num,
                external_vcf=True,
            )

        # indel annotation
        df = rl.indel_annotator(df, genome, exons, chr_prefixed)

        # feature calculation
        df, df_filtered_premerge = rl.indel_sequence_processor(
            df, genome, alignments, args.uniq_mapq, chr_prefixed
        )

        df = rl.indel_protein_processor(df, refgene, protein)

        # merging equivalent indels
        df, df_filtered_postmerge = rl.indel_equivalence_solver(
            df, genome, refgene, chr_prefixed
        )

        # dbSNP annotation
        df = rl.indel_snp_annotator(df, genome, dbsnp, clinvar, chr_prefixed)

        # command "feature" exits here
        if command == "feature":
            df = rl.indel_feature_reporter(df, genome, args.output_tab, chr_prefixed)
            print("rnaindel feature completed successfully.", file=sys.stdout)
            sys.exit(0)

        # prediction
        df = rl.indel_classifier(df, model_dir, args.process_num)

        # concatenating invalid(filtered) entries
        df_filtered = pd.concat(
            [df_filtered_premerge, df_filtered_postmerge],
            axis=0,
            ignore_index=True,
            sort=True,
        )

        # panel of non somatic
        pons = (
            os.path.join(args.data_dir, "non_somatic/non_somatic.vcf.gz")
            if not args.non_somatic_panel
            else args.non_somatic_panel
        )
        pons = pysam.TabixFile(pons)

        df = rl.indel_reclassifier(df, genome, pons, chr_prefixed)

        # postProcessing & VCF formatting
        df, df_filtered = rl.indel_postprocessor(
            df, df_filtered, genome, exons, chr_prefixed
        )
        rl.indel_vcf_writer(
            df,
            df_filtered,
            args.fasta,
            genome,
            alignments,
            chr_prefixed,
            args.output_vcf,
            model_dir,
            __version__,
        )

        print("rnaindel analysis completed successfully.", file=sys.stdout)


def get_args(command):
    prog = "rnaindel " + command
    parser = argparse.ArgumentParser(prog=prog)

    if command != "training":
        parser.add_argument(
            "-b",
            "--bam",
            metavar="FILE",
            required=True,
            type=partial(check_file, file_name="BAM file (.bam)"),
            help="input tumor RNA-Seq BAM file (must be STAR-mapped).",
        )

    if command == "training":
        parser.add_argument(
            "-t",
            "--training-data",
            metavar="FILE",
            required=True,
            type=partial(check_file, file_name="data file (.tab)"),
            help="input training data file (tab delimited).",
        )

    if command != "training":
        parser.add_argument(
            "-f",
            "--fasta",
            metavar="FILE",
            required=True,
            type=partial(check_file, file_name="FASTA file"),
            help="reference genome FASTA file.",
        )

    if command != "training":
        parser.add_argument(
            "-d",
            "--data-dir",
            metavar="DIR",
            required=True,
            help="data directory contains databases and models",
            type=check_folder_existence,
        )

    if command == "training":
        parser.add_argument(
            "-d",
            "--data-dir",
            metavar="DIR",
            required=True,
            help="data directory contains databases and models. training will update the models in the directory",
            type=check_folder_existence,
        )

    if command == "analysis":
        parser.add_argument(
            "-o", "--output-vcf", metavar="FILE", required=True, help="output VCF file"
        )
    elif command == "feature":
        parser.add_argument(
            "-o",
            "--output-tab",
            metavar="FILE",
            required=True,
            help="output tab-delimited file",
        )
    else:
        pass

    if command == "training":
        parser.add_argument(
            "-c",
            "--indel-class",
            metavar="STR",
            required=True,
            help="indel class to be trained: s for single-nucleotide indel or m for multi-nucleotide indels",
            type=check_indel_class,
        )

    # input VCF from other callers (optional)
    if command != "training":
        parser.add_argument(
            "-v",
            "--input-vcf",
            metavar="FILE",
            type=partial(check_file, file_name="VCF (.vcf) file"),
            help="input VCF file from other callers",
        )

    if command != "training":
        parser.add_argument(
            "-q",
            "--uniq-mapq",
            metavar="INT",
            default=255,
            type=check_mapq,
            help="STAR mapping quality MAPQ for unique mappers (default: 255)",
        )

    if command == "training":
        parser.add_argument(
            "-k",
            "--k-fold",
            metavar="INT",
            default=5,
            type=check_k_fold,
            help="number of folds in k-fold cross-validation (default: 5)",
        )

    parser.add_argument(
        "-p",
        "--process-num",
        metavar="INT",
        default=1,
        type=check_pos_int,
        help="number of processes (default: 1)",
    )

    if command == "analysis":
        parser.add_argument(
            "-n",
            "--non-somatic-panel",
            metavar="FILE",
            type=partial(check_file, file_name="Panel of non-somatic (.vcf)"),
            help="user-defined panel of non-somatic indels in VCF format",
        )

    if command != "training":
        parser.add_argument(
            "-m",
            "--heap-memory",
            metavar="STR",
            default="6000m",
            help="maximum heap space (default: 6000m)",
        )

    parser.add_argument(
        "-l",
        "--log-dir",
        metavar="DIR",
        default=os.getcwd(),
        type=check_folder_existence,
        help="directory to ouput log files (default: current)",
    )

    if command == "training":
        parser.add_argument(
            "--ds-beta",
            metavar="INT",
            default="10",
            type=check_beta,
            help="F_beta to be optimized in down_sampling step. optimized for TPR when beta >100 given. (default: 10)",
        )

        parser.add_argument(
            "--fs-beta",
            metavar="INT",
            default="10",
            type=check_beta,
            help="F_beta to be optimized in feature selection step. optimized for TPR when beta >100 given. (default: 10)",
        )

        parser.add_argument(
            "--pt-beta",
            metavar="INT",
            default="10",
            type=check_beta,
            help="F_beta to be optimized in parameter_tuning step. optimized for TPR when beta >100 given. (default: 10)",
        )

    args = parser.parse_args(sys.argv[2:])
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
        sys.exit("Error: the number of processes must be a positive integer.")
    return val


def check_k_fold(val):
    val = int(val)
    if val <= 1:
        sys.exit(
            "Error: the number of folds must be a positive integer greater than 1."
        )
    return val


def check_mapq(val):
    val = int(val)
    if not 0 <= val <= 255:
        sys.exit("Error: the MAPQ value must be between 0 and 255.")
    return val


def check_indel_class(val):
    if val != "s" and val != "m":
        sys.exit(
            "Error: indel class must be s for single-nucleotide indels (1-nt) or m for multi-nucleotide indels (>1-nt) indels"
        )
    return val


def check_beta(val):
    val = int(val)
    if val < 1:
        sys.exit("Error: beta must be an interger larger than 0.")
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
