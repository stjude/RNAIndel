#!/usr/bin/env python3

import os
import re
import sys
import pysam
import uuid
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
import rnaindel.nonsomatic_lib as nl


warnings.filterwarnings("ignore", category=UserWarning)


def main():
    Commands()


class Commands(object):
    def __init__(self):
        parser = argparse.ArgumentParser(
            prog="rnaindel",
            usage="""rnaindel <subcommand> [<args>]

subcommands are:
    analysis              Predict somatic indels from tumor RNA-Seq data
    feature               Calculate and report features for training
    nonsomatic            Compile non-somatic indel panel
    reclassification      Reclassify false positives by non-somatic panel
    recurrence            Annotate false positives by recurrence
    training              Train models""",
        )

        parser.add_argument(
            "subcommand",
            help="analysis, feature, nonsomatic, reclassification, recurrence, training",
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

    def analysis(self):
        run("analysis")

    def feature(self):
        run("feature")

    def training(self):
        run("training")

    def nonsomatic(self):
        run("nonsomatic")

    def reclassification(self):
        run("reclassification")

    def recurrence(self):
        run("recurrence")


def run(subcommand):
    args = get_args(subcommand)

    if subcommand == "reclassification":
        nl.filterate_by_panel(
            args.input_vcf,
            args.output_vcf,
            pysam.FastaFile(args.fasta),
            args.non_somatic_panel,
        )
        print("rnaindel reclassification completed successfully.", file=sys.stdout)
        sys.exit(0)

    data_dir = args.data_dir.rstrip("/")
    model_dir = "{}/models".format(data_dir)
    # database check
    path2cosmic = pathlib.Path("{}/cosmic".format(data_dir))

    if not path2cosmic.exists():
        print(
            "Please download the latest database: http://ftp.stjude.org/pub/software/RNAIndel/"
        )
        sys.exit(1)

    if subcommand == "nonsomatic" or subcommand == "recurrence":
        cosmic = pysam.TabixFile(
            "{}/cosmic/CosmicCodingMuts.indel.vcf.gz".format(data_dir)
        )
        if subcommand == "nonsomatic":
            nl.make_non_somatic_panel(
                args.vcf_list,
                args.output_vcf,
                pysam.FastaFile(args.fasta),
                cosmic,
                args.count,
            )
            print("rnaindel nonsomaic completed successfully.", file=sys.stdout)
            sys.exit(0)
        else:
            nl.annotate_recurrence(
                args.vcf_list, pysam.FastaFile(args.fasta), cosmic, args.out_dir
            )
            print("rnaindel recurrence completed successfully.", file=sys.stdout)
            sys.exit(0)

    log_dir = args.log_dir.rstrip("/")

    if subcommand == "training":
        df = tl.input_validator(args.training_data, args.indel_class)

        # downsampling
        artifact_ratio, ds_f_beta, ds_precision = tl.downsampler(
            df,
            args.k_fold,
            args.indel_class,
            args.ds_beta,
            args.process_num,
            args.downsample_ratio,
        )

        # feature_selection
        selected_features, fs_f_beta, fs_precision = tl.selector(
            df,
            args.k_fold,
            args.indel_class,
            artifact_ratio,
            args.fs_beta,
            args.process_num,
            args.feature_names,
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
            args.auto_param,
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

        alignments = pysam.AlignmentFile(args.bam)
        genome = pysam.FastaFile(args.fasta)
        refgene = "{}/refgene/refCodingExon.bed.gz".format(data_dir)
        exons = pysam.TabixFile(refgene)
        protein = "{}/protein/proteinConservedDomains.txt".format(data_dir)
        dbsnp = pysam.TabixFile("{}/dbsnp/dbsnp.indel.vcf.gz".format(data_dir))
        clinvar = pysam.TabixFile("{}/clinvar/clinvar.indel.vcf.gz".format(data_dir))
        cosmic = pysam.TabixFile(
            "{}/cosmic/CosmicCodingMuts.indel.vcf.gz".format(data_dir)
        )

        germline_db = pysam.TabixFile(args.germline_db) if args.germline_db else None

        # input validation
        rl.input_validator(alignments, genome, args.uniq_mapq)

        # region analysis
        region = args.region if subcommand == "analysis" else None

        # preprocessing
        # variant calling will be performed if no external VCF is supplied
        if not args.input_vcf:

            with tempfile.TemporaryDirectory() as tmp_dir:
                # indel calling
                bambino_output = os.path.join(tmp_dir, "bambino.txt")

                bl.bambino(
                    args.bam, args.fasta, bambino_output, args.heap_memory, region
                )

                # preprocess indels from the built-in caller
                df, chr_prefixed = rl.indel_preprocessor(
                    bambino_output, genome, alignments, exons
                )

                df = rl.indel_rescuer(
                    df, args.fasta, args.bam, chr_prefixed, args.process_num
                )

        else:
            # preprocess indels from external VCF
            df, chr_prefixed = rl.indel_vcf_preprocessor(
                args.input_vcf, genome, alignments, exons, region
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
        if subcommand == "feature":
            df, df_filtered_premerge = rl.indel_sequence_processor(
                df,
                genome,
                alignments,
                args.uniq_mapq,
                chr_prefixed,
                softclip_analysis=args.softclip_analysis,
            )
        else:
            coverage_in_trainingset = "{}/models/coverage.txt".format(data_dir)
            downsample_thresholds = {}
            with open(coverage_in_trainingset) as f:
                for line in f:
                    if line.startswith("s"):
                        downsample_thresholds["single_nuleotide_indels"] = int(
                            line.rstrip().split("\t")[1]
                        )
                    else:
                        downsample_thresholds["multi_nuleotide_indels"] = int(
                            line.rstrip().split("\t")[1]
                        )

            df, df_filtered_premerge = rl.indel_sequence_processor(
                df,
                genome,
                alignments,
                args.uniq_mapq,
                chr_prefixed,
                softclip_analysis=args.softclip_analysis,
                downsample_thresholds=downsample_thresholds,
            )

        df = rl.indel_protein_processor(df, refgene, protein)

        # merging equivalent indels
        df, df_filtered_postmerge = rl.indel_equivalence_solver(
            df, genome, refgene, chr_prefixed
        )

        # SNP annotation
        df = rl.indel_snp_annotator(
            df, genome, dbsnp, clinvar, cosmic, germline_db, chr_prefixed
        )

        # subcommand "feature" exits here
        if subcommand == "feature":
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
        default_pons = pysam.TabixFile(
            os.path.join(args.data_dir, "non_somatic/non_somatic.vcf.gz")
        )
        user_pons = (
            pysam.TabixFile(args.non_somatic_panel) if args.non_somatic_panel else None
        )

        df = rl.indel_reclassifier(
            df, genome, default_pons, user_pons, cosmic, chr_prefixed
        )

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


def get_args(subcommand):
    prog = "rnaindel " + subcommand
    parser = argparse.ArgumentParser(prog=prog)

    if subcommand == "analysis" or subcommand == "feature" or subcommand == "training":
        parser.add_argument(
            "-p",
            "--process-num",
            metavar="INT",
            default=1,
            type=check_int,
            help="number of processes (default: 1)",
        )

        parser.add_argument(
            "-l",
            "--log-dir",
            metavar="DIR",
            default=os.getcwd(),
            type=check_folder_existence,
            help="directory to ouput log files (default: current)",
        )

        if subcommand == "analysis" or subcommand == "feature":
            parser.add_argument(
                "-i",
                "--bam",
                metavar="FILE",
                required=True,
                type=partial(check_file, file_name="BAM file (.bam)"),
                help="input tumor RNA-Seq BAM file (must be STAR-mapped).",
            )

            parser.add_argument(
                "-r",
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
                help="data directory contains databases and models",
                type=check_folder_existence,
            )

            parser.add_argument(
                "-v",
                "--input-vcf",
                metavar="FILE",
                type=partial(check_file, file_name="VCF (.vcf) file"),
                help="input VCF file from other callers",
            )

            parser.add_argument(
                "-q",
                "--uniq-mapq",
                metavar="INT",
                default=255,
                type=partial(check_int, preset="mapq"),
                help="STAR mapping quality MAPQ for unique mappers (default: 255)",
            )

            parser.add_argument(
                "-m",
                "--heap-memory",
                metavar="STR",
                default="6000m",
                help="maximum heap space (default: 6000m)",
            )

            parser.add_argument(
                "-g",
                "--germline-db",
                metavar="FILE",
                type=check_file,
                help="user-provided germline database in bgzip-compressed VCF format",
            )

            parser.add_argument(
                "--exclude-softclipped-alignments",
                dest="softclip_analysis",
                action="store_false",
                help="do not use softclipped alignments for feature calculation",
            )
            parser.set_defaults(softclip_analysis=True)

            if subcommand == "analysis":

                parser.add_argument(
                    "-o",
                    "--output-vcf",
                    metavar="FILE",
                    required=True,
                    help="output VCF file",
                )

                parser.add_argument(
                    "-n",
                    "--non-somatic-panel",
                    metavar="FILE",
                    default=None,
                    type=partial(check_file, file_name="nonsomatic panel (.vcf.gz)"),
                    help="user-defined panel of non-somatic indels in bgzip-compressed VCF format (default: None)",
                )

                parser.add_argument(
                    "--region",
                    metavar="STR",
                    default=None,
                    type=check_region,
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

        else:
            parser.add_argument(
                "-t",
                "--training-data",
                metavar="FILE",
                required=True,
                type=partial(check_file, file_name="data file (.tab)"),
                help="input training data file (tab delimited).",
            )

            parser.add_argument(
                "-d",
                "--data-dir",
                metavar="DIR",
                required=True,
                help="data directory contains databases and models. training will update the models in the directory",
                type=check_folder_existence,
            )

            parser.add_argument(
                "-c",
                "--indel-class",
                metavar="STR",
                required=True,
                help="indel class to be trained: s for single-nucleotide indel or m for multi-nucleotide indels",
                type=check_indel_class,
            )

            parser.add_argument(
                "-k",
                "--k-fold",
                metavar="INT",
                default=5,
                type=partial(check_int, preset="k_fold"),
                help="number of folds in k-fold cross-validation (default: 5)",
            )

            parser.add_argument(
                "--ds-beta",
                metavar="INT",
                default="10",
                type=check_int,
                help="F_beta to be optimized in down_sampling step. optimized for TPR when beta >100 given. (default: 10)",
            )

            parser.add_argument(
                "--fs-beta",
                metavar="INT",
                default="10",
                type=check_int,
                help="F_beta to be optimized in feature selection step. optimized for TPR when beta >100 given. (default: 10)",
            )

            parser.add_argument(
                "--pt-beta",
                metavar="INT",
                default="10",
                type=check_int,
                help="F_beta to be optimized in parameter_tuning step. optimized for TPR when beta >100 given. (default: 10)",
            )

            parser.add_argument(
                "--downsample-ratio",
                metavar="INT",
                default=None,
                type=partial(check_int, preset="downsample"),
                help="train with specified downsample ratio in [1, 20]. (default: None)",
            )

            parser.add_argument(
                "--feature-names",
                metavar="FILE",
                default=None,
                type=check_file,
                help="train with specified subset of features. Supply as file containing a feature name per line (default: None)",
            )

            parser.add_argument(
                "--auto-param",
                action="store_true",
                help='train with sklearn.RandomForestClassifer\'s max_features="auto"',
            )

    elif (
        subcommand == "nonsomatic"
        or subcommand == "reclassification"
        or subcommand == "recurrence"
    ):
        parser.add_argument(
            "-r",
            "--fasta",
            metavar="FILE",
            required=True,
            type=partial(check_file, file_name="FASTA file"),
            help="reference genome FASTA file.",
        )

        if subcommand == "nonsomatic" or subcommand == "recurrence":
            parser.add_argument(
                "-d",
                "--data-dir",
                metavar="DIR",
                required=True,
                help="data directory contains databases and models",
                type=check_folder_existence,
            )
            if subcommand == "nonsomatic":
                parser.add_argument(
                    "--vcf-list",
                    metavar="FILE",
                    required=True,
                    type=check_file,
                    help="file containing paths to normal VCF files",
                )

                parser.add_argument(
                    "--count",
                    metavar="INT",
                    required=True,
                    type=check_int,
                    help="use indels observed >= count times for panel compilation",
                )

                parser.add_argument(
                    "-o",
                    "--output-vcf",
                    metavar="FILE",
                    required=True,
                    help="nonsomatic VCF file",
                )
            else:
                parser.add_argument(
                    "--vcf-list",
                    metavar="FILE",
                    required=True,
                    type=check_file,
                    help="file containing paths to RNAIndel output VCF files to be annotated",
                )

                parser.add_argument(
                    "--out-dir",
                    metavar="DIR",
                    default=None,
                    type=check_folder_existence,
                    help="directory to ouput (default: input file directory)",
                )
        else:
            parser.add_argument(
                "-i",
                "--input-vcf",
                metavar="FILE",
                required=True,
                type=partial(check_file, file_name="VCF (.vcf) file"),
                help="RNAIndel ouput VCF to be reclassified",
            )

            parser.add_argument(
                "-o",
                "--output-vcf",
                metavar="FILE",
                required=True,
                help="reclassified VCF file",
            )

            parser.add_argument(
                "-n",
                "--non-somatic-panel",
                metavar="FILE",
                required=True,
                type=partial(check_file, file_name="nonsomatic panel (.vcf.gz)"),
                help="user-defined panel of non-somatic indels in bgzip-compressed VCF format",
            )

    args = parser.parse_args(sys.argv[2:])
    return args


def create_logger(log_dir):
    logger = logging.getLogger("")
    logger.setLevel(logging.INFO)

    if log_dir:
        logfilename = "rnaindel_" + str(uuid.uuid4()) + ".log"
        fh = logging.FileHandler(os.path.join(log_dir, logfilename), delay=True)
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


def check_int(val, preset=None):
    val = int(val)
    if val <= 0 and not preset:
        sys.exit("Error: the input must be a positive integer.")
    elif val <= 1 and preset == "k_fold":
        sys.exit(
            "Error: the number of folds must be a positive integer greater than 1."
        )
    elif not 0 <= val <= 255 and prest == "mapq":
        sys.exit("Error: the MAPQ value must be an integer between 0 and 255.")
    elif not 1 <= val <= 20 and preset == "downsample":
        sys.exit("Error: downsample ratio must be an integer between 1 and 20")
    else:
        return val


def check_indel_class(val):
    if val != "s" and val != "m":
        sys.exit(
            "Error: indel class must be s for single-nucleotide indels (1-nt) or m for multi-nucleotide indels (>1-nt) indels"
        )
    return val


def check_folder_existence(folder):
    p = pathlib.Path(folder)
    if not p.exists():
        sys.exit("Error: {} directory Not Found.".format(folder))
    return folder


def check_file(file_path, file_name=None):
    if not os.path.isfile(file_path):
        if file_name:
            sys.exit("Error: {} Not Found.".format(file_name))
        else:
            sys.exit("Error: {} Not Found.".format(file_path))
    return file_path


def check_region(region):
    region = region.replace("chr", "")
    is_valid = bool(re.match(r"[0-9XY]+:[0-9]+-[0-9]+", region))

    if is_valid:
        roi_lst = region.split(":")
        chr, span = roi_lst[0], roi_lst[1].split("-")
        start, stop = int(span[0]), int(span[1])
    else:
        sys.exit("Check the region format: chrN:start-stop")

    canonicals = [str(i) for i in range(1, 23)] + ["X", "Y"]
    if not chr in canonicals:
        sys.exit("RNAIndel only supports human canonical chrmosomes: chr1-22,X,Y")

    if start >= stop:
        sys.exit(
            "Check the stop is larger than the start in your region: chrN:start-stop"
        )

    return ("chr{}".format(chr), start, stop)


if __name__ == "__main__":
    main()
