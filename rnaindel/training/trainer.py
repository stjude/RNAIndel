import os
import sys
import tempfile
import argparse
from functools import partial

from .model_selection import input_validator
from .downsampler import downsampler
from .feature_selector import selector
from .parameter_tuner import tuner
from .model_updater import updater
from .result_reporter import reporter
from .homopolyer_trainer import train_homolopolymer


def train():

    subcommand = "Train"
    args = get_args(subcommand)

    indel_class = args.indel_class
    data_dir = args.data_dir.rstrip("/")

    df = input_validator(args.training_data, indel_class)

    if indel_class in ["s", "m"]:

        # downsampling
        artifact_ratio, ds_f_beta, ds_precision = downsampler(
            df,
            args.k_fold,
            args.indel_class,
            args.ds_beta,
            args.process_num,
            args.downsample_ratio,
        )

        # feature_selection
        selected_features, fs_f_beta, fs_precision = selector(
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
        max_features, pt_f_beta, pt_precision = tuner(
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
        model_dir = "{}/models".format(data_dir)
        updater(
            df, args.indel_class, artifact_ratio, feature_lst, max_features, model_dir
        )

        # make report
        reporter(
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
        model_dir = "{}/outliers".format(data_dir)
        train_homolopolymer(df, model_dir)


def get_args(subcommand):
    prog = "rnaindel " + subcommand
    parser = argparse.ArgumentParser(prog=prog)

    parser.add_argument(
        "-t",
        "--training-data",
        metavar="FILE",
        required=True,
        type=validate_file_input,
        help="input training data file (tab delimited).",
    )

    parser.add_argument(
        "-c",
        "--indel-class",
        metavar="STR",
        required=True,
        type=validate_indel_class,
        help="indel class to be trained: s for single-nucleotide indels, m for multi-nucleotide indels, h for homopolymer indels",
    )

    parser.add_argument(
        "-d",
        "--data-dir",
        metavar="DIR",
        required=True,
        type=validate_dir_input,
        help="data directory contains databases and models. training will update the models in the directory",
    )

    parser.add_argument(
        "-k",
        "--k-fold",
        metavar="INT",
        default=5,
        help="number of folds in k-fold cross-validation (default: 5)",
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
        "--ds-beta",
        metavar="INT",
        default="10",
        type=validate_int_inputs,
        help="F_beta to be optimized in down_sampling step. optimized for TPR when beta >100 given. (default: 10)",
    )

    parser.add_argument(
        "--fs-beta",
        metavar="INT",
        default="10",
        type=validate_int_inputs,
        help="F_beta to be optimized in feature selection step. optimized for TPR when beta >100 given. (default: 10)",
    )

    parser.add_argument(
        "--pt-beta",
        metavar="INT",
        default="10",
        type=validate_int_inputs,
        help="F_beta to be optimized in parameter_tuning step. optimized for TPR when beta >100 given. (default: 10)",
    )

    parser.add_argument(
        "--downsample-ratio",
        metavar="INT",
        default=None,
        type=validate_int_inputs,
        help="train with specified downsample ratio in [1, 20]. (default: None)",
    )

    parser.add_argument(
        "--feature-names",
        metavar="FILE",
        default=None,
        type=validate_file_input,
        help="train with specified subset of features. Supply as file containing a feature name per line (default: None)",
    )

    parser.add_argument(
        "--auto-param",
        action="store_true",
        help='train with sklearn.RandomForestClassifer\'s max_features="auto"',
    )

    parser.add_argument(
        "-l",
        "--log-dir",
        metavar="DIR",
        default=os.getcwd(),
        type=validate_dir_input,
        help="directory to ouput training results (default: current)",
    )

    return parser.parse_args(sys.argv[2:])


def validate_file_input(file_path):
    if os.path.isfile(file_path):
        return file_path
    else:
        sys.exit("Error: {} not found.".format(file_path))


def validate_dir_input(dir_path):
    if os.path.isdir(dir_path):
        return dir_path
    else:
        sys.exit("Error: {} not found.".format(dir_path))


def validate_indel_class(val):
    if val not in ["s", "m", "h"]:
        sys.exit(
            "Error: indel class must be s for single-nucleotide indels (1-nt) or m for multi-nucleotide indels (>1-nt) indels or h for homopolyer"
        )

    return val


def validate_int_inputs(val):
    val = int(val)
    if val <= 0:
        sys.exit("Error: the input must be a positive integer.")

    return val
