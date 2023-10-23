import os
import sys
import argparse
from .model_updater import updater
from .model_selection import input_validator


def setup():
    """Serialize pretrained models with user's sklearn version"""

    subcommand = "SetUp"
    args = get_args(subcommand)
    data_dir = args.data_dir.rstrip("/")
    model_dir = "{}/models".format(data_dir)
    training_data = "{}/trainingset/training_set.txt.gz".format(data_dir)

    for indel_class in ["s", "m"]:
        df = input_validator(training_data, indel_class)
        params = pretrained_params(indel_class)

        updater(df, indel_class, params[0], params[1], params[2], model_dir)


def pretrained_params(indel_class):
    """Hard-coded training params used in the pretrained models
    input:
     indel_class (str): "s" for 1-nt indels, "m" for longer indels
    """
    if indel_class == "s":
        artifact_ratio = 17
        features = [
            "repeat",
            "is_on_db",
            "alt_count",
            "ref_count",
            "indel_location",
            "is_at_ins",
            "is_near_boundary",
            "is_multiallelic",
            "indel_complexity",
            "local_lc",
            "strength",
            "is_uniq_mapped",
        ]
        max_features = 3
    else:
        artifact_ratio = 2
        features = [
            "dissimilarity",
            "cds_length",
            "local_gc",
            "is_on_db",
            "alt_count",
            "ref_count",
            "local_lc",
            "is_truncating",
            "strength",
            "is_bidirectional",
            "lc",
            "is_nmd_insensitive",
            "is_uniq_mapped",
            "gc",
            "is_splice",
            "is_near_boundary",
            "is_multiallelic",
            "indel_size",
            "is_ins",
            "repeat",
            "indel_complexity",
            "is_in_cdd",
            "indel_location",
        ]
        max_features = 4

    return artifact_ratio, features, max_features


def get_args(subcommand):
    prog = "rnaindel " + subcommand
    parser = argparse.ArgumentParser(prog=prog)

    parser.add_argument(
        "-d",
        "--data-dir",
        metavar="DIR",
        required=True,
        type=validate_dir_input,
        help="data directory contains databases and models. training will update the models in the directory",
    )

    return parser.parse_args(sys.argv[2:])


def validate_dir_input(dir_path):
    if os.path.isdir(dir_path):
        return dir_path
    else:
        sys.exit("Error: {} not found.".format(dir_path))
