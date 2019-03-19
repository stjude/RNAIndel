#!/usr/bin/env python3

import os
import gzip
import pickle
from .model_selection import train_classifiers


def updater(df, indel_class, artifact_ratio, features, max_features, model_dir):
    """Update trained models

    Args:
        df (pandas.DataFrame)
        indel_class (str): s for single-nucleotide indels, m for multi-nucleotide indels
        artifact_ratio (int): downsampling ratio for artifact class
        features (list): list of feature names
        max_features (int): maximum num of features considered in sklearn random forest
        model_dir (str): path to dir where trained models are stored
    Returns:
        None
    """
    models = train_classifiers(df, artifact_ratio, features, max_features)

    prefix = "sni." if indel_class == "s" else "mni."

    i = 0
    for model in models:
        path = os.path.join(model_dir, prefix + str(i) + ".pkl.gz")
        model_pkl = gzip.open(path, "wb")
        pickle.dump(model, model_pkl)
        model_pkl.close()
        i += 1

    update_featurefile(indel_class, features, model_dir)


def update_featurefile(indel_class, features, model_dir):
    """Update feature list used in trained model

    Args:
        indel_class (str): s for single-nucleotide indels, m for multi-nucleotide indels
        features (list) list of feature names
        model_dir (str): path to dir where "features.txt" is located
    Returns:
        None
    """
    featurefile = os.path.join(model_dir, "features.txt")
    fr = open(featurefile, "r")
    data = [line.rstrip() for line in fr.readlines()]
    fr.close()

    fw = open(featurefile, "w")
    for line in data:
        if line.startswith(indel_class):
            class_to_be_updated = line.split("\t")[0]
            newline = class_to_be_updated + "\t" + ";".join(features)
            fw.write(newline + "\n")
        else:
            fw.write(line + "\n")
    fw.close()
