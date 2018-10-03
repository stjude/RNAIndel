#!/usr/bin/env python3
"""7th step of analysis

Make prediction for 1-nt (mono) and >1-nt (non-mono) indels

'indel_classifier' is the main routine of this module
"""

import os
import gzip
import pickle
import logging
import numpy as np
import pandas as pd
from functools import partial
from multiprocessing import Pool

logger = logging.getLogger(__name__)


def indel_classifier(df, model_dir, **kwargs):
    """ Makes prediction
    Args:
        df (pandas.DataFrame)
        model_dir (str): path to dir where models are locaded
        processes (int): a positive interger for the number of processes 
    Returns:
       df (pandas.DataFrame) : with prediction
    """
    num_of_processes = kwargs.pop("num_of_processes", 1)

    df = calculate_proba(df, model_dir, num_of_processes)
    df["predicted_class"] = df.apply(predict_class, axis=1)

    # used in later step
    df["reclassified"] = "-"

    return df


def calculate_proba(df, model_dir, num_of_processes):
    """ Calculates prediction probability for 1-nt (mono) and >1-mt (non-mono) indels
    Args:
        df (pandas.DataFrame): with features calculated 
        model_dir (str): path to dir where model pickle files are located
        num_of_processes (int): a kwarg to specify number of processes for multiprocessing.pool
                                Default = 1
    Returns:
        df (pandas.DataFrame): with prediction probabaility for somatic, germline, artifact
    """
    # DO NOT CHANGE THE FEATURE ORDER
    mono_features = [
        "repeat",
        "is_at_del",
        "is_on_dbsnp",
        "alt_count",
        "ref_count",
        "is_at_ins",
        "is_nmd_insensitive",
        "is_near_boundary",
        "indel_complexity",
        "ipg",
        "is_uniq_mapped",
    ]

    non_mono_features = [
        "indel_size",
        "ipg",
        "dissimilarity",
        "alt_count",
        "is_on_dbsnp",
        "ref_count",
        "is_near_boundary",
        "is_truncating",
        "local_strength",
        "indel_complexity",
        "is_uniq_mapped",
        "is_ins",
        "is_multiallelic",
        "is_bidirectional",
    ]

    # to keep the original row order
    df["order"] = df.index
    df_mono, df_non_mono = split_by_indel_size(df)

    pool = Pool(num_of_processes)
    header = ["prob_a", "prob_g", "prob_s"]

    # prediction for mono indels
    if len(df_mono) > 0:
        mono_models = [
            os.path.join(model_dir, "mono." + str(i) + ".pkl.gz") for i in range(20)
        ]
        mono_pred = partial(predict, data=df_mono, features=mono_features)
        mono_proba = np.average(pool.map(mono_pred, mono_models), axis=0)
        dfp_mono = pd.DataFrame(data=mono_proba)
        dfp_mono.columns = header
    else:
        dfp_mono = pd.DataFrame(columns=header)

    df_mono = pd.concat([df_mono, dfp_mono], axis=1)

    # prediction for non mono indels
    if len(df_non_mono) > 0:
        non_mono_models = [
            os.path.join(model_dir, "non_mono." + str(i) + ".pkl.gz") for i in range(20)
        ]
        non_mono_pred = partial(predict, data=df_non_mono, features=non_mono_features)
        non_mono_proba = np.average(pool.map(non_mono_pred, non_mono_models), axis=0)
        dfp_non_mono = pd.DataFrame(data=non_mono_proba)
        dfp_non_mono.columns = header
    else:
        dfp_non_mono = pd.DataFrame(columns=header)

    df_non_mono = pd.concat([df_non_mono, dfp_non_mono], axis=1)

    # format output
    df = pd.concat([df_mono, df_non_mono], axis=0)
    df.sort_values("order", inplace=True)
    df.drop("order", axis=1, inplace=True)
    df.reset_index(drop=True, inplace=True)

    return df


def split_by_indel_size(df):
    """ Sort 1-nt and >1-nt indels
    Args:
        df (pandas.DataFrame)
    Returns:
        df_mono (pandas.DataFrame): df for 1-nt indels
        df_non_mono (pandas.DataFrame): df for >1-nt indels
    """
    df_mono = df[df["indel_size"] == 1]
    df_non_mono = df[df["indel_size"] > 1]

    df_mono.reset_index(drop=True, inplace=True)
    df_non_mono.reset_index(drop=True, inplace=True)

    return df_mono, df_non_mono


def predict(model, data, features):
    """ Calculate prediction probabaility
    Args:
        model (file): trained model stored in .pkl.gz 
        data (pandas.DataFrame): df_mono or df_non_mono
        features (list): a subset of features used for prediction
    Returns:
        prob (tuple): (artifact_prob, germline_prob, somatic_prob) 
    """
    X = data[features]
    model_pkl = gzip.open(model, "rb")
    rf = pickle.load(model_pkl)
    prob = rf.predict_proba(X)
    return prob


def predict_class(row):
    """ Assign class based on the highest probability
    Args:
        row (pandas.Series)
    Returns:
        predicted class (str): 'artifact', 'germline', or 'somatic'
    """
    maxp = max(row["prob_a"], row["prob_g"], row["prob_s"])
    if maxp == row["prob_a"]:
        return "artifact"
    elif maxp == row["prob_g"]:
        return "germline"
    else:
        return "somatic"
