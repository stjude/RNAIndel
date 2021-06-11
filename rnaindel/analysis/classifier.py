#!/usr/bin/env python3

import os
import gzip
import pickle
import warnings
import numpy as np
import pandas as pd
from functools import partial
from multiprocessing import Pool

warnings.filterwarnings("ignore")

def classify(df, model_dir, num_of_processes):
    """ Makes prediction
    Args:
        df (pandas.DataFrame)
        model_dir (str): path to dir where models are locaded
        num_of_processes (int): the number of processes 
    Returns:
       df (pandas.DataFrame) : with prediction
    """
    df = calculate_proba(df, model_dir, num_of_processes)
    df["predicted_class"] = df.apply(predict_class, axis=1)

    # used in later step
    df["reclassified"] = "-"

    return df


def calculate_proba(df, model_dir, num_of_processes):
    """ Calculates prediction probability for 1-nt (single-nucleotide indels (sni)) 
        and >1-mt (multi-nucleotide indels (mni)) indels
    Args:
        df (pandas.DataFrame): with features calculated 
        model_dir (str): path to dir where model pickle files are located
        num_of_processes (int): a kwarg to specify number of processes for multiprocessing.pool
                                Default = 1
    Returns:
        df (pandas.DataFrame): with prediction probabaility for somatic, germline, artifact
    """
    feature_dict = make_feature_dict(model_dir)

    sni_features = feature_dict["single_nucleotide_indels"]
    mni_features = feature_dict["multi_nucleotide_indels"]
    
    # to keep the original row order
    df["order"] = df.index
    df_sni, df_mni = split_by_indel_size(df)

    pool = Pool(num_of_processes)
    header = ["prob_a", "prob_g", "prob_s"]

    # prediction for 1-nt (sni) indels
    if len(df_sni) > 0:
        sni_models = [
            os.path.join(model_dir, "sni." + str(i) + ".pkl.gz") for i in range(20)
        ]
        sni_pred = partial(predict, data=df_sni, features=sni_features)
        sni_proba = np.average(pool.map(sni_pred, sni_models), axis=0)
        dfp_sni = pd.DataFrame(data=sni_proba)
        dfp_sni.columns = header
    else:
        dfp_sni = pd.DataFrame(columns=header)

    df_sni = pd.concat([df_sni, dfp_sni], axis=1)

    # prediction for >1-nt (mni) indels
    if len(df_mni) > 0:
        mni_models = [
            os.path.join(model_dir, "mni." + str(i) + ".pkl.gz") for i in range(20)
        ]
        mni_pred = partial(predict, data=df_mni, features=mni_features)
        mni_proba = np.average(pool.map(mni_pred, mni_models), axis=0)
        dfp_mni = pd.DataFrame(data=mni_proba)
        dfp_mni.columns = header
    else:
        dfp_mni= pd.DataFrame(columns=header)

    df_mni = pd.concat([df_mni, dfp_mni], axis=1)

    # format output
    df = pd.concat([df_sni, df_mni], axis=0)
    df.sort_values("order", inplace=True)
    df.drop("order", axis=1, inplace=True)
    df.reset_index(drop=True, inplace=True)

    return df


def split_by_indel_size(df):
    """ Sort 1-nt and >1-nt indels
    Args:
        df (pandas.DataFrame)
    Returns:
        df_sni (pandas.DataFrame): df for 1-nt indels
        df_mni (pandas.DataFrame): df for >1-nt indels
    """
    df_sni = df[df["indel_size"] == 1]
    df_mni = df[df["indel_size"] > 1]

    df_sni.reset_index(drop=True, inplace=True)
    df_mni.reset_index(drop=True, inplace=True)

    return df_sni, df_mni


def make_feature_dict(model_dir):
    """Parse features.txt to dict
    Args:
        model_dir (str): path to data directory where "features.txt" is located.
    Returns:
        feature_dict (dict): {indel_class : [feture names]}
    """ 
    features = os.path.join(model_dir, "features.txt")
    f = open(features)
    feature_dict = {line.split("\t")[0]: line.rstrip().split("\t")[1].split(";") for line in f}
    f.close()
    return feature_dict

def predict(model, data, features):
    """ Calculate prediction probabaility
    Args:
        model (file): trained model stored in .pkl.gz 
        data (pandas.DataFrame): df_sni or df_mni
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
