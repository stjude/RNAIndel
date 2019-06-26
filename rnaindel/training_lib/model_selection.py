#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
from functools import partial
from multiprocessing import Pool
from sklearn.ensemble import RandomForestClassifier


def input_validator(filename, indel_class):
    """Validate and shuffle data
    Args:
        filename (str): path to input training data
        indel_class (str): "s" for 1-nt, "m" for >1-nt indels
    Returns
        df (pandas.DataFrame)
    """
    df = pd.read_csv(filename, sep="\t")

    if not "truth" in df.columns:
        sys.exit('Error: column label "truth" not found.')

    if indel_class == "s":
        df = df[df["indel_size"] == 1]
    else:
        df = df[df["indel_size"] > 1]

    truth_set = set(df["truth"].values.tolist())
    if truth_set != {"somatic", "germline", "artifact"}:
        which_class = (
            "single-nucleotide indels."
            if indel_class == "s"
            else "multi-nucleotide indels."
        )
        sys.exit('Error: invalid values in column "truth" for ' + which_class)

    # shuffle
    df = df.sample(frac=1, random_state=111).reset_index(drop=True)

    return df


def make_k_folds(df, k, indel_class):
    """Make training and validation sets for k-fold CV
    Args:
        df (pandas.DataFrame)
        k (int): a positive int for the num of folds
        indel_class (str): "s" for 1-nt, "m" for >1-nt indels
    Returns:
        folds (list): a list of [training_df, validation_df]

                      [[df_t_1, df_v_1],
                             ...
                       [df_t_k, df_v_k]]              
    """
    # df with 1-nt or >1-nt
    df = split_by_indel_size(df, indel_class)

    np.random.seed(111)
    df["fold"] = np.random.randint(0, k, df.shape[0])

    folds = []
    for i in range(k):
        df_t = df[df["fold"] != i].reset_index(drop=True)
        df_v = df[df["fold"] == i].reset_index(drop=True)

        df_t.drop("fold", axis=1, inplace=True)
        df_v.drop("fold", axis=1, inplace=True)

        # df_t = preprocess(df_t)
        folds.append([df_t, df_v])

    return folds


def split_by_indel_size(df, indel_class):
    """Split dataframe by indel size
    Args:
        df (pandas.DataFrame)
        indel_class (str): "s" for 1-nt, "m" for >1-nt indels
    Returns:
        df (pandas.DataFrame): only contains 1-nt 
                               or >1-nt indels
    """
    if indel_class == "s":
        return df[df["indel_size"] == 1].reset_index(drop=True)
    else:
        return df[df["indel_size"] > 1].reset_index(drop=True)


def perform_k_fold_cv(folds, features, artifact_ratio, num_of_processes, max_features):
    """Cross Validate in k folds

    Args:
        folds (list): a lisf of [training df, validation df]
        features (list): a list of feature names
        num_of_processes (int): positive int for multiprocessing
    Returns:
        stats (np.array): aggreated from k folds
                           arr([true positive, 
                               false positive,
                               false negative])
    """
    rf = partial(
        run_a_fold,
        features=features,
        artifact_ratio=artifact_ratio,
        max_features=max_features,
    )
    stats = np.zeros(3)
    with Pool(num_of_processes) as p:
        ans = p.map(rf, folds)
    for stat in ans:
        stats = np.add(stats, stat)

    return stats


def run_a_fold(fold, features, artifact_ratio, max_features):
    """Train and validate in a signle fold

    Args:
        fold (list): [training df, validation df]
        features (list): a list of feature names
        artifact_ratio (int): positive int
    Returns:
        stat (np.array): arr([true positive,
                              false positive,
                              false negative])
    """
    df_t, df_v = fold[0], fold[1]
    models = train_classifiers(df_t, artifact_ratio, features, max_features)
    df = predict(df_v, models, features)
    stat = get_stat(df)

    return stat


def train_classifiers(df, artifact_ratio, features, max_features):
    """Train random forest classifiers

    Args:
        df (pandas.DataFrame): training dataframe
        artifact_ratio (int): positive int
        features (list): feature set used for training
    Returns:
        models (list): a list of 20 trained RandomForest obj 
    """
    models = []
    for i in range(20):
        df = downsample(df, artifact_ratio, i)
        X, y = df[features], df["truth"]

        models.append(
            RandomForestClassifier(
                n_estimators=1000, max_features=max_features, random_state=i
            ).fit(X, y)
        )

    return models


def predict(df, models, features):
    """Make prediction

    Args:
        df (pandas.DataFrame): validation dataframe
        models (list): a list of trained classifiers
        features (list): a list of feature names
    Returns:
        df (pandas.DataFrame): dataframe with prediction
    """
    header = ["prob_artifact", "prob_germline", "prob_somatic"]
    prob = partial(predict_prob, df=df, features=features)
    avg_prob = np.average(list(map(prob, models)), axis=0)
    dfp = pd.DataFrame(data=avg_prob, columns=header)
    df = pd.concat([df, dfp], axis=1)
    df["prediction"] = df.apply(annotate_prediction, axis=1)

    return df


def get_stat(df):
    """Calculate performance meterics

    Args:
        df (pandas.DataFrame): dataframe with prediction
    Returns:
        stat (np.array): arr([true positive,
                              false positive,
                              false negative])
    """
    tp = len(df[(df["truth"] == "somatic") & (df["prediction"] == "somatic")])
    fp = len(df[(df["truth"] != "somatic") & (df["prediction"] == "somatic")])
    fn = len(df[(df["truth"] == "somatic") & (df["prediction"] != "somatic")])
    stat = np.array([tp, fp, fn])

    return stat


def downsample(df, artifact_ratio, i):
    """Downsample with specified artifact ratio

    Args:
        df (pandas.DataFrame)
        artifact_ratio (int): positive int x in s:g:a = 1:1:x
        i (int): i-th classifier. used to set seed 
    Return:
        df (pandas.DataFrame): s:g:a = 1:1:x
    """
    df_s = df[df["truth"] == "somatic"]
    df_g = df[df["truth"] == "germline"]
    df_a = df[df["truth"] == "artifact"]

    n = len(df_s)

    # using all somatic indels
    # somatic: germline: artifact = 1:1:x
    df_g = df_g.sample(n=n, random_state=i)
    df_a = df_a.sample(n=n * artifact_ratio, random_state=i)
    df = pd.concat([df_s, df_g, df_a])

    return df


def predict_prob(model, df, features):
    """Calculate prediciton probability

    Args:
        model (RandomForestClassifier): trained classifier
        df (pandas.DataFrame): validation dataframe
        features (list): a list of feature names
    Returns:
        prob (np.array): len(df)x3
    """
    X = df[features]
    prob = model.predict_proba(X)

    return prob


def annotate_prediction(row):
    """Annotate predicated class based 
    on the highest prediction probability
    
    Args:
        row (pandas.Series)
    Returns
        predicted class (str)
    """
    maxp = max(row["prob_artifact"], row["prob_germline"], row["prob_somatic"])

    if maxp == row["prob_artifact"]:
        predicted_class = "artifact"
    elif maxp == row["prob_germline"]:
        predicted_class = "germline"
    else:
        predicted_class = "somatic"

    return predicted_class


def make_score_dict(param_to_be_optimized, stat_arr, beta):
    """Associate performance scores and param being optimized 

    Args:
        param_to_be_optimized (int or str): int for artifact ratio
                                            str for feature subset
        stat_arr (np.array): arr([tp, fp, fn])
        beta (int): int (>=1) for F_beta score
    Returns:
        d (dict)
    """
    f_beta = calculate_performance_score(stat_arr, beta)
    f_label = "tpr" if beta > 100 else "f" + str(beta)

    # precision (beta = 0)
    precision = calculate_performance_score(stat_arr, 0)

    d = {"param": param_to_be_optimized, f_label: f_beta, "precision": precision}

    return d


def calculate_performance_score(stat_arr, beta):
    """Calculate F beta, or precision or TPR

    Args:
        stat_arr (np.array): arr([tp, fp, fn])
        beta (int): non-negative int for F beta score
    Returns:
        Precision (float): for beta = 0
        F_beta (float): for 1 <= beta <= 100
        TPR (float): for > 100
    """
    tp, fp, fn = stat_arr[0], stat_arr[1], stat_arr[2]

    if beta > 100:
        numerator = tp
        denominator = tp + fn
    else:
        numerator = (1 + beta * beta) * tp
        denominator = (1 + beta * beta) * tp + beta * beta * fn + fp
    if denominator == 0:
        return 0.0
    else:
        return numerator / denominator


def report_result(df, beta):
    """Report the best result

    Args:
        df (pandas.DataFrame): df of param and performance score
        beta (int): int (>=1) for F beta score
    Returns:
        best (tuple): (optimized param name, primary score, secondary score)
    """
    f_label = "tpr" if beta > 100 else "f" + str(beta)

    # check if entry with highest f_beta is uniq
    df = df[df[f_label] == df[f_label].max()]

    if len(df) == 1:
        reportable = df

    # check if entry with highest secondary score is uniq
    df = df[df["precision"] == df["precision"].max()]

    if len(df) == 1:
        reportable = df
    else:
        if ("ratio" in df.columns) or ("feature_subset" in df.columns):
            reportable = df.iloc[0]
        else:
            reportable = df.sample(n=1, random_state=123)

    reportable = reportable.reset_index(drop=True)

    labels = reportable.columns.tolist()
    param = [i for i in labels if i not in {f_label, "precision"}][0]

    best = (
        reportable[param].iloc[0],
        reportable[f_label].iloc[0],
        reportable["precision"].iloc[0],
    )

    return best


def features(indel_class):
    """Manage feature names 

    Args:
        indel_class (str): "s" for single-nucleotide indels (1-nt indels)
                           "m" for multi-nucleotide indels (>1-nt indels)
    Returns:
        features (list): a list of feature names
    """
    assert indel_class == "s" or indel_class == "m"

    features_s = [
        "is_ins",
        "is_gc_ins",
        "is_gc_del",
        "is_at_ins",
        "is_at_del",
        "is_splice",
        "is_truncating",
        "is_nmd_insensitive",
        "gc",
        "local_gc",
        "lc",
        "local_lc",
        "strength",
        "local_strength",
        "repeat",
        "indel_complexity",
        "ref_count",
        "alt_count",
        "is_multiallelic",
        "is_near_boundary",
        "is_bidirectional",
        "is_uniq_mapped",
        "cds_length",
        "indel_location",
        "equivalence_exists",
        "ipg",
        "is_on_db",
    ]

    features_m = [
        "is_ins",
        "indel_size",
        "is_inframe",
        "is_truncating",
        "is_splice",
        "is_nmd_insensitive",
        "gc",
        "local_gc",
        "lc",
        "local_lc",
        "strength",
        "local_strength",
        "repeat",
        "dissimilarity",
        "indel_complexity",
        "ref_count",
        "alt_count",
        "is_multiallelic",
        "is_near_boundary",
        "is_bidirectional",
        "is_uniq_mapped",
        "cds_length",
        "indel_location",
        "is_in_cdd",
        "equivalence_exists",
        "ipg",
        "is_on_db",
    ]

    if indel_class == "s":
        return features_s
    elif indel_class == "m":
        return features_m
    else:
        return None
