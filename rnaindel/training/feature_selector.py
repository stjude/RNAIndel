#!/usr/bin/env python3

import numpy as np
import pandas as pd
from functools import partial
from .model_selection import features
from .model_selection import make_k_folds
from .model_selection import perform_k_fold_cv
from .model_selection import make_score_dict
from .model_selection import report_result


def selector(
    df,
    k,
    indel_class,
    artifact_ratio,
    beta,
    num_of_processes,
    feature_names,
    max_features="auto",
):
    """Select s subset of features optimizing F beta
    
    Args:
        df (pandas,DataFrame)
        k (int): num of folds in cross validation 
        indel_class (str): s for single-nucleotide indels, m for multi-nucleotide indels
        artifact_ratio (int): downsampling ratio for artifact class
        beta (int): specify F beta score to be optimized
        num_of_processes (int): num of processes in parallelism
        feature_names (str): filename specifying a subset of features to be selected
        max_features (str or int): maximum num of features considered in sklearn random forest. default to 'auto'
    Returns:
        report_result (tuple): (selected_features (str), fs_f_beta (float), fs_precision (float))
                              selected_features: subset of features optimizing fs_f_beta. features are semicolon-delimited(;)
                              fs_f_beta: F beta score optimized in feature selection (fs) step
                              fs_precision: associated precision at the F beta optimum
    """
    folds = make_k_folds(df, k, indel_class)

    if feature_names:
        feature_list = [line.rstrip() for line in open(feature_names)]
        selected_features = feature_list[:-1]
        remaining_features = [feature_list[-1]]
    else:
        selected_features = []
        remaining_features = features(indel_class)

    result_dict_lst = []
    while remaining_features:
        d = {}
        best, f_beta, precision = greedy_search(
            selected_features,
            remaining_features,
            folds,
            artifact_ratio,
            beta,
            num_of_processes,
            max_features,
        )

        selected_features.append(best)
        remaining_features.remove(best)

        f_label = "tpr" if beta > 100 else "f" + str(beta)

        d["param"] = ";".join(selected_features)
        d[f_label] = f_beta
        d["precision"] = precision
        result_dict_lst.append(d)

    return report_result(pd.DataFrame(result_dict_lst), beta)


def greedy_search(
    selected_features,
    remaining_features,
    folds,
    artifact_ratio,
    beta,
    num_of_processes,
    max_features,
):
    """Pick up a feature with the greatest increment in F beta in a greedy manner
    
    Args:
        selected_features (list): selected features
        remaining_features (list): features to be examined
        folds (list): folds (list): a k-element list of list [training_df, validation_df]
        artifact_ratio (int): downsampling ratio for artifact class
        beta (int): specify F beta score to be optimized
        num_of_processes (int): num of processes in parallelism
        max_features (str or int): maximum num of features considered in sklearn random forest model
    Returns:
        report_result (tuple): (best, f_beta, precision)
                              best (str): feature name with maximum increment in f_beta
                              f_beta (float): F beta score for model with selected_features + best
                              precision (float): associated precision
    """
    scores = []

    for feature in remaining_features:

        selected_features.append(feature)

        # do k-fold CV for each feature
        stats = perform_k_fold_cv(
            folds, selected_features, artifact_ratio, num_of_processes, max_features
        )
        scores.append(make_score_dict(feature, stats, beta))

        del selected_features[-1]

    return report_result(pd.DataFrame(scores), beta)
