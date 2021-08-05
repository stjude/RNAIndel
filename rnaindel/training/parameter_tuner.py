#!/usr/bin/env python3

import numpy as np
import pandas as pd
from .model_selection import make_k_folds
from .model_selection import perform_k_fold_cv
from .model_selection import make_score_dict
from .model_selection import report_result
from sklearn.model_selection import ParameterGrid


def tuner(
    df, k, indel_class, artifact_ratio, features, beta, num_of_processes, auto_param,
):
    """Optimize the maximum number of features considered in sklearn random forest model

    Args: 
        df (pandas.DataFrame)
        k (int): num of folds in cross validation
        indel_class (str): s for single-nucleotide indels, m for multi-nucleotide indels
        artifact_ratio (int): downsampling ratio for artifact class
        features (list): a list of feature names
        beta (int): specify F beta score to be optimized
        num_of_processes (int): num of processes in parallelism
        auto_param (bool): Train with sklearn's default params if True 
    Returns:
        report_result (tuple): (max_features (int), pt_f_beta (float), pt_precision)
                                max_features: maximum num of features 
                                pt_f_beta: F beta score optimized in parameter tuning (pt) step
                                pt_precision: associated precision
    """
    n_features = len(features)

    folds = make_k_folds(df, k, indel_class)

    # somatic:germline:artifact = 1:1:x for 1 <= x <= 20
    result_dict_lst = []

    search_space = (
        [{"max_features": "auto"}]
        if auto_param
        else prepare_grid(indel_class, n_features)
    )
    for param in search_space:
        # k-fold CrossValidation
        max_features = param["max_features"]
        stats = perform_k_fold_cv(
            folds, features, artifact_ratio, num_of_processes, max_features=max_features
        )
        d = make_score_dict(max_features, stats, beta)

        result_dict_lst.append(d)

    return report_result(pd.DataFrame(result_dict_lst), beta)


def prepare_grid(indel_class, n_features):
    """Make grid for parameter optimization

    Args:
        indel_class (str): s for single-nucleotide indels, m for multi-nucleotide indels
        n_features (int): maximum num of features considered in sklearn random forest
    Returns:
        ParameterGrid (class)
    """
    param_grid = {"max_features": range(1, n_features + 1)}

    return ParameterGrid(param_grid)
