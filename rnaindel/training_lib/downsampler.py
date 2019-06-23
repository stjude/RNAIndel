#!/usr/bin/env python3

import numpy as np
import pandas as pd
from .model_selection import features
from .model_selection import make_k_folds
from .model_selection import perform_k_fold_cv
from .model_selection import make_score_dict
from .model_selection import report_result


def downsampler(
    df,
    k,
    indel_class,
    beta,
    num_of_processes,
    downsample_ratio,
    max_features="auto",
):
    """Optimize downsampling ratio optimizing F beta: somatic:germline:artifact = 1:1:x

    Args:
        df (pandas.DataFrame)
        k (int): num of folds in cross validation
        indel_class (str): s for single-nucleotide indels, m for multi-nucleotide indels
        beta (int): specify F beta score to be optimized
        num_of_processes (int): num of processes in parallelism
        downsample_ratio (int): specify a downsample ration (default None)
        max_features (str or int): maximum num of features considered in sklearn random forest. default to 'auto'
    Returns:
        report_result (tuple): (artifact_ratio (int), ds_f_beta (float), ds_precision (float)
                               artifact_ratio: ratio optimizing ds_f_beta
                               ds_f_beta: F beta score optimized in downsampling (ds) step. 
                               ds_precision: associated precision at the F beta optimum 
    """

    all_features = features(indel_class)

    folds = make_k_folds(df, k, indel_class)

    # somatic:germline:artifact = 1:1:x for 1 <= x <= 20
    result_dict_lst = []
    
    search_range = (
        range(downsample_ratio, downsample_ratio + 1)
        if downsample_ratio
        else range(1, 21)
    )
    
    for x in search_range:
        # k-fold CrossValidation
        stats = perform_k_fold_cv(
            folds, all_features, x, num_of_processes, max_features
        )
        d = make_score_dict(x, stats, beta)

        result_dict_lst.append(d)

    return report_result(pd.DataFrame(result_dict_lst), beta)
