#!/usr/bin/env python3

import numpy as np
import pandas as pd
from functools import partial
from .model_selection import features
from .model_selection import make_k_folds
from .model_selection import perform_k_fold_cv
from .model_selection import make_score_dict
from .model_selection import report_result


def selector(df, k, indel_class, artifact_ratio, beta, num_of_processes, max_features="auto"):
    remaining_features = features(indel_class)

    folds = make_k_folds(df, k, indel_class)
    
    selected_features = []
    
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
### Remove later !!!!!!!!!!!!!!!!!!!!!!!
    df = pd.DataFrame(result_dict_lst)
    df.to_csv("feature_selection.m.txt", sep="\t", index=False)

    return report_result(pd.DataFrame(result_dict_lst), beta)


def greedy_search(
    selected_features, remaining_features, folds, artifact_ratio, beta, num_of_processes, max_features):

    scores = []
    for feature in remaining_features:

        selected_features.append(feature)

        # do k-fold CV for each feature
        stats = perform_k_fold_cv(
            folds, selected_features, artifact_ratio, num_of_processes, max_features
        )
        scores.append(make_score_dict(feature, stats, beta))

        del selected_features[-1]
        print(feature, stats)
    return report_result(pd.DataFrame(scores), beta)
