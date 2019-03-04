#!/usr/bin/env python3

import numpy as np
import pandas as pd
from .model_selection import make_k_folds
from .model_selection import perform_k_fold_cv
from .model_selection import make_score_dict
from .model_selection import report_result
from sklearn.model_selection import ParameterGrid


def tuner(df, k, indel_class, artifact_ratio, features, pt_beta, num_of_processes):
  
    n_features = len(features)
    
    folds = make_k_folds(df, k, indel_class)
    
    # somatic:germline:artifact = 1:1:x for 1 <= x <= 20
    result_dict_lst = []
    for param in prepare_grid(indel_class, n_features):
        # k-fold CrossValidation
        max_features = param["max_features"]
        stats = perform_k_fold_cv(folds, features, artifact_ratio, num_of_processes, max_features=max_features)
        d = make_score_dict(max_features, stats, pt_beta)
        
        result_dict_lst.append(d)
     
#### Remove later!!!!!!!!!!!!!!!!!!!!!!!!!!!
    df = pd.DataFrame(result_dict_lst)
    df.to_csv("test_param.m.txt", sep="\t", index=False)
    return report_result(pd.DataFrame(result_dict_lst))

def prepare_grid(indel_class, n_features):
    param_grid ={"max_features":range(1, n_features+1)}
    
    return ParameterGrid(param_grid)
    






if __name__ == "__main__":
    file = "all.txt"
    #file = "whole_training.txt"
    df = pd.read_csv(file, sep="\t")
    features = tl.features("m")
    tunuer(df, 5, "m", 3, features, 10, 0, 12)
