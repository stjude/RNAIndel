#!/usr/bin/env python3

import os
import sys
import gzip
import pickle
import logging
import numpy as np
import pandas as pd
from functools import partial
from multiprocessing import Pool
from sklearn.ensemble import RandomForestClassifier

logger = logging.getLogger(__name__)

def main(df, machine_dir, processes):
    
    # requires at least 2 indel fragments
    df = df[df['alt_count'] > 1]      
    
    df = calculate_proba(df, machine_dir)
    df['predicted_class'] = df.apply(predict_class, axis=1)

    return df


def calculate_proba(df, machine_dir, **kwargs):
    mono_features= ['repeat', 'is_at_del', 'is_on_dbsnp', 'alt_count',\
                    'ref_count', 'is_at_ins', 'is_nmd_insensitive', 'is_near_boundary',\
                    'indel_complexity', 'ipkc', 'is_uniq_mapped']

    non_mono_features = ['indel_size', 'ipkc', 'dissimilarity', 'alt_count',\
                         'is_on_dbsnp', 'ref_count', 'is_near_boundary',\
                         'is_truncating', 'local_strength', 'indel_complexity',\
                         'is_uniq_mapped', 'is_ins', 'is_multiallelic', 'is_bidirectional']
    
    # to keep the original row order
    df['order'] = df.index
    df_mono, df_non_mono = split_by_indel_size(df)
    
    pool = Pool(kwargs.pop('num_of_processes', 1))
    header = ['prob_a', 'prob_g', 'prob_s']
    
    # prediction for mono indels
    if len(df_mono) > 0:
        mono_machines\
        = [os.path.join(machine_dir, 'mono.'+str(i)+'.pkl.gz') for i in range(20)]  
        mono_pred = partial(predict, data=df_mono, features=mono_features)
        mono_proba = np.average(pool.map(mono_pred, mono_machines), axis=0)
        dfp_mono = pd.DataFrame(data=mono_proba)
        dfp_mono.columns = header
    else:
        dfp_mono = pd.DataFrame(columns=header)
     
    df_mono = pd.concat([df_mono, dfp_mono], axis=1)
    
    # prediction for non mono indels
    if len(df_non_mono) > 0 :
        non_mono_machines\
        = [os.path.join(machine_dir, 'non_mono.'+str(i)+'.pkl.gz') for i in range(20)]
        non_mono_pred = partial(predict, data=df_non_mono, features=non_mono_features)
        non_mono_proba = np.average(pool.map(non_mono_pred, non_mono_machines), axis=0)
        dfp_non_mono = pd.DataFrame(data=non_mono_proba)
        dfp_non_mono.columns = header
    else:
        dfp_non_mono = pd.DataFrame(columns=header)
    
    df_non_mono = pd.concat([df_non_mono, dfp_non_mono], axis=1)
    
    # format output
    df = pd.concat([df_mono, df_non_mono], axis=0)
    df.sort_values('order', inplace=True)
    df.drop('order', axis=1, inplace=True)
   
    return df

    
def split_by_indel_size(df):
    df_mono = df[df['indel_size'] == 1]
    df_non_mono = df[df['indel_size'] > 1]
   
    df_mono.reset_index(drop=True, inplace=True)
    df_non_mono.reset_index(drop=True, inplace=True)
     
    return df_mono, df_non_mono


def predict(machine, data, features):
    X = data[features]
    model_pkl = gzip.open(machine, 'rb')
    model = pickle.load(model_pkl)
    prob = model.predict_proba(X)
    return prob


def predict_class(row):
    maxp = max(row['prob_a'], row['prob_g'], row['prob_s'])
    
    if maxp == row['prob_a']:
        return 'artifact'
    elif maxp == row['prob_g']:
        return 'germline'
    else:
        return 'somatic'

