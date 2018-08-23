#!/usr/bin/env python3

import pysam
import pandas as pd
from functools import partial

def main(df, reclf_lst=None):
    df['reclassified'] = df['predicted_class']
    df['comment'] = '-'

    # reclassification by common snp (first priority)
    df['reclassified'],\
    df['comment'] = zip(*df.apply(reclassify_by_common, axis=1))
    
    # OPTIONAL reclassification by bad list (lowest priority)
    if reclf_lst:
        db_reclf = pysam.TabixFile(reclf_lst)
        reclf = partial(relassify_by_list, db_reclf=db_reclf)
        df['reclassified'], df['comment'] = zip(*df.apply(reclf, axis=1))

    return df
     

def reclassify_by_common(row):
    if row['predicted_class'] == 'somatic' and row['is_common'] == 1:
        return 'germline', 'to_germline_by_common_snp'
    else:
        return row['reclassified'], row['comment']


def relassify_by_list(row, db_reclf):
    chr = row['chr']
    pos = row['pos']
    ref = row['ref']
    alt = row['alt']
    prob_a = row['prob_a']
    prob_g = row['prob_g']
    
    if row['reclassified'] == 'somatic' and row['comment'] == '-':
        try:
            to_be_reclf = db_reclf.fetch(chr, pos, pos+1, parser=pysam.asTuple())
        except:
            return row['reclassified'], row['comment']

        if to_be_reclf:
            for instance in to_be_reclf:
                if ref == instance[3] and alt == instance[4]:
                    if prob_a <= prob_g:
                        return 'germline', 'to_gernline_by_reclassification_list'
                    else:
                        return 'artifact', 'to_artifact_by_reclassification_list'
         
    return row['reclassified'], row['comment'] 
