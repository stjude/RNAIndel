#!/usr/bin/env python3
"""Optional step for reclassification.

Knowledge-based reclassification is performed
for indels predicted somatic.

'indel_reclassifier' is the main routine in this module
"""

import pysam
import pandas as pd
from functools import partial

def indel_reclassifier(df, reclf_bed_file=None):
    """Main module to reclassify based on common SNP 
    and, if user-defined list is provided, reclassify based on
    the list.

    Args:
        df (pandas.DataFrame): df with prediction made
        reclf_bed_file (file): user-defined .bed. Tabix-indexed. 
                               Default=None (not provided)
    Returns:
        df (pandas.DataFrame): df reclassified
    """
    df['reclassified'] = df['predicted_class']
    df['comment'] = '-'

    # reclassification by common snp (first priority)
    df['reclassified'], df['comment'] = zip(*df.apply(reclassify_by_common, axis=1))
    
    # OPTIONAL reclassification by bad list (lowest priority)
    if reclf_bed_file:
        db_reclf = pysam.TabixFile(reclf_bed_file)
        reclf = partial(relassify_by_list, db_reclf=db_reclf)
        df['reclassified'], df['comment'] = zip(*df.apply(reclf, axis=1))

    return df
     

def reclassify_by_common(row):
    """Reclassifies common SNP indels predicted somatic to germline
   
    Args:
        row (pandas.Series)
    Returns:
        row['reclassified'] (str): reclassifed class
        row['comment'] (str): 'to_germline_by_common_snp' if reclassified 
    """
    if row['predicted_class'] == 'somatic' and row['is_common'] == 1:
        return 'germline', 'to_germline_by_common_snp'
    else:
        return row['reclassified'], row['comment']


def relassify_by_list(row, db_reclf):
    """Reclassifies indels predicted somatic using user-defined reclassification list

    Args:
        row (pandas.Series)
        db_reclf (pysam.TabixFile obj): object storing indexed .bed file 
    Returns:
        row['reclassified'] (str): reclassifed class
        row['comment'] (str): 'to_***_by_reclassification_list' if reclassified
    """
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
                        return 'germline', 'to_germline_by_reclassification_list'
                    else:
                        return 'artifact', 'to_artifact_by_reclassification_list'
         
    return row['reclassified'], row['comment'] 
