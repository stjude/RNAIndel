#!/usr/bin/env python3 
"""5th step of analysis

Find equivalent indels and merge them.

'indel_equivalence_solver' is the main routine of this module
"""

import re
import numpy as np
import pandas as pd
import subprocess as sp
from .indel_curator_dev import curate_indel_in_genome
from functools import partial
from .indel_protein_processor import acc_len_dict
from .indel_sequence_dev import SequenceWithIndel

mrna = re.compile(r'NM_[0-9]+')

def indel_equivalence_solver(df, fasta, refgene):
    """Solve indel equivalence and calculates 
    indels_per_kilo_coding_sequence (ipkc)
    
    Args:
        df (pandas.DataFrame)
        fasta (str): path to .fa
        refgene (str): path to refCodingExon.bed.gz
    Returns:
        df (pandas.DataFrame)
    """
    # finds and merge equivalent indels
    df = solve_equivalence(df, fasta)
    dfe = df.groupby('equivalence_id')
    df = dfe.apply(merge_equivalents)
    
    # counts indels per transcript in an equivalent aware way
    acc_len = acc_len_dict(refgene)
    dfg = df.groupby('gene_symbol')
    df  = dfg.apply(partial(indels_per_kilo_cds, d=acc_len))
        
    df.drop(['gene_symbol', 'equivalence_id'], axis=1, inplace=True)

    return df


def solve_equivalence(df, fasta):
    """Finds equivalent indels and assigns IDs based on equivalnece
       Equivalent indels are assigned the same ID numbers.

    Args: 
       df (pandas dataframe)
       fasta (str): complete path to FASTA file

    Returns:
       df (pandas datagrame): 'equivalence_id' column added
    """
    # generates indel objects
    df['indel_obj'] = df.apply(partial(generate_indel, fasta=fasta), axis=1)
    
    # finds equivalent indels and outputs an intermediate file
    df['eq'] = df.apply(partial(check_equivalence, df=df), axis=1)
    dfe = df[['eq']].drop_duplicates(['eq'])
    dfe.to_csv('equivalent_indels.txt', sep='\t', index=False, header=False)
    
    # generates dict and assigns IDs 
    d = equivalence_id_dict()
    df['equivalence_id'] = df.apply(partial(assign_id, d=d), axis=1)
    df.drop(['indel_obj', 'eq'], axis=1, inplace=True)
   
    # cleans the intermediate file 
    sp.call(['rm', 'equivalent_indels.txt'])
    
    return df


def generate_indel(row, fasta):
    """Generates indel objects

    Args:
       row (pandas.Series): 4 Series objects respectively
                            column indexed 'chr', 'pos', 
                            'idl_type', 'idl_seq'  
       fasta (str): complete path to FASTA file
    Returns:
       SequenceWithIndel object
    """
    chr = row['chr']
    pos = row['pos']
    idl_type = row['is_ins']
    idl_seq = row['indel_seq']
    return curate_indel_in_genome(fasta, chr, pos, idl_type, idl_seq)


def check_equivalence(row, df):
    """Finds equivalent indels stored in dataframe column
    
    Args:
       row (pandas.Series): A Series object column labeled 'indel_obj'
       df (pandas.DataFrame)
    Returns:
       str: see example below

    Example:
     indel objects are stored in 'indel_obj' column.
     Suppose they are specificed by: 
     (chromosome, position, indel type (ins or del), indel sequence) 
     input  'indel_obj'
              indel 1  (chr1, 100, 1, G)
              indel 2  (chr1, 200, 0, T)
              indel 3  (chr1, 204, 0, T) 
              indel 4  (chr2, 90, 1, ATT)
              indel 5  (chr3, 300, 1, GC)
              ...
              indel n  (chrn, nnn, 1, 'ACTG')
     
     This script finds equivalent indel objects for each object.
     The same indel will be returned as a trivial case
     (for indel i, the same indel i is obviously equivalent).

     finds equivalence (~)
              indel 1 ~ indel 1
              indel 2 ~ indel 2 indel 3
              indel 3 ~ indel 2 indel 3
              indel 4 ~ indel 4
              indel 5 ~ indel 5 indel 9 indel 11
              ...
              indel n ~ indel n
              
     Outputs is a string (chr:pos:indel_type:indel_seq), comma-delimited 
     if multiple equivalent indels exist.
     outout
              chr1:100:1:G
              chr1:200:0:T,chr1:204:0:T
              chr1:200:0:T,chr1:204:0:T
              chr2:90:1:ATT
              chr3:300:1:GC,chr3:302:1:CG,chr3:303:1:GC
              ...
              chrn:nnn:1:ACTG             
    """                  
    idl1 = row['indel_obj']
    chr1 = idl1.chr
     
    # limit search space onto the same chromosome
    search_space = [idl2 for idl2 in df['indel_obj'].values if idl2.chr == chr1]
    
    res = []
    for idl2 in search_space:
        if are_equivalent(idl1, idl2):
            msg = idl2.chr + ':' + str(idl2.pos) + ':'\
                + str(idl2.idl_type) + ':' + idl2.idl_seq
            res.append(msg)
    return ','.join(res)


def equivalence_id_dict():
    """Make dict {chr:pos:indel_type:indel_seq : id}
       the ID number is same for equivalents

    Args:
       None
          The hard-coded filename 'equivalent_indels.txt'
          is an intermediate file generated within this module, 

    Returns:
         dict 
    """
    with open('equivalent_indels.txt') as f:
        d = {}
        i = 1
        for line in f:
            lst = line.rstrip().split(',')
            for key in lst:
                d[key] = i
            i += 1
    return d 


def assign_id(row, d):
    """Assigns equivalence ID to each indel object
    """
    chr = row['chr']
    pos = str(row['pos'])
    idl_type = str(row['is_ins'])
    idl_seq = row['indel_seq']
    key = chr + ':' + pos + ':' + idl_type + ':' + idl_seq
    return d[key]


def merge_equivalents(df):
    """Merges indel counts over equivalent indels.
       
    Args:
       df (pandas.DataFrame): grouped by 'equivalence_id'
    Returns:
       df (pandas.DataFrame): column 'equivalents_exist' added.

    Example:
        indel 1, indel 2 and indel 3 are equivalent
          
                  ref  alt         ref            alt
        indel 1   50     6     ->  33 (50 - 17)   23 (6 + 17)
        indel 2   40     5     ->  22 (40 - 18)   23 (5 + 18)
        indel 3   80    12     ->  69 (80 - 11)   23 (11 + 12)
    """
    pd.options.mode.chained_assignment = None
    
    # if no equivalent indels, return the input
    if len(df) == 1:
        df.loc[:, 'equivalents_exist'] = 0

        return df
    
    # merge all equivalent indel counts
    merged_indel_count = df.loc[:, 'alt_count'].sum()
    
    # adjust ref counts
    diff = merged_indel_count - df['alt_count']
    df.loc[:, 'ref_count'] = df['ref_count'] - diff
    # to ascertain the non-negativity
    df.loc[df['ref_count'] < 0, 'ref_count'] = 0  

    # assign the merged indel count
    df.loc[:, 'alt_count'] = merged_indel_count

    # homoginizes 'is_multiallelic' for equivalents
    if df['is_multiallelic'].sum() > 0:
        df.loc[:, 'is_multiallelic'] = 1
    
    # homoginizes 'is_near_boundary' for equivalents
    if df['is_near_boundary'].sum() > 0:
        df.loc[:, 'is_near_boundary'] = 1

    # homoginizes 'is_bidirectional for equivalents
    if df['is_bidirectional'].sum() > 0:
        df.loc[:, 'is_bidirectional'] = 1
        
    # homoginizes 'is_uniq_mapped' for equivalents
    if df['is_uniq_mapped'].sum() > 0:
        df.loc[:, 'is_uniq_mapped'] = 1
    
    # flags the presence of equivalents
    df.loc[:, 'equivalents_exist'] = 1
     
    return df


def indels_per_kilo_cds(df, d):
    """Counts the number of Indels Per Kilo Coding sequence (IPKC)
    
    Args:
       df (pandas.DataFrame): grouped by 'gene_symbol'
       d (dict): acc_len dict
    Returns:
       df (pandas.DataFrame): 'ipkc' column added
    """
    equivalence_corrected_num_of_indels = len(df['equivalence_id'].unique())
    anno_str = ','.join(df['annotation'].values)
    
    try:
        acc_lst = re.findall(mrna, anno_str)
        ipkc = [equivalence_corrected_num_of_indels*1000/d[acc]\
                for acc in acc_lst]
        df['ipkc'] =  np.median(ipkc)
    except:
        median_cds_len = 1323
        df['ipkc'] = equivalence_corrected_num_of_indels / median_cds_len
    
    return df


def are_equivalent(idl1, idl2):
    """Checks if two indels are equivalent.
       A python implementation of Steve Rive's 
       algorithm.

    Args:
       idl1 (SequenceWithIndel obj)
       idl2 (SequenceWithIndel obj)
    
    Returns:
       bool: True for idl1 and idl2 are equivalent
             False otherwise 
    """
    
    # assume idl2 is on the left side   
    if idl1.pos > idl2.pos:
        idl1, idl2 = idl2, idl1 
    
    chr1 = idl1.chr
    chr2 = idl2.chr
    idl_type1 = idl1.idl_type
    idl_type2 = idl2.idl_type
    idl_seq1 = idl1.idl_seq
    idl_seq2 = idl2.idl_seq
     
    # rejects trivial cases
    if idl1.chr != idl2.chr:
        return False
    if idl_type1 != idl_type2:
        return False
    if len(idl_seq1) != len(idl_seq2):
        return False

    # here after the two indels are 
    # of same type (ins or del) with
    # the same indel length, and 
    # indel 2 pos >= indel 1 pos.

    n = len(idl_seq1)
    m = idl2.pos - idl1.pos
    
    # insertion cases
    if idl_type1 == 1:
        s = idl1.rt_seq[0:m]
        
        if m > n:
            if idl_seq1 == s[:n] and s[:(m-n)] == s[n:] and idl_seq2 == s[-n:]:
                return True
            else:
                return False
        
        elif m == n:
            if idl_seq1 == s == idl_seq2:
                return True
            else:
                return False

        elif m > 0 and m < n:
            if idl_seq1[:m] == s == idl_seq2[-m:] and \
               idl_seq1[-(n-m):] == idl_seq2[:(n-m)]:
                return True
            else:
                return False

        else:
            if idl_seq1 == idl_seq2:
                return True
            else:
                return False

    # deletion cases
    else:
        if m == 0:
            if idl_seq1 == idl_seq2:
                return True
            else:
                return False
        else:
            s = (idl_seq1 + idl1.rt_seq)[:m] + idl_seq2
            if s[:m] == s[n:(m+n)]:
                return True
            else:
                return False
