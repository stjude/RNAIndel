#!/usr/bin/env python3
"""5th step of analysis

Find equivalent indels and merge them.

'indel_equivalence_solver' is the main routine of this module
"""

import os
import re
import pathlib
import numpy as np
import pandas as pd
import subprocess as sp
from .indel_curator import curate_indel_in_genome
from .indel_protein_processor import acc_len_dict

mrna = re.compile(r"NM_[0-9]+")


def indel_equivalence_solver(df, genome, refgene, chr_prefixed):
    """Solve indel equivalence and calculates 
    indels_per_gene (ipg)
    
    Args:
        df (pandas.DataFrame)
        genome (pysam.FastaFile): reference genome
        refgene (str): path to refCodingExon.bed.gz
        chr_prefixed (bool): True if chromosome names are "chr"-prefixed
    Returns:
        df (pandas.DataFrame): dataframe with valid entries
        df_filtered_postmerge (pandas.DataFrame):
                               dataframe with entries with alt_count == 1
                               after the merge of equivalen indels
    """
    # finds and merge equivalent indels
    df = solve_equivalence(df, genome, chr_prefixed)
    dfe = df.groupby("equivalence_id")
    df = dfe.apply(merge_equivalents)

    # counts indels per transcript in an equivalent aware way
    acc_len = acc_len_dict(refgene)
    dfg = df.groupby("gene_symbol")
    df = dfg.apply(indels_per_gene, d=acc_len)

    df.drop(["gene_symbol", "equivalence_id"], axis=1, inplace=True)

    df["filtered"] = df.apply(flag_entry_with_one_read, axis=1)

    df, df_filtered_postmerge = df[df["filtered"] == "-"], df[df["filtered"] != "-"]

    return df, df_filtered_postmerge


def flag_entry_with_one_read(row):
    filtered = "lt2count" if row["alt_count"] < 2 else "-"
    return filtered


def solve_equivalence(df, genome, chr_prefixed):
    """Find equivalent indels and assigns IDs based on equivalnece
       Equivalent indels are assigned the same ID numbers.

    Args: 
       df (pandas dataframe)
       genome (pysam.FastaFile): reference genome
       chr_prefixed (bool): True is chromosome names in BAM are "chr"-prefixed
    Returns:
       df (pandas datagrame): 'equivalence_id' column added
    """
    # generate indel objects
    df["indel_obj"] = df.apply(
        generate_indel, genome=genome, chr_prefixed=chr_prefixed, axis=1
    )

    # find equivalent indels and assign id
    df["eq"] = df.apply(check_equivalence, df=df, axis=1)
    data_array = df[["eq"]].drop_duplicates(["eq"]).values
    lst_of_equivalent_indels = [i[0] for i in data_array]
    d = assign_id(lst_of_equivalent_indels)

    # annotate equivalence by id
    df["equivalence_id"] = df.apply(annotate_by_id, d=d, axis=1)
    df.drop(["indel_obj", "eq"], axis=1, inplace=True)

    return df


def generate_indel(row, genome, chr_prefixed):
    """Generate indel objects

    Args:
       row (pandas.Series): 4 Series objects respectively
                            column indexed 'chr', 'pos', 
                            'idl_type', 'idl_seq'  
       genome (pysam.FastaFile): referene genome
       chr_prefixed (bool): True if chromosome names in BAM are "chr"-prefixed
    Returns:
       SequenceWithIndel (class)
    """
    chr = row["chr"]
    pos = row["pos"]
    idl_type = row["is_ins"]
    idl_seq = row["indel_seq"]

    return curate_indel_in_genome(genome, chr, pos, idl_type, idl_seq, chr_prefixed)


def check_equivalence(row, df):
    """Find equivalent indels stored in dataframe column
    
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
              
     Outputs is a string (chr:pos:indel_type:indel_seq), underscore-delimited 
     if multiple equivalent indels exist.
     outout
              chr1:100:1:G
              chr1:200:0:T_chr1:204:0:T
              chr1:200:0:T_chr1:204:0:T
              chr2:90:1:ATT
              chr3:300:1:GC_chr3:302:1:CG_chr3:303:1:GC
              ...
              chrn:nnn:1:ACTG             
    """
    idl1 = row["indel_obj"]
    chr1 = idl1.chr

    # limit search space onto the same chromosome
    search_space = [idl2 for idl2 in df["indel_obj"].values if idl2.chr == chr1]

    res = []
    for idl2 in search_space:
        # equality by equivalence
        if idl1 == idl2:
            msg = (
                idl2.chr
                + ":"
                + str(idl2.pos)
                + ":"
                + str(idl2.idl_type)
                + ":"
                + idl2.idl_seq
            )
            res.append(msg)
    return "_".join(res)


def assign_id(lst_of_equivalent_indels):
    """Assign ID to equivalent indels. 
    Equivalent ones share the same ID.

    Args:
        lst_of_equivalent_indels (list): deduplicated list of 
                                         indels with equivalents grouped 
                                         into the same element
    Returns:
        d (dict): {chr:pos:ins/del:indel_seq: ID(int)}

    Example:
       Equivalent indels are in the same list element delimited by '_'.
       
       lst_of_equivalent_indels = [chr1:100:1:G,
                                   chr1:200:0:T_chr1:204:0:T
                                   chr2:90:1:ATT,
                                   ...]
       Return the dist:

       d = {
            chr1:100:1:G: 1, 
            chr1:200:0:T: 2 (same ID), 
            chr1:204:0:T: 2 (same ID),
            chr2:90:1:ATT: 3,
            ...
           } 
    """
    d = {}
    i = 1
    for idls in lst_of_equivalent_indels:
        for key in idls.split("_"):
            d[key] = i
        i += 1
    return d


def annotate_by_id(row, d):
    """Annotate indels by equivalence id
    
    Args:
        row (pandas.Series)
        d (dict): {chr:pos:ins/del:indel_seq: ID(int)}
    Returns:
        ID (int)
    """
    chr = row["chr"]
    pos = str(row["pos"])
    idl_type = str(row["is_ins"])
    idl_seq = row["indel_seq"]
    key = chr + ":" + pos + ":" + idl_type + ":" + idl_seq
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
        df.loc[:, "equivalence_exists"] = 0

        return df

    # merge all equivalent indel counts
    merged_indel_count = df.loc[:, "alt_count"].sum()

    # adjust ref counts
    diff = merged_indel_count - df["alt_count"]
    df.loc[:, "ref_count"] = df["ref_count"] - diff
    # to ascertain the non-negativity
    df.loc[df["ref_count"] < 0, "ref_count"] = 0

    # assign the merged indel count
    df.loc[:, "alt_count"] = merged_indel_count

    # homoginizes 'is_multiallelic' for equivalents
    if df["is_multiallelic"].sum() > 0:
        df.loc[:, "is_multiallelic"] = 1

    # homoginizes 'is_near_boundary' for equivalents
    if df["is_near_boundary"].sum() > 0:
        df.loc[:, "is_near_boundary"] = 1

    # homoginizes 'is_bidirectional for equivalents
    if df["is_bidirectional"].sum() > 0:
        df.loc[:, "is_bidirectional"] = 1

    # homoginizes 'is_uniq_mapped' for equivalents
    if df["is_uniq_mapped"].sum() > 0:
        df.loc[:, "is_uniq_mapped"] = 1

    # flags the presence of equivalents
    df.loc[:, "equivalence_exists"] = 1

    return df


def indels_per_gene(df, d):
    """Counts the number of indels per gene (ipg)
    
    Args:
       df (pandas.DataFrame): grouped by 'gene_symbol'
       d (dict): acc_len dict
    Returns:
       df (pandas.DataFrame): 'ipg' column added
    """
    equivalence_corrected_num_of_indels = len(df["equivalence_id"].unique())
    anno_str = ",".join(df["annotation"].values)

    try:
        acc_lst = re.findall(mrna, anno_str)
        ipg = [equivalence_corrected_num_of_indels * 1000 / d[acc] for acc in acc_lst]
        df["ipg"] = np.median(ipg)
    except:
        median_cds_len = 1323
        df["ipg"] = equivalence_corrected_num_of_indels / median_cds_len

    return df
