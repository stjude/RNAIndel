#!/usr/bin/env python3

import re
import pysam
import argparse
import operator
import pandas as pd
from functools import partial
from indel_curator_dev import curate_indel_in_genome
from indel_sequence_dev import SequenceWithIndel
from indel_features import IndelReport
from indel_equivalence_solver_dev import are_equivalent


def main(df, fasta, dbsnp, clnvr):
    
    dbsnp = pysam.TabixFile(dbsnp)
    clnvr = pysam.TabixFile(clnvr)

    db_anno = partial(annotate_indel_on_db,\
                      fasta=fasta, dbsnp=dbsnp, clnvr=clnvr)
    df['db'] = df.apply(db_anno, axis=1)
    df['dbsnp'] = df.apply(lambda x: x['db'].report_dbsnp_id(), axis=1)
    df['is_on_dbsnp'] = df.apply(is_on_dbsnp, axis=1)
    df['max_maf'] = df.apply(lambda x: x['db'].report_freq(), axis=1)
    df['is_common'] = df.apply(lambda x: x['db'].is_common(), axis=1)
    #df['is_not_pathogenic'] = df.apply(lambda x: x['db'].is_not_pathogenic(), axis=1)
    #df['with_germline_reports'] = df.apply(lambda x: x['db'].with_germline_reports(), axis=1)
    df['clin_info'] = df.apply(lambda x: x['db'].report_clnvr_info(), axis=1)
    df['is_on_dbsnp'] = df.apply(dbsnp_delister, axis=1)
    
    df.drop('db', axis=1, inplace=True)
    
    return df


def annotate_indel_on_db(row, fasta, dbsnp, clnvr):
    chr = row['chr']
    pos = row['pos']
    idl_type = row['is_ins']
    idl_seq = row['indel_seq']
    
    # obj representing indel to be reported
    this_idl = IndelReport(chr, pos, idl_type, idl_seq)
    # obj representing 'this indel' in reference genome
    idl = curate_indel_in_genome(fasta, chr, pos, idl_type, idl_seq)
    
    # search for equivalent indels over pos +/- search_window nt
    search_window = 50
    start, end = pos - search_window, pos + search_window
    chr_vcf = row['chr'].replace('chr', '')

    for record in dbsnp.fetch(chr_vcf, start, end, parser=pysam.asTuple()):    
        bambinos = vcf2bambino(record)
        for bb in bambinos:
            if idl_type == bb.idl_type and len(idl_seq) == len(bb.idl_seq):
                # indel on db representing in reference genome
                db_idl = curate_indel_in_genome(fasta, chr, bb.pos,\
                                                bb.idl_type, bb.idl_seq)
                if are_equivalent(idl, db_idl):
                    rs = record[2]
                    this_idl.list_dbsnp_id(rs)
                    this_idl.list_dbsnp_freq(dbsnp_freq(record))
                    this_idl.list_dbsnp_origin(dbsnp_origin(record))
                    this_idl.list_dbsnp_common(dbsnp_common(record))

    for record in clnvr.fetch(chr_vcf, start, end, parser=pysam.asTuple()):
        bambinos = vcf2bambino(record)
        for bb in bambinos:
            if idl_type == bb.idl_type and len(idl_seq) == len(bb.idl_seq):
                db_idl = curate_indel_in_genome(fasta, chr, bb.pos,\
                                                bb.idl_type, bb.idl_seq)
                if are_equivalent(idl, db_idl):
                    id = record[2]
                    this_idl.list_clnvr_id(id)
                    this_idl.list_clnvr_freq(clnvr_freq(record))
                    this_idl.list_clnvr_origin(clnvr_origin(record))
                    this_idl.list_clnvr_info(cln_info(record))
    
    return this_idl 


def is_on_dbsnp(row):
    if row['dbsnp'] == '-':
        return 0
    else:
        return 1


def dbsnp_delister(row):
    """Delist dbsnp status if pathogenic
    """
    if 'Pathogenic' in row['clin_info'] or\
        'Likely_pathogenic' in row['clin_info']:
        return 0
    else:
        return row['is_on_dbsnp']


def count_padding_bases(seq1, seq2):
    """
    Args:
       seq1 (str)
       seq2 (str)
    
    Returns:
       int

    Examples:
        a deletion is represented in .vcf as follows
        
        ref   alt
        
        ATATC  AT
        
        By 'left-alignment', ref and alt are alignmed: 
             AT  
             ||      
             ATATC
       
       The first 2 bases are left-aligned.
       In this case, 2 will be returned
    """
    if len(seq1) <= len(seq2):
        shorter, longer  = seq1, seq2
    else:
        shorter, longer  = seq2, seq1
    
    n = 0
    for base1, base2  in zip(shorter, longer[:len(shorter)]):
        if base1 == base2:
            n += 1
        else:
            break

    return n


def vcf2bambino(record):
    """Converts .vcf format to Bambino compatible format
    
    Args:
       record (tuple): vcf line with fields separated in tuple
    Returns:
       list: each element IndelReport obj

    Example:
                 pos 123456789012
           reference ATTAGTAGATGT
           deletion  ATTA---GATGT
           
           del in vcf  N 4 AGTA(ref)  A(alt)
       del in bambino  chr_N 5  GTA(ref)  -(alt)

                 pos 1234***56789012
           reference ATTA***GTAGATGT
           insertion ATTAGTAGTAGATGT
           
          ins in vcf N 4 A(ref) AGTA(alt)
      ins in bambino chr_N 5 -(ref) GTA(alt)    
              
    """ 
    chr = 'chr' + record[0]
    # add 1 to pos (see Example in doc)
    pos = int(record[1]) + 1 
    ref = record[3]
    # alt sequences may be multiple
    alts = record[4].split(',')
    
    parsed = []
    for alt in alts:
        n = count_padding_bases(ref, alt)
        # insertion
        if len(ref) < len(alt):
            idl_type = 1
            idl_seq = alt[n:]
            idl = IndelReport(chr, pos, idl_type, idl_seq)
            parsed.append(idl)
        # deletion
        elif len(ref) > len(alt):
            idl_type = 0
            idl_seq = ref[n:]
            idl = IndelReport(chr, pos, idl_type, idl_seq)
            parsed.append(idl)
        else:
            pass
        
    return parsed


def dbsnp_freq(record):
    """Collects frequency value and common SNP annotation
   
    Args:
       record (tuple): vcf line with fields separated in tuple
    Returns:
       allele_freq (float): the larger value of 1000G and TOPMED
                            -1 if freq info not avail
       is_common (int): 1 if annotated as COMMON
                        0 if not annotated as COMMON
                       -1 if annoation not avail
       origin (int): 0 for unspecified origin
                     1 for germline
                     2 for somatic
                     3 for both
                    -1 for annotation not avail
   """ 
    try:
        kg = re.search(r'(CAF=)([0-9,.e-]+)', record[7]).group(2)
        kg_af = float(kg.split(',')[1])
    except:
        kg_af = -1

    try:
        topmed = re.search(r'(TOPMED=)([0-9,.e-]+)', record[7]).group(2)
        topmed_af = float(topmed.split(',')[1])
    except:
        topmed_af = -1

    return max(kg_af, topmed_af)


def clnvr_freq(record):
    try:
        esp = float(re.search(r'(AF_ESP=)([0-9,.e-]+)', record[7]).group(2))
    except:
        esp = -1

    try:
        exac = float(re.search(r'(AF_EXAC=)([0-9,.e-]+)', record[7]).group(2))
    except:
        exac = -1

    try:
        tgp = float(re.search(r'(AF_TGP=)([0-9,.e-]+)', record[7]).group(2))
    except:
        tgp = -1

    return max(esp, exac, tgp)


def dbsnp_common(record):
    try:
        common = int(re.search(r'(COMMON=)([01])', record[7]).group(2))
    except:
        common = -1
    
    return common


def dbsnp_origin(record):
    try: 
        origin = int(re.search(r'(SAO=)([0123])', record[7]).group(2))
    except:
        origin = -1
    
    return origin


def clnvr_origin(record):
    try:
        origin = int(re.search(r'(ORIGIN=)([0-9]+)', record[7]).group(2))
    except:
        origin = -1

    return origin


def cln_info(record):
    try:
       significance = re.search(r'(CLNSIG=)([A-Za-z_]+)', record[7]).group(2)
    except:
       significance = 'clnSigNA'
    
    try:
        disease = re.search(r'(CLNDN=)([A-Za-z0-9_,-]+)', record[7]).group(2)
    except:
        disease = 'diseaseNA'
    
    return significance + '|' + disease
