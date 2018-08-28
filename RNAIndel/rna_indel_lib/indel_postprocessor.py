#!/usr/bin/env python3
"""8th (last) step of analysis 

Left-align, unify equivalents and format the result

'indel_postprocessor' is the main routine of this module
"""

import sys
import pysam
import logging
import pandas as pd
from functools import partial
from .indel_sequence_dev import Indel
from .left_aligner import lt_aln
from .indel_annotator_dev import annotate_indels

logger = logging.getLogger(__name__)

def indel_postprocessor(df, refgene, fasta, reclf=False):
    """Main routine to perform left-alingment, unification, and formatting
     
    Args:
        df (pandas.DataFrame): df with prediction made
        refgene (str): path to refCodingExon.bed.gz
        fasta (str): path to .fa
        reclf (bool): True if reclassification is performed. Default=False 
    Returns:
        df (pandas.DataFrame): df with all post-processing done
    """
    fa = pysam.FastaFile(fasta)
    
    # left-alignment
    lt_aln_indel_generator = partial(generate_lt_aln_indel, fa=fa)
    df['lt'] = df.apply(lt_aln_indel_generator, axis=1)
    df['pos'], df['ref'], df['alt'] = zip(*df.apply(left_align_report, axis=1))
     
    # reannotate afer left-alignment
    exon_data = pysam.TabixFile(refgene)
    anno = partial(annotate_indels, exon_data=exon_data, fasta=fasta, postprocess=True)
    df['annotation'] = df.apply(anno, axis=1)
    df = df[df['annotation'] != '-']

    if len(df) == 0:
        logging.warning('No indels annotated in coding region after left-alignment. Analysis done.')
        sys.exit(0)

    df = unify_equivalent_indels(df)
    df = format_header(df, reclf)
    
    return df


def generate_lt_aln_indel(row, fa):
    """Generates a left-aligned Indel object

    Args:
        row (pandas.Series): with 'chr', 'pos', 'is_ins', 'indel_seq'
                            specifies original (not lt-aligned) indel  
        fa (pysam.FastaFile): obj storing reference seq 
    Returns:
        idl (Indel obj): Indel obj left-aligned against reference
    """
    idl = Indel(row['chr'], row['pos'], row['is_ins'], row['indel_seq'])
    idl = lt_aln(idl, fa)
     
    return idl


def left_align_report(row):
    """Gets and formats info from left-aligned indel
    
    Args:
        row (pandas.DataFrame): with 'lt', 'is_ins' labels
                                in 'lt', left-aligned indels are stored
    Returns:
        pos (int): 1-based coordinate
        alt, ref (str): alt or ref allele 
    """
    pos = row['lt'].pos
    if row['is_ins'] == 1:
        ref = '-'
        alt = row['lt'].idl_seq
    else:
        ref = row['lt'].idl_seq
        alt = '-'
       
    return pos, alt, ref 


def unify_equivalent_indels(df):
    """Unify equivalent indels with highest somatic probability

    Args:
        df (pandas.DataFrame): df after left-alignment
                               -> all equivalent indels have same 
                                  'chr', 'pos', 'ref', 'alt'
    Returns:
        df (pandas.Dataframe): de-duplicated dataframe
    """
    # to keep original order
    df['order'] = df.index
   
    # select one with highest somatic probability
    df = df.sort_values('prob_s', ascending=False)
    df = df.drop_duplicates(['chr', 'pos', 'ref', 'alt'])
    df = df.sort_values('order')

    return df


def format_header(df, reclf=False):
    """Formats header for readability
   
    Args:
        df (pandas.DataFrame): df with left-alinged and de-duplicated
        reclf (bool): True if reclassification is performed, Default=False
    Returns:
        df (pandas.DataFrame): formatted df
    """
    header1 = [
               'chr', 
               'pos', 
               'ref', 
               'alt', 
               'annotation',
               'dbsnp', 
               'max_maf',
               'is_common',
               'clin_info',
               'prob_s',
               'prob_g',
               'prob_a',
               'predicted_class'
              ]

    header2 = [
               'ref_count',
               'alt_count',
               'indel_complexity',
               'dissimilarity',
               'indel_size',
               'repeat',
               'is_uniq_mapped',
               'is_near_boundary',
               'is_bidirectional',
               'is_multiallelic',
               'is_truncating',
               'is_nmd_insensitive',
               'ipkc',
               'local_strength',
               'is_at_ins',
               'is_at_del'
               ]
    if reclf:
        header1 = header1 + ['reclassified', 'comment']

    header = header1 + header2

    df = df[header]

    return df
