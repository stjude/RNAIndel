#!/usr/bin/env python3

import sys
import pysam
import logging
import pandas as pd
from functools import partial
from indel_sequence_dev import Indel
from left_aligner import lt_aln
from indel_annotator_dev import indel_annotator 

logger = logging.getLogger(__name__)

def main(df, refgene, fasta, reclf=False):
    fa = pysam.FastaFile(fasta)
    
    # left-alignment
    lt_aln_indel_generator = partial(generate_lt_aln_indel, fa=fa)
    df['lt'] = df.apply(lt_aln_indel_generator, axis=1)
    df['pos'], df['ref'], df['alt'] = zip(*df.apply(left_align_report, axis=1))
     
    # reannotate afer left-alignment
    exon_data = pysam.TabixFile(refgene)
    anno = partial(indel_annotator, exon_data=exon_data, fasta=fasta)
    df['annotation'] = df.apply(anno, axis=1)
    df = df[df['annotation'] != '-']

    if len(df) == 0:
        logging.warning('No indels annotated in coding region\
                         after left-alignment. Analysis done.')
        sys.exit(0)

    df = unify_equivalent_indels(df)
    df = format_header(df, reclf)
    
    return df


def generate_lt_aln_indel(row, fa):
    idl = Indel(row['chr'], row['pos'], row['is_ins'], row['indel_seq'])
    idl = lt_aln(idl, fa)
     
    return idl


def left_align_report(row):
    pos = row['lt'].pos
    if row['is_ins'] == 1:
        ref = '-'
        alt = row['lt'].idl_seq
    else:
        ref = row['lt'].idl_seq
        alt = '-'
       
    return pos, alt, ref 


def unify_equivalent_indels(df):
    # to keep original order
    df['order'] = df.index
   
    # select one with highest somatic probability
    df = df.sort_values('prob_s', ascending=False)
    df = df.drop_duplicates(['chr', 'pos'])
    df = df.sort_values('order')

    return df


def format_header(df, reclf=False):
    header1 = ['chr', 
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

    header2 = ['ref_count',
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
