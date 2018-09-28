#!/usr/bin/env python3

import os
import sys
import pathlib
import logging
import argparse
import pandas as pd
import rna_indel_lib as rna


def main():
    args = get_args()
    
    create_logger(args.output)
    
    if args.input_bambino:
        df = rna.indel_preprocessor(args.input_bambino, args.refgene, args.fasta)
        df = rna.indel_rescuer(df, args.fasta, args.bam, 
                               num_of_processes=args.num_of_processes)
    else:
        df = rna.indel_vcf_preprocessor(args.input_vcf, args.refgene, args.fasta)
        df = rna.indel_rescuer(df, args.fasta, args.bam,
                               num_of_processes=args.num_of_processes,
                               left_aligned=True,
                               external_vcf=True)
              
    df = rna.indel_annotator(df, args.refgene, args.fasta)
    df, df_filtered = rna.indel_sequence_processor(df, args.fasta, args.bam, args.uniq_mapq)
    df = rna.indel_protein_processor(df, args.refgene) 
    df = rna.indel_equivalence_solver(df, args.fasta, args.refgene)
    df = rna.indel_snp_annotator(df, args.fasta, args.dbsnp, args.clinvar)
    df = rna.indel_classifier(df, args.dir_for_models, num_of_processes=args.num_of_processes)
    
    if args.panel_of_non_somatic:
        df = rna.indel_reclassifier(df, args.fasta, args.panel_of_non_somatic)
    
    df = rna.indel_postprocessor(df, args.refgene, args.fasta, args.panel_of_non_somatic)
    rna.indel_vcf_writer(df, df_filtered, args.bam, args.fasta, args.output)
    

def create_logger(output):
    p = pathlib.Path(output)
    out_dir = p.parents[0].as_posix()
 
    logger = logging.getLogger('')
    logger.setLevel(logging.INFO)

    fh = logging.FileHandler(os.path.join(out_dir, 'run.log'), delay=True)
    fh.setLevel(logging.INFO)
    fh_formatter = logging.Formatter('%(asctime)s %(module)-12s %(levelname)-8s %(message)s')
    fh.setFormatter(fh_formatter)
    logger.addHandler(fh)

    sh = logging.StreamHandler()
    sh.setLevel(logging.WARNING)
    logger.addHandler(sh)
    
    return logger


def check_pos_int(val):
    val = int(val)
    if val <=0:
        sys.exit('The number of processes must be a positve integer')

    return val


def check_output(output):
    p = pathlib.Path(output)
    out_dir = p.parents[0]

    if not out_dir.exists():
        sys.exit('Output directory Not Found.')
    
    return output


def check_panel_of_non_somatic(filepath):
    if not os.path.isfile(filepath):
        sys.exit('Panel of non somatic (.vcf) Not Found.')
    
    return filepath


def get_args():
        
    REFGENE = './refgene/refCodingExon.bed.gz'
    DBSNP = './dbsnp/00-All.151.indel.vcf.gz'
    CLINVAR = './clnvr/clinvar.indel.vcf.gz'
    MODELS_DIR = './models'

    parser = argparse.ArgumentParser()
   
    # input required: bam
    parser.add_argument('-b', '--bam', required=True) 
    # input required either: bambino output or vcf
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-ib', '--input-bambino')
    group.add_argument('-iv', '--input-vcf') 
    # output required
    parser.add_argument('-o', '--output', required=True, type=check_output)
    # reference required: fasta
    parser.add_argument('-r', '--fasta', required=True)
    # optional: MAPQ
    parser.add_argument('-q', '--uniq-mapq', default=255, type=int) 
    # optional: number of processes 
    parser.add_argument('-p', '--num-of-processes', type=check_pos_int)
    # optional: reclassification by non-somatic list
    parser.add_argument('-pons', '--panel-of-non-somatic', type=check_panel_of_non_somatic)
    # configurations: with defaut value
    parser.add_argument('-refgene', '--refgene', default=REFGENE)
    parser.add_argument('-dbsnp', '--dbsnp', default=DBSNP)
    parser.add_argument('-clinvar', '--clinvar', default=CLINVAR)
    parser.add_argument('-models', '--dir-for-models', default=MODELS_DIR)
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    main()
