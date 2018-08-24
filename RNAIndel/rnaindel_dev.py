#!/usr/bin/env python3

import os
import logging
import argparse
import pandas as pd
import rna_indel_lib as rna

logger = logging.getLogger('')
logger.setLevel(logging.INFO)

fh = logging.FileHandler('./run.log', delay=True)
fh.setLevel(logging.INFO)
fh_formatter = logging.Formatter('%(asctime)s %(module)-12s %(levelname)-8s %(message)s')
fh.setFormatter(fh_formatter)
logger.addHandler(fh)

sh = logging.StreamHandler()
sh.setLevel(logging.WARNING)
logger.addHandler(sh)


def main():
    args = get_args()

    if args.bambino_output:
        df = rna.indel_preprocessor(args.bambino_output)
    else:
        df = rna.indel_vcf_processor(args.vcf)

    df = rna.indel_annotator(df, args.refgene, args.fasta)
    df = rna.indel_sequence_processor(df, args.fasta, args.bam)
    df = rna.indel_protein_processor(df, args.refgene) 
    df = rna.indel_equivalence_solver(df, args.fasta, args.refgene)
    df = rna.indel_snp_annotator(df, args.fasta, args.dbsnp, args.clinvar)
    df = rna.indel_classifier(df, args.dir_for_models, processes=args.num_of_processes)
    
    if args.reclassification:
        try:
            if os.path.isfile(args.reclassification_list): 
                df = rna.reclassify_indels(df, args.reclassification_list)
        except:
            df = rna.indel_reclassifier(df)
    
    df = rna.indel_postprocessor(df, args.refgene, args.fasta, args.reclassification)

    df.to_csv('rna_indels.txt', index=False, sep='\t')


def check_pos_int(val):
    val = int(val)
    if val <=0:
        raise argparse.ArgumentTypeError('The number of processes must be a positve integer')

    return val


def get_args():
        
    REFGENE = './refgene/refCodingExon.bed.gz'
    DBSNP = './dbsnp/00-All.151.indel.vcf.gz'
    CLINVAR = './clnvr/clinvar.indel.vcf.gz'
    MODELS_DIR = './models'

    parser = argparse.ArgumentParser()
   
    # required: bam
    parser.add_argument('-b', '--bam', required=True) 
    # required either: bambino output or vcf
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-f', '--bambino-output')
    group.add_argument('-v', '--vcf') 
    # required: fasta
    parser.add_argument('-r', '--fasta', required=True)
    # optional: number of processes 
    parser.add_argument('-p', '--num-of-processes', type=check_pos_int)
    # optional: reclassification by common SNP
    parser.add_argument('-re-clf', '--reclassification', action='store_true')
    # optional: reclassification by common SNP + user's list
    parser.add_argument('-re-clf-w-lst', '--reclassification-list')
    # configurations: with defaut value
    parser.add_argument('-refgene', '--refgene', default=REFGENE)
    parser.add_argument('-dbsnp', '--dbsnp', default=DBSNP)
    parser.add_argument('-clinvar', '--clinvar', default=CLINVAR)
    parser.add_argument('-models', '--dir-for-models', default=MODELS_DIR)
    args = parser.parse_args()

    return args

if __name__ == '__main__':
    main()
