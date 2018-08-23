#!/usr/bin/env python3

import os
import logging
import argparse
import pandas as pd
import indel_preprocessor as indel_preprocessor
import indel_vcf_processor as indel_vcf_processor
import indel_annotator_dev as indel_annotator_dev
import indel_processor_sequence_dev as indel_processor_sequence_dev
import indel_processor_protein as indel_processor_protein
import indel_equivalence_solver_dev as indel_equivalence_solver_dev
import indel_snp_annotator_dev as indel_snp_annotator_dev 
import indel_classifier as indel_classifier
import indel_reclassifier as indel_reclassifier
import indel_postprocessor as indel_postprocessor
 

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
        df = indel_preprocessor.main(args.bambino_output)
    else:
        df = indel_vcf_processor.main(args.vcf)

    df = indel_annotator_dev.main(df, args.refgene, args.fasta)
    df = indel_processor_sequence_dev.main(df, args.fasta, args.bam)
    df = indel_processor_protein.main(df, args.refgene) 
    df = indel_equivalence_solver_dev.main(df, args.fasta, args.refgene)
    df = indel_snp_annotator_dev.main(df, args.fasta, args.dbsnp, args.clinvar)
    df = indel_classifier.main(df, 
                               args.machine, 
                               processes=args.num_of_processes
                               )
    
    if args.reclassification:
        try:
            if os.path.isfile(args.reclassification_list): 
                df = indel_reclassifier.main(df, args.reclassification_list)
        except:
            df = indel_reclassifier.main(df)
    
    df = indel_postprocessor.main(df, args.refgene, args.fasta, args.reclassification)

    df.to_csv('rna_indels.txt', index=False, sep='\t')


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
    parser.add_argument('-p', '--num-of-processes')
    # optional: reclassification by common SNP
    parser.add_argument('-re-clf', '--reclassification', action='store_true')
    # optional: reclassification by common SNP + user's list
    parser.add_argument('-re-clf-w-lst', '--reclassification-list')
    # configurations: with defaut value
    parser.add_argument('-refgene', '--refgene', default=REFGENE)
    parser.add_argument('-dbsnp', '--dbsnp', default=DBSNP)
    parser.add_argument('-clinvar', '--clinvar', default=CLINVAR)
    parser.add_argument('-machine', '--machine', default=MODELS_DIR)
    args = parser.parse_args()

    return args

if __name__ == '__main__':
    main()
