#!/usr/bin/env python3
# Somatic indel detector for tumor RNA-Seq data.

import os
import sys
import pathlib
import logging
import argparse
import RNAIndel


def main():
    args = get_args()
    create_logger(args.output)
    
    if args.input_bambino:
        df = RNAIndel.indel_preprocessor(args.input_bambino, args.refgene, args.fasta)
        df = RNAIndel.indel_rescue(df, args.fasta, args.bam, 
                                   num_of_processes=args.num_of_processes)
    else:
        df = RNAIndel.indel_vcf_preprocessor(args.input_vcf, args.refgene, args.fasta)
        df = RNAIndel.indel_rescuer(df, args.fasta, args.bam,
                                    num_of_processes=args.num_of_processes,
                                    left_aligned=True, external_vcf=True)
              
    df = RNAIndel.indel_annotator(df, args.refgene, args.fasta)
    df = RNAIndel.indel_sequence_processor(df, args.fasta, args.bam, args.uniq_mapq)
    df = RNAIndel.indel_protein_processor(df, args.refgene) 
    df = RNAIndel.indel_equivalence_solver(df, args.fasta, args.refgene, args.output)
    df = RNAIndel.indel_snp_annotator(df, args.fasta, args.dbsnp, args.clinvar)
    df = RNAIndel.indel_classifier(df, args.dir_for_models, num_of_processes=args.num_of_processes)
    
    if args.panel_of_non_somatic:
        df = RNAIndel.indel_reclassifier(df, args.fasta, args.panel_of_non_somatic)
    
    df = RNAIndel.indel_postprocessor(df, args.refgene, args.fasta, args.panel_of_non_somatic)
    RNAIndel.indel_vcf_writer(df, args.bam, args.fasta, args.output)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--input-bam', metavar='FILE', required=True, help='input tumor bam file')

    # input indel calls required either: bambino output or a vcf file
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-i', '--input-bambino', metavar='FILE', help='input file with Bambino indel calls')
    group.add_argument('-c', '--input-vcf', metavar='FILE', help='input vcf file with indel calls')

    parser.add_argument('-o', '--output-vcf', metavar='FILE', required=True, type=check_output,
                        help='output vcf file')
    parser.add_argument('-f', '--fasta', metavar='FILE', required=True,
                        help='reference genome (GRCh38) FASTA file')
    parser.add_argument('-r', '--refgene', metavar='FILE', default=REFGENE, required=True,
                        help='refgene coding exon database')
    parser.add_argument('-d', '--dbsnp', metavar='FILE', default=DBSNP, required=True,
                        help='indels on dbSNP database in vcf format')
    parser.add_argument('-l', '--clinvar', metavar='FILE', default=CLINVAR, required=True,
                        help='ClinVar database')
    parser.add_argument('-m', '--model-dir', metavar='DIR', default=MODELS_DIR,
                        help='directory with trained random forest models')
    parser.add_argument('-q', '--uniq-mapq', metavar='INT', default=255, type=int,
                        help='STAR mapping quality MAPQ for unique mappers')
    parser.add_argument('-p', '--process-num', metavar='INT', default=1, type=check_pos_int,
                        help='number of processes (default is 1)')
    parser.add_argument('-n', '--non-somatic-panel', metavar='FILE', type=check_panel_of_non_somatic,
                        help='user-defined panel of non-somatic indel list in vcf format')
    args = parser.parse_args()
    return args


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
    if val <= 0:
        sys.exit('The number of processes must be a positve integer')
    return val


def check_output(output):
    p = pathlib.Path(output)
    out_dir = p.parents[0]
    if not out_dir.exists():
        sys.exit('Output directory Not Found.')
    return output


def check_panel_of_non_somatic(file_path):
    if not os.path.isfile(file_path):
        sys.exit('Panel of non somatic (.vcf) Not Found.')
    return file_path


if __name__ == '__main__':
    main()
