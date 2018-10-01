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
    data_dir = args.data_dir.rstrip('/')
    create_logger(data_dir)
    refgene = '{}/refgene/refCodingExon.bed.gz'.format(data_dir)
    dbsnp = '{}/dbsnp/00-All.151.indel.vcf.gz'.format(data_dir)
    clinvar = '{}/clinvar/clinvar.indel.vcf.gz'.format(data_dir)

    if args.input_bambino:
        df = RNAIndel.indel_preprocessor(args.input_bambino, refgene, args.fasta)
        df = RNAIndel.indel_rescue(df, args.fasta, args.bam, num_of_processes=args.num_of_processes)
    else:
        df = RNAIndel.indel_vcf_preprocessor(args.input_vcf, args.refgene, args.fasta)
        df = RNAIndel.indel_rescuer(df, args.fasta, args.bam, num_of_processes=args.num_of_processes,
                                    left_aligned=True, external_vcf=True)
              
    df = RNAIndel.indel_annotator(df, refgene, args.fasta)
    df = RNAIndel.indel_sequence_processor(df, args.fasta, args.bam, args.uniq_mapq)
    df = RNAIndel.indel_protein_processor(df, refgene)
    df = RNAIndel.indel_equivalence_solver(df, args.fasta, refgene, args.output)
    df = RNAIndel.indel_snp_annotator(df, args.fasta, dbsnp, clinvar)
    df = RNAIndel.indel_classifier(df, args.dir_for_models, num_of_processes=args.num_of_processes)
    
    if args.panel_of_non_somatic:
        df = RNAIndel.indel_reclassifier(df, args.fasta, args.panel_of_non_somatic)
    
    df = RNAIndel.indel_postprocessor(df, refgene, args.fasta, args.panel_of_non_somatic)
    RNAIndel.indel_vcf_writer(df, args.bam, args.fasta, args.output)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--input-bam', metavar='FILE', required=True, help='input tumor bam file')

    # input indel calls required either: bambino output or a vcf file
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-i', '--input-bambino', metavar='FILE', help='input file with Bambino indel calls')
    group.add_argument('-c', '--input-vcf', metavar='FILE', help='input vcf file with indel calls')

    parser.add_argument('-o', '--output-vcf', metavar='FILE', required=True, type=check_folder_existence,
                        help='output vcf file')
    parser.add_argument('-f', '--fasta', metavar='FILE', required=True,
                        help='reference genome (GRCh38) FASTA file')
    parser.add_argument('-d', '--data-dir', metavar='DIR', required=True, type=check_folder_existence,
                        help='data directory contains refgene, dbsnp and clivar databases')
    parser.add_argument('-q', '--uniq-mapq', metavar='INT', default=255, type=int,
                        help='STAR mapping quality MAPQ for unique mappers')
    parser.add_argument('-p', '--process-num', metavar='INT', default=1, type=check_pos_int,
                        help='number of processes (default is 1)')
    parser.add_argument('-n', '--non-somatic-panel', metavar='FILE', type=check_panel_of_non_somatic,
                        help='user-defined panel of non-somatic indel list in vcf format')
    parser.add_argument('-l', '--log-dir', metavar='DIR', type=check_folder_existence,
                        help='directory for storing log files.')
    args = parser.parse_args()
    return args


def create_logger(log_dir):
    logger = logging.getLogger('')
    logger.setLevel(logging.INFO)

    print('log_dir: ' + log_dir)
    if log_dir:
        fh = logging.FileHandler(os.path.join(log_dir, 'rna_indel.log'), delay=True)
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
        sys.exit('Error: The number of processes must be a positive integer.')
    return val


def check_folder_existence(folder):
    p = pathlib.Path(dir)
    out_dir = p.parents[0]
    if not out_dir.exists():
        sys.exit('Error: {} directory Not Found.'.format(folder))
    return folder


def check_panel_of_non_somatic(file_path):
    if not os.path.isfile(file_path):
        sys.exit('Error: Panel of non somatic (.vcf) Not Found.')
    return file_path


if __name__ == '__main__':
    main()
