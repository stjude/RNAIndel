#!/usr/bin/env python3
# Somatic indel detector for tumor RNA-Seq data.

import os
import sys
import pathlib
import logging
import argparse
import pandas as pd
from functools import partial

#try:
#    import RNAIndel.rna_indel_lib as ri
#except ImportError:
    # try import rna_indel_lib package directly
import rna_indel_lib as ri


def main():
    args = get_args()
    create_logger(args.log_dir)
    data_dir = args.data_dir.rstrip("/")
    refgene = "{}/refgene/refCodingExon.bed.gz".format(data_dir)
    dbsnp = "{}/dbsnp/00-All.151.indel.vcf.gz".format(data_dir)
    clinvar = "{}/clinvar/clinvar.indel.vcf.gz".format(data_dir)
    model_dir = "{}/models".format(data_dir)
    
    # Preprocessing
    if args.input_bambino:
        df, chr_prefixed = ri.indel_preprocessor(args.input_bambino, args.bam, refgene, args.fasta)
        df = ri.indel_rescuer(
            df, args.fasta, args.bam, chr_prefixed, num_of_processes=args.process_num
        )
        df.to_csv("after_rescue.txt", sep="\t", index=False)
    else:
        df = ri.indel_vcf_preprocessor(args.input_vcf, refgene, args.fasta)
        df = ri.indel_rescuer(
            df,
            args.fasta,
            args.bam,
            num_of_processes=args.process_num,
            left_aligned=True,
            external_vcf=True
        )

    # Analysis 1: indel annotation
    df = ri.indel_annotator(df, refgene, args.fasta)
    
    # Analysis 2: feature calculation using  
    df, df_filtered_premerge = ri.indel_sequence_processor(
        df, args.fasta, args.bam, args.uniq_mapq
    )
    df = ri.indel_protein_processor(df, refgene)
    
    # Analysis 3: mergeing equivalent indels
    df, df_filtered_postmerge = ri.indel_equivalence_solver(
       df, args.fasta, refgene
    )
    
    # Analysis 4: dbSNP annotation
    df = ri.indel_snp_annotator(df, args.fasta, dbsnp, clinvar)
    
    # Analysis 5: prediction
    df = ri.indel_classifier(df, model_dir, num_of_processes=args.process_num)
    
    # Analysis 6: concatenating invalid(filtered) entries
    df_filtered = pd.concat(
        [df_filtered_premerge, df_filtered_postmerge], 
        axis=0, 
        ignore_index=True,
        sort=True
    )
    
    # Analysis 7(Optional): custom refinement of somatic prediction
    if args.non_somatic_panel:
        df = ri.indel_reclassifier(df, args.fasta, args.non_somatic_panel)

    # PostProcessing & VCF formatting
    df, df_filtered = ri.indel_postprocessor(
        df, df_filtered, refgene, args.fasta, args.non_somatic_panel
    )
    ri.indel_vcf_writer(df, df_filtered, args.bam, args.fasta, args.output_vcf)
    
    print("rna_indel completed successfully", file=sys.stderr)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-b",
        "--bam",
        metavar="FILE",
        required=True,
        type=partial(check_file, file_name="BAM file (.bam)"),
        help="input tumor bam file"
    )

    # input indel calls required either: bambino output or a vcf file
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-i",
        "--input-bambino",
        metavar="FILE",
        type=partial(check_file, file_name="Bambino Call Format file"),
        help="input file with calls from Bambino"
    )
    group.add_argument(
        "-c",
        "--input-vcf",
        metavar="FILE",
        type=partial(check_file, file_name="VCF (.vcf) file"),
        help="input vcf file from other callers"
    )
    parser.add_argument(
        "-o", "--output-vcf", metavar="FILE", required=True, help="output vcf file"
    )
    parser.add_argument(
        "-f",
        "--fasta",
        metavar="FILE",
        required=True,
        type=partial(check_file, file_name="FASTA file"),
        help="reference genome (GRCh38) FASTA file. Use the same FASTA file used for mapping"
    )
    parser.add_argument(
        "-d",
        "--data-dir",
        metavar="DIR",
        required=True,
        help="data directory contains refgene, dbsnp and clinvar databases and models",
        type=check_folder_existence
    )
    parser.add_argument(
        "-q",
        "--uniq-mapq",
        metavar="INT",
        default=255,
        type=check_mapq,
        help="STAR mapping quality MAPQ for unique mappers (default: 255)",
    )
    parser.add_argument(
        "-p",
        "--process-num",
        metavar="INT",
        default=1,
        type=check_pos_int,
        help="number of processes (default: 1)",
    )
    parser.add_argument(
        "-n",
        "--non-somatic-panel",
        metavar="FILE",
        type=partial(check_file, file_name="Panel of non-somatic (.vcf)"),
        help="user-defined panel of non-somatic indels in VCF format",
    )
    parser.add_argument(
        "-l",
        "--log-dir",
        metavar="DIR",
        type=check_folder_existence,
        help="directory for storing log files",
    )
    args = parser.parse_args()
    return args


def create_logger(log_dir):
    logger = logging.getLogger("")
    logger.setLevel(logging.INFO)

    if log_dir:
        fh = logging.FileHandler(os.path.join(log_dir, "rna_indel.log"), delay=True)
        fh.setLevel(logging.INFO)
        fh_formatter = logging.Formatter(
            "%(asctime)s %(module)-12s %(levelname)-8s %(message)s"
        )
        fh.setFormatter(fh_formatter)
        logger.addHandler(fh)

    sh = logging.StreamHandler()
    sh.setLevel(logging.WARNING)
    logger.addHandler(sh)
    return logger


def check_pos_int(val):
    val = int(val)
    if val <= 0:
        sys.exit("Error: The number of processes must be a positive integer.")
    return val


def check_mapq(val):
    val = int(val)
    if not 0 <= val <= 255:
        sys.exit("Error: the MAPQ value must be between 0 and 255.")
    return val

def check_folder_existence(folder):
    p = pathlib.Path(folder)
    if not p.exists():
        sys.exit("Error: {} directory Not Found.".format(folder))
    return folder


def check_file(file_path, file_name):
    if not os.path.isfile(file_path):
        sys.exit("Error: {} Not Found.".format(file_name))
    return file_path


if __name__ == "__main__":
    main()
