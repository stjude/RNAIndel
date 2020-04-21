import os
import sys
import pysam
import tempfile
import argparse
from functools import partial
from defaultcaller.defaultcaller import callindel
from .preprocessor import preprocess


dbsnp = "/research/rgs01/project_space/zhanggrp/MethodDevelopment/common/mosaicism/test_data_set/dbsnp_v151_grch37.vcf.gz"
#dbsnp = "/home/khagiwar/git/OutLier/data/00-All.38.vcf.gz"

dbsnp = pysam.VariantFile(dbsnp)

def analyze(subcommand):
    
    args = get_args(subcommand) 
 
    with tempfile.TemporaryDirectory() as tmp_dir:
        outfile = os.path.join(tmp_dir, "outfile.txt")
        
        callindel(args.bam, args.reference, outfile, args.heap_memory, args.region) 
        df = preprocess(outfile, pysam.FastaFile(args.reference), pysam.AlignmentFile(args.bam), dbsnp)
    
    
    df = df[["CHROM", "POS", "REF", "ALT", "REF_CNT", "ALT_CNT", "CPOS",  "CREF",  "CALT"]]
    
    df.to_csv(args.output_vcf, sep="\t", index=False)


def get_args(subcommand):
    prog = "rnaindel " + subcommand
    parser = argparse.ArgumentParser(prog=prog)

    parser.add_argument("-i", "--bam", metavar="FILE", required=True, help="input tumor RNA-Seq BAM file (must be STAR-mapped).")

    parser.add_argument("-r", "--reference", metavar="FILE", required=True, help="reference genome FASTA file.")
    
    parser.add_argument("-m", "--heap-memory", metavar="STR", default="6000m", help="maximum heap space (default: 6000m)")

    if subcommand == "analysis":
        parser.add_argument("-o", "--output-vcf", metavar="FILE", required=True, help="output VCF file")

        parser.add_argument("--region", metavar="STR", default=None, help="specify region for target analysis: chrN:start-stop (default: None)")
    else:
        parser.add_argument("-o", "--output-tab", metavar="FILE", required=True, help="output tab-delimited file")
   
    return parser.parse_args(sys.argv[2:])
