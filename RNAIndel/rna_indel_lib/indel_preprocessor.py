#!/usr/bin/env python3
"""1st step of analysis.

Performs sanity check of Bambino output and extracts 
indel calls on canonical chromosomes.

'indel_preprocessor' is the main routine of this module. 
"""

import os
import sys
import pysam
import logging
import pandas as pd
from functools import partial
from .indel_annotator import generate_coding_indels

logger = logging.getLogger(__name__)


def indel_preprocessor(bambinofile, bam, refgene, fasta):
    """ Validate, extract and format indel calls from Bambino output
    Args:
        bambinofile (str): Bambino output filename (contains SNVs + indels)
        refgene (bed file): refCodingExon.bed.gz (contained in data_dir)
        fasta (str): path to reference
    Returns:
        df (pandas.DataFrame): Contains coding indels
                               4 columns formatted as: 

                               chr  pos  ref  alt
                        
                               chr1 123  -    ATC         
                               chr2 456  GG   -
                                      ....
                               chrY 987  CCT  -
        chr_prefixed (bool): True if chromosome names are prefixed with "chr" in BAM
    """
    exon_data = pysam.TabixFile(refgene)
    bam_data = pysam.AlignmentFile(bam)

    if not exists_bambino_output(bambinofile):
        sys.exit(1)

    df = pd.read_csv(bambinofile, sep="\t", dtype={"dbSNP": str})

    if not is_non_trivial_result(df):
        sys.exit(1)

    df = extract_necessary_info(df)
    df = extract_indel_calls(df)

    if len(df) == 0:
        logging.warning("No indels detected in variant calling. Analysis done.")
        sys.exit(0)

    df = rename_header(df)
    df = format_indel_report(df)

    chr_prefixed = is_chr_prefixed(bam_data)
    coding = partial(
        flag_coding_indels, exon_data=exon_data, fasta=fasta, chr_prefixed=chr_prefixed
    )
    df["is_coding"] = df.apply(coding, axis=1)
    df = df[df["is_coding"] == True]

    if len(df) == 0:
        logging.warning("No coding indels annotated. Analysis done.")
        sys.exit(0)

    df.drop("is_coding", axis=1, inplace=True)
    df = df.reset_index(drop=True)

    return df, chr_prefixed


def is_chr_prefixed(bam_data):
    """Check if chromosome names are prefixed with "chr"

    Args:
        bam_data (pysam.AlignmeentFile obj)
    Returns:
        is_prefixed (bool): True if prefixed.
    """
    header_dict = bam_data.header
    chromosome_names = header_dict["SQ"]

    is_prefixed = False
    if chromosome_names[0]["SN"].startswith("chr"):
        is_prefixed = True

    return is_prefixed


def flag_coding_indels(row, exon_data, fasta, chr_prefixed):
    """Flag indels if they are coding indels
    Args:
        row (panda.Series)
        exon_data (pysam.TabixFile obj): coding exon database obj
        fasta (str): path to Fasta
    Return:
        is_coding (bool): True for coding indels
    """
    is_coding = False

    if row["ref"] == "-":
        idl_type, idl_seq = 1, row["alt"]
    else:
        idl_type, idl_seq = 0, row["ref"]

    res = generate_coding_indels(
        row["chr"], row["pos"], idl_type, idl_seq, exon_data, fasta, chr_prefixed
    )
    if res != []:
        is_coding = True

    return is_coding


def exists_bambino_output(filename):
    """Assert if Bambino output file exists
       
    Args:
        filename (str): Bambino output filename
    Returns:
        it_exists (bool): True if exists
    """
    it_exists = False

    if not os.path.exists(filename):
        logging.critical("Error: Bambino output file not found.")
    else:
        it_exists = True
    return it_exists


def is_non_trivial_result(df):
    """Assert if Bambino outputs contains data other than the header line

    Args:
        df (pandas.DataFrame): Bambino output as pd.DataFrame
    Returns:
        is_non_trivial (bool): True if it has more than header line
    """
    is_non_trivial = False

    if len(df) == 0:
        logging.critical("Bambino output only contains the header line.")
    else:
        is_non_trivial = True

    return is_non_trivial


def extract_necessary_info(df):
    """Reduce data size by selecting columns and chromosomes

    Args:
        df (pandas.DataFrame): 45 columns
    Retuns:
        df (pandas.DataFrame): 5 columns with canonical chromosomes
    """
    pd.options.mode.chained_assignment = None

    df = df[["Chr", "Pos", "Type", "Chr_Allele", "Alternative_Allele"]]

    df["is_canonical"] = df.apply(lambda x: is_canonical_chromosome(x["Chr"]), axis=1)
    df = df[df["is_canonical"] == True]
    df.drop("is_canonical", axis=1, inplace=True)

    return df


def is_canonical_chromosome(chr):
    """Check if chr is 1-22, X or Y (M not included)

    Args:
        chr (str): chromosome name
    Returns:
        is_canonical (bool): True if chr is 1-22, X or Y
    """
    is_canonical = False

    if chr.startswith("chr"):
        chr = chr.replace("chr", "")

    if chr == "X" or chr == "Y":
        is_canonical = True
    else:
        try:
            if 1 <= int(chr) <= 22:
                is_canonical = True
        except:
            pass

    return is_canonical


def extract_indel_calls(df):
    """Extract indels calls and sort by position

    Args: 
        df (pandas.DataFrame)
    Returns:
        df (pandas.DataFrame)
    """
    df["original_order"] = df.index

    df_d = df[df["Type"] == "deletion"]
    df_i = df[df["Type"] == "insertion"]
    df = pd.concat([df_d, df_i])

    df.sort_values("original_order", inplace=True)
    df.drop("original_order", axis=1, inplace=True)

    return df


def rename_header(df):
    """Rename as follows
       Chr -> chr
       Pos -> pos
       Chr_Allele -> ref
       Alternative_Allele -> alt

    Args:
        df (pandas.DataFrame)
    Returns:
        df (pandas.DataFrame)
    """
    df.drop("Type", axis=1, inplace=True)
    df = df.rename(
        columns={
            "Chr": "chr",
            "Pos": "pos",
            "Chr_Allele": "ref",
            "Alternative_Allele": "alt",
        }
    )

    return df


def format_indel_report(df):
    """Format as follows
    For insertion
           ref          alt
           -            inserted_seq
    For deletion
           ref          alt
           deleted_seq  -
    Args:
        df (pandas.DataFrame)
    Returns:
        df (pandas.DataFrame)
    """
    df["ref"] = df.apply(lambda x: "-" if x["ref"] != x["ref"] else x["ref"], axis=1)
    df["alt"] = df.apply(lambda x: "-" if x["alt"] != x["alt"] else x["alt"], axis=1)

    return df
