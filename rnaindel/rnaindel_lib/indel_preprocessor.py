#!/usr/bin/env python3
"""1st step of analysis.

Performs sanity check of Bambino output and extracts 
indel calls on canonical chromosomes.

'indel_preprocessor' is the main routine of this module. 
"""

import os
import sys
import logging
import pandas as pd
from .left_aligner import lt_aln
from .indel_sequence import Indel
from .indel_annotator import generate_coding_indels


logger = logging.getLogger(__name__)


def indel_preprocessor(bambinofile, genome, alignments, exons):
    """ Validate, extract and format indel calls from Bambino output
    Args:
        bambinofile (str): Bambino output filename (contains SNVs + indels)
        genome (pysam.FastaFile): reference genome
        alignments (pysam.AlignmentFile): bam data
        exons (pysam.TabixFile): coding exon data
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
    if not os.path.exists(bambinofile):
        logging.critical("Error: call sets from the built-in caller not found.")
        sys.exit(1)

    df = pd.read_csv(bambinofile, sep="\t", dtype={"dbSNP": str})

    if len(df) == 0:
        logging.critical("Error: no variants called by the built-in caller.")
        sys.exit(1)

    df = extract_necessary_info(df)
    df = extract_indel_calls(df)

    if len(df) == 0:
        logging.warning("No indels detected in variant calling. Analysis done.")
        sys.exit(0)

    df = rename_header(df)

    # check chromosome name format
    chr_prefixed = is_chr_prefixed(alignments)

    # left-alignment
    df = perform_left_alignment(df, genome, chr_prefixed)

    # filter non coding indels
    df["is_coding"] = df.apply(
        flag_coding_indels,
        genome=genome,
        exons=exons,
        chr_prefixed=chr_prefixed,
        axis=1,
    )
    df = df[df["is_coding"]]

    if len(df) == 0:
        logging.warning("No coding indels annotated. Analysis done.")
        sys.exit(0)

    df.drop("is_coding", axis=1, inplace=True)
    df = df.reset_index(drop=True)

    return df, chr_prefixed


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
        bool: True if chr is 1-22, X or Y
    """
    chr = chr.replace("chr", "")

    if chr == "X" or chr == "Y":
        return True
    
    try:
        return (1 <= int(chr) <= 22)
    except:
        return False


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


def is_chr_prefixed(alignments):
    """Check if chromosome names are prefixed with "chr"

    Args:
        alignments (pysam.AlignmeentFile)
    Returns:
        is_prefixed (bool): True if prefixed.
    """
    header_dict = alignments.header
    chromosome_names = header_dict["SQ"]

    is_prefixed = True if chromosome_names[0]["SN"].startswith("chr") else False

    return is_prefixed


def perform_left_alignment(df, genome, chr_prefixed):
    """Perform left alignemnt 
    Args:
        df (pandas.DataFrame)
        genome (pysam.FastaFile): reference genome
        chr_prefixed (bool): True if chromosome names are prefixed with "chr" in BAM
    Returns:
        df (pandas.DataFrame): left aligned
    """
    df = format_indel_report(df)

    df["lt"] = df.apply(generate_lt_aln_indel, genome=genome, chr_prefixed=chr_prefixed, axis=1)

    df["pos"], df["ref"], df["alt"] = zip(*df.apply(left_align_report, axis=1))

    return df


def generate_lt_aln_indel(row, genome, chr_prefixed):
    """Generates a left-aligned Indel object

    Args:
        row (pandas.Series): with 'chr', 'pos', 'is_ins', 'indel_seq'
                            specifies original (not lt-aligned) indel  
        genome (pysam.FastaFile): reference genome
        chr_prefixed (bool): True if chromosome names in BAM are "chr"-prefixed
    Returns:
        idl (Indel obj): Indel obj left-aligned against reference
    """
    idl = Indel(row["chr"], row["pos"], row["is_ins"], row["indel_seq"])
    idl = lt_aln(idl, genome, chr_prefixed)

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
    pos = row["lt"].pos
    if row["is_ins"] == 1:
        ref = "-"
        alt = row["lt"].idl_seq
    else:
        ref = row["lt"].idl_seq
        alt = "-"

    return pos, ref, alt


def format_indel_report(df):
    """Format as follows
    For insertion
           ref          alt               is_ins      indel_seq
           -            inserted_seq      1           inserted_seq
    For deletion
           ref          alt               is_ins      indel_seq
           deleted_seq  -                 0           deleted_seq
    Args:
        df (pandas.DataFrame)
    Returns:
        df (pandas.DataFrame)
    """
    df["ref"] = df.apply(lambda x: "-" if x["ref"] != x["ref"] else x["ref"], axis=1)
    df["alt"] = df.apply(lambda x: "-" if x["alt"] != x["alt"] else x["alt"], axis=1)
    df["is_ins"] = df.apply(lambda x: 1 if x["ref"] == "-" else 0, axis=1)
    df["indel_seq"] = df.apply(lambda x: x["alt"] if x["is_ins"] else x["ref"], axis=1)

    return df


def flag_coding_indels(row, genome, exons, chr_prefixed):
    """Flag indels if they are coding indels
    Args:
        row (panda.Series)
        genome (pysam.FastaFile): reference genome
        exons (pysam.TabixFile): coding exon data
        chr_prefixed (bool): True if chromosome names are prefixed with "chr" in BAM
    Return:
        is_coding (bool): True for coding indels
    """
    res = generate_coding_indels(
        row["chr"], row["pos"], row["is_ins"], row["indel_seq"], genome, exons, chr_prefixed
    )
    is_coding = True if res else False

    return is_coding
