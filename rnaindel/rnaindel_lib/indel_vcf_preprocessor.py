#!/usr/bin/env python3
"""Optional step for .vcf input

Convert .vcf to a format compatible with Bambino

'indel_vcf_preprocessor' is the main routine of this module
"""

import sys
import logging
import pandas as pd
from .indel_preprocessor import flag_coding_indels
from .indel_preprocessor import is_chr_prefixed
from .indel_preprocessor import is_canonical_chromosome
from .indel_preprocessor import perform_left_alignment

logger = logging.getLogger(__name__)


def indel_vcf_preprocessor(vcffile, genome, alignments, exons):
    """Convert input VCF to Bambino format and check chromosome name format
    
    Args:
        vcffile (str): path to input vcf
        genome (pysam.FastaFile): reference genome
        alignments (pysam.AlignmentFile): bam data
        exons (pysam.TabixFile): coding exon data
    Returns:
        df (pandas.DataFrame): df with indels reported as in Bambino output
    """
    vcf_data = open(vcffile)
    df = pd.DataFrame(make_data_list(vcf_data))
    vcf_data.close()

    if len(df) == 0:
        logging.warning("No indels detected in input vcf. Analysis done.")
        sys.exit(0)

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


def make_data_list(vcf_data):
    """
    Args:
        vcf_data (file): VCF content 
    Returns:
        data (list): with dict element 
    """
    data = []
    for line in vcf_data:
        parsed_line_lst = parse_vcf_line(line)
        if parsed_line_lst:
            for parsed_line in parsed_line_lst:
                d = {
                    "chr": parsed_line[0],
                    "pos": parsed_line[1],
                    "ref": parsed_line[2],
                    "alt": parsed_line[3],
                }

                data.append(d)

    return data


def parse_vcf_line(line):
    """Parse VCF to Bambino format

    Args:
        line (str): a line in VCF file.
    Returns:
        parsed_line_lst (lst): with tuple elem (chr, pos, ref, alt)
    Example:
      deletion
              pos 123456789012
        reference ATTAGTAGATGT
        deletion  ATTA---GATGT
         
        VCF:
            CHROM POS REF  ALT
            N     4   AGTA A
        Bambino:
            chr   pos ref alt
            chr_N 5   GTA -
         
      insertion
              pos 1234***56789012
        reference ATTA***GTAGATGT
        insertion ATTAGTAGTAGATGT
         
        VCF:
            CHROM POS REF ALT
            N     4   A   AGTA
        Bambino:
            chr   pos ref alt
            chr_N 5   -   GTA   
    """
    parsed_line_lst = []

    # skip header lines
    if line.startswith("#"):
        return parsed_line_lst

    lst = line.rstrip().split("\t")
    chr = lst[0]
    vcf_pos = int(lst[1])
    vcf_ref = lst[3]
    vcf_alts = lst[4].split(",")  # possibly multi-allelic

    if not chr.startswith("chr"):
        chr = "chr" + chr

    # skip non canonical chrmosomes
    if not is_canonical_chromosome(chr):
        return parsed_line_lst

    for vcf_alt in vcf_alts:
        trimmed_ref, trimmed_alt = right_trim(vcf_ref, vcf_alt)
        n = count_padding_bases(trimmed_ref, trimmed_alt)
        pos = vcf_pos + n

        if len(trimmed_ref) < len(trimmed_alt) and len(trimmed_alt) <= 300:
            ref = "-"
            alt = trimmed_alt[n:]
            parsed_line_lst.append((chr, pos, ref, alt)) 
        elif len(trimmed_ref) > len(trimmed_alt) and len(trimmed_ref) <= 300:
            ref = trimmed_ref[n:]
            alt = "-"
            parsed_line_lst.append((chr, pos, ref, alt))
        else:
            pass  # not indel

    return parsed_line_lst


def right_trim(seq1, seq2):
    """Trim 3'bases if they are identical
    Args:
      seq1, seq2 (str): alleles in .vcf
    Returns:
      seq1, seq2 (str)
    """
    if len(seq1) == 1 or len(seq2) == 1:
        return seq1, seq2

    while seq1[-1] == seq2[-1] and len(seq1) > 1 and len(seq2) > 1:
        seq1, seq2 = seq1[:-1], seq2[:-1]

    return seq1, seq2


def count_padding_bases(seq1, seq2):
    """Count the number of bases padded to
    report indels in .vcf

    Args:
       seq1, seq2 (str): alleles in .vcf
    Returns:
       n (int): 

    Examples:
        REF    ALT
        
        GCG    GCGCG
        
        By 'left-alignment', REF and ATL are alignmed: 
             GCG 
             |||      
             GCGCG
       
       The first 3 bases are left-aligned.
       In this case, 3 will be returned
    """
    if len(seq2) < len(seq1):
        return count_padding_bases(seq2, seq1)

    n = 0
    for base1, base2 in zip(seq1, seq2[: len(seq1)]):
        if base1 == base2:
            n += 1
        else:
            break

    return n
