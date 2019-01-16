#!/usr/bin/env python3
"""Optional step for .vcf input

Convert .vcf to a format compatible with Bambino

'indel_vcf_preprocessor' is the main routine of this module
"""

import sys
import pysam
import logging
import pandas as pd
from functools import partial
from .indel_preprocessor import flag_coding_indels
from .indel_snp_annotator import count_padding_bases
from .indel_preprocessor import is_chr_prefixed
from .indel_preprocessor import is_canonical_chromosome

logger = logging.getLogger(__name__)


def indel_vcf_preprocessor(vcffile, bam, refgene, fasta):
    """Convert input VCF to Bambino format and check chromosome name format
    
    Args:
        vcffile (str): path to input vcf
        bam (str): path to bam
        refgene (str): path to refCodingExon.bed.gz
        fasta (str): path to fasta
    Returns:
        df (pandas.DataFrame): df with indels reported as in Bambino output
    """
    vcf_data = open(vcffile)
    df = pd.DataFrame(make_data_list(vcf_data))
    vcf_data.close()

    bam_data = pysam.AlignmentFile(bam)
    exon_data = pysam.TabixFile(refgene)

    chr_prefixed = is_chr_prefixed(bam_data)

    datasize = len(df)
    if len(df) == 0:
        logging.warning("No indels detected in input vcf. Analysis done.")
        sys.exit(0)

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


def make_data_list(vcf_data):
    """
    Args:
        vcf_data (File obj.): VCF content 
    Returns:
        data (list): with dict element 
    """
    data = []
    for line in vcf_data:
        parsed_line = parse_vcf_line(line)
        if parsed_line:
            d = {
                "chr": parsed_line[0],
                "pos": parsed_line[1],
                "ref": parsed_line[2],
                "alt": parsed_line[3],
            }

            data.append(d)

    return data


def parse_vcf_line(line):
    """
    Args:
        line (str): line in VCF file obj.
    Returns:
        parsed_line (tuple): (chr, pos, ref, alt) 
                             Bambino compatible format
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
    parsed_line = None

    # skip header lines
    if line.startswith("#"):
        return parsed_line

    lst = line.rstrip().split("\t")
    chr = lst[0]
    pos = int(lst[1])
    ref = lst[3]
    alts = lst[4].split(",")  # possibly multi-allelic

    if not chr.startswith("chr"):
        chr = "chr" + chr

    # skip non canonical chrmosomes
    if not is_canonical_chromosome(chr):
        return parsed_line

    for alt in alts:
        n = count_padding_bases(ref, alt)
        pos += n

        if len(ref) < len(alt):
            ref = "-"
            alt = alt[n:]
        elif len(ref) > len(alt):
            ref = ref[n:]
            alt = "-"
        else:
            return parsed_line  # not indel

        parsed_line = (chr, pos, ref, alt)

    return parsed_line
