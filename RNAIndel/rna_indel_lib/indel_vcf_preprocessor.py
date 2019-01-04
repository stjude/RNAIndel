#!/usr/bin/env python3
"""Optional step for .vcf input

Convert .vcf to a format compatible with Bambino

'indel_vcf_preprocessor' is the main routine of this module
"""

import sys
import vcf
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
    vcf_data = vcf.Reader(open(vcffile, "r"))
    bam_data = pysam.AlignmentFile(bam)
    exon_data = pysam.TabixFile(refgene)

    chr_prefixed = is_chr_prefixed(bam_data)
    df = pd.DataFrame(make_data_list(vcf_data))
    
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
    data = []
    for record in vcf_data:
        parsed_record = parse_vcf_record(record)
        if parsed_record != None:
            d = {
                "chr": parsed_record[0],
                "pos": parsed_record[1],
                "ref": parsed_record[2],
                "alt": parsed_record[3],
            }

            data.append(d)

    return data


def parse_vcf_record(record):
    """
    Args:
        record (vcf.Record): record in .vcf 
    Returns:
        parsed_record (tuple): (chr, pos, ref, alt) 
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
    parsed_record = None

    chr = record.CHROM  # from input VCF, may or may not be "chr"-prefixed
    if not chr.startswith("chr"):
        chr = "chr" + chr

    if record.is_indel and is_canonical_chromosome(chr):
        ref = record.REF
        alts = record.ALT
        for alt in alts:
            alt = str(alt)  # cast to str from _Substitution obj
            n = count_padding_bases(ref, alt)
            pos = record.POS + n  # converts to Bambino coordinate

            if len(ref) < len(alt):
                ref = "-"
                alt = alt[n:]
            else:
                ref = ref[n:]
                alt = "-"

            parsed_record = (chr, pos, ref, alt)

    return parsed_record
