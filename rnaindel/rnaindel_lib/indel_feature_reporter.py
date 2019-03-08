#!/usr/bin/env python3
import pandas as pd
from .indel_vcf import IndelVcfReport


def indel_feature_reporter(df, genome, reportname, chr_prefixed):
    """Report calculated feature values for training
    
    Args:
        df (pandas.DataFrame)
        genome (pysam.FastaFile): reference genome
        reportname (str): filename of the ouput report
        chr_prefixed (bool): True if chromosome names in BAM are chr-prefixed
    Returns:
        None (a tab-delimited file will be generated)
    """
    df["vcf"] = df.apply(to_vcfreport, genome=genome, chr_prefixed=chr_prefixed, axis=1)

    # reformat to VCF style
    df["chr"] = df.apply(lambda x: x["vcf"].CHROM, axis=1)
    df["pos"] = df.apply(lambda x: x["vcf"].POS, axis=1)
    df["ref"] = df.apply(lambda x: x["vcf"].REF, axis=1)
    df["alt"] = df.apply(lambda x: x["vcf"].ALT, axis=1)
    
    # add truth column. this is to be filled by user
    df["truth"] = ""


    df.drop(
        [
            "rescued",
            "indel_seq",
            "filtered",
            "dbsnp",
            "max_maf",
            "is_common",
            "clin_info",
            "vcf",
        ],
        inplace=True,
        axis=1,
    )

    df.to_csv(reportname, sep="\t", index=False)


def to_vcfreport(row, genome, chr_prefixed):
    """Convert bamabino call to VCF record

    Args:
        row (pandas.Series)
        genome (pysam.FastaFile): reference genome
        chr_prefixed (bool): True if chromosome names in BAM are chr-prefixed
    Returns
        IndelVcfReport (class): indel call in VCF format
    """
    idl_vcf = IndelVcfReport(
        genome, row["chr"], row["pos"], row["ref"], row["alt"], chr_prefixed
    )

    return idl_vcf
