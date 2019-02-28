#!/usr/bin/env python3
import pysam
import pandas as pd
from .indel_vcf import IndelVcfReport


def report_features(df, fasta, reportname, chr_prefixed):
    fa = pysam.FastaFile(fasta)
    df["vcf"] = df.apply(to_vcfreport, fa=fa, chr_prefixed=chr_prefixed, axis=1)

    # reformat to VCF style
    df["chr"] = df.apply(lambda x: x["vcf"].CHROM, axis=1)
    df["pos"] = df.apply(lambda x: x["vcf"].POS, axis=1)
    df["ref"] = df.apply(lambda x: x["vcf"].REF, axis=1)
    df["alt"] = df.apply(lambda x: x["vcf"].ALT, axis=1)

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


def to_vcfreport(row, fa, chr_prefixed):
    idl_vcf = IndelVcfReport(
        fa, row["chr"], row["pos"], row["ref"], row["alt"], chr_prefixed
    )

    return idl_vcf
