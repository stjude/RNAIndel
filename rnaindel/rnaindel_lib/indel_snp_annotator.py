#!/usr/bin/env python3
"""6th step of analysis

Annotate indels if they are on dbSNP

'indel_snp_annotator' is the main routine of this module
"""
import re
import pysam
from .indel_features import IndelSnpFeatures
from .indel_curator import curate_indel_in_genome
from .indel_vcf_preprocessor import right_trim
from .indel_vcf_preprocessor import parse_vcf_line
from .indel_vcf_preprocessor import count_padding_bases


def indel_snp_annotator(df, genome, dbsnp, clinvar, chr_prefixed):
    """Annotates indels with dbSNP and ClinVar info

    Args:
        df (pandas.DataFrame)
        genome (pysam.FastaFile): reference genome
        dbsnp (pysam.TabixFile): dbSNP database
        clinvar (pysam.TabixFile): ClinVar database
        chr_prefixed (bool): True if chromosome names in BAM are "chr"-prefixed
    Returns:
        df (pandas.DataFrame): with SNP info
    """
    df["db"] = df.apply(
        annotate_indel_on_db,
        genome=genome,
        dbsnp=dbsnp,
        clinvar=clinvar,
        chr_prefixed=chr_prefixed,
        axis=1,
    )
    df["dbsnp"] = df.apply(lambda x: x["db"].report_dbsnp_id(), axis=1)
    df["is_on_dbsnp"] = df.apply(is_on_dbsnp, axis=1)
    df["max_maf"] = df.apply(lambda x: x["db"].report_freq(), axis=1)
    df["is_common"] = df.apply(lambda x: x["db"].is_common(), axis=1)
    # df['is_not_pathogenic'] = df.apply(lambda x: x['db'].is_not_pathogenic(), axis=1)
    # df['with_germline_reports'] = df.apply(lambda x: x['db'].with_germline_reports(), axis=1)
    df["clin_info"] = df.apply(lambda x: x["db"].report_clnvr_info(), axis=1)
    df["is_on_dbsnp"] = df.apply(negate_on_dbsnp_if_pathogenic, axis=1)

    df.drop("db", axis=1, inplace=True)

    return df


def annotate_indel_on_db(row, genome, dbsnp, clinvar, chr_prefixed):
    """Check dbSNP and ClinVar membership by equivalence and annotate.

    Args:
        row (pandas.Series): with 'chr', 'pos', 'is_ins', 'indel_seq' lables
        genome (pysam.FastaFile): reference genome
        dbsnp (pysam.TabixFile): dbSNP database
        clinvar (pysam.TabixFile): ClinVar database
        chr_prefixed (bool): True if chromosome names in BAM are "chr"-prefixed
    Returns:
        report (IndelSnpFeatures): idl object reporting SNP info
    """
    chr = row["chr"]
    pos = row["pos"]
    idl_type = row["is_ins"]
    idl_seq = row["indel_seq"]

    # obj representing the indel in reference genome
    idl = curate_indel_in_genome(genome, chr, pos, idl_type, idl_seq, chr_prefixed)
    # obj representing report of the indel
    report = IndelSnpFeatures(chr, pos, idl_type, idl_seq)

    # search for equivalent indels over pos +/- search_window nt
    search_window = 50
    start, end = pos - search_window, pos + search_window
    chr_vcf = row["chr"].replace("chr", "")

    for record in dbsnp.fetch(chr_vcf, start, end, parser=pysam.asTuple()):
        bambinos = vcf2bambino(record)
        for bb in bambinos:
            if idl_type == bb.idl_type and len(idl_seq) == len(bb.idl_seq):
                # indel on db representing in reference genome
                db_idl = curate_indel_in_genome(
                    genome, chr, bb.pos, bb.idl_type, bb.idl_seq, chr_prefixed
                )
                if idl == db_idl:
                    rs = record[2]
                    report.add_dbsnp_id(rs)
                    report.add_dbsnp_freq(dbsnp_freq(record))
                    # report.add_dbsnp_origin(dbsnp_origin(record))
                    report.add_dbsnp_common(dbsnp_common(record))

    for record in clinvar.fetch(chr_vcf, start, end, parser=pysam.asTuple()):
        bambinos = vcf2bambino(record)
        for bb in bambinos:
            if idl_type == bb.idl_type and len(idl_seq) == len(bb.idl_seq):
                db_idl = curate_indel_in_genome(
                    genome, chr, bb.pos, bb.idl_type, bb.idl_seq, chr_prefixed
                )
                if idl == db_idl:
                    id = record[2]
                    report.add_clnvr_id(id)
                    report.add_clnvr_freq(clnvr_freq(record))
                    # report.add_clnvr_origin(clnvr_origin(record))
                    report.add_clnvr_info(cln_info(record))

    return report


def is_on_dbsnp(row):
    """Encodes if the indel is found on dbSNP 
    
    Args:
        row (pandas.Series): with 'dbsnp' label
    Returns:
        is_on_dbsnp (int): 1 if yes 0 othewise
    """
    is_on_dbsnp = 0 if row["dbsnp"] == "-" else 1

    return is_on_dbsnp


def negate_on_dbsnp_if_pathogenic(row):
    """Returns 0 if the indel is on dbSNP but rated
    'Pathogenic' or 'Likely pathogenic'. This is to
    prevent pathogenic (possibly tumorigenic) indels
    from being classified into germline.
    
    Args:
        row (pandas.Series): with 'clin_info' and 'is_on_dbsnp' labels 
    Returns:
        is_on_dbsnp (int): 0 if indel is pathogenic, 
                           otherwise as annotated in row['is_on_dbsnp']
    """
    is_on_dbsnp = row["is_on_dbsnp"]

    if "Pathogenic" in row["clin_info"] or "Likely_pathogenic" in row["clin_info"]:
        is_on_dbsnp = 0

    return is_on_dbsnp


def vcf2bambino(record):
    """Generate IndelSnpFeatures obj. with bambibo coordinate from VCF record
    
    Args:
       record (tuple): vcf line with fields separated in tuple
    Returns:
       parsed (list): a list of IndelSnpFeatures obj
    """
    vcf_line = "\t".join(
        ["chr" + record[0], record[1], record[2], record[3], record[4]]
    )

    parsed = parse_vcf_line(vcf_line)

    indels_in_bambino_format = []
    for elem in parsed:
        chr, pos, ref, alt = elem[0], elem[1], elem[2], elem[3]
        idl_type = 1 if ref == "-" else 0
        idl_seq = alt if idl_type else ref

        indels_in_bambino_format.append(IndelSnpFeatures(chr, pos, idl_type, idl_seq))

    return indels_in_bambino_format


def dbsnp_freq(record):
    """Collects minor allele frequency and common SNP annotation
   
    Args:
       record (tuple): vcf line with fields separated in tuple
    Returns:
       allele_freq (float): the larger value of 1000G and TOPMED
                            -1 if freq info not avail
       is_common (int): 1 if annotated as COMMON
                        0 if not annotated as COMMON
                       -1 if annoation not avail
       origin (int): 0 for unspecified origin
                     1 for germline
                     2 for somatic
                     3 for both
                    -1 for annotation not avail
   """
    try:
        kg = re.search(r"(CAF=)([0-9,.e-]+)", record[7]).group(2)
        kg_af = float(kg.split(",")[1])
    except:
        kg_af = -1

    try:
        topmed = re.search(r"(TOPMED=)([0-9,.e-]+)", record[7]).group(2)
        topmed_af = float(topmed.split(",")[1])
    except:
        topmed_af = -1

    return max(kg_af, topmed_af)


def clnvr_freq(record):
    """Collects minor allele frequency (MAF) reported in ClinVar

    Args:
        record (tuple): vcf line with fields separated in tuple
    Returns:
        max(esp, exac, tgp) (float): return -1 if MAF is not reported
    """
    try:
        esp = float(re.search(r"(AF_ESP=)([0-9,.e-]+)", record[7]).group(2))
    except:
        esp = -1

    try:
        exac = float(re.search(r"(AF_EXAC=)([0-9,.e-]+)", record[7]).group(2))
    except:
        exac = -1

    try:
        tgp = float(re.search(r"(AF_TGP=)([0-9,.e-]+)", record[7]).group(2))
    except:
        tgp = -1

    return max(esp, exac, tgp)


def dbsnp_common(record):
    """Annotates if the indel is a commmon SNP

    Args:
         record (tuple): vcf line with fields separated in tuple
    Returns:
         common (int): 1 if common
                       0 if uncommon
                      -1 if info is not available
    """
    try:
        common = int(re.search(r"(COMMON=)([01])", record[7]).group(2))
    except:
        common = -1

    return common


def dbsnp_origin(record):
    """Annotates the origin of SNP
   
    Args:
        record (tuple): vcf line with fields separated in tuple
    Returns:
        origin (int): 0 if unspecified
                      1 if germline
                      2 if somatic
                      3 if both
                     -1 if info is not available
    """

    try:
        origin = int(re.search(r"(SAO=)([0123])", record[7]).group(2))
    except:
        origin = -1

    return origin


def clnvr_origin(record):
    """Annotate the origin of SNP in ClinVar

    Args:
        record (tuple): vcf line with fields separated in tuple
    Returns:
        origin (int): 0 if unknown
                      1 if germline
                      2 if somatic
                      4 if inherited
                      8 if parternal
                      16 if maternal
                      32 if de-novo
                      64 if biparental
                      128 if uniparental
                      256 if not-tested
                      512 if tested-inconclusive
                      1073741824 if other
                      -1 if info is not available
    """
    try:
        origin = int(re.search(r"(ORIGIN=)([0-9]+)", record[7]).group(2))
    except:
        origin = -1

    return origin


def cln_info(record):
    """Annotates pathogenicity using ClinVar
    
    Args:
        record (tuple): vcf line with fields separated in tuple
    Returns:
        significance + '|' + disease (str): pathogenicity and disease name implicated
    """
    try:
        significance = re.search(r"(CLNSIG=)([A-Za-z_]+)", record[7]).group(2)
    except:
        significance = "clnSigNA"

    try:
        disease = re.search(r"(CLNDN=)([A-Za-z0-9_,-]+)", record[7]).group(2)
    except:
        disease = "diseaseNA"

    return significance + "|" + disease
