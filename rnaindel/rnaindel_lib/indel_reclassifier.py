#!/usr/bin/env python3
"""Optional step for reclassification.

Knowledge-based reclassification is performed
for indels predicted somatic.

'indel_reclassifier' is the main routine in this module
"""
import pysam
from functools import partial
from .indel_snp_annotator import vcf2bambino
from .indel_curator import curate_indel_in_genome


def indel_reclassifier(df, genome, pons, chr_prefixed):
    """Main module to reclassify based on user-defined 
    panel of non somatic (PONS). 
    
    Args:
        df (pandas.DataFrame): df with prediction made
        genome (pysam.FastaFile): reference genome 
        pons (pysam.TabixFile): user-defined VCF non-somatic indels
        chr_prefixed (bool): True if chromosome names in BAM are "chr"-prefixed
    Returns:
        df (pandas.DataFrame): df reclassified
    """
    reclf = partial(
        wrap_reclassify_by_pons, genome=genome, pons=pons, chr_prefixed=chr_prefixed
    )
    df["predicted_class"], df["reclassified"] = zip(*df.apply(reclf, axis=1))

    return df


def wrap_reclassify_by_pons(row, genome, pons, chr_prefixed):
    """Wrap 'relassify_by_panel_of_non_somatic' so that this function is
    only applied to instances predicted 'somatic'

    Args: see 'relassify_by_panel_of_non_somatic'
    Returns: see 'relassify_by_panel_of_non_somatic'
    """
    if row["predicted_class"] == "somatic" and row["is_common"] != 1:
        return relassify_by_panel_of_non_somatic(row, genome, pons, chr_prefixed)
    else:
        return row["predicted_class"], row["reclassified"]


def relassify_by_panel_of_non_somatic(row, genome, pons, chr_prefixed):
    """Reclassifies indels predicted somatic using user-defined
    panel of non somatic (PONS).

    Args:
        row (pandas.Series)
        genome (pysam.FastaFile): reference genome
        pons (pysam.TabixFile obj): user-defined database for non-somatic indels
        chr_prefixed (bool): True if chromosome names in BAM are "chr"-prefixed
    Returns:
        'predicted_class' (str): reclassifed class if applicable
        'comment' (str): 'reclassified' if reclassified, '-' otherwise 
    """
    search_window = 50

    chr = row["chr"]
    pos = row["pos"]
    idl_type = row["is_ins"]
    idl_seq = row["indel_seq"]

    idl = curate_indel_in_genome(genome, chr, pos, idl_type, idl_seq, chr_prefixed)

    # check if contif names in vcf are prefixed with 'chr'
    sample_contig = pons.contigs[0]
    if not sample_contig.startswith("chr"):
        chr_vcf = chr.replace("chr", "")
    else:
        chr_vcf = chr

    start, end = pos - search_window, pos + search_window

    # check if the indel is equivalent to indel on the panel of non somatic (PONS)
    # reclassify based on the 2nd highest probability if equivalent PONS indel found
    for record in pons.fetch(chr_vcf, start, end, parser=pysam.asTuple()):
        bambinos = vcf2bambino(record)
        for bb in bambinos:
            if idl_type == bb.idl_type and len(idl_seq) == len(bb.idl_seq):
                pons_idl = curate_indel_in_genome(
                    genome, chr, bb.pos, bb.idl_type, bb.idl_seq, chr_prefixed
                )

                if idl == pons_idl:
                    if row["prob_a"] >= row["prob_g"]:

                        return "artifact", "reclassified"
                    else:
                        return "germline", "reclassified"

    return row["predicted_class"], row["reclassified"]
