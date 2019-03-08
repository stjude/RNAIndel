#!/usr/bin/env python3
"""8th  step of analysis 

Left-align, unify equivalents and format the result

'indel_postprocessor' is the main routine of this module
"""

import sys
import logging
from .indel_preprocessor import perform_left_alignment
from .indel_annotator import annotate_indels


logger = logging.getLogger(__name__)


def indel_postprocessor(df, df_filtered, genome, exons, chr_prefixed):
    """Main routine to perform left-alingment, unification, and formatting
     
    Args:
        df (pandas.DataFrame): df with successful entries
        df_filtered (pandas.DataFrame): df with filtered entries
        genome (pysam.FastaFile): reference genome
        exons (pysam.TabixFile): coding exon data 
        chr_prefixed (bool): True if chromosome names in BAM are "chr"-prefixed
    Returns:
        df (pandas.DataFrame): df with all post-processing done
    """

    # left-alignment
    df = perform_left_alignment(df, genome, chr_prefixed)

    if not df_filtered.empty:
        df_filtered = perform_left_alignment(df_filtered, genome, chr_prefixed)
        df_filtered = df_filtered.drop_duplicates(["chr", "pos", "ref", "alt"])

    # re-classify common indels to germline
    df["predicted_class"], df["reclassified"] = zip(
        *df.apply(reclassify_common_indels, axis=1)
    )

    # reannotate afer left-alignment
    df["annotation"] = df.apply(
        annotate_indels, genome=genome, exons=exons, chr_prefixed=chr_prefixed, axis=1
    )
    df = df[df["annotation"] != "-"]

    if len(df) == 0:
        logging.warning(
            "No indels annotated in coding region after left-alignment. Analysis done."
        )
        sys.exit(0)

    df = unify_equivalent_indels(df)

    return df, df_filtered


def unify_equivalent_indels(df):
    """Unify equivalent indels with highest somatic probability

    Args:
        df (pandas.DataFrame): df after left-alignment
                               -> all equivalent indels have same 
                                  'chr', 'pos', 'ref', 'alt'
    Returns:
        df (pandas.Dataframe): de-duplicated dataframe
    """
    # to keep original order
    df["order"] = df.index

    # select one with highest somatic probability
    df = df.sort_values("prob_s", ascending=False)
    df = df.drop_duplicates(["chr", "pos", "ref", "alt"])
    df = df.sort_values("order")

    return df


def reclassify_common_indels(row):
    """Reclassify common indels predicted 'somatic' to germline

    Args:
        row (pandas.Series)
    Returns:
    """
    pred = row["predicted_class"]
    proba_somatic = row["prob_s"]
    proba_germline = row["prob_g"]
    proba_artifact = row["prob_a"]
    msg = row["reclassified"]

    if pred == "somatic" and row["is_common"] == 1:
        if (
            "Pathogenic" not in row["clin_info"]
            and "Likely_pathogenic" not in row["clin_info"]
        ):
            pred, msg = "germline", "reclassified"

    return pred, msg
