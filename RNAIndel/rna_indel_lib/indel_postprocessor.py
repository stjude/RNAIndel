#!/usr/bin/env python3
"""8th  step of analysis 

Left-align, unify equivalents and format the result

'indel_postprocessor' is the main routine of this module
"""

import sys
import pysam
import logging
from functools import partial
from .left_aligner import lt_aln
from .indel_sequence import Indel
from .indel_annotator import annotate_indels

logger = logging.getLogger(__name__)


def indel_postprocessor(df, df_filtered, refgene, fasta, chr_prefixed):
    """Main routine to perform left-alingment, unification, and formatting
     
    Args:
        df (pandas.DataFrame): df with successful entries
        df_filtered (pandas.DataFrame): df with filtered entries 
        refgene (str): path to refCodingExon.bed.gz
        fasta (str): path to .fa
        chr_prefixed (bool): True if chromosome names in BAM are "chr"-prefixed
    Returns:
        df (pandas.DataFrame): df with all post-processing done
    """
    fa = pysam.FastaFile(fasta)

    # left-alignment
    lt_aln_indel_generator = partial(
        generate_lt_aln_indel, fa=fa, chr_prefixed=chr_prefixed
    )
    df["lt"] = df.apply(lt_aln_indel_generator, axis=1)
    df["pos"], df["ref"], df["alt"] = zip(*df.apply(left_align_report, axis=1))

    if not df_filtered.empty:
        df_filtered["lt"] = df_filtered.apply(lt_aln_indel_generator, axis=1)
        df_filtered["pos"], df_filtered["ref"], df_filtered["alt"] = zip(
            *df_filtered.apply(left_align_report, axis=1)
        )
        df_filtered = df_filtered.drop_duplicates(["chr", "pos", "ref", "alt"])

    # re-classify common indels to germline
    df["predicted_class"], df["reclassified"] = zip(
        *df.apply(reclassify_common_indels, axis=1)
    )

    # reannotate afer left-alignment
    exon_data = pysam.TabixFile(refgene)
    anno = partial(
        annotate_indels,
        exon_data=exon_data,
        fasta=fasta,
        chr_prefixed=chr_prefixed,
        postprocess=True,
    )
    df["annotation"] = df.apply(anno, axis=1)
    df = df[df["annotation"] != "-"]

    if len(df) == 0:
        logging.warning(
            "No indels annotated in coding region after left-alignment. Analysis done."
        )
        sys.exit(0)

    df = unify_equivalent_indels(df)

    return df, df_filtered


def generate_lt_aln_indel(row, fa, chr_prefixed):
    """Generates a left-aligned Indel object

    Args:
        row (pandas.Series): with 'chr', 'pos', 'is_ins', 'indel_seq'
                            specifies original (not lt-aligned) indel  
        fa (pysam.FastaFile): obj storing reference seq
        chr_prefixed (bool): True if chromosome names in BAM are "chr"-prefixed
    Returns:
        idl (Indel obj): Indel obj left-aligned against reference
    """
    idl = Indel(row["chr"], row["pos"], row["is_ins"], row["indel_seq"])
    idl = lt_aln(idl, fa, chr_prefixed)

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
