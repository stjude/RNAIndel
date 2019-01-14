#!/usr/bin/env python3
import pysam
import numpy as np
import pandas as pd
from functools import partial
from multiprocessing import Pool
from .most_common import most_common
from .indel_vcf import IndelVcfReport
from .indel_curator import extract_indel_reads
from .indel_curator import decompose_indel_read
from .indel_curator import curate_indel_in_genome
from .indel_curator import extract_all_valid_reads


def indel_rescuer(df, fasta, bam, chr_prefixed, **kwargs):

    num_of_processes = kwargs.pop("num_of_processes", 1)
    left_aligned = kwargs.pop("left_aligned", False)
    external_vcf = kwargs.pop("external_vcf", False)

    pool = Pool(num_of_processes)

    df["rescued"] = "-"

    # rescue by equivalence
    rqxeq = partial(
        rescue_by_equivalence,
        fasta=fasta,
        bam=bam,
        search_window=50,
        pool=pool,
        left_aligned=left_aligned,
        chr_prefixed=chr_prefixed,
    )

    df["rescued_indels"] = df.apply(rqxeq, axis=1)
    df["rescued"] = df.apply(flag_indel_rescued_by_equivalence, axis=1)

    # rescue by nearest
    if external_vcf:
        rqxnr = partial(
            rescue_by_nearest,
            fasta=fasta,
            bam=bam,
            search_window=10,
            chr_prefixed=chr_prefixed,
        )
        df["rescued_indels"] = df.apply(
            lambda x: rqxnr(x) if x["rescued_indels"] == [] else x["rescued_indels"],
            axis=1,
        )
    df["rescued"] = df.apply(flag_indel_rescued_by_nearest, axis=1)

    list_of_data_dict = df["rescued_indels"].sum()
    df_rescued = pd.DataFrame(list_of_data_dict)

    df = df[["chr", "pos", "ref", "alt", "rescued"]]

    df = pd.concat([df, df_rescued], axis=0, sort=True)

    df = sort_positionally(df)
    df = df.drop_duplicates(["chr", "pos", "ref", "alt"])
    df.reset_index(drop=True, inplace=True)

    return df


def rescue_by_equivalence(
    row, fasta, bam, search_window, pool, left_aligned, chr_prefixed
):
    """Recover equivalent indels from left-aligned indel report
    
    Args:
        row (pandas.Series)
        fasta (str): path to fasta
        bam (str): path to bam
        search_window (int): to define search range
        pool (Multiprocessing.Pool obj)
        left_aligned (bool): True if input indels are left-aligned.
        chr_prefixed (bool): True if chromosome names are "chr" prefixed
    Returns
        equivalents (list): dict element
                            {'chr':chromosome,
                             'pos':1-based pos,
                             'ref':ref allele,
                             'alt':alt allele}

                            empty list if no equivalent indels found
    """
    chr = row["chr"]  # this is "chr"-prefixed
    pos = row["pos"]

    if row["ref"] == "-":
        idl_type, idl_seq = 1, row["alt"]
    else:
        idl_type, idl_seq = 0, row["ref"]

    called_idl = curate_indel_in_genome(
        fasta, chr, pos, idl_type, idl_seq, chr_prefixed
    )

    if left_aligned:
        rt_window, lt_window = search_window, search_window - search_window
    else:
        rt_window, lt_window = int(search_window / 2), int(search_window / 2)

    rescue = partial(
        extract_indel,
        fasta=fasta,
        bam=bam,
        chr=chr,
        idl_type=idl_type,
        chr_prefixed=chr_prefixed,
        equivalent_to=called_idl,
    )

    rt_range = [pos + i for i in range(rt_window)]
    rt_equivalents = pool.map(rescue, rt_range)

    lt_range = [pos - i for i in range(lt_window)]
    lt_equivalents = pool.map(rescue, lt_range)

    equivalents = rt_equivalents + lt_equivalents

    equivalents = [
        {
            "chr": eq.chr,
            "pos": eq.pos,
            "ref": eq.ref,
            "alt": eq.alt,
            "rescued": "by_equivalence",
        }
        for eq in equivalents
        if eq
    ]

    return equivalents


def flag_indel_rescued_by_equivalence(row):
    flag = row["rescued"]

    if row["rescued_indels"] != []:
        flag = "by_equivalence"

    return flag


def rescue_by_nearest(row, fasta, bam, search_window, chr_prefixed):

    chr = row["chr"]
    pos = row["pos"]
    idl_type = 0

    if row["ref"] == "-":
        idl_type = 1

    # make [0, 1, -1, 2, -2, ...]
    pos_move = [
        -int(i) if int(i) == i else int(i + 1)
        for i in np.array(range(search_window)) / 2
    ]
    i = 0
    idl_found = None
    while i < len(pos_move) and not idl_found:
        idl_found = extract_indel(
            pos + pos_move[i], fasta, bam, chr, idl_type, chr_prefixed
        )
        i += 1

    if not idl_found:
        return []
    else:
        chr = idl_found.chr
        pos = idl_found.pos
        ref = idl_found.ref
        alt = idl_found.alt
        fa = pysam.FastaFile(fasta)
        idl_vcf = IndelVcfReport(fa, chr, pos, ref, alt, chr_prefixed)
        in_vcf_style = (
            idl_vcf.CHROM
            + ":"
            + str(idl_vcf.POS)
            + ":"
            + idl_vcf.REF
            + ":"
            + idl_vcf.ALT
        )
        return [
            {
                "chr": chr,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "rescued": "rescued_by:" + in_vcf_style,
            }
        ]


def flag_indel_rescued_by_nearest(row):
    flag = row["rescued"]
    if row["rescued_indels"] != [] and flag != "by_equivalence":
        flag = row["rescued_indels"][0]["rescued"]

    return flag


def extract_indel(pos, fasta, bam, chr, idl_type, chr_prefixed, **kwargs):
    """Extract equivalent indel if exists at the locus (chr, pos)

    Args:
        pos (int): 1-based coordinate
                   Placed as 1st arg for map()
        fasta (str): path to fasta
        bam (str): path to bam
        chr (str): chr1-22, chrX or chrY. Note "chr"-prefixed.
        idl_type (int): 1 for insertion, 0 for deletion
        chr_prefixed (bool): True if chrosomome names in BAM is "chr"-prefixed
    KwArgs:
        equivalent_to (SequenceWithIndel obj): indel obj to be compared for equivalence
    Returns:
        equivalent_indel(SequenceWithIndel or None) :SequenceWithIndel if found
                                                     None if not found
    """
    idl_at_this_locus = None
    bam_data = pysam.AlignmentFile(bam, "rb")
    idl_to_compare = kwargs.pop("equivalent_to", None)

    inferred_idl_seq = get_most_common_indel_seq(
        bam_data, chr, pos, idl_type, chr_prefixed
    )

    if inferred_idl_seq:
        idl_at_this_locus = curate_indel_in_genome(
            fasta, chr, pos, idl_type, inferred_idl_seq, chr_prefixed
        )

        if idl_to_compare:
            if idl_to_compare != idl_at_this_locus:
                idl_at_this_locus = None

    return idl_at_this_locus


def get_most_common_indel_seq(bam_data, chr, pos, idl_type, chr_prefixed):
    """Extract most frequent indel sequnece from bam data

    Args:
        bam_data (pysam.pysam.AlignmentFile)
        chr (str): chr1-22, chrX or chrY. Note "chr"-prefixed. 
        pos (int): 1-based coordinate
        idl_type (int): 1 for insertion, 0 for deletion
        chr_prefixed (bool): True if chromosome names in BAM are "chr"-prefixed
    Returns:
        idl_seq (str or None): None type if no indels found
    """
    idl_seq = None

    if idl_type == 1:
        ins_or_del = "I"
    else:
        ins_or_del = "D"

    # convert 0-based coordinate
    pos = pos - 1

    try:
        valid_reads = extract_all_valid_reads(bam_data, chr, pos, chr_prefixed)
    except:
        return idl_seq

    try:
        parsed_indel_reads = extract_indel_reads(valid_reads, pos, ins_or_del)
    except:
        return idl_seq

    if not parsed_indel_reads:
        return idl_seq
    else:
        decomposed = [
            decompose_indel_read(parsed_read) for parsed_read in parsed_indel_reads
        ]

    idl_seq = most_common([decomp[1] for decomp in decomposed])

    return idl_seq


def sort_positionally(df):
    df["chr"] = df.apply(lambda x: x["chr"].replace("chr", ""), axis=1)
    df["chr"] = df.apply(lambda x: 23 if x["chr"] == "X" else x["chr"], axis=1)
    df["chr"] = df.apply(lambda x: 24 if x["chr"] == "Y" else x["chr"], axis=1)
    df["chr"] = df.apply(lambda x: int(x["chr"]), axis=1)

    df.sort_values(["chr", "pos"], inplace=True)

    df["chr"] = df.apply(lambda x: "Y" if x["chr"] == 24 else x["chr"], axis=1)
    df["chr"] = df.apply(lambda x: "X" if x["chr"] == 23 else x["chr"], axis=1)
    df["chr"] = df.apply(lambda x: "chr" + str(x["chr"]), axis=1)

    return df
