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


def indel_rescuer(df, fasta, bam, chr_prefixed, num_of_processes, **kwargs):
    """Search indels using locus info from caller and rescue them by equivalence 
    or by nearest

    Args:
        df (pandas.DataFrame): left-aligned indels
        fasta (str): path to reference fasta (cannot use pysam.FastaFile here)
        bam (str): path to bam (cannot use pysam.AlignmentFile here)
        chr_prefixed (bool): True if chromosome names in BAM are "chr" prefixed
        num_of_processes (int): the num of processes in parallelism
    Returns:
        df (pandas.DataFrame): indels rescued 
    """
    external_vcf = kwargs.pop("external_vcf", False)
    pool = Pool(num_of_processes)
    
    # rescue by equivalence
    df["rescued_indels"] = df.apply(
        rescue_by_equivalence,
        fasta=fasta,
        bam=bam,
        search_window=50,
        pool=pool,
        chr_prefixed=chr_prefixed,
        axis=1,
    )
    df["rescued"] = df.apply(flag_indel_rescued_by_equivalence, axis=1)

    # rescue by nearest
    if external_vcf:
        rqxnr = partial(
            rescue_by_nearest, fasta=fasta, bam=bam, chr_prefixed=chr_prefixed
        )
        # not applied if already rescued by equivalence
        df["rescued_indels"] = df.apply(
            lambda x: x["rescued_indels"] if x["rescued_indels"] else rqxnr(x), axis=1
        )
    df["rescued"] = df.apply(flag_indel_rescued_by_nearest, axis=1)

    # create and format dataframe
    list_of_data_dict = df["rescued_indels"].sum()
    df_rescued = pd.DataFrame(list_of_data_dict)
    df = df[["chr", "pos", "ref", "alt", "rescued"]]
    df = pd.concat([df, df_rescued], axis=0, sort=True)
    df = sort_positionally(df)
    df = df.drop_duplicates(["chr", "pos", "ref", "alt"])
    df.reset_index(drop=True, inplace=True)

    return df


def rescue_by_equivalence(row, fasta, bam, search_window, pool, chr_prefixed):
    """Recover equivalent indels from left-aligned indel report
    
    Args:
        row (pandas.Series)
        fasta (str): path to reference fasta 
                     (CANNOT pass pysam.FastaFile(fasta), which is not pickleble)
        bam (str): path to bam 
                   (CANNOT pass pysam.AlignmentFile(bam), which is not pickleble)
        search_window (int): to define search range
        pool (multiprocessing.Pool)
        chr_prefixed (bool): True if chromosome names in BAM are "chr" prefixed
    Returns
        equivalents (list): dict element
                            {
                              'chr':chromosome,
                              'pos':1-based pos,
                              'ref':ref allele,
                              'alt':alt allele,
                              'rescued':rescue status
                             }

                            empty list if no equivalent indels found
    """
    chr, pos = row["chr"], row["pos"]

    if row["ref"] == "-":
        idl_type, idl_seq = 1, row["alt"]
    else:
        idl_type, idl_seq = 0, row["ref"]
    
    called_idl = curate_indel_in_genome(
        pysam.FastaFile(fasta), chr, pos, idl_type, idl_seq, chr_prefixed
    )
    
    rescue = partial(
        extract_indel,
        fasta=fasta,
        bam=bam,
        chr=chr,
        idl_type=idl_type,
        chr_prefixed=chr_prefixed,
        equivalent_to=called_idl,
    )
    
    # equivalent indels only exist on 5' side as left-aligned
    equivalents = pool.map(rescue,  [pos + i for i in range(search_window)])

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
    flag = "by_equivalence" if row["rescued_indels"] else "-"
    return flag


def rescue_by_nearest(row, fasta, bam, chr_prefixed):
    """Recover nearest non-equivalent indel
       within an internally defined search window
       (this is done only if rescue by equivalence failed)

    Args:
        row (pandas.Series)
        fasta (str): path to reference fasta 
        bam (str): path to bam
        chr_prefixed (bool): True if chromosome names are "chr" prefixed
    Returns
        nearests (list): dict element
                         {
                            'chr':chromosome,
                             'pos':1-based pos,
                             'ref':ref allele,
                             'alt':alt allele,
                             'rescued': rescue status
                         }             
                         empty list if no equivalent indels found
    """
    chr, pos = row["chr"], row["pos"]

    idl_type = 1 if row["ref"] == "-" else 0
    idl_size = len(row["alt"]) if idl_type else len(row["ref"])
    search_window = 10 if idl_size < 4 else 20
    
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
        genome = pysam.FastaFile(fasta)

        chr = idl_found.chr
        pos = idl_found.pos
        ref = idl_found.ref
        alt = idl_found.alt
        idl_vcf = IndelVcfReport(genome, chr, pos, ref, alt, chr_prefixed)
        idl_vcf.left_align()
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
    flag = (
        row["rescued_indels"][0]["rescued"]
        if row["rescued_indels"] and row["rescued"] != "by_equivalence"
        else row["rescued"]
    )
    return flag


def extract_indel(pos, fasta, bam, chr, idl_type, chr_prefixed, **kwargs):
    """Extract equivalent indel if exists at the locus (chr, pos)

    Args:
        pos (int): 1-based coordinate
                   placed as 1st param for map()
        fasta (str): path to reference fasta
        bam (str): path to bam
        chr (str): chr1-22, chrX or chrY. Note "chr"-prefixed.
        idl_type (int): 1 for insertion, 0 for deletion
        chr_prefixed (bool): True if chrosomome names in BAM is "chr"-prefixed
    KwArgs:
        equivalent_to (SequenceWithIndel obj): indel obj to be compared for equivalence
    Returns:
        target_indel(SequenceWithIndel or None) :SequenceWithIndel if found
                                                     None if not found
    """

    target_indel = None
    
    genome = pysam.FastaFile(fasta)
    alignments = pysam.AlignmentFile(bam)

    idl_to_compare = kwargs.pop("equivalent_to", None)

    collected_seqs = collect_indel_seqs(alignments, chr, pos, idl_type, chr_prefixed)
    # equivalence search
    if idl_to_compare and collected_seqs:
        indels_at_this_locus = [
            curate_indel_in_genome(genome, chr, pos, idl_type, idl_seq, chr_prefixed)
            for idl_seq in set(collected_seqs)
        ]
        for idl in indels_at_this_locus:
            if idl_to_compare == idl:
                target_indel = idl

    # nearest indel search
    if not idl_to_compare and collected_seqs:
        target_indel = curate_indel_in_genome(
            genome, chr, pos, idl_type, most_common(collected_seqs), chr_prefixed
        )
    return target_indel


def collect_indel_seqs(alignments, chr, pos, idl_type, chr_prefixed):
    """Extract most frequent indel sequnece from bam data

    Args:
        alignments (pysam.AlignmentFile): bam data
        chr (str): chr1-22, chrX or chrY. Note "chr"-prefixed. 
        pos (int): 1-based coordinate
        idl_type (int): 1 for insertion, 0 for deletion
        chr_prefixed (bool): True if chromosome names in BAM are "chr"-prefixed
    Returns:
        idl_seqs (list): empty list if no indels found
    """
    idl_seqs = []

    ins_or_del = "I" if idl_type == 1 else "D"

    # convert 0-based coordinate
    pos = pos - 1

    try:
        valid_reads = extract_all_valid_reads(alignments, chr, pos, chr_prefixed)
    except:
        return idl_seqs

    try:
        parsed_indel_reads = extract_indel_reads(valid_reads, pos, ins_or_del)
    except:
        return idl_seqs

    if not parsed_indel_reads:
        return idl_seqs
    else:
        decomposed = [
            decompose_indel_read(parsed_read) for parsed_read in parsed_indel_reads
        ]

    idl_seqs = [decomp[1] for decomp in decomposed]
    
    return idl_seqs


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
