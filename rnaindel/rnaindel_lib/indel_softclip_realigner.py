#!/usr/bin/env python3

import re
import pysam
from .most_common import most_common
from .indel_curator import decompose_non_indel_read

cigar_ptn = re.compile(r"[0-9]+[MIDNSHPX=]")


def realn_softclips(reads, pos, ins_or_del, idl_seq, idl_flanks):

    candidate_reads = [
        classify_softclip_read(read, pos)
        for read in reads
        if classify_softclipped_read(read, pos)
    ]

    if not candidate_reads:
        return []

    decomposed_candidates = [
        decompose_softclip(read, softclip_ptrn, pos, ins_or_del, idl_seq)
        for read, softclip_ptrn in candidate_reads
    ]

    template = make_indel_template(idl_seq, idl_flanks)
    compatible_softclip_reads = [
        decom[0]
        for decom in decomposed_candidates
        if is_compatible(decom[1:], template, ins_or_del)
    ]

    return compatible_softclip_reads


def decompose_softclip(read, softclip_ptrn, pos, ins_or_del, idl_seq):

    decom = decompose_non_indel_read(read, pos, ins_or_del, idl_seq)
    lt_flank, mid_seq, rt_flank = decom[2][0], decom[1], decom[2][1]

    if ins_or_del == "I" and softclip_ptrn == "leading":
        mid_seq = lt_flank[-len(idl_seq) :]
        lt_flank = lt_flank[: -len(idl_seq)]
    elif ins_or_del == "I" and softclip_ptrn == "trailing":
        mid_seq = rt_flank[: len(idl_seq)]
        rt_flank = rt_flank[len(idl_seq) :]
    elif ins_or_del == "D" and softclip_ptrn == "leading":
        lt_flank = lt_flank + mid_seq
        mid_seq = idl_seq
    else:
        rt_flank = mid_seq + rt_flank
        mid_seq = idl_seq

    return (read, lt_flank, mid_seq, rt_flank)


def classify_softclip_read(read, pos):
    cigarstring = read.cigarstring

    if not "S" in cigarstring:
        return None

    cigarlst = cigar_ptn.findall(read.cigarstring)
    if "N" in cigarstring:
        idx_at_splicesite = [
            i for i, cigartoken in enumerate(cigarlst) if cigartoken.endswith("N")
        ]
        exonic_cigarlst = split_lst_by_index(cigarlst, idx_at_splicesite)
        idx_at_this_exon = [
            i
            for i, block in enumerate(read.get_blocks())
            if block[0] <= pos <= block[1]
        ]
        this_exon_cigarstring = exonic_cigarlst[idx_at_this_exon[0]]
    else:
        this_exon_cigarstring = cigarlst

    first, last = this_exon_cigarstring[0][-1], this_exon_cigarstring[-1][-1]

    if first == "S" and last != "S":
        return (read, "leading")
    elif first != "S" and last == "S":
        return (read, "trailing")
    elif first == "S" and last == "S":
        # give up this pattern for now
        return None
    else:
        return None


def split_lst_by_index(lst, idx):

    if idx:
        idx = (0,) + tuple(data + 1 for data in idx) + (len(lst) + 1,)

    my_lst = []
    for start, end in zip(idx, idx[1:]):
        my_lst.append(lst[start : end - 1])
    return my_lst


def make_indel_template(idl_seq, idl_flanks):
    lt_flanks = [flank[0][::-1] for flank in idl_flanks]
    rt_flanks = [flank[1] for flank in idl_flanks]
    lt_template = find_consensus_seq(lt_flanks)[::-1]
    rt_template = find_consensus_seq(rt_flanks)

    return lt_template, idl_seq, rt_template


def get_ith_char(seq, i):
    try:
        return seq[i]
    except:
        return None


def find_consensus_seq(seq_lst):
    """
    lst of str
    """
    consensus = ""
    for i in range(len(max(seq_lst, key=len))):
        ith_chars = [get_ith_char(seq, i) for seq in seq_lst if get_ith_char(seq, i)]
        if most_common(ith_chars) == "N":
            break
        else:
            consensus += most_common(ith_chars)

    return consensus.upper()


def is_compatible(read_tuple, template_tuple, ins_or_del):

    read_lt_flank, read_indel, read_rt_flank = (
        read_tuple[0],
        read_tuple[1],
        read_tuple[2],
    )

    template_lt_flank, template_indel, template_rt_flank = (
        template_tuple[0],
        template_tuple[1],
        template_tuple[2],
    )

    lt_len = min(len(read_lt_flank), len(template_lt_flank))
    rt_len = min(len(read_rt_flank), len(template_rt_flank))

    if lt_len > 0:
        lt_read, lt_template = read_lt_flank[-lt_len:], template_lt_flank[-lt_len:]
    else:
        lt_read, lt_template = "", ""

    rt_read, rt_template = read_rt_flank[:rt_len], template_rt_flank[:rt_len]

    if not is_almost_same(lt_read, lt_template) or not is_almost_same(
        rt_read, rt_template
    ):
        return False

    if read_indel and ins_or_del == "I":
        subject_len = len(template_indel)
        query_len = len(read_indel)
        if subject_len < query_len:
            return False
        elif query_indel == subject_indel:
            return True
        elif 4 <= subject_len <= 5:
            return identical_for_end_n_bases(query_indel, subject_indel, 2)
        elif 6 <= subject_len <= 7:
            return identical_for_end_n_bases(query_indel, subject_indel, 3)
        else:
            return identical_for_end_n_bases(query_indel, subject_indel, 4)
    elif ins_or_del == "D":
        return True
    else:
        return False


def identical_for_end_n_bases(query_str, subject_str, n):
    return (query_str[:n] == subject_str[:n]) or (query_str[-n:] == subject_str[-n:])


def is_almost_same(seq1, seq2, len_lim=10, mismatch_lim=1):
    seq_len = len(seq1)
    hamming = sum([seq1[i] != seq2[i] for i in range(seq_len)])
    if seq_len >= len_lim:
        return hamming <= mismatch_lim
    else:
        return hamming == 0
