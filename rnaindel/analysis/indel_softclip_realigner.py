#!/usr/bin/env python3

import re
import pysam
from .most_common import most_common
from .sequence_properties import repeat

cigar_ptn = re.compile(r"[0-9]+[MIDNSHPX=]")


def realn_softclips(
    reads, pos, ins_or_del, idl_seq, idl_flanks, decompose_non_indel_read
):
    template = make_indel_template(idl_seq, idl_flanks)

    candidate_reads = [
        classify_softclip_read(read, pos)
        for read in reads
        if classify_softclip_read(read, pos)
    ]
    if not candidate_reads:
        return []

    fw_decomposed = [
        forward_decomposition(
            read, softclip_ptrn, pos, ins_or_del, idl_seq, decompose_non_indel_read
        )
        for read, softclip_ptrn in candidate_reads
    ]
    rv_decomposed = [
        reverse_decomposition(read, pos, ins_or_del, idl_seq)
        for read, softclip_ptrn in candidate_reads
    ]
    decomposed_candidates = fw_decomposed + rv_decomposed

    compatible_softclip_reads = [
        decom[0]
        for decom in decomposed_candidates
        if is_compatible(decom, template, ins_or_del)
    ]

    return compatible_softclip_reads


def make_indel_template(idl_seq, idl_flanks):
    """Make consensus contig
    """
    lt_flanks = [flank[0][::-1] for flank in idl_flanks if flank[0][-1] != "N"]
    rt_flanks = [flank[1] for flank in idl_flanks if flank[1][0] != "N"]
    lt_template = find_consensus_seq(lt_flanks)[::-1]
    rt_template = find_consensus_seq(rt_flanks)

    return lt_template, idl_seq, rt_template


def get_ith_char(seq, i):
    try:
        return seq[i]
    except:
        return None


def find_consensus_seq(seq_lst):
    consensus = ""
    if not seq_lst:
        return consensus

    for i in range(len(max(seq_lst, key=len))):
        ith_chars = [get_ith_char(seq, i) for seq in seq_lst if get_ith_char(seq, i)]
        if most_common(ith_chars) == "N":
            break
        else:
            consensus += most_common(ith_chars)

    return consensus.upper()


def classify_softclip_read(read, pos):
    """Check softclip pattern and the clipped alignment is in the exon of interest
    """
    cigarstring = read.cigarstring

    if not "S" in cigarstring:
        return None

    cigarlst = cigar_ptn.findall(read.cigarstring)

    start_adjust = int(cigarlst[0][:-1]) if cigarlst[0].endswith("S") else 0
    read_start = read.reference_start - start_adjust
    end_adjust = int(cigarlst[-1][:-1]) if cigarlst[-1].endswith("S") else 0
    read_end = read.reference_end + end_adjust

    if "N" in cigarstring:
        idx_at_splicesite = [
            i for i, cigartoken in enumerate(cigarlst) if cigartoken.endswith("N")
        ]
        exonic_cigarlst = split_lst_by_index(cigarlst, idx_at_splicesite)

        # merge blocks separated by insertion/deletions
        deletion_lengths = [
            int(token[:-1]) for token in cigarlst if token.endswith("D")
        ]
        d = max(deletion_lengths) if deletion_lengths else 0
        blocks = merge_blocks(read.get_blocks(), d)

        idx_at_this_exon = []
        for i, block in enumerate(blocks):
            if i == 0 and read_start <= pos <= block[1]:
                idx_at_this_exon.append(i)
            elif i == len(blocks) - 1 and block[0] <= pos <= read_end:
                idx_at_this_exon.append(i)
            elif block[0] <= pos <= block[1]:
                idx_at_this_exon.append(i)
            else:
                pass

        if idx_at_this_exon:
            this_exon_cigarstring = exonic_cigarlst[idx_at_this_exon[0]]
        else:
            return None
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


def merge_blocks(lst, d):
    merged = []
    for i, b in enumerate(lst):
        if i <= len(lst) - 2:
            if lst[i + 1][0] <= b[1] + d:
                merged.append((b[0], lst[i + 1][1]))
                del lst[i + 1]
            else:
                merged.append((b[0], b[1]))
        else:
            if lst[i - 1][1] < b[0]:
                merged.append((b[0], b[1]))
            else:
                pass

    if merged == lst or len(merged) == 1:
        return merged
    else:
        return merge_blocks(merged, d)


def split_lst_by_index(lst, idx):

    if idx:
        idx = (0,) + tuple(data + 1 for data in idx) + (len(lst) + 1,)

    my_lst = []
    for start, end in zip(idx, idx[1:]):
        my_lst.append(lst[start : end - 1])
    return my_lst


def forward_decomposition(
    read, softclip_ptrn, pos, ins_or_del, idl_seq, decompose_non_indel_read
):
    """Decompose softclipped read from 5'side
    """
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


def reverse_decomposition(read, pos, ins_or_del, idl_seq):
    """Decompose softclipped read from 3'side
    """
    read_seq = read.query_sequence[::-1]
    cigarstring = read.cigarstring
    cigarlst = cigar_ptn.findall(read.cigarstring)[::-1]

    adjust = int(cigarlst[0][:-1]) if cigarlst[0].endswith("S") else 0
    read_pos = read.reference_end + adjust

    idl_len = len(idl_seq)

    read_idx = 0
    pos = pos if ins_or_del == "I" else pos + idl_len
    for token in cigarlst:
        event, event_len = token[-1], int(token[:-1])

        if pos < read_pos:
            read_pos = read_pos if event == "I" else (read_pos - event_len)
            read_idx = (
                read_idx if event == "D" or event == "N" else read_idx + event_len
            )
        else:
            break

    diff = read_pos - pos

    if ins_or_del == "D":
        lt_flank, mid_seq, rt_flank = (
            read_seq[read_idx + diff :],
            "",
            read_seq[: read_idx + diff],
        )
    else:
        rt_flank, mid_seq, lt_flank = (
            read_seq[: read_idx + diff],
            read_seq[read_idx + diff : read_idx + diff + idl_len],
            read_seq[read_idx + diff + idl_len :],
        )

    lt_flank, mid_seq, rt_flank = lt_flank[::-1], mid_seq[::-1], rt_flank[::-1]

    return (read, lt_flank, mid_seq, rt_flank)


def is_compatible(read_tuple, template_tuple, ins_or_del):

    read = read_tuple[0]
    read_lt_flank, read_indel, read_rt_flank = (
        read_tuple[1],
        read_tuple[2],
        read_tuple[3],
    )

    template_lt_flank, template_indel, template_rt_flank = (
        template_tuple[0],
        template_tuple[1],
        template_tuple[2],
    )

    lt_len = min(len(read_lt_flank), len(template_lt_flank))
    rt_len = min(len(read_rt_flank), len(template_rt_flank))

    # count repeat in template
    idl_type = 1 if ins_or_del == "I" else 0
    if template_lt_flank and template_rt_flank:
        template_repeat = repeat(
            idl_type, template_lt_flank, template_indel, template_rt_flank
        )
    else:
        return None

    if template_repeat > 0:
        if lt_len == 0 or rt_len == 0:
            return False
        else:
            read_repeat = repeat(idl_type, read_lt_flank, template_indel, read_rt_flank)
            if template_repeat != read_repeat:
                return False

    if ins_or_del == "D" and (lt_len == 0 or rt_len == 0):
        return False

    if lt_len > 0:
        lt_read, lt_template = read_lt_flank[-lt_len:], template_lt_flank[-lt_len:]
    else:
        lt_read, lt_template = "", ""

    rt_read, rt_template = read_rt_flank[:rt_len], template_rt_flank[:rt_len]

    if not is_almost_same(lt_read[::-1], lt_template[::-1]) or not is_almost_same(
        rt_read, rt_template
    ):
        return False

    if read_indel and ins_or_del == "I":
        template_indel_len = len(template_indel)
        read_indel_len = len(read_indel)
        if template_indel_len < read_indel_len:
            return False
        elif read_indel == template_indel:
            return True
        elif 4 <= template_indel_len <= 5:
            return identical_for_end_n_bases(read_indel, template_indel, 2)
        elif 6 <= template_indel_len <= 7:
            return identical_for_end_n_bases(read_indel, template_indel, 3)
        else:
            return identical_for_end_n_bases(read_indel, template_indel, 4)
    elif ins_or_del == "D":
        return True
    else:
        return False


def identical_for_end_n_bases(query_str, subject_str, n):
    return (query_str[:n] == subject_str[:n]) or (query_str[-n:] == subject_str[-n:])


def is_almost_same(seq1, seq2, len_lim=10, mismatch_lim=1):
    seq_len = len(seq1)
    if seq_len > 0 and seq1[0] != seq2[0]:
        return False

    hamming = sum([seq1[i] != seq2[i] for i in range(seq_len)])
    if seq_len >= len_lim:
        return hamming <= mismatch_lim
    else:
        return hamming == 0
