import re
import os
import sys
import numpy as np

cigar_ptrn = re.compile(r"[0-9]+[MIDNSHPX=]")


def most_common(lst):
    return max(set(lst), key=lst.count)


def flatten_list_of_list(lst_of_lst):
    return [i for lst in lst_of_lst for i in lst]


def validate_int_inputs(val, preset=None):
    val = int(val)
    if val <= 0 and not preset:
        sys.exit("Error: the input must be a positive integer.")
    elif 0 <= val <= 255 and preset == "mapq":
        sys.exit("Error: the MAPQ value must be an integer between 0 and 255.")

    return val


def validate_file_input(file_path):
    if os.path.isfile(file_path):
        return file_path
    else:
        sys.exit("Error: {} not found.".format(file_path))


def validate_dir_input(dir_path):
    if os.path.isdir(dir_path):
        return dir_path
    else:
        sys.exit("Error: {} not found.".format(dir_path))


def adjust_start_pos(aligned_segment, is_for_ref):
    start = aligned_segment.reference_start + 1  # to 1-based

    if is_for_ref:
        return start
    else:
        cigar_list = cigar_ptrn.findall(aligned_segment.cigarstring)
        start_offset = int(cigar_list[0][:-1]) if cigar_list[0].endswith("S") else 0
        start = start - start_offset
        return start


def get_ref_seq(aligned_segment, target_indel, cigar_string, cigar_list):

    aln_start, aln_end = (
        aligned_segment.reference_start + 1,
        aligned_segment.reference_end,
    )
    chrom, reference = target_indel.chrom, target_indel.reference

    current_pos = aln_start - 1

    if not "N" in cigar_string:
        return reference.fetch(chrom, current_pos, aln_end)

    ref_seq = ""

    for cigar in cigar_list:
        event, event_len = cigar[-1], int(cigar[:-1])
        if event == "M" or event == "D":
            ref_seq += reference.fetch(chrom, current_pos, current_pos + event_len)
            current_pos += event_len
        elif event in ("I", "S", "H", "P"):
            pass
        else:
            current_pos += event_len

    return ref_seq


def split(aligned_segment, target_indel, is_for_ref):

    target_pos = target_indel.pos
    string_pos = adjust_start_pos(aligned_segment, is_for_ref)

    cigarstring = aligned_segment.cigarstring
    cigar_lst = cigar_ptrn.findall(cigarstring)

    data = (
        get_ref_seq(aligned_segment, target_indel, cigarstring, cigar_lst)
        if is_for_ref
        else aligned_segment.query_sequence
    )

    _size = len(cigar_lst)
    data_moves = np.zeros((_size,))
    genome_moves = np.zeros((_size,))

    i, j = 0, 0
    for cigar in cigar_lst:
        event, event_len = cigar[-1], int(cigar[:-1])

        if event == "N":
            d_move = 0
            g_move = event_len
        elif event == "I":
            g_move = 0
            d_move = 0 if is_for_ref else event_len
        elif event == "D":
            g_move = event_len
            d_move = event_len if is_for_ref else 0
        elif event == "H" or event == "P":
            d_move = 0
            g_move = 0
        else:
            g_move, d_move = event_len, event_len

        data_moves[i] = d_move
        genome_moves[i] = g_move
        i += 1

    string_pos -= 1

    for d_move, g_move in zip(data_moves, genome_moves):
        if string_pos < target_pos:
            string_pos += g_move
        else:
            break
        j += d_move

    diff = target_pos - string_pos
    idx = int(j + diff)

    lt = data[:idx]
    rt = data[idx:]

    return lt, rt
