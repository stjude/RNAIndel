#!/usr/bin/env python3
"""Library for nucleotide sequence properties
"""

import numpy as np
from operator import mul
from functools import reduce


def editdistance(seq1, seq2):
    """Calculates edit distance
    
    Args:
        seq1, seq2 (str): may be empty
    Returns:
        edit distance (int)
    
    The original source of this implementation:
    https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Python
    """
    if len(seq1) < len(seq2):
        return editdistance(seq2, seq1)

    if len(seq2) == 0:
        return len(seq1)

    seq1 = np.array(tuple(seq1))
    seq2 = np.array(tuple(seq2))

    prev = np.arange(seq2.size + 1)
    for s in seq1:
        current = prev + 1

        current[1:] = np.minimum(current[1:], np.add(prev[:-1], seq2 != s))

        current[1:] = np.minimum(current[1:], current[0:-1] + 1)

        prev = current

    return prev[-1]


def linguistic_complexity(seq):
    """Quantifies the vocabulary usage of sequence.
    
    Args:
        seq (str): sequnece consists of 'A', 'G', 'T', 'C'.
    Returns:
        linguistic complexity (float)
    Raises:
        ValueError: if input string contains 'N'
       
    Example: vocabulary of ACACCA:

        1-letter vocabulary = 2; 'A', 'C'
        2-letter vocabulary = 3; 'AC', 'CA', 'CC'
        3-letter vocabulary = 4; 'ACA, 'CAC', 'ACC', 'CCA'
        ...

    Theory:
        seq: length n
        vi: i-letter vocabulary in seq

        For a DNA sequnce of length n, the maximun number of
        i-letter vocabulary is min(4**i, n-i+1).
                 
        ui = vi / min(4**i, n-i+1)

        linguistic_complexity(seq) = u1*u2*...*un

    """
    if "N" in seq:
        raise ValueError("input string may not contain 'N'.")

    n = len(seq)
    if n == 0:
        return 0
    elif n == 1:
        return 1
    else:
        usage = []
        for i in range(1, n):
            max_vocabulary = min(4 ** i, n - i + 1)

            i_mer = []
            for j in range(n - i + 1):
                i_mer.append(seq[j : j + i])
            vocabulary = len(set(i_mer))  # set() to count uniq

            usage.append(vocabulary / max_vocabulary)

        return reduce(mul, usage)


def reduce_repetitive_sequence(seq):
    """Reduces repetitive sequence to shortest non-overlapping repeat pattern
       
    Args:
        seq (str)
    Returns:
        min_unit (str)
    
    Example:
        'AGAGAGAG': repetitive patters are 'AG' and 'AGAG'.
                    this function returns 'AG' 
    """
    mid = int(len(seq) / 2)
    min_unit = seq

    j = 1
    found = False

    while j <= mid and not found:
        tandems = [seq[i : i + j] for i in range(0, len(seq), j)]
        if len(set(tandems)) == 1:
            found = True
            min_unit = list(set(tandems))[0]
        j += 1

    return min_unit


def repeat(idl_type, lt_seq, idl_seq, rt_seq):
    """Counts the repeat of the indel sequence
    
    Args:
        idl_type (int): 1 for insertion, 0 for deletion
        lt_seq (str): 5' flanking seq
        idl_seq (str): inserted or deleted seq
        rt_seq (str): 3' flanking seq
    Returns:
        repeat (int): see example below
    
    Raises:
        ValueError: if idl_type is not either 1 or 0
        ValueError: if input is None or empty str
    Example:
        Reference ATCAGAGAGAGAGAGAGCATCA
                     
        del(AGAG) ATCAG----AGAGAGAGCATCA
        
        'AGAG' reduced to 'AG'            
        Repeat of 'AG' in lt sequence: 1
        Repeat of 'AG' in rt sequence: 4
        Repeat of 'AG' in deleted seq: 2

        Return 7 
        
        For deletion, the repeat of the original seq is
        considered.
            
    """
    if idl_type != 1 and idl_type != 0:
        raise ValueError("indel type must be 1 for insertion or 0 for deletion")

    if not lt_seq or not idl_seq or not rt_seq:
        raise ValueError("Input sequences must be non-empty string.")

    idl_len = len(idl_seq)
    lt_len = len(lt_seq)
    rt_len = len(rt_seq)

    min_unit = reduce_repetitive_sequence(idl_seq)
    rep_len = len(min_unit)

    if min(lt_len, rt_len) < rep_len:
        return 0

    # replace idl_seq with min_unit if found
    original_idl_seq = ""
    if idl_seq != min_unit:
        original_idl_seq, idl_seq = idl_seq, min_unit
        idl_len = len(idl_seq)

    lt_repeat = 0
    lt_seq = lt_seq[::-1]  # reverse for easy counting
    for i in range(0, lt_len, idl_len):
        if lt_seq[i : i + idl_len] == idl_seq[::-1]:
            lt_repeat += 1
        else:
            break  # break at the 1st occurrce of non-repetitive pattern

    rt_repeat = 0
    for i in range(0, rt_len, idl_len):
        if rt_seq[i : i + idl_len] == idl_seq:
            rt_repeat += 1
        else:
            break

    # correct if deleted in repetitive region
    if idl_type == 0 and (lt_repeat + rt_repeat) > 0:
        correction = 1
    else:
        correction = 0

    # Check if indel_seq consists of a signle letter.
    # correction = 1 x len if del in repetitive region
    #            = 0  otherwise
    if original_idl_seq:
        correction = correction * len(original_idl_seq) / len(idl_seq)

    return lt_repeat + rt_repeat + correction


def dna_strength(seq):
    """Calculates DNA strength
      
    Args:
        seq (str): sequence of A, C, T, G, N with len > 1
    Returns:
        DNA strength (float): normalized by length
    Raises:
        ValueError: if input is None or str shorter than 2 
    
    Citation:   
       Khandelwal et al. 2010. 'A Phenomenological Model for 
       Predicting Melting Temperature of DNA sequences', PLOS ONE
    """
    # Adapted from Table 1 Khandelwal et al (2010).
    # Modified for N-containing di-nucleotide (averaged over possible cases)
    strength_values = {
        "GC": 13,
        "CC": 11,
        "GG": 11,
        "CG": 10,
        "AC": 10,
        "TC": 8,
        "AG": 8,
        "TG": 7,
        "GT": 10,
        "CT": 8,
        "GA": 8,
        "CA": 7,
        "AT": 7,
        "TT": 5,
        "AA": 5,
        "TA": 4,
        "AN": 7.5,
        "CN": 9,
        "GN": 10.5,
        "TN": 6,
        "NA": 6,
        "NC": 10.5,
        "NG": 9,
        "NT": 7.5,
        "NN": 8.25,
    }

    if not seq or len(seq) < 2:
        raise ValueError("Input must be string with 2-nt or longer")

    strength = 0
    for i in range(len(seq) - 1):
        # return 8.25 ('NN'-value) for 2-mer not in the dictionary
        strength = strength + strength_values.get(seq[i : i + 2], 8.25)

    return strength / len(seq)


def gc(seq):
    """ Calculates GC content

    Args:
        seq (str): non empty
    Returns:
        GC content (float)
    Raises:
        ValueError: if the sequence is None or empty str
    """
    if not seq:
        raise ValueError("Input must be non-empty string")

    g = seq.count("G")
    c = seq.count("C")

    return (g + c) / len(seq)


def dissimilarity(lt_seq, idl_seq, rt_seq):
    """Calculate how dissimilar between indel and flanking sequences.
    
    Args:
        lt_seq (str): 5' flanking seq
        indel_seq (str): inserted or deleted seq
        rt_seq (srt): 3' flanking seq
    Returns:
        dissimilarity (float)
    Raises:
        ValueError: if input is None or empty str   
    
    Example:  
       
       ref:  ATAGAAG****ATGCGGA
       data: ATAGAAGATCGATGCGGA

       The inserted seq 'ATCG' is compared with the flanking seqs 
       with the same length:GAAG and ATGC.
       
       Returns the smaller length normalized editdistance
    """
    if not lt_seq or not idl_seq or not rt_seq:
        raise ValueError("input sequences must be non-empty strings")

    n = len(idl_seq)

    lt_seq_n = lt_seq[-n:]
    rt_seq_n = rt_seq[:n]

    # in case where flanking is shorter than n,
    # normalize with len(lt(rt)_seq_n) rather than n
    lt_dissim = editdistance(lt_seq_n, idl_seq) / len(lt_seq_n)
    rt_dissim = editdistance(idl_seq, rt_seq_n) / len(rt_seq_n)

    return min(lt_dissim, rt_dissim)


def reverse_complement(seq):
    """Takes the reverse-complement of DNA seq
    
    Args:
        seq (str)
    Returns:
        seq (str): reverse complemented
    Raises:
        ValueError: if input sequence is None or empty str
    """
    if not seq:
        raise ValueError("input sequences mut be non-empty string")

    return seq[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def exists_stop_codon(strand, seq):
    """Checks if stop codon exists in sequence with 0-frame
       (0-frame: seq[i:i+3] for i = 0, 3, 6, ...)
    
    Args:
        strand (str): '+' for positive, '-' for negative
    Returns:
        stop_codon_exists (bool): True if stop codon exists, False othewise
        For strand = '-', the reverser-complement
        of the sequnce is considered.
    Raises:
        ValueError: if strand is None or invalid
        ValueError: if input sequence is None or empty string

    Example:
        For stand='+', sequence=TACTGCTATCGTCAACCATA
        frames = TAC TGC TAT CGT CAA CCA TA
        -> no stop codon
           return False
                  
        For strand = '-', sequence=same_as_above
        frames = TAT GGT TGA CGA TAG CAG TA
        -> 'TGA' exists
           return True
    """
    if strand != "+" and strand != "-":
        raise ValueError("strand must be specified either '+' or '-'")

    if not seq:
        raise ValueError("input sequence must be non-empty string")

    if strand == "-":
        seq = reverse_complement(seq)

    stop_codons = ["TAA", "TAG", "TGA"]

    stop_codon_exists = False
    for i in range(0, len(seq) - 2, 3):
        if seq[i : i + 3] in stop_codons:
            stop_codon_exists = True
            return stop_codon_exists

    return stop_codon_exists
