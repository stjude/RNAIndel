#!/usr/bin/env python3
"""A python implementation of left-alignment algorithm
by Tan et al 2015 Bioinformatics, 31:2202-2204
"""

import pysam
from .indel_sequence import Indel


def lt_aln(idl, fa):
    """Perfoms left alignment 
    
    Args:
        idl (Indel obj)
        fa (pysam.FastaFile obj)
    Returns:
        idl (Indel obj)
    """
    while idl.idl_seq[-1] == peek_left_base(idl, fa):
        idl = shift_to_left(idl, fa)

    return idl


def shift_to_left(idl, fa):
    """Move indel to the left by one nucleotide

    pos  012345   6     01234   56
    ref: ACTCCT   G  -> ACTCC   TG
    alt: ACTCCTCCTG     ACTCCTCCTG
      ins(CCT) at 6  -> ins(TCC) at 5 
     
    Args:
        idl (Indel obj)
        fa (pysam.FastaFile obj)
    Returns:
        idl (Indel obj)
    """
    left_base = peek_left_base(idl, fa)
    # shift to the left by 1 nucleotide 
    idl.idl_seq = left_base + idl.idl_seq[:-1]
    idl.pos = idl.pos - 1
    
    return idl 


def peek_left_base(idl, fa):
    """Get the nucleotide immediately left to the indel
         
         pos:  0123  4567      
         ref:  ATGC  GCTA
         alt:  ATGCTAGCTA
         
         'TA' is inserted at 4.
         the left base is 'C' at 3

    Args:
        idl (Indel obj)
        fa (pysam.FastaFile obj)
    Returns:
        left_base (str)
    """
    left_base = fa.fetch(idl.chr, idl.pos - 2, idl.pos - 1)
    
    return left_base
