#!/usr/bin/env python3
"""A python implementation of left-alignment algorithm
by Tan et al 2015 Bioinformatics, 31:2202-2204
"""


def lt_aln(idl, genome, chr_prefixed):
    """Perfoms left alignment 
    
    Args:
        idl (Indel obj)
        genome (pysam.FastaFile): reference genome
        chr_prefixed (bool): True if chromosome names are "chr"-prefixed
    Returns:
        idl (Indel obj)
    """
    while idl.idl_seq[-1] == peek_left_base(idl, genome, chr_prefixed):
        idl = shift_to_left(idl, genome, chr_prefixed)

    return idl


def shift_to_left(idl, genome, chr_prefixed):
    """Move indel to the left by one nucleotide
   
    Args:
        idl (Indel obj)
        genome (pysam.FastaFile): reference genome 
        chr_prefixed (bool): True if chromosome names are "chr"-prefixed 
    Returns:
        idl (Indel obj)
    Example:
        pos  012345   6     01234   56
        ref: ACTCCT   G  -> ACTCC   TG
        alt: ACTCCTCCTG     ACTCCTCCTG
        ins(CCT) at 6  -> ins(TCC) at 5 
    """
    left_base = peek_left_base(idl, genome, chr_prefixed)
    # shift to the left by 1 nucleotide
    idl.idl_seq = left_base + idl.idl_seq[:-1]
    idl.pos = idl.pos - 1

    return idl


def peek_left_base(idl, genome, chr_prefixed):
    """Get the nucleotide immediately left to the indel
         
    Args:
        idl (Indel obj)
        genome (pysam.FastaFile): reference genome
        chr_prefixed (bool): True if chromosome names in BAM are "chr"-prefixed
    Returns:
        left_base (str)
    Example:
         pos:  0123  4567      
         ref:  ATGC  GCTA
         alt:  ATGCTAGCTA
         
         'TA' is inserted at 4.
         the left base is 'C' at 3
    """
    chr = idl.chr
    if not chr_prefixed:
        chr = chr.replace("chr", "")
    left_base = genome.fetch(chr, idl.pos - 2, idl.pos - 1)

    return left_base
