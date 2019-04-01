#!/usr/bin/env python3

import re
import random
import numpy as np
from .most_common import most_common
from .indel_sequence import SequenceWithIndel
from .indel_sequence import PileupWithIndel
from .indel_sequence import PileupWithIndelNotFound

random.seed(123)
cigar_ptn = re.compile(r"[0-9]+[MIDNSHPX=]")


def curate_indel_in_genome(genome, chr, pos, idl_type, idl_seq, chr_prefixed):
    """Gerenates an indel object with reference flanking sequences.
       Splicing will NOT be considered.
       
    Args:
        genome (pysam.FastaFile): reference genome
        chr (str): chr1-22, chrX, chrY. Note "chr"-prefixed. 
        pos (int): 1-based 
        idl_type (int): 1 for insertion 0 for deletion
        idl_seq (str): inserted or deleted sequence
        chr_prefixed (bool): True if chromosome names in BAM or FASTA is prefixed with "chr"
    Returns:
        SequenceWithIndel (obj)

    Example:
        The read below is spliced (>) and has an 'A' deletion (*)
        at chr1:13-13.
      
        Pos:       12345678901234567890123456789012345
        Reference: CGTATGATGTTCAGGTATGCGTATATAGAAAATCGA
                               *
        Read:      CGTATGATGTTC-G>>>>>>>>>>>>>>AAAATCGA

        SequenceWithIndel obj. represents the deleted 'A' with
        flanking sequence like this:
              
                    TATGATGTTC del(A) GGTATGCGTA
        
        Note that the right flanking sequence is not spliced.             
    """
    # extract flanking seq +- window(nt)
    window = 50

    if not chr_prefixed:
        chr = chr.replace("chr", "")

    # left flank seq
    start, end = pos - 1 - window, pos - 1
    lt_seq = genome.fetch(chr, start, end)

    # right flank seq
    # for insertion
    if idl_type:
        start, end = pos - 1, pos - 1 + window
        rt_seq = genome.fetch(chr, start, end)

    # for deletion
    else:
        size = len(idl_seq)
        start, end = pos - 1 + size, pos - 1 + size + window
        rt_seq = genome.fetch(chr, start, end)

    # adding back "chr" if removed above
    if not chr_prefixed:
        chr = "chr" + chr

    return SequenceWithIndel(chr, pos, idl_type, lt_seq, idl_seq, rt_seq)


def is_close_to_exon_boundary(cigarstring, idx):
    """Checks if indel is within 2-nt to the exon boundary.
    
    Args:
        cigarstring (str): see Example below
        idx (int): list index to point which CIGAR token 
                   specifies the indel of interest
    Returns:
        is_close (int): 1 for true othewise 0

    Example:
         The distance between indel and the splice site
         is 1-nt. This is 'close'
         
    Reference: CGTATGATGTTCAGGTATGCGTATATAGAAAATCGA
           
    Read:         ATGACGTTC-G>>>>>>>>>>>>>>AAAATCGA
    cigarstring: 9M1D1M14N8M         
    cigarlst: ['9M, '1D', '1M', '14N', '8M']
    idx: 1 
         del('A') is specified by '1D' which is 
         indexed 1 in CIGAR_LST
    """
    cigarlst = cigar_ptn.findall(read.cigarstring)

    is_close = 0
    # dist to 5' exon boundary
    if idx >= 2 and "N" in cigarlst[idx - 2]:
        dist_to_lt_boundary = int(cigarlst[idx - 1].replace("M", ""))
        if dist_to_lt_boundary <= 2:
            is_close = 1

    # dist to 3' exon boundary
    elif (idx + 2) <= len(cigarlst) - 1 and "N" in cigarlst[idx + 2]:
        dist_to_rt_boundary = int(cigarlst[idx + 1].replace("M", ""))
        if dist_to_rt_boundary <= 2:
            is_close = 1
    else:
        pass

    return is_close


def extract_all_valid_reads(alignments, chr, pos, chr_prefixed):
    """Extracts reads that are
        1. non-duplicate
        2. primary alignment
        3. covering the locus specified by chr and pos (non-skipping)
    
    Args:
        alignments (pysam.AlignmentFile): bam data
        chr (str): chr1-22, chrX or chrY. Note "chr"-prefixed
        pos (int): 0-based coordinate
        chr_prefixed (bool): True if chromosome names are "chr"-prefixed
    Returns:
        valid_reads (list): a list of pysam.AlignedSegment
    
    Example:
           locus of interest: chr1:13-13
           
           Extract all reads except for   
               Read_1_dup (duplicate)  
               Read_4 (non-primary)
               Read_5 (not covering, skipping)
 
           Chr: 1
           Pos:       12345678901234567890123456789012345
           Reference: CGTATGATGTTCAGGTATGCGTATATAGAAAATCGA
                             ^    *              
           Read_1        ATGACGTTC-G>>>>>>>>>>>>>>AAAATCGA
           Read_1_dup    ATGACGTTC-G>>>>>>>>>>>>>>AAAATCGA
           Read_2      GTATGACGTTCAG>>>>>>>>>>>>>>AAAAT
           Read_3     CGTATGACGTTC-G>>>>>>>>>>>>>>AAA
           Read_4            CGTTC-G>>>>>>>>>>>>>>AAATCGA (non-primary)   
           Read_5     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    """
    chr = chr if chr_prefixed else chr.replace("chr", "")

    all_reads = alignments.fetch(chr, pos, pos + 1, until_eof=True)

    valid_reads = []
    for read in all_reads:
        # excludes duplicate or non-primary alignments
        if not read.is_duplicate and not read.is_secondary:
            blocks = read.get_blocks()
            for block in blocks:
                # excludes skipping reads
                if block[0] <= pos <= block[1]:
                    valid_reads.append(read)

    return valid_reads


def extract_indel_reads(reads, pos, ins_or_del):
    """Extract reads with indel at locus specified by chr and pos

    Args:
        reads (list): [pysam.AlignedSegment]
        pos (int): 0-based coordinate
        ins_or_del (str): 'I' for insertion or 'D' insertion
    Returns:
        parsed_indel_reads (list): a list of (pysam.AligedSegment obj, idx, adjust)
                            idx: the index of cigar token specifying the indel
                            abjust: the number of 5' soft-clipped bases
    Example:
            pos          012345678901234567
            reference:   AATGATAGAAGGATGATG
            read:        aaT--AAG-AGGATGATG
            Cigar:       ['2S', '1M', '2D', '3M', '1D', '9M']
            
            The 'A' deletion at 8 is specified by '1D', whose idx = 4
            The first 2 bases are soft-clipped. adjust = 2.  
    """
    parsed_indel_reads = []
    for read in reads:
        ref_pos = read.reference_start
        cigarstring = read.cigarstring

        if ins_or_del in cigarstring:
            cigarlst = cigar_ptn.findall(cigarstring)

            # adjust ref_start if the read starts with soft-clipping
            # Pos on Genome:   0 1 2 3 4 5 6 7 8
            # reference:           A A A C C C G
            # read:                a a A C C - G (first 2 small 'a' soft-clipped)
            # ref_pos = 4 (before adjustment)
            # adjusted ref_pos = 2
            adjust = 0
            if cigarlst[0].endswith("S"):
                adjust = int(cigarlst[0].replace("S", ""))
                ref_pos = ref_pos - adjust

            for idx, cigartoken in enumerate(cigarlst):
                # cigartoken example :'34M' (34 bases mapped)
                # cigartoken value = 34
                # cigartoken operation = 'M' (mapped)
                ope = cigartoken[-1] # cigartoken operation
                val = int(cigartoken.replace(ope, ""))  # cigartoken value

                if ref_pos == pos and ope == ins_or_del:
                    parsed_indel_reads.append((read, idx, adjust))
                elif ope == "I":  # if ins, no move on reference
                    ref_pos = ref_pos
                else:
                    ref_pos = ref_pos + val

    return parsed_indel_reads


def decompose_indel_read(parsed_indel_read):
    """Decompose read and ref sequences 
    into flanking and inserted/deleted sequences

    Args:
        parsed_indel_read (tuple): (pysam.AlignedSegment, idx, adjust)
    Returns:
        decomposed_reads (tuple): (
                                   pysam.AlignedSegment, 
                                   idl_seq (str), 
                                   read_flanks (list),
                                   ref_flanks(list)
                                  ) 
    Example:
        This sample is carrying a G>A SNP (*).
         
                         *
        reference   ATGAGGTAGATAGAT
        read        ATGAGAT---TAGAT            
        ref         ATGAGGT---TAGAT
     
        idl_seq = 'AGA'
        read_flanks = ['ATGAGAT', 'TAGAT']
        ref_flanks = ['ATGAGGT', 'TAGAT']
    """
    read = parsed_indel_read[0]
    idx = parsed_indel_read[1]
    adjust = parsed_indel_read[2]

    cigarlst = cigar_ptn.findall(read.cigarstring)
    cigarlst_up_to_the_indel = cigarlst[:idx]
    ins_or_del = cigarlst[idx][-1]
    last_move = int(cigarlst[idx].replace(ins_or_del, ""))

    # read (actual sequence)
    read_seq = read.query_sequence
    # reference sequence
    ref_seq = read.get_reference_sequence()

    i = 0  # pos on read_seq
    j = 0 - adjust  # pos on ref_seq
    for cigartoken in cigarlst_up_to_the_indel:
        ope = cigartoken[-1]
        val = int(cigartoken.replace(ope, ""))

        if ope == "N":  # if spliced, no move
            i = i
            j = j
        elif ope == "D":  # if del, move on ref_seq
            i = i
            j = j + val
        elif ope == "I":
            i = i + val  # if ins, move on read_seq
            j = j
        else:
            i = i + val
            j = j + val

    if ins_or_del == "I":
        lt_read = read_seq[:i].upper()
        rt_read = read_seq[i + last_move :].upper()

        lt_ref = ref_seq[:j].upper()
        rt_ref = ref_seq[j:].upper()

        idl_seq = read_seq[i : i + last_move]
    else:
        lt_read = read_seq[:i].upper()
        rt_read = read_seq[i:].upper()

        lt_ref = ref_seq[:j].upper()
        rt_ref = ref_seq[j + last_move :].upper()

        idl_seq = ref_seq[j : j + last_move]

    read_flank = [lt_read, rt_read]
    ref_flank = [lt_ref, rt_ref]

    return (read, idl_seq, read_flank, ref_flank)


def decompose_non_indel_read(read, pos, ins_or_del, idl_seq):
    """Decompose non-indel reads into flanking and
    sequence observed at the indel event locus.

    Args:
        read (pysam.AligedSegment)
        pos (int): 0-based
        ins_or_del (str): 'I' or 'D'
        idl_seq (str)
    Returns:
        decomposed read (tuple): (
                                  pysam.AligedSegment,
                                  idl_seq (str),
                                  non_indel_flanks (list)
                                 )
    Example: 
                   01234567890123
        reference: CAGCAGCATCAGCA
     non_idl_read: CAGCAGCAGCAGCA 
     
     input: read = non_idl_read
            pos = 6
            ins_or_del = 'D'
            idl_seq = 'CAT'
     
     output: pysam.AligedSegment = read
             idl_seq = 'CAG'
             non_indel_flanks = ['CAGCAG', 'CAGCA']
      
     For insertion, idl_seq = '-'
    """
    ref_pos = read.reference_start
    read_seq = read.query_sequence
    cigarstring = read.cigarstring
    cigarlst = cigar_ptn.findall(cigarstring)

    # adjust if the read starts with softclipping
    adjust = 0
    if cigarlst[0].endswith("S"):
        adjust = int(cigarlst[0].replace("S", ""))
        ref_pos = ref_pos - adjust

    i = 0
    ref_pos
    for cigartoken in cigarlst:
        val = int(cigartoken.replace(cigartoken[-1], ""))
        ope = cigartoken[-1]

        if ref_pos < pos:
            if ope == "I":
                ref_pos = ref_pos
            else:
                ref_pos = ref_pos + val

            if ope == "N":
                i = i
            elif ope == "D":
                i = i
            elif ope == "I":
                i = i + val
            else:
                i = i + val
        else:
            break

    diff = pos - ref_pos
    if ins_or_del == "I":
        lt = read_seq[: i + diff]
        rt = read_seq[i + diff :]
        idl_seq = "-"
    else:
        size = len(idl_seq)
        lt = read_seq[: i + diff]
        rt = read_seq[i + diff + size :]
        # recover deleted sequnece from the data
        idl_seq = read_seq[i + diff : i + diff + size]

    non_idl_flanks = [lt, rt]

    return read, idl_seq, non_idl_flanks


def is_near_exon_boundary(parsed_indel_read):
    """Checks if the dist to the exon boundary is
    within threshold.
    
    Args:
        cigarstring (str): see Example below
        idx (int): list index to point which CIGAR token 
                   specifies the indel of interest
    Returns:
        is_close (int): 1 for true othewise 0

    Example:
         The distance between indel and the splice site
         is 1-nt. This is 'close'
         
    Reference: CGTATGATGTTCAGGTATGCGTATATAGAAAATCGA
           
    Read:         ATGACGTTC-G>>>>>>>>>>>>>>AAAATCGA
    cigarstring: 9M1D1M14N8M         
    cigarlst: ['9M, '1D', '1M', '14N', '8M']
    idx: 1 
         del('A') is specified by '1D' which is 
         indexed 1 in CIGAR_LST
    """
    read = parsed_indel_read[0]
    idx = parsed_indel_read[1]
    cigarlst = cigar_ptn.findall(read.cigarstring)
    idl_size = int(cigarlst[idx].replace(cigarlst[idx][-1], ""))

    if idl_size <= 2:
        threshold = 2
    else:
        threshold = 3

    is_near = 0
    # dist to 5' exon boundary
    if idx >= 2 and "N" in cigarlst[idx - 2]:
        dist_to_lt_boundary = int(cigarlst[idx - 1].replace("M", ""))
        if dist_to_lt_boundary <= threshold:
            is_near = 1

    # dist to 3' exon boundary
    elif (idx + 2) <= len(cigarlst) - 1 and "N" in cigarlst[idx + 2]:
        dist_to_rt_boundary = int(cigarlst[idx + 1].replace("M", ""))
        if dist_to_rt_boundary <= threshold:
            is_near = 1
    else:
        pass

    return is_near


def infer_del_seq_from_data(decomposed_non_idl_reads, idl_flanks, del_seq):
    """Incorporate individual difference in deleted sequence.
    
    Args:
        decomposed_non_idl_reads (list): a list of non indel read decomposition
        idl_read_flanks (list): a list of [lt_flank, rt_flank]
        del_seq (str): deleted sequence derived from ref.
    Return:
        inferred_seq (str)

    Example:                   *
        reference      CAGCAGCATCAGCAGCAG
        non_idl_read   CAGCAGCAGCAGCAGCAG       
        idl_read       CAGCAG---CAGCAGCAG
       
    Deleted sequence from the reference is 'CAT'. But this sample
    is carrying synonymous variant CAT>CAG (*). From the data
    the deleted seq is more likely to be 'CAG'.  
    """
    inferred_seq = del_seq
    size = len(del_seq)
    reads = [decomp[0] for decomp in decomposed_non_idl_reads]
    recovered_del_seq = [decomp[1] for decomp in decomposed_non_idl_reads]
    non_idl_flanks = [decomp[2] for decomp in decomposed_non_idl_reads]

    # comparison with reference (del_seq)
    lt_comparison_with_ref = [del_seq == flank[0][-size:] for flank in idl_flanks]
    rt_comparison_with_ref = [del_seq == flank[1][:size] for flank in idl_flanks]

    if True in lt_comparison_with_ref or True in rt_comparison_with_ref:
        return inferred_seq

    inferred_ptn = []
    for recovered, flank, read in zip(recovered_del_seq, non_idl_flanks, reads):

        lt_comparison_with_data = [
            recovered == flank[0][-size:] for flank in idl_flanks
        ]
        lt_flank_len = len(flank[0])

        rt_comparison_with_data = [recovered == flank[0][:size] for flank in idl_flanks]
        rt_flank_len = len(flank[1])

        if (
            True in lt_comparison_with_data
            and lt_flank_len > 9
            and not "S" in read.cigarstring
        ):
            inferred_ptn.append(recovered)

        elif (
            True in rt_comparison_with_data
            and rt_flank_len > 9
            and not "S" in read.cigarstring
        ):
            inferred_ptn.append(recovered)

        else:
            pass

    if inferred_ptn:
        inferred_seq = most_common(inferred_ptn)

    return inferred_seq


def curate_indel_in_pileup(alignments, chr, pos, idl_type, idl_seq, mapq, chr_prefixed):
    """Abstract indels in the alignment pileup
    
    Args:    
        alignments (pysam.AlignmentFile): bam data
        chr (str): chr1-22, chrX or chrY. Note "chr"-prefixed
        pos (int): 1-based position of indel on the reference
        idl_type (int): 1 for insertion 0 for deletion
        idl_seq (str): inserted or deleted sequence
        mapq (int): MAPQ for uniquely mapped reads
        chr_prefixed (bool): True if chromosome names in BAM are "chr"-prefixed
    Returns:
        PileupWithIndel object: if indels found as specified with 
                                chr, pos, idl_type and idl_seq 
        
        PileupWithIndelNotFound object: otherwise
    """
    # convert to 0-based coordinate
    pos = pos - 1

    # convert indel type to CIGAR token
    if idl_type == 1:
        ins_or_del = "I"  # as specified by input
        del_or_ins = "D"  # the other pattern
    else:
        ins_or_del = "D"
        del_or_ins = "I"

    # extract all good reads covering the locus of interest
    all_reads = extract_all_valid_reads(alignments, chr, pos, chr_prefixed)

    ###########################
    # Analysis of indel reads #
    ###########################

    # extract all indel reads at the specified pos
    # this may contain indel reads with indel sequences
    # other than specified by the input 'idl_seq'
    idl_reads = extract_indel_reads(all_reads, pos, ins_or_del)

    #
    # sanity check by position and pattern
    #
    if idl_reads == []:
        return PileupWithIndelNotFound(chr, pos, idl_type, idl_seq)

    # decompose indel read into indel sequence and flanking sequences
    decomposed_idl_reads = [decompose_indel_read(idl_read) for idl_read in idl_reads]

    # filter decomposed indels by 'idl_seq'
    filtered_decomposed_idl_reads = [
        decomp for decomp in decomposed_idl_reads if decomp[1] == idl_seq
    ]

    #
    # sanity check by indel sequence
    #
    if filtered_decomposed_idl_reads == []:
        return PileupWithIndelNotFound(chr, pos, idl_type, idl_seq)

    # extract indel read flanking
    idl_flanks = [decomp[2] for decomp in filtered_decomposed_idl_reads]

    # extract reference flanking
    ref_flanks = [decomp[3] for decomp in filtered_decomposed_idl_reads]

    is_multiallelic = 0
    # check for multiallelic by heterogeneous indel sequences
    if len(decomposed_idl_reads) > len(filtered_decomposed_idl_reads):
        is_multiallelic = 1

    # check for multiallelic by co-occurrence of insertion and deletion
    idl_reads_with_the_other_pattern = extract_indel_reads(
        all_reads, pos, del_or_ins
    )  # note the last arg!
    if idl_reads_with_the_other_pattern:
        is_multiallelic = 1

    # check for nearness to the exon boundary
    exon_boundary = [is_near_exon_boundary(idl_read) for idl_read in idl_reads]

    # collect mapping quality
    map_qual = [idl_read[0].mapping_quality for idl_read in idl_reads]

    # check bidirectionality
    bidirectional = [idl_read[0].is_reverse for idl_read in idl_reads]

    ###############################
    # Analysis of non-indel reads #
    ###############################

    # make dict to access read sequenc by name
    reads_by_name = {read.query_name: read for read in all_reads}

    # collect non-indel read by name
    all_read_names = [read.query_name for read in all_reads]
    idl_read_names = [decomp[0].query_name for decomp in filtered_decomposed_idl_reads]

    non_idl_read_names = list(set(all_read_names) - set(idl_read_names))
    non_idl_read_names.sort()

    # sample 10 non-indel reads if too many
    if len(non_idl_read_names) > 10:
        non_idl_read_names = random.sample(non_idl_read_names, 10)

    # decompose non-indel reads
    decomposed_non_idl_reads = [
        decompose_non_indel_read(reads_by_name[name], pos, ins_or_del, idl_seq)
        for name in non_idl_read_names
    ]

    # collect non indel flankings
    if non_idl_read_names == [] and len(ref_flanks) > 10:
        non_idl_flanks = random.sample(ref_flanks, 10)
    elif non_idl_read_names == []:
        non_idl_flanks = ref_flanks
    else:
        non_idl_flanks = [decomp[2] for decomp in decomposed_non_idl_reads]

    # infer deleted seq from data
    if ins_or_del == "D":
        del_seq = idl_seq
        idl_seq = infer_del_seq_from_data(decomposed_non_idl_reads, idl_flanks, del_seq)

    ########################
    # Summarize the result #
    ########################

    # fragment count by unifiying the read name
    total_count, alt_count = len(set(all_read_names)), len(set(idl_read_names))
    ref_count = total_count - alt_count

    # decide if the indel is close to exon boundary
    is_near_boundary = most_common(exon_boundary)

    # decide if uniquely mapped
    is_uniq_mapped = 0
    if most_common(map_qual) == mapq:
        is_uniq_mapped = 1

    # decide if bidirectionally supported
    is_bidirectional = 0
    if len(set(bidirectional)) == 2:
        is_bidirectional = 1

    # back to 1-based coordinate
    pos = pos + 1

    return PileupWithIndel(
        chr,
        pos,
        idl_type,
        idl_seq,
        ref_flanks,
        idl_flanks,
        ref_count,
        alt_count,
        is_multiallelic,
        is_near_boundary,
        is_bidirectional,
        is_uniq_mapped,
        non_idl_flanks,
    )
