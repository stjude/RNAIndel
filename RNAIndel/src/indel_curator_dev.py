#!/usr/bin/env python3

import re
import pysam
import random
import logging
import numpy as np
from most_common import most_common
from indel_sequence_dev import SequenceWithIndel
from indel_sequence_dev import PileupWithIndel
from indel_sequence_dev import PileupWithIndelNotFound


random.seed(123)
logger = logging.getLogger(__name__)


def curate_indel_in_genome(fasta, chr, pos, idl_type, idl_seq):
    """Gerenates an indel object with reference genome flanking sequences.
       Splicing will NOT be considered as it is genome.
       
       Example:
           The read below is spliced (>) and has an 'A' deletion (*)
           at chr1:13-13.
           
           Chr: 1
           Pos:       12345678901234567890123456789012345
           Reference: CGTATGATGTTCAGGTATGCGTATATAGAAAATCGA
                                  *
           Read:      CGTATGATGTTC-G>>>>>>>>>>>>>>AAAATCGA

           For window=10, this function returns an indel object represented by:
               chr: 1
               pos: 13
               idl_type: 0 (0 for del, 1 for ins)
               lt_seq: 'TATGATGTTC' (left flanking reference with len=10)
               idl_seq: 'A'
               rt_seq: 'GGTATGCGTA' (right flanking reference with len=10)    

       Args:
           fasta: str, complete path to .fa
           chr: str, chr1-22, chrX, chrY
           pos: int, 1-based position of indel
           idl_type: boolean, 1 for insertion 0 for deletion
           idl_seq: str, non-empty string of indel sequence
           (window: int, range for flank region # hard-coded as 50nt)

       Returns:
           GenomeWithIndel: class

       Raises:
           AssertionError: if idl_type is not 1 or 0
           AssertionError: if idl_seq is empty
    """
    window = 50
    
    # left flank seq
    start, end  = pos - window, pos - 1
    # retrieve left flank sequence in FASTA format 
    lt_fasta = pysam.faidx(fasta, chr+':'+str(start)+'-'+str(end))
    # extract the sequence string
    lt_seq = lt_fasta.split('\n')[1]
     
    # right flank seq
    # for insertion
    if idl_type == 1:  
        start, end = pos, pos - 1 + window

        # retrieve right flank sequence in FASTA format
        rt_fasta = pysam.faidx(fasta, chr+':'+str(start)+'-'+str(end))
        # extract the seqence string
        rt_seq  = rt_fasta.split('\n')[1]
    # for deletion
    else:
        size = len(idl_seq)  # needed to adjust by the indel length
        start, end = pos + size, pos - 1 + size + window  
        
        # retrieve right flank sequence in FASTA format
        rt_fasta = pysam.faidx(fasta, chr+':'+str(start)+'-'+str(end))
        # extract the seqence string
        rt_seq = rt_fasta.split('\n')[1]
    
    return SequenceWithIndel(chr, pos, idl_type, lt_seq, idl_seq, rt_seq)


def is_close_to_exon_boundary(cigar_lst, idx):
    """Checks if indel is close to the exon boundary.
       'close' means within 2-nt to the boundary.

       Example:
         The distance between indel and the splice site
         is 1-nt. This is 'close'
         
         Reference: CGTATGATGTTCAGGTATGCGTATATAGAAAATCGA
           
           Read:       ATGACGTTC-G>>>>>>>>>>>>>>AAAATCGA
           CIGAR: 9M1D1M14N8M         
           CIGAR_LST: ['9M, '1D', '1M', '14N', '8M']
                 IDX:  1 
       Args:
           cigar_lst: list, see the example above
           idx: int, the index of cigar_lst to specify 
                the indel of interest.
       Returns:
           is_close: boolean, 1 for yes 0 for no.   
    """
    is_close = 0

    if idx >= 2 and 'N' in cigar_lst[idx-2]:
        dist_to_lt_boundary = int(cigar_lst[idx-1].replace('M', ''))
        if dist_to_lt_boundary <= 2:
            is_close = 1
                        
    elif (idx+2) <= len(cigar_lst)-1 and 'N' in cigar_lst[idx+2]:
        dist_to_rt_boundary = int(cigar_lst[idx+1].replace('M', ''))      
        if dist_to_rt_boundary <= 2:
            is_close = 1
    else:
        pass

    return is_close


def extract_all_valid_reads(bam_data, chr, pos):
    """Extracts reads that are
        1. non-duplicate
        2. primary alignment
        3. covering the locus of interest (non-skipping)

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
    all_reads = bam_data.fetch(chr, pos, pos+1, until_eof=True)
    
    valid_reads = []
    for read in all_reads:
        # excludes duplicate or non-primary alignments
        if read.is_duplicate == False and read.is_secondary == False:
            blocks = read.get_blocks()  
            for block in blocks:                 
                # excludes skipping reads
                if block[0] <= pos <= block[1]: 
                    valid_reads.append(read)

    return valid_reads


def curate_indel_in_pileup(bam_data, chr, pos, idl_type, idl_seq):
    """Generates an object describing what indel looks like
       in the pileup view.
    
       Example:
           In the alignment below, the reads are spliced (>).
           
           Read_1 and Read_3 are indel-reads ('A' deletion (*)).
           Read_2 is a non-indel read.
           Read_4 skips this locus.
           The T>C SNP is found in Read_1-3. (^).
           

           Chr: 1
           Pos:       12345678901234567890123456789012345
           Reference: CGTATGATGTTCAGGTATGCGTATATAGAAAATCGA
                             ^    *              
           Read_1        ATGACGTTC-G>>>>>>>>>>>>>>AAAATCGA
           Read_2      GTATGACGTTCAG>>>>>>>>>>>>>>AAAAT
           Read_3     CGTATGACGTTC-G>>>>>>>>>>>>>>AAA
           Read_4     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   

           This function does:
               1. extraction of spliced reference sequence flanking the indel

                  For Read_1, 'ATGATGTTC'(lt), 'GAAAATCGA'(rt)            
                  For Read_3, 'CGTATGATGTTC'(lt), 'GAAA'(rt)

                  No individual gene variations (e.g. SNP) contained. 
               
               2. extraction of indel read flanking seqence
                   
                  For Read_1, 'ATGACGTTC'(lt), 'GAAAATCGA'(rt)
                  For Read_3, 'CGTATGACGTTC'(lt), 'GAAA'(rt)

                  May contain individual gene variations.

               3. count of non-indel fragments (NOT reads).
               4. count of indel fragments (NOT reads).
                   
                  Non-indel count: 1 (Read_2)
                  Indel count: 1 if Read_1 and Read_3 are mate
                               2 otherwise 
                   
                  Fragment count:

                    Fw ----------->   
                       GACGTTC-GAAA
                            <---------- Rv
                            
                  The indel is overlapped by Fw and Rv reads (mate reads).
                  In this case, the Fw and Rv reads belong to the same
                  fragment. Counted as 1. 
                
               5. check if indel is variable in length
                    
                  The insertions below are variable in length.
                  Read_a, Read_c, Read_d: 'A' inserted
                                  Read_b: 'AA' inserted 
                    
                  Reference: CTGGAAAA**TCACAA
                    
                  Read_a     CTGGAAAAA*TCACAA  
                  Read_b     CTGGAAAAAATCACAA
                  Read_c     CTGGAAAAA*TCACAA
                  Read_d     CTGGAAAAA*TCACAA

               6. check if indel is within exon but close to the exon boundary.
                  'close' means within 2-nt from the boundary.
                  In Example, the boundary is chr1:14-14
                              the indel is chr1:13-13
                  This case is 'close'.                 
                
               7. check if indel is bidirectionally supported.
                  'bidirectionally supported' if Fw and Rv reads support the 
                  indel.
               
               8. check if indel reads are uniquely mapped.
                  This is based on mapping quality score by STAR.
                  Collect scores for all indel reads.
                  Uniquely mapped if the most frequent score is 255.
               
               9. random sampling non-indel reads (up to 10 reads) 
                  and extraction of the seence flnaking the indel
                  
                  For Read_2, 'GTATGACGTTC'(lt), 'GAAAA'(rt)   
                   
       Args:    
           bam: str, path to .bam
           chr: str, chr1-22, chrX, chrY
           pos: int, 1-based position of indel on the reference
           idl_type: boolean, 1 for insertion 0 for deletion
           idl_seq: str, non-empty string of indel sequence

       Returns:
           PileupWithIndel object: if indels found as specified with 
                                   chr, pos, idl_type and idl_seq 
           PileupWithIndelNotFound object: otherwise

       Raises:
           AssertionError: if idl_type is not 1 or 0
           AssertionError: if idl_seq is an empty string
    """
    # convert to 0-based coordinate
    pos = pos - 1
    # convert indel type to CIGAR token
    if idl_type == 1:
        idl_type = 'I'
    else:
        idl_type = 'D'
    # compile regex to parse CIGAR string 
    cigar_ptn = re.compile(r'[0-9]+[MIDNSHPX=]')
    
    size = len(idl_seq) 
   
    # extract all reads covering the locus of interest
    all_reads = extract_all_valid_reads(bam_data, chr, pos)
    
    read_names = []
    idl_read_names = [] 
    read_flanks = []  
    ref_flanks =[]  
    is_multiallelic = 0
    exon_boundary = []
    map_qual = []
    bidirectional = []
    for read in all_reads:
        name = read.query_name
        read_names.append(name)  # for total coverage by uniq read
                
        ref_start = read.reference_start
        cigarstring = read.cigarstring
        read_seq = read.query_sequence
        ref_seq = read.get_reference_sequence() 
        
        if idl_type in cigarstring and read.is_duplicate == False:
            cigars = cigar_ptn.findall(cigarstring)
            
            # adjust ref_start if the read starts with soft-clipping
            adjust = 0 
            if cigars[0].endswith('S'):
                adjust = int(cigars[0].replace('S', ''))
                ref_start = ref_start - adjust
            
            # Pos on Genome:       1 2 3 4 5 6 7 8 9 
            # read_seq                 a a A C - C G (first 2 small 'a' soft-clipped)
            # ref_seq                  A A A C C C G     
            # Pos on read_seq (i)      0 1 2 3   4 5       
            # Pos on ref_seq (j)      -2-1 0 1 2 3 4 (count from the first mapped base) 
            # ref_start = 5 (before adjust)
            # ref_start = 3 (after adjust)

            i = 0  # position on read_seq
            j = 0 - adjust  # position on ref_seq
            ref_pos = ref_start   # postion on genome
            for idx, cigar in enumerate(cigars):
                val = int(cigar.replace(cigar[-1], ''))  # cigar value  
                ope = cigar[-1]  # cigar operation
                
                if ref_pos == pos:
                    
                    # check for variable indels
                    if val != size and ope == 'I' and idl_type == 'I':
                        is_multiallelic = 1
                    if val != size and ope == 'D' and idl_type == 'D':
                        is_multiallelic = 1
                    
                    # extract flanks for indel reads 
                    if val == size and ope == 'I' and idl_type == 'I':
                        lt_read = read_seq[:i].upper()
                        lt_ref = ref_seq[:j].upper()
                        rt_read = read_seq[i+val:].upper()
                        rt_ref = ref_seq[j:].upper()
                        
                        read_flanks.append([lt_read, rt_read])
                        
                        ref_flanks.append([lt_ref, rt_ref])
                        
                        idl_read_names.append(name)
                        
                        exon_boundary.append(is_close_to_exon_boundary(cigars, idx))
                        
                        map_qual.append(read.mapping_quality)
                         
                        bidirectional.append(read.is_reverse)

                    if val == size and ope == 'D' and idl_type == 'D':
                        lt_read = read_seq[:i].upper()
                        lt_ref = ref_seq[:j].upper()
                        rt_read = read_seq[i:].upper()
                        rt_ref =ref_seq[j+val:].upper()
                        
                        read_flanks.append([lt_read, rt_read])
                        
                        ref_flanks.append([lt_ref, rt_ref])
                        
                        idl_read_names.append(name)
                        
                        exon_boundary.append(is_close_to_exon_boundary(cigars, idx))

                        map_qual.append(read.mapping_quality)

                        bidirectional.append(read.is_reverse)
                
                if ope == 'I':   # if ins, no move on reference
                    ref_pos = ref_pos
                else:
                    ref_pos = ref_pos + val

                if ope == 'N':  # if spliced, no move
                    i = i
                    j = j  
                elif ope == 'D':  # if del, move on ref_read
                    i = i
                    j = j + val 
                elif ope == 'I':
                    i = i + val  # if ins, move on read
                    j = j
                else:
                    i = i + val
                    j = j + val 

    # when indels not found as specified...
    if idl_read_names == []:
        
        notfound = 'NotFoundAsSpecified: '
        if idl_type == 'I':
            msg = notfound + chr + '|' + str(pos+1) + '|' + '-' + '|' + idl_seq
        else:
            msg = notfound + chr + '|' + str(pos+1) + '|' + idl_seq + '|' + '-'
        
        logging.info(msg)
        
        # retuns NotFound obj
        return PileupWithIndelNotFound(chr, pos, idl_type, idl_seq)
    
    # collect non-indel read names
    non_idl_read_names = list(set(read_names) - set(idl_read_names))
    non_idl_read_names.sort()
    
    # sample 10 non-indel reads if too many
    if len(non_idl_read_names) > 10:
        non_idl_read_names = random.sample(non_idl_read_names, 10)
  
    # make dict to access read seqence by name
    reads_by_name = {read.query_name: read for read in all_reads}
    
    # get flanking seq for non-indel reads
    non_idl_flanks = []
    for name in non_idl_read_names:
       read = reads_by_name[name]
         
       ref_start = read.reference_start
       sequence = read.query_sequence
       cigarstring = read.cigarstring
       cigars = cigar_ptn.findall(cigarstring)
       
       # adjust if the read starts with softclipping
       adjust = 0
       if cigars[0].endswith('S'):
           adjust = int(cigars[0].replace('S', '')) 
           ref_start = ref_start - adjust
        
       i = 0 
       ref_pos = ref_start
       for cigar in cigars:
           val = int(cigar.replace(cigar[-1], '')) 
           ope = cigar[-1]
           
           if ref_pos < pos:
               if ope == 'I':
                   ref_pos = ref_pos
               else:
                   ref_pos = ref_pos + val

               if ope == 'N':
                   i = i
               elif ope == 'D':
                   i = i
               elif ope == 'I':
                   i = i + val
               else:
                   i = i + val
           else:
               break
        
       diff = pos - ref_pos
       if idl_type == 'I':
           lt = sequence[:i+diff]
           rt = sequence[i+diff:]
           non_idl_flanks.append([lt, rt])
       else:
           lt = sequence[:i+diff]
           rt = sequence[i+diff+size:]
           non_idl_flanks.append([lt, rt])
    
    # when no no-indel reads found (i.e., indel VAF = 1.0)
    if non_idl_flanks == []:
        if len(ref_flanks) > 10:
            non_idl_flanks = random.sample(ref_flanks, 10)
        else:
            non_idl_flanks = ref_flanks
    
    # fragment count by unifiying the read name
    total_count, alt_count  = len(set(read_names)), len(set(idl_read_names))
    ref_count = total_count - alt_count
    
    # decide if the indel is close to exon boundary
    is_near_boundary = most_common(exon_boundary)

    # decide if uniquely mapped
    most_common_map_qual = most_common(map_qual)
    is_uniq_mapped = 0
    if most_common_map_qual == 255:
        is_uniq_mapped = 1
   
    # decide if bidirectionally supported
    is_bidirectional = 0
    if len(set(bidirectional)) == 2:
        is_bidirectional = 1

    # back to 1-based coordinate
    pos = pos + 1
    
    # convert idl_type 
    if idl_type == 'I':
        idl_type = 1
    else:
        idl_type = 0 

    return PileupWithIndel(chr, pos, idl_type, idl_seq,\
                           ref_flanks, read_flanks,\
                           ref_count, alt_count,\
                           is_multiallelic, is_near_boundary,\
                           is_bidirectional, is_uniq_mapped,\
                           non_idl_flanks) 
