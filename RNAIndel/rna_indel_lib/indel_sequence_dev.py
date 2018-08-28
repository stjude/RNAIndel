#!/usr/bin/env python3

import pysam
import numpy as np
import rna_indel_lib.sequence_properties as sp
from rna_indel_lib.most_common import most_common


class Indel(object):
   """Represents indel by chr, pos, ins/del, seq

    Attributes:
        chr (str): chromosomes chr1-11, chrX, or chrY
        pos (int): 1-based position specifying the first base 
                   affected by the indel
        idl_type (int): 1 for insertion, 0 for deletion
        idl_seq (str): inserted or deleted sequence
   """
   def __init__(self, chr, pos, idl_type, idl_seq):
       self.chr = chr
       self.pos = pos
       self.idl_type = idl_type
       self.idl_seq = idl_seq
   
   @property
   def chr(self):
       return self.__chr

   @chr.setter
   def chr(self, chr):
       if not chr.startswith('chr'):
           self.__chr = 'chr' + chr
       else:
           self.__chr = chr
   
  
class SequenceWithIndel(Indel):
    """ Represents indel with its flanking sequence
    
    Attributes:
        lt_seq (str): 5' flanking seq
        rt_seq (str): 3' flanking seq
    """ 
  
    def __init__(self, chr, pos, idl_type, lt_seq, idl_seq, rt_seq):
        Indel.__init__(self, chr, pos, idl_type, idl_seq)
        self.lt_seq = lt_seq  
        self.rt_seq = rt_seq 


    def gc(self, n):
        """GC content
         
        Args:
            n (int)
        Returns:
            gc (float)
        """       
        if self.idl_type == 1:
            seq = self.lt_seq[-n:] + self.rt_seq[:n]
        # include the deleted sequence to recover the original sequence
        else:
            seq = self.lt_seq[-n:] + self.idl_seq + self.rt_seq[:n]
        
        return sp.gc(seq)

    
    def lc(self, n):
        """Average linguistic complexity in n-nt flanks
        
        Args:
            n (int): window in flanking seq
        Returns:
            lc (float or None): None if linguistic complexity is not defined
                                (ie., seq contains 'N')
        """
        lt_seq_n = self.lt_seq[-n:]
        rt_seq_n = self.rt_seq[:n]
        
        try:
            lt_lc = sp.linguistic_complexity(lt_seq_n)
            rt_lc = sp.linguistic_complexity(rt_seq_n)
            lc = (lt_lc + rt_lc) / 2
        except:
            lc = None

        return lc


    def local_lc(self, n):
        """The smaller linguistic complexity of n-nt flanks
        
        Args:
            n (int): window in flanking seq
        Returns:
            local_lc (float or None): None if linguistic complexity is not defined
        """
        lt_seq_n = self.lt_seq[-n:]
        rt_seq_n = self.rt_seq[:n]

        try:
            lt_lc = sp.linguistic_complexity(lt_seq_n)
            rt_lc = sp.linguistic_complexity(rt_seq_n)
            local_lc = min(lt_lc, rt_lc)
        except:
            local_lc = None

        return local_lc


    def strength(self, n):
        """DNA strength
        
        Args:
            n (int)
        Returns:
            dna_strength (float)
        """
        if self.idl_type == 1:
            seq = self.lt_seq[-n:] + self.rt_seq[:n]
        # include the deleted sequence to recover the original sequence
        else:
            seq = self.lt_seq[-n:] + self.idl_seq + self.rt_seq[:n]
        
        return sp.dna_strength(seq)


    def repeat(self):
        """Repeat
        
        Args:
            None
        Returns:
            repeat (int)
        """

        return sp.repeat(self.idl_type, self.lt_seq, self.idl_seq, self.rt_seq)


    def dissimilarity(self):
        """Dissimilarity

        Args:
            None
        Returns:
            dissimilarity (float)
        """
       
        return sp.dissimilarity(self.lt_seq, self.idl_seq, self.rt_seq)


class PileupWithIndelNotFound(Indel):
    """Represents indels not aligned/found as specified by caller
    """
    
    def __init__(self, chr, pos, idl_type, idl_seq):
        Indel.__init__(self, chr, pos, idl_type, idl_seq)
        self.ref_count = None
        self.alt_count = None
        self.is_multiallelic = None
        self.is_near_boundary = None
        self.is_bidirectional = None
        self.is_uniq_mapped = None

    def repeat(self):
        return None
        
    def local_gc(self):
        return None
     
    def local_lc(self):
        return None

    def local_strength(self):
        return None

    def dissimilarity(self):
        return None

    def indel_complexity(self):
        return None


class PileupWithIndel(Indel):
    """Represents pileup with indel 

    Attributes
        ref_flanks (list): list of [5' reference flank, 3' reference flank]
        idl_flanks (list): list of [5' indel_read flank, 3' indel_read flank]
        ref_count (int)
        alt_count (int)
        is_multiallelic (bool): 1 if true 0 otherwise
        is_near_boundary (bool) 1 if true 0 otherwise
        is_bidirectional (bool): 1 if true 0 otherwise
        is_uniq_mapped (bool): 1 if true 0 otherwise
        non_idl_flanks (list): list of [5' non_indel_read flank, 3' non_indel_read_flank]
    """

    def __init__(self, chr, pos, idl_type, idl_seq,
                 ref_flanks, idl_flanks,
                 ref_count, alt_count,
                 is_multiallelic, is_near_boundary,
                 is_bidirectional, is_uniq_mapped,
                 non_idl_flanks):
       Indel.__init__(self, chr, pos, idl_type, idl_seq)
       self.ref_flanks = ref_flanks
       self.idl_flanks = idl_flanks
       self.ref_count = ref_count
       self.alt_count = alt_count
       self.is_multiallelic = is_multiallelic
       self.is_near_boundary = is_near_boundary
       self.is_bidirectional = is_bidirectional
       self.is_uniq_mapped = is_uniq_mapped
       self.non_idl_flanks = non_idl_flanks 
    

    def generate_ref_reads(self):
        """Generates reference read
        
        Args: 
            None
        Returns:
            SequenceWithIndel (obj): representing reference sequence 
                                     with indel and splicing (if spliced)
        """
        ref_reads = [SequenceWithIndel(self.chr, self.pos, self.idl_type,\
                                       flank[0], self.idl_seq, flank[1]) \
                                       for flank in self.ref_flanks]
        return ref_reads


    def generate_indel_reads(self):
        """Generates indel read
        
        Args:
            None
        Returns:
            SequenceWithIndel (obj): representing indel read as aligned in bam
        """     
        indel_reads = [SequenceWithIndel(self.chr, self.pos, self.idl_type,\
                                         flank[0], self.idl_seq, flank[1])\
                                         for flank in self.idl_flanks]  
        return indel_reads
    

    def generate_non_indel_reads(self):
        """Generates non-indel read
        
        Args:
            None
        Returns:
            SequenceWithIndel (obj): representing non-indel reads as aligned in bam
                                     this read may contain polymorphisms.
        """
        non_indel_reads = [SequenceWithIndel(self.chr, self.pos, self.idl_type,\
                                             flank[0], self.idl_seq, flank[1])\
                                             for flank in self.non_idl_flanks]

        return non_indel_reads


    def repeat(self):
        """Most freqent number of repeats
        """
        repeats = []
        for indel in self.generate_indel_reads():
            repeat = indel.repeat()
            repeats.append(repeat)
        
        return most_common(repeats)


    def local_gc(self, n):
        """
        """
        local_vals = []
        for indel in self.generate_indel_reads():
            if len(indel.lt_seq) >= n and len(indel.rt_seq) >= n:
                local_vals.append(indel.gc(n))
            else:
                pass

        return np.mean(local_vals)
                 

    def local_lc(self, n):
        """Average local linguistic complexity values over indel reads
        
        Args:
            n (int): window in flanking seq
                     to make sure the flanking seq is n-nt or longer
        Returns:
            mean local lc (float): averaged over indel read found in the pileup
        """
        local_vals = []
        for indel in self.generate_indel_reads():
            if len(indel.lt_seq) >= n and len(indel.rt_seq) >= n:
                local_vals.append(indel.local_lc(n))
            else:
                pass
        
        return np.mean(local_vals)


    def local_strength(self, n):
        """
        """
        local_strengths = []
        for indel in self.generate_indel_reads():
            if len(indel.lt_seq) >= n and len(indel.rt_seq) >= n:
                local_strengths.append(indel.strength(n))
            else:
                pass
        
        return np.mean(local_strengths)

   
    def dissimilarity(self):
        """
        """
        dissimilarities = []
        idl_size = len(self.idl_seq)
        for indel in self.generate_indel_reads():
            if len(indel.lt_seq) >= idl_size and len(indel.rt_seq) >= idl_size:
                dissimilarities.append(indel.dissimilarity())
            else:
                pass
        
        return np.mean(dissimilarities)

    
    def indel_complexity(self, n):
        """
        """
        complexities = []
        indel_reads = self.generate_indel_reads()
        ref_reads = self.generate_ref_reads()
        for idl, ref in zip(indel_reads, ref_reads):
            if len(idl.lt_seq) >= n and len(idl.rt_seq) >= n:
                lt = idl.lt_seq[-n:]
                lt_ref = ref.lt_seq[-n:]
                lt_edit_dist = sp.editdistance(lt, lt_ref)

                rt = idl.rt_seq[:n]
                rt_ref = ref.rt_seq[:n]
                rt_edit_dist = sp.editdistance(rt, rt_ref)

                complexity = lt_edit_dist + rt_edit_dist
                        
                complexities.append(complexity)
            else:
                pass
        
        
        if complexities == []:
            return 0

        raw_value = min(complexities)
        if raw_value == 0:
            return 0

        # raw_value > 0, check for SNP-induced compleixty
        complexities = []
        indel_reads = self.generate_indel_reads()
        non_reads = self.generate_non_indel_reads()
        for idl in indel_reads:
            if len(idl.lt_seq) >= n and len(idl.rt_seq) >= n:
                for non in non_reads:
                    if len(non.lt_seq) >= n and len(non.rt_seq) >= n:
                        lt = idl.lt_seq[-n:]
                        lt_non = non.lt_seq[-n:]
                        lt_edit_dist = sp.editdistance(lt, lt_non)

                        rt = idl.rt_seq[:n]
                        rt_non = non.rt_seq[:n]
                        rt_edit_dist = sp.editdistance(rt, rt_non)

                        complexity = lt_edit_dist + rt_edit_dist
                        complexities.append(complexity)
        if complexities == []:
            return raw_value
        else:
            refined_value = min(complexities)
            return min(raw_value, refined_value)


class CodingSequenceWithIndel(SequenceWithIndel):
    def __init__(self, chr, pos, idl_type, lt_seq, idl_seq, rt_seq,
                 strand, accession, gene_symbol,
                 exon, exon_start, exon_end, last_exon,
                 cds_start,
                 prev_exon_start, prev_exon_end,
                 next_exon_start, next_exon_end):
        SequenceWithIndel.__init__(self, chr, pos,
                                   idl_type, lt_seq, idl_seq, rt_seq)

        self.strand = strand
        self.accession = accession
        self.gene_symbol = gene_symbol
        self.exon = exon
        self.exon_start = exon_start
        self.exon_end = exon_end
        self.last_exon = last_exon
        self.cds_start = cds_start
        self.prev_exon_start = prev_exon_start
        self.prev_exon_end = prev_exon_end
        self.next_exon_start = next_exon_start
        self.next_exon_end = next_exon_end
    
    def is_nmd_insensitive(self):
        is_insensitive = 0
        
        if self.exon == 1 or self.exon == self.last_exon:
            is_insensitive = 1

        return is_insensitive

    def effect(self):
        if self.strand == '+':
            if self.exon_start <= self.pos <= self.exon_end:
                return self.exonic_on_pos_strand()
            elif 0 < self.exon_start - self.pos <= 2 or \
            0 < self.pos - self.exon_end <= 2:
                return self.splice_site_on_pos_strand()
            elif 2 < self.exon_start - self.pos <= 11 or \
            2 < self.pos - self.exon_end <= 11:
                return self.splice_region_on_pos_strand()
            else:
                pass
        else:
            if self.exon_start <= self.pos <= self.exon_end:
                return self.exonic_on_neg_strand()
            elif  0 < self.exon_start - self.pos <= 2 or \
            0 < self.pos - self.exon_end <= 2:
                return self.splice_site_on_neg_strand()
            elif 2 < self.exon_start - self.pos <= 11 or \
            2 < self.pos - self.exon_end <= 11:
                return self.splice_region_on_neg_strand()
            else:
                pass
    

    def cds_pos_in_exonic_indels(self):
        """the position of the first coding sequence (CDS) base affected
           by the indel. 'cds_pos' - 1 gives the last unaffected base position.
           
           As name suggests, this method is for indels within exon.

           Example:          1234567890123
                      CDS  : ATGCTACGACTGA
                       del : ATGCTA-GACTGA  -> cds_pos = 7
                             
                             123456 7890123
                      CDS  : ATGCTA*CGACTGA             
                       ins : ATGCTATCGACTGA  -> cds_pos = 7
               Note that the sequences are unaffected upto 
               first 6 (i.e., cds_pos-1) bases in both cases.  
        """
        # insertion/deletion on positive strand
        if self.strand == '+':
                cds_pos = self.cds_start + self.pos - self.exon_start
        else: 
            # insertion on negative strand
            if self.idl_type == 1:
                cds_pos = self.cds_start + self.exon_end - (self.pos - 1)
            # deletion on negative strand
            else:
                cds_pos = self.cds_start + self.exon_end - self.pos\
                          - (len(self.idl_seq) - 1)
        return cds_pos


    def exonic_on_pos_strand(self):
        """return like aapos_effect
        """
        
        # insertion at 5'exon_start
        if self.idl_type == 1 and self.pos == self.exon_start:
            cds_pos = self.cds_start - 1
            codon_pos = int(cds_pos/3) + 1        
            if len(self.idl_seq) > 1 and self.idl_seq[-2:] == 'AG':
                return codon_pos, 'splicePreserving'
            else:
                return codon_pos, 'spliceTruncating'
            
        # indels within exon
        else:
            cds_pos = self.cds_pos_in_exonic_indels()
          
            frame = (cds_pos - 1) % 3
            if frame == 2:
                codon_pos = int(cds_pos/3)
            else:
                codon_pos = int(cds_pos/3) + 1

            # insertion
            if self.idl_type == 1:
                if frame == 0:
                    seq = self.idl_seq + self.rt_seq[:2]
                elif frame == 1:
                    seq = self.lt_seq[-1:] + self.idl_seq + self.rt_seq[:1] 
                else:
                    seq = self.lt_seq[-2:] + self.idl_seq + self.rt_seq[:3]   
            # deletion
            else:
                if frame == 0:
                    seq = self.rt_seq[:3]
                elif frame == 1:
                    seq = self.lt_seq[-1:] + self.rt_seq[:2]
                else:
                    seq = self.lt_seq[-2:] + self.rt_seq[:1]
            # check for stop codon
            if sp.exists_stop_codon(self.strand, seq):
                return codon_pos, 'nonsenseTruncating'
            else:
                if len(self.idl_seq) % 3 == 0 and self.idl_type == 1:
                    return codon_pos, 'inframeIns'
                elif len(self.idl_seq) % 3 == 0 and self.idl_type == 0:    
                    return codon_pos, 'inframeDel'
                else:
                    return codon_pos, 'frameshiftTruncating'

    
    def splice_site_on_pos_strand(self):
        
        # splicing motif + at least 1 base
        # ith_exon_end + GT + (at least one) + AG + (i+1)th_exon_start
        min_motif_len = 6

        # 5'splice
        if self.exon_start > self.pos:
            cds_pos = self.cds_start - 1
            codon_pos = int(cds_pos/3) + 1

            if self.exon_start - self.prev_exon_end <= min_motif_len:
                return codon_pos, 'spliceShortIntron'
            
            else:
                # insertion at 1-nt upstream of exon start
                if self.idl_type == 1 and self.exon_start - self.pos == 1:
                    if self.idl_seq[-1] == 'A':
                        return codon_pos, 'splicePreserving'
                    else:
                        return codon_pos, 'spliceTruncating'
                
                # insertion at 2-nt upstream of exon start
                elif self.idl_type == 1 and self.exon_start - self.pos == 2:
                    return codon_pos, 'spliceRegion'
                
                # deletion at 1-nt upstream of exon start
                elif self.idl_type == 0 and self.exon_start - self.pos == 1:
                    return codon_pos, 'spliceTruncating'
                
                # deletion at 2-nt upstream of exon start
                elif self.idl_type == 0 and self.exon_start - self.pos == 2:
                    if len(self.idl_seq) == 1 and self.lt_seq[-1] == 'A':
                        return codon_pos, 'splicePreserving'
                    elif  len(self.idl_seq) == 2 and self.lt_seq[-2:] == 'AG':
                        return codon_pos, 'splicePreserving'
                    else:
                        return codon_pos, 'spliceTruncating'
                else:
                    pass
        
        # 3'splice
        else:
            cds_pos = self.cds_start + self.exon_end - self.exon_start
            codon_pos = int(cds_pos/3) + 1
            
            if self.next_exon_start - self.exon_end <= min_motif_len:
                return codon_pos, 'spliceShortIntron'

            else:
                # insertion 1-nt downstream of exon end
                if self.idl_type == 1 and self.pos - self.exon_end == 1:
                    if len(self.idl_seq) > 1 and self.idl_seq[:2] == 'GT':
                        return codon_pos, 'splicePreserving'
                    else:
                        return codon_pos, 'spliceTruncating'
                
                # insertion 2-nt downstream of exon end
                elif self.idl_type == 1 and self.pos - self.exon_end == 2:
                    if self.idl_seq[0] == 'T':
                        return codon_pos, 'splicePreserving'
                    else:
                        return codon_pos, 'spliceTruncating'
                
                # deletion 1-nt downstream of exon end
                elif self.idl_type == 0 and self.pos - self.exon_end == 1:
                    if len(self.idl_seq) > 1 and self.rt_seq[:2] == 'GT':
                        return codon_pos, 'splicePreserving'
                    else:
                        return codon_pos, 'spliceTruncating'
                
                # deletion 2-nt downstream of exon end
                elif self.idl_type == 0 and self.pos - self.exon_end == 2:
                    if self.rt_seq[0] == 'T':
                        return codon_pos, 'splicePreserving'
                    else:
                         return codon_pos, 'spliceTruncating'
                else:
                    pass


    def splice_region_on_pos_strand(self):
        # 5'splice region
        if self.exon_start > self.pos:
            cds_pos = self.cds_start - 1 
            codon_pos = int(cds_pos/3) + 1
            return codon_pos, 'spliceRegion'
        
        # 3'splice region
        else:
            cds_pos = self.cds_start + self.exon_end - self.exon_start
            codon_pos = int(cds_pos/3) + 1
            return codon_pos, 'spliceRegion'
            

    def exonic_on_neg_strand(self):
        
        # insertion at 3'exon_start
        if self.idl_type == 1 and self.pos == self.exon_start:
            cds_pos = self.cds_start + self.exon_end - self.pos
            codon_pos = int(cds_pos/3) + 1
            if len(self.idl_seq) > 1 and self.idl_seq[-2:] == 'AC':
                return codon_pos, 'splicePreserving'
            else:
                return codon_pos, 'spliceTruncating'
        
        # indels within exon
        else:
            cds_pos = self.cds_pos_in_exonic_indels()
            
            frame = (cds_pos - 1) % 3
            if frame == 2:
                codon_pos = int(cds_pos/3)
            else:
                codon_pos = int(cds_pos/3) + 1
                 
            # insertion
            if self.idl_type == 1:
                if frame == 0:
                    seq = self.lt_seq[-2:] + self.idl_seq
                elif frame == 1:
                    seq = self.lt_seq[-1:] + self.idl_seq + self.rt_seq[:1]
                else:
                    seq = self.lt_seq[-3:] + self.idl_seq + self.rt_seq[:2]
            # deletion
            else:
                if frame == 0:
                    seq = self.lt_seq[-3:]
                elif frame == 1:
                    seq = self.lt_seq[-2:] + self.rt_seq[:1]
                else:
                    seq = self.lt_seq[-1:] + self.rt_seq[:2]
            # check for stop codon
            if sp.exists_stop_codon(self.strand, seq):
                return codon_pos, 'nonsenseTruncating'
            else:
                if len(self.idl_seq) % 3 == 0 and self.idl_type == 1:
                    return codon_pos, 'inframeIns'
                elif len(self.idl_seq) % 3 == 0 and self.idl_type == 0:
                    return codon_pos, 'inframeDel'
                else:              
                    return codon_pos, 'frameshiftTruncating'     
            
   
    def splice_site_on_neg_strand(self):
        # splicing motif + at least 1 base
        # ith_exon_end + GT + (at least one) + AG + (i+1)th_exon_start
        min_motif_len = 6

        # 5'splice
        if self.pos > self.exon_end:
            cds_pos = self.cds_start - 1
            codon_pos = int(cds_pos/3) + 1

            if self.prev_exon_start - self.exon_end <= min_motif_len:
                return codon_pos, 'spliceShortIntron'
            
            else:
                # insertion at 1-nt downstream of exon end
                if self.idl_type == 1 and self.pos - self.exon_end == 1:
                    if len(self.idl_seq) > 1 and self.idl_seq[:2] == 'CT':
                        return codon_pos, 'splicePreserving'
                    else:
                        return codon_pos, 'spliceTruncating'
                
                # insertion at 2-nt downstream of exon end
                elif self.idl_type == 1 and self.pos - self.exon_end == 2:
                    if self.idl_seq[0] == 'T':
                        return codon_pos, 'splicePreserving'
                    else:
                        return codon_pos, 'spliceTruncating'
                
                # deletion at 1-nt downstream of exon end
                elif self.idl_type == 0 and self.pos - self.exon_end == 1:
                    if len(self.idl_seq) > 1 and self.rt_seq[:2] == 'CT':
                        return codon_pos, 'splicePreserving'
                    else:
                        return codon_pos, 'spliceTruncating'
                
                # deletion at 2-nt downstream of exon end
                elif self.idl_type == 0 and self.pos - self.exon_end == 2:
                    if self.rt_seq[0] == 'T':
                        return codon_pos, 'splicePreserving'
                    else:
                        return codon_pos, 'spliceTruncating'
                else:
                    pass
        
        # 3'splice
        else:
            cds_pos = self.cds_start + self.exon_end - self.exon_start
            codon_pos = int(cds_pos/3) + 1
            
            if self.exon_start - self.next_exon_end <= min_motif_len:
                return codon_pos, 'spliceShortIntron'

            else:
                # insertion 1-nt upstream of exon start
                if self.idl_type == 1 and self.exon_start - self.pos == 1:
                    if self.idl_seq[-1] == 'A':
                        return codon_pos, 'splicePreserving'
                    else:
                        return codon_pos, 'spliceTruncating'
                
                # insertion 2-nt upstream of exon start
                elif self.idl_type == 1 and self.exon_start - self.pos == 2:
                    return codon_pos, 'spliceRegion'
                
                # deletion 1-nt upstream of exon start
                elif self.idl_type == 0 and self.exon_start - self.pos == 1:
                    return codon_pos, 'spliceTruncating'
                
                # deletion 2-nt upstream of exon start
                elif self.idl_type == 0 and self.exon_start - self.pos == 2:
                    if len(self.idl_seq) == 2 and self.lt_seq[-2:] == 'AC':
                        return codon_pos, 'splicePreserving'
                    else:
                        return codon_pos, 'spliceTruncating'
                else:
                    pass
       
        
    def splice_region_on_neg_strand(self):
        # 5'splice region
        if self.pos > self.exon_end:
            cds_pos = self.cds_start - 1 
            codon_pos = int(cds_pos/3) + 1
            return codon_pos, 'spliceRegion'
        
        # 3'splice region
        else:
            cds_pos = self.cds_start + self.exon_end - self.exon_start
            codon_pos = int(cds_pos/3) + 1
            return codon_pos, 'spliceRegion'


