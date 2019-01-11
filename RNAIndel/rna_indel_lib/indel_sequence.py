#!/usr/bin/env python3

from .most_common import most_common
from .sequence_properties import *


class Indel(object):
    """ Represents indel by chr, pos, ins/del, seq
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
    def ref(self):
        if self.idl_type == 1:
            return "-"
        else:
            return self.idl_seq

    @property
    def alt(self):
        if self.idl_type == 1:
            return self.idl_seq
        else:
            return "-"


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

        return gc(seq)

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
            lt_lc = linguistic_complexity(lt_seq_n)
            rt_lc = linguistic_complexity(rt_seq_n)
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
            lt_lc = linguistic_complexity(lt_seq_n)
            rt_lc = linguistic_complexity(rt_seq_n)
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

        return dna_strength(seq)

    def repeat(self):
        """Repeat
        
        Args:
            None
        Returns:
            repeat (int)
        """

        return repeat(self.idl_type, self.lt_seq, self.idl_seq, self.rt_seq)

    def dissimilarity(self):
        """Dissimilarity

        Args:
            None
        Returns:
            dissimilarity (float)
        """

        return dissimilarity(self.lt_seq, self.idl_seq, self.rt_seq)

    def __eq__(self, other):
        """Equality for equivalent indels

        Args: 
            self, other: SequenceWithIndel objects
        Returns:
            True/False (bool): True if they are equivalent
            
        Example:
            Pos       12345678
            Refernce  ATGATACC
            Indel 1   ATG--ACC
            Indel 2   ATGA--CC
            
            Indel 1 ('AT' del at 4) and Indel 2 ('TA' del at 5)
            are alternaive alignments of the same sequence.
            The objects created from Indel 1 and 2 should
            be identical.
        """
        idl1, idl2 = self, other

        # assume idl2 is on the left side
        if idl1.pos > idl2.pos:
            idl1, idl2 = idl2, idl1

        chr1 = idl1.chr
        chr2 = idl2.chr
        idl_type1 = idl1.idl_type
        idl_type2 = idl2.idl_type
        idl_seq1 = idl1.idl_seq
        idl_seq2 = idl2.idl_seq

        # rejects trivial cases
        if idl1.chr != idl2.chr:
            return False
        if idl_type1 != idl_type2:
            return False
        if len(idl_seq1) != len(idl_seq2):
            return False

        # here after the two indels are
        # of same type (ins or del) with
        # the same indel length, and
        # indel 2 pos >= indel 1 pos.

        n = len(idl_seq1)
        m = idl2.pos - idl1.pos

        # insertion cases
        if idl_type1 == 1:
            s = idl1.rt_seq[0:m]

            if m > n:
                if idl_seq1 == s[:n] and s[: (m - n)] == s[n:] and idl_seq2 == s[-n:]:
                    return True
                else:
                    return False

            elif m == n:
                if idl_seq1 == s == idl_seq2:
                    return True
                else:
                    return False

            elif m > 0 and m < n:
                if (
                    idl_seq1[:m] == s == idl_seq2[-m:]
                    and idl_seq1[-(n - m) :] == idl_seq2[: (n - m)]
                ):
                    return True
                else:
                    return False

            else:
                if idl_seq1 == idl_seq2:
                    return True
                else:
                    return False

        # deletion cases
        else:
            if m == 0:
                if idl_seq1 == idl_seq2:
                    return True
                else:
                    return False
            else:
                s = (idl_seq1 + idl1.rt_seq)[:m] + idl_seq2
                if s[:m] == s[n : (m + n)]:
                    return True
                else:
                    return False


class PileupWithIndelNotFound(Indel):
    """Represents indels not aligned/found as specified by caller

    Attributes:
        chr (str): as specified by caller
        pos (int): as specified by caller
        idl_type (int): as specified by caller
        idl_seq (str): as specified by caller
        
        other attributes are set None      
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

    For each read, SequenceWithIndel obj. is created
    and sequence properties are calculated. The summary
    of the properties of all read represents the pileup.

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

    def __init__(
        self,
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
    ):

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
        ref_reads = [
            SequenceWithIndel(
                self.chr, self.pos, self.idl_type, flank[0], self.idl_seq, flank[1]
            )
            for flank in self.ref_flanks
        ]

        return ref_reads

    def generate_indel_reads(self):
        """Generates indel read
        
        Args:
            None
        Returns:
            SequenceWithIndel (obj): representing indel read as aligned in bam
        """
        indel_reads = [
            SequenceWithIndel(
                self.chr, self.pos, self.idl_type, flank[0], self.idl_seq, flank[1]
            )
            for flank in self.idl_flanks
        ]

        return indel_reads

    def generate_non_indel_reads(self):
        """Generates non-indel read
        
        Args:
            None
        Returns:
            SequenceWithIndel (obj): representing non-indel reads as aligned in bam
                                     this read may contain polymorphisms.
        """
        non_indel_reads = [
            SequenceWithIndel(
                self.chr, self.pos, self.idl_type, flank[0], self.idl_seq, flank[1]
            )
            for flank in self.non_idl_flanks
        ]

        return non_indel_reads

    def repeat(self):
        """Most frequent number of repeats in pileup

        Args:
            None
        Returns:
            most frequent repeat number (int)
        """
        repeats = []
        for indel in self.generate_indel_reads():
            repeat = indel.repeat()
            repeats.append(repeat)

        return most_common(repeats)

    def local_gc(self, n):
        """Average GC content
        
        Args:
            n (int): length of flanking sequence considered
                     6 for RNA read
                     50 for reference genome
        Returns:
            average local GC content (float) 
        """
        local_vals = []
        for indel in self.generate_indel_reads():
            if len(indel.lt_seq) >= n and len(indel.rt_seq) >= n:
                local_vals.append(indel.gc(n))
            else:
                pass

        return np.mean(local_vals)

    def local_lc(self, n):
        """Average local Linguistic Complexity 
        
        Args:
            n (int): length of flanking sequence considered
                     6 for RNA read
                     50 for reference genome
        Returns:
            average local lc (float)
        """
        local_vals = []
        for indel in self.generate_indel_reads():
            if len(indel.lt_seq) >= n and len(indel.rt_seq) >= n:
                local_vals.append(indel.local_lc(n))
            else:
                pass

        return np.mean(local_vals)

    def local_strength(self, n):
        """Average local DNA-strength

        Args:
            n (int): length of flanking sequence considered
                     6 for RNA read
                     50 for reference genome
        Returns:
            average DNA-strength (float)  
        """
        local_strengths = []
        for indel in self.generate_indel_reads():
            if len(indel.lt_seq) >= n and len(indel.rt_seq) >= n:
                local_strengths.append(indel.strength(n))
            else:
                pass

        return np.mean(local_strengths)

    def dissimilarity(self):
        """Average Dissimilarity

        Args:
            None
        Returns:
            average Dissimilarity (float)
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
        """Mininum edit distance between indel and non-indel flanking sequences
        
        Args:
            n (int): length of flanking sequence considered
                     n = 6 for RNA read
                     indel complexity is not defined for reference.
        Returns:   
           indel complexity (int) 
        
        Example:           * ***      mismatch
                 CGTAGTAT  GAAGCAAAGT (non-indel_read 1)
                 CGTAGTATAGAAGCAAAAGT (indel_read 1)
                 CGTAGTATAGAAGCAAAAGT (indel_read 2)
                 CGTAGTAT  GAAGCAAAGT (non-indel_read_2)
                 
                 For indel_read_1 and non-indel_read_1, 

                 left edit_dist(TAGTAT, TAGTAT) = 0
                 right edit_dist(AAGCAA, GAAGCA) = 2
                
                 indel complexity = left edit_dist + right edit_dist
                                  = 2 
                                      
                 Calculate for indel_read_i nad non-indel_read_j 
                 and take mininum.
        """
        # first compare with reference
        # if indel complexity against ref is 0, return 0.
        # -> O(N) (N: num of indel reads)
        complexities = []
        indel_reads = self.generate_indel_reads()
        ref_reads = self.generate_ref_reads()
        for idl, ref in zip(indel_reads, ref_reads):
            if len(idl.lt_seq) >= n and len(idl.rt_seq) >= n:
                lt = idl.lt_seq[-n:]
                lt_ref = ref.lt_seq[-n:]
                lt_edit_dist = editdistance(lt, lt_ref)

                rt = idl.rt_seq[:n]
                rt_ref = ref.rt_seq[:n]
                rt_edit_dist = editdistance(rt, rt_ref)

                complexity = lt_edit_dist + rt_edit_dist

                complexities.append(complexity)
            else:
                pass

        if complexities == []:
            return 0

        indel_complexity_against_ref = min(complexities)
        if indel_complexity_against_ref == 0:
            return 0

        # indel_complexity_against_ref > 0, check for SNP-induced compleixty
        # -> O(NxM) (M: num of non-indel reads)
        complexities = []
        indel_reads = self.generate_indel_reads()
        non_reads = self.generate_non_indel_reads()
        for idl in indel_reads:
            if len(idl.lt_seq) >= n and len(idl.rt_seq) >= n:
                for non in non_reads:
                    if len(non.lt_seq) >= n and len(non.rt_seq) >= n:
                        lt = idl.lt_seq[-n:]
                        lt_non = non.lt_seq[-n:]
                        lt_edit_dist = editdistance(lt, lt_non)

                        rt = idl.rt_seq[:n]
                        rt_non = non.rt_seq[:n]
                        rt_edit_dist = editdistance(rt, rt_non)

                        complexity = lt_edit_dist + rt_edit_dist
                        complexities.append(complexity)
        if complexities == []:
            return indel_complexity_against_ref
        else:
            refined_value = min(complexities)
            return min(indel_complexity_against_ref, refined_value)


class CodingSequenceWithIndel(SequenceWithIndel):
    """Represents indel annotated with gene info
    
    Attributes:
        strand (str): '+' for positive strand '-' for negative
        accession (str): RefSeq accession number (e.g. NM_****)
        gene_symbol (str): gene name
        exon (int): exon number. 1 is the first exon 
        exon_start (int): the exon start pos on genome coordinate
        exon_end (int): the exon end pos on genome coordinate
        last_exon (int): 1 if the current exon is the last exon, 0 otherwise
        cds_start (int): the pos of coding sequence (cds) starting at the exon_start
        prev_exon_start (int): (current - 1) exon start pos on genome coordinate
                              -1 if current = 1 (first exon)
        prev_exon_end (int): (current - 1) exon end pos on genome coordinate
                              -1 if current = 1
        next_exon_start (int): (current + 1) exon start pos on genome coordinate
                              -1 if current = last exon
        next_exon_end (int): (current + 1) exon end pos on genome coordinate
                              -1 if current = last exon
    """

    def __init__(
        self,
        chr,
        pos,
        idl_type,
        lt_seq,
        idl_seq,
        rt_seq,
        strand,
        accession,
        gene_symbol,
        exon,
        exon_start,
        exon_end,
        last_exon,
        cds_start,
        prev_exon_start,
        prev_exon_end,
        next_exon_start,
        next_exon_end,
    ):

        SequenceWithIndel.__init__(self, chr, pos, idl_type, lt_seq, idl_seq, rt_seq)

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
        """Nonsense-mediatate decay (NMD) insensitivity
         
        Args:
            None
        Returns:
            is_insensitive (int): 1 if insensitive 0 otherwise
        """
        is_insensitive = 0

        if self.exon == 1 or self.exon == self.last_exon:
            is_insensitive = 1

        return is_insensitive

    def effect(self):
        """Report indel annotation based on the region where
        indel is annotated. 
        
        Possible regions:
            Exon, 
            Splice site (0 < dist.to exon boundary < 3)
            Splice region (2 < dist.to exon boundary < 11)

        Args:
            None
        Returns:
            indel annotation (str): see Example
        
        Example:
                           
            SDF4|NM_016547|167|frameshiftTruncating|0
           
        Pipe-delimited string reports GeneName, Accession, 
        Codon pos, Effect and NMD-insensitivity. 
        """
        if self.strand == "+":
            if self.exon_start <= self.pos <= self.exon_end:
                return self.exonic_on_pos_strand()
            elif (
                0 < self.exon_start - self.pos <= 2 or 0 < self.pos - self.exon_end <= 2
            ):
                return self.splice_site_on_pos_strand()
            elif (
                2 < self.exon_start - self.pos <= 11
                or 2 < self.pos - self.exon_end <= 11
            ):
                return self.splice_region_on_pos_strand()
            else:
                pass
        else:
            if self.exon_start <= self.pos <= self.exon_end:
                return self.exonic_on_neg_strand()
            elif (
                0 < self.exon_start - self.pos <= 2 or 0 < self.pos - self.exon_end <= 2
            ):
                return self.splice_site_on_neg_strand()
            elif (
                2 < self.exon_start - self.pos <= 11
                or 2 < self.pos - self.exon_end <= 11
            ):
                return self.splice_region_on_neg_strand()
            else:
                pass

    def cds_pos_in_exonic_indels(self):
        """Report coding sequence (CDS) pos affected by indel
        
        Args:
            None
        Returns:
            cds pos (int): The first coding sequence base affected by the indel
       
        Example:        1234567890123
                 CDS  : ATGCTACGACTGA
                  del : ATGCTA---CTGA  -> cds_pos = 7
                             
                        123456   7890123
                 CDS  : ATGCTA   CGACTGA             
                  ins : ATGCTATAGCGACTGA  -> cds_pos = 7

        Note that the sequences are unaffected upto first 6 bases. 
        """
        # insertion/deletion on positive strand
        if self.strand == "+":
            cds_pos = self.cds_start + self.pos - self.exon_start
        else:
            # insertion on negative strand
            if self.idl_type == 1:
                cds_pos = self.cds_start + self.exon_end - (self.pos - 1)
            # deletion on negative strand
            else:
                cds_pos = (
                    self.cds_start + self.exon_end - self.pos - (len(self.idl_seq) - 1)
                )
        return cds_pos

    def exonic_on_pos_strand(self):
        """Annotate coding exon indel on positve strand

        Args:
            None
        Returns:
            indel annotation (str): gene|acc|codon_pos|effect|nmd_insensitivity
            
            possible effect: frameshiftTruncating
                             inframeDel
                             inframeIns
                             nonsenseTruncating
                             spliceTruncating (the GT-AG motif broken)
                             splicePreserving (the GT-AG motif preserved)
           
           The splice effect is possible when insertion occurs at the 5'exon
           boundary.  
        """

        # insertion at 5'exon_start
        if self.idl_type == 1 and self.pos == self.exon_start:
            cds_pos = self.cds_start - 1
            codon_pos = int(cds_pos / 3) + 1
            if len(self.idl_seq) > 1 and self.idl_seq[-2:] == "AG":
                return codon_pos, "splicePreserving"
            else:
                return codon_pos, "spliceTruncating"

        # indels within exon
        else:
            cds_pos = self.cds_pos_in_exonic_indels()

            frame = (cds_pos - 1) % 3
            if frame == 2:
                codon_pos = int(cds_pos / 3)
            else:
                codon_pos = int(cds_pos / 3) + 1

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
            if exists_stop_codon(self.strand, seq):
                return codon_pos, "nonsenseTruncating"
            else:
                if len(self.idl_seq) % 3 == 0 and self.idl_type == 1:
                    return codon_pos, "inframeIns"
                elif len(self.idl_seq) % 3 == 0 and self.idl_type == 0:
                    return codon_pos, "inframeDel"
                else:
                    return codon_pos, "frameshiftTruncating"

    def splice_site_on_pos_strand(self):
        """Annotate indel within 2-nt to exon boundary on positve strand

        Args:
            None
        Returns:
           indel annotation (str): gene|acc|codon_pos|effect|nmd_insensitivity 
                                    
           possible effect: 
                           spliceShortIntron (for intron <= 5-nt)
                           splicePreserving (the GT-AG motif preserved)
                           spliceTruncating (the GT-AG motif broken)
                           spliceRegion (for ins at 2-nt upstream of 5'splice site)
        """
        # splicing motif + at least 1 base
        # GT + (at least one) + AG
        min_motif_len = 5

        # 5'splice
        if self.exon_start > self.pos:
            cds_pos = self.cds_start - 1
            codon_pos = int(cds_pos / 3) + 1

            if (self.exon_start - 1) - self.prev_exon_end <= min_motif_len:
                return codon_pos, "spliceShortIntron"

            else:
                # insertion at 1-nt upstream of exon start
                if self.idl_type == 1 and self.exon_start - self.pos == 1:
                    if self.idl_seq[-1] == "A":
                        return codon_pos, "splicePreserving"
                    else:
                        return codon_pos, "spliceTruncating"

                # insertion at 2-nt upstream of exon start
                elif self.idl_type == 1 and self.exon_start - self.pos == 2:
                    return codon_pos, "spliceRegion"

                # deletion at 1-nt upstream of exon start
                elif self.idl_type == 0 and self.exon_start - self.pos == 1:
                    return codon_pos, "spliceTruncating"

                # deletion at 2-nt upstream of exon start
                elif self.idl_type == 0 and self.exon_start - self.pos == 2:
                    if len(self.idl_seq) == 1 and self.lt_seq[-1] == "A":
                        return codon_pos, "splicePreserving"
                    elif len(self.idl_seq) == 2 and self.lt_seq[-2:] == "AG":
                        return codon_pos, "splicePreserving"
                    else:
                        return codon_pos, "spliceTruncating"
                else:
                    pass

        # 3'splice
        else:
            cds_pos = self.cds_start + self.exon_end - self.exon_start
            codon_pos = int(cds_pos / 3) + 1

            if self.next_exon_start - self.exon_end <= min_motif_len:
                return codon_pos, "spliceShortIntron"

            else:
                # insertion 1-nt downstream of exon end
                if self.idl_type == 1 and self.pos - self.exon_end == 1:
                    if len(self.idl_seq) > 1 and self.idl_seq[:2] == "GT":
                        return codon_pos, "splicePreserving"
                    else:
                        return codon_pos, "spliceTruncating"

                # insertion 2-nt downstream of exon end
                elif self.idl_type == 1 and self.pos - self.exon_end == 2:
                    if self.idl_seq[0] == "T":
                        return codon_pos, "splicePreserving"
                    else:
                        return codon_pos, "spliceTruncating"

                # deletion 1-nt downstream of exon end
                elif self.idl_type == 0 and self.pos - self.exon_end == 1:
                    if len(self.idl_seq) > 1 and self.rt_seq[:2] == "GT":
                        return codon_pos, "splicePreserving"
                    else:
                        return codon_pos, "spliceTruncating"

                # deletion 2-nt downstream of exon end
                elif self.idl_type == 0 and self.pos - self.exon_end == 2:
                    if self.rt_seq[0] == "T":
                        return codon_pos, "splicePreserving"
                    else:
                        return codon_pos, "spliceTruncating"
                else:
                    pass

    def splice_region_on_pos_strand(self):
        """Annotate indel in splice region on positive strand
        
        Splice region is defined intronic region where
             2 < distance to the exon boundary < 11
        
        Args:
            None
        Returns:
            indel annotation (str): gene|acc|codon_pos|effect|nmd_insensitivity
            
            possible effect: spliceRegion
        """
        # 5'splice region
        if self.exon_start > self.pos:
            cds_pos = self.cds_start - 1
            codon_pos = int(cds_pos / 3) + 1
            return codon_pos, "spliceRegion"

        # 3'splice region
        else:
            cds_pos = self.cds_start + self.exon_end - self.exon_start
            codon_pos = int(cds_pos / 3) + 1
            return codon_pos, "spliceRegion"

    def exonic_on_neg_strand(self):
        """Annotate coding indel on negative strand
        
         Args:
             None
         Returns:
             indel annotation (str): gene|acc|codon_pos|effect|nmd_insensitivity
         
             possible effect: frameshiftTruncating
                              inframeDel
                              inframeIns
                              nonsenseTruncating
                              spliceTruncating (the GT-AG motif broken)
                              splicePreserving (the GT-AG motif preserved)
         
             The splice effect is possible when insertion occurs at the 3'exon
             boundary.
        """

        # insertion at 3'exon_start
        if self.idl_type == 1 and self.pos == self.exon_start:
            cds_pos = self.cds_start + self.exon_end - self.pos
            codon_pos = int(cds_pos / 3) + 1
            if len(self.idl_seq) > 1 and self.idl_seq[-2:] == "AC":
                return codon_pos, "splicePreserving"
            else:
                return codon_pos, "spliceTruncating"

        # indels within exon
        else:
            cds_pos = self.cds_pos_in_exonic_indels()

            frame = (cds_pos - 1) % 3
            if frame == 2:
                codon_pos = int(cds_pos / 3)
            else:
                codon_pos = int(cds_pos / 3) + 1

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
            if exists_stop_codon(self.strand, seq):
                return codon_pos, "nonsenseTruncating"
            else:
                if len(self.idl_seq) % 3 == 0 and self.idl_type == 1:
                    return codon_pos, "inframeIns"
                elif len(self.idl_seq) % 3 == 0 and self.idl_type == 0:
                    return codon_pos, "inframeDel"
                else:
                    return codon_pos, "frameshiftTruncating"

    def splice_site_on_neg_strand(self):
        """Annotate indel within 2-nt to exon boundary on negative strand

        Args:
            None
        Returns:
            indel annotation (str): gene|acc|codon_pos|effect|nmd_insensitivity

            possible effect:
                            spliceShortIntron (for intron <= 5-nt)
                            splicePreserving (the GT-AG motif preserved)
                            spliceTruncating (the GT-AG motif broken)
                            spliceRegion (for ins at 2-nt upstream of 3'splice site)
        """

        # splicing motif + at least 1 base
        # GT + (at least one) + AG
        min_motif_len = 5

        # 5'splice
        if self.pos > self.exon_end:
            cds_pos = self.cds_start - 1
            codon_pos = int(cds_pos / 3) + 1

            if (self.prev_exon_start - 1) - self.exon_end <= min_motif_len:
                return codon_pos, "spliceShortIntron"

            else:
                # insertion at 1-nt downstream of exon end
                if self.idl_type == 1 and self.pos - self.exon_end == 1:
                    if len(self.idl_seq) > 1 and self.idl_seq[:2] == "CT":
                        return codon_pos, "splicePreserving"
                    else:
                        return codon_pos, "spliceTruncating"

                # insertion at 2-nt downstream of exon end
                elif self.idl_type == 1 and self.pos - self.exon_end == 2:
                    if self.idl_seq[0] == "T":
                        return codon_pos, "splicePreserving"
                    else:
                        return codon_pos, "spliceTruncating"

                # deletion at 1-nt downstream of exon end
                elif self.idl_type == 0 and self.pos - self.exon_end == 1:
                    if len(self.idl_seq) > 1 and self.rt_seq[:2] == "CT":
                        return codon_pos, "splicePreserving"
                    else:
                        return codon_pos, "spliceTruncating"

                # deletion at 2-nt downstream of exon end
                elif self.idl_type == 0 and self.pos - self.exon_end == 2:
                    if self.rt_seq[0] == "T":
                        return codon_pos, "splicePreserving"
                    else:
                        return codon_pos, "spliceTruncating"
                else:
                    pass

        # 3'splice
        else:
            cds_pos = self.cds_start + self.exon_end - self.exon_start
            codon_pos = int(cds_pos / 3) + 1

            if self.exon_start - self.next_exon_end <= min_motif_len:
                return codon_pos, "spliceShortIntron"

            else:
                # insertion 1-nt upstream of exon start
                if self.idl_type == 1 and self.exon_start - self.pos == 1:
                    if self.idl_seq[-1] == "A":
                        return codon_pos, "splicePreserving"
                    else:
                        return codon_pos, "spliceTruncating"

                # insertion 2-nt upstream of exon start
                elif self.idl_type == 1 and self.exon_start - self.pos == 2:
                    return codon_pos, "spliceRegion"

                # deletion 1-nt upstream of exon start
                elif self.idl_type == 0 and self.exon_start - self.pos == 1:
                    return codon_pos, "spliceTruncating"

                # deletion 2-nt upstream of exon start
                elif self.idl_type == 0 and self.exon_start - self.pos == 2:
                    if len(self.idl_seq) == 2 and self.lt_seq[-2:] == "AC":
                        return codon_pos, "splicePreserving"
                    else:
                        return codon_pos, "spliceTruncating"
                else:
                    pass

    def splice_region_on_neg_strand(self):
        """Annotate indel in splice region on negative strand

        Splice region is defined intronic region where
             2 < distance to the exon boundary < 11
         
        Args:
            None
        Returns:
            indel annotation (str): gene|acc|codon_pos|effect|nmd_insensitivity
         
            possible effect: spliceRegion
        """
        # 5'splice region
        if self.pos > self.exon_end:
            cds_pos = self.cds_start - 1
            codon_pos = int(cds_pos / 3) + 1
            return codon_pos, "spliceRegion"

        # 3'splice region
        else:
            cds_pos = self.cds_start + self.exon_end - self.exon_start
            codon_pos = int(cds_pos / 3) + 1
            return codon_pos, "spliceRegion"
