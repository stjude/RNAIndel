from .sequence_properties import exists_stop_codon

def annotate_coding_info(indel, coding_gene_db):
    """Generate coding indel objects
    
    Args:
        chr (str): chr1-22, chrX or chrY. Note "chr"-prefixed.
        pos (int): 1-based genomic position
        idl_type (int): 1 for insertion, 0 for deletion
        idl_seq (str): inserted or deleted sequence
        genome (pysam.FastaFile): reference genome
        exons (pysam.TabixFile): coding exon data
        chr_prefixed (bool): True if chromosome names in BAM are "chr"-prefixed

    Returns:
        coding_idl_lst (list): a list of CodingSequenceWithIndel obj
                               empty list if non-coding indel  
    """
    coding_annots = []
    chrom, pos, indel_type, indel_seq = (
        indel.chrom,
        indel.pos,
        indel.variant_type,
        indel.indel_seq,
    )

    indel_type_db = 1 if indel_type == "I" else 0  # conversion for coding_gene_db style
    pos = pos + 1  # conversion for bambino style
    margin = 11  # allow a margin of 11bp for splice region
    try:
        candidate_genes = coding_gene_db.fetch(chrom, pos - margin, pos + margin)
    except:
        candidate_genes = None

    # check for UTR
    if candidate_genes:
        for line in candidate_genes:
            lst = line.split("\t")

            # parsing exon info
            info = lst[3].split("|")
            exon = int(info[2])
            last_exon = int(info[3])

            # exon start and end
            exon_start, exon_end = int(lst[1]), int(lst[2])

            # strand
            strand = lst[4]

            # 5'UTR on positive strand (insertion)
            if strand == "+" and exon == 1 and indel_type_db == 1 and exon_start >= pos:
                pass
            # 5'UTR on positive strand (deletion)
            elif (
                strand == "+" and exon == 1 and indel_type_db == 0 and exon_start > pos
            ):
                pass
            # 3'UTR on positive strand
            elif strand == "+" and exon == last_exon and pos > exon_end:
                pass
            # 5'UTR on negative strand
            elif strand == "-" and exon == 1 and pos > exon_end:
                pass
            # 3'UTR on negative strand (insertion)
            elif (
                strand == "-"
                and exon == last_exon
                and indel_type_db == 1
                and exon_start >= pos
            ):
                pass
            # 3'UTR on negative strand (deletion)
            elif (
                strand == "-"
                and exon == last_exon
                and indel_type_db == 0
                and exon_start > pos
            ):
                pass
            else:
                accession = info[0]
                gene_symbol = info[1]
                cds_start = int(info[4])
                cds_len = int(info[5])
                prev_exon = lst[5].split("|")
                prev_exon_start, prev_exon_end = int(prev_exon[0]), int(prev_exon[1])
                next_exon = lst[6].split("|")
                next_exon_start, next_exon_end = int(next_exon[0]), int(next_exon[1])

                coding_annot = CodingAnnotation(
                    indel,
                    strand,
                    accession,
                    gene_symbol,
                    exon,
                    exon_start,
                    exon_end,
                    last_exon,
                    cds_start,
                    cds_len,
                    prev_exon_start,
                    prev_exon_end,
                    next_exon_start,
                    next_exon_end,
                )

                coding_annots.append(coding_annot)

        return coding_annots


def get_gene_symbol(row):
    """Extracts gene name from annotation

    Args:
        row (pandas.Series): annotation info (str) at 'annotation' index
    Returns:
        gene_symbol (str): gene name(s)
    """
    pd.options.mode.chained_assignment = None

    lst = row["annotation"].split(",")
    genes = [token.split("|")[0] for token in lst]

    gene_symbol = ",".join(set(genes))

    return gene_symbol


class CodingAnnotation(object):
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
        indel,
        strand,
        accession,
        gene_symbol,
        exon,
        exon_start,
        exon_end,
        last_exon,
        cds_start,
        cds_len,
        prev_exon_start,
        prev_exon_end,
        next_exon_start,
        next_exon_end,
    ):
        self.indel = indel
        self.indel_type = indel.variant_type
        self.indel_seq = indel.indel_seq
        self.lt_reference_flank = indel.left_flank()
        self.rt_reference_flank = indel.right_flank()

        self.pos = indel.pos + 1  # conversion for bambino style

        self.strand = strand
        self.accession = accession
        self.gene_symbol = gene_symbol
        self.exon = exon
        self.exon_start = exon_start
        self.exon_end = exon_end
        self.last_exon = last_exon
        self.cds_start = cds_start
        self.cds_len = cds_len
        self.prev_exon_start = prev_exon_start
        self.prev_exon_end = prev_exon_end
        self.next_exon_start = next_exon_start
        self.next_exon_end = next_exon_end

        self.codon_pos, self.variant_effect = self.effect()

    def to_str(self):
        annotation_string = "{}|{}|{}|{}|{}".format(
            self.gene_symbol,
            self.accession,
            self.codon_pos,
            self.variant_effect,
            self.is_nmd_insensitive(),
        )
        return annotation_string

    def is_nmd_insensitive(self):
        if self.exon == 1 or self.exon == self.last_exon:
            is_insensitive = 1
        else:
            is_insensitive = 0

        return is_insensitive

    def get_relative_location(self):
        return (self.codon_pos * 3) / self.cds_len

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
            if self.indel_type == "I":
                cds_pos = self.cds_start + self.exon_end - (self.pos - 1)
            # deletion on negative strand
            else:
                cds_pos = (
                    self.cds_start
                    + self.exon_end
                    - self.pos
                    - (len(self.indel_seq) - 1)
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
        if self.indel_type == "I" and self.pos == self.exon_start:
            cds_pos = self.cds_start - 1
            codon_pos = int(cds_pos / 3) + 1
            if len(self.indel_seq) > 1 and self.indel_seq[-2:] == "AG":
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
            if self.indel_type == "I":
                if frame == 0:
                    seq = self.indel_seq + self.rt_reference_flank[:2]
                elif frame == 1:
                    seq = (
                        self.lt_reference_flank[-1:]
                        + self.indel_seq
                        + self.rt_reference_flank[:1]
                    )
                else:
                    seq = (
                        self.lt_reference_flank[-2:]
                        + self.indel_seq
                        + self.rt_reference_flank[:3]
                    )
            # deletion
            else:
                if frame == 0:
                    seq = self.rt_reference_flank[:3]
                elif frame == 1:
                    seq = self.lt_reference_flank[-1:] + self.rt_reference_flank[:2]
                else:
                    seq = self.lt_reference_flank[-2:] + self.rt_reference_flank[:1]
            # check for stop codon
            if exists_stop_codon(self.strand, seq):
                return codon_pos, "nonsenseTruncating"
            else:
                if len(self.indel_seq) % 3 == 0 and self.indel_type == "I":
                    return codon_pos, "inframeIns"
                elif len(self.indel_seq) % 3 == 0 and self.indel_type == "D":
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
                if self.indel_type == "I" and self.exon_start - self.pos == 1:
                    if self.indel_seq[-1] == "A":
                        return codon_pos, "splicePreserving"
                    else:
                        return codon_pos, "spliceTruncating"

                # insertion at 2-nt upstream of exon start
                elif self.indel_type == "I" and self.exon_start - self.pos == 2:
                    return codon_pos, "spliceRegion"

                # deletion at 1-nt upstream of exon start
                elif self.indel_type == "D" and self.exon_start - self.pos == 1:
                    return codon_pos, "spliceTruncating"

                # deletion at 2-nt upstream of exon start
                elif self.indel_type == "D" and self.exon_start - self.pos == 2:
                    if len(self.indel_seq) == 1 and self.lt_reference_flank[-1] == "A":
                        return codon_pos, "splicePreserving"
                    elif (
                        len(self.indel_seq) == 2
                        and self.lt_reference_flank[-2:] == "AG"
                    ):
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
                if self.indel_type == "I" and self.pos - self.exon_end == 1:
                    if len(self.indel_seq) > 1 and self.indel_seq[:2] == "GT":
                        return codon_pos, "splicePreserving"
                    else:
                        return codon_pos, "spliceTruncating"

                # insertion 2-nt downstream of exon end
                elif self.indel_type == "I" and self.pos - self.exon_end == 2:
                    if self.indel_seq[0] == "T":
                        return codon_pos, "splicePreserving"
                    else:
                        return codon_pos, "spliceTruncating"

                # deletion 1-nt downstream :of exon end
                elif self.indel_type == "D" and self.pos - self.exon_end == 1:
                    if len(self.indel_seq) > 1 and self.rt_reference_flank[:2] == "GT":
                        return codon_pos, "splicePreserving"
                    else:
                        return codon_pos, "spliceTruncating"

                # deletion 2-nt downstream of exon end
                elif self.indel_type == "D" and self.pos - self.exon_end == 2:
                    if self.rt_reference_flank[0] == "T":
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
        if self.indel_type == "I" and self.pos == self.exon_start:
            cds_pos = self.cds_start + self.exon_end - self.pos
            codon_pos = int(cds_pos / 3) + 1
            if len(self.indel_seq) > 1 and self.indel_seq[-2:] == "AC":
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
            if self.indel_type == "I":
                if frame == 0:
                    seq = self.lt_reference_flank[-2:] + self.indel_seq
                elif frame == 1:
                    seq = (
                        self.lt_reference_flank[-1:]
                        + self.indel_seq
                        + self.rt_reference_flank[:1]
                    )
                else:
                    seq = (
                        self.lt_reference_flank[-3:]
                        + self.indel_seq
                        + self.rt_reference_flank[:2]
                    )
            # deletion
            else:
                if frame == 0:
                    seq = self.lt_reference_flank[-3:]
                elif frame == 1:
                    seq = self.lt_reference_flank[-2:] + self.rt_reference_flank[:1]
                else:
                    seq = self.lt_reference_flank[-1:] + self.rt_reference_flank[:2]
            # check for stop codon
            if exists_stop_codon(self.strand, seq):
                return codon_pos, "nonsenseTruncating"
            else:
                if len(self.indel_seq) % 3 == 0 and self.indel_type == "I":
                    return codon_pos, "inframeIns"
                elif len(self.indel_seq) % 3 == 0 and self.indel_type == "D":
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
                if self.indel_type == "I" and self.pos - self.exon_end == 1:
                    if len(self.indel_seq) > 1 and self.indel_seq[:2] == "CT":
                        return codon_pos, "splicePreserving"
                    else:
                        return codon_pos, "spliceTruncating"

                # insertion at 2-nt downstream of exon end
                elif self.indel_type == "I" and self.pos - self.exon_end == 2:
                    if self.indel_seq[0] == "T":
                        return codon_pos, "splicePreserving"
                    else:
                        return codon_pos, "spliceTruncating"

                # deletion at 1-nt downstream of exon end
                elif self.indel_type == "D" and self.pos - self.exon_end == 1:
                    if len(self.indel_seq) > 1 and self.rt_reference_flank[:2] == "CT":
                        return codon_pos, "splicePreserving"
                    else:
                        return codon_pos, "spliceTruncating"

                # deletion at 2-nt downstream of exon end
                elif self.indel_type == "D" and self.pos - self.exon_end == 2:
                    if self.rt_reference_flank[0] == "T":
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
                if self.indel_type == "I" and self.exon_start - self.pos == 1:
                    if self.indel_seq[-1] == "A":
                        return codon_pos, "splicePreserving"
                    else:
                        return codon_pos, "spliceTruncating"

                # insertion 2-nt upstream of exon start
                elif self.indel_type == "I" and self.exon_start - self.pos == 2:
                    return codon_pos, "spliceRegion"

                # deletion 1-nt upstream of exon start
                elif self.indel_type == "D" and self.exon_start - self.pos == 1:
                    return codon_pos, "spliceTruncating"

                # deletion 2-nt upstream of exon start
                elif self.indel_type == "D" and self.exon_start - self.pos == 2:
                    if (
                        len(self.indel_seq) == 2
                        and self.lt_reference_flank[-2:] == "AC"
                    ):
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
