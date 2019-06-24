#!/usr/bin/env python3

class Variant(object):
    """Represents genomic variants based on VCF record
    """

    def __init__(self, chrom, pos, ref, alt, genome):
        """Constructor 
        
        Args:
            chr (str): chromosome names
            pos (int): 1-based position
            ref (str): reference allele
            alt (str): alternative allele
            genome (pysam.FastaFile): reference genome
        """
        self.chrom = self.format_chrom_name(chrom, genome)
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.genome = genome

        # check if variant is in valid genomic locus
        self.validate_locus()

        # to store dbSNP IDs (a single variant may have multiple IDs)
        self._ids = []
        # to store if the dbSNP IDs are annotated as COMMON
        self._common = []

    def format_chrom_name(self, chrom, genome):
        chrom_names = genome.references
        is_prefixed = True if chrom_names[0].startswith("chr") else False
        is_mt = True if "chrMT" in chrom_names or "MT" in chrom_names else False

        chrom = chrom.replace("chr", "")
        if chrom == "M" and is_mt:
            chrom = "MT"
        elif chrom == "MT" and not is_mt:
            chrom = "M"

        if is_prefixed:
            chrom = "chr" + chrom

        return chrom

    def validate_allele(self):
        """Check if allele consists of A,C,T,G

        Args:
            None
        Returns:
            None
        Raises:
            Exceptions: if None, empty str, str with char other than ACTG 
        """
        if not self.ref or not self.alt:
            raise Exception("Allele may not be None")

        reflst, altlst = list(self.ref), list(self.alt)
        if not reflst or not altlst:
            raise Exception("Allele may not be empty")

        bases = {"A", "C", "T", "G"}
        if not set(reflst) <= base or not set(altlst) <= base:
            raise Exception("Allele contains char other than A, C, T, G")

    def validate_locus(self):
        """Check if the locus is curated in the reference
       
        Args:
            None
        Returns:
            None
        Raises:
            Exception: if invalid
        """
        # check if contig is valide
        try:
            # note pos is converted into 0-based
            seq = self.genome.fetch(self.chrom, (self.pos - 1) - 1, self.pos - 1)
        except:
            raise Exception("This contig is not contained in the reference")

        # check if the base is available and non-N
        if len(seq) < 1 or "N" in seq:
            raise Exception("This locus is not curated by the reference")

    def normalize(self):
        """Nomalize to a left-aligned and minimal VCF representation
           
        A Python implementation of Tan et al. Bioinformatics, 2015, 2202-2204
        
        Args:
            None
        Returns:
            Variant: normalized copy
                     (I prefer to return a copy to keep the original representiaon)
        """
        i = Variant(self.chrom, self.pos, self.ref, self.alt, self.genome)

        condition_1 = i.ref[-1] == i.alt[-1]
        while condition_1:
            # note 0-based
            left_base = i.genome.fetch(i.chrom, (i.pos - 1) - 1, i.pos - 1)
            i.ref = left_base + i.ref[:-1]
            i.alt = left_base + i.alt[:-1]
            i.pos -= 1
            condition_1 = i.ref[-1] == i.alt[-1]

        condition_2 = i.ref[0] == i.alt[0]
        condition_3 = len(i.ref) > 1 and len(i.alt) > 1
        while condition_2 and condition_3:
            i.ref = i.ref[1:]
            i.alt = i.alt[1:]
            i.pos += 1
            condition_2 = i.ref[0] == i.alt[0]
            condition_3 = len(i.ref) > 1 and len(i.alt) > 1
        return i

    @property
    def is_indel(self):
        """Ask if this is an indel

        Args:
            None
        Returns:
            bool: True if indel
        """
        return False if len(self.ref) == len(self.alt) else True

    def __eq__(self, other):
        """Ask if this is equivalent to another variant

        Args:
            None
        Returns:
            bool: True if they are equivalent to each other
        """
        i, j = self.normalize(), other.normalize()

        equivalent = (
            i.chrom.replace("chr", "") == j.chrom.replace("chr", "")
            and i.pos == j.pos
            and j.ref == i.ref
            and i.alt == j.alt
        )

        return equivalent

    def __hash__(self):
        """Make this hashable

        Args:
            None
        Returns: 
            hash: hashed value
        """
        i = self.normalize() if self.is_indel else self
        hashable = (i.chrom, i.pos, i.ref, i.alt)

        return hash(hashable)


def make_indel_from_vcf_line(line, genome):
    """Make a list of indel obj from a vcf line

    Args: 
        line (str): a VCF line
        genome (pysam.FastaFile)
    Returns:
        indels (list):[indel_obj_1, indel_obj_2]
                None: if no indel created (e.g., a VCF line for SNV)
    """
    if line.startswith("#"):
        return None

    lst = line.rstrip().split("\t")
    chrom, pos = lst[0], int(lst[1])
    ref, alts = lst[3], lst[4].split(",")

    # SNVs
    if all(len(ref) == len(alt) == 1 for alt in alts):
        return None

    indels = [
        Variant(chrom, pos, ref, alt, genome)
        for alt in alts
        if Variant(chrom, pos, ref, alt, genome).is_indel
    ]

    if indels:
        return indels
    else:
        return None

