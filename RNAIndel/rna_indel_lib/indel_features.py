#!/usr/bin/env python3

from .indel_sequence import Indel


class AnnotationFeatures(object):
    """ Store annotation summary 

    Attributes:
        is_inframe (int): 1 if true 0 otherwise
        is_truncation (int): 1 if true 0 otherwise
        is_splice (int): 1 if true 0 otherwise
        is_nmd_insensitive: 1 if true 0 otherwise
    """

    def __init__(self, is_inframe, is_truncating, is_splice, is_nmd_insensitive):

        self.is_inframe = is_inframe
        self.is_truncating = is_truncating
        self.is_splice = is_splice
        self.is_nmd_insensitive = is_nmd_insensitive


class SamFeatures(object):
    """Store Sequence/Alignment (SAM) features
    """

    def __init__(
        self,
        gc,
        lc,
        strength,
        local_gc,
        local_lc,
        local_strength,
        repeat,
        dissimilarity,
        indel_complexity,
        ref_count,
        alt_count,
        is_multiallelic,
        is_near_boundary,
        is_bidirectional,
        is_uniq_mapped,
    ):

        self.gc = gc
        self.lc = lc
        self.strength = strength
        self.local_gc = local_gc
        self.local_lc = local_lc
        self.local_strength = local_strength
        self.repeat = repeat
        self.dissimilarity = dissimilarity
        self.indel_complexity = indel_complexity
        self.ref_count = ref_count
        self.alt_count = alt_count
        self.is_multiallelic = is_multiallelic
        self.is_near_boundary = is_near_boundary
        self.is_bidirectional = is_bidirectional
        self.is_uniq_mapped = is_uniq_mapped


class IndelSnpFeatures(Indel):
    """Store and summarize dbSNP and ClinVar info
    
    Attributes:
        chr (str): chr1-22, chrX, chrY
        pos (int): 1-based
        idl_type (int): 1 for insertion 0 for deletion
        idl_seq (str): inserted or deleted sequence
    """

    def __init__(self, chr, pos, idl_type, idl_seq):
        Indel.__init__(self, chr, pos, idl_type, idl_seq)

        self.dbsnp_id = []
        self.dbsnp_freq = []
        self.dbsnp_common = []
        self.dbsnp_origin = []
        self.clnvr_id = []
        self.clnvr_freq = []
        self.clnvr_info = []
        self.clnvr_origin = []

    def add_dbsnp_id(self, rs):
        """Add dbSNP ID
         
        Args:
            rs (str): dbSNP ID (starting with 'rs')
        Returns:
            None
        """
        self.dbsnp_id.append(rs)

    def report_dbsnp_id(self):
        """Report dbSNP ID 

        Args:
            None
        Returns:
            dbSNP ID (str): '-' if the indel is not on dbSNP
                            delimited with ',' if multiple IDs found
        """
        if self.dbsnp_id == []:
            return "-"
        else:
            return ",".join(self.dbsnp_id)

    def add_clnvr_id(self, id):
        """Add ClinVar ID

        Args:
            id (str): ClinVar ID
        Returns:
            None
        """
        self.clnvr_id.append(id)

    def report_clnvr_id(self):
        """Report ClinVar ID

        Args:
            None
        Returns:
            ClinVar (str): '-' if the indel is not on ClinVar
                           delimited with ',' if multiple IDs found
        """
        if self.clnvr_id == []:
            return "-"
        else:
            return ",".join(self.clnvr_id)

    def add_dbsnp_freq(self, freq):
        """Add MAF on dbSNP

        Args:
            freq (float): minor allele frequency reported on dbSNP
        Returns:
            None
        """
        self.dbsnp_freq.append(freq)

    def add_clnvr_freq(self, freq):
        """Add MAF on ClinVar

        Args:
            freq (float): minor allele frequency reported on ClinVar
        Returns:
            None
        """
        self.clnvr_freq.append(freq)

    def report_freq(self):
        """Report maximun MAF

        Args:
            None
        Returns:
            max MAF: max(dbSNP_MAF, ClinVAR_MAF)
                     -1 if no MAF info is available
        """
        freqs = self.dbsnp_freq + self.clnvr_freq

        if freqs == []:
            return -1
        else:
            return max(freqs)

    def add_dbsnp_common(self, common):
        """Add 'Common' annotation on dbSNP

        Args:
            common (int): 1 if common
                          0 if uncommon
                         -1 if no info
        Returns:
            None
        """
        self.dbsnp_common.append(common)

    def is_common(self):
        """Encodes if the indel is common

        Args:
            None
        Returns:
            is_common: 1 if annotated 'Common' on dbSNP
                            or MAF >= 0.01
                      -1 if no info available
                       0 otherwise
        """
        is_common = 0

        if self.dbsnp_common == [] and self.report_freq() == "-":
            is_common = -1
        elif 1 in self.dbsnp_common:
            is_common = 1
        elif self.report_freq() != "-":
            if self.report_freq() >= 0.01:
                is_common = 1
        else:
            pass

        return is_common

    def add_dbsnp_origin(self, origin):
        """Add the tissue origin of dbSNP report

        Args:
            origin (int)
        Returns:
            None
        """
        if origin == 0:
            ori = "unspecified"
        elif origin == 1:
            ori = "germline"
        elif origin == 2:
            ori = "somatic"
        elif origin == 3:
            ori = "both"
        else:
            ori = "-"

        self.dbsnp_origin.append(ori)

    def add_clnvr_origin(self, origin):
        """Add the genetic origin of ClinVar report

        Args:
            origin (int)
        Returns:
            None
        """
        if origin == 0:
            ori = "unknown"
        elif origin == 1:
            ori = "germline"
        elif origin == 2:
            ori = "somatic"
        # 4 = 'inherited'
        elif origin == 4:
            ori = "germline"
        # 8 = 'paternal'
        elif origin == 8:
            ori = "germline"
        # 16 = 'maternal'
        elif origin == 16:
            ori = "germline"
        elif origin == 32:
            ori = "de-novo"
        # 64 = 'biparental'
        elif origin == 64:
            ori = "germline"
        # 128 = 'uniparental'
        elif origin == 128:
            ori = "germline"
        elif origin == 256:
            ori = "not-tested"
        elif origin == 512:
            ori = "tested-inconclusive"
        elif origin == 1073741824:
            ori = "others"
        else:
            ori = "-"

        self.clnvr_origin.append(ori)

    def with_germline_reports(self):
        """Summarise origin report
        
        Args:
            None
        Returns:
            with_germline_reports: 1 if germline report found
                                   0 otherwise
                                  -1 for no info
        """
        origins = self.dbsnp_origin + self.clnvr_origin
        if origins == []:
            return -1
        elif "germline" in origins and not "somatic" in origins:
            return 1
        elif "germline" in origins and not "both" in origins:
            return 1
        else:
            return 0

    def add_clnvr_info(self, clninfo):
        """Add pathogenecity info

        Args:
            clninfo (str): significance|disease name
        Returns:
            None
        """
        self.clnvr_info.append(clninfo)

    def report_clnvr_info(self):
        """Report pathogenecity info

        Args:
            None
        Returns:
            clnvr_info (str): significance|disease name
                              '-' if no info available
        """
        if self.clnvr_info == []:
            return "-"
        else:
            return ",".join(self.clnvr_info)

    def is_not_pathogenic(self):
        """Encode if the indel is not pathogenic

        Args:
            None
        Returns:
            is_not_pathogenic: 1 if not pathogenic
                               0 otherwise
                              -1 if no info
        """
        info = self.report_clnvr_info()

        if info == "-":
            return -1

        info = info.lower()
        if "pathogenic" not in info:
            if "benign" in info:
                return 1
            elif "protective" in info:
                return 1
            elif "affects" in info:
                return 1
            else:
                return 0
        else:
            return 0
