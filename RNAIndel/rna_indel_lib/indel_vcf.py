#!/usr/bin/env python3

import re
from .indel_sequence import Indel
from .left_aligner import peek_left_base

snp_ptn = re.compile(r"rs[0-9]+")


class IndelVcfReport(object):
    """Represent VCF record and meta info of indel 
    specified by Bambino coordinate

    Attributes:
        fa (pysam.FastaFile): representing the refernce
        chr (str): chr1-22, chrX, chrY
        pos (int): 1-based indel pos (Bambino coordinate)
        ref (str): Bambino style ref allele
        alt (str): Bambino style alt allele
    """

    def __init__(self, fa, chr, pos, ref, alt):
        self.fa = fa
        self.chr = chr
        self.pos = pos
        self.ref = ref
        self.alt = alt

    def generate_indel(self):
        if self.ref == "-":
            idl = Indel(self.chr, self.pos, 1, self.alt)
        else:
            idl = Indel(self.chr, self.pos, 0, self.ref)

        return idl

    @property
    def CHROM(self):
        return self.chr

    @property
    def POS(self):
        return self.pos - 1

    @property
    def ID(self):
        """Return ID info (typically dbSNP ID).
        If info is none or not defined (not set), return the 
        defined missing value: a dot '.'

        Returns:
            self.__ID (str) or . (str)
        """
        try:
            if self.__ID == None or self.__ID == "-":
                return "."
            else:
                return self.__ID
        except:
            return "."

    @ID.setter
    def ID(self, ID):
        """Set dbSNP IDs as semi-colon delimited list
        
        Args:
            ID (str): dbSNP IDs with any delimiter
        """
        rs_lst = snp_ptn.findall(ID)
        self.__ID = ";".join(rs_lst)

    @property
    def REF(self):
        """Convert ref allele Bambino to VCF
        
        Returns:
            for insertion, base at pos - 1 (str)
            for deletion, bsse at pos - 1 + Bambino's ref base (str)
        Example:
            pos       12345678 9012
            referece  ATGATGAT TAGA
            ins       ATGATGATCTAGA
            del       ATG-TGAT TAGA
            
            for insertion
                Bambino: ref = '-', alt = 'C' at 9
                VCF: REF = 'T', ALT = 'TC' at 8
            for deletion
                Bambino: ref = 'A', alt = '-' at 4
                VCF: REF = 'GA', ALT = 'G' at 3
        """
        if self.ref == "-":
            return peek_left_base(self.generate_indel(), self.fa)
        else:
            return peek_left_base(self.generate_indel(), self.fa) + self.ref

    @property
    def ALT(self):
        """Convert alt allele Bambino to VCF

        Returs:
            for insertion, base at pos - 1 + Bambino's alt (str)
            for deletion, base at pos - 1 
        """
        if self.ref == "-":
            return self.REF + self.alt
        else:
            return self.REF.replace(self.ref, "")

    @property
    def QUAL(self):
        """Return QUAL info (typically Phred-scaled quality score)
        If info is none or not defined (not set), return the defined 
        missing value: a dot '.'

        Returns:
            self.__QUAL (int) or . (str)
        """
        try:
            if self.__QUAL == None:
                return "."
            else:
                return self.__QUAL
        except:
            return "."

    @QUAL.setter
    def QUAL(self, QUAL):
        """Set QUAL info
        
        Args:
            QUAL (int)
        """
        self.__QUAL = QUAL

    @property
    def FILTER(self):
        """Return FILTER info 
        If info is none or not defined (not set), return 'PASS'

        Returns:
            self.__FILTER (str) or 'PASS' (str)
        """
        try:
            if self.__FILTER == None or self.__FILTER == "-":
                return "PASS"
            else:
                return self.__FILTER
        except:
            return "PASS"

    @FILTER.setter
    def FILTER(self, FILTER):
        """Set FILTER info
        
        Args:
            FILTER (str)
        """
        self.__FILTER = FILTER

    @property
    def vcf_record(self):

        record = [
            self.CHROM,
            str(self.POS),
            self.ID,
            self.REF,
            self.ALT,
            str(self.QUAL),
            self.FILTER,
            self.INFO,
            self.FORMAT,
        ]
        return "\t".join(record)

    ################
    #  INFO fields #
    ################
    @property
    def INFO(self):
        pred = "PRED=" + self.PRED

        prob = "PROB=" + ",".join([str(p) for p in self.PROB])

        db = "DB"
        if self.DB == "-":
            db = ""

        anno = "ANNO=" + self.ANNO

        maxmaf = "MAXMAF=" + str(self.MAXMAF)
        if self.MAXMAF == -1:
            maxmaf = ""

        common = "COMMON"
        if self.COMMON != 1:
            common = ""

        clin = "CLIN=" + self.CLIN
        if self.CLIN == "-":
            clin = ""

        icp = "ICP=" + str(self.ICP)

        dsm = "DSM=" + str(self.DSM)

        isz = "ISZ=" + str(self.ISZ)

        rep = "REP=" + str(self.REP)

        uqm = "UQM"
        if self.UQM == 0:
            uqm = ""

        neb = "NEB"
        if self.NEB == 0:
            neb = ""

        bid = "BID"
        if self.BID == 0:
            bid = ""

        mta = "MTA"
        if self.MTA == 0:
            mta = ""

        nmd = "NMD"
        if self.NMD == 0:
            nmd = ""

        ipg = "IPG=" + str(self.IPG)

        lsg = "LSG=" + str(self.LSG)

        ati = "ATI"
        if self.ATI == 0:
            ati = ""

        atd = "ATD"
        if self.ATD == 0:
            atd = ""

        info_lst = [
            pred,
            prob,
            db,
            anno,
            maxmaf,
            common,
            clin,
            icp,
            dsm,
            isz,
            rep,
            uqm,
            neb,
            bid,
            mta,
            nmd,
            ipg,
            lsg,
            ati,
            atd,
        ]

        return ";".join([i for i in info_lst if i != ""])

    @INFO.setter
    def INFO(self, INFO):
        self.PRED = INFO["PRED"]
        self.PROB = INFO["PROB"]
        self.DB = INFO["DB"]
        self.ANNO = INFO["ANNO"]
        self.MAXMAF = INFO["MAXMAF"]
        self.COMMON = INFO["COMMON"]
        self.CLIN = INFO["CLIN"]
        self.ICP = INFO["ICP"]
        self.DSM = INFO["DSM"]
        self.ISZ = INFO["ISZ"]
        self.REP = INFO["REP"]
        self.UQM = INFO["UQM"]
        self.NEB = INFO["NEB"]
        self.BID = INFO["BID"]
        self.MTA = INFO["MTA"]
        self.TRC = INFO["TRC"]
        self.NMD = INFO["NMD"]
        self.IPG = INFO["IPG"]
        self.LSG = INFO["LSG"]
        self.ATI = INFO["ATI"]
        self.ATD = INFO["ATD"]

    ################
    # FORMAT field #
    ################
    @property
    def FORMAT(self):
        ad = ",".join([str(int(i)) for i in self.AD])

        return "AD\t" + ad

    @FORMAT.setter
    def FORMAT(self, FORMAT):
        self.AD = FORMAT["AD"]
