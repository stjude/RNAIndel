#!/usr/bin/env python3

import re
from .indel_sequence import Indel
from .left_aligner import peek_left_base
from .left_aligner import lt_aln

snp_ptn = re.compile(r"rs[0-9]+")


class IndelVcfReport(object):
    """Make VCF record and meta info of indel from Bambino coordinate

    Attributes:
        genome (pysam.FastaFile): reference genome
        chr (str): chr1-22, chrX, chrY
        pos (int): 1-based indel pos (Bambino coordinate)
        ref (str): Bambino style ref allele
        alt (str): Bambino style alt allele
        chr_prefixed (bool): True if chromosome names in BAM are "chr"-prefixed
    """

    def __init__(self, genome, chr, pos, ref, alt, chr_prefixed):
        self.genome = genome
        self.chr = chr
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.chr_prefixed = chr_prefixed

    def generate_indel(self):
        if self.ref == "-":
            idl = Indel(self.chr, self.pos, 1, self.alt)
        else:
            idl = Indel(self.chr, self.pos, 0, self.ref)

        return idl

    def left_align(self):
        idl = lt_aln(self.generate_indel(), self.genome, self.chr_prefixed)
        self.pos = idl.pos
        if self.ref == "-":
            self.alt = idl.idl_seq
        else:
            self.ref = idl.idl_seq
    
    @property
    def CHROM(self):
        return self.chr

    @property
    def POS(self):
        return self.pos - 1

    @property
    def ID(self):
        """Return ID info (dbSNP ID).
        If info is none or not defined (not set), return the 
        defined missing value: a dot '.'

        Returns:
            self.__ID (str) or . (str)
        """
        if self.__ID and self.__ID.startswith("rs"):
            return self.__ID
        else:
            return "."

    @ID.setter
    def ID(self, ID):
        """Set dbSNP IDs as semi-colon delimited list
        
        Args:
            ID (str): dbSNP IDs with any delimiter
        """
        self.__ID = ID
        if self.__ID:
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
            return peek_left_base(self.generate_indel(), self.genome, self.chr_prefixed)
        else:
            return (
                peek_left_base(self.generate_indel(), self.genome, self.chr_prefixed)
                + self.ref
            )

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
            self.QUAL (int) or . (str)
        """
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
        if self.__FILTER == "-":
            return "PASS"
        elif self.__FILTER == "notfound":
            return "NtF"
        elif self.__FILTER == "lt2count":
            return "Lt2"
        elif self.__FILTER == "by_nearest":
            return "RqN"
        else:
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
        if self.chr_prefixed:
            chr = self.CHROM
        else:
            chr = self.CHROM.replace("chr", "")

        record = [
            chr,
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
        if not self.PRED:
            pred = ""
        else:
            pred = "PRED=" + self.PRED

        if not self.PROB:
            prob = ""
        else:
            prob = "PROB=" + ",".join([str(p) for p in self.PROB])

        anno = "ANNO=" + self.ANNO

        if not self.MAXMAF or self.MAXMAF == -1:
            maxmaf = ""
        else:
            maxmaf = "MAXMAF=" + str(self.MAXMAF)

        common = ""
        if self.COMMON == 1:
            common = "COMMON"

        if not self.CLIN or self.CLIN == "-":
            clin = ""
        else:
            clin = "CLIN=" + self.CLIN

        if not self.REP:
            rep = ""
        else:
            rep = "REP=" + str(self.REP)
        
        if not self.LC:
            lc = ""
        else:
            lc = "LC=" + str(self.LC)
        
        if not self.LLC:
            llc = ""
        else:
            llc = "LLC=" + str(self.LLC)
        
        if not self.GC:
            gc = ""
        else:
            gc = "GC=" + str(self.GC)

        if not self.LGC:
            lgc = ""
        else:
            lgc = "LGC=" + str(self.LGC)
        
        if not self.SG:
            sg = ""
        else:
            sg = "SG=" + str(self.SG)
       
        if not self.LSG:
            lsg = ""
        else:
            lsg = "LSG=" + str(self.LSG)
       
        if not self.DSM:
            dsm = ""
        else:
            dsm = "DSM=" + str(self.DSM)
        
        if not self.ICP:
            icp = ""
        else:
            icp = "ICP=" + str(self.ICP)
        
        if not self.ISZ:
            isz = ""
        else:
            isz = "ISZ=" + str(self.ISZ)
        
        ins = "INS"
        if not self.INS or self.INS == 0:
            ins = ""

        ati = "ATI"
        if not self.ATI or self.ATI == 0:
            ati = ""

        atd = "ATD"
        if not self.ATD or self.ATD == 0:
            atd = ""

        gci = "GCI"
        if not self.GCI or self.GCI == 0:
            gci = ""

        gcd = "GCD"
        if not self.GCD or self.GCD == 0:
            gcd = ""

        if not self.REFC:
            refc = ""
        else:
            refc = "REFC=" + str(self.REFC)

        if not self.ALTC:
            altc = ""
        else:
            altc = "ALTC=" + str(self.ALTC)
        
        bid = "BID"
        if not self.BID or self.BID == 0:
            bid = ""

        uqm = "UQM"
        if not self.UQM or self.UQM == 0:
            uqm = ""

        neb = "NEB"
        if not self.NEB or self.NEB == 0:
            neb = ""

        eqx = "EQX"
        if not self.EQX or self.EQX == 0:
            eqx = ""

        mta = "MTA"
        if not self.MTA or self.MTA == 0:
            mta = ""
        
        frm = "FRM"
        if not self.FRM or self.FRM == 0:
            frm = ""
        
        spl = "SPL"
        if not self.SPL or self.SPL == 0:
            spl = ""
        
        trn = "TRN"
        if not self.TRN or self.TRN == 0:
            trn = ""
        
        cdd = "CDD"
        if not self.CDD or self.CDD == 0:
            cdd = ""

        if not self.LOC:
            loc = ""
        else:
            loc = "LOC=" + str(self.LOC)

        nmd = "NMD"
        if not self.NMD or self.NMD == 0:
            nmd = ""

        if not self.IPG:
            ipg = ""
        else:
            ipg = "IPG=" + str(self.IPG)

        if not self.LEN:
            len = ""
        else:
            len = "LEN=" + str(self.LEN)

        dbsnp = "DBSNP"
        if not self.DBSNP or self.DBSNP == 0:
            dbsnp = ""

        rcf = "RCF"
        if not self.RCF or self.RCF == "-":
            rcf = ""

        if self.RQB and self.RQB[0] == "by_nearest":
            rqb = "RQB=" + self.RQB[1].replace("rescued_by:", "")
        else:
            rqb = ""

        info_lst = [
            pred,
            prob,
            anno,
            maxmaf,
            common,
            clin,
            rep,
            lc,
            llc,
            gc,
            lgc,
            sg,
            lsg,
            dsm,
            icp,
            isz,
            ins,
            ati,
            atd,
            gci,
            gcd,
            refc,
            altc,
            bid,
            uqm,
            neb,
            eqx,
            mta,
            frm,
            spl,
            trn,
            cdd,
            loc,
            nmd,
            ipg,
            len,
            dbsnp,
            rcf,
            rqb,
        ]

        return ";".join([i for i in info_lst if i != ""])

    @INFO.setter
    def INFO(self, INFO):
        self.PRED = INFO["PRED"]
        self.PROB = INFO["PROB"]
        self.ANNO = INFO["ANNO"]
        self.MAXMAF = INFO["MAXMAF"]
        self.COMMON = INFO["COMMON"]
        self.CLIN = INFO["CLIN"]
        self.REP = INFO["REP"]
        self.LC = INFO["LC"]
        self.LLC = INFO["LLC"]
        self.GC = INFO["GC"]
        self.LGC = INFO["LGC"]
        self.SG = INFO["SG"]
        self.LSG = INFO["LSG"]
        self.DSM = INFO["DSM"]
        self.ICP = INFO["ICP"]
        self.ISZ = INFO["ISZ"]
        self.INS = INFO["INS"]
        self.ATI = INFO["ATI"]
        self.ATD = INFO["ATD"]
        self.GCI = INFO["GCI"]
        self.GCD = INFO["GCD"]
        self.REFC = INFO["REFC"]
        self.ALTC = INFO["ALTC"]
        self.BID = INFO["BID"]
        self.UQM = INFO["UQM"]
        self.NEB = INFO["NEB"]
        self.EQX = INFO["EQX"]
        self.MTA = INFO["MTA"]
        self.FRM = INFO["FRM"]
        self.SPL = INFO["SPL"]
        self.TRN = INFO["TRN"]
        self.CDD = INFO["CDD"]
        self.LOC = INFO["LOC"]
        self.NMD = INFO["NMD"]
        self.IPG = INFO["IPG"]
        self.LEN = INFO["LEN"]
        self.DBSNP= INFO["DBSNP"]
        self.RCF = INFO["RCF"]
        self.RQB = INFO["RQB"]

    ################
    # FORMAT field #
    ################
    @property
    def FORMAT(self):
        if self.AD:
            ad = ",".join([str(int(i)) for i in self.AD])
            return "AD\t" + ad
        else:
            return ""

    @FORMAT.setter
    def FORMAT(self, FORMAT):
        self.AD = FORMAT["AD"]
