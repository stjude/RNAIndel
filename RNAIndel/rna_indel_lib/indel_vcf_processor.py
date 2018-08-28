#!/usr/bin/env python3
"""Optional step for .vcf input

Convert .vcf to a format compatible with Bambino

'vcf_processor' is the main routine of this module
"""

import sys
import vcf
import logging
import pandas as pd
from .indel_snp_annotator_dev import count_padding_bases

logger = logging.getLogger(__name__)

def vcf_processor(vcffile):
    """
    
    Args:
        vcffile (str): vcf file name
    Returns:
        df (pandas.DataFrame): df with indels reported as in Bambino output
    """
    df = pd.DataFrame(columns=['chr', 'pos', 'ref', 'alt'], index=[0])
    
    vcf_reader = vcf.Reader(open(vcffile, 'r'))
    
    i = 0
    for record in vcf_reader:
        parsed_tuple = parse_vcf_record(record)
        if parsed_tuple:
            df.ix[i,'chr'] = parsed_tuple[0]
            df.ix[i,'pos'] = parsed_tuple[1]
            df.ix[i,'ref'] = parsed_tuple[2]
            df.ix[i,'alt'] = parsed_tuple[3] 
            i += 1
     
    datasize = len(df)
    if len(df) ==  0:
        logging.warning('No indels detected in your vcf. Analysis done.')
        sys.exit(0)

    return df


def select_canonical_chromosome(chr):
    """Extract indels on chr1-22, X and Y

    Args:
        chr (str): may not be prefixed with 'chr'
    Returns:
        chromosome (str): prefixed with 'chr'
    """
    chromosome = None

    if chr.startswith('chr'):
        chr = chr.replace('chr', '')
    
    if chr == 'X' or chr == 'Y':
        chromosome = 'chr' + chr
    else:
        try:
            if 1 <= int(chr) <= 22:
                chromosome = 'chr' + chr
        except:
            pass
    
    return chromosome


def parse_vcf_record(record):
    """
    Args:
        record (vcf.Record): record in .vcf 
    Returns:
        parsed_record (tuple): (chr, pos, ref, alt) 
                               Bambino compatible format
    Example:
      deletion
              pos 123456789012
        reference ATTAGTAGATGT
        deletion  ATTA---GATGT
         
        VCF:
            CHROM POS REF  ALT
            N     4   AGTA A
        Bambino:
            chr   pos ref alt
            chr_N 5   GTA -
         
      insertion
           pos 1234***56789012
        reference ATTA***GTAGATGT
        insertion ATTAGTAGTAGATGT
         
        VCF:
            CHROM POS REF ALT
            N     4   A   AGTA
        Bambino:
            chr   pos ref alt
            chr_N 5   -   GTA   
    """
    parsed_record = None

    if record.is_indel:
        chr = record.CHROM
        chr = select_canonical_chromosome(chr)
        if chr: 
            pos = record.POS + 1  # converts to Bambino coordinate
            ref = record.REF
            alts = record.ALT
            for alt in alts:
                alt = str(alt)  # cast to str from _Substitution obj
                n = count_padding_bases(ref, alt)
                if len(ref) < len(alt):
                    ref = '-'
                    alt = alt[n:]
                else:
                    ref = ref[n:]
                    alt = '-'
           
                parsed_record = (chr, pos, ref, alt)

                return parsed_record
