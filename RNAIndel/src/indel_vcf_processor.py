#!/usr/bin/env python3

import sys
import vcf
import logging
import pandas as pd
from indel_snp_annotator_dev import count_padding_bases

logger = logging.getLogger(__name__)

def main(vcffile):
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
