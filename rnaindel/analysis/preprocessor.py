import sys
sys.path.append("/research/rgs01/project_space/zhanggrp/MethodDevelopment/common/indeltools")

import csv
import pysam
import pandas as pd

from indeltools import Variant, VariantAlignment

CANONICALS = [str(i) for i in range(1, 23)] + ["X", "Y"]



def preprocess(callset, reference, bam, dbsnp, from_default_caller=True):
    
    if from_default_caller:
        f = open(callset)
        records = csv.DictReader(f, delimiter="\t")
       
        data = [instantiate_defalt_caller_record(record, reference, bam, dbsnp) for record in records]
    
    data = [i for i in data if i]      
    return pd.DataFrame(data)


def instantiate_defalt_caller_record(record, reference, bam, dbsnp):

    #TODO check fasta BAM chrnames
    chrom = record["Chr"].replace("chr", "")
    
    if not chrom.replace("chr", "") in CANONICALS: 
        return None

    pos = int(record["Pos"])
    ref = record["Chr_Allele"]
    alt = record["Alternative_Allele"]
    var_type = record["Type"]
    
    if var_type == "SNP":
        return None

    pos -= 1
    padding_base = reference.fetch(chrom, pos - 1, pos)  
    
    if var_type == "deletion":
        alt = padding_base
        ref = alt + ref
    else:
        ref = padding_base
        alt = ref + alt
    
    var = Variant(chrom, pos, ref, alt, reference).normalize()
    
    try:
        varaln = VariantAlignment(var, bam)
        cnt = varaln.count_alleles(by_fragment=True)
        ref_cnt, alt_cnt = cnt[0], cnt[1]
        su = varaln.to_complex(dbsnp=dbsnp)
        if su:
            su_pos, su_ref, su_alt = su.pos, su.ref, su.alt
        else:
            su_pos, su_ref, su_alt = -1, -1, -1
    except:
        ref_cnt, alt_cnt = -99, -99
        su_pos, su_ref, su_alt = -99, -99, -99

    return {"CHROM": chrom, "POS": pos, "ALT": alt, "REF": ref, "REF_CNT": ref_cnt, "ALT_CNT": alt_cnt, "CPOS": su_pos, "CREF": su_ref, "CALT": su_alt} 
   
def instantiate_vcf_record(chrom, pos, ref, alt, reference, allele_len_thresh=100):
    if len(ref) == len(alt):
        return None

    if not chrom.replace("chr", "") in CANONICALS:
        return None
    
    if len(ref) > allele_len_thresh or len(alt) > allele_len_thresh:
        return None
    
    return Variant(chrom, pos, ref, alt, reference).normalize()

    

def is_canonical_indel(record):
    is_indel = (len(record["REF"]) != len(record["ALT"]))  
    
    chrom_name = record["CHROM"].replace("chr", "")
    
    is_canonical = chrom_name in CANONICALS



    if is_default:
        var_type = record["Type"]
        is_indel = (var_type == "insertion") or (var_type == "deletion")
    
        chrom_name = record["Chr"].replace("chr", "")
    else:
        is_indel = (len(record["REF"]) != len(record["ALT"]))
        
        chrom_name = record["CHROM"].replace("chr", "")
        
    is_canonical = chrom_name in CANONICALS

    


def to_vcf(record_from_default_caller, reference):
    pass    
    
    
