#!/usr/bin/env python3

from indel_sequence_dev import Indel


class AnnotationFeatures(object):

    def __init__(self, is_inframe, is_truncating,\
                       is_splice, is_nmd_insensitive):

        self.is_inframe = is_inframe
        self.is_truncating = is_truncating
        self.is_splice = is_splice
        self.is_nmd_insensitive = is_nmd_insensitive


class SamFeatures(object):

    def __init__(self, gc, lc, strength,\
                 local_gc, local_lc, local_strength,\
                 repeat, dissimilarity, indel_complexity,\
                 ref_count, alt_count,\
                 is_multiallelic, is_near_boundary,\
                 is_bidirectional, is_uniq_mapped):

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


class IndelReport(Indel):
    
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


    def list_dbsnp_id(self, rs):
        self.dbsnp_id.append(rs)


    def report_dbsnp_id(self):
        if self.dbsnp_id == []:
            return '-'
        else:
            return ','.join(self.dbsnp_id)

    
    def list_clnvr_id(self, id):
        self.clnvr_id.append(id)


    def report_clnvr_id(self):
        if self.clnvr_id == []:
            return '-'
        else:
            return ','.join(self.clnvr_id)
    
        
    def list_dbsnp_freq(self, freq):
        self.dbsnp_freq.append(freq)


    def list_clnvr_freq(self, freq):
        self.clnvr_freq.append(freq)


    def report_freq(self):
        freqs = self.dbsnp_freq + self.clnvr_freq
       
        if freqs == []: 
            return -1
        else:
            return max(freqs)   
        
     
    def list_dbsnp_common(self, common):
        self.dbsnp_common.append(common)


    def is_common(self):
        if self.dbsnp_common == []:
            return -1
        elif 1 in self.dbsnp_common:
            return 1
        elif self.report_freq() != '-':
            if self.report_freq() >= 0.01:
                return 1
            else:
                return 0
        else:
            return 0
        
       
    def list_dbsnp_origin(self, origin):
        if origin == 0:
            ori = 'unspecified'
        elif origin == 1:
            ori = 'germline'
        elif origin == 2:
            ori = 'somatic'
        elif origin == 3:
            ori = 'both'
        else:
            ori = '-'

        self.dbsnp_origin.append(ori)
        

    def list_clnvr_origin(self, origin):
        if origin == 0:
            ori = 'unknown'
        elif origin == 1:
            ori = 'germline'
        elif origin == 2:
            ori = 'somatic'
        # 4 = 'inherited'
        elif origin == 4:
            ori = 'germline'
        # 8 = 'paternal'
        elif origin == 8:
            ori = 'germline'
        # 16 = 'maternal'
        elif origin == 16:
            ori = 'germline'
        elif origin == 32:
            ori = 'de-novo'
        # 64 = 'biparental'
        elif origin == 64:
            ori = 'germline'
        # 128 = 'uniparental'
        elif origin == 128:
            ori = 'germline'
        elif origin == 256:
            ori = 'not-tested'
        elif origin == 512:
            ori = 'tested-inconclusive'
        elif origin == 1073741824:
            ori = 'others'
        else:
            ori = '-'
             
        self.clnvr_origin.append(ori)        


    def with_germline_reports(self):
        origins = self.dbsnp_origin + self.clnvr_origin
        if origins == []:
            return -1
        elif 'germline' in origins and not 'somatic' in origins:
            return 1
        elif 'germline' in origins and not 'both' in origins:
            return 1
        else:
            return 0


    def list_clnvr_info(self, clninfo):
        self.clnvr_info.append(clninfo)
    

    def report_clnvr_info(self):
        if self.clnvr_info == []:
            return '-'
        else:
            return ','.join(self.clnvr_info)
    

    def is_not_pathogenic(self):
        info = self.report_clnvr_info()
        
        if info == '-':
            return -1
        
        info = info.lower()
        if 'pathogenic' not in info:
            if 'benign' in info:
                return 1
            elif 'protective' in info:
                return 1
            elif 'affects' in info:
                return 1
            else:
                return 0
        else:
            return 0  

