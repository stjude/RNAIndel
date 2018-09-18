import pysam
import pandas as pd
from functools import partial
from multiprocessing import Pool
from .most_common import most_common
from .indel_curator import extract_indel_reads
from .indel_curator import decompose_indel_read
from .indel_curator import curate_indel_in_genome
from .indel_curator import extract_all_valid_reads
from .indel_equivalence_solver import are_equivalent


def equivalence_collector(df, fasta, bam, num_of_process, left_alinged): 
    #bam_data = pysam.AlignmentFile(bam, 'rb')
    
    pool = Pool(num_of_process)

    # wrap recover_equivalent_indels
    recover = partial(recover_equivalent_indels,
                      fasta=fasta,
                      bam=bam,
                      search_window=50,
                      pool=pool,
                      left_aligned=left_alinged)
    
    df['recovered_indels'] = df.apply(recover, axis=1)
    
    data_in_list_of_dict = df['recovered_indels'].sum()
    df_recovered = pd.DataFrame(data_in_list_of_dict)
    df = df[['chr', 'pos', 'ref', 'alt']]
    df = pd.concat([df, df_recovered], axis=0, sort=True)
    
    df = sort_positionally(df)
    df = df.drop_duplicates(['chr', 'pos', 'ref', 'alt'])     
    df.reset_index(drop=True, inplace=True) 
    
    return df
                          
   
def recover_equivalent_indels(row, fasta, bam, search_window, pool, left_aligned=True):
    """Recorver equivalent indels from left-aligned indel report
    
    Args:
        row (pandas.Series)
        fasta (str): path to fasta
        bam_data (pysam.AlignmentFile)
        search_window (int): to define search range
    Returns
        equivalents (list): dict element
                            {'chr':chromosome,
                             'pos':1-based pos,
                             'ref':ref allele,
                             'alt':alt allele}

                            empty list if no equivalent indels found
    """
    chr = row['chr']         
    pos = row['pos']
    
    if row['ref'] == '-':
        idl_type, idl_seq = 1, row['alt']
    else:
        idl_type, idl_seq = 0, row['ref']
    
    called_idl = curate_indel_in_genome(fasta, chr, pos, idl_type, idl_seq)
    
    if left_aligned:
        rt_window, lt_window = search_window, search_window - search_window
    else:
        rt_window, lt_window = int(search_window/2), int(search_window/2)
     

    eq_idl = partial(extract_equivalent_indel,
                     called_idl=called_idl,
                     fasta=fasta,
                     bam=bam,
                     chr=chr,
                     idl_type=idl_type)
    
    
    rt_range = [pos+i for i in range(rt_window)]
    rt_equivalents = pool.map(eq_idl, rt_range)
      
    lt_range = [pos-i for i in range(lt_window)] 
    lt_equivalents = pool.map(eq_idl, lt_range)
    
    equivalents = list(rt_equivalents) + list(lt_equivalents)

    equivalents = [{'chr':eq.chr,
                    'pos':eq.pos,
                    'ref':eq.ref,
                    'alt':eq.alt} for eq in equivalents if eq != None]
    
    return equivalents


def extract_equivalent_indel(pos, called_idl, fasta, bam, chr, idl_type):
    """Extract equivalent indel if exists at the locus (chr, pos)

    Args:
        called_idl (SequenceWithIndel): indel object to compare
        fasta (str): path to fasta 
        bam_data (pysam.AlignmentFile)
        pos (int): 1-based coordinate
        idl_type (int): 1 for insertion, 0 for deletion
    Returns:
        equivalent_indel(SequenceWithIndel or None) :SequenceWithIndel if found
                                                     None if not found
    """
    equivalent_indel = None

    bam_data = bam_data = pysam.AlignmentFile(bam, 'rb')
    inferred_idl_seq = get_most_common_indel_seq(bam_data, chr, pos, idl_type)
    
    if inferred_idl_seq:
        idl_at_this_locus = curate_indel_in_genome(fasta, chr, pos, idl_type, inferred_idl_seq)
        
        if are_equivalent(called_idl, idl_at_this_locus):
            equivalent_indel = idl_at_this_locus
    
    return equivalent_indel

    
def get_most_common_indel_seq(bam_data, chr, pos, idl_type):
    """Extract most frequent indel sequnece from bam data

    Args:
        bam_data (pysam.pysam.AlignmentFile)
        chr (str): chr1-22, chrX or chrY
        pos (int): 1-based coordinate
        idl_type (int): 1 for insertion, 0 for deletion
    Returns:
        idl_seq (str or None): None type if no indels found
    """
    idl_seq = None
    
    if idl_type == 1:
        ins_or_del = 'I'
    else:
        ins_or_del = 'D'
    
    # convert 0-based coordinate
    pos = pos - 1 
    
    try:
        valid_reads = extract_all_valid_reads(bam_data, chr, pos)
    except:
        return idl_seq

    try:
        parsed_indel_reads = extract_indel_reads(valid_reads, pos, ins_or_del)
    except:
        return idl_seq

    if parsed_indel_reads == []:
        return idl_seq
    else:
        decomposed = [decompose_indel_read(parsed_read) \
                      for parsed_read in parsed_indel_reads] 
    
    idl_seq = most_common([decomp[1] for decomp in decomposed]) 
   
    return idl_seq

    
def sort_positionally(df):
    df['chr'] = df.apply(lambda x:x['chr'].replace('chr', ''), axis=1)
    df['chr'] = df.apply(lambda x:23 if x['chr'] == 'X' else x['chr'], axis=1)
    df['chr'] = df.apply(lambda x:24 if x['chr'] == 'Y' else x['chr'], axis=1)
    df['chr'] = df.apply(lambda x:int(x['chr']), axis=1)
    
    df.sort_values(['chr', 'pos'], inplace=True)

    df['chr'] = df.apply(lambda x:'Y' if x['chr'] == 24 else x['chr'], axis=1)
    df['chr'] = df.apply(lambda x:'X' if x['chr'] == 23 else x['chr'], axis=1)
    df['chr'] = df.apply(lambda x:'chr'+str(x['chr']), axis=1)
    
    return df


