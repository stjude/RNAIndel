#!/usr/bin/env python3
"""3rd step of analysis.

Calculates features at sequence and alignment level

'indel_sequence_processor' is the main routine of this module
"""

from .most_common import most_common
from .indel_features import SamFeatures
from .indel_features import AnnotationFeatures
from .indel_curator import curate_indel_in_genome
from .indel_curator import curate_indel_in_pileup


def indel_sequence_processor(df, genome, alignments, mapq, chr_prefixed):
    """Calculate features from Bambino output, annotation, and .bam
    
    Features not used for final model are commented out '#'
   
    Args:
        df (pandas.DataFrame)
        genome (pysam.FastaFile): reference genome
        alignments (pysam.AlignmentFile): bam data
        mapq (int): MAPQ score for uniquely mapped reads
        chr_prefixed (bool): True if chromosome names are "chr"-prefixed
    Returns:
        df (pandas.DataFrame): dataframe with valid entries
        df_filtered_premerge (pandas.DataFrame): dataframe with invalid entries
    """
    # features derived from call set
    df["is_gc_ins"] = df.apply(is_gc_ins, axis=1)
    df["is_gc_del"] = df.apply(is_gc_del, axis=1)
    df["is_at_ins"] = df.apply(is_at_ins, axis=1)
    df["is_at_del"] = df.apply(is_at_del, axis=1)
    df["indel_size"] = df.apply(indel_size, axis=1)

    # features derived from annotation
    df["a"] = df.apply(anno_features, axis=1)
    df["is_inframe"] = df.apply(lambda x: x["a"].is_inframe, axis=1)
    df["is_truncating"] = df.apply(lambda x: x["a"].is_truncating, axis=1)
    df["is_splice"] = df.apply(lambda x: x["a"].is_splice, axis=1)
    df["is_nmd_insensitive"] = df.apply(lambda x: x["a"].is_nmd_insensitive, axis=1)

    # features derived from sequence and alingment
    df["s"] = df.apply(
        sam_features,
        genome=genome,
        alignments=alignments,
        mapq=mapq,
        chr_prefixed=chr_prefixed,
        axis=1,
    )
    df["gc"] = df.apply(lambda x: x["s"].gc, axis=1)
    df["local_gc"] = df.apply(lambda x: x["s"].local_gc, axis=1)
    df["lc"] = df.apply(lambda x: x["s"].lc, axis=1)
    df["local_lc"] = df.apply(lambda x: x["s"].local_lc, axis=1)
    df["strength"] = df.apply(lambda x: x["s"].strength, axis=1)
    df["local_strength"] = df.apply(lambda x: x["s"].local_strength, axis=1)
    df["repeat"] = df.apply(lambda x: x["s"].repeat, axis=1)
    df["dissimilarity"] = df.apply(lambda x: x["s"].dissimilarity, axis=1)
    df["indel_complexity"] = df.apply(lambda x: x["s"].indel_complexity, axis=1)
    df["ref_count"] = df.apply(lambda x: x["s"].ref_count, axis=1)
    df["alt_count"] = df.apply(lambda x: x["s"].alt_count, axis=1)
    df["is_multiallelic"] = df.apply(lambda x: x["s"].is_multiallelic, axis=1)
    df["is_near_boundary"] = df.apply(lambda x: x["s"].is_near_boundary, axis=1)
    df["is_bidirectional"] = df.apply(lambda x: x["s"].is_bidirectional, axis=1)
    df["is_uniq_mapped"] = df.apply(lambda x: x["s"].is_uniq_mapped, axis=1)

    df.drop(["a", "s"], axis=1, inplace=True)

    df["filtered"] = df.apply(flag_invalid_entry, axis=1)

    df, df_filtered_premerge = df[df["filtered"] == "-"], df[df["filtered"] != "-"]

    # drop original calls rescued by equivalence
    df.dropna(inplace=True)

    return df, df_filtered_premerge


def is_gc_ins(row):
    """Encodes if the indel is an insertion of 'G' or 'C'
    
    Args:
        row (pandas.Series): a Series with 'is_ins'(bool) 
                             and 'indel_seq'(str) indexes 
    Returns:
        is_gc_ins (bool): 1 for 'G' or 'C' insertion
                          0 otherwise
    """
    if row["is_ins"] == 1 and row["indel_seq"] == "G":
        is_gc_ins = 1
    elif row["is_ins"] == 1 and row["indel_seq"] == "C":
        is_gc_ins = 1
    else:
        is_gc_ins = 0

    return is_gc_ins


def is_gc_del(row):
    """Encodes if the indel is a deletion of 'G' or 'C'
    
    Args:
        row (pandas.Series): a Series with 'is_ins'(bool)
                             and 'indel_seq'(str) indexes
    Returns:
        is_gc_del (bool): 1 for 'G' or 'C' deletion
                          0 otherwise
    """
    if row["is_ins"] == 0 and row["indel_seq"] == "G":
        is_gc_del = 1
    elif row["is_ins"] == 0 and row["indel_seq"] == "C":
        is_gc_del = 1
    else:
        is_gc_del = 0

    return is_gc_del


def is_at_ins(row):
    """Encodes if the indel is an insertion of 'A' or 'T'
    
    Args:
        row (pandas.Series): a Series with 'is_ins'(bool)
                             and 'indel_seq'(str) indexes

    Returns:
        is_at_ins (bool): 1 for 'A' or 'T' insertion
                          0 otherwise
    """
    if row["is_ins"] == 1 and row["indel_seq"] == "A":
        is_at_ins = 1
    elif row["is_ins"] == 1 and row["indel_seq"] == "T":
        is_at_ins = 1
    else:
        is_at_ins = 0

    return is_at_ins


def is_at_del(row):
    """Encodes if the indel is an deletion with 'A' or 'T'
    
    Args:
        row (pandas.Series): a Series with 'is_ins'(bool)
                             and 'indel_seq'(str) indexes

    Returns:
        is_at_del (bool): 1 for 'A' or 'T' deletion
                          0 otherwise
    """
    if row["is_ins"] == 0 and row["indel_seq"] == "A":
        is_at_del = 1
    elif row["is_ins"] == 0 and row["indel_seq"] == "T":
        is_at_del = 1
    else:
        is_at_del = 0

    return is_at_del


def indel_size(row):
    """Calculate the length of indel sequence.

    Args:
        row (pandas.Series): a Series with 'indel_seq'(str) index
    Returns:
        indel_size (int): a positive int
    """
    indel_size = len(row["indel_seq"])

    return indel_size


def anno_features(row):
    """Encodes features derived from variant annotaion:
        
    1. inframe 
    2. truncating (frameshift, nonsense, or splice-motif affected)
    3. whether it is in splice site/region
    4. whether it is in the first or last exon
         
    In multiple-isoform case, the common pattern is returned 
    except for is_inframe and is_splice

    Args:
        row (pandas.Series): a Series with 'annotation' index
    Returns:
        AnnotationFeatures (class): class to store variant annotation features 
    """
    lst = row["annotation"].split(",")

    inframe = 0
    truncates = []
    splice = 0
    insensitivities = []
    for anno in lst:

        if "inframe" in anno:
            inframe += 1
        else:
            pass

        if "Truncating" in anno:
            truncates.append(1)
        else:
            truncates.append(0)

        if "splice" in anno:
            splice += 1
        else:
            pass

        is_insensitive = int(anno.split("|")[-1])
        insensitivities.append(is_insensitive)

    is_inframe = 0
    if inframe > 0:
        is_inframe = 1

    is_splice = 0
    if splice > 0:
        is_splice = 1

    return AnnotationFeatures(
        is_inframe, most_common(truncates), is_splice, most_common(insensitivities)
    )


def sam_features(row, genome, alignments, mapq, chr_prefixed):
    """Encodes features derived from sequence alignment/map(SAM)
    
    Args:
        row (pandas.Series): a Series with 'chr', 'pos', 
                             'is_ins', 'indel_seq' indexes
        genome (pysam.FastaFile): reference genome 
        alignments (pysam.AlignmentFile): bam data
        mapq (int): MAPQ score for unique mappers
        chr_prefixed (bool): True if chromosome names in BAM are "chr"-prefixed
    Returns:
        SamFeatures (class): class to store sequence and alignment features            
    """
    dna_window = 50
    rna_window = 6

    chr = row["chr"]  # this is "chr"-prefixed
    pos = row["pos"]
    idl_type = row["is_ins"]
    idl_seq = row["indel_seq"]

    # SequenceWithIndel obj in refrence genome
    idl_ref_genome = curate_indel_in_genome(
        genome, chr, pos, idl_type, idl_seq, chr_prefixed
    )
    # PileupWithIndel obj in bam
    idl_bam = curate_indel_in_pileup(
        alignments, chr, pos, idl_type, idl_seq, mapq, chr_prefixed
    )

    # global sequence properties
    # derived from reference genome
    gc = idl_ref_genome.gc(dna_window)
    lc = idl_ref_genome.lc(dna_window)
    strength = idl_ref_genome.strength(dna_window)

    # local sequence properties derived from bam
    # these consider individual variations such SNPs
    # replace with info from fasta if failed to retrieve
    # info from bam (this may happen if the reads are too short)
    try:
        local_gc = idl_bam.gc(rna_window)
    except:
        local_gc = idl_ref_genome.gc(rna_window)

    try:
        local_lc = idl_bam.local_lc(rna_window)
    except:
        local_lc = idl_ref_genome.local_lc(rna_window)

    try:
        local_strength = idl_bam.strength(rna_window)
    except:
        local_strength = idl_ref_genome.strength(rna_window)

    try:
        repeat = idl_bam.repeat()
    except:
        repeat = idl_ref_genome.repeat()

    try:
        dissimilarity = idl_bam.dissimilarity()
    except:
        dissimilarity = idl_ref_genome.dissimilarity()

    try:
        indel_complexity = idl_bam.indel_complexity(rna_window)
    except:
        indel_complexity = 0

    # alignment/map properities
    try:
        ref_count = idl_bam.ref_count
    except:
        ref_count = None

    try:
        alt_count = idl_bam.alt_count
    except:
        alt_count = None

    try:
        is_multiallelic = idl_bam.is_multiallelic
    except:
        is_multiallelic = 0

    try:
        is_near_boundary = idl_bam.is_near_boundary
    except:
        is_near_boundary = 0

    try:
        is_bidirectional = idl_bam.is_bidirectional
    except:
        is_bidirectional = 1

    try:
        is_uniq_mapped = idl_bam.is_uniq_mapped
    except:
        is_uniq_mapped = 0

    return SamFeatures(
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
    )


def flag_invalid_entry(row):
    """Flag entries to be filtered
   
    Args:
        row (pandas.Series)
    Returns:
        filtered (str): '-' for valid entrr, otherwise invalid.

    Example:
        
       Case 1:
            chr  pos  rescued         ref_count alt_count
            chrN 123  by_equivalence  30        1
         ->      
            chr  pos  rescued         ref_count alt_count filtered
            chrN 123  by_equivalence  30        1         lt2count

       Case 2:
            chr  pos  rescued ref_count alt_count
            chrN 123  - 
       ->
            chr  pos  rescued ref_count alt_count filtered
            chrN 123  -                           notfound
            
       Case 3:
       The indel at chrN:123 was rescued by the nearest indel chrN:125.
            chr  pos  rescued               ref_count  alt_count
            chrN 123  rescued_by:chrN:125
            chrN 125                        10         5
       
       ->   chr  pos  rescued               ref_count  alt_count filtereed
            chrN 123  rescued_by:chrN:125                        by_nearest
            chrN 125                        10         5         -   
    """
    filtered = "-"

    # Case 1
    if row["alt_count"] < 2 and row["rescued"] != "by_equivalence":
        filtered = "lt2count"
    # Case 2
    # not rescued and not found in the bam
    elif row["rescued"] == "-" and row["alt_count"] != row["alt_count"]:
        filtered = "notfound"
    # Case 3
    # original call that is rescued by nearest indel
    # (they are not found as specified)
    elif (
        row["rescued"].startswith("rescued_by:")
        and row["alt_count"] != row["alt_count"]
    ):
        filtered = "by_nearest"
    else:
        pass

    return filtered
