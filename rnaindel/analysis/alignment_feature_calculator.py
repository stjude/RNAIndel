import re
import pysam
import random
from indelpost import Variant, VariantAlignment

from .utils import most_common, split, flatten_list_of_list
from .sequence_properties import (
    repeat,
    dissimilarity,
    editdistance,
    linguistic_complexity,
    gc,
    dna_strength,
)

cigar_ptrn = re.compile(r"[0-9]+[MIDNSHPX=]")


def alignment_features(df, bam, mapq, downsample_threshold=1000):
    bam = pysam.AlignmentFile(bam)

    (
        df["repeat"],
        df["lc"],
        df["local_lc"],
        df["gc"],
        df["local_gc"],
        df["strength"],
        df["local_strength"],
        df["dissimilarity"],
        df["indel_complexity"],
        df["indel_size"],
        df["is_ins"],
        df["is_at_ins"],
        df["is_at_del"],
        df["is_gc_ins"],
        df["is_gc_del"],
        df["ref_count"],
        df["alt_count"],
        df["orig_ref_cnt"],
        df["orig_alt_cnt"],
        df["is_bidirectional"],
        df["is_uniq_mapped"],
        df["uniq_mapping_rate"],
        df["is_near_boundary"],
        df["equivalence_exists"],
        df["is_multiallelic"],
        df["cplx_variant"]
    ) = zip(*df.apply(_wrapper, bam=bam, mapq=mapq, downsample_threshold=downsample_threshold, axis=1))

    df = df[df["alt_count"] > 1]

    return df


def _wrapper(row, bam, mapq, downsample_threshold):
    variant = row["indel"]
    (
        indel_size,
        is_ins,
        is_at_ins,
        is_at_del,
        is_gc_ins,
        is_gc_del,
    ) = indel_type_features(variant)

    (
        n_repeats,
        lc,
        loc_lc,
        gc,
        loc_gc,
        strength,
        loc_strength,
        dissim,
        indel_complexity,
        ref_cnt,
        alt_cnt,
        orig_ref_cnt,
        orig_alt_cnt,
        is_bidirectional,
        is_uniq_mapped,
        uniq_mapping_rate,
        is_near_exon_boundaray,
        equivalent_exists,
        is_multiallelic,
        cplx_variant
    ) = (-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, variant)

    res = make_indel_alignment(variant, bam, downsample_threshold)
    if res:
        try:
    #    if res:
            valn, contig = res[0], res[1]

            (
                n_repeats,
                lc,
                loc_lc,
                gc,
                loc_gc,
                strength,
                loc_strength,
                dissim,
                indel_complexity,
                cplx_variant,
            ) = sequence_features(variant, valn, contig)
            ref_cnt, alt_cnt, orig_ref_cnt, orig_alt_cnt, is_bidirectional = read_support_features(valn, downsample_threshold)

            (
                is_uniq_mapped,
                uniq_mapping_rate,
                is_near_exon_boundaray,
                equivalent_exists,
                is_multiallelic,
            ) = mapping_features(variant, valn, bam, mapq)
        except:
            pass
            #print(variant.pos, variant.chrom, variant.ref, variant.alt)

    return (
        n_repeats,
        lc,
        loc_lc,
        gc,
        loc_gc,
        strength,
        loc_strength,
        dissim,
        indel_complexity,
        indel_size,
        is_ins,
        is_at_ins,
        is_at_del,
        is_gc_ins,
        is_gc_del,
        ref_cnt,
        alt_cnt,
        orig_ref_cnt,
        orig_alt_cnt,
        is_bidirectional,
        is_uniq_mapped,
        uniq_mapping_rate,
        is_near_exon_boundaray,
        equivalent_exists,
        is_multiallelic,
        cplx_variant,
    )


def indel_type_features(variant):
    indel_seq = variant.indel_seq
    indel_size = len(indel_seq)

    is_ins = variant.is_ins
    is_at_ins, is_at_del, is_gc_ins, is_gc_del = 0, 0, 0, 0
    if indel_size == 1:
        if indel_seq in ["A", "T"]:
            if is_ins:
                is_at_ins = 1
            else:
                is_at_del = 1
        else:
            if is_ins:
                is_gc_ins = 1
            else:
                is_gc_del = 1

    return indel_size, is_ins, is_at_ins, is_at_del, is_gc_ins, is_gc_del


def make_indel_alignment(variant, bam, downsample_threshold=1500):
    
    valn = VariantAlignment(variant, bam, downsample_threshold=downsample_threshold)
    
    contig = valn.get_contig()
    if contig:
        return valn, contig
    else:
        return None


def read_support_features(valn, downsample_threshold=1500):
    
    orig_ref_cnt, orig_alt_cnt = valn.count_alleles(by_fragment=True)
    cov = orig_ref_cnt + orig_alt_cnt

    if cov > downsample_threshold:
        ref_cnt, alt_cnt = orig_ref_cnt * downsample_threshold / cov, orig_alt_cnt * downsample_threshold / cov 
    else:
        ref_cnt, alt_cnt = orig_ref_cnt, orig_alt_cnt

    alt_fw_rv = valn.count_alleles(fwrv=True)[1]
    is_bidirectional = all(alt_fw_rv)

    return int(ref_cnt), int(alt_cnt), int(orig_ref_cnt), int(orig_alt_cnt), is_bidirectional


def sequence_features(target_indel, valn, contig):

    contig_seq_tuple = contig.get_contig_seq(split=True)
    lt_seq, rt_seq = contig_seq_tuple[0], contig_seq_tuple[-1]

    contig_ref_seq_tuple = contig.get_reference_seq(split=True)
    lt_ref, rt_ref = contig_ref_seq_tuple[0], contig_ref_seq_tuple[-1]

    # repeat
    if target_indel.is_del:
        non_target_reads = valn.fetch_reads(how="non_target")
        indel_seq = infer_del_seq_from_data(non_target_reads, target_indel)
        indel_type = 0
    else:
        indel_seq = contig_seq_tuple[1]
        indel_type = 1

    n_repeats = repeat(indel_type, lt_seq, indel_seq, rt_seq)

    regional, local = 50, 6
    lt_regional, rt_regional = lt_seq[-regional:], rt_seq[:regional]
    lt_local, rt_local = lt_seq[-local:], rt_seq[:local]

    if "N" in lt_regional:
        lt_regional = lt_ref[-regional:]

    if "N" in rt_regional:
        rt_regional = rt_ref[:regional]

    if "N" in lt_local:
        lt_local = lt_ref[-local:]

    if "N" in rt_local:
        rt_local = rt_ref[:local]
    
    # linguistic complexity
    lt_lc = linguistic_complexity(lt_regional)
    rt_lc = linguistic_complexity(rt_regional)
    lc = (lt_lc + rt_lc) / 2

    # local lc
    loc_lt_lc = linguistic_complexity(lt_local)
    loc_rt_lc = linguistic_complexity(rt_local)
    loc_lc = min(loc_lt_lc, loc_rt_lc)

    # gc and local gc
    if target_indel.is_ins:
        gc_seq = lt_regional + indel_seq + rt_regional
        loc_gc_seq = lt_local + indel_seq + rt_local
    else:
        gc_seq = lt_regional + rt_regional
        loc_gc_seq = lt_local + rt_local

    _gc, loc_gc = gc(gc_seq), gc(loc_gc_seq)

    # strength and local strength
    if target_indel.is_ins:
        strength_seq = lt_regional + rt_regional  # orignal seq
        loc_strength_seq = lt_local + rt_local
    else:
        strength_seq = lt_regional + target_indel.indel_seq + rt_regional
        loc_strength_seq = lt_local + target_indel.indel_seq + rt_local

    strength, loc_strength = dna_strength(strength_seq), dna_strength(loc_strength_seq)

    # dissimilarity
    dissim = dissimilarity(lt_seq, indel_seq, rt_seq)
    
    # indel complexity
    cplx_var = valn.phase(how="complex")
    if not cplx_var.is_non_complex_indel():
        
        lt_len = min(len(lt_seq), len(lt_ref), local)
        rt_len = min(len(rt_seq), len(rt_ref), local)

        lt_complexity = editdistance(lt_ref[-lt_len:], lt_seq[-lt_len:])
        rt_complexity = editdistance(rt_ref[:rt_len], rt_seq[:rt_len])

        indel_complexity = lt_complexity + rt_complexity
    else:
        indel_complexity = 0

    return (
        n_repeats,
        lc,
        loc_lc,
        _gc,
        loc_gc,
        strength,
        loc_strength,
        dissim,
        indel_complexity,
        cplx_var
    )


def infer_del_seq_from_data(non_target_reads, target_deletion):
    if len(non_target_reads) > 20:
        random.seed(123)
        non_target_reads = random.sample(non_target_reads, 20)

    non_ref_del_seq = []

    del_seq = target_deletion.indel_seq
    del_len = len(del_seq)
    for non_target_read in non_target_reads:
        lt_seq, rt_seq = split(non_target_read, target_deletion, is_for_ref=False)
        if len(rt_seq) > del_len:
            inferred = rt_seq[:del_len]
            if inferred != del_seq:
                non_ref_del_seq.append(inferred)

    if non_ref_del_seq:
        return most_common(non_ref_del_seq)
    else:
        return del_seq


def mapping_features(target_indel, valn, bam, mapq):
    target_reads = valn.fetch_reads()
    if len(target_reads) > 20:
        random.seed(123)
        target_reads = random.sample(target_reads, 20)

    non_target_reads = valn.fetch_reads(how="non_target")
    if len(non_target_reads) > 20:
        random.seed(123)
        non_target_reads = random.sample(non_target_reads, 20)

    # mapping
    mapq_from_bam_header = get_star_uniq_mapq(bam)
    uniq_map = mapq_from_bam_header if mapq_from_bam_header else mapq
    is_uniq_mapped = (
        1
        if most_common([read.mapping_quality for read in target_reads]) == uniq_map
        else 0
    )
    
    # uniq_mapping_rate
    covering_reads = valn.fetch_reads(how="covering")
    n_uniq_reads = sum([read.mapping_quality == uniq_map for read in covering_reads])
    uniq_mapping_rate = n_uniq_reads / (len(covering_reads) + 0.0001)


    # near exon
    is_near_exon_boundaray = (
        1
        if most_common([is_near_boundary(read, target_indel) for read in target_reads])
        else 0
    )

    reference = target_indel.reference
    # equivalence
    indels = flatten_list_of_list(
        [get_indels_from_read(read, reference) for read in target_reads]
    )
    equivalent_positions = set([indel.pos for indel in indels if indel == target_indel])
    equivalent_exists = 1 if len(equivalent_positions) > 1 else 0

    # multialleleic
    non_target_indels = flatten_list_of_list(
        [get_indels_from_read(read, reference) for read in non_target_reads]
    )
    non_target_indels_at_target_pos = [
        i for i in non_target_indels if i.normalize().pos == target_indel.pos
    ]
    is_multiallelic = 1 if len(non_target_indels_at_target_pos) > 0 else 0

    return is_uniq_mapped, uniq_mapping_rate, is_near_exon_boundaray, equivalent_exists, is_multiallelic


def get_star_uniq_mapq(bam):
    star_cmd_line = bam.header["PG"][0]["CL"].split("--")
    mapq_option = [option for option in star_cmd_line if "outSAMmapqUnique" in option]

    if mapq_option:
        return int(mapq_option[0].replace("outSAMmapqUnique", ""))


class MockVariant(object):
    def __init__(self, chrom, pos, reference):
        self.chrom = str(chrom)
        self.pos = pos
        self.reference = reference


def get_indels_from_read(aligned_segment, reference, return_index=False, target=None):

    indels = []
    cigarstring = aligned_segment.cigarstring
    if not "I" in cigarstring and not "D" in cigarstring:
        return indels

    aln_start_pos = aligned_segment.reference_start  # 0-based

    is_prefixed = True if reference.references[0].startswith("chr") else False
    chrom = str(aligned_segment.reference_name).replace("chr", "")
    if is_prefixed:
        chrom = "chr" + chrom

    for idx, cigar in enumerate(cigar_ptrn.findall(cigarstring)):
        event, event_len = cigar[-1], int(cigar[:-1])
        if event == "I":
            m = MockVariant(chrom, aln_start_pos, reference)
            lt_seq, rt_seq = split(aligned_segment, m, False)
            lt_ref, rt_ref = split(aligned_segment, m, True)

            i = Variant(
                chrom,
                aln_start_pos,
                lt_ref[-1],
                lt_ref[-1] + rt_seq[:event_len],
                reference,
            )
            indels.append(i)

            if return_index and (i == target):
                return idx

        elif event == "D":
            m = MockVariant(chrom, aln_start_pos, reference)
            lt_ref, rt_ref = split(aligned_segment, m, True)

            d = Variant(
                chrom,
                aln_start_pos,
                lt_ref[-1] + rt_ref[:event_len],
                lt_ref[-1],
                reference,
            )
            indels.append(d)

            aln_start_pos += event_len

            if return_index and (d == target):
                return idx

        elif event in ["S", "H", "P"]:
            pass
        else:
            aln_start_pos += event_len

    return indels


def is_near_boundary(aligned_segment, target):
    cigarstring = aligned_segment.cigarstring
    if not "N" in cigarstring:
        return False

    idx = get_indels_from_read(
        aligned_segment, target.reference, return_index=True, target=target
    )

    cigar_lst = cigar_ptrn.findall(cigarstring)
    if not isinstance(idx, int):
        return False

    threshold = 2 if len(target.indel_seq) <= 2 else 3

    is_near = 0
    if idx >= 2 and "N" in cigar_lst[idx - 2]:
        dist_to_lt_boundary = int(cigar_lst[idx - 1].replace("M", ""))
        if dist_to_lt_boundary <= threshold:
            is_near = 1
    elif (idx + 2) <= len(cigar_lst) - 1 and "N" in cigar_lst[idx + 2]:
        dist_to_rt_boundary = int(cigar_lst[idx + 1].replace("M", ""))
        if dist_to_rt_boundary <= threshold:
            is_near = 1

    return is_near
