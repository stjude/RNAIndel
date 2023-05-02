import os
import pysam
import pandas as pd
from ssw import AlignmentMgr
from indelpost import Variant
from functools import partial
from multiprocessing import Pool

CANONICALS = [str(i) for i in range(1, 23)] + ["X", "Y"]


def realn_softclips(
    bam, fasta, tmp_dir, data_dir, region, num_of_processes, safety_mode
):
    if safety_mode:
        num_of_processes = 1

    if region:
        softclip(region, bam, fasta, data_dir, tmp_dir)
    elif num_of_processes == 1:
        for chromosome in CANONICALS:
            softclip(chromosome, bam, fasta, data_dir, tmp_dir)
        merge_outfiles(tmp_dir)
    else:
        pool = Pool(num_of_processes)

        pool.map(
            partial(softclip, bam=bam, fasta=fasta, data_dir=data_dir, tmp_dir=tmp_dir),
            CANONICALS,
        )


def merge_outfiles(tmp_dir):
    outfilename = os.path.join(tmp_dir, "outfile.sftclp.txt")
    files_by_chrom = [
        os.path.join(tmp_dir, "chr{}.sftclp.txt".format(chrom))
        for chrom in CANONICALS
        if os.path.isfile(os.path.join(tmp_dir, "chr{}.sftclp.txt".format(chrom)))
    ]

    _dfs = [pd.read_csv(_file, sep="\t") for _file in files_by_chrom]
    df = pd.concat(_dfs, ignore_index=True, sort=False)
    df.to_csv(outfilename, sep="\t", index=False)


def softclip(
    region,
    bam,
    fasta,
    data_dir,
    tmp_dir,
    min_clip_len=10,
    base_qual_thresh=12,
    dirty_rate_thresh=0.2,
    mapq_thresh=1,
):
    _parsed_region = region.split(":")
    chrom = _parsed_region[0]
    start, end = 0, 0
    if len(_parsed_region) > 1:
        span = _parsed_region[1].split("-")
        start, end = int(span[0]), int(span[1])

    bam = pysam.AlignmentFile(bam)

    if bam.header.references[0].startswith("chr"):
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
    else:
        if chrom.startswith("chr"):
            chrom = chrom[3:]

    fasta = pysam.FastaFile(fasta)
    if fasta.references[0].startswith("chr"):
        if not chrom.startswith("chr"):
            fasta_chrom = "chr" + chrom
        else:
            fasta_chrom = chrom
    else:
        if chrom.startswith("chr"):
            fasta_chrom = chrom[3:]
        else:
            fasta_chrom = chrom

    if not chrom.startswith("chr"):
        db_chrom = "chr" + chrom
    else:
        db_chrom = chrom

    refseq_db = pysam.TabixFile("{}/refgene/refCodingExon.bed.gz".format(data_dir))

    if end:
        reads = bam.fetch(chrom, start, end)
    else:
        reads = bam.fetch(chrom)

    clipped_reads = []

    is_first_read = True
    window_width = 20
    window_start_pos = start
    window_end_pos = window_start_pos + window_width
    window_volume = 0

    if ":" in region:
        outfilename = os.path.join(tmp_dir, "outfile.sftclp.txt")
    else:
        no_chr = chrom.replace("chr", "")
        outfilename = os.path.join(tmp_dir, "chr{}.sftclp.txt".format(no_chr))

    f = open(outfilename, "w")
    f.write("Chr\tPos\tChr_Allele\tAlternative_Allele\tType\n")

    indels = []
    for read in reads:
        if is_analyzable(
            read, window_start_pos, base_qual_thresh, dirty_rate_thresh, mapq_thresh
        ):
            # use first analyzable read to start window
            if is_first_read:
                window_start_pos = read.reference_start
                is_first_read = False

            window_volume += 1

            if is_hq_clip(
                read,
                base_qual_thresh,
                dirty_rate_thresh,
                min_clip_len,
                homopoly_junction_len=5,
            ):
                clipped_reads.append(read)

        if is_beyond_window(read, window_end_pos):
            try:
                tmp = process_clipped_reads(
                    fasta_chrom,
                    db_chrom,
                    clipped_reads,
                    window_volume,
                    window_start_pos,
                    window_end_pos,
                    refseq_db,
                    fasta,
                    min_clip_len,
                    density_thresh=0.05,
                    margin=50,
                )
            except:
                tmp = []

            if tmp:
                indels += tmp

            clipped_reads.clear()
            window_start_pos = read.reference_start
            window_end_pos = window_start_pos + window_width
            window_volume = 0

    indels = filter_indels(indels, min_occurrence=2)

    for v in indels:
        _chrom = v.chrom
        if not _chrom.startswith("chr"):
            _chrom = "chr" + _chrom

        if v.is_del:
            _ref = v.ref[1:]
            _alt = "-"
            _type = "deletion"
        else:
            _ref = "-"
            _alt = v.alt[1:]
            _type = "insertion"

        f.write("{}\t{}\t{}\t{}\t{}\n".format(_chrom, v.pos + 1, _ref, _alt, _type))

    f.close()


def filter_indels(indel_lst, min_occurrence):
    indel_dict = {}
    for _idl in indel_lst:
        if indel_dict.get(_idl[0], None):
            indel_stat = indel_dict[_idl[0]]
            indel_stat[0].append(_idl[1])
            indel_stat[1].append(_idl[2])
        else:
            indel_dict[_idl[0]] = [[_idl[1]], [_idl[2]]]

    lst = []
    for k, v in indel_dict.items():
        if len(set(v[0])) >= min_occurrence:
            lst.append(k)
    return lst


def process_clipped_reads(
    fasta_chrom,
    db_chrom,
    clipped_reads,
    window_volume,
    window_start_pos,
    window_end_pos,
    refseq_db,
    fasta,
    min_clip_len,
    density_thresh=0.1,
    margin=50,
):
    if len(clipped_reads) / (window_volume + 1) < density_thresh:
        return None

    try:
        is_coding = len(
            list(refseq_db.fetch(db_chrom, window_start_pos - 1, window_end_pos + 1))
        )
    except:
        is_coding = 0
    
    if not is_coding:
        return None

    lt_clips, rt_clips = sort_by_clip_ends(clipped_reads)

    indels_lt = process_clip_dict(
        fasta_chrom, lt_clips, fasta, margin, min_clip_len, lt_sided=True
    )

    indels_rt = process_clip_dict(
        fasta_chrom, rt_clips, fasta, margin, min_clip_len, lt_sided=False
    )

    return indels_lt + indels_rt


def process_clip_dict(
    fasta_chrom, clips, fasta, margin, min_clip_len, lt_sided, mut_num_thresh=5
):
    indels = []
    for clip_pos, clip_data in clips.items():
        for _read in clip_data:
            read = _read[1]
            ref_seq, pos_dict = preprocess_realn(
                fasta_chrom, read, fasta, lt_sided, margin
            )
            aln = realn(
                read,
                ref_seq,
                match_score=3,
                mismatch_penalty=2,
                gap_open=3,
                gap_extension=1,
            )

            if is_hq_realn(read, aln, ref_seq, min_clip_len):
                variants = postprocess_realn(
                    fasta_chrom, fasta, ref_seq, pos_dict, read, aln, n_gap_thresh=2
                )

                if len(variants) <= mut_num_thresh:
                    tmp = [
                        (v, read.query_name, read.is_reverse)
                        for v in variants
                        if v.is_indel
                    ]
                    indels += tmp

    return indels


def preprocess_realn(fasta_chrom, read, fasta, lt_sided, margin=50):
    ref_seq, pos_dict = get_spliced_reference_seq(
        fasta_chrom, read.reference_start, read.cigartuples, fasta, lt_sided, margin
    )

    return ref_seq, pos_dict


def realn(
    read, ref_seq, match_score=3, mismatch_penalty=2, gap_open=3, gap_extension=0
):
    aln_mgr = AlignmentMgr(match_score=match_score, mismatch_penalty=mismatch_penalty)
    aln_mgr.set_read(read.query_sequence)
    aln_mgr.set_reference(ref_seq)

    return aln_mgr.align(gap_open=gap_open, gap_extension=gap_extension)


def postprocess_realn(chrom, fasta, ref_seq, pos_dict, read, aln, n_gap_thresh):
    candidate_indels = []

    n_gaps = aln.CIGAR.count("D") + aln.CIGAR.count("I")

    if n_gaps < 1 or n_gaps > n_gap_thresh:
        return candidate_indels

    n_mapped_bases_orig = sum([c[1] for c in read.cigartuples if c[0] in (0, 1)])
    n_mapped_bases_realn = sum(
        [c[0] for c in aln.cigar_pair_list if c[1] in ("M", "I")]
    )

    if n_mapped_bases_realn <= n_mapped_bases_orig:
        return candidate_indels

    return find_all_variants(chrom, aln, pos_dict, ref_seq, read, fasta)


def get_mapped_ends_n(aln):
    cigar_list = aln.cigar_pair_list
    lt_n, rt_n = 0, 0
    if cigar_list[0][1] == "M":
        lt_n = cigar_list[0][0]

    if cigar_list[-1][1] == "M":
        rt_n = cigar_list[-1][0]

    return lt_n, rt_n


def is_hq_realn(read, aln, ref_seq, min_clip_len, end_match_thresh=6):
    read_seq = read_seq = read.query_sequence

    read_start_idx = aln.read_start
    read_end_idx = aln.read_end
    ref_start_idx = aln.reference_start
    ref_end_idx = aln.reference_end

    n_lt_clip = read_start_idx
    n_rt_clip = len(read_seq) - read_end_idx

    if (n_lt_clip + n_rt_clip) > 5:
        return False

    lt_n_mapped, rt_n_mapped = get_mapped_ends_n(aln)

    if min(lt_n_mapped, rt_n_mapped) < min_clip_len:
        return False

    lt_mapped_read = read_seq[read_start_idx : read_start_idx + lt_n_mapped]
    lt_mapped_ref = ref_seq[ref_start_idx : ref_start_idx + lt_n_mapped]

    if lt_mapped_read[:end_match_thresh] != lt_mapped_ref[:end_match_thresh]:
        return False

    rt_mapped_read = read_seq[read_end_idx + 1 - rt_n_mapped : read_end_idx + 1]
    rt_mapped_ref = ref_seq[ref_end_idx + 1 - rt_n_mapped : ref_end_idx + 1]

    if rt_mapped_read[-end_match_thresh:] != rt_mapped_ref[-end_match_thresh:]:
        return False

    return True


def find_all_variants(chrom, aln, pos_dict, ref_seq, read, fasta):
    variants = []
    genome_pos = pos_dict[0]
    ref_idx = aln.reference_start
    read_idx = aln.read_start

    read_seq = read.query_sequence
    cigar_list = aln.cigar_pair_list

    for c in cigar_list:
        op, op_len = c[1], c[0]

        if op == "M":
            i = 0
            while i < op_len:
                _ref = ref_seq[ref_idx + i : ref_idx + i + 1]
                _alt = read_seq[read_idx + i : read_idx + i + 1]

                if _ref != _alt:
                    snv = Variant(chrom, genome_pos + i + 1, _ref, _alt, fasta)

                    variants.append(snv)

                i += 1

            ref_idx += op_len
            read_idx += op_len
            genome_pos = pos_dict.get(ref_idx, -1)
        elif op == "I":
            _ref = ref_seq[ref_idx - 1 : ref_idx]
            _alt = read_seq[read_idx : read_idx + op_len]
            ins = Variant(chrom, genome_pos, _ref, _ref + _alt, fasta)
            variants.append(ins)
            read_idx += op_len
        elif op == "D":
            _alt = ref_seq[ref_idx - 1 : ref_idx]
            _ref = ref_seq[ref_idx : ref_idx + op_len]
            del_ = Variant(chrom, genome_pos, _alt + _ref, _alt, fasta)

            ref_idx += op_len

            _prev_pos = genome_pos
            genome_pos = pos_dict.get(ref_idx, -1)

            # skip involved
            if genome_pos - _prev_pos != op_len:
                pass
            else:
                variants.append(del_)

    return variants


def is_hq_clip(
    read,
    base_qual_thresh,
    dirty_rate_thresh,
    min_clip_len,
    homopoly_junction_len=5,
):
    if "S" not in read.cigarstring:
        return False

    read_seq = read.query_sequence
    quals = read.query_qualities
    cigar_tuples = read.cigartuples

    first_cigar = cigar_tuples[0]
    last_cigar = cigar_tuples[-1]
    len_lt, len_rt = 0, 0
    n_lt_dirty_bases, n_rt_dirty_bases = 0, 0

    if first_cigar[0] == 4:
        len_lt = first_cigar[1]
        lt_clip = read_seq[:len_lt]

        if "N" in lt_clip:
            return False

        # clipped-homopolymer
        if lt_clip == len_lt * lt_clip[0]:
            return False

        # homopolymer junction
        lt_junc = read_seq[len_lt : len_lt + homopoly_junction_len]
        if lt_clip == homopoly_junction_len * lt_junc[0]:
            return False

        n_lt_dirty_bases = sum([q <= base_qual_thresh for q in quals[:len_lt]])

    if last_cigar[0] == 4:
        len_rt = last_cigar[1]
        rt_clip = read_seq[-len_rt:]

        if "N" in rt_clip:
            return False

        if rt_clip == len_rt * rt_clip[-1]:
            return False

        rt_junc = read_seq[-(len_rt + homopoly_junction_len) : -len_rt]
        if rt_junc == homopoly_junction_len * rt_junc[-1]:
            return False

        n_rt_dirty_bases = sum([q <= base_qual_thresh for q in quals[-len_rt:]])

    if max(len_lt, len_rt) < min_clip_len:
        return False

    dirty_clip_rates = (n_lt_dirty_bases + n_rt_dirty_bases) / (len_lt + len_rt)

    return dirty_clip_rates <= dirty_rate_thresh / 2


def update_pos_dict(pos_dict, start_pos, span):
    i = 0
    offset = len(pos_dict)
    for j in range(offset, offset + span):
        pos_dict[j] = start_pos + i
        i += 1


def get_spliced_reference_seq(chrom, aln_start, cigar_tuples, fasta, lt_sided, margin):
    __pos = aln_start  # 0-based

    consuming_operations = (0, 2, 8)

    pos_dict = {}

    if lt_sided:
        ref_seq = fasta.fetch(chrom, __pos - margin, __pos)
        update_pos_dict(pos_dict, __pos - margin, margin)
    else:
        ref_seq = ""

    for c in cigar_tuples:
        op, op_len = c[0], c[1]
        if op in consuming_operations:
            ref_seq += fasta.fetch(chrom, __pos, __pos + op_len)

            update_pos_dict(pos_dict, __pos, op_len)

            __pos += op_len
        elif op == 3:
            __pos += op_len
        else:
            pass

    if not lt_sided:
        ref_seq += fasta.fetch(chrom, __pos, __pos + margin)
        update_pos_dict(pos_dict, __pos, margin)

    return ref_seq, pos_dict


def dirty_base_rate(read, base_qual_thresh):
    try:
        quals = read.query_qualities
        return sum([q <= base_qual_thresh for q in quals]) / len(quals)
    except:
        return 1.0


def is_analyzable(
    read, window_start_pos, base_qual_thresh, dirty_rate_thresh, mapq_thresh
):
    cigar_str = read.cigarstring

    if (
        cigar_str
        and not read.is_duplicate
        and not read.is_secondary
        and not read.is_unmapped
    ):
        aln_start = read.reference_start

        if dirty_base_rate(read, base_qual_thresh) >= dirty_rate_thresh:
            return False

        if read.mapping_quality < mapq_thresh:
            return False

        return True

    return False


def is_beyond_window(read, window_end_pos):
    cigar_str = read.cigarstring

    if cigar_str:
        cigar_tuples = read.cigartuples
        curr_pos = read.reference_start + 1

        return window_end_pos < curr_pos

    return False


def sort_by_clip_ends(clipped_reads):
    lt_clipps, rt_clipps = {}, {}

    for read in clipped_reads:
        cigar_tupes = read.cigartuples
        read_seq = read.query_sequence
        first_cigar, last_cigar = cigar_tupes[0], cigar_tupes[-1]
        if first_cigar[0] == 4:
            clip_end = read.reference_start
            if lt_clipps.get(clip_end, None):
                lt_clipps[clip_end].append((read_seq[: first_cigar[1]], read))
            else:
                lt_clipps[clip_end] = [(read_seq[: first_cigar[1]], read)]

        if last_cigar[0] == 4:
            clip_start = read.reference_end
            if rt_clipps.get(clip_start, None):
                rt_clipps[clip_start].append((read_seq[-last_cigar[1] :], read))
            else:
                rt_clipps[clip_start] = [(read_seq[-last_cigar[1] :], read)]

    return lt_clipps, rt_clipps
