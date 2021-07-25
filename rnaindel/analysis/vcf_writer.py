import os
import pysam
import datetime
from .classifier import make_feature_dict


def write_vcf(df, version, arguments):

    header = (
        header_1(version, arguments.reference)
        + filter_fields()
        + info_fields()
        + format_fields()
        + get_cmd_line(arguments)
        + get_features_used(arguments.data_dir)
        + bottom_header(arguments.bam)
    )

    f = open(arguments.output_vcf, "w")

    f.write("\n".join(header))

    for idx, row in df.iterrows():
        vcf_record = parse_row_to_vcf_record(row)
        f.write(vcf_record + "\n")

    f.close()


def header_1(version, fasta):
    h = [
        "##fileformat=VCFv4.2",
        "##filedate={}".format(get_today()),
        "##source=RNAIndel",
        "##source_version={}".format(version),
        "##reference={}".format(fasta),
    ]

    return h


def get_today():
    dt = datetime.datetime.now()
    return str(dt.year) + str(dt.month) + str(dt.day)


def get_cmd_line(args):
    cmd = "rnaindel PredictSomaticIndels -i {} -r {} -o {} -d {} -m {} -p {} -q {} --region {}".format(
        os.path.abspath(args.bam),
        os.path.abspath(args.reference),
        os.path.abspath(args.output_vcf),
        os.path.abspath(args.data_dir),
        args.heap_memory,
        args.process_num,
        args.uniq_mapq,
        args.region,
    )

    h = ["##cmdline={}".format(cmd)]
    return h


def get_features_used(data_dir):
    feature_dict = make_feature_dict("{}/models".format(data_dir))
    used_for_single = ",".join(feature_dict["single_nucleotide_indels"])
    used_for_multi = ",".join(feature_dict["multi_nucleotide_indels"])

    h = [
        "##features_used_for_single_nucleotide_indel_prediction={}".format(
            used_for_single
        ),
        "##features_used_for_multi_nucleotide_indel_prediction={}".format(
            used_for_multi
        ),
    ]

    return h


def filter_fields():
    h = [
        '##FILTER=<ID=ProbableArtifact,Description="Matches to known non-somatic indel with the predicted artifact probability higher than the germline probability.">',
        '##FILTER=<ID=ProbableGermline,Description="Matches to known common indel or non-somatic indel with the predicted germline probability higher than the artifact probability">',
        '##FILTER=<ID=LowMappabilityRegion,Description="More than half of reads covering the indel locus having MAPQ < unique mapping quality score.">',
    ]
    return h


def info_fields():
    h = [
        '##INFO=<ID=predicted_class,Number=1,Type=String,Description="Either of somatic, germline, or, artifact.">',
        '##INFO=<ID=probabilities,Number=3,Type=Float,Description="Probabality of being somatic, germline, artifact in this order.">',
        '##INFO=<ID=annotation,Number=.,Type=String,Description="Indel annotation formatted as GeneSymbol|RefSeqAccession|CodonPos|VariantEffect. Delimited by comma for multiple isoforms.">',
        '##INFO=<ID=cosmic_cnt,Number=1,Type=Integer,Description="COSMIC count found in COSMICv89.">',
        '##INFO=<ID=maximum_population_frequency,Number=1,Type=Float,Description="Maximum of the allele frequencies reported in dbSNPv151, ClinVar20180603, and gnomAD_r2.1. non-cancer population.">',
        '##INFO=<ID=COMMON,Number=0,Type=Flag,Description="Annotated COMMON in dbsnpv151 or maximum_population_frequency > 0.01.">',
        '##INFO=<ID=repeat,Number=1,Type=Integer,Description="occurrence of the indel-sequence repeat units inflanking region.">',
        '##INFO=<ID=lc,Number=1,Type=Float,Description="Linguistic complexity in flanking 50-bp region.">',
        '##INFO=<ID=local_lc,Number=1,Type=Float,Description="Linguistic complexity in flanking 6-bp region.">',
        '##INFO=<ID=gc,Number=1,Type=Float,Description="GC content in flanking 50-bp region.">',
        '##INFO=<ID=local_gc,Number=1,Type=Float,Description="GC content in flanking 6-bp region.">',
        '##INFO=<ID=strength,Number=1,Type=Float,Description="DNA pair-bond strength of 2-mer in flanking 50-bp region.">',
        '##INFO=<ID=local_strength,Number=1,Type=Float,Description="DNA pair-bond strength of 2-mer in flanking 6-bp region.">',
        '##INFO=<ID=dissimilarity,Number=1,Type=Float,Description="Edit distance between indel and flanking sequences.">',
        '##INFO=<ID=indel_complexity,Number=1,Type=Integer,Description="Mismatches in cis configuration to indel.">',
        '##INFO=<ID=indel_size,Number=1,Type=Integer,Description="Length of the indel sequence.">',
        '##INFO=<ID=INS,Number=0,Type=Flag,Description="This indel is an insertion.">',
        '##INFO=<ID=AT_INS,Number=0,Type=Flag,Description="This indel is an  A/T single-insertion.">',
        '##INFO=<ID=AT_DEL,Number=0,Type=Flag,Description="This indel is an  A/T single-deletion.">',
        '##INFO=<ID=GC_INS,Number=0,Type=Flag,Description="This indel is an  G/C single-insertion.">',
        '##INFO=<ID=GC_DEL,Number=0,Type=Flag,Description="This indel is an  G/C single-deletion.">',
        '##INFO=<ID=ref_count,Number=1,Type=Integer,Description="Read-fragments representing the reference allele.">',
        '##INFO=<ID=alt_count,Number=1,Type=Integer,Description="Read-fragments representing the indel allele.">',
        '##INFO=<ID=BIDIRECTIONAL,Number=0,Type=Flag,Description="This indel is supported by both forward and reverse reads.">',
        '##INFO=<ID=UNIQ_MAPPED,Number=0,Type=Flag,Description="This indel is supported by at least one uniquely mapped read.">',
        '##INFO=<ID=NEAR_EXON_BOUNDARY,Number=0,Type=Flag,Description="This indel is within exon but near the exon boundary.">',
        '##INFO=<ID=EQUIVALENCE_EXISTS,Number=0,Type=Flag,Description="Alternative indel alignments are observed for this indel.">',
        '##INFO=<ID=MULTIALLELIC,Number=0,Type=Flag,Description="Multiple different indels are observed at this locus.">',
        '##INFO=<ID=INFRAME,Number=0,Type=Flag,Description="This indel is inframe.">',
        '##INFO=<ID=SPLICE,Number=0,Type=Flag,Description="This indel is intronic and within 10-bp to exon (splice region).">',
        '##INFO=<ID=TRUNCATING,Number=0,Type=Flag,Description="This indel is frameshift or stopgain or destroys splice motif.">',
        '##INFO=<ID=CONSERVED_DOMAIN,Number=0,Type=Flag,Description="This indel is located in conserved domain.">',
        '##INFO=<ID=indel_location,Number=1,Type=Float,Description="Relative indel location in coding region.">',
        '##INFO=<ID=NMD_INSENSITIVE,Number=0,Type=Flag,Description="This indel is insensitive to nonsense-mediated decay.">',
        '##INFO=<ID=indel_per_gene,Number=1,Type=Float,Description="Indels detected in the gene in this sample. Normalized by coding sequence length.">',
        '##INFO=<ID=coding_sequence_length,Number=1,Type=Float,Description="Coding sequence length. Median value if multiple isoforms exist.">',
        '##INFO=<ID=GERMLINE_DB,Number=0,Type=Flag,Description="This indel is in the (default or user-provided) germline database.">',
        '##INFO=<ID=RECLASSIFIED,Number=0,Type=Flag,Description="This indel is reclassified from the initial prediction.">',
        '##INFO=<ID=OUTLYING,Number=0,Type=Flag,Description="This homopolymer indel is outlying from artifact homoplymer indels of the same kind.">',
        '##INFO=<ID=CALLER,Number=1,Type=String,Description="This indel is called by built_in, external, or, both.">',
    ]

    return h


def format_fields():
    h = ['##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allele depth.">']
    return h


def bottom_header(bam):
    h = [
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n".format(
            get_sample_name(bam)
        )
    ]

    return h


def get_sample_name(bam):
    try:
        bamheader = pysam.AlignmentFile(bam).header
        samplename = bamheader["RG"][0]["SM"]
    except:
        samplename = "SampleNameNotAvailable"

    return samplename


def parse_row_to_vcf_record(row):
    chrom = str(row["chrom"])
    pos = str(row["cpos"])
    snpid = row["dbsnp"]
    ref = row["cref"]
    alt = row["calt"]
    qual = "."
    fltr = row["filter"]
    info_str = generate_info_str(row)
    fmt = "AD"
    sample = "{},{}".format(row["ref_count"], row["alt_count"])

    return "\t".join([chrom, pos, snpid, ref, alt, qual, fltr, info_str, fmt, sample])


def generate_info_str(row):
    info_str = "predicted_class={};probabilities={},{},{};annotation={};cosmic_cnt={};".format(
        row["predicted_class"],
        row["prob_s"],
        row["prob_g"],
        row["prob_a"],
        row["annotation"],
        row["cosmic_cnt"],
    )

    if row["pop_freq"] > 0:
        info_str += "maximum_population_frequency={};".format(row["pop_freq"])

    info_str = extend_flag_str(info_str, "COMMON", row["is_common"])

    info_str += "repeat={};lc={};local_lc={};gc={};local_gc={};strength={};local_strength={};dissimilarity={};indel_complexity={};indel_size={};".format(
        int(row["repeat"]),
        row["lc"],
        row["local_lc"],
        row["gc"],
        row["local_gc"],
        row["strength"],
        row["local_strength"],
        row["dissimilarity"],
        row["indel_complexity"],
        row["indel_size"],
    )

    info_str = extend_flag_str(info_str, "INS", row["is_ins"])
    info_str = extend_flag_str(info_str, "AT_INS", row["is_at_ins"])
    info_str = extend_flag_str(info_str, "AT_DEL", row["is_at_del"])
    info_str = extend_flag_str(info_str, "GC_INS", row["is_gc_ins"])
    info_str = extend_flag_str(info_str, "GC_DEL", row["is_gc_del"])

    info_str += "ref_count={};alt_count={};".format(
        int(row["ref_count"]), int(row["alt_count"])
    )

    info_str = extend_flag_str(info_str, "BIDIRECTIONAL", row["is_bidirectional"])
    info_str = extend_flag_str(info_str, "UNIQ_MAPPED", row["is_uniq_mapped"])
    info_str = extend_flag_str(info_str, "NEAR_EXON_BOUNDARY", row["is_near_boundary"])
    info_str = extend_flag_str(
        info_str, "EQUIVALENCE_EXISTS", row["equivalence_exists"]
    )
    info_str = extend_flag_str(info_str, "MULTIALLELIC", row["is_multiallelic"])
    info_str = extend_flag_str(info_str, "INFRAME", row["is_inframe"])
    info_str = extend_flag_str(info_str, "SPLICE", row["is_splice"])
    info_str = extend_flag_str(info_str, "TRUNCATING", row["is_truncating"])
    info_str = extend_flag_str(info_str, "CONSERVED_DOMAIN", row["is_in_cdd"])

    info_str += "indel_location={};".format(row["indel_location"])

    info_str = extend_flag_str(info_str, "NMD_INSENSITIVE", row["is_nmd_insensitive"])

    info_str += "indel_per_gene={};coding_sequence_length={};".format(
        row["ipg"], row["cds_length"]
    )

    info_str = extend_flag_str(info_str, "GERMLINE_DB", row["is_on_db"])

    if row["reclassified"] != "-":
        info_str += "RECLASSIFIED;"

        if "outlier" in row["reclassified"]:
            info_str += "OUTLYING;"

    info_str += "CALLER={};".format(row["origin"])

    return info_str.rstrip(";")


def extend_flag_str(info_str, flag, flag_boolean):
    if flag_boolean:
        info_str += "{};".format(flag)

    return info_str
