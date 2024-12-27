import os
import csv
import pysam
import pandas as pd
from functools import partial
from multiprocessing import Pool

from indelpost import Variant, VariantAlignment

from .callset_formatter import format_callset
from .coding_indel import annotate_coding_info
from .transcript_feature_calculator import transcript_features
from .alignment_feature_calculator import alignment_features
from .database_feature_calculator import database_features


CANONICALS = [str(i) for i in range(1, 23)] + ["X", "Y"]


def preprocess(
    tmp_dir,
    fasta_file,
    bam_file,
    data_dir,
    mapq,
    num_of_processes,
    region,
    external_vcf,
    pass_only,
    safety_mode,
):
    if num_of_processes == 1:

        callset = format_callset(tmp_dir, external_vcf, pass_only, region)
        df = calculate_features(
            callset, fasta_file, bam_file, data_dir, mapq, external_vcf
        )
    else:
        callsets_by_chrom = format_callset(tmp_dir, external_vcf, pass_only, region)

        if safety_mode:
            dfs = [
                calculate_features(
                    call_set_by_chrom,
                    fasta_file=fasta_file,
                    bam_file=bam_file,
                    data_dir=data_dir,
                    mapq=mapq,
                    external_vcf=external_vcf,
                )
                for call_set_by_chrom in callsets_by_chrom
            ]

            df = pd.concat(dfs)
        else:
            pool = Pool(num_of_processes)

            # workaround for higher coverage datasets (eg., cell-free RNA-Seq)
            _data_files = pool.map(
                partial(
                    calculate_features,
                    fasta_file=fasta_file,
                    bam_file=bam_file,
                    data_dir=data_dir,
                    mapq=mapq,
                    external_vcf=external_vcf,
                    out_as_df=False,
                ),
                callsets_by_chrom,
            )

            pool.close()
            pool.join()

            dfs = []
            for _ in _data_files:
                if ".done.txt" in _:
                    dfs.append(
                        pd.read_csv(
                            _,
                            sep="\t",
                            dtype={
                                "is_in_cdd": bool,
                                "is_ins": bool,
                                "is_bidirectional": bool,
                            },
                        )
                    )

            df = pd.concat(dfs)

            fasta = pysam.FastaFile(fasta_file)
            df["indel"] = df.apply(_remake, fasta=fasta, axis=1)

    return df


def _remake(row, fasta):
    return Variant(row["chrom"], row["pos"], row["ref"], row["alt"], fasta)


def calculate_features(
    callset, fasta_file, bam_file, data_dir, mapq, external_vcf, out_as_df=True
):
    path_to_coding_gene_db = "{}/refgene/refCodingExon.bed.gz".format(data_dir)
    path_to_proteindb = "{}/protein/proteinConservedDomains.txt".format(data_dir)
    path_to_dbsnp = "{}/dbsnp/dbsnp.indel.vcf.gz".format(data_dir)
    path_to_clinvar = "{}/clinvar/clinvar.indel.vcf.gz".format(data_dir)
    path_to_cosmic = "{}/cosmic/CosmicCodingMuts.indel.vcf.gz".format(data_dir)

    df = filter_non_coding_indels(
        callset, fasta_file, path_to_coding_gene_db, external_vcf
    )

    _df_filename = callset.replace(".txt", ".failed.txt")
    if len(df) > 0:
        df = transcript_features(df, path_to_proteindb)
        if len(df) > 0:
            df = alignment_features(df, bam_file, mapq)

        if len(df) > 0:
            df = database_features(df, path_to_dbsnp, path_to_clinvar, path_to_cosmic)
            if out_as_df:
                return df
            else:
                _df_filename = _df_filename.replace(".failed.txt", ".done.txt")
                df.to_csv(_df_filename, sep="\t", index=False)
            return _df_filename

    if out_as_df:
        return make_empty_df()
    else:
        return _df_filename


def filter_non_coding_indels(callset, fasta_file, path_to_coding_gene_db, external_vcf):

    reference = pysam.FastaFile(fasta_file)
    coding_gene_db = pysam.TabixFile(path_to_coding_gene_db)

    coding_indels = []
    is_prefixed = reference.references[0].startswith("chr")
    with open(callset) as f:
        records = csv.DictReader(f, delimiter="\t")
        for record in records:
            try:
                indel, origin = bambino2variant(record, reference, is_prefixed)
                update_coding_indels(coding_indels, indel, origin, coding_gene_db)
            except:
                pass

    if coding_indels:
        df = pd.DataFrame(coding_indels)

        if external_vcf:
            dfg = df.groupby(["chrom", "pos", "ref", "alt"])
            df = dfg.apply(summarize_caller_origin)

        df = df.drop_duplicates(subset=["chrom", "pos", "ref", "alt", "origin"])
        return df
    else:
        header = ["empty"]
        return pd.DataFrame(columns=header)


def update_coding_indels(coding_indels, indel, origin, coding_gene_db):

    coding_annotations = annotate_coding_info(indel, coding_gene_db)
    if coding_annotations:
        d = {
            "indel": indel,
            "chrom": indel.chrom,
            "pos": indel.pos,
            "ref": indel.ref,
            "alt": indel.alt,
            "coding_indel_isoforms": coding_annotations,
            "origin": origin,
        }
        coding_indels.append(d)


def summarize_caller_origin(df_groupedby_indel):
    origins = set(df_groupedby_indel["origin"].to_list())

    if len(origins) > 1:
        df_groupedby_indel["origin"] = "both"

    return df_groupedby_indel


def bambino2variant(record, reference, is_prefixed):
    chrom = record["Chr"].replace("chr", "")

    if not chrom in CANONICALS:
        return None

    chrom = "chr" + chrom if is_prefixed else chrom

    pos = int(record["Pos"])
    ref = record["Chr_Allele"]
    alt = record["Alternative_Allele"]
    var_type = record["Type"]

    origin = "external"
    if var_type in ["deletion", "insertion"]:
        origin = "built_in"
        pos -= 1
        padding_base = reference.fetch(chrom, pos - 1, pos)
        if var_type == "deletion":
            alt = padding_base
            ref = alt + ref
        else:
            ref = padding_base
            alt = ref + alt

    return Variant(chrom, pos, ref, alt, reference).normalize(), origin


def make_empty_df():
    header = [
        "indel",
        "origin",
        "chrom",
        "pos",
        "ref",
        "alt",
        "annotation",
        "cds_length",
        "indel_location",
        "is_inframe",
        "is_splice",
        "is_truncating",
        "is_nmd_insensitive",
        "is_in_cdd",
        "gene_symbol",
        "ipg",
        "repeat",
        "lc",
        "local_lc",
        "gc",
        "local_gc",
        "strength",
        "local_strength",
        "dissimilarity",
        "indel_complexity",
        "indel_size",
        "is_ins",
        "is_at_ins",
        "is_at_del",
        "is_gc_ins",
        "is_gc_del",
        "ref_count",
        "alt_count",
        "orig_ref_cnt",
        "orig_alt_cnt",
        "is_bidirectional",
        "is_uniq_mapped",
        "uniq_mapping_rate",
        "is_near_boundary",
        "equivalence_exists",
        "is_multiallelic",
        "cpos",
        "cref",
        "calt",
        "dbsnp",
        "pop_freq",
        "is_common",
        "is_on_db",
        "is_pathogenic",
        "cosmic_cnt",
    ]

    return pd.DataFrame(columns=header)
