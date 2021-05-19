import os
import csv
import pysam
import pandas as pd
from functools import partial
from multiprocessing import Pool

from indelpost import Variant, VariantAlignment

from .coding_indel import annotate_coding_info
from .transcript_feature_calculator import transcript_features
from .alignment_feature_calculator import alignment_features
from .database_feature_calculator import database_features

CANONICALS = [str(i) for i in range(1, 23)] + ["X", "Y"]


def preprocess(
    tmp_dir, fasta_file, bam_file, data_dir, mapq, num_of_processes, from_default_caller=True
):
    if num_of_processes == 1:
        callset = os.path.join(tmp_dir, "outfile.txt")
        df = calculate_features(callset, fasta_file, bam_file, data_dir, mapq, from_default_caller)

    else:
        callsets_by_chrom = [
            os.path.join(tmp_dir, "chr{}.txt".format(chrom)) for chrom in CANONICALS
        ]
        pool = Pool(num_of_processes)

        dfs = pool.map(
            partial(
                calculate_features,
                fasta_file=fasta_file,
                bam_file=bam_file,
                data_dir=data_dir,
                mapq=mapq,
                from_default_caller=from_default_caller,
            ),
            callsets_by_chrom,
        )
        
        df = pd.concat(dfs)

    return df


def calculate_features(callset, fasta_file, bam_file, data_dir, mapq, from_default_caller):

    path_to_coding_gene_db = "{}/refgene/refCodingExon.bed.gz".format(data_dir)
    path_to_proteindb = "{}/protein/proteinConservedDomains.txt".format(data_dir)
    path_to_dbsnp = "{}/dbsnp/dbsnp.indel.vcf.gz".format(data_dir)
    path_to_clinvar = "{}/clinvar/clinvar.indel.vcf.gz".format(data_dir)
    path_to_cosmic = "{}/cosmic/CosmicCodingMuts.indel.vcf.gz".format(data_dir)
    
    df = filter_non_coding_indels(
        callset, fasta_file, path_to_coding_gene_db, from_default_caller
    )
    
    if len(df) > 0:
        df = transcript_features(df, path_to_proteindb)
        df = alignment_features(df, bam_file, mapq)
        df = database_features(df, path_to_dbsnp, path_to_clinvar, path_to_cosmic)

        return df


def filter_non_coding_indels(
    callset, fasta_file, path_to_coding_gene_db, from_default_caller
):
    reference = pysam.FastaFile(fasta_file)
    coding_gene_db = pysam.TabixFile(path_to_coding_gene_db)

    coding_indels = []
    if from_default_caller:

        # check if chromosome name is "chr"-prefixed
        is_prefixed = reference.references[0].startswith("chr")

        f = open(callset)
        records = csv.DictReader(f, delimiter="\t")
        for record in records:
            indel = bambino2variant(record, reference, is_prefixed)
            update_coding_indels(coding_indels, indel, coding_gene_db)
        f.close()

    return pd.DataFrame(coding_indels)


def update_coding_indels(coding_indels, indel, coding_gene_db):
    if indel:
        coding_annotations = annotate_coding_info(indel, coding_gene_db)
        if coding_annotations:
            d = {"indel": indel, "coding_indel_isoforms": coding_annotations}
            coding_indels.append(d)


def bambino2variant(record, reference, is_prefixed):
    chrom = record["Chr"].replace("chr", "")

    if not chrom.replace("chr", "") in CANONICALS:
        return None

    chrom = "chr" + chrom if is_prefixed else chrom

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

    return Variant(chrom, pos, ref, alt, reference).normalize()


def instantiate_vcf_record(chrom, pos, ref, alt, reference, allele_len_thresh=100):
    if len(ref) == len(alt):
        return None

    if not chrom.replace("chr", "") in CANONICALS:
        return None

    if len(ref) > allele_len_thresh or len(alt) > allele_len_thresh:
        return None

    return Variant(chrom, pos, ref, alt, reference).normalize()


def is_canonical_indel(record):
    is_indel = len(record["REF"]) != len(record["ALT"])

    chrom_name = record["CHROM"].replace("chr", "")

    is_canonical = chrom_name in CANONICALS

    if is_default:
        var_type = record["Type"]
        is_indel = (var_type == "insertion") or (var_type == "deletion")

        chrom_name = record["Chr"].replace("chr", "")
    else:
        is_indel = len(record["REF"]) != len(record["ALT"])

        chrom_name = record["CHROM"].replace("chr", "")

    is_canonical = chrom_name in CANONICALS
