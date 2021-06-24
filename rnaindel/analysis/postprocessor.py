import os
import pysam
from .outlier import outlier_analysis


def postprocess(df, data_dir, perform_outlier_analysis):
    path_to_cosmic = "{}/cosmic/CosmicCodingMuts.indel.vcf.gz".format(data_dir)
    path_to_non_somatic = "{}/non_somatic/non_somatic.vcf.gz".format(data_dir)

    non_somatic = pysam.VariantFile(path_to_non_somatic)
    cosmic = pysam.VariantFile(path_to_cosmic)

    df["filter"], df["reclassified"], df["predicted_class"] = zip(
        *df.apply(_wrapper, non_somatic_db=non_somatic, cosmic=cosmic, axis=1)
    )

    df["is_rescurable_homopolymer"] = df.apply(is_rescurable_homopolymer, axis=1)

    if perform_outlier_analysis:
        df = outlier_analysis(df, os.path.join(data_dir, "outliers"))

    return sort_positionally(df)


def _wrapper(row, non_somatic_db, cosmic, mapping_thresh=0.5):
    fltr_str = filter_str(row, non_somatic_db, mapping_thresh)
    is_rescued = reclassify_by_knowledge(row, cosmic)
    if is_rescued:
        return fltr_str, "reclassified_by_knowledge", "somatic"
    else:
        return fltr_str, row["reclassified"], row["predicted_class"]


def is_rescurable_homopolymer(row):
    if row["is_common"]:
        return False

    if row["filter"] != "PASS":
        return False

    if row["reclassified"] != "-":
        return False

    if row["predicted_class"] == "somatic":
        return False
    
    if row["indel_size"] == 1 and row["repeat"] >= 5 and row["prob_s"] >= 0.2:
        vaf = row["alt_count"] / (row["ref_count"] + row["alt_count"])

        if (vaf > 0.3 and row["alt_count"] > 7) or vaf > 0.6:
            return True

    return False


def filter_str(row, non_somatic_db, mapping_thresh):

    pred = row["predicted_class"]

    if pred == "somatic":
        s = [
            filter_by_db(row, non_somatic_db),
            filter_by_mappability(row, mapping_thresh),
        ]
        if any(s):
            return ",".join(s).strip(",")
        else:
            return "PASS"
    else:
        # TODO from other caller
        return "PASS"


def filter_by_db(row, non_somatic_db):

    non_somatic_hits = row["indel"].query_vcf(non_somatic_db)

    if non_somatic_hits:
        if row["prob_a"] >= row["prob_g"]:
            return "ProbableArtifact"
        else:
            return "ProbableGermline"

    if row["is_common"] and not row["is_pathogenic"]:
        return "ProbableGermline"

    return ""


def filter_by_mappability(row, mapping_thresh):
    if row["uniq_mapping_rate"] < mapping_thresh:
        return "LowMappabilityRegion"
    else:
        return ""


def reclassify_by_knowledge(row, cosmic):
    if row["is_common"]:
        return False

    # known pathogenic event with high cosmic count
    cosmic_cnts = row["cosmic_cnt"]
    if row["prob_s"] > 0.1:
        if cosmic_cnts >= 10 and row["is_pathogenic"]:
            return True
        elif cosmic_cnts >= 30:
            return True

    return False


def sort_positionally(df):
    df["chrom"] = df.apply(lambda x: x["chrom"].replace("chr", ""), axis=1)
    df["chrom"] = df.apply(lambda x: 23 if x["chrom"] == "X" else x["chrom"], axis=1)
    df["chrom"] = df.apply(lambda x: 24 if x["chrom"] == "Y" else x["chrom"], axis=1)
    df["chrom"] = df.apply(lambda x: int(x["chrom"]), axis=1)

    df.sort_values(["chrom", "pos"], inplace=True)

    df["chrom"] = df.apply(lambda x: "Y" if x["chrom"] == 24 else x["chrom"], axis=1)
    df["chrom"] = df.apply(lambda x: "X" if x["chrom"] == 23 else x["chrom"], axis=1)
    df["chrom"] = df.apply(lambda x: "chr" + str(x["chrom"]), axis=1)

    return df
