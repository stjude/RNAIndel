import pysam


def postprocess(df, data_dir, user_db=None):
    path_to_cosmic = "{}/cosmic/CosmicCodingMuts.indel.vcf.gz".format(data_dir)
    path_to_non_somatic = "{}/non_somatic/non_somatic.vcf.gz".format(data_dir)

    non_somatic = pysam.VariantFile(path_to_non_somatic)
    cosmic = pysam.VariantFile(path_to_cosmic)

    df["filter"], df["reclassified"] = zip(
        *df.apply(_wrapper, non_somatic_db=non_somatic, cosmic=cosmic, axis=1)
    )

    return sort_positionally(df)


def _wrapper(row, non_somatic_db, cosmic, mapping_thresh=0.5):
    fltr_str = filter_str(row, non_somatic_db, mapping_thresh)
    is_rescued = rescue_possible_misclassified(row, cosmic)
    if is_rescued:
        return fltr_str, "reclassified"
    else:
        return fltr_str, row["reclassified"]


def filter_str(row, non_somatic_db, mapping_thresh):
    
    pred = row["predicted_class"]
    
    if pred == "somatic":
        s = [filter_by_db(row, non_somatic_db), filter_by_mappability(row, mapping_thresh)]
        if any(s):
            return ",".join(s).strip(",")
        else:
            return "PASS"
    else:
        #TODO from other caller
        return "PASS"


def filter_by_db(row, non_somatic_db):

    non_somatic_hits = row["indel"].query_vcf(non_somatic_db)
    
    if non_somatic_hits:
        if row["prob_a"] >= row["prob_g"]:
            return "ProbableArtifact"
        else:
            return "ProbableGermline"
    
    if row["is_common"] and not row["is_pathogenic"]:
        return  "ProbableGermline"

    return ""


def filter_by_mappability(row, mapping_thresh):
    if row["uniq_mapping_rate"] < mapping_thresh:
        return "LowMappabilityRegion"
    else:
        return ""


def rescue_possible_misclassified(row, cosmic):

    # high_quality homopolymer with near cosmic records
    if row["indel_size"] == 1 and row["repeat"] >= 6 and row["prob_s"] >= 0.33:
        near_by_cosmic_indels = cosmic.fetch(
            row["chrom"].replace("chr", ""), row["pos"] - 30, row["pos"] + 30
        )
        if near_by_cosmic_indels:
            return True

    # known cosmics and pathogenics
    cosmic_cnts = row["cosmic_cnt"]
    is_on_db = row["is_on_db"]
    if cosmic_cnts >= 3 and row["is_pathogenic"]:
        return True
    elif cosmic_cnts >= 10 and not is_on_db:
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
