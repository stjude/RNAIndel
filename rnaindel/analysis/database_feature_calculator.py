import pysam


def database_features(df, dbsnp, clinvar, cosmic):
    dbsnp = pysam.VariantFile(dbsnp)
    clinvar = pysam.VariantFile(clinvar)
    cosmic = pysam.VariantFile(cosmic)

    (
        df["dbsnp"],
        df["pop_freq"],
        df["is_common"],
        df["is_on_db"],
        df["is_pathogenic"],
        df["cosmic_cnt"],
    ) = zip(*df.apply(_wrapper, dbsnp=dbsnp, clinvar=clinvar, cosmic=cosmic, axis=1))

    return df


def _wrapper(row, dbsnp, clinvar, cosmic):
    variant = row["indel"]

    dbsnp_annotation = dbsnp_annot(variant, dbsnp)
    try:
        clin_annotation = clinvar_annot(variant, clinvar)
    except:
        clin_annotation = -1, 0

    # dbSNP ID
    rs_id = dbsnp_annotation[0]

    pop_freq = max(dbsnp_annotation[1], clin_annotation[0])

    is_common = _is_common(pop_freq, dbsnp_annotation)

    is_on_db = 0
    if is_common:
        is_on_db = 1
    elif dbsnp_annotation[-1]:
        is_on_db = 1

    cosmic_cnt = cosmic_annot(variant, cosmic)

    return rs_id, pop_freq, is_common, is_on_db, clin_annotation[1], cosmic_cnt


def get_population_freq(dbsnp_annotation, clin_annotation):

    max_dbsnp_freq = dbsnp_annotation[1] if dbsnp_annotation else -1
    max_clinvar_freq = clin_annotation[0] if clin_annotation else -1

    return max(max_dbsnp_freq, max_clinvar_freq)


def _is_common(pop_freq, dbsnp_annotation):
    if 1 in dbsnp_annotation[2]:
        return 1
    else:
        if pop_freq >= 0.01:
            return 1
        else:
            return 0


def dbsnp_annot(variant, db):
    hits = variant.query_vcf(db)

    if hits:
        ids = ",".join([hit["ID"] for hit in hits])

        max_thousand_genome_freq = max(
            [get_allele_freq(hit["INFO"], preset="CAF") for hit in hits]
        )
        max_topmed_freq = max(
            [get_allele_freq(hit["INFO"], preset="TOPMED") for hit in hits]
        )
        max_non_cancer_freq = max(
            [get_allele_freq(hit["INFO"], preset="non_cancer_AF") for hit in hits]
        )

        max_freq = max(max_thousand_genome_freq, max_topmed_freq, max_non_cancer_freq)

        dbsnp_common = [int(hit["INFO"].get("COMMON", -1)) for hit in hits]
        is_common_in_non_cancer_pop = max_non_cancer_freq > 0.0001

        return ids, max_freq, dbsnp_common, is_common_in_non_cancer_pop
    else:
        return ".", -1, [], False


def clinvar_annot(variant, db):
    hits = variant.query_vcf(db)

    if hits:
        max_esp_freq = max(
            [get_allele_freq(hit["INFO"], preset="AF_ESP") for hit in hits]
        )
        max_exac_freq = max(
            [get_allele_freq(hit["INFO"], preset="AF_EXAC") for hit in hits]
        )
        max_tgp_freq = max(
            [get_allele_freq(hit["INFO"], preset="AF_TGP") for hit in hits]
        )

        max_freq = max(max_esp_freq, max_exac_freq, max_tgp_freq)

        clinical_signigicance = ",".join(
            [hit["INFO"].get("CLNSIG", "xxx").upper() for hit in hits]
        )
        is_pathogenic = 1 if "PATHOGENIC" in clinical_signigicance else 0

        return max_freq, is_pathogenic
    else:
        return -1, 0


def cosmic_annot(variant, db):
    hits = variant.query_vcf(db)

    if hits:
        cnts = sum([int(hit["INFO"]["CNT"]) for hit in hits])

        return cnts
    else:
        return 0


def get_allele_freq(vcf_info, preset):
    freq_str = str(vcf_info.get(preset, "-1.0,-1.0"))

    freq = -1
    try:
        if preset == "non_cancer_AF":
            freq = max([float(freq) for freq in freq_str.split(",")])
        elif preset in ["CAF", "TOPMED"]:
            major_freq = float(freq_str.split(",")[0])
            if major_freq < 0:
                return major_freq
            else:
                return 1.0 - major_freq
        elif preset in ["AF_ESP", "AF_EXAC", "AF_TGP"]:
            freq = max([float(freq) for freq in freq_str.split(",")])
    except:
        pass

    return freq
