#!/usr/bin/env python3
"""Optional step for reclassification.

Knowledge-based reclassification is performed
for indels predicted somatic.

'indel_reclassifier' is the main routine in this module
"""
import pysam
from functools import partial
from .indel_snp_annotator import vcf2bambino
from .indel_curator import curate_indel_in_genome


def indel_reclassifier(df, genome, default_pons, user_pons, cosmic, chr_prefixed):
    """Knowledge-based reclassification
    
    Args:
        df (pandas.DataFrame): df with prediction made
        genome (pysam.FastaFile): reference genome 
        default_pons (pysam.TabixFile): default black list indels in .vcf.gz
        user_pons (pysam.TabixFile): user's panel of non-somatic indels in .vcf.gz. may be None
        cosmic (pysam.TabixFile): COSMIC db 
        chr_prefixed (bool): True if chromosome names in BAM are "chr"-prefixed
    Returns:
        df (pandas.DataFrame): df reclassified
    """

    df["predicted_class"], df["reclassified"] = zip(
        *df.apply(
            reclassify_by_pons,
            genome=genome,
            pons=default_pons,
            chr_prefixed=chr_prefixed,
            axis=1,
        )
    )

    if user_pons:
        df["predicted_class"], df["reclassified"] = zip(
            *df.apply(
                reclassify_by_pons,
                genome=genome,
                pons=user_pons,
                chr_prefixed=chr_prefixed,
                axis=1,
            )
        )

    df["predicted_class"], df["reclassified"] = zip(
        *df.apply(
            reclassify_by_cosmic,
            genome=genome,
            cosmic=cosmic,
            chr_prefixed=chr_prefixed,
            axis=1,
        )
    )
    df["predicted_class"], df["reclassified"] = zip(
        *df.apply(rescue_highquality_homopolymer, cosmic=cosmic, axis=1)
    )

    df["predicted_class"], df["reclassified"] = zip(
        *df.apply(reclassify_by_mappability, axis=1)
    )

    return df


def reclassify_by_pons(row, genome, pons, chr_prefixed, preset="pons"):
    """Reclassify by blacklist

    Args: see 'relassify_by_db'
    Returns: see 'relassify_by_db'
    """
    if (
        row["predicted_class"] == "somatic"
        and row["is_common"] != 1
        and row["reclassified"] == "-"
    ):
        return relassify_by_db(row, genome, pons, chr_prefixed, preset)
    else:
        return row["predicted_class"], row["reclassified"]


def reclassify_by_cosmic(row, genome, cosmic, chr_prefixed, preset="cosmic"):
    """Reclassify by whitelist

    Args: see 'relassify_by_db'
    Returns: see 'relassify_by_db'
    """
    if (
        row["predicted_class"] != "somatic"
        and row["prob_a"] < 0.6
        and row["reclassified"] == "-"
    ):
        return relassify_by_db(row, genome, cosmic, chr_prefixed, preset)
    else:
        return row["predicted_class"], row["reclassified"]


def rescue_highquality_homopolymer(row, cosmic):
    """Lower threshold for high-quality homopolymer with COSMIC support
    """
    if (
        row["indel_size"] > 1
        or row["repeat"] < 6
        or row["prob_s"] < 0.33
        or row["is_on_db"] == 1
        or row["reclassified"] == "reclassified"
    ):
        return row["predicted_class"], row["reclassified"]

    if (
        row["alt_count"] < 5
        or row["alt_count"] / (row["ref_count"] + row["alt_count"]) < 0.5
    ):
        return row["predicted_class"], row["reclassified"]

    near_20_aminoacids = cosmic.fetch(
        row["chr"].replace("chr", ""), row["pos"] - 60, row["pos"] + 60
    )
    try:
        cosmic_cnt_info = [
            int(str(record).split("\t")[7].split("CNT=")[1])
            for record in near_20_aminoacids
            if "CNT=" in str(record)
        ]
        if sum(cosmic_cnt_info) >= 3:
            return "somatic", "reclassified"
        else:
            return row["predicted_class"], row["reclassified"]
    except:
        return row["predicted_class"], row["reclassified"]


def relassify_by_db(row, genome, db, chr_prefixed, preset):
    """Knowledge-based reclassification of indels predicted somatic

    Args:
        row (pandas.Series)
        genome (pysam.FastaFile): reference genome
        db (pysam.TabixFile obj): indel db in tabixed VCF
        chr_prefixed (bool): True if chromosome names in BAM are "chr"-prefixed
    Returns:
        'predicted_class' (str): reclassifed class if applicable
        'comment' (str): 'reclassified' if reclassified, '-' otherwise 
    """
    search_window = 50
    chr = row["chr"]
    pos = row["pos"]
    idl_type = row["is_ins"]
    idl_seq = row["indel_seq"]

    idl = curate_indel_in_genome(genome, chr, pos, idl_type, idl_seq, chr_prefixed)

    # check if contif names in vcf are prefixed with 'chr'
    sample_contig = db.contigs[0]
    if not sample_contig.startswith("chr"):
        chr_vcf = chr.replace("chr", "")
    else:
        chr_vcf = chr

    start, end = pos - search_window, pos + search_window

    # check if the indel is equivalent to indel on database
    for record in db.fetch(chr_vcf, start, end, parser=pysam.asTuple()):
        bambinos = vcf2bambino(record)
        for bb in bambinos:
            if idl_type == bb.idl_type and len(idl_seq) == len(bb.idl_seq):
                db_idl = curate_indel_in_genome(
                    genome, chr, bb.pos, bb.idl_type, bb.idl_seq, chr_prefixed
                )

                if idl == db_idl:
                    if preset == "pons":
                        if row["prob_a"] >= row["prob_g"]:
                            return "artifact", "reclassified"
                        else:
                            return "germline", "reclassified"
                    elif preset == "cosmic":
                        cosmic_info = str(record).split("\t")[7].split(";")
                        cosmic_cnt_data = [
                            i.replace("CNT=", "") for i in cosmic_info if "CNT=" in i
                        ]
                        if cosmic_cnt_data:
                            cosmic_cnt = int(cosmic_cnt_data[0])
                            if (
                                "Pathogenic" in row["clin_info"]
                                or "Likely Pathogenic" in row["clin_info"]
                            ) and cosmic_cnt >= 3:
                                return "somatic", "reclassified"
                            elif row["is_on_db"] == 0 and cosmic_cnt >= 20:
                                return "somatic", "reclassified"
                            else:
                                pass
                    else:
                        pass

    return row["predicted_class"], row["reclassified"]


def reclassify_by_mappability(row):
    if (
        row["predicted_class"] == "somatic"
        and row["reclassified"] != "reclassified"
        and row["mappability"] < 0.5
    ):
        return "artifact", "reclassified"
    else:
        return row["predicted_class"], row["reclassified"]
