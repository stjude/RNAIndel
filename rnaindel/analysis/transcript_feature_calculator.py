import re
import numpy as np

from .utils import most_common, flatten_list_of_list

numeric = re.compile(r"[0-9]+")


def transcript_features(df, proteindb):

    domain_dict = make_conserved_domain_dict(proteindb)

    (
        df["annotation"],
        df["cds_length"],
        df["indel_location"],
        df["is_inframe"],
        df["is_splice"],
        df["is_truncating"],
        df["is_nmd_insensitive"],
        df["is_in_cdd"],
        df["gene_symbol"],
    ) = zip(*df.apply(_wrapper, domain_dict=domain_dict, axis=1))

    dfg = df.groupby("gene_symbol")
    df = dfg.apply(indels_per_gene)

    df.drop(columns=["coding_indel_isoforms"], inplace=True)

    return df


def _wrapper(row, domain_dict):

    coding_indel_isoforms = row["coding_indel_isoforms"]

    annotations = get_annot_str(coding_indel_isoforms)
    cds_length = get_median_cds_length(coding_indel_isoforms)
    indel_location = get_median_indel_location(coding_indel_isoforms)
    is_frame = is_inframe_for_at_least_one_isoform(coding_indel_isoforms)
    is_splice = is_splice_for_at_least_one_isoform(coding_indel_isoforms)
    is_truncating = truncating_is_most_common(coding_indel_isoforms)
    is_nmd_insensitive = nmd_insensitive_is_most_common(coding_indel_isoforms)

    is_in_cdd = in_conserved_domain_is_most_common(coding_indel_isoforms, domain_dict)

    gene_symbol = get_gene_symbol(coding_indel_isoforms)

    return (
        annotations,
        cds_length,
        indel_location,
        is_frame,
        is_splice,
        is_truncating,
        is_nmd_insensitive,
        is_in_cdd,
        gene_symbol,
    )


def get_annot_str(coding_indel_isoforms):
    return ",".join([isoform.to_str() for isoform in coding_indel_isoforms])


def get_median_cds_length(coding_indel_isoforms):
    return np.median([isoform.cds_len for isoform in coding_indel_isoforms])


def get_median_indel_location(coding_indel_isoforms):
    return np.median(
        [isoform.get_relative_location() for isoform in coding_indel_isoforms]
    )


def is_inframe_for_at_least_one_isoform(coding_indel_isoforms):
    res = any(
        ["inframe" in isoform.variant_effect for isoform in coding_indel_isoforms]
    )
    return int(res)


def is_splice_for_at_least_one_isoform(coding_indel_isoforms):
    res = any(["splice" in isoform.variant_effect for isoform in coding_indel_isoforms])
    return int(res)


def truncating_is_most_common(coding_indel_isoforms):
    res = most_common(
        ["Truncating" in isoform.variant_effect for isoform in coding_indel_isoforms]
    )
    return int(res)


def nmd_insensitive_is_most_common(coding_indel_isoforms):
    return most_common(
        [isoform.is_nmd_insensitive() for isoform in coding_indel_isoforms]
    )


def make_conserved_domain_dict(proteindb):
    d = {}
    with open(proteindb) as f:
        for line in f:
            lst = line.split("\t")
            mrna_accession = lst[0]
            domains = [int(n) for n in re.findall(numeric, lst[1])]
            d[mrna_accession] = domains

    return d


def in_conserved_domain_is_most_common(coding_indel_isoforms, domain_dict):
    res = []
    for isoform in coding_indel_isoforms:
        res.append(is_within(isoform, domain_dict))
    return most_common(res)


def is_within(isoform, domain_dict):
    domain_position_lst = domain_dict.get(isoform.accession, [])
    if not domain_position_lst:
        return 0

    codon_pos = isoform.codon_pos
    domain_position_lst.append(codon_pos)
    domain_position_lst.sort()

    # if domain start <= codon_pos <= domain end
    # the index for aa_pos is odd, otherwise even
    idx = domain_position_lst.index(codon_pos)

    return idx % 2 == 1


def get_gene_symbol(coding_indel_isoforms):
    gene_symbols = [isoform.gene_symbol for isoform in coding_indel_isoforms]
    return most_common(gene_symbols)


def indels_per_gene(df_grouped_by_gene_symbol):
    num_of_indels = len(df_grouped_by_gene_symbol)

    list_of_coding_indel_isoforms = df_grouped_by_gene_symbol[
        "coding_indel_isoforms"
    ].to_list()

    # normalize by coding region len
    normalized = [
        (num_of_indels * 1000) / isoform.cds_len
        for isoform in flatten_list_of_list(list_of_coding_indel_isoforms)
    ]

    df_grouped_by_gene_symbol["ipg"] = np.median(normalized)

    return df_grouped_by_gene_symbol
