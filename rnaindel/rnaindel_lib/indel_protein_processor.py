#!/usr/bin/env python3
"""4th step of analysis

Calculate features at protein level

'indel_protein_processor' is the main routine of this module
"""

import re
import gzip
import numpy as np
from functools import partial

numeric = re.compile(r"[0-9]+")


def indel_protein_processor(df, refgene, proteincdd):
    """Calculate protein features
     
    Features not used in the final model are commented out

    Args:
        df (pandas.DataFrame)
        refgene (str): path to refCodingExon.bed.gz
        proteincdd (str): optional, path to proteinConservedDomains.txt
    Returns:
        df (pandas.DataFrame)
    """
    # cds length & indel location
    acc_len = acc_len_dict(refgene)
    df["cds_length"], df["indel_location"] = zip(
        *df.apply(partial(len_loc, d=acc_len), axis=1)
    )

    # check if the indel is in conserved domain (CDD)
    acc_dom = acc_domain_dict(proteincdd)
    df["is_in_cdd"] = df.apply(partial(is_in_conserved_domain, d=acc_dom), axis=1)

    return df


def acc_len_dict(refgene):
    """Making a dictionary {accession: CDS_length}

    Args:
        refgene (str): path to refCodingExon.bed.gz
    Returns:
        d (dict): {accession (str): coding_seq_length (int)}
    """
    d = {}
    with gzip.open(refgene, "rb") as f:
        for line in f:
            lst = line.split(b"\t")
            info = lst[3].split(b"|")
            acc = info[0].decode("utf-8")
            cds_len = info[5].decode("utf-8")

            d[acc] = int(cds_len)

    return d


def len_loc(row, d):
    """Calculating median CDS length and indel location

    Args:
        row (pandas.Series): a Series with 'annotation' index
        d (dict): {accession (str): CDS_length (int)}
    Returns:
        tuple (float) : 1st element = median CDS len over all isoforms
                        2nd element = median indel location over all isofroms
    """
    median_cds_len = 1323
    tokens = row["annotation"].split(",")

    lengths = []
    locations = []
    for token in tokens:
        info = token.split("|")
        acc = info[1]
        cds_pos = int(info[2]) * 3
        cds_len = d.get(acc, median_cds_len)

        lengths.append(cds_len)
        locations.append(cds_pos / cds_len)

    return np.median(lengths), np.median(locations)


def acc_domain_dict(proteincdd):
    """Making a dictionary {accession: domain positions}

    Args:
        proteincdd (str): path to proteinConservedDomains.txt
    Returns:
        d (dict): {accession (str): list[domain_start(int), domain_end(int)]}
    """
    d = {}
    with open(proteincdd) as f:
        for line in f:
            lst = line.split("\t")
            acc = lst[0]
            domains = [int(n) for n in re.findall(numeric, lst[1])]
            d[acc] = domains

    return d


def is_in_conserved_domain(row, d):
    """Annotating whether indel is in conserved domain

    Args:
        row (pandas.Series): a Series with 'annotation' index
        d (dict): {accession (str): list[domain_start(int), domain_end(int)]}
    Returns:
        is_in_cdd (bool): 1 if in conserved domain, 0 otherwise
    """
    tokens = row["annotation"].split(",")

    is_in_domain = []
    for token in tokens:
        info = token.split("|")
        acc = info[1]
        aa_pos = int(info[2])
        try:
            domain_pos = d[acc]

            # add aa_pos and sort
            domain_pos.append(aa_pos)
            domain_pos.sort()

            # if domain start <= aa_pos <= domain end
            # the index for aa_pos is odd, otherwise even
            idx = domain_pos.index(aa_pos)
            if idx % 2 == 1:
                is_in_domain.append(1)
            else:
                is_in_domain.append(0)
        except:
            is_in_domain.append(0)

    # most frequent pattern
    is_in_cdd = round(np.mean(is_in_domain))

    return is_in_cdd
