#!/usr/bin/env python3

import os
import sys
import pysam
import datetime
import pandas as pd
from collections import Counter

from .variant import Variant
from .variant import make_indel_from_vcf_line


def make_non_somatic_panel(file_lst, panelname, genome, cosmic_db, cnt):

    indel_lst = filter_indels(file_lst, genome, cosmic_db, cnt)
    vcf_data = to_vcf_data(indel_lst)

    with open(panelname, "w") as f:
        f.write(vcf_header() + "\n")
        f.write("\n".join(vcf_data))
   
    pysam.tabix_index(panelname, preset="vcf")

def vcf_header():
    dt = datetime.datetime.now()
    today = str(dt.year) + str(dt.month) + str(dt.day)
    header = [
        "##fileformat=VCFv4.2",
        "##filedate=" + today,
        '##INFO=<ID=OCCURRENCE,Number=1,Type=Integer,Description="Indel count in the dataset">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    return "\n".join(header)


def to_vcf_data(indel_lst):
    df = pd.DataFrame(
        [{"indel": indel, "occurrence": occurrence} for indel, occurrence in indel_lst]
    )
    df["CHROM"], df["POS"], df["ID"], df["REF"], df["ALT"], df["QUAL"], df[
        "FILTER"
    ], df["INFO"] = zip(*df.apply(populate_data, axis=1))

    df["CHROM"] = df.apply(lambda x: x["CHROM"].replace("chr", ""), axis=1)
    df["CHROM"] = df.apply(lambda x: 23 if x["CHROM"] == "X" else x["CHROM"], axis=1)
    df["CHROM"] = df.apply(lambda x: 24 if x["CHROM"] == "Y" else x["CHROM"], axis=1)
    df["CHROM"] = df.apply(lambda x: int(x["CHROM"]), axis=1)

    df.sort_values(["CHROM", "POS"], inplace=True)

    df["CHROM"] = df.apply(lambda x: "Y" if x["CHROM"] == 24 else x["CHROM"], axis=1)
    df["CHROM"] = df.apply(lambda x: "X" if x["CHROM"] == 23 else x["CHROM"], axis=1)
    df["CHROM"] = df.apply(lambda x: "chr" + str(x["CHROM"]), axis=1)

    df["VCF"] = df.apply(
        lambda x: "\t".join(
            [
                x["CHROM"],
                str(x["POS"]),
                x["ID"],
                x["REF"],
                x["ALT"],
                x["QUAL"],
                x["FILTER"],
                x["INFO"],
            ]
        ),
        axis=1,
    )

    return df["VCF"].values


def populate_data(row):
    var, occurrence = row["indel"], row["occurrence"]
    return (
        var.chrom,
        var.pos,
        ".",
        var.ref,
        var.alt,
        ".",
        "PASS",
        "OCCURRENCE=" + str(occurrence),
    )


def filter_indels(file_lst, genome, cosmic_db, cnt):
    vcf_lst = validate_file_lst(file_lst)
    processed_vcfs = [process_vcf_file(vcf, genome) for vcf in vcf_lst]
    indel_lst = to_flat_lst(processed_vcfs)
    frequent_indels = [(k, v) for k, v in Counter(indel_lst).items() if v >= cnt]
    non_cosmic_frequent_indels = [
        (k, v) for k, v in frequent_indels if is_absent_in_cosmic(k, genome, cosmic_db)
    ]

    return non_cosmic_frequent_indels


def to_flat_lst(lst_of_lst):
    return [elem for lst in lst_of_lst for elem in lst]


def validate_file_lst(filename):
    path_lst = []
    with open(filename) as f:
        for line in f:
            line = line.rstrip()
            if os.path.isfile(line):
                path_lst.append(line)
            else:
                print("Not found: " + line)

    if path_lst:
        return path_lst
    else:
        print("Input files not found.")
        sys.exit(1)


def process_vcf_file(vcf, genome):
    f = open(vcf)
    parsed_lines = [
        make_indel_from_vcf_line(line, genome)
        for line in f
        if make_indel_from_vcf_line(line, genome)
    ]
    f.close()
    return to_flat_lst(parsed_lines)


def is_absent_in_cosmic(var, genome, cosmic_db, fuzzy_match=False):
    window = 100
    search_space = cosmic_db.fetch(
        var.chrom.replace("chr", ""), var.pos - window / 2, var.pos + window / 2
    )

    for line in search_space:
        for db_indel in make_indel_from_vcf_line(line, genome):
            # fuzzy allows positional match
            if fuzzy_match:
                var_norm = var.normalize()
                db_indel_norm = db_indel.normalize()
                if var_norm.pos == db_indel_norm.pos:
                    return False
            else:
                if var == db_indel:
                    return False
    return True
