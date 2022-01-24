#!/usr/bin/env python3

import os
import sys
import pysam
import argparse
from functools import partial
from collections import Counter
from indelpost import Variant


# def count(file_lst, genome, out_dir):
def count():

    subcommand = "CountOccurrence"
    args = get_args(subcommand)

    file_lst = args.vcf_list
    genome = pysam.FastaFile(args.reference)
    out_dir = args.out_dir

    validated_vcfs = validate_vcfs(file_lst)
    occurrence_dict = count_somatic_preditions(validated_vcfs, genome)

    for vcf in validated_vcfs:
        annotate_vcf_with_recurrence(vcf, genome, occurrence_dict, out_dir)


def count_somatic_preditions(validated_vcfs, genome):
    processed_vcfs = [
        collect_somatic_predictions_from_vcf(vcf, genome) for vcf in validated_vcfs
    ]
    somatic_indels = to_flat_lst(processed_vcfs)
    return Counter(somatic_indels)


def to_flat_lst(lst_of_lst):
    return [elem for lst in lst_of_lst for elem in lst]


def validate_vcfs(file_lst):
    vcf_lst = validate_file_lst(file_lst)
    return [vcf for vcf in vcf_lst if is_rnaindel_output_vcf(vcf)]


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


def is_rnaindel_output_vcf(vcf):
    headers = [line for line in open(vcf) if line.startswith("##source=RNAIndel")]
    if headers:
        return True
    else:
        print(vcf + ": validation failed. Check if this is a RNAIndel output VCF.")
        return False


def make_indel_from_vcf_line(line, genome):
    if line.startswith("#"):
        return None

    lst = line.rstrip().split("\t")
    chrom, pos = lst[0], int(lst[1])
    ref, alts = lst[3], lst[4].split(",")

    if all(len(ref) == len(alt) == 1 for alt in alts):
        return None

    indels = [
        Variant(chrom, pos, ref, alt, genome)
        for alt in alts
        if Variant(chrom, pos, ref, alt, genome).variant_type != "S"
    ]

    if indels:
        return indels
    else:
        return None


def collect_somatic_predictions_from_vcf(vcf, genome):
    somatic_lines = [line for line in open(vcf) if "predicted_class=somatic" in line]
    parsed_somatic_lines = [
        make_indel_from_vcf_line(line, genome)
        for line in somatic_lines
        if make_indel_from_vcf_line(line, genome)
    ]
    return to_flat_lst(parsed_somatic_lines)


def annotate_vcf_with_recurrence(vcf, genome, occurrence_dict, outdir):
    new_vcf = edit_header(vcf)

    fi = open(vcf)
    for line in fi:
        if line.startswith("#"):
            pass
        elif "predicted_class=somatic" in line:
            indels = make_indel_from_vcf_line(line, genome)
            if indels:
                putative_somatic = indels[0]
                occurrence = occurrence_dict[putative_somatic]
                new_vcf.append(append_recurrence(line, occurrence) + "\n")
        else:
            new_vcf.append(line)
    fi.close()

    if outdir:
        filename = os.path.basename(vcf).replace(".vcf", ".occurrence.vcf")
        fo = open(os.path.join(outdir, filename), "w")
    else:
        fo = open(vcf, "w")
    fo.write("".join(new_vcf))
    fo.close()


def edit_header(vcf):
    header_lines = [line for line in open(vcf) if line.startswith("##")]
    new_header_line = '##INFO=<ID=occurrence,Number=1,Type=Integer,Description="Occurrence count in the input VCF files. Only counted for indels predicted as somatic">\n'
    bottom = [line for line in open(vcf) if line.startswith("#CHROM")]

    if new_header_line not in header_lines:
        return header_lines + [new_header_line] + bottom
    else:
        return header_lines + bottom


def append_recurrence(line, occurrence):
    lst = line.rstrip().split("\t")
    info = (
        lst[7] if "occurrence=" in lst[7] else lst[7] + ";occurrence=" + str(occurrence)
    )
    new_lst = lst[0:7] + [info] + lst[8:]
    return "\t".join(new_lst)


def get_args(subcommand):
    prog = "rnaindel " + subcommand
    parser = argparse.ArgumentParser(prog=prog)

    parser.add_argument(
        "-r",
        "--reference",
        metavar="FILE",
        required=True,
        type=validate_file_input,
        help="reference genome FASTA file.",
    )

    parser.add_argument(
        "--vcf-list",
        metavar="FILE",
        required=True,
        type=validate_file_input,
        help="file containing paths to RNAIndel output VCF files to be annotated",
    )

    parser.add_argument(
        "--out-dir",
        metavar="DIR",
        default=None,
        type=validate_dir_input,
        help="directory to ouput (default: input file directory)",
    )

    return parser.parse_args(sys.argv[2:])


def validate_file_input(file_path):
    if os.path.isfile(file_path):
        return file_path
    else:
        sys.exit("Error: {} not found.".format(file_path))


def validate_dir_input(dir_path):
    if os.path.isdir(dir_path):
        return dir_path
    else:
        sys.exit("Error: {} not found.".format(dir_path))
