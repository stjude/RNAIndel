#!/usr/bin/env python3

import re
import os
from functools import partial
from collections import Counter
from .variant import make_indel_from_vcf_line
from .filterate_by_panel import is_somatic_prediction
from .make_non_somatic_panel import to_flat_lst
from .make_non_somatic_panel import validate_file_lst
from .make_non_somatic_panel import is_absent_in_cosmic


def annotate_recurrence(file_lst, genome, cosmic_db, out_dir):
    validated_vcfs = validate_vcfs(file_lst)
    occurrence_dict = count_somatic_preditions(validated_vcfs, genome)

    for vcf in validated_vcfs:
        annotate_vcf_with_recurrence(vcf, genome, cosmic_db, occurrence_dict, out_dir)


def count_somatic_preditions(validated_vcfs, genome):
    processed_vcfs = [
        collect_somatic_predictions_from_vcf(vcf, genome) for vcf in validated_vcfs
    ]
    somatic_indels = to_flat_lst(processed_vcfs)
    return Counter(somatic_indels)


def validate_vcfs(file_lst):
    vcf_lst = validate_file_lst(file_lst)
    return [vcf for vcf in vcf_lst if is_rnaindel_output_vcf(vcf)]


def is_rnaindel_output_vcf(vcf):
    headers = [line for line in open(vcf) if line.startswith("##source=RNAIndel")]
    if headers:
        return True
    else:
        print(vcf + ": validation failed. Check if this is a RNAIndel output VCF.")
        return False


def collect_somatic_predictions_from_vcf(vcf, genome):
    somatic_lines = [line for line in open(vcf) if "PRED=somatic" in line]
    parsed_somatic_lines = [
        make_indel_from_vcf_line(line, genome)
        for line in somatic_lines
        if make_indel_from_vcf_line(line, genome)
    ]
    return to_flat_lst(parsed_somatic_lines)


def annotate_vcf_with_recurrence(vcf, genome, cosmic_db, occurrence_dict, outdir):
    new_vcf = edit_header(vcf)

    fi = open(vcf)
    for line in fi:
        if line.startswith("#"):
            pass
        elif is_somatic_prediction(line):
            line = re.sub(r"REC=[0-9]+", "", line) 
            putative_somatic = make_indel_from_vcf_line(line, genome)[0]
            occurrence = occurrence_dict[putative_somatic]
            if occurrence > 1:
                new_vcf.append(append_recurrence(line, occurrence) + "\n")
            else:
                new_vcf.append(line)
        else:
            new_vcf.append(line)
    fi.close()

    if outdir:
        filename = os.path.basename(vcf)
        fo = open(os.path.join(outdir, filename), "w")
    else:
        fo = open(vcf, "w")
    fo.write("".join(new_vcf))
    fo.close()


def edit_header(vcf):
    header_lines = [line for line in open(vcf) if line.startswith("##")]
    new_header_line = '##INFO=<ID=REC,Number=1,Type=Integer,Description="Recurrent somatic predictions. Annotated as the number of recurrent samples in which the indel is predicted as somatic in the input cohort. Not annotated for artifact, germline and non-recurrent somatic predictions">\n'
    bottom = [line for line in open(vcf) if line.startswith("#CHROM")]
    return header_lines + [new_header_line] + bottom


def append_recurrence(line, occurrence):
    lst = line.rstrip().split("\t")
    info = lst[7] if "REC=" in lst[7] else lst[7] + ";REC=" + str(occurrence)
    new_lst = lst[0:7] + [info] + lst[8:]
    return "\t".join(new_lst)
