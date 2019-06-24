#!/usr/bin/env python3

import pysam
from .make_non_somatic_panel import to_flat_lst
from .variant import make_indel_from_vcf_line


def filterate_by_panel(vcf, outvcf, genome, nonsomatic_panel):
    new_vcf = edit_header(vcf, nonsomatic_panel)

    nonsomatic_db = pysam.TabixFile(nonsomatic_panel)
    chr_prefixed = is_chr_prefixed(nonsomatic_db)

    with open(vcf) as f:
        for line in f:
            if line.startswith("#"):
                pass
            elif is_somatic_prediction(line):
                putative_somatic = make_indel_from_vcf_line(line, genome)[0]
                candidates = extract_candidate_indels(line, genome, nonsomatic_db, chr_prefixed)
                if candidates and putative_somatic in candidates:
                    new_vcf.append(reclassify(line))
                else:
                    new_vcf.append(line)
            else:
                new_vcf.append(line)

    output_vcf = open(outvcf, "w")
    output_vcf.write("".join(new_vcf))
    output_vcf.close()


def edit_header(vcf, nonsomatic_panel):
    header_lines = [line for line in open(vcf) if line.startswith("##")]
    new_header_line = "##reclassified_by=" + nonsomatic_panel + "\n"
    bottom = [line for line in open(vcf) if line.startswith("#CHROM")]

    return header_lines + [new_header_line] + bottom


def extract_candidate_indels(line, genome, nonsomatic_db, is_chr_prefixed):
    chrom, pos = line.split("\t")[0].replace("chr", ""), int(line.split("\t")[1])
    chrom = "chr" + chrom if is_chr_prefixed else chrom

    sliced_db = nonsomatic_db.fetch(chrom, pos - 50, pos + 50)
    lst_of_indel_lst = [
        make_indel_from_vcf_line(line, genome) for line in sliced_db
    ]
    
    return to_flat_lst(lst_of_indel_lst)


def is_chr_prefixed(nonsomatic_db):
    return nonsomatic_db.contigs[0].startswith("chr")


def is_somatic_prediction(line):
    return "PRED=somatic" in line


def append_rcf_to_info(line):
    lst = line.rstrip().split("\t")
    info = lst[7] + ";RCF"
    new_lst = lst[0:7] + [info] + lst[8:]
    return "\t".join(new_lst)


def reclassify(line):
    line = line.rstrip()
    probabilities = [
        i.replace("PROB=", "") for i in line.split("\t")[7].split(";") if "PROB=" in i
    ][0].split(",")
    
    germline_prob, artifact_porb = float(probabilities[1]), float(probabilities[2])
    if artifact_porb > germline_prob:
        line = line.replace("PRED=somatic", "PRED=artifact")
        line = append_rcf_to_info(line) if not "RCF" in line else line        

    else:
        line = line.replace("PRED=somatic", "PRED=germline")
        line = append_rcf_to_info(line) if not "RCF" in line else line
    return line + "\n"
