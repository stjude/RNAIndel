#!/usr/bin/env python3

import sys
import pysam

canonicals = [str(i) for i in range(1, 22)] + ["X", "Y"]


def input_validator(alignments, genome, uniq_mapq):
    bam_canonical_chromosome_names = [
        i for i in alignments.references if is_canonical(i)
    ]
    fasta_canonical_chromosome_names = [i for i in genome.references if is_canonical(i)]

    if len(bam_canonical_chromosome_names) <= len(fasta_canonical_chromosome_names):
        shorter, longer = (
            bam_canonical_chromosome_names,
            fasta_canonical_chromosome_names,
        )
    else:
        shorter, longer = (
            fasta_canonical_chromosome_names,
            bam_canonical_chromosome_names,
        )

    if not set(shorter) <= set(longer):
        print(
            "Please input the same reference FASTA file used for mapping.",
            file=sys.stderr,
        )
        sys.exit(1)

    # select only good reads
    reads = sample_reads(alignments, 300)

    if reads:
        max_mapq = max([read.mapping_quality for read in reads])
        if max_mapq != uniq_mapq:
            reads = sample_reads(alignments, 30000)
            max_mapq = max([read.mapping_quality for read in reads])
            if max_mapq != uniq_mapq:
                print("The MAPQ for unique mapping looks like {}.".format(max_mapq))
                print("Please specify this MAPQ by -q (default: 255).", file=sys.stderr)
                sys.exit(1)

        has_md = sum([read.has_tag("MD") for read in reads])
        if has_md == 0:
            print("MD tag is missing in input BAM.", file=sys.stderr)
            sys.exit(1)
    else:
        print("Please check if the input BAM contains alignments.", file=sys.stderr)
        sys.exit(1)


def is_canonical(chrom_name):
    chrom_name_nochr = chrom_name.replace("chr", "")

    if chrom_name_nochr in canonicals:
        return True
    else:
        return False


def sample_reads(alignments, n):
    try:
        i = 0
        reads = []
        for read in alignments.fetch():
            if not read.is_secondary and not read.is_duplicate:
                reads.append(read)
                i += 1
            
            if i > n:
                break

        return reads
    except:
        return None
