#!/usr/bin/env python3
"""2nd step of the analysis

Checks if indels are coding or non-coding and annotates
coding indels with variant effect 

indel_annotator is the main routine of this module
"""

import sys
import logging
import pandas as pd
from .indel_curator import curate_indel_in_genome
from .indel_sequence import CodingSequenceWithIndel

logger = logging.getLogger(__name__)


def indel_annotator(df, genome, exons, chr_prefixed):
    """Sort coding indels and annotate coding indels with variant effect

    Args:
        df (pandas.DataFrame): with a header:'chr', 'pos', 'ref', 'alt' 
        genome (pysam.FastaFile): reference genome
        exons (pysam.TabixFile): coding exon data
    Returns:
        df (pandas.DataFrame): with indels annotated
    """
    df["is_ins"] = df.apply(lambda x: 1 if x["ref"] == "-" else 0, axis=1)
    df["indel_seq"] = df.apply(lambda x: x["alt"] if x["is_ins"] else x["ref"], axis=1)

    # annotate coding indels
    df["annotation"] = df.apply(
        annotate_indels, genome=genome, exons=exons, chr_prefixed=chr_prefixed, axis=1
    )

    # annotate non-coding indels equivalent to a coding indel (occurs in 3' splice and 3'UTR)
    df["annotation"] = df.apply(
        annotate_equivalent_non_coding_indels,
        df=df,
        genome=genome,
        chr_prefixed=chr_prefixed,
        axis=1,
    )

    # removes unannotated calls
    df = df[df["annotation"] != "-"]

    if len(df) == 0:
        logging.warning("No indels annotated in coding region. Analysis done.")
        sys.exit(0)

    # gene symbols
    df["gene_symbol"] = df.apply(get_gene_symbol, axis=1)

    # formats the header
    df = df[
        [
            "chr",
            "pos",
            "ref",
            "alt",
            "rescued",
            "indel_seq",
            "annotation",
            "gene_symbol",
            "is_ins",
        ]
    ]
    
    return df


def annotate_indels(row, genome, exons, chr_prefixed, postprocess=False):
    """Annotates indels for all RefSeq isoforms

    Args: 
        row (pandas.Series): a Series with indices 
                             'chr', 'pos', 'is_ins', 'indel_seq'
        genome (pysam.FastaFile): reference genome
        exons (pysam.TabixFile): coding exon data
        chr_prefixed (bool): True if chromosome names in BAM are "chr"-prefixed
        postprocess (bool): True if used in indel_postprocessor. Default to False 
    
    Returns:
        annotation (str): Each token represents an annotation for one 
                          of the isoforms and is formatted as:
           
            GeneSymbol|RefSeqAccession|AminoAcidPostion|Effect|IsInsensitive
       
            GeneSymbol: RefSeq gene name
            RefSeqAccession: RefSeq mRNA accession number
            CodonPostion: the position of codon (not amino acid) affected in 
                              the isoform specified in RefSeqAccession
            Effect: consequences of the indel.
                    See CodingSequenceWithIndel for detail
            IsInsensitive: 1 if the indel is nonsense-mediated-decay insensitive, 
                             0 otherwise

            '-' for non-coding indels
    """
    chr = row["chr"]
    pos = row["pos"]
    idl_type = row["is_ins"]
    idl_seq = row["indel_seq"]

    # generates CodingSequenceWithIndel instances
    idls = generate_coding_indels(
        chr, pos, idl_type, idl_seq, genome, exons, chr_prefixed
    )

    # annotates for all RefSeq isoforms
    annots = []
    if idls:
        for idl in idls:
            gene = idl.gene_symbol
            refseq_acc = idl.accession
            codon_pos, effect = idl.effect()
            is_insensitive = idl.is_nmd_insensitive()

            if not postprocess:
                anno = (
                    gene
                    + "|"
                    + refseq_acc
                    + "|"
                    + str(codon_pos)
                    + "|"
                    + effect
                    + "|"
                    + str(is_insensitive)
                )
            else:
                anno = gene + "|" + refseq_acc + "|" + str(codon_pos) + "|" + effect

            annots.append(anno)

    annotation = "-" if len(annots) == 0 else ",".join(annots)

    return annotation


def generate_coding_indels(chr, pos, idl_type, idl_seq, genome, exons, chr_prefixed):
    """Generate coding indel objects
    
    Args:
        chr (str): chr1-22, chrX or chrY. Note "chr"-prefixed.
        pos (int): 1-based genomic position
        idl_type (int): 1 for insertion, 0 for deletion
        idl_seq (str): inserted or deleted sequence
        genome (pysam.FastaFile): reference genome
        exons (pysam.TabixFile): coding exon data
        chr_prefixed (bool): True if chromosome names in BAM are "chr"-prefixed

    Returns:
        coding_idl_lst (list): a list of CodingSequenceWithIndel obj
                               empty list if non-coding indel  
    """
    coding_idl_lst = []

    try:
        candidate_genes = exons.fetch(chr, pos - 11, pos + 11)
    except:
        candidate_genes = None

    # check for UTR
    if candidate_genes:
        for line in candidate_genes:
            lst = line.split("\t")

            # parsing exon info
            info = lst[3].split("|")
            exon = int(info[2])
            last_exon = int(info[3])

            # exon start and end
            exon_start, exon_end = int(lst[1]), int(lst[2])

            # strand
            strand = lst[4]

            # 5'UTR on positive strand (insertion)
            if strand == "+" and exon == 1 and idl_type == 1 and exon_start >= pos:
                pass
            # 5'UTR on positive strand (deletion)
            elif strand == "+" and exon == 1 and idl_type == 0 and exon_start > pos:
                pass
            # 3'UTR on positive strand
            elif strand == "+" and exon == last_exon and pos > exon_end:
                pass
            # 5'UTR on negative strand
            elif strand == "-" and exon == 1 and pos > exon_end:
                pass
            # 3'UTR on negative strand (insertion)
            elif (
                strand == "-"
                and exon == last_exon
                and idl_type == 1
                and exon_start >= pos
            ):
                pass
            # 3'UTR on negative strand (deletion)
            elif (
                strand == "-"
                and exon == last_exon
                and idl_type == 0
                and exon_start > pos
            ):
                pass
            else:
                indel_in_reference_genome = curate_indel_in_genome(
                    genome, chr, pos, idl_type, idl_seq, chr_prefixed
                )
                lt_seq = indel_in_reference_genome.lt_seq
                rt_seq = indel_in_reference_genome.rt_seq

                accession = info[0]
                gene_symbol = info[1]
                cds_start = int(info[4])
                prev_exon = lst[5].split("|")
                prev_exon_start, prev_exon_end = int(prev_exon[0]), int(prev_exon[1])
                next_exon = lst[6].split("|")
                next_exon_start, next_exon_end = int(next_exon[0]), int(next_exon[1])

                indel = CodingSequenceWithIndel(
                    chr,
                    pos,
                    idl_type,
                    lt_seq,
                    idl_seq,
                    rt_seq,
                    strand,
                    accession,
                    gene_symbol,
                    exon,
                    exon_start,
                    exon_end,
                    last_exon,
                    cds_start,
                    prev_exon_start,
                    prev_exon_end,
                    next_exon_start,
                    next_exon_end,
                )
                coding_idl_lst.append(indel)

        return coding_idl_lst


def annotate_equivalent_non_coding_indels(row, df, genome, chr_prefixed):
    """Annotate 3'UTR or 3'splice region indels equivalent to 5'coding indels
       
    Args: 
        row (pandas.Series)
        df (pandas.DataFrame): Positionally sorted   
        genome (pysam.FastaFile): reference genome
        chr_prefixed (bool): True if chromosome names in BAM are "chr"-prefixed                        
    Returns:
        annotation (str)
    """
    if row["annotation"] == "-" and row["rescued"] == "by_equivalence":
        query_idl = curate_indel_in_genome(
            genome,
            row["chr"],
            row["pos"],
            row["is_ins"],
            row["indel_seq"],
            chr_prefixed,
        )
        idx = row.name
        found = None
        while not found:
            idx = idx - 1  # look up previous rows
            pos = df.iloc[idx]["pos"]
            idl_seq = df.iloc[idx]["indel_seq"]
            subject_idl = curate_indel_in_genome(
                genome, row["chr"], pos, row["is_ins"], idl_seq, chr_prefixed
            )
            if query_idl == subject_idl and df.iloc[idx]["annotation"] != "-":
                found = True
        return df.iloc[idx]["annotation"]
    else:
        return row["annotation"]


def get_gene_symbol(row):
    """Extracts gene name from annotation

    Args:
        row (pandas.Series): annotation info (str) at 'annotation' index
    Returns:
        gene_symbol (str): gene name(s)
    """
    pd.options.mode.chained_assignment = None

    lst = row["annotation"].split(",")
    genes = [token.split("|")[0] for token in lst]

    gene_symbol = ",".join(set(genes))

    return gene_symbol
