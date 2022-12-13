import pysam
from indelpost import Variant, VariantAlignment


def cross_check(df, reference, tdna, ndna):
    reference = pysam.FastaFile(reference)
    if tdna:
        tdna = pysam.AlignmentFile(tdna)

    if ndna:
        ndna = pysam.AlignmentFile(ndna)

    (
        df["t_dna_ref_count"],
        df["t_dna_alt_count"],
        df["n_dna_ref_count"],
        df["n_dna_alt_count"],
    ) = zip(*df.apply(wrap_count, reference=reference, tdna=tdna, ndna=ndna, axis=1))

    return df


def wrap_count(row, reference, tdna, ndna):
    tcnt, ncnt = (-1, -1), (-1, -1)
    variant = Variant(row["chrom"], row["pos"], row["ref"], row["alt"], reference)

    if tdna:
        if row["predicted_class"] == "somatic":
            tcnt = count(variant, tdna)

    if ndna:
        if row["predicted_class"] == "somatic":
            ncnt = count(variant, ndna)

    return tcnt[0], tcnt[1], ncnt[0], ncnt[1]


def count(variant, bam):
    try:
        valn = VariantAlignment(variant, bam)
        return valn.count_alleles(by_fragment=True)
    except:
        return (-1, -1)
