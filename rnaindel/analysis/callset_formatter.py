import os
import pysam
import pandas as pd

CANONICALS = [str(i) for i in range(1, 23)] + ["X", "Y"]


def format_callset(tmp_dir, external_vcf, pass_only, region):
    allchroms = get_outfile(tmp_dir)
    if allchroms:
        append_external_indels(
            allchroms, external_vcf, pass_only, which_chrom=0, region=region
        )
        return allchroms

    by_chroms = get_chrom_files(tmp_dir)
    available_chroms = []
    for each_chrom in by_chroms:
        append_external_indels(
            each_chrom[0], external_vcf, pass_only, which_chrom=each_chrom[1]
        )
        available_chroms.append(each_chrom[0])

    return available_chroms


def get_outfile(tmp_dir):
    if os.path.isfile("{}/outfile.txt".format(tmp_dir)):
        return "{}/outfile.txt".format(tmp_dir)
    else:
        return None


def get_chrom_files(tmp_dir):
    return (
        (os.path.join(tmp_dir, "chr{}.txt".format(chrom)), chrom)
        for chrom in CANONICALS
        if os.path.isfile(os.path.join(tmp_dir, "chr{}.txt".format(chrom)))
    )


def append_external_indels(filepath, external_vcf, pass_only, which_chrom, region=None):

    if os.stat(filepath).st_size == 0:
        df1 = supply_empty_df()
    else:
        df1 = pd.read_csv(filepath, sep="\t")
        df1 = df1[(df1["Type"] == "deletion") | (df1["Type"] == "insertion")]
        df1 = df1[["Chr", "Pos", "Chr_Allele", "Alternative_Allele", "Type"]]

    if not external_vcf:
        df1.to_csv(filepath, sep="\t", index=False)
    else:
        external_vcf = pysam.VariantFile(external_vcf)

        if which_chrom and not region:
            if is_chr_prefixed_vcf(external_vcf):
                this_chrom = "chr{}".format(which_chrom)
            else:
                this_chrom = "{}".format(which_chrom)

            records = external_vcf.fetch(contig=this_chrom)
        elif not which_chrom and region:
            chrom, start_stop = (
                region.split(":")[0].replace("chr", ""),
                region.split(":")[1],
            )
            start, stop = int(start_stop.split("-")[0]), int(start_stop.split("-")[1])
            if is_chr_prefixed_vcf(external_vcf):
                this_chrom = "chr{}".format(chrom)
            else:
                this_chrom = "{}".format(chrom)

            records = external_vcf.fetch(
                contig=this_chrom, start=(start - 1), stop=stop
            )
        elif not which_chrom and not region:
            records = external_vcf.fetch()

        data_lst = []
        for record in records:
            update_data(data_lst, record, pass_only)

        if data_lst:
            df2 = pd.DataFrame(data_lst)
            df = pd.concat([df1, df2], ignore_index=True, sort=False)
        else:
            df = df1

        df.to_csv(filepath, sep="\t", index=False)


def supply_empty_df():
    header = ["Chr", "Pos", "Chr_Allele", "Alternative_Allele", "Type"]
    return pd.DataFrame(columns=header)


def is_chr_prefixed_vcf(external_vcf):
    return list(external_vcf.header.contigs)[0].startswith("chr")


def update_data(data_lst, record, pass_only, max_indel_len=50):
 
    if pass_only:
        try:
            if record.filter.__getitem__("PASS").record:
                pass
        except:
            return None
    else:
        pass

    d = {}
    d["Chr"] = "chr" + record.chrom.replace("chr", "")
    d["Pos"] = record.pos
    d["Chr_Allele"] = record.ref
    for alt in record.alts:
        if (
            len(alt) != len(record.ref)
            and max(len(alt), len(record.ref)) < max_indel_len
        ):
            d["Alternative_Allele"] = alt
            d["Type"] = "external"

            data_lst.append(d)
