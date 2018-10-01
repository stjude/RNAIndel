#!/usr/bin/env bash

# set -x

print_usage() {
    >&2 echo "This wrapper runs Bambino and RNAIndel. Installation of bambino and rna_indel executables is required."
	>&2 echo "Usage: bambino_rna_indel.sh bam fasta data_dir output_vcf"
	>&2 echo
	>&2 echo "Required arguments:"
	>&2 echo "\$1 bam: input bam file"
	>&2 echo "\$2 fasta: reference genome (GRCh38) FASTA file"
	>&2 echo "\$3 data_dir: path to RNAIndel data directory"
	>&2 echo "\$4 output_vcf: output vcf file with predicted indels"
}

if [ "$#" == 0 ]; then print_usage ; exit 0; fi


# Initialize parameters
bam=$1
fasta=$2
data_dir=${3%/}
output_vcf=$4

if [ ! -f $bam ]; then
    >&2 echo "Error: bam file not found."
    exit 1
fi

if [ ! -f $fasta ]; then
    >&2 echo "Error: fasta file not found."
    exit 1
fi

if [ ! -d $data_dir ]; then
    >&2 echo "Error: data directory not found."
    exit 1
fi

db_snp="$data_dir/dbsnp/00-All.151.indel.vcf.gz"
if [ ! -f $db_snp ]; then
    >&2 echo "Error: db_snp file $db_snp not found in data directory."
    exit 1
fi

bambino_output=$(mktemp)
bambino -b $bam -f $fasta -d $db_snp -o $bambino_output
rna_indel -b $bam -i $bambino_output -o $output_vcf -f $fasta -d $data_dir
