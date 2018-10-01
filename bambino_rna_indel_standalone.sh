#!/usr/bin/env bash

# set -x

print_usage() {
    >&2 echo "Run Bambino and RNAIndel without the need of an installation."
	>&2 echo "Usage: bambino_rna_indel_standalone.sh bam fasta data_dir output_vcf"
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
    >&2 echo "Error: db_snp file not found in data directory."
    exit 1
fi

bambino_output=$(mktemp)

# Export Bambino JAR files to CLASSPATH
RNAIndel_PACKAGE_HOME=$(dirname $(readlink -f $0))
export CLASSPATH=$RNAIndel_PACKAGE_HOME/Bambino/*:$CLASSPATH

# Run unpaired Bambino
java -Xmx6000m Ace2.SAMStreamingSNPFinder -of $bambino_output -fasta $fasta -min-mapq 1 -optional-tags XT!=R \
-bam $bam -tn N -dbsnp-file $db_snp -min-quality 20 -min-flanking-quality 20 -min-alt-allele-count 3 \
-min-minor-frequency 0 -broad-min-quality 10 -mmf-max-hq-mismatches 6 -mmf-max-hq-mismatches-xt-u 10 \
-mmf-min-quality 15 -mmf-max-any-mismatches 6 -unique-filter-coverage 2 -no-strand-skew-filter -illumina-q2 1

# Export root directory of RNAIndel to PATH
export PATH=$RNAIndel_PACKAGE_HOME/RNAIndel:$PATH

# Setting RNAIndel_HOME and cd into it, which is NOT the best practice
RNAIndel_HOME="$RNAIndel_PACKAGE_HOME/RNAIndel"
cd $RNAIndel_HOME

# Run RNAIndel
rna_indel.py -b $bam -i $bambino_output -o $output_vcf -f $fasta -d $data_dir

>&2 echo "Completed successfully"
