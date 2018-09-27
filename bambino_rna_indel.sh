#!/usr/bin/env bash
# This wrapper runs Bambino and RNAIndel without the need of the installation.

# set -x

print_usage() {
    >&2 echo "Runs Bambino first and supplies the indel calls to RNAIndel, then runs RNAIndel."
	>&2 echo "Usage: bambino_rna_indel.sh tumor_bam fasta db_snp"
	>&2 echo "Required parameters:"
	>&2 echo "\$1 tumor_bam: input tumor bam file"
	>&2 echo "\$2 fasta: reference genome (GRCh38) FASTA file"
	>&2 echo "\$3 db_snp: indels on dbSNP database in vcf format"
}

if [ "$#" == 0 ]; then print_usage ; exit 0; fi

# Initialize parameters
tumor_bam=$1
fasta=$2
db_snp=$3

if [ ! -f $tumor_bam ]; then
    >&2 echo "Error: tumor bam file not found."
    exit 1
fi

if [ ! -f $fasta ]; then
    >&2 echo "Error: fasta file not found."
    exit 1
fi

if [ ! -f $db_snp ]; then
    >&2 echo "Error: db_snp file not found."
    exit 1
fi

high_20_file=$(mktemp)

# Export Bambino JAR files to CLASSPATH
RNAIndel_PACKAGE_HOME=$(dirname $(readlink -f $0))
export CLASSPATH=$RNAIndel_PACKAGE_HOME/Bambino/*:$CLASSPATH

# Run unpaired Bambino
java -Xmx6000m Ace2.SAMStreamingSNPFinder -of $high_20_file -fasta $fasta -min-mapq 1 -optional-tags XT!=R \
-bam $tumor_bam -tn N -dbsnp-file $db_snp -min-quality 20 -min-flanking-quality 20 -min-alt-allele-count 3 \
-min-minor-frequency 0 -broad-min-quality 10 -mmf-max-hq-mismatches 6 -mmf-max-hq-mismatches-xt-u 10 \
-mmf-min-quality 15 -mmf-max-any-mismatches 6 -unique-filter-coverage 2 -no-strand-skew-filter -illumina-q2 1

# Export root directory of RNAIndel to PATH
export PATH=$RNAIndel_PACKAGE_HOME/RNAIndel:$PATH

# Setting RNAIndel_HOME and cd into it, which is NOT the best practice
RNAIndel_HOME="$RNAIndel_PACKAGE_HOME/RNAIndel"
cd $RNAIndel_HOME

# Run RNAIndel
rna_indel.py -b $tumor_bam -f $high_20_file -r $fasta -dbsnp $db_snp
