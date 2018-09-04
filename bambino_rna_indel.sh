#!/usr/bin/env bash

set -x

module load java/1.8.0_66

RNAIndel_HOME=$(dirname $(readlink -f $0))
export CLASSPATH=$RNAIndel_HOME/Bambino/*:$CLASSPATH

input_bam=$1
output_high_20_file=$2
fasta="/research/rgs01/resgen/prod/tartan/runs/tartan_import/FASTA-sQ8GE1r7/output/FASTA/GRCh38_no_alt.fa"
db_snp="/research/rgs01/resgen/prod/tartan/runs/tartan_import/snp142_binary.blob-15yHdPKI/output/snp142_binary.blo"

# Run unpaired Bambino
java -Xmx6000m Ace2.SAMStreamingSNPFinder -of $output_high_20_file -fasta $fasta -min-mapq 1 -optional-tags XT!=R -bam $input_bam -tn N -dbsnp-file $db_snp -min-quality 20 -min-flanking-quality 20 -min-alt-allele-count 3 -min-minor-frequency 0 -broad-min-quality 10 -mmf-max-hq-mismatches 6 -mmf-max-hq-mismatches-xt-u 10 -mmf-min-quality 15 -mmf-max-any-mismatches 6 -unique-filter-coverage 2 -no-strand-skew-filter -illumina-q2 1

export PATH=$RNAIndel_HOME/RNAIndel:$PATH

# Run RNAIndel
rna_indel.py
