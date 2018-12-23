#!/usr/bin/env bash
# Test bambino executable on testdata

module load java/1.8.0_66

# User should provide a tumor bam
TUMOR_BAM="./testdata/inputs/test.bam"

# User should provide a fasta file. Note the same fasta file should be used in mapping, variance calling and indel prediction.
#FASTA="/research/rgs01/resgen/prod/tartan/runs/tartan_import/FASTA-sQ8GE1r7/output/FASTA/GRCh38_no_alt.fa"
FASTA="/research/rgs01/project_space/zhanggrp/MethodDevelopment/common/RNAIndel/fasta/GRCh38_no_alt_ERCC92_SpikeIn.fa"

# We should provide this file and its index file for user to download
DB_SNP="/research/rgs01/project_space/zhanggrp/MethodDevelopment/common/RNAIndel/dbsnp/00-All.151.indel.vcf.gz"

BAMBINO_OUTPUT='./testdata/outputs/bambino_calls.txt'

bambino -b $TUMOR_BAM -f $FASTA -d $DB_SNP -o $BAMBINO_OUTPUT
