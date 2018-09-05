#!/usr/bin/env bash
# Test bambino_rna_indel.sh on testdata

module load java/1.8.0_66

# User should provide a tumor bam
TUMOR_BAM="/research/rgs01/home/clusterHome/lding/Git/RNAIndel/testdata/inputs/test.bam"

# User should provide a fasta file. Note the same fasta file should be used in mapping, variance calling and indel prediction.
#FASTA="/research/rgs01/resgen/prod/tartan/runs/tartan_import/FASTA-sQ8GE1r7/output/FASTA/GRCh38_no_alt.fa"
FASTA="/research/rgs01/project_space/zhanggrp/MethodDevelopment/common/rna_seq_indels/fasta/GRCh38_no_alt_ERCC92_SpikeIn.fa"

# We should provide this file and its index file for user to download
DB_SNP="/research/rgs01/project_space/zhanggrp/MethodDevelopment/common/rna_seq_indels/dbsnp/00-All.151.indel.vcf.gz"

./bambino_rna_indel.sh $TUMOR_BAM $FASTA $DB_SNP