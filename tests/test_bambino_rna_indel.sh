#!/usr/bin/env bash
# Test bambino_rna_indel.sh on testdata

module load java/1.8.0_66

# User should provide a tumor bam
BAM="/research/rgs01/home/clusterHome/lding/Git/RNAIndel/testdata/inputs/test.bam"

# User should provide a fasta file.
# Note the same fasta file should be used in mapping, variance calling and indel prediction.
#FASTA="/research/rgs01/resgen/prod/tartan/runs/tartan_import/FASTA-sQ8GE1r7/output/FASTA/GRCh38_no_alt.fa"
FASTA="/research/rgs01/project_space/zhanggrp/MethodDevelopment/common/RNAIndel/fasta/GRCh38_no_alt_ERCC92_SpikeIn.fa"

DATA_DIR="/research/projects/zhanggrp/MethodDevelopment/common/RNAIndel/data_dir"
OUTPUT_VCF="./testdata/outputs/predicted_indels.vcf"


#./bambino_rna_indel.sh $BAM $FASTA $DATA_DIR $OUTPUT_VCF
./bambino_rna_indel_standalone.sh $BAM $FASTA $DATA_DIR $OUTPUT_VCF