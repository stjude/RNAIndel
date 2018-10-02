#!/usr/bin/env bash
# Test bambino_rna_indel.sh on testdata

BAM="../testdata/inputs/test.bam"
BAMBINO_OUTPUT='./testdata/inputs/bambino_calls.txt'
FASTA="/research/rgs01/project_space/zhanggrp/MethodDevelopment/common/RNAIndel/fasta/GRCh38_no_alt_ERCC92_SpikeIn.fa"
DATA_DIR="/research/projects/zhanggrp/MethodDevelopment/common/RNAIndel/data_dir"
OUTPUT="./testdata/outputs/predicted_indels.vcf"
LOG_DIR="/research/rgs01/home/clusterHome/lding/test/rna_indel"

#rna_indel -b $BAM -i $BAMBINO_OUTPUT -o $OUTPUT -f $FASTA -d $DATA_DIR
./rna_indel.py -b $BAM -i $BAMBINO_OUTPUT -o $OUTPUT -f $FASTA -d $DATA_DIR -l $LOG_DIR
