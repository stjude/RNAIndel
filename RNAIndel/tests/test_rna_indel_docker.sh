#!/usr/bin/env bash

docker run \
--mount type=bind,source=/Users/lding/Git/RNAIndel/RNAIndel/testdata/inputs,target=/inputs,readonly \
--mount type=bind,source=/Users/lding/Git/RNAIndel/fasta,target=/refs,readonly \
--mount type=bind,source=/Users/lding/Git/RNAIndel/RNAIndel/testdata/outputs,target=/outputs \
-i rna_indel:0.1.0 \
rna_indel -b /inputs/test.bam -i /inputs/bambino_calls.txt \
-f /refs/GRCh38_no_alt_ERCC92_SpikeIn.fa \
-d /RNAIndel/data_dir \
-o /outputs/predicted_indels_docker.vcf