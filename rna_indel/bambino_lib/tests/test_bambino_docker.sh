#!/usr/bin/env bash

docker run --mount type=bind,source=/Users/lding/Git/RNAIndel/testdata/inputs,target=/inputs,readonly \
--mount type=bind,source=/Users/lding/Git/RNAIndel/fasta,target=/refs,readonly \
--mount type=bind,source=/Users/lding/Git/RNAIndel/Bambino/testdata/outputs,target=/outputs \
-i rna_indel:0.1.0 \
bambino -b /inputs/test.bam -f /refs/GRCh38_no_alt_ERCC92_SpikeIn.fa \
-d /RNAIndel/data_dir/dbsnp/00-All.151.indel.vcf.gz \
-o /outputs/bambino_calls_docker.txt
