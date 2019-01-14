#!/usr/bin/env bash
#cwl-runner --preserve-environment CLASSPATH --outdir ./sample_data/outputs ./bambino.cwl ./sample_data/inputs/bambino.yml
#cwl-runner --outdir ./sample_data/outputs ./rna_indel.cwl ./sample_data/inputs/rna_indel.yml
cwl-runner --outdir ./sample_data/outputs ./bambino_rna_indel.cwl ./sample_data/inputs/bambino_rna_indel.yml
