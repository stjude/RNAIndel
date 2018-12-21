#!/usr/bin/env bash
#cwl-runner --preserve-environment CLASSPATH --outdir ./testdata/outputs ./bambino.cwl ./testdata/inputs/bambino.yml
#cwl-runner --outdir ./testdata/outputs ./rna_indel.cwl ./testdata/inputs/rna_indel.yml
cwl-runner --outdir ./testdata/outputs ./bambino_rna_indel.cwl ./testdata/inputs/bambino_rna_indel.yml