#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

inputs:
  bam:
    type: File
    secondaryFiles: .bai
    inputBinding:
      prefix: -b
    doc: input bam file

  reference_fasta:
    type: File
    secondaryFiles: .fai
    inputBinding:
      prefix: -f
    doc: reference genome (GRCh38) FASTA file

  bambino_output:
    type: string
    inputBinding:
      prefix: -o
    doc: Bambino output file

  heap_memory:
    type: string
    default: 6000m
    inputBinding:
      prefix: -m
    doc: maximum heap space

  input_vcf:
    type: File?
    inputBinding:
      prefix: -c
    doc: vcf file with indel calls (required for using other callers, e.g. GATK)

  output_vcf:
    type: string
    inputBinding:
      prefix: -o
    doc: output vcf file name

  data_dir:
    type: Directory
    inputBinding:
      prefix: -d
    doc: data directory contains refgene, dbsnp and clivar databases

  star_mapq:
    type: int
    default: 255
    inputBinding:
      prefix: -q
    doc: STAR mapping quality MAPQ for unique mappers (default=255)

  process_num:
    type: int?
    default: 1
    inputBinding:
      prefix: -p
    doc: number of cores (default=1)

  non-somatic_list:
    type: File?
    inputBinding:
      prefix: -n
    doc: user-defined panel of non-somatic indel list in vcf format

outputs:
  out:
    type: File
    outputSource: rna_indel/out_vcf

steps:
  bambino:
    run: bambino.cwl
    in:
      bam: bam
      reference_fasta: reference_fasta
      bambino_output: bambino_output
      heap_memory: heap_memory
    out: [out_calls]

  rna_indel:
    run: rna_indel.cwl
    in:
      bam: bam
      bambino_output: bambino/out_calls
      data_dir: data_dir
      input_vcf: input_vcf
      output_vcf: output_vcf
      reference_fasta: reference_fasta
    out: [out_vcf]

s:author:
  class: s:Person
  s:name: Liang Ding
  s:email: mailto:Liang.Ding@stjude.org
  s:worksFor:
  - class: s:Organization
    s:name: St.Jude Children's Research Hospital
    s:location: 262 Danny Thomas Place Memphis, TN 38105
    s:department:
    - class: s:Organization
      s:name: Computational Biology

doc: |
    Somatic indel detector for tumor RNA-Seq data

$namespaces:
  s: http://schema.org/

$schemas:
 - https://schema.org/docs/schema_org_rdfa.html