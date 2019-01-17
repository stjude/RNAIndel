#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: rnaindel

hints:
  DockerRequirement:
    dockerPull: adamdingliang/rnaindel:0.3.0

inputs:
  bam:
    type: File
    secondaryFiles: .bai
    inputBinding:
      prefix: -b
    doc: input BAM file

  reference_fasta:
    type: File
    secondaryFiles: .fai
    inputBinding:
      prefix: -f
    doc: reference genome FASTA file

  data_dir:
    type: Directory
    inputBinding:
      prefix: -d
    doc: data directory contains databases and models
  
  output_vcf:
    type: string
    inputBinding:
      prefix: -o
    doc: output vcf file name

  input_vcf:
    type: File?
    inputBinding:
      prefix: -c
    doc: VCF file from other callers

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

  non-somatic_panel:
    type: File?
    inputBinding:
      prefix: -n
    doc: user-defined panel of non-somatic indels in VCF format

  heap_memory:
    type: string?
    inputBinding:
      prefix: -m
    doc: maximum heap space

  log_dir:
    type: Directory?
    inputBinding:
      prefix: -l
    doc: directory to store log files

outputs:
outputs:
  out_vcf:
    type: File
    outputBinding:
      glob: $(inputs.output_vcf)


s:mainEntity:
  class: s:SoftWareSourceCode
  s:name: RNAIndel
  s:url: https://github.com/stjude/RNAIndel
  s:codeRepository: https://github.com/stjude/RNAIndel.git
  s:author:
  - class: s:Person
    s:name: Kohei Hagiwara
    s:email: mailto:Kohei.Hagiwara@stjude.org 

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
    Somatic indel discovery tool for tumor RNA-Seq data

$namespaces:
  s: http://schema.org/

$schemas:
 - https://schema.org/docs/schema_org_rdfa.html
