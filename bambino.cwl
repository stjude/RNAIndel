#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: bambino

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

outputs:
  out_calls:
    type: File
    outputBinding:
      glob: $(inputs.bambino_output)

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
    Bambino is a Variant Detector for next-Generation Sequencing Data.
    This CWL script is for a python wrapper of Bambino for RNAIndel usage.

$namespaces:
  s: http://schema.org/

$schemas:
 - https://schema.org/docs/schema_org_rdfa.html