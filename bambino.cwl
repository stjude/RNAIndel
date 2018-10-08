#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: bambino

inputs:

outputs:


s:author:
  class: s:Person
  s:name: Kohei Hagiwara, Liang Ding
  s:email: mailto:Liang.Ding@stjude.org, kohei.hagiwara@stjude.org
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