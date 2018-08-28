#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: rnaindel_dev.py

inputs:
  bambinoOutput:
    type: File
    inputBinding:
      prefix: -f

  bamFile:
    type: File
    secondaryFiles: .bai
    inputBinding:
      prefix: -b
 
  numOfProcesses:
    type: int?
    inputBinding:
      prefix: -p

  reclassification:
    type: boolean?
    inputBinding:
      prefix: -re-clf
outputs:
  predictedresult:
    type: File
    outputBinding:
      glob: rna_indels.txt
