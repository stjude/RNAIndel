class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://sevenbridges.com'
baseCommand: []
inputs:
  - id: subcommand
    type:
      type: enum
      symbols:
        - analysis
        - feature
        - nonsomatic
        - reclassification
        - recurrence
        - training
      name: subcommand
    inputBinding:
      position: 0
      shellQuote: false
  - id: input_bam
    type: File
    inputBinding:
      position: 1
      prefix: '-i'
      shellQuote: false
    label: Input BAM file
    'sbg:fileTypes': BAM
    secondaryFiles:
      - .bai
  - id: fasta
    type: File
    inputBinding:
      position: 3
      prefix: '-r'
      shellQuote: false
    label: Input FASTA
    'sbg:fileTypes': FA
    secondaryFiles:
      - .fai
  - id: data
    type: Directory
    inputBinding:
      position: 4
      prefix: '-d'
      shellQuote: false
    label: Data Directory
  - id: user_caller_vcf
    type: File?
    inputBinding:
      position: 10
      prefix: '-v'
      shellQuote: false
    'sbg:fileTypes': VCF
  - id: star_mapping
    type: int?
    inputBinding:
      position: 10
      prefix: '-q'
      shellQuote: false
  - id: num_cores
    type: int?
    inputBinding:
      position: 10
      prefix: '-p'
      shellQuote: false
  - id: max_heap_space
    type: string?
    inputBinding:
      position: 10
      prefix: '-m'
      shellQuote: false
  - id: log_directory
    type: string?
    inputBinding:
      position: 10
      prefix: '-l'
      shellQuote: false
  - id: nonsomatic_indels_vcf
    type: File?
    inputBinding:
      position: 10
      prefix: '-n'
      shellQuote: false
    'sbg:fileTypes': VCF
  - id: germline_indels_vcf
    type: File?
    inputBinding:
      position: 10
      prefix: '-g'
      shellQuote: false
    'sbg:fileTypes': VCF
  - id: input
    type: string?
    inputBinding:
      position: 10
      prefix: '--region'
      shellQuote: false
outputs:
  - id: '#output_file'
    type: File
    outputBinding:
      glob: '*.vcf'
label: rnaindel
arguments:
  - position: 101
    prefix: '-o'
    shellQuote: false
    valueFrom: $(inputs.input_bam.nameroot).vcf
requirements:
  - class: ShellCommandRequirement
  - class: EnvVarRequirement
    envDef:
      TMPDIR: /data
  - class: InlineJavascriptRequirement
hints:
  - class: DockerRequirement
    dockerPull: 'ghcr.io/stjude/rnaindel:latest'
'sbg:projectName': rnaindel
