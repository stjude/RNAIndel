# RNAIndel
Somatic indel detector for tumor RNA-Seq data.

RNAIndel only works with hg38 (GRCh38).

## Citations
1. Edmonson, M.N., Zhang, J., Yan, C., Finney, R.P., Meerzaman, D.M., and Buetow, K.H.. Bambino: A Variant Detector and Alignment Viewer for next-Generation Sequencing Data in 
the SAM/BAM Format. Bioinformatics 27.6 (2011): 865–866. 
DOI: [10.1093/bioinformatics/btr032](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3051333/)
2. RNAIndel. Submitted.

## Prerequisites
* python>=3.5.2
    * pandas>=0.22.0
    * numpy>=1.12.0
    * sklearn>=0.18.1
    * pysam=0.15.1
    * pyvcf=0.6.8
* java=1.8.0_66 (required for Bambino only)

## Download
```
git clone https://github.com/adamdingliang/RNAIndel/tree/master  # Clone the repo
```

## Installation
It is highly recommended to setup a virtual python environment using [conda](https://conda.io/docs/) and install 
the python dependencies in the virtual environment:
```
conda create -n py36 python=3.6 anaconda    # Create a python3.6 virtual environment
source activate py36                        # Activate the virtual environment
pip install -r requirements.txt             # Install python dependencies
```

You can install RNAIndel from source directly:
```
cd RNAIndel                 # Switch to source directory
python setup.py install     # Install bambino and rna_indel from source
bambino -h                  # Check if bambino works correctly
rna_indel -h                # Check if rna_indel works correctly
```

## Run on the command line
### Indel calling using Bambino
```
bambino -i BAM -f REF_FASTA -o BAMBINO_OUTPUT
```
#### Bambino options
* ```-b``` input bam file (required)
* ```-f``` reference genome FASTA file (required)
* ```-d``` indels on [dbSNP database](https://www.ncbi.nlm.nih.gov/snp) (required)
* ```-o``` Bambino output file (required)

### Run RNAIndel with Bambino calls (recommended)
```
rna_indel -b BAM -i BAMBINO_OUTPUT -o OUTPUT_VCF -f REF_FASTA -d DATA_DIR [other options]
```
### Run RNAIndel with indels from other callers
```
rna_indel -b BAM -c INDEL_CALL_VCF -o OUTPUT_VCF -f REF_FASTA -d DATA_DIR [other options]
```
#### RNAIndel options
* ```-b``` input bam file (required)
* ```-i``` Bambino output file (required for using Bambino as the indel caller)
* ```-c``` vcf file with indel calls (required for using other callers, e.g. [GATK](https://software.broadinstitute.org/gatk/))
* ```-o``` output vcf file (required)
* ```-f``` reference genome (GRCh38) FASTA file (required)
* ```-d``` data directory contains refgene, dbsnp and clivar databases
* ```-q``` STAR mapping quality MAPQ for unique mappers (default=255)
* ```-p``` number of cores (default=1)
* ```-n``` user-defined panel of non-somatic indel list in vcf format
<!--
* ```-r``` [refgene](https://www.ncbi.nlm.nih.gov/refseq/) coding exon database
* ```-d``` indels on [dbSNP database](https://www.ncbi.nlm.nih.gov/snp) in vcf format
* ```-l``` [ClinVar database](https://www.ncbi.nlm.nih.gov/clinvar/)
* ```-m``` directory with trained random forest models -->

## Run as a Workflow
Run Bambino and RNAIndel as a workflow
