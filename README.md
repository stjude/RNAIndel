# RNAIndel
Somatic indel detector for tumor RNA-Seq data.

RNAIndel only works with hg38 (GRCh38).

## References
1. Edmonson, Michael N. et al. “Bambino: A Variant Detector and Alignment Viewer for next-Generation Sequencing Data in 
the SAM/BAM Format.” Bioinformatics 27.6 (2011): 865–866. PMC. Web. 27 Sept. 2018.
2. RNAIndel

## Prerequisites
* python>=3.5.2
    * pandas>=0.22.0
    * numpy>=1.12.0
    * sklearn>=0.18.1
    * pysam=0.15.1
    * pyvcf=0.6.8
* java=1.8.0_66 (for [Bambino](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3051333/) only)

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
python setup.py install     # Install `bambino` and `rna_indel` from source
bambino -h                  # Check if `bambino` works correctly
rna_indel -h                # Check if `rna_indel` works correctly
```

## Run on the command line
Call indels using Bambino:
```
bambino -i BAM -f REF_FASTA -o BAMBINO_RAW_INDEL
```
### bambino Options
* ```-b``` input tumor bam file
* ```-f``` reference genome FASTA file
* ```-o``` output file with Bambino indel calls

Use Bambino indel calls as an input to RNAIndel (highly recommended):
```
rna_indel -b BAM -i BAMBINO_INDEL_CALL -f REF_FASTA -o OUTPUT_VCF [other options]
```

Use indels from other indel caller (e.g., GATK) as an input to RNAIndel:
```
rna_indel -b BAM -c RAW_INDEL_VCF -r REF_FASTA -o OUTPUT_VCF [other options]
```
### rna_indel Options
* ```-b``` input tumor bam file (required)
* ```-i``` input file with Bambino indel calls (required for using Bambino indel calls)
* ```-c``` input vcf file with indel calls (required for using other indel callers, e.g. [GATK](https://software.broadinstitute.org/gatk/))
* ```-o``` output vcf file (required)
* ```-f``` reference genome (GRCh38) FASTA file (required)
* ```-q``` STAR mapping quality MAPQ for unique mappers (default=255)
* ```-p``` number of cores (default=1)
* ```-n``` user-defined panel of non-somatic indel list in vcf format
* ```-r``` [refgene](https://www.ncbi.nlm.nih.gov/refseq/) coding exon database
* ```-d``` indels on [dbSNP database](https://www.ncbi.nlm.nih.gov/snp) in vcf format
* ```-l``` [ClinVar database](https://www.ncbi.nlm.nih.gov/clinvar/)
* ```-m``` directory with trained random forest models

# Run as a Workflow
Run Bambino and RNAIndel as a workflow
