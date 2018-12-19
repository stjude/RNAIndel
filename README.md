# RNAIndel
RNAIndel calls coding indels and classifies them into 
somatic, germline, and artifact from a tumor RNA-Seq data.
Users can also classify indels called by their own callers by
supplying a VCF file.

## Table of Contents
**[Citations](#citations)**<br>
**[Prerequisites](#prerequisites)**<br>
**[Download](#download)**<br>
**[Installation](#installation)**<br>
**[Input BAM file](#input-bam-file)**<br>
**[Run on the command line](#run-on-the-command-line)**<br>
**[Run Bambino and RNAIndel as a Workflow](#run-bambino-and-rnaindel-as-a-workflow)**<br>
**[Preparation of non-somatic indel panel](#preparation-of-non-somatic-indel-panel)**<br>

## Citations
1. RNAIndel (in preparation).

2. Edmonson, M.N., Zhang, J., Yan, C., Finney, R.P., Meerzaman, D.M., and Buetow, K.H. Bambino: A Variant Detector 
and Alignment Viewer for next-Generation Sequencing Data in 
the SAM/BAM Format. Bioinformatics 27.6 (2011): 865–866. 
DOI: [10.1093/bioinformatics/btr032](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3051333/)

## Prerequisites
* [python>=3.5.2](https://www.python.org/downloads/)
    * [pandas>=0.23.0](https://pandas.pydata.org/)
    * [numpy>=1.12.0](https://www.scipy.org/scipylib/download.html)
    * [scikit-learn=0.18.1](http://scikit-learn.org/stable/install.html#)
    * [pysam=0.15.1](https://pysam.readthedocs.io/en/latest/index.html)
    * [pyvcf=0.6.8](https://pyvcf.readthedocs.io/en/latest/index.html)
* [java=1.8.0_66](https://www.java.com/en/download/) (required for Bambino only)


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

## Data directory set up
Download [data_dir.tar.gz](http://ftp.stjude.org/pub/software/RNAIndel/data_dir.tar.gz) <br>
Place the gziped file under your working directory and unpack it.<br>
```
tar xzvf data_dir.tar.gz
```

## Input BAM file
Currently, RNAIndel only supports STAR-mapped BAM files (GRCh38).<br>
Please prepare your input as follows:<br>

Step 1. Map your reads with the STAR 2-pass mode to GRCh38.<br>
Step 2. Add read groups, sort, mark duplicates, and index the BAM file with Picard.<br>

Please input the BAM file from Step 2 without additional caller specific preprocessing such as indel realignment.<br>
Additional processing steps may prevent desired behavior.

## Run on the command line

### Indel calling using Bambino
A separate Bambino executable is provided with parameters optimized for RNA-Seq variant calling.
The output is a flat tab-delimited file contains SNVs and indels. 
```
bambino -i BAM -f REF_FASTA -o BAMBINO_OUTPUT [optional argument (-m)]
```

#### Bambino options
* ```-m``` maximum heap space (default 6000m)
* ```-b``` input BAM file (required)
* ```-f``` reference genome FASTA file (required)
* ```-o``` Bambino output file (required)

### Run RNAIndel for classification
#### Classification of Bambino calls
```
rna_indel -b BAM -i BAMBINO_OUTPUT -o OUTPUT_VCF -f REF_FASTA -d DATA_DIR [optional arguments]
```
#### Classification of calls from other callers
The input VCF file may contain SNVs.
```
rna_indel -b BAM -c INPUT_VCF -o OUTPUT_VCF -f REF_FASTA -d DATA_DIR [optional arguments]
```

#### RNAIndel options
* ```-b``` input BAM file (required)
* ```-i``` Bambino output file (required for using Bambino as the indel caller)
* ```-c``` VCF file from other caller (required for using other callers, e.g., [GATK](https://software.broadinstitute.org/gatk/))
* ```-o``` output VCF file (required)
* ```-f``` reference genome (GRCh38) FASTA file (required)
* ```-d``` path to data directory contains refgene, dbsnp and clinvar databases [Data directory set up](#data-direcotry-set-up) 
* ```-q``` STAR mapping quality MAPQ for unique mappers (default=255)
* ```-p``` number of cores (default=1)
* ```-n``` user-defined panel of non-somatic indels in VCF format
<!--
* ```-r``` [refgene](https://www.ncbi.nlm.nih.gov/refseq/) coding exon database
* ```-d``` indels on [dbSNP database](https://www.ncbi.nlm.nih.gov/snp) in vcf format
* ```-l``` [ClinVar database](https://www.ncbi.nlm.nih.gov/clinvar/)
* ```-m``` directory with trained random forest models -->


## Run Bambino and RNAIndel as a workflow
### Use CWl scripts (recommended)
To do

### Use BASH wrapper
This requires the [installation](#installation) of `bambino` and `rna_indel` executables.<br>
This pipeline calls indels by Bambino and classifies them.
```
rna_indel_pipeline.sh -b BAM -o OUTPUT_VCF -f REF_FASTA -d DATA_DIR [other options]
```
When a VCF file is supplied by -c,  indel entries in the VCF file are used for classification (variant calling by Bambino will not be performed).
```
rna_indel_piepline.sh -b BAM -c INPUT_VCF -o OUTPUT_VCF -f REF_FASTA -d DATA_DIR [other options]
```
See [Bambino options](#bambino-options) and [RNAIndel options](#rnaindel-options) for the explanations of the options.

## Preparation of non-somatic indel panel
Somatic prediction can be refined by applying a user-defined exclusion panel. RNAIndel reclassifies indels predicted as somatic if <br>
they are found on the panel to either germline or artifact class with a higher probability. This panel is especially useful to <br>
remove recurrent artifacts. Users can compile their panel as follows: <br>
1. Prepare a dataset with RNA-Seq and DNA-Seq performed. <br>
2. Apply RNAIndel and collect indels predicted somatic. <br>
3. Validate the indels with the DNA-Seq data. <br>
4. Collect indels which are validated as germline or artifact in N samples or more. <br>   
5. Format the recurrent non-somatic indels in VCF format.<br>
6. Index the VCF file with Tabix.<br>     
7. Apply the panel to your new sample. <br>
<br>
The panel described in the RNAIndel manuscript can be download here.<br> 
