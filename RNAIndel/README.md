# Indel classification executable 
The classification executable (rna_indel) classifies indel calls from the built-in Bambino caller or other callers.   

## Prerequisites
* [python>=3.5.2](https://www.python.org/downloads/)
    * [pandas>=0.23.0](https://pandas.pydata.org/)
    * [numpy>=1.12.0](https://www.scipy.org/scipylib/download.html)
    * [scikit-learn=0.18.1](http://scikit-learn.org/stable/install.html#)
    * [pysam=0.15.1](https://pysam.readthedocs.io/en/latest/index.html)
    * [pyvcf=0.6.8](https://pyvcf.readthedocs.io/en/latest/index.html)

## Classify Bambino calls
Specify the input Bambino calls by -i.
```
rna_indel -b BAM \
          -i BAMBINO_CALLS \
          -o OUTPUT_VCF \
          -f FASTA \
          -d DATA_DIR \
          [optional arguments]
```

## Classify calls from other callers
Specify the input VCF by -c.  
```
rna_indel -b BAM \
          -c INPUT_VCF \
          -o OUTPUT_VCF \
          -f FASTA \
          -d DATA_DIR \ 
          [optional arguments]
```

### RNAIndel options
* ```-b``` input BAM file (required)
* ```-i``` [Bambino output file](../Bambino) (required for using Bambino as the indel caller)
* ```-c``` VCF file from other caller (required for using other callers, e.g., [GATK](https://software.broadinstitute.org/gatk/))
* ```-o``` output VCF file (required)
* ```-f``` reference genome (GRCh38) FASTA file (required)
* ```-d``` path to data directory contains trained models and databases. See [Data directory set up](../README.md#data-directory-set-up). 
* ```-q``` STAR mapping quality MAPQ for unique mappers (default=255)
* ```-p``` number of cores (default=1)
* ```-n``` user-defined panel of non-somatic indels in VCF format
* ```-l``` direcotry to store log files
* ```-h``` print help
* ```--version``` print RNAIndel version

