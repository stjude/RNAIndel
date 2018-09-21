# RNAIndel
Somatic indel detector for tumor RNA-Seq data.

# Prerequisites
* python>=3.5.2
    * pandas>=0.22.0
    * numpy>=1.12.0
    * sklearn>=0.18.1
    * pysam=0.15.1
    * pyvcf=0.6.8
* java=1.8.0_66 (for [Bambino](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3051333/) only)

# Download
```
git clone https://github.com/adamdingliang/RNAIndel/tree/master  # Clone the repo
cd RNAIndel                                                      # Switch to source directory
```

# Installation
It is highly recommended to setup a virtual python environment using [conda](https://conda.io/docs/) and install 
the python dependencies in the virtual environment:
```
conda create -n py36 python=3.6 anaconda    # Create a python3.6 virtual environment
source activate py36                        # Activate the virtual environment
pip install -r requirements.txt             # Install python dependencies
```

You can install RNAIndel from source directly:
```
python setup.py install     # Install `bambino` and `rna_indel` from source
bambino -h                  # Check if `bambino` works correctly
rna_indel -h                # Check if `rna_indel` works correctly
```

# Run on the command line
Call indels using Bambino:
```
bambino -i BAM -f REF_FASTA -o BAMBINO_RAW_INDEL
```

Use Bambino indel calls as an input:
```
rna_indel -b BAM -ib BAMBINO_RAW_INDEL -r REF_FASTA -dbsnp DBSNP [other options]
```
Use indels from other indel caller (e.g., GATK) as an input:
```
rna_indel -b BAM -iv RAW_INDEL_VCF -r REF_FASTA -dbsnp DBSNP [other options]
```

# Run as a Workflow
Run Bambino and RNAIndel as a workflow
