# RNAIndel 

[RNAIndel](https://doi.org/10.1093/bioinformatics/btz753) calls coding indels from tumor RNA-Seq data and predicts them as somatic, germline, and artifactual. 


## What's new in Version 3:
New implementation with [indelpost](https://github.com/stjude/indelPost), an indel realigner/phaser. 
* faster analysis (typically < 20 min with 8 cores)
* somatic complex indel calling in RNA-Seq
* ensemble calling with your own caller (e.g., GATK HaplotypeCaller/MuTect2)  
* improved sensitivity for homopolymer indels by error-profile outlier analysis  

## Quick Start
### PyPI
RNAIndel depends on [python>=3.6.0](https://www.python.org/downloads/) and [java>=1.8.0](https://www.java.com/en/download/).<br> 
Installing via the pip command will install the following packages:
* [indelpost>=0.0.4](https://github.com/stjude/indelPost)
* [pandas>=0.23.0](https://pandas.pydata.org/)
* [scikit-learn>=0.22.0](http://scikit-learn.org/stable/install.html#)

```
> pip install rnaindel
```

Test the installation.
```
> rnaindel -h
usage: rnaindel <subcommand> [<args>]

subcommands are:
    PredictIndels             Predict somatic/germline/artifact indels from tumor RNA-Seq data
    CalculateFeatures         Calculate and report features for training
    Train                     Perform model training
    CountOccurrence           Count occurrence within cohort to filter false somatic predictions
positional arguments:
  subcommand  PredictIndels, CalculateFeatures, Train, CountOccurrence

optional arguments:
  -h, --help  show this help message and exit
  --version   show program's version number and exit
```

## Usage
RNAIndel has 4 subcommands:
* ```PredictIndels``` analyze RNA-Seq data for indel discovery
* ```CalculateFeatures``` calculate features for training
* ```Train``` train models with user's dataset
* ```CountOccurrence``` annotate over-represented somatic predictions

Subcommands are invoked:
```
> rnaindel subcommand [subcommand-specific options]
```

### Discover somatic indels

#### Input BAM file
RNAIndel expects [STAR](https://academic.oup.com/bioinformatics/article/29/1/15/272537) 2-pass mapped BAM file with sorted by coordinate 
and [MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates). Further preprocessing such as 
indel realignment may prevent desired behavior.

#### Standard calling
This mode uses the built-in caller to analyze simple and complex indels.
```
> rnaindel PredictIndels -i input.bam -o output.vcf -r ref.fa -d data_dir -p 8 (default 1) 
```

#### Ensemble calling 
Indels in the exernal VCF (supplied by -v) are integrated to the callset by the built-in caller to boost performance.<br> 
See [demo](./docs/walkthrough/README.md)
```
> rnaindel PredictIndels -i input.bam -o output.vcf -r ref.fa -d data_dir -v gatk.vcf.gz -p 8
```
#### Options
* ```-i``` input [STAR](https://academic.oup.com/bioinformatics/article/29/1/15/272537)-mapped BAM file (required)
* ```-o``` output VCF file (required)
* ```-r``` reference genome FASTA file (required)
* ```-d``` [data directory](#setup) contains trained models and databases (required)
* ```-v``` VCF file (must be .vcf.gz + index) from user's caller. (default: None)
* ```-p``` number of cores (default: 1)
* <details>
    <summary>other options (click to open)</summary><p>
        
    * ```-q``` STAR mapping quality MAPQ for unique mappers (default: 255)
    * ```-m``` maximum heap space (default: 6000m)
    * ```--region``` target genomic region. specify by chrN:start-stop (default: None)
    * ```--pon``` user's defined list of non-somatic calls such as PanelOfNormals. Supply as .vcf.gz with index (default: None)
    * ```--include-all-external-calls``` set to include all indels in VCF file supplied by -v. (default: False. Use only calls with PASS in FILTER) 
    * ```--skip-homopolyer-outlier-analysis``` no outlier analysis for homopolymer indels (repeat > 4) performed if set. (default: False)  

</p></details>

### Train RNAIndel

### Annotate over-represented putative somatic indels

## Citation
aaa

This version needs [indelPost](https://github.com/stjude/indelPost).
Please install this first.

Also, depends on [scikit-learn>=0.22](http://scikit-learn.org/stable/install.html#).


To install
```
pip install indelpost

git clone https://github.com/rawagiha/RNAIndel.git
cd RNAIndel
python setup.py install
```

Download the datapackage(135MB). The current version requires the new package:
```
curl -LO http://ftp.stjude.org/pub/software/RNAIndel/data_dir_grch38.v3.tar.gz
tar -zxf data_dir_grch38.v3.tar.gz
```

To run somatic analysis:
```
rnaindel PredictIndels -i RNA-Seq BAM -d data_dir -r ref.fa -o output.vcf -p 8 (default 1)
``` 
