<p align="center">

  <h1 align="center">
    RNAIndel
  </h1>

  <p align="center">
   <a href="https://github.com/stjude/RNAIndel" target="_blank">
     <img alt="Status"
          src="https://img.shields.io/badge/status-active-success.svg" />
   </a>
   <a href="https://github.com/stjude/RNAIndel/issues" target="_blank">
     <img alt="Github Issues"
          src="https://img.shields.io/github/issues/stjude/RNAIndel" />
   </a>
   <a href="https://github.com/stjude/RNAIndel/pulls" target="_blank">
     <img alt="Pull Requests"
          src="https://img.shields.io/github/issues-pr/stjude/RNAIndel" />
   </a>
   <a href="https://github.com/stjude/RNAIndel/blob/master/LICENSE.md" target="_blank">
     <img alt="License: MIT"
          src="https://img.shields.io/badge/License-MIT-blue.svg" />
   </a>
   <a href="https://badge.fury.io/py/rnaindel" target="_blank">
     <img alt="PyPI version"
          src="https://badge.fury.io/py/rnaindel.png" />
   </a>
   <br />
   <a href="https://github.com/stjude/RNAIndel/actions?query=workflow%3ADocumentation" target="_blank">
     <img alt="Actions: Documentation Status"
          src="https://github.com/stjude/RNAIndel/workflows/Documentation/badge.svg" />
   </a>
   <a href="https://github.com/stjude/RNAIndel/actions?query=workflow%3APackage" target="_blank">
     <img alt="Actions: Package Status"
          src="https://github.com/stjude/RNAIndel/workflows/Package/badge.svg" />
   </a>
  </p>


  <p align="center">
   RNAIndel calls coding indels from tumor RNA-Seq data and classifies them as somatic, germline, and artifactual. RNAIndel supports GRCh38 and 37. <br> 
   <br />
   <a href="https://stjude.github.io/RNAIndel/"><strong>Explore the docs »</strong></a>
   <br />
   <a href="https://doi.org/10.1093/bioinformatics/btz753"><strong>Read the paper »</strong></a>
   <br />
   <br />
   <a href="https://github.com/stjude/RNAIndel/issues/new?assignees=&labels=&template=feature_request.md&title=Descriptive%20Title&labels=enhancement">Request Feature</a>
    | 
   <a href="https://github.com/stjude/RNAIndel/issues/new?assignees=&labels=&template=bug_report.md&title=Descriptive%20Title&labels=bug">Report Bug</a>
   <br />
    ⭐ Consider starring the repo! ⭐
   <br />
  </p>
</p>

## What's new in Version 3
New implementation with [indelpost](https://github.com/stjude/indelPost), an indel realigner/phaser. 
* [faster analysis](#benchmarking) (typically < 20 min with 8 cores)
* somatic complex indel calling in RNA-Seq
* ensemble calling with your own caller (e.g., GATK HaplotypeCaller/MuTect2)  
* improved sensitivity for homopolymer indels by error-profile outlier analysis  

## Quick Start
RNAIndel can be executed via Docker or run locally, downloadable via PyPI.

### Docker
We publish our latest docker builds on GitHub.  You can run the latest code base by running the following command
```
> docker run --rm -v ${pwd}:/data ghcr.io/stjude/rnaindel:latest
```

If you want to have a more native feel, you can add an alias to your shell's rc file.
```
> alias rnaindel="docker run --rm -v ${pwd}:/data ghcr.io/stjude/rnaindel:latest"
```
Note: if its the first time you are executing the `docker run` command, you will see the output of docker downloading the image

### PyPI
RNAIndel depends on [python>=3.6.0](https://www.python.org/downloads/) and [java>=1.8.0](https://www.java.com/en/download/).<br> 
Installing via the pip command will install the following packages:
* [indelpost>=0.0.4](https://github.com/stjude/indelPost)
* [pysam>=0.15.0](https://github.com/pysam-developers)
* [cython>=0.29.12](https://cython.org/)
* [numpy>=1.16.0](https://numpy.org/)
* [ssw-py==0.2.6](https://github.com/Wyss/ssw-py)
* [pandas>=0.23.0](https://pandas.pydata.org/)
* [scikit-learn>=0.22.0](http://scikit-learn.org/stable/install.html#)

```
> pip install indelpost --no-binary indelpost --no-build-isolation  
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

### DataPackage
Download data package (version 3 is not compatible with the previous data package). 
```
#GRCh38
curl -LO http://ftp.stjude.org/pub/software/RNAIndel/data_dir_grch38.v3.tar.gz
tar -zxf data_dir_grch38.v3.tar.gz

#GRCh37
curl -LO http://ftp.stjude.org/pub/software/RNAIndel/data_dir_grch37.v3.tar.gz
tar -zxf data_dir_grch37.v3.tar.gz
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
See [demo](./docs/walkthrough/README.md).
```
> rnaindel PredictIndels -i input.bam -o output.vcf -r ref.fa -d data_dir -v gatk.vcf.gz -p 8
```
#### Options
* ```-i``` input [STAR](https://academic.oup.com/bioinformatics/article/29/1/15/272537)-mapped BAM file (required)
* ```-o``` output VCF file (required)
* ```-r``` reference genome FASTA file (required)
* ```-d``` [data directory](#datapackage) contains trained models and databases (required)
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

#### Benchmarking
Using pediatric tumor RNA-Seq samples ([SJC-DS-1003](https://platform.stjude.cloud/data/cohorts#), n=77), 
the time and memory consumption was benchmarked for ensemble calling with 8 cores (i.e., -p 8) 
on a server with 32-core AMD EPYC 7542 CPU @2.90 GHz.

|       | Run time (wall) | Max memory | 
|------ | -------------   | ---------- |     
|median | 374 sec         | 18.6 GB    |
|max    | 1388 sec        | 23.5 GB    |

### Train RNAIndel
Users can [train](./docs/training) RNAIndel with their own training set. 

### Annotate over-represented putative somatic indels
Check [occurrence](./docs/filtering) to filter probable false positives.

## Contact
* kohei.hagiwara[AT]stjude.org   
Please let me know what your experience with RNAIndel was like (even bad comments are welcome)!

## Citation
Published in [Bioinformatics](https://doi.org/10.1093/bioinformatics/btz753)
