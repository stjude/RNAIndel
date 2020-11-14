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
   RNAIndel calls small coding indels from tumor RNA-Seq data and classifies them as somatic, germline, and artifactual.  You can also use RNAIndel as a postprocessor to classify indels called by your own caller. RNAIndel supports GRCh38 and 37. <br> 
   <br />
   <a href="https://stjudecloud.github.io/RNAIndel/"><strong>Explore the docs »</strong></a>
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

---
## Quick Start
RNAIndel can be executed via Docker or ran locally, downloadable via PyPI.

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
Installing RNAIndel via the pip command will install the dependencies except for Java.  

#### Dependencies
* [python>=3.6.0](https://www.python.org/downloads/)
    * [pandas>=0.23.0](https://pandas.pydata.org/) 
    * [scikit-learn>=0.20.0](http://scikit-learn.org/stable/install.html#)
    * [pysam>=0.13.0](https://pysam.readthedocs.io/en/latest/index.html)
* [java>=1.8.0](https://www.java.com/en/download/) 

```
> pip install rnaindel
```

Test the installation.
```
> rnaindel
usage: rnaindel <subcommand> [<args>]

subcommands are:
    analysis              Predict somatic indels from tumor RNA-Seq data
    feature               Calculate and report features for training
    nonsomatic            Compile non-somatic indel panel
    reclassification      Reclassify false positives by non-somatic panel
    recurrence            Annotate false positives by recurrence
    training              Train models

positional arguments:
  subcommand  analysis, feature, nonsomatic, reclassification, recurrence,
              training

optional arguments:
  -h, --help  show this help message and exit
  --version   show program's version number and exit
```
You can download data package (GRCh38, GRCh37) and unpack it in a convenient directory on your system. 
```
# GRCh38
> curl -LO http://ftp.stjude.org/pub/software/RNAIndel/data_dir_38.tar.gz
> tar -xzf data_dir_38.tar.gz
# GRCh37
> curl -LO http://ftp.stjude.org/pub/software/RNAIndel/data_dir_37.tar.gz
> tar -xzf data_dir_37.tar.gz
```


Test it out!
```
❯ rnaindel analysis -i ./sample_data/inputs/sample.bam \
                      -o output.vcf \
                      -r ./sample_data/inputs/chr10.fa \
                      -d ./data_dir_38
indel calling completed successfully.
rnaindel analysis completed successfully.
```
## Usage 
RNAIndel has 6 subcommands:   
* ```analysis``` analyze RNA-Seq data for indel discovery   
* ```feature``` calculate features for training   
* ```nonsomatic``` make a non-somatic indel panel    
* ```reclassification``` reclassify false positives by non-somatic panel    
* ```recurrence``` annotate false positives by recurrence   
* ```training``` train and update the models   

Subcommands are invoked:
```
rnaindel subcommand [subcommand-specific options]
```

### Discover somatic indels ([demo](./sample_data))

#### Input BAM file
RNAIndel expects [STAR](https://academic.oup.com/bioinformatics/article/29/1/15/272537) 2-pass mapped BAM file with sorted by coordinate 
and [MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates). Further preprocessing such as 
indel realignment may prevent desired behavior (RNAIndel internally realigns indels to correct allele count).  

#### Use the built-in caller
RNAIndel calls indels by the [built-in caller](https://academic.oup.com/bioinformatics/article/27/6/865/236751), which is optimized 
for RNA-Seq indel calling, and classifies detected indels as somatic, germline, and artifactual. 
```
rnaindel analysis -i BAM -o OUTPUT_VCF -r REFERENCE -d DATA_DIR [other options]
```
#### Use your caller 
RNAIndel can be used as a postprocessor for indel calls generated by your caller such as 
[GATK-HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.8.0/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php), 
[Mutect2](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.8.0/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php)
and [freebayes](https://github.com/ekg/freebayes). Specify the input VCF file with ```-v```.
```
rnaindel analysis -i BAM -v INPUT_VCF -o OUTPUT_VCF -r REFERENCE -d DATA_DIR [other options]
```
#### Options
* ```-i``` input [STAR](https://academic.oup.com/bioinformatics/article/29/1/15/272537)-mapped BAM file (required)
* ```-o``` output VCF file (required)
* ```-r``` reference genome FASTA file (required)
* ```-d``` [data directory](#setup) contains trained models and databases (required)
* ```-v``` VCF file from user's caller (default: None)
* <details>
    <summary>other options (click to open)</summary><p>
    
    * ```-q``` STAR mapping quality MAPQ for unique mappers (default: 255)
    * ```-p``` number of cores (default: 1)
    * ```-m``` maximum heap space (default: 6000m)
    * ```-l``` directory to store log files (default: current)
    * ```-n``` user-defined panel of non-somatic indels in tabixed VCF format (default: built-in reviewed indel set)
    * ```-g``` user-provided germline indel database in tabixed VCF format (default: built-in database in data dir) <br>
    &nbsp;   &nbsp;   &nbsp;   &nbsp;use only if the model is trained with the user-provided database ([more](./docs/training)).      
    * ```--region``` target genomic region. specify by chrN:start-stop (default: None)
    * ```--exclude-softclipped-alignments``` softclipped indels will not be used for analysis if added (default: False)

</p></details>

### Train RNAIndel
Users can [train](./docs/training) RNAIndel with their own training set. 

### Filter false positives
RNAIndel supports [custom filtering](./docs/filtering) to refine the predicted results.

## Contact
* kohei.hagiwara[AT]stjude.org   
Please let me know what your experience with RNAIndel was like (even bad comments are welcome)!

## Limitations
1. RNAIndel does not perform well for samples with microsatellite instability such as colon adenocarcinoma hypermutated subtype. 
