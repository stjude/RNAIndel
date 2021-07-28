# RNAIndel

[RNAIndel](https://doi.org/10.1093/bioinformatics/btz753) calls coding indels from tumor RNA-Seq data and predicts them as somatic, germline, and artifactual. Key features of version 3 include:

* faster analysis (typically < 20 min with 8 cores)
* somatic complex indel calling in RNA-Seq
* ensemble calling with your own caller (e.g., MuTect2)  
 


## Usage
RNAIndel has 4 subcommands:
* ```PredictIndels``` analyze RNA-Seq data for indel discovery
* ```CalculateFeatures``` calculate features for training
* ```Train``` train models with user's dataset
* ```CountOccurrence``` annotate over-represented somatic predictions

Subcommands are invoked:
```
rnaindel subcommand [subcommand-specific options]
```

### Discover somatic indels

#### Input BAM file
RNAIndel expects [STAR](https://academic.oup.com/bioinformatics/article/29/1/15/272537) 2-pass mapped BAM file with sorted by coordinate 
and [MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates). Further preprocessing such as 
indel realignment may prevent desired behavior (RNAIndel internally realigns indels using [indelPost](https://github.com/stjude/indelPost)).

#### Standard calling
```
rnaindel PredictIndels -i input.bam -o output.vcf -r ref.fa -d data_dir -p 8 (default 1) 
```

#### Ensemble calling 
```
rnaindel PredictIndels -i input.bam -o output.vcf -r ref.fa -d data_dir -v mutect2.vcf.gz -p 8 (default 1)
```


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
