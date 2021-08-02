## Walking through with sample dataset
Using reads around the *PTEN* tumor suppresor gene locus (chr10:87,863,113-87,971,930 on GRCh38), the sample dataset was prepared by [GATK Best Practice for RNA-Seq](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-). The dataset contains the BAM and VCF files (not FASTQ)<br>


<p align="center">
  <img src="files.JPG"/>
</p>


To download:
```
curl -LO http://ftp.stjude.org/pub/software/RNAIndel/sampledataset.v3.tar.gz
tar -zxf sampledataset.v3.tar.gz
```

RNAIndel integrates the callset from the built-in caller and that from GAKT (pre-called).<br>
To do this, input the original BAM file (i.e., sample.bam, not the pre-processed one) and the GATK VCF file:
```
rnaindel PredictIndels -i ./sampledataset/sample.bam \ # do not use sample.gatk.bam
                       -o test.vcf \
                       -r ./sampledataset/chr10.fa \
                       -d ./data_dir_grch38 \
                       -v ./sampledataset/sample.gatk.vcf.gz \
```

The output VCF reports the expressed coding indels called from the 2 methods. Each indel is annoated in INFO field for prdicted class (somatic, germline, artifact) and by which caller it was detected (built-in, external, both). In this example, indels predicted as somatic are:
```
chr10   87957917    AC   GGCCCATGG  CALLER=both
chr10   87957955    C    CCTGGGTT   CALLER=both
```
These *PTEN* indels are called by the built-in and GATK HaplotypeCaller. The first indel is reported as a complex indel AC>GGCCCATGG. 

