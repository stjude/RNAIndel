## Walking through with sample dataset
The sample RNA-Seq dataset was prepared using reads around the *PTEN* tumor suppresor gene locus (chr10:87,863,113-87,971,930 on GRCh38). As shown below, variants were called by [GATK Best Practice for RNA-Seq](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-).<br>






To download:
```
curl -LO http://ftp.stjude.org/pub/software/RNAIndel/sampledataset.v3.tar.gz
tar -zxf sampledataset.v3.tar.gz
```

RNAIndel integrates the callset from the built-in caller and that from GAKT (pre-called):
```
rnaindel PredictIndels -i ./sampledataset/sample.bam \ # do not use sample.gatk.bam
                       -o test.vcf \
                       -r ./sampledataset/chr10.fa \
                       -d ./data_dir_grch38 \
                       -v ./sampledataset/sample.gatk.vcf.gz \
```

The output VCF reports the expressed coding indels called from the 2 callers. Each indel is annoated in INFO field for prdicted class (somatic, germline, artifact) and by which caller is detected (built-in, external, both). Also, 

