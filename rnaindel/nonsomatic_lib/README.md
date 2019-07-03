# Indel filtration 

RNAIndel has subcommands to filter artifact and germline indels predicted as somatic.

### Compile non-somatic indel panel
RNAIndel compiles a non-somatic indel panel using common non-COSMIC indels found in matched or public normal samples.
```
rnaindel nonsomatic --vcf-list FILE --count INT -o NONSOMATIC_VCF -r REFERENCE -d DATA_DIR
bgzip NONSOMATIC_VCF
tabix -p vcf NONSOMATIC_VCF.gz
```
#### Options
* ```--vcf-list``` file containing paths to normal VCF files ([example](../../sample_data/inputs/normals.txt)) (required)
* ```--count``` indel occurrence count to be defined as common in the supplied VCF files (required)
* ```-o``` output non-somatic VCF file (required)
* ```-r``` reference genome FASTA file (required)
* ```-d``` [data directory](../../README.md/#setup) contains the COSMIC database (required)
 
### Reclassify indels found in the panel
```
rnaindel reclassification -i RNAIndel_OUTPUT_VCF -o RECLASSIFIED_VCF -r REFERENCE -n NON_SOMATIC_VCF.gz 
```
#### Options
* ```-i``` RNAIndel output VCF file to be reclassified (required).
* ```-o``` reclassifed VCF file (required)
* ```-r``` reference genome FASTA file (required)
* ```-n``` non-somatic indel panel in bgzip-compressed VCF file (required)

### Annotate non-COSMIC recurrent indels
When multiple RNAIndel output VCF files are available as in a case of cohort study, 
RNAIndel annotates indels recurring in the samples that are not present in the 
COSMIC databse. Such recurrent indels are possibly common artifact/germline rather 
than a somatic hotspot. Annotation is made in INFO field. 
```
ranindel recurrence --vcf-list FILE -r REFERENCE -d DATADIR --out-dir
```
#### Options
* ```--vcf-list``` file containing paths to RNAIndel output VCF files ([example](../../sample_data/inputs/rnaindel_vcfs.txt)) (required)
* ```-r``` reference genome FASTA file (required)
* ```-d``` [data directory](../../README.md/#setup) contains the COSMIC database (required)
* ```--out-dir``` output direcory for annotated VCF file. The input file dir is used if not specified.
                 (file name will not be changed after annotation)
