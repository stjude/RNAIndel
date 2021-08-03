# Indel filtration 
Except for a few known hotspots, true somatic indels rarely recurr.   

### Recurrent somatic indels
When multiple RNAIndel output VCF files are generated from the same cohort, occurrence with in the cohort is annotated in INFO field. 

```
rnaindel CountOccurrence --vcf-list vcf_paths.txt -r reference.fa
```

#### Options
* ```--vcf-list``` file containing paths to RNAIndel output VCF files (path/line) (required)
* ```-r``` reference genome FASTA file (required)
* ```--out-dir``` output directory for annotated VCF file. The input file dir is used if not specified. <br>
&nbsp;   &nbsp;   &nbsp;   &nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp;  &nbsp; (file name will not be changed after annotation)
