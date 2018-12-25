# Indel calling executable
The Bambino variant caller is optimized to call SNVs and indels from RNA-Seq data.

## Prerequisites
* [java=1.8.0_66](https://www.java.com/en/download/)

## Variant calling  
```
bambino -i input.bam \
        -f reference.fa \
        -o bambino_call.txt \ 
        [optional argument (-m)]
```
### Bambino options
* ```-m``` maximum heap space (default 6000m)
* ```-b``` input BAM file (required)
* ```-f``` reference genome FASTA file (required)
* ```-o``` Bambino output file. Tab-delimited flat file (required)

