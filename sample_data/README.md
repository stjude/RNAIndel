# Run Example
Here, we demonstrate an analysis example using a sample data prepared from the Jurkat T-cell leukemia cell line.<br>
This cell line harbors two known indels in the *PTEN* tumor suppressor gene: a 2-nt deletion followed by a 9-nt insertion at codon 234 and a 39-nt insertion at codon 246 (Figure 2 in [Reference](#reference)). 
We apply the RNAIndel pipeline to the sample BAM file (sample.bam), which contains the GRCh38 region chr10:80,000,000-90,000,000 (the *PTEN* locus is chr10:87,863,113-87,971,930). 

## Setup
**Step 1:** [Clone](../README.md#download) the RNAIndel repository and [install](../README.md#installation) RNAIndel. <br>
Your direcotry will be structured:
```
RNAIndel
     |_ Bambino
     ...
     |_ RNAIndel
     ...
     |_ sample_data
     ...
```
**Step 2:**  Setup Data Directory.<br> 
Unpack [data_dir.tar.gz](http://ftp.stjude.org/pub/software/RNAIndel/data_dir.tar.gz) under the RNAIndel root.  
```
RNAIndel
    |_ Bambino
    ...
    |_ data_dir
    ...
    |_ RNAIndel
    ...
    |_ sample_data
    ...
```

## Perform analysis
We are assumed to be in the RNAIndel root directory.

```
RNAIndel <----------- We are here.     
    |_ Bambino
    ...
    |_ data_dir
    ...
    |_ RNAIndel
    ...
    |_ sample_data
    ...
```

### Working with the built-in Bambino caller
```
$ ./rna_indel_pipeline.sh -b ./sample_data/inputs/sample.bam \
                          -o ./sample_data/outputs/sample.vcf \
                          -f path/to/your_GRCh38.fa \
                          -d ./data_dir

bambino completed successfully.
rna_indel completed successfully
```
Fifteen coding indels are reported in the ouput [VCF](./outputs/sample.vcf) file: 2 somatic, 1 germline, and 12 artifact indels.
The two *PTEN* indels are predicted as somatic. The first indel is a complex indel in which a 2-nt deletion and a 9-nt insertion 
are involved. This indel is detected as a 7-nt insertion at codon 233 with indel complexity = 2 (See INFO field in [VCF](./outputs/sample.vcf)). The second indel
is detected as a 7-nt insertion at codon 246, not as a 39-nt insertion, due to the soft-clipped alignment. 

### Working with GATK-HaplotypeCaller
The sample BAM file was preprocessed following GATK RNA-Seq Variant Calling [BestPractice](https://software.broadinstitute.org/gatk/documentation/article.php?id=3891). 
GATK-HC (ver 4.0.2.1) called variants in the preprocessed BAM file (./inputs/sample.gatk.bam) and generated a VCF file (./inputs/sample_gatk.vcf).
Now, the indels in the GATK VCF file are classified. **Please input the original BAM file (sample.bam), not the preprocessed one (sample.gatk.bam)**.
```
$ ./rna_indel_pipeline.sh -b ./sample_data/inputs/sample.bam \
                          -c ./sample_data/inputs/sample_gatk.vcf \
                          -o ./sample_data/outputs/sample_gatk_classified.vcf
                          -f path/to/your_GRCh38.fa \
                          -d ./data_dir

rna_indel completed successfully
```
Four coding indels are reported in the output [VCF](./outputs/sample_gatk_classified.vcf) file: 2 somatic, 1 germline, and 1 artifact. 
Both *PTEN* indels are predicated as somatic. The 39-insertion at codon 246 is also detected as a 7-nt insertion by GATK-HC. 
However, the indel at codon 234 is detected as a combination of two insertions: 
the GGCCC insertion at chr10:87957916 and the TG insertion at chr10:87957917. In the output [VCF](./outputs/sample_gatk_classified.vcf), 
these two are annotated as "RQB=chr10:87957916:G:GGGCCCAT". This means that RNAIndel could not find these two indels in the BAM 
file as specified by the input [VCF](./inputs/sample_gatk.vcf) file and, instead, found the GGCCCAT insertion at chr10:87957916 in the BAM file and 
used this for prediction. 

## Reference
Shan, X, Czar, J.C., Bunnell, S.C., Liu, P., Liu, Y., Schwartzberg, P.L., Wange, R.L. (2000) Deficiency of PTEN in Jurkat T Cells Causes Constitutive Localization of Itk to the Plasma Membrane and Hyperresponsiveness to CD3 Stimulation. Mol Cell Biol., 20: 6945–6957. PMID [10958690](https://www.ncbi.nlm.nih.gov/pubmed/10958690)      
