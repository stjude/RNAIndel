# Run Example
Here, we demonstrate an analysis example using a sample data prepared from the Jurkat T-cell leukemia cell line.<br>
This cell line harbors two known indels in the PTEN tumor suppressor gene: codon 234 and codon 246 (Figure 2 in [Reference](#reference)). We apply the RNAIndel pipeline to the sample BAM file (sample.bam), whcih contains the GRCh38 region chr10:80,000,000-90,000,000 (the PTEN locus is chr10:87,863,113-87,971,930). 

## Setup
**Step 1:** [Clone](../README.md) the RNAIndel repository and [install](../README.md) RNAIndel. <br>
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
Unpack [data_dir.tar.gz](http//ftp.stjude.org/pub/software/RNAIndel/data_dir.tar.gz) under the RNAIndel root.  
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
### Working with the built-in Bambino caller
```
$ ./rna_indel_pipeline.sh -b sample_data/sample.bam \
                          -o sample_data/sample.vcf \
                          -f path/to/your_GRCh38.fa \
                          -d ./data_dir

bambino completed successfully.
rna_indel completed successfully
```
Fifteen indels are reported in the ouput [VCF](sample.vcf) file: 2 somatic, 1 germline, and 12 artifact indels.
The two PTEN indels are predicted as somatic. The first indel is annotated in codon 233 (not 234) for NM_000314.
This is due to the 2-nt deletion (AC) at chr10:87957916, which is also described in [Reference](#reference). <br>

### Working with GATK-HaplotypeCaller
The sample BAM file was preprocessed following GATK RNA-Seq Variant Calling [BestPractice](https://software.broadinstitute.org/gatk/documentation/article.php?id=3891). 
GATK-HC (ver 4.0.2.1) called variants in the preprocessed BAM file (sample.gatk.bam) and generated a VCF file (sample_gatk.vcf). <br>
Now, the indels in the GATK VCF file are classified. **Please note that the original BAM file (sample.bam), not the preprocessed one, is used**.
```
$ ./rna_indel_pipeline.sh -b sample_data/sample.bam \ **NOT sample.gatk.bam !!**  
                          -c sample_gatk.vcf \
                          -o sample_data/sample_gatk_classified.vcf
                          -f path/to/your_GRCh38.fa \
                          -d ./data_dir

rna_indel completed successfully
```
Four indels are reported in the output [VCF](sample_gatk_classified.vcf) file: 2 somatic, 1 germline, and 1 artifact. 
Both PTEN indels are predicated as somatic. However, the first indel was detected by GATK as a combination of two insertions: 
the GGCCC insertion at chr10:87957916 and the TG insertion at chr10:87957917. In the output, these two are annotated as "RQB=chr10:87957916:G:GGGCCCAT". 
This means that RNAIndel could not find these indels in the BAM file as represented by the input VCF file and, instead, used 
the GGCCCAT insertion at chr10:87957916 for prediction, which is complexed with the AC>GG substitution at chr10:87957917-18. 

## Reference
Shan, X, Czar, J.C., Bunnell, S.C., Liu, P., Liu, Y., Schwartzberg, P.L., Wange, R.L. (2000) Deficiency of PTEN in Jurkat T Cells Causes Constitutive Localization of Itk to the Plasma Membrane and Hyperresponsiveness to CD3 Stimulation. Mol Cell Biol., 20: 6945–6957. PMID [10958690](https://www.ncbi.nlm.nih.gov/pubmed/10958690)      
