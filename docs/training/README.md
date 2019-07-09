# Train RNAIndel
RNAIndel trains the model using your training set.

### Step 1 (feature calculation)
Features are calculated for each indel and reported in a tab-delimited file.<br>

By default, RNAIndel calculates features based on a call set by the built-in caller. 
```
rnaindel feature -i BAM -o OUTPUT_TAB -r REFERENCE -d DATA_DIR [other options]
```
To calculate for your caller's call set, specify the input VCF from your caller with ```-v```.
```
rnaindel feature -i BAM -v INPUT_VCF -o OUTPUT_TAB -r REFERENCE -d DATA_DIR [other options]
```
To use your germline indel database for training, specify the database with ```-g```. <br>
The database is expected to be a bgzip-compressed VCF file with no missing values in the ID field.
```
rnaindel feature -i BAM -o OUTPUT_TAB -r REFERENCE -d DATA_DIR -g YOUR_DB [other options]
```

#### Options
* ```-i``` input [STAR](https://academic.oup.com/bioinformatics/article/29/1/15/272537)-mapped BAM file (required)
* ```-o``` output tab-delimited file (required)
* ```-r``` reference genome (GRCh37 or 38) FASTA file (required)
* ```-d``` [data directory](../../README.md/#setup) contains trained models and databases (required)
* ```-v``` VCF file from user's caller (default: None)
* ```-g``` user-provided germline indel database in tabixed VCF format (default: None)
* <details>
    <summary>other options (click to open)</summary><p>
    
    * ```-q``` STAR mapping quality MAPQ for unique mappers (default: 255)
    * ```-p``` number of cores (default: 1)
    * ```-m``` maximum heap space (default: 6000m)
    * ```-l``` direcotry to store log files (default: current)
    * ```-n``` user-defined panel of non-somatic indels in tabixed VCF format (default: built-in reviewed indel set)
    * ```--exclude-softclipped-alignments``` softclipped indels will not be used for analysis if used (default: False)

</p></details>

### Step 2 (annotation)
The output tab-delimited file has \"truth\" column. Users annotate each indel by filling the column.
Possible values are:
```
somatic, germline, artifact 
```
<br>
Repeat Step 1 and 2 for N samples.
<br>

### Step 3 (update models)
Users concatenate the annotated files. Here, assuming the files are named \"sample.i.tab\" (i = 1,...,N), 
```
head -1 sample.1.tab > training_set.tab           # keep the header line
```
```
tail -n +2 -q sample.*.tab > training_set.tab     # concatenate files without header
```
The concatenated file is used as a training set to update the models.
Specify the indel class to be trained by ```-c```. 
```
rnaindel training -t TRAINING_SET -d DATA_DIR -c INDEL_CLASS [other options]
```
#### Options
* ```-t``` training set with annotation (required)
* ```-d``` [data directory](../../README.md/#setup) contains trained models and databases (required) 
* ```-c``` indel class to be trained. "s" for single-nucleotide indel and "m" for multi-nucleotide indel (required)
* <details>
    <summary>other options (click to open)</summary><p>
    
    * ```-k``` number of folds in k-fold cross-validation (default: 5)
    * ```-p``` number of processes (default: 1)
    * ```-l``` directory to ouput log files (default: current)
    * ```-ds-beta``` F beta to be optimized in down sampling step. Optimized for TPR if beta > 100. (default: 10)
    * ```-fs-beta``` F beta to be optimized in feature selection step. Optimized for TPR if beta > 100. (default: 10)
    * ```-pt-beta``` F beta to be optimized in parameter tuning step. Optimized for TPR if beta > 100. (default: 10)
    * ```--downsample-ratio``` train with a user-specified downsample ratio: integer between 1 and 20. (default: None)
    * ```--feature-names``` train with a user-specified subset of features: [input example](../../sample_data/inputs/feature_names.txt) (default: None)
    * ```--auto-param``` train with sklearn.RandomForestClassifer's max_features="auto" (default: False)

</p></details>
