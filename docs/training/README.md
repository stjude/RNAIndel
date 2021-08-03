# Train RNAIndel
RNAIndel trains the model using your training set.

### Step 1 (feature calculation)

Features are calculated for each indel and reported in a tab-delimited file.<br>
Suppose we have N samples. For i-th sample:
```
rnaindel CalculateFeatures  -i sample.i.bam \
                            -o sample.i.tab \
                            -r reference.fa \
                            -d ./data_dir_grch38\
                            [-v sample.i.external.vcf.gz]
```

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
Concateate the annotated files .  
```
head -1 sample.1.tab > training_set.tab           # keep the header line
```
```
tail -n +2 -q sample.*.tab > training_set.tab     # concatenate files without header
```
The concatenated file is used as a training set to update the models.
Specify the indel class to be trained by ```-c```. 
```
rnaindel Train -t training_set.tab -d ./data_dir_grch38 -c indel_class_to_train [other options]
```
#### Options
* ```-t``` training set with annotation (required)
* ```-d``` [data directory](../../README.md/#setup) contains trained models and databases (required) 
* ```-c``` indel class to be trained. "s" for single-nucleotide indel, "m" for multi-nucleotide indel, "h" for homopolymer indel(required)
* <details>
    <summary>other options (click to open)</summary><p>
    
    * ```-k``` number of folds in k-fold cross-validation (default: 5)
    * ```-p``` number of processes (default: 1)
    * ```-l``` directory to ouput log files (default: current)
    * ```--ds-beta``` F beta to be optimized in down sampling step. Optimized for TPR if beta > 100. (default: 10)
    * ```--fs-beta``` F beta to be optimized in feature selection step. Optimized for TPR if beta > 100. (default: 10)
    * ```--pt-beta``` F beta to be optimized in parameter tuning step. Optimized for TPR if beta > 100. (default: 10)
    * ```--downsample-ratio``` train with a user-specified downsample ratio: integer between 1 and 20. (default: None)
    * ```--feature-names``` train with a user-specified subset of features: [input example](../../sample_data/inputs/feature_names.txt) (default: None)
    * ```--auto-param``` train with sklearn.RandomForestClassifer's max_features="auto" (default: False)

</p></details>
