### Train RNAIndel
RNAIndel trains the model based on the truth annotated by the user.

#### Step 1 (feature calculation)
Features are calculated for each indel and reported in a tab-delimited file.<br>
Using a callset by the built-in caller, 
```
rnaindel feature -b BAM -o OUTPUT_TAB -f FASTA -d DATA_DIR [other options]
```
To train based on a callset from your caller, specify the input VCF by ```-v```.
```
rnaindel feature -b BAM -v INPUT_VCF -o OUTPUT_TAB -f FASTA -d DATA_DIR [other options]
```
To train based on a germline database of your choice, specify the database by ```-g```.
```
rnaindel feature -b BAM -o OUTPUT_TAB -f FASTA -d DATA_DIR -g YOUR_DB [other options]
```
The germline databse is expected to be a tabixed VCF file with no missing value in ID field.

See [options](#../../#Options) for detail.
 
#### Step 2 (annotation)
The output tab-delimited file has a column \"truth\". Users annotate each indel
by filling the column with either of <br> 
\"somatic\", \"germline\", or \"artifact\". 

#### Step 3 (update models)
Repeat Step 1 and 2 for N samples.<br>
Users concatenate the annotated files. Here, assuming the files are \"sample.i.tab\" (i = 1,...,N), 
```
head -1 sample.1.tab > training_set.tab           # keep the header line
```
```
tail -n +2 -q sample.*.tab > training_set.tab     # concatenate files without header
```
The concatenated file is used as a training set to update the models.
Specify the indel class to be trained by -c. 
```
rnaindel training -t TRAINING_SET -d DATA_DIR -c INDEL_CLASS [other options]
```
#### Options
* ```-t``` training set with annotation (required)
* ```-d``` [data directory](#setup) contains trained models and databases (required) 
* ```-c``` indel class to be trained. s for single-nucleotide indel and m for multi-nucleotide indel (required)
*<details>
     <summary>other options (click to open)</summary>
    * ```-k``` number of folds in k-fold cross-validation (default: 5)
    * ```-p``` number of processes (default: 1)
    * ```-l``` directory to ouput log files (default: current)
    * ```-ds-beta``` F beta to be optimized in down sampling step. Optimized for TPR if beta > 100. (default: 10)
    * ```-fs-beta``` F beta to be optimized in feature selection step. Optimized for TPR if beta > 100. (default: 10)
    * ```-pt-beta``` F beta to be optimized in parameter tuning step. Optimized for TPR if beta > 100. (default: 10)
    * ```--downsample-ratio``` Train with a user-specified downsample ratio: integer between 1 and 20. (default: None)
    * ```--feature-names``` Train with a user-specified subset of features. Supply feature names in a text file containing a feature name per line (default: None)
    * ```--auto-param``` Train with sklearn.RandomForestClassifer's max_features="auto" (default: False)
*</details>
