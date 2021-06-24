RNAIndel version 3 


This version needs [indelPost](https://github.com/stjude/indelPost).
Please install this first.

Also, depends on [scikit-learn>=0.22](http://scikit-learn.org/stable/install.html#).


To install
```
pip install indelpost

git clone https://github.com/rawagiha/RNAIndel.git
cd RNAIndel
python setup.py install
```

To run somatic analysis:
```
rnaindel PredictIndels -i RNA-Seq BAM -d data_dir -r ref.fa -o output.vcf -p 8 (default 1)
``` 
