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

Download the datapackage(135MB). The current version requires the new package:
```
curl -LO http://ftp.stjude.org/pub/software/RNAIndel/data_dir_grch38.v3.tar.gz
tar -zxf data_dir_grch38.v3.tar.gz
```

To run somatic analysis:
```
rnaindel PredictIndels -i RNA-Seq BAM -d data_dir -r ref.fa -o output.vcf -p 8 (default 1)
``` 
