############################################################
# Purpose   : Dockerize RNAIndel and Bambino
# OS        : Ubuntu 18.04
# VERSION   : 0.1.0
############################################################
FROM adamdingliang/rnaindel:base
MAINTAINER Liang.Ding@stjude.org

WORKDIR /RNAIndel
COPY . .
RUN python3 setup.py install

CMD ["echo", "rna_indel image created"]