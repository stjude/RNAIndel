############################################################
# Purpose   : Dockerize RNAIndel
# OS        : Ubuntu 18.04
# VERSION   : 0.3.0
############################################################
FROM adamdingliang/rnaindel:base
MAINTAINER Liang.Ding@stjude.org

WORKDIR /RNAIndel
COPY . .
RUN python3 setup.py install

CMD ["echo", "rnaindel image created"]
