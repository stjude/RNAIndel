FROM openjdk:8

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    libbz2-dev \
    liblzma-dev \
    && rm -rf /var/lib/apt/lists/* 

RUN mkdir /data
VOLUME /data
WORKDIR /tmp

ADD . .

RUN python3 -m pip install --upgrade pip setuptools wheel

RUN pip install cython numpy scikit-learn==1.3.2

RUN git clone https://github.com/libnano/ssw-py.git && cd ssw-py && python3 setup.py install && cd /tmp

ENV PYTHONWARNINGS="ignore"
RUN python3 -m pip install . && rm -rf /tmp/*
COPY rnaindel/RNAIndel.sh /usr/local/bin/RNAIndel.sh

WORKDIR /data
ENTRYPOINT ["RNAIndel.sh"]
CMD ["-h"]
