FROM openjdk:8

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/* 

RUN mkdir /data
VOLUME /data
WORKDIR /tmp

ADD . .
RUN curl -LO http://ftp.stjude.org/pub/software/RNAIndel/sample_data.tar.gz && tar xzf sample_data.tar.gz -C /data
RUN curl -LO http://ftp.stjude.org/pub/software/RNAIndel/data_dir_38.tar.gz && tar xzf data_dir_38.tar.gz -C /data

ENV PYTHONWARNINGS="ignore"
RUN python3 -m pip install . && rm -rf /tmp/*

WORKDIR /data
RUN chmod 777 -R /data
ENTRYPOINT ["rnaindel"]
CMD ["-h"]
