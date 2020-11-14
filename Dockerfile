FROM openjdk:8

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/* 

VOLUME /data
WORKDIR /tmp

ADD . .

ENV PYTHONWARNINGS="ignore"
RUN python3 -m pip install . && rm -rf /tmp/*

WORKDIR /data
ENTRYPOINT ["rnaindel"]
CMD ["-h"]
