FROM python:3.6

WORKDIR /data
ADD . .

RUN pip install .

ENTRYPOINT ["rnaindel"]
CMD ["-h"]
