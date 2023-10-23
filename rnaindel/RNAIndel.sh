#!/usr/bin/env bash

BAMFILE=""
DATA_DIR=
NEXT=0

for var in "$@"
do
  if [[ $var == *.bam ]]
  then
    BAMFILE=$var
    if [[ -f "${BAMFILE}" ]]
    then 
      bam_dir=$(dirname ${BAMFILE})
      bam_name=$(basename ${BAMFILE} ".bam")
      if [[ -f "${bam_dir}/${bam_name}.bai" ]]
      then
        ln -s ${bam_dir}/${bam_name}.bai ${BAMFILE}.bai
      fi
     fi
  fi
  if [[ $NEXT == 1 ]]
  then
    NEXT=0
    DATA_DIR=$var
  fi
  if [[ $var == '-d' ]]
  then
    NEXT=1
  fi
done

if [[ -d $DATA_DIR ]]
then
  rnaindel SetUp -d $DATA_DIR
fi
rnaindel "$@"
