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
done

rnaindel "$@"
