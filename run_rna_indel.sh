#!/usr/bin/env bash

# Wrapper script to chain bambino and rna_indel executables
# Only rna_indel will be run when an input VCF file is supplied
#
# Author: Kohei Hagiwara
#

usage() {
    >&2 echo ""
    >&2 echo ""
    >&2 echo "This wrapper requires bambino and rna_indel executables installed."
    >&2 echo ""
    >&2 echo "Usage:"
    >&2 echo "  To run the Bambino/RNAIndel pipeline:"
    >&2 echo "      ./run_rna_indel.sh -b yourBAM.bam -o yourOUT.vcf -f yourFASTA.fa -d dataDir [optional parameters]"
    >&2 echo ""
    >&2 echo "  To classify indel entries in your VCF file:"
    >&2 echo "      ./run_rna_indel.sh -b yourBAM.bam -c yourINPUT.vcf -o yourOUT.vcf -f yourFASTA.fa -d dataDir [optional parameters]"
    >&2 echo ""
    >&2 echo "Arguments:"
    >&2 echo "-b BAM file (required)"
    >&2 echo "-c input VCF file (only required for using other callers)"
    >&2 echo "-o oupput VCF file (required)"
    >&2 echo "-f reference FASTA file (GRCh38) (required)"
    >&2 echo "-d data directory contains refgene, dbsnp and clinvar databases (required)"
    >&2 echo "-q STAR mapping quality MAPQ for unique mappers (default 255)"
    >&2 echo "-p number of cores (default 1)"
    >&2 echo "-m maximum heap space. Java parameter used for Bambino (default 6000m)"
    >&2 echo "-n user-defined panel of non-somatic indels in VCF format"
    >&2 echo "-l directory to store log files"
    >&2 echo "-h show this message"
    >&2 echo ""
    >&2 echo ""
}

# show usage when run w/o args
if [ $# -eq 0 ]; then 
    usage
    exit 0
fi

# parse args
while getopts "b:c:o:f:d:q:p:m:n:l:" opt
do
    case $opt in
        b) BAM=$OPTARG ;;
        c) INVCF=$OPTARG ;;
        o) OUTVCF=$OPTARG ;;
        f) FASTA=$OPTARG ;;
        d) DATADIR=$OPTARG ;;
        q) MAPQ=$OPTARG ;;
        p) PROCESS=$OPTARG ;;
        m) HEAP=$OPTARG ;;
        n) PONS=$OPTARG ;;
        l) LOG=$OPTARG ;;
        h) usage ;;
    esac
done
shift $(($OPTIND - 1))

# check params
if [ ! $BAM ]; then
    >&2 echo "Error: input BAM file is required."
    exit 1
fi

if [ ! -f $BAM ]; then
    >&2 echo "Error: input BAM file is not found."
    exit 1
fi

if [ ! -f $INVCF ]; then
    >&2 echo "Error: input VCF file is not found."
    exit 1
fi

if [ ! $OUTVCF ]; then
    >&2 echo "Error: output VCF file must be specified."
    exit 1
fi

if [ ! $FASTA ]; then
    >&2 echo "Error: reference FASTA file is required."
    exit 1
fi

if [ ! -f $FASTA ]; then
    >&2 echo "Error: reference FASTA file is not found."
    exit 1
fi

if [ ! $DATADIR ]; then
    >&2 echo "Error: Data directory is required."
    exit 1
fi

if [ ! -d $DATADIR ]; then
    >&2 echo "Error: Data directory is not found."
    exit 1
fi

if [ $MAPQ ]; then
    if [[ ! $MAPQ =~ ^-?[0-9]+$ ]] || [ $MAPQ -le 0 ] || [ $MAPQ -gt 255 ]; then
        >&2 echo "Error: MAPQ must be an integer between 1 and 255."
        exit 1
    fi
else
    MAPQ=255
fi 


if [ $PROCESS ]; then
    if [[ ! $PROCESS =~ ^-?[0-9]+$ ]] || [ $PROCESS -le 0 ]; then
        >&2 echo "Error: the number of cores must be a positive integer."
        exit 1
    fi
else
    PROCESS=1
fi

if [ ! $HEAP ];then
    HEAP=6000m
fi

if [ ! -f $PONS ]; then
    >&2 echo "Error: panel of non-somatic indels is not found."
    exit 1
fi

if [ $LOG ]; then
    if [ ! -d $LOG ]; then
        >&2 echo "Error: Log directory is not found."
        exit 1
    fi
else
    LOG=`pwd`
fi

# when input VCF is not supplied -> run Bambino & RNAIndel
if [ ! $INVCF ]; then
    BAMBINO_OUT=$(mktemp)
    bambino -b $BAM -f $FASTA -o $BAMBINO_OUT -m $HEAP
    if [ ! $PONS ]; then
        rna_indel -b $BAM -i $BAMBINO_OUT -o $OUTVCF -f $FASTA -d $DATADIR -q $MAPQ -p $PROCESS -l $LOG
    else
        rna_indel -b $BAM -i $BAMBINO_OUT -o $OUTVCF -f $FASTA -d $DATADIR -q $MAPQ -p $PROCESS -l $LOG -n $PONS
    fi
    rm $BAMBINO_OUT
# input VCF is supplied skip Bambino
else
   if [ ! $PONS ]; then
       rna_indel -b $BAM -c $INVCF -o $OUTVCF -f $FASTA -d $DATADIR -q $MAPQ -p $PROCESS -l $LOG
   else
       rna_indel -b $BAM -c $INVCF -o $OUTVCF -f $FASTA -d $DATADIR -q $MAPQ -p $PROCESS -l $LOG -n $PONS
   fi
fi
