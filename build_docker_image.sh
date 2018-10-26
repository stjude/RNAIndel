#!/usr/bin/env bash

version_num=$(grep -o '".*"' version.py | sed 's/"//g')
echo "Version: $version_num" >&2

usage() {
   echo "Usage:"
   echo "    build_docker_image <SUBCOMMAND> [args...]"
   echo "Subcommands:"
   echo "    base           build a base image with all the dependencies installed for Bambino and RNAIndel"
   echo "    rna_indel      build an image for rna_indel"
   exit 1
}


if [ "$#" == 0 ]; then
    usage
    return 1
    exit 1
fi


build_base() {
    echo "Building a base image ..." >&2
    docker build -t rna_indel_base:$version_num -f Dockerfile.base .
    echo "Done" >&2
}


build() {
    echo "Building an image for rna_indel ..." >&2
    docker build -t rna_indel:$version_num .
    echo "Done" >&2
}


SUBCOMMAND=$1

case $SUBCOMMAND in
    base) build_base ;;
    rna_indel) build ;;
    *) usage ;;
esac
