#!/usr/bin/env bash 
# Modified from Nathan S. Watson-Haigh's script https://gist.github.com/nathanhaigh/3521724#file-deinterleave_fastq-sh
# Added check for gzipped input


# Deinterleaves a FASTQ file of paired reads into two FASTQ
# files specified on the command line. Optionally GZip compresses the output
# FASTQ files using pigz if the 3rd command line argument is the word "compress"
# 
# Can deinterleave 100 million paired reads (200 million total
# reads; a 43Gbyte file), in memory (/dev/shm), in 4m15s (255s)
# 
# Latest code: https://gist.github.com/3521724
# Also see my interleaving script: https://gist.github.com/4544979
# 
# Inspired by Torsten Seemann's blog post:
# http://thegenomefactory.blogspot.com.au/2012/05/cool-use-of-unix-paste-with-ngs.html

# Set up some defaults
GZIP_OUTPUT=0

if [ "$1" == "" ] || [ "$2" == "" ] || [ "$3" == "" ] ; then
    echo "Usage: deinterleave_fastq.sh <interleaved.fastq[.gz]> <f.fastq> <r.fastq> [compress]"
    exit 0
fi

zipper=`which pigz`

if [ "$pigz_path" == "" ] ; then
  zipper=`which gzip`
fi

GZIPPED_INPUT=`file -L $1 | grep -c "gzip compressed data"`

# If the third argument is the word "compress" then we'll compress the output using pigz
if [[ $4 == "compress" ]]; then
    GZIP_OUTPUT=1
fi

if [[ ${GZIP_OUTPUT} == 0 ]]; then
    if [ $GZIPPED_INPUT == 1 ] ; then
        gunzip -c $1 | paste - - - - - - - -  | tee >(cut -f 1-4 | tr "\t" "\n" > $2) | cut -f 5-8 | tr "\t" "\n" > $3
    else
        cat $1 | paste - - - - - - - -  | tee >(cut -f 1-4 | tr "\t" "\n" > $2) | cut -f 5-8 | tr "\t" "\n" > $3
    fi
else
    if [ $GZIPPED_INPUT == 1 ] ; then
        gunzip -c $1 |  paste - - - - - - - -  | tee >(cut -f 1-4 | tr "\t" "\n" | $zipper --best > $2) | cut -f 5-8 | tr "\t" "\n" | $zipper --best > $3
    else
        cat $1 | paste - - - - - - - -  | tee >(cut -f 1-4 | tr "\t" "\n" | $zipper --best > $2) | cut -f 5-8 | tr "\t" "\n" | $zipper --best > $3
    fi
fi
