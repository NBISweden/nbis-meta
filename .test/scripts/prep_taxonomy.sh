#!/bin/bash

mkdir -p resources/uniref100
gzip -c .test/data/uniref100.fasta > resources/uniref100/uniref100.fasta.gz
