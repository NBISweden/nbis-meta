#!/bin/bash

base="https://github.com/eggnogdb/eggnog-mapper/blob/master/tests/fixtures"

mkdir -p resources/eggnog-mapper
curl -L -o resources/eggnog-mapper/eggnog.db "${base}/eggnog.db?raw=true"
curl -L -o resources/eggnog-mapper/eggnog_proteins.dmnd "${base}/eggnog_proteins.dmnd?raw=true"
touch resources/eggnog-mapper/eggnog.version
