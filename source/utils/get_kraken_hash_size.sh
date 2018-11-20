#!/usr/bin/env bash

lib=$1
kernel=`uname -s`

if [ "$kernel" == "Linux" ] ; then
  KRAKEN_HASH_SIZE=$(find $lib/ '(' -name '*.fna' -o -name '*.fa' -o -name '*.ffn' ')' -printf '%s\n' | perl -nle '$sum += $_; END {print int(1.15 * $sum)}')
elif [ "$kernel" == "Darwin" ] ; then
  KRAKEN_HASH_SIZE=$(find $lib/ '(' -name '*.fna' -o -name '*.fa' -o -name '*.ffn' ')' -print0 | xargs -0 stat -f '%i '| perl -nle '$sum += $_; END {{print int(1.15 * $sum)}}')
fi

echo $KRAKEN_HASH_SIZE