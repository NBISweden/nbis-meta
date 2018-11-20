#!/usr/bin/env bash

bgcolor=$2

if [ "$bgcolor" == "" ] ; then
  bgcolor='transparent'
fi

cat $1 | sed "s/bgcolor=white/bgcolor=$bgcolor/g" | sed 's/rounded"/rounded", fontcolor=white/g' | sed 's/rounded,dashed"/rounded,dashed", fontcolor=white/g' | dot -Tpng -Gdpi=300