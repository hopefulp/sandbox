#!/usr/bin/bash

for i in {6..20..4}; do
    ptname=hl${i}${i}${i}${i}
    echo "mllorenz.py tr -hl $i $i $i $i -ms ${ptname}.pt -mi 1000000 -nd 200000 -dbp 2 > ${ptname}.out"
    mllorenz.py tr -hl $i $i $i $i -ms ${ptname}.pt -mi 1000000 -nd 200000 -dbp 2 > ${ptname}.out
    done


