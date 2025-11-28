#!/bin/bash

for d in */ ; do
    if [[ -f "${d}CONTCAR" ]]; then
        base=$(basename "$d")
        cp "${d}CONTCAR" "CONTCAR.${base}"
        echo "Copied ${d}CONTCAR → CONTCAR.${base}"
    fi
done

