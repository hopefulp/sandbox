#!/bin/bash

grep molecules $1

if [ $? -eq 0 ]; then
    sed -i  -e '/molecules/ i\#include \"amber03.ff/ions.itp\"\n' $1
    echo "\"amber03.ff/ions.itp\" was included"
else
    echo "nothing has been done"
fi

