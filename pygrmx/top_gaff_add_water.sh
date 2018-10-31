#!/bin/bash

# replace [ defaults ] with include
grep defaults $1
if [ $? -eq 0 ]; then
    sed -i -e '3,6d' $1
    echo " [ defaults ] is deleted"
    sed -i -e '/atomtypes/ i\#include \"amber03.ff/forcefield.itp\"\n' $1
    echo "\"amber03.ff/forcefield.itp\" is included"
else
    echo " nothing happend "
fi    

grep water $1
if [ $? -eq 0 ]; then
    sed -i -e '/water/ a\#include \"amber03.ff/tip3p.itp\"' $1
    echo "\"amber03.ff/tip3p.itp\" is included"
else
    echo " nothing happend "
fi    
