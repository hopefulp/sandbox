#!/bin/bash



sed -i -e '3,6d'  -e '/atomtypes/ i\#include \"amber03.ff/forcefield.itp\"\n' -e '/water/ a\#include \"amber03.ff/tip3p.itp\"' $1
