#!/bin/bash



sed -i  -e '/molecules/ i\#include \"amber03.ff/ions.itp\"\n' $1
