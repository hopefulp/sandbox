#!/usr/bin/bash

com=$1
shift

echo python $(which $com) $@
python $(which $com) $@
