#!/bin/bash

fn=$1

mkdir $fn
cp $fn.inp $fn
cd $fn
~/bin/bsse_2mol.pl $fn $2 $3
