#!/bin/sh

o_model=$1
n_model=$2

V_mkdir.sh $o_model $n_model 
cp $n_model.pos $n_model/POSCAR

 
