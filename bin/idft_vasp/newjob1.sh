#!/bin/bash

old=Relax/$1
new=$2

if [ ! -e $old ]; then
    echo "There is no $old dir"
    exit
fi

if [ -e $new ]; then
    echo "There is already $new dir"
    exit
fi

mkdir $new

cp $old/POTCAR  $new
cp $old/POSCAR  $new
cp ./incar.dos  $new/INCAR
cp ./kp_dos     $new/KPOINTS
cp $old/CHGCAR  $new
cp $old/WAVECAR $new
