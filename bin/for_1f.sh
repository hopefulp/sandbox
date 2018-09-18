#!/bin/bash

for Me in `cat metal.dir`
  do
    ./insert_mol2mof.pl CONTCAR.MOF74.$Me CONTCAR outf=MOF74-$Me-Ace.pos "C 2 H 2"
  done
