#!/bin/bash

prefix=mol
suffix=

INCAR=incar.520.fine.d2.mol


for model in `cat mol2.txt`
    do
        fname=$prefix-$model.pos
        new=$prefix-$model
        echo "file: $fname  dir: $new"
#        ./job_mkdir_sub_mol.sh $fname $new Fe  1 $model

	if [ -e $new ]; then
    	    echo "There is already $new dir"
    	    exit
	fi

	mkdir $new

   	cp ./Pot/pawpbe_h.pot            $new/POTCAR
	cp $fname                         $new/POSCAR
	cp ./$INCAR                       $new/INCAR
	cp ./kp1.gamma                     $new/KPOINTS

	~/bin/changeline_pbs.pl pbsfold-idft.d2.csh $new nodes=1
    done

