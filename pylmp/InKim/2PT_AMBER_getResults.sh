#!/bin/bash

curr_dir=`pwd`

infile=`ls -1 *grps | head -n 1`
suffix=${infile%.grps}
filelist=`ls -1 *thermo`
thermofile=`ls -1 *thermo | head -n 1`
title1=`grep property $thermofile`
title2=`grep nmolecules $thermofile`
title3=`grep natom $thermofile`

Sq=${curr_dir}/${suffix}.Sq.dat
Aq=${curr_dir}/${suffix}.Aq.dat
Eq=${curr_dir}/${suffix}.Eq.dat

# Write header
for file in $Sq $Aq $Eq
do
	echo $curr_dir > $file
	echo $title1 >> $file
	echo $title2 >> $file
	echo $title3 >> $file
	echo >> $file
done

echo "Sq_(J/mol_K/SimBox)" >> $Sq
echo "Aq__(kJ/mol/SimBox)" >> $Aq
echo "Eq__(kJ/mol/SimBox)" >> $Eq

# Grep data
echo "Getting 2PT results.."
for file in $filelist
do
	fnumber=${file#restrt_}
	fnumber=${fnumber%.2pt.thermo}

	echo $fnumber `grep Sq_ $file | sed -e "s/${suffix}.npt.//" -e "s/.2pt.mol.grps.thermo://" -e "s/Sq_(J\/mol_K\/SimBox)//"` >> $Sq
	echo $fnumber `grep Aq_ $file | sed -e "s/${suffix}.npt.//" -e "s/.2pt.mol.grps.thermo://" -e "s/Aq__(kJ\/mol\/SimBox)//"` >> $Aq
	echo $fnumber `grep Eq_ $file | sed -e "s/${suffix}.npt.//" -e "s/.2pt.mol.grps.thermo://" -e "s/Eq__(kJ\/mol\/SimBox)//"` >> $Eq
done

echo "Done.. check your *dat file."

### end of 2PT_getResults.sh
