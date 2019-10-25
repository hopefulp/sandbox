#!/bin/bash

curr_dir=`pwd`
usage="2PT_getResult.sh (result_dir)"

# Check whether result directory is specified
if [ $# -eq 0 ]
then
	result_dir="results"
else
	result_dir=$1
fi

# Find results directory
if [ -d $result_dir ]
then
	infile=`ls -1 in* | head -n 1`
	suffix=${infile#in.}
	cd $result_dir
	echo 2PT_getResult.sh: successfully found the 2PT thermo directory $result_dir
else
	cd ..
	if [ -d $result_dir ]; then
		infile=`ls -1 in* | head -n 1`
		suffix=${infile#in.}
		cd $result_dir
		echo 2PT_getResult.sh: successfully found the 2PT thermo directory $result_dir
	else
		echo 'ERROR: cannot find result directory'
		exit
	fi
fi

ls -1 *thermo > /dev/null 2>&1
if [ "$?" = "0" ]; then
	echo "2PT_getResult.sh: ready to run the scripts"
else
	echo "ERROR: did not find any *thermo file"
	exit
fi

filelist=`ls -1 *thermo | sort -n -k 3 -t '.'`
thermofile=`ls -1 *thermo | sort -n -k 3 -t '.' | head -n 1`
outfile=`ls -1 *out | sort -n -k 3 -t '.'`
title1=`grep property $thermofile`
title2=`grep nmolecules $thermofile`
title3=`grep natom $thermofile`

Sq=${curr_dir}/${suffix}.Sq.dat
Aq=${curr_dir}/${suffix}.Aq.dat
Eq=${curr_dir}/${suffix}.Eq.dat
Df=${curr_dir}/${suffix}.Df.dat
Temp=${curr_dir}/${suffix}.Temp.dat
ZPE=${curr_dir}/${suffix}.ZPE.dat
vol=${curr_dir}/${suffix}.volume.dat

# Write header
for file in $Sq $Aq $Eq $Df $ZPE $vol $Temp
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
echo "Diffus(cm2\/s_\,_1\/s)" >> $Df
echo "ZPE_(kJ/mol/SimBox)" >> $ZPE
echo "volume_(A^3)" >> $vol
echo "temperature_____(K)" >> $Temp

# Grep data
echo "Getting 2PT results.." Sq Aq Eq Df ZPE Temp
for file in $filelist
do
	fnumber=${file#${suffix}.npt.}
	fnumber=${fnumber%.2pt.mol.grps.thermo}

	echo $fnumber `grep Sq_ $file | sed -e "s/${suffix}.npt.//" -e "s/.2pt.mol.grps.thermo://" -e "s/Sq_(J\/mol_K\/SimBox)//"` >> $Sq
	echo $fnumber `grep Aq_ $file | sed -e "s/${suffix}.npt.//" -e "s/.2pt.mol.grps.thermo://" -e "s/Aq__(kJ\/mol\/SimBox)//"` >> $Aq
	echo $fnumber `grep Eq_ $file | sed -e "s/${suffix}.npt.//" -e "s/.2pt.mol.grps.thermo://" -e "s/Eq__(kJ\/mol\/SimBox)//"` >> $Eq
	echo $fnumber `grep Diffus $file | sed -e "s/${suffix}.npt.//" -e "s/.2pt.mol.grps.thermo://" -e "s/Diffus(cm2\/s_\,_1\/s)//"` >> $Df
	echo $fnumber `grep ZPE $file | sed -e "s/${suffix}.npt.//" -e "s/.2pt.mol.grps.thermo://" -e "s/ZPE_(kJ\/mol\/SimBox)//"` >> $ZPE
	echo $fnumber `grep vol $file | sed -e "s/${suffix}.npt.//" -e "s/.2pt.mol.grps.thermo://" -e "s/volume_       (A^3)//"` >> $vol
	echo $fnumber `grep temperature $file | sed -e "s/${suffix}.npt.//" -e "s/.2pt.mol.grps.thermo://" -e "s/temperature_____(K)//"` >> $Temp
done

# get Averaged temperatures -- deprecated
#echo 'Getting Averaged Temperatures from 2PT log..'
#grep "(K)" $outfile > $Temp
#sed -e "s/${suffix}.npt.//" -e "s/.2pt.mol.grps.screen.out:  Temperature ://" -e "s/ (K)//" -i $Temp

# replace space as tab
for file in $Sq $Aq $Eq $Df $ZPE $vol $Temp
do
	sed -e "s/ /\t/g" -i $file
done

echo "Done.. check your *dat file."
echo "* Sq: "$Sq
echo "* Aq: "$Aq
echo "* Eq: "$Eq
echo "* Df: "$Df
echo "* ZPE: "$ZPE
echo "* volume: "$vol
echo "* Temp: "$Temp

### end of 2PT_getResults.sh
