#!/bin/tcsh

if ( $# < 1 ) then
    echo "input vasp directory"
endif

set pd = `pwd`
set bpd = `basename $pd`
echo "/$bpd"

foreach x ($*)
    if ( ! -d $x ) then
	echo "There is no $x directory; then step"
	exit
    endif

    cd $x
	set pd = `pwd`
	set bpd1 = `basename $pd`
 	echo "./$bpd1"
	mkdir tmp
	mv * tmp
	cd tmp
	    set pd = `pwd`
	    set bpd2 = `basename $pd`
 	    echo "./$bpd1/$bpd2"
	    cp POSCAR POTCAR INCAR KPOINTS ..
	    cd ..
	set pd = `pwd`
	set bpd = `basename $pd`
 	echo "./$bpd"
	rm -r tmp
    	cd ..
    set pd = `pwd`
    set bpd = `basename $pd`
    echo "/$bpd"
    echo "initialized $x"
end


