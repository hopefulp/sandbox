#!/bin/tcsh
#!/bin/csh

#foreach f (../../ffs/orig/*.par)
if ($#argv < 4) then
    echo "usage: $0 [forcefields] [bgf files] [lammps file] [save name]"
    exit(1)
endif

set fflist = ($1)
set bgflist = ($2)
set lmp = ($3)
set save = ($4)

if ($#fflist == 0) then
    echo "Error: No valid forcefiles found while searching \"${fflist}\""
    exit(1)
endif
if ($#bgflist == 0) then
    echo "Error: No valid bgf files found while searching \"${bgflist}\""
    exit(1)
endif

foreach f ($fflist)
    set fle = `basename $f .par`
    set dir = `dirname $f`
    set input_file = ${dir}/in.${fle}
    set fle = `basename $fle _opt`

    foreach b ($bgflist)
	set molname = `basename $b .bgf`
        echo running calc on $molname using $fle using $input_file
	/home/yjn1818/scripts/createLammpsInput.pl -f ${f} -b $b -s ${molname}_${fle} > /dev/null
	rm -fr in.${molname}_${fle}* ${molname}_${fle}_lammps.script
	set datafile = data.${molname}_${fle}
	$lmp -in $input_file -log ${molname}_${fle}_log.lammps -var datafile $datafile > /dev/null
        rm -fr $datafile
    end
end
/home/yjn1818/scripts/parseLammpsVDWResults.pl -l "*log.lammps" -s ../${save}
