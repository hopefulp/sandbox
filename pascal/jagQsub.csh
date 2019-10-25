#!/bin/tcsh
#!/bin/csh

if ($#argv < 1) then
echo "usage: $0 [jaguar prefix| input file]"
exit(1)
endif

set mol = `/home/yjn1818/scripts/getFileAbsPath.pl $argv[1]`
set molname = `basename $mol .in`
set inputfile = ${molname}.in
set maefile = ${molname}.mae
set curr_dir = `echo $PWD`
set jagType = "run"

if ($#argv > 1) then
    set jagType = $argv[2]
endif

if (((-e $inputfile) == 0) || ((-e $maefile) == 0)) then
echo "ERROR: Cannot find $inputfile or $maefile"
exit(1)
endif

cat > ${molname}_jaguar.script <<EOF
#PBS -l nodes=1,walltime=960:00:00
#PBS -q workq
#PBS -e ${curr_dir}/${molname}_jag.err
#PBS -o ${curr_dir}/${molname}_jag.out
#PBS -N ${molname}_jag
#!/bin/csh
#!/bin/tcsh

echo "Nodes:"
cat \$PBS_NODEFILE
echo "Jaguar run of ${molname}"
echo

# For jaguar on cluster uses a special license server
setenv LM_LICENSE_FILE @10.0.0.1

# borg's nodes aren't in jaguar.hosts file so define temp dir
setenv JAGUAR_TEMP /temp1/$user

cd $curr_dir

jaguar $jagType $inputfile -WAIT || goto error

echo " No errors detected - All Tasks Completed "
exit(0)

error:
echo "Error occurred. See ${molname}.out"
exit(1)

EOF
