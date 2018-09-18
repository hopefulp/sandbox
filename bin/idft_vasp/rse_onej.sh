fn=$1

mkdir $fn
cp $fn.inp $fn
cd $fn
qchem $fn.inp $fn.out &
