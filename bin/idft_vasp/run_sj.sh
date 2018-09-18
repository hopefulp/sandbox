fdname=$1

mkdir $fdname
cp $fdname.inp $fdname
cd $fdname 
qchem $fdname.inp $fdname.out &

