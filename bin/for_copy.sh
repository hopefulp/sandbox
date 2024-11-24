#!/bin/bash

cwd=$(pwd)
dir={cwd:-$cwd}

base1=d201559
base=${1:-$base1}

i=a
f=j
f1="kp.gamma"
f2=KPOINTS
a=${2:-$i}
b=${3:-$f}
fname1=${2:-$f1}
fname2=${3:-$f2}
for x in {a..j}
    do
        dir=$base$x
        #echo $dir
    	#echo cp $dir/kp.gamma $dir/KPOINTS
    	echo cp $dir/${fname1} $dir/${fname2}
    done

