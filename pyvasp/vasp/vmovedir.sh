#!/bin/tcsh
set src = `basename $1 /`
set tgt = `basename $2 /`
mv $src $tgt
cd $tgt
if (-e $src.1st.out) then
	mv $src.1st.out $tgt.1st.out
endif
mv $src.out $tgt.out
