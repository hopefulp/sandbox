#!/bin/tcsh
#!/bin/csh

set seed = `awk 'BEGIN {srand();myRandomSeed=int(rand()*systime());print myRandomSeed}'`
set slen = `expr length $seed`
if($slen > 8) then
    set seed = `echo $seed | cut -c1-8`
endif
echo $seed
