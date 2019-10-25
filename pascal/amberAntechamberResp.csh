#!/bin/tcsh
#!/bin/csh

antechamber -i ttha_protonated.mol2 -fi mol2 -o ttha_protonated.ac -fo ac -nc -5 -m 1 -rn TTH -rn ttha -at gaff -s 2 -pf y
respgen -i ttha_protonated.ac -o ttha_protonated.respin1 -f resp1
respgen -i ttha_protonated.ac -o ttha_protonated.respin2 -f resp2
resp -O -i ttha_protonated.respin1 -o ttha_protonated.respout1 -e ttha_protonated.resp -t qout_stage1 
resp -O -i ttha_protonated.respin2 -o ttha_protonated.respout2 -e ttha_protonated.resp -q qout_stage1 -t qout_stage2
antechamber -i ttha_protonated.ac -fi ac -o ttha_protonated.prepi -fo prepi -c rc -cf qout_stage2 -pf y -s 2 -at gaff -rf TTHA -rn tth
