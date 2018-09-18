#!/bin/bash

#dir1=DOScpo27
#dir2a=_cell_anal
#dir2b=_CO2_anal

#dir1a=cpo27_1co2_
#dir1b=_new
#dir1b=_new.d
#dir2a=cpo27_1co2_ads_
#dir2b=_mag
#dir2b=_dos.extg

#dir1=../vasp_relax/cpo27
#dir2a=_cell_done
#dir2b=_CO2_done
dir1=cpo27
dir2=_CO2_rlda

file=conv.chg
for Me in `cat metal.dir`
  do
     echo  $Me
#    dir_cell=$dir1a$Me
#    dir_co2=$dir1$Me$dir2b
#     file=$dir1$Me$dir2.out
#     dir_des=$dir1a$Me$dir1b
#     dir_ads=$dir2a$Me$dir2b

     dir=$dir1$Me$dir2
     ls -l $dir/WAVECAR
#     dir=$dir1$Me$dir2a
#      head -7 $dir_des/conv.chg | awk 'BEGIN {i=0;} {if(i!=0)  print $3 ; i++; }' | perl -ne 'END {print "\n"} chomp($_); print "$_\t";'
#      head -7 $dir_des/conv.chg | awk 'BEGIN {i=0;} {if(i!=0)  print $3 ; i++; }' | perl -ne 'END {print "\n"} chomp($_); print "$_\t";'
#      head -2 $dir_des/conv.chg | awk 'BEGIN {i=0;} {if(i!=0)  print $2 ; i++; }' ; head -2 $dir_ads/conv.chg | awk 'BEGIN {i=0;} {if(i!=0)  print $2 ; i++; }' | perl -ne 'END {print "\n"} chomp($_); print "$_\t";'
#     dir=$dir1$Me$dir2b
#     head  $dir/POSCAR
#     head -2 $dir/POTCAR
#    cd $dir_des
#    	/qcfs/biduri/scripts/vasp/convasp -chgint CHGCAR > conv.chg &
#    cd $dir_ads
#      ../get_procar.pl PROCAR -8 -5
#	/qcfs/biduri/scripts/vasp/convasp -chgint CHGCAR > conv.chg &
#    cd ..
#    ./get_vcharge.pl $dir_des/ACF.dat 25 26 51 | perl -ne 'chomp($_); print "$_\t";' ; ./get_vcharge.pl $dir_ads/ACF.dat 25 26 51
#    ./get_vcharge.pl $dir_des/ACF.dat 26 | perl -ne 'chomp($_); print "$_\t";' ; ./get_vcharge.pl $dir_des/ACF.dat 25 | perl -ne 'chomp($_); print "$_\t";'  ; ./get_vcharge.pl $dir_des/ACF.dat 51 | perl -ne 'chomp($_); print "$_\t";' ;  ./get_vcharge.pl $dir_ads/ACF.dat 26 | perl -ne 'chomp($_); print "$_\t";' ; ./get_vcharge.pl $dir_ads/ACF.dat 25 | perl -ne 'chomp($_); print "$_\t";'  ; ./get_vcharge.pl $dir_ads/ACF.dat 51

#      grep F=  $file
#      dir=$dir_ads
#     more $dir_des/$file | perl -ne '@list=split(/\s+/,$_); if($list[0] eq ""){ shift(@list);}  if($list[0]==26) {print $list[1],"\t";}'; more $dir_des/$file | perl -ne '@list=split(/\s+/,$_); if($list[0] eq ""){ shift(@list);}  if($list[0]==25) {print $list[1],"\t";}' ; more $dir_des/$file | perl -ne '@list=split(/\s+/,$_); if($list[0] eq ""){ shift(@list);}  if($list[0]==51) {print $list[1],"\n";}'
#     more $dir/$file | perl -ne '@list=split(/\s+/,$_); if($list[0] eq ""){ shift(@list);}  if($list[0]==26) {print $list[1],"\t";}'; more $dir/$file | perl -ne '@list=split(/\s+/,$_); if($list[0] eq ""){ shift(@list);}  if($list[0]==25) {print $list[1],"\t";}' ; more $dir/$file | perl -ne '@list=split(/\s+/,$_); if($list[0] eq ""){ shift(@list);}  if($list[0]==51) {print $list[1],"\n";}'
#      dir=$dir2a$Me
#      grep ZVAL $dir/POTCAR   

	
  done



