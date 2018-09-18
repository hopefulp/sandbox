#!/bin/bash

dir1=cpo27
dir2=_CO2_r1_ch
for Me in `cat metal.dir`
  do
     dir=$dir1$Me$dir2
#     echo  $Me
#     ./get_vcharge.pl $dir/ACF.dat 1 #| perl -ne 'chomp($_); print "$_\t";' ; ./get_vcharge.pl $dir/ACF.dat 25 26 51
#    ./get_vcharge.pl $dir/ACF.dat 25 26 51 #| perl -ne 'chomp($_); print "$_\t";' ; ./get_vcharge.pl $dir/ACF.dat 25 26 51

     ./get_pcharge.pl $dir/ACF.dat 1 2 3 4 5 6 $Me $Me $Me $Me $Me $Me
#     ./get_pcharge.pl $dir/ACF.dat  1     7     9    11    12    13 $Me O O O O O
#    ./get_pcharge.pl $dir/ACF.dat  2     7    10    12    14    15 $Me O O O O O
#    ./get_pcharge.pl $dir/ACF.dat  3     8     9    10    13    15 $Me O O O O O
#    ./get_pcharge.pl $dir/ACF.dat  4    16    18    20    21    22 $Me O O O O O
#    ./get_pcharge.pl $dir/ACF.dat  5    16    19    21    23    24 $Me O O O O O
#    ./get_pcharge.pl $dir/ACF.dat  6    17    18    19    22    24 $Me O O O O O
  done



