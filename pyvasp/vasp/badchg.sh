#!/bin/tcsh
chgsum.pl AECCAR0 AECCAR2
bader -c voronoi CHGCAR -ref CHGCAR_sum > ACF.voronoi
#bader CHGCAR -ref CHGCAR_sum
