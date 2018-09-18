###### 	MOFbp-74 + CO2
######	DOSCAR
"Find overlap energy"
### run dos
V_run.sh dos old_dir new_dir pbs-file[check exe in advance]
### extract dos
vasp_anal_dos.sh dir1 dir2 ...
dosall.pl DOSCAR
	: for Tdos and Fermi.dat
	: for CO2, 9, 49
### draw dos data in the directory
gnu_multi.sh Pdos_d_a11.dat Ldos_a49-135.dat Fermi.dat

gnu_2f.sh Tdos.dat Fermi.dat
gnu_4f.sh Tdos.dat Ldos_a1-6.dat Ldos_a10-79.dat Fermi.dat
# y = [][25]
gnu_3f.sh Ldos_a1-6.dat Ldos_a10-79.dat Fermi.dat
# y = [][7]
$qcvin/vasp_anal_gnu3f.sh Mg1u-1Ace-dos Ldos_a10.dat Ldos_a1.dat Ldos_a2.dat Ldos_a3.dat Ldos_a4.dat Ldos_a5.dat Ldos_a6.dat

####### 	PROCAR
"Find partial charge at the overlap energy
get_Procar.pl e=-5 "3  6 10"
get_Procar.pl e=-5 "6 10"

### Ni
get_Procar.pl d=Ni2u-1co2-dos e=-9:-7 "1 2 3 4 5 6 7 8 9 10 11 12 49"
   : atom list: (13 atoms) 1  2  3  4  5  6  7  8  9  10  11  12  49
   : density criteria: 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.015  0.  0.015
IBAND = 210 219 220 221 222 223 224 225 226 227 228 229
### Mn
get_Procar.pl d=Mn2u-1co2-dos e=-9:-7 "1 2 3 4 5 6 7 8 9 10 11 12 49"
IBAND = 211
### Co
get_Procar.pl d=Co2u-1co2-dos e=-8:-6.5 "1 2 3 4 5 6 7 8 9 10 11 12 49"
IBAND = 217 218 219 220 221 226 227 228 229 230 233 239
### Mg
get_Procar.pl d=Mg2u-1co2-dos e=-8:-7 "1 2 3 4 5 6 7 8 9 10 11 12 49"
### Zn
get_Procar.pl d=Zn2u-1co2-dos e=-8:-7 "1 2 3 4 5 6 7 8 9 10 11 12 49"
IBAND = 265 271 275

last line is sum of all the k and b normalized by num of K-points

####### Partial CHG
V_run.sh pchg Co1u-1Ace-Nupd-dos Co1u-1Ace-Nupd-band
task_mv.sh dir

##### MOFbp-74 + CO2
### U correction
forU.sh job dir1 dir2 start_U interval_U Num_U
    	: job= U_run OUTCAR gref
read_outcar.pl dir ini_atom fin_atom
	: to see magnetization

### Electrostatic energy calculation
extract_index poscar index_file chg_file
	: index_file is made manually by visual program

extract_anchor poscar "anchor [metal] atom index" "Co O C O[series of atoms]"

