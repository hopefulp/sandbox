### after optimization, do osp (sp after opt) calculation with the same params
### in all the systems.

### another functional

### U correction
forU.sh Job dir1 dir2 start_U interval_U n_points n_TM Mag
	: U_run
		V_mkdir.tch dir1 dir2_u
    		changeline_incar_U.pl U=[LDAUU]
		LDALL=1.0
	: OUTCAR
		dir1 dir1[dummy] start_U interval_U n_points n_TM Mag
	   e.g. OUTCAR Mn_AFM Mn_AFM  1.8 0.0 1 12 AFM
	: gref
read_outcar.pl dir 1 n_tm Mag
	: ARGV[0]=dirname which have OUTCAR
	: ARGV[1]=initial TM index
	: ARGV[2]=number of TM atoms
	: ARGV[3]=Magnetization

### Atom rearrangement for selective dynamics

vasp_1d_arr_ave_dyn.pl inputfile dynamic_option Ni a[bc] ave1 ave2 ave3 
	:
vasp_1d_arr_div_dyn.pl inputfile dynamic_option Ni a[bc] division_limit_1 2 ...
	: by division limit, the group is divided by n+1 group 

vasp_x_arr_ave_dyn.pl inputfile  Ni a[bc] ave1 ave2 ave3 ..
	: group is divided by number of average distances

### Design MOF + Mol in MS
### MSI to POS

extract_natoms.pl poscar "N 10 C 5 etc"
	: direct coordinate
	: extrace mol in mol+MOF system
insert_mol2mof_add.pl bare_mof mol_lattice "N 10 C 5 etc"
	: after "extract_natoms.pl"
	: 


insert_mag.pl INCAR Me nmetal FMAll[FM|AFM|read file] 
	: insert MAGMOM for 3d transition metals
	: modify AFM magnetism
	option: alternative, half

insert_mol2mof.pl mof1 mof2+mol [new_filename] "C 2 H 1"
	: mof1-CONTCAR.MOF74.Me mof2-A.pls 
	: it should have blank in molecule type

vasp_task_loop.sh
	: for file with dirname
vasp_task.sh 
	: cp grep grepa grepf mv qsub rename rm rmr dos_fermi msi2pos.pl ls
	: grep - 3 argument prior to dirname
		{1:-jobname} {2:-filename} {3:-string}

vasp_anal.sh job[bader bader2]
    vasp_chgcar_del.pl dir
	: delete atom symbol line in CHGCAR for bader analysis

V_mkdir.sh old_dir new_dir 0|1|2
	: copy old_dir new_dir
	: copy POSCAR(1) or CONTCAR(2|no arg)
 
V_run.sh job[cont|cont2|dirrun] dir vasp_exe sub[version]
	: dirrun dir vasp_exe
	: cont   dir vasp_exe version
	: cont2  dir vasp_exe version
	:	make sub directory for save magnetization

V_run_mof74.sh
    echo "Usage :: $0  job w_dir  metal mag incar-cont [ini:poscar|dos:o_dir|cont:version]"
    echo "      :: job == rerun ini cont cpdir dos"
    echo "      ::        rerun w_dir    "
    echo "      ::        ini   w_dir metal mag ini ac"
    echo "      ::        cont  w_dir metal mag wav 1 "
    echo "      ::        cpdir n_dir metal mag ini o_dir "
    echo "      :: for job==ini POSCAR from POSCAR.Fe.mol mol= \"zz, ac, py, pa etc\" "
    echo "      :: mag : incar.AFM or FM or NM"
    echo "      :: INCAR continous == ini wav chg"
    exit

msi2pos.pl A.msi Co O C H
	: makes A.pos
./job_mkdir_ff.sh 1
	: read A.pos make A-dir and run
./job_mol.sh
	: make job dir but is supposed to be modified in INCAr, Mag
copy_dir.sh
	: make new dir and copy four files


job_mkdir_co2.sh nnode
    insert_co2opt.pl MOF.pos CONTCAR.6co2.dirt
	: "CONTCAR.6co2.dirt" was obtained from PE calculation already
	: insert 6 co2 to relaxed cell structure (cp $dir/CONTCAR MOF.pos)
job_mkdir_1co2.sh nnode		# initial job for new system
    insert_co2.pl MOF.pos CONTCAR.1co2.dirt  	# insert 1CO2 in the MOF
	: insert 1 CO2 from CONTCAR.1co2.dirt to MOF, MOF.pos
    ext_rm_nco2a.pl POSCAR	# remove 5 CO2 from MOF-CO2 system
	: makes POSCAR.ads1co2
mk_dir1.sh old_dir new_dir Me
	: for just one job
job_mkdir_cont.sh 1
	:
    job_mkdir_cont_sub.sh old_dir new_dir
job_mkdir.sh 1
	: loop for all metals
	: "modify" input directory - dir1=cpo27 dir2=_CO2_done or _cell_done
	: "modify" new directory   - newdir1=gga newdir2=sp or CO2 sp
	: "$1" is overloaded for "changeline.pl" for number of nodes
	: define "sys" "func" "w_type"
	    - sys: cell or CO2, func: gga or lda, w_type: (cell)relax, CO2opt,
	      CO2sp etc
    job_mkdir_sub.sh $old_dir $new_dir $Me $func
		: mkdir
		: copy files
		  - pot file, kp2.monk, poscar, incar-ini file
     	  ~/bin/changeline_tmpl1.pl 'incar.520.CO2.gga' Me
     	  ~/bin/changeline_tmpl.pl 'incar.520.CO2.gga' Me cell[CO2|1CO2]
			: INCAR file should be changed for MAGMOM for each Me
			: modify "$mag_later" for number of atoms in MAGMOM
			: make "t.incar"
		: copy t.incar to $new_dir/INCAR
      	~/bin/changeline.pl pbsfold-idft.csh $new_dir nodes=$1
		: run in "job_mkdir_sub" for exit with error
		: read pbsfold-idft.csh
		: $new_dir is jobname for #PBS -N $new_dir	
		: $1 for number of nodes in #PBS -l nodes=$1:ppn=8
## wrap up 
./task_metal_F.sh



## Make potential
task_pot.sh
	: cp vasp/pot and concatenate to POTCAR
task_gnu.sh Me
	: plot gnu with three figures
### 
task_metal.sh
	: draw band structure for all metals
./dos_multif.sh

	
### DOS run and analysis 
V_run.sh dos old_dir new_dir [pbsfile]
	: copy CHGCAR
	: incar_dos.pl
		: make INCAR.dos and copy into new_dir
incar_dos.pl INCAR-scf
	: make INCAR.dos using SCF-INCAR file

vasp_anal_dos.sh


dos in more than 2 nodes is not working for cpo_1co2 now working
usable file for INCAR: "incar.dos.algo" then "incar.dos.alog.gauss" ->
"incar.dos"
# Dos run and pcharge run also.
job_mkdir_dos.sh nnode
    # for 1 job
    job_mkdir_dos_sub.sh old_dir new_dir nnode dos 	for DOS
    job_mkdir_dos_sub.sh old_dir new_dir nnode ch  $Me 	for pcharge
	: mkdir new
	: copy all 6 files
	: qsub pbs

./task_dos_argdir.sh dir1 dir2 ...
	: dir list is transfered in $@ in shell
	: ../dosall.pl in the dir
	
./dosall.pl [a=1[:6]] [l=]3 [m=]1
    ./dosall1.pl [1=]3 [m=]1 [dir| ] [a=]3:[10]
	: Usage
	    dosall.pl DOSCAR for Tdos.dat and Fermi.dat
	: control Fermi_shift='F' or 'T'
	: if there is dir, it read dir/DOSCAR; if not find in the dir
	: ARGV[0]=DOSCAR for "Tdos.dat" and "Fermi.dat"
	: a=1 or a=1:5 "Ldos_a1.dat"
	: a=1 and l=2 [m=2] Pdos "Ldos_dxz_a1.dat"
	: the first number for l(angular momentum) can be # or l=#
	: the second is same    just # or m=#
	: if you use delimiter "=" or ":", the order of arg.'s are not
	: sum pdos for atoms a=1:6 or 1:6 or just 1: 2: 3: etc
	: atom index start from 1 and end at natom [! 0, natom-1 ]
	: if just number - first l and the second m
	: if alphabet, it's dir name and load $dir/DOSCAR, default= ./DOSCAR
	: if a= or num: it's the chosen atoms for series n:m or single n:
   ../gnu_multi.sh Fermi.dat dos1.dat dos2.dat ... 

./task_metal_dos.sh | grep Sum | awk '{print $4}' | perl -ne 'BEGIN {$i=0;}
    chomp($_); print "$_\t"; if($i%2==1) { print "\n";} $i++;'
### GNU plot
task_gnu_1me.sh Me
	: 
    "./gnu_multi.sh" $1 $2 $3 $4 
	: just file name with the number of files should be $#
	: Title is given by dir=`basename pwd`
	: first is black for Fermi.dat
./gnu_for.sh
    ./gnu2.sh dir1 dir2
	: plot dir1/dos.dat dir2/dos.dat
    ./gnu2f_pdos.sh dir1 dir2 fname

./gnu.sh [$1]
	: plot "Tdos.dat"
	: if $1=dir plot $1/Tdos.dat
./gnu1f.sh $1 [$2]
	: if not $2, plot $1=filename
	: if $2, plot $2/$1=filename
	:
### partial charge analysis 
get_procar.pl is now revising
get_Procar.pl PROCAR e=Emin:Emax "atom lists"
	: 2014 new edition
	: print k, band when pdos is over a certain value
	: ARGV[1]= -3.5 or  e=-3.5 or e=-4.0:-5.0
	: output: energy kpoint band pdos[1st atom|Me] pdos[2nd atom|CO2]
job_pch_density.sh old_dir new_dir nnodes
	: old_dir is pre-scf dir
	: make new dir for partial charge density calculation prom PROCAR
job_mkdir_pcharge.sh
### Bader charge analysis
task_metal_mv.sh
  ./vasp_ext_meco2.pl Me.xyz
	: make qchem input file for Me-CO2 electrostatic interaction
	: Me_co2.mol with external chg of Me and Me__co2.mol for only co2
  ./vasp_ext_6meco2.pl Me.xyz
	: for 6CO2+MOF system, first make "Me.xyz", "metal_pcharges.dat"
	: read "metal_pcharge.dat" for charge of 6 metals
	: select 6 charges using metal index in ARG[0]
1. pos2xyz_park.pl $dir/CONTCAR
	: CONTCAR.xyz as xyz format in $dir
	: mv $dir/CONTCAR.xyz "Me.xyz"
    pos2xyz_park.pl $dir/CONTCAR 28 30 34
	   : to fix up the weird coordinate of oxygen
2 get_pcharge.pl ACF.dat [index] 1 2 3 4 5 6 9 [atom name] Mg Mg Mg Mg Mg Mg
	: to get partial charge + or - for electrostatic calculation
	: save > "metal_pcharge.pjh" for 12 metals (it should be same as script
	    ./task_metal_pcahrges.sh > metal_pcharges.dat
    #2.1 get_vcharge.pl ACF.dat 25 26 51 [atom_list]
    #	: sum all the valence charges of listed atoms
    #  get_vcharge.pl	ACF.dat
    #	: average of charges for several atoms (from 1 to 6 in atom index)
# electrostatic energy by coordinates and partial charge
3. coulomb_mof72.pl CONTCAR.xyz ACF.dat Ca 
	: get partial charges from ACF.dat <IN1>
	: get distance from coordinates in POSCAR <IN2>
	: $Me to get referenc charges in the %hash
xyz_dist.pl CONTCAR.xyz
	: to measure distance in 6 pairs	
xyz_angle.pl
	: to measure angle of CO2 in 6 molecules.
task_metal_ene.sh 
	: make ienergy.dat
    anal_energy.pl ienergy.dat Me
	: write energy Me
get_lattice.pl CONTCAR option[0|any]
	: get volume lattice constants
	: with option==1, write in 1 row
grep volume cpo27Sc_CO2_r/OUTCAR
get_fermi.pl dir1 dir2
	: it reads dir1/DOSCAR; DOSCAR is read from "dosall1.pl" also
	: if there are two arguments, it reads both and write them in a line
######### split dos peak
ldos_int.pl fi ref_sep
	: integration of dos peak
	: everage dos into the seperate single peak
ldos_int2.pl fi1 fi2 ref_sep [Me|me|CO2]
	: compare two dosfile fi1 and fi2
	: reference seperation - seperate dos peak if the values are lower
ldos_int_2fermi.pl cpo27_1co2_V_dos.extg/Ldos_a25-51.dat 1.e-5 V [ads|des]
	: input ldos, don't use 0 for reference dos
	: option for "ads" and "des" to see the Fermi level
	: t_dos to the fermi level and t_anti for upper level dos
	    t_dos integrates all the dos than reference in ldos.dat 
	: sys==
### CO2 extract for quadrupole interaction in vasp and qchem
1. Conversion between direct and cartesian coordintes: CONTCAR(d) to POSCAR(c)
  cont2pos_c.pl CONTCAR(d)
	: makes POSCAR(c)
  con22pos.pl POSCAR(c)
	: makes CONTCAR(d)
    cf. posd2c.pl

2. Extract molecule from POSCAR or CONTCAR
2.1 based on atom index (from 0)
extract_index.pl POSCAR index_file charge_file
e.g. extract_index.pl Zn.pos Zn12CO2.index Zn.chg
	: arrange in the order of @atom_extract=(M O C O)
	: write mol file
	: calculate distance between atoms
2.2 from anchored atom
extract_anchor.pl poscar metal_index atom_series
	: extract_anchor.pl Co-chain1.pos "6 7 8 9 10 11" "Co O C O"

2. ./extract_co2.pl CONTCAR D
        : "extract_cartco2_pore.pl" was changed into "extract_co2.pl"
        : to extract 6 CO2 from MOF+6CO2 as direct format with the same cell
        : make "CONTCAR.6co2.dirt"
	: this is not working if C and O has specific order in POSCAR
2. ./pos2cart.pl CONTCAR.6co2.dirt
        : in poscar change direct to cartesian with the same cell
        : make "CONTCAR.cart"
        : it was "posd2c.pl"
3. ./shift_6co2_poscart.pl CONTCAR.cart  1.5
        : $cell_expansion = 1.5 while atoms are fixed in that position in
Cartesian format
        : make "CONTCAR.shift.1.5"
        : molecules are gathered in the pore out of unitcell
  3.1 ./shift_co2_poscart.pl CO2.cart  1.0
        : to gether single molecule in the same sites
4. ./poscart2xyzmol.pl CONTCAR.shift
        : xyz format where molecules are gathered in the xyz file to be used
for BSSE correction
5. xyz2mol.pl 6CO2.xyz
        : make "6CO2.mol" $molecule block in the qchem input file
6. add method part by hand
        : to make "6CO2.inp" qchem input file
7.
        : script for BSSE for 6 cases of 6 CO2

2. extract CO2 coordinate from CONTCAR.cart to CONTCAR.6co2.pos
        "extract_co2.pl"
./extract_nco2a.pl CONTCAR 1 D
	: extract 1 CO2 as direct coord from CONTCAR
# for CO2 calculation
shift_outofcell_poscart.pl CONTCAR.6co2.cart
        : "CONTCAR.6co2.cart" poscar in cartesian coordinates of only 6 CO2
        : "CONTCAR.shift" 6CO2 position in the pore of MOF
### extract CO2 and copy into the MOF file for relaxation of CO2 in the MOF
./insert_co2opt.pl MOF.pos 6CO2.pos
        : write combination of bare-MOF and 6CO2(made from MgMOF-6CO2)
	: to include selective dynamics
./insert_co2.pl MOF.pos CO2.pos
	: at the moment the order in CO2.pos [O C] should be same as MOF.pos
./insert_co2b.pl


### Clean VASP directory
vasp_cleanf.sh dir option[all|W|CW]
	: clean vasp directory upto 2 sub-dirs
