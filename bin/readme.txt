##### VASP ##########################################################
##### Analysis
dosall.pl
	: old type of DOSCAR
dosall_newf.pl
	: new type of DOSCAR
##### Change file
cp_vasp_run.sh dir1 dir2 ...
	: locate old home directory
	: modify new home directory

changeline_pbs.pl pbsfold-idft.csh new_dir nodes=#
	: make new pbs script with # of nodes and make work-directory
	: read pbsfold-idft.csh and write t.csh
	: t.csh include new_dir name (job directory) and number of nodes
changeline_pbs_kdft.pl pbsfold-kdft.csh  new_dir nnodes g[1|2]
	:
##### Converting file: msi2pos
msi2pos.pl A.msi [atom series]
	: make A.pos
	: stdout is different from A.pos
	: no md T/F
msi2pos_mdnew.pl
	: extend MD T/F, default all F F F
msi2pos_molec.pl A.msi [atom series]
	: shift molecule from origine to center by shifting half of lattice
constants
##### Electric potential
Acc-mol.sh M
    	: charge file should be M.chg
	: change Accelry mol file to mol w/wo external charge
	: "file.mol" and "file_ext_charges.mol"


##### GNU #####################


##### Q-CHEM ########################################################
### PSI : Parallel
#   PRIMP2: Parallel for SP but serial for Opt
# single job

# parallel single job
qsub pbs_qchem1_psii[kdft].sh
	: modify jobname[.in]
	: correct rem
pbs_qchem1_in_psi.sh
	:modify jobname[.in]

### KDFT : Serial
# single job
qc_ser_kdft.sh(run_sj_kdft.sh) mol [method] [job]
	: $1="mol".in
	: $2=[rimp2]
	: $3=opt
qc_ser_bsse_kdft.sh(run_bsse_kdft.sh) fn[.in] nA nB
	: make fn.in for SP
	: call "bsse_2mol.pl"

##### BSSE job program
qbsse_psi.sh mol[.in] nA nB pbs_qchem1_in_psi.sh [16]  [5]
	: basename of .in file
	: make new dir "mol"
	: /qcfs/joonho/bin/"bsse_2mol_psi.pl" $mol $na $nb
	: if $jobnumber==4 (in case a sp was done) , remove a..in (not sp for original mol)
	: /qcfs/joonho/bin/"changeline_pbs_bsse.pl" $pbs_file $filename $ppn 

bsse_2mol.pl	fname	na	nb
	: fname without extension (e.g. co2.inp => co2)
	: na is the number of atoms of the first molecule
	: nb is the number of atoms of the second molecule
	: output afname, bfname, cfname, dfname, efname
	: bsse correction = a - (b - c + d - e)
##### ANALYSIS
find_key.pl a.out 2
	: find 2nd keyword FREQ] and show 10 lines
qgeo_mo.pl  file.out [any]
	: get last geometry in MOLDEN format, find "GEOMETRIES"
	: with option [any], get last geometry in "opt" in multiple job outfile
qnon_geo.pl file.out 
	: select last unconverged geometry, find last "Optimization Cycle"
qgeo_o.pl file
	: old file
	: get last geometry in "opt"


task_case_qchem1.sh bsse dir_name mp2

grep "RIMP2" *out | grep "total" | awk '{print $6}'
grep "Total energy" phenol_co2b/*out | awk '{ print $10 }' | perl -ne '{END
{print "\n";} chomp($_); print "$_\t"; }' | awk '{ print $1-$2+$3-$4 }'

get_basis.pl Mg
	: get basis of cc-pvtz from "/home/joonho/basis/cc-pvtz.bas"
	: standard output

#### Qchem job file
get_freq_all.pl qchem.out[jobtype=freq]
	: makes qchem.freq
./qcbatch-1f.sh .outfile
        : for a "continuous job" from Opt ot sp of 1 optimized file
        : get mol file with (-augdz.mol)
        : make in file with rem.pbed_spaugth
        : run pbs with "pbs_qchem1.sh"

./qcbatch-tailor.sh | sh
        : change file name with tailing
./qcbatch_out.sh name
        : get optimized geometry from *.out file and add tail name with "arg"
./qcbatch-run.sh 
        : predefine pbsfile and prepare rem file
	: run all the mol file
./qc1f.sh file.mol
	: modify rem
./qc1f_getgeo.sh outfile base
	: for 1 outfile


$qcbin/changeline_pbs_kdft_QC.pl pbs_qchem_inp.sh azocop2 b3lyp_opt g2
	: "pbs_qchem_inp.sh" pbs input file
	: jobname = input filename = mol file name
	: change "rem" 
	: class in cluster
	: does not change mol name
changeline_pbs_kdft.pl


get_remgeo.pl geo=a.inp rem=b.inp inp=c.inp
	: read two files and make c.inp
	: e.g. get_remgeo.pl geo=geo.inp rem=rem.inp inp=model.inp
make_mol_atomchange.pl datfile molfile model_name
	: datfile has two column for metal name and multiplicity
	: molfile should have qchem format for molecule
	: check charge
qchem_rem_arg.pl	a=b c=d ...
	: read "dir.txt" which include input filename one by one per line
qchem_change_rem.pl filename jobname=new a=b c=d
	: read input filename
	: make input file with the same geometry
qchem_rem_arg_1f_readgeo.pl  filename jobname=new a=b c=d
	: actually it doesnot need to read an input file
	: outfile=qchem input file with read molecular file
	:
qchem_rem_arg_1f_2job.pl  filename a=b c=d
	: read input filename and add consecutive job
	: call  "qchem_addjob.pl" to add one more job
	: new job is listed in "qchem_addjob.pl"

common	: default
		sp hf rimp2/cc-pvtz
	: modify default or use command line argument by variable=option
		e.g. EXCHANGE=l3lyp

qchem_addjob.pl a.inp
	: add new job in a.inp
	: default is the same as above "common"
qchem_finalgeo.pl name.out
	: read qchem geometry optimization out file
	: std output xyz format
	: mol.inp in qchem format

##### SCAN: orientation manipulation
plane_vector.pl 
	: input	: number of line for three atom for plane
	: output: directional unit vector
	: theory: Cramer theory to get coefficients of plane equation,
Ax+By+Cz+D=0

##### SHELL
cp_dir_sub2 . obj_dir filename
	: at present directory find filename; then copy with dirname in
obj_dir
cp_tree_dir2_for.sh [$1]
	: copy trees from old directory of $1 to new directory of pwd
	: pwd is level 2 parent directory
	: scans "POSCAR CONTCAR KPOINTS INCAR *out"
	: cp_tree_dir1.sh is for level 1 parent directory
	
