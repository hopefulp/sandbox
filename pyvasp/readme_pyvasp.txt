### Env
sandbox/pypbs/Env_msg.py

### lib connection
make_ini -> myvasp : make_incar (for only read)
make_incar -> myvasp : make_incar (read or write)

### how to run
vasp_drun.py dir [dir [ ..] ] -r
    : run outside of directory
    : all for all directory runs
    : sed to replace jobname in pbs file
    : if dir == 'all': runs all the directories


#### Run
v_prepare.py [new,cont,cont4] new_dir -o old_dir
    : prepare new_dir with INCAR, POSCAR, POTCAR, KPOINTS
v_cp_pbs.sh dir_job type[std,ncl]
    : prepare pbs file from /qcfs/joonho/bin/template/vsp_all.pbs
    qsub vsp_all.pbs
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
##### Clean directory
vclean_rec.py
    : clean directories recursively

