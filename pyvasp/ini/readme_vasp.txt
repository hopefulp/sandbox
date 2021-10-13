####################################################################################
##################### VASP ##########################################################

########## RUN
V_run_ini.sh POT INC KPOINT dirsuff
   "USE": vrun_poscar.pl
	: vrun_incar.pl
	: changeline_pbs_$HOST.pl

V_run_cont.sh o_dir n_dir [cont]
	: default: continous job (1 1 2)
V_run_again.sh dir # [ of version ]
	: all the same directory, run again for geometry optimization


changeline_pbs.pl pbsfold-idft.csh new_dir nodes=#
        : make new pbs script with # of nodes and make work-directory
        : read pbsfold-idft.csh and write t.csh
        : t.csh include new_dir name (job directory) and number of nodes
changeline_pbs_kdft.pl pbsfold-kdft.csh  new_dir nnodes g[1|2]





##### Analysis
dosall.pl
        : old type of DOSCAR
dosall_newf.pl
        : new type of DOSCAR


##### Change file
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


