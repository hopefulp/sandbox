#!/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_all, MyClass, dir_classify_n, whereami


amp_collection = {
    'file_conv':    'im2extxyz.py'  ,
    'amp_run':      'amp_ene.py'    ,
    'amp_valid':    'amp_validation.sh',
    'amp_scan':     'amp_loop.sh'}

models= {
    'ethylene': 'Ethylene.extxyz',
    'Diss_CHO': 'Diss_H2COH.extxyz',
    'water'   : 'water128.extxyz'}
### these are included as key in globals()
make    = MyClass('make')
run     = MyClass('run')
clean   = MyClass('clean')
poscar  = MyClass('poscar')
incar   = MyClass('incar')
potcar  = MyClass('potcar')
dosband = MyClass('dosband')
procar  = MyClass('procar')
outcar  = MyClass('outcar')
md      = MyClass('md')
ase     = MyClass('ase')
charge  = MyClass('charge')
convert = MyClass('convert')
analysis = MyClass('analysis')
make.vas_make_ini   ="==================== VASP INI START =======================================\
                    \n\t===== Make VASP initial directory ====\
                    \n\t===== Prepare POSCAR, POTCAR, KPOINTS, INCAR\
                    \n\tOPTIONS:\
                    \n\t    -s  poscar\
                    \n\t\tPOSCAR.dirname <- convention\
                    \n\t    -j/sj [sp,opt,kp,fake,md,mdnve,nnff,...]/[sp,kp,cool,heat,quench]\
                    \n\t\tkp moved to subjob: kp/sp(default) -> sp/kp\
                    \n\t\tfake:: -sj for -j is required to use INCAR.job, KPOINTS.job\
                    \n\t    -i  can designate INCAR file explicitly\
                    \n\t    -io can modify the INCAR by k-v pairs\
                    \n\t    -ra == -r a: yes for all\
                    \n\t    -r  on [mkdir,qsub] overwrite/nosub\
                    \n\t\tk kill but show poscars and dirs\
                    \n\t    -o  qopt to use different qscript\
                    \n\t\t long queue=long in KISTI\
                    \n\t    -al stop questioning except -s poscar\
                    \n\t    -hpp pseudo hydrogen in /TGM/Apps/VASP/POTCAR\
                    \n\t    Group k-sampling\
                    \n\t\t-kd --kdim [1,2,3], 1: 1 1 k; 2: k[0] k[1] 1;  3: k[0], [1], [2]\
                    \n\t\t-kps 4 4 6 4 8 4 10 4 for -kdim 2\
                    \n\tPOTCAR: \
                    \n\t    will be made by 'genpotcar.py -pp pbe' after cd and reading POSCAR\
                    \n\tKPOINTS:\
                    \n\t    constructed by raw: if fail M, try G\
                    \n\tINCAR:\
                    \n\t    need to be prepared in advance\
                    \n\te.g. (KISTI)\
                    \n\t    kpy vas_make_ini.py -s POSCAR.IntMoS2 -j opt\
                    \n\t    kpy vas_make_ini.py -s POSCAR.IntHfSe2md -j md -al -o long (queue long option)\
                    \n\t    vas_make_ini.py -s POSCAR.MoS2T2Kb -j md -io TEBEG 2000 TEEND 2000 -k g -d MoS2nvtT2K\
                    \n\te.g. (Pseudo Hydrogen)\
                    \n\t    vas_make_ini.py -s POSCAR.sc11sHHf -j sp -al -hpp .66 1.33\
                    \n\t    kpy vas_make_ini.py -s POSCAR.HfOHSe2MLf1 -j opt -d t1736a -hpp .5 .66 1.5 1.33\
                    \n\te.g. (kpoint sampling in Pt)\
                    \n\t    vas_make_ini.py -s POSCAR -kdim 2 -kps 4 4 6 4 8 4 -x 2 -N 1 -n 4\
                    \n\te.g. (fake job for kisti)\
                    \n\t    vas_make_ini.py -j fake -sj opt -d dname -nd ndirs\
                    \n\te.g. (NEB)\
                    \n\t    qsub -N HfSe2mO2sglnebafix -l select=32:ncpus=40:mpiprocs=40:ompthreads=1 $SB/pypbs/pbs_vasp_kisti_skl.sh\
                    \n\t    kpy vas_make_ini.py -s POSCAR.HfSe2mO2tplini -j neb -i INCAR.nebmag -al -ja POSCAR.HfSe2mO2tplfin\
                    \n\t===== D2D w. NEB  ===================\
                    "
make.vas_make_d2d   =" Make Vasp dir from the existing old dir\
                    \n\tvas_make_d2d.py old_dir new_dir job [options]\
                    \n\tjob = {cont, md, neb,... post process}\
                    \n\toptions:\
                    \n\t    -j [ini,cont,neb,nebcont]\
                    \n\t    -o qsub command options\
                    \n\t\tlong: long queue\
                    \n\t\tmem: memory issue - run half cpu to double memory\
                    \n\tNEB::\
                    \n\t    kpy vas_make_d2d.py HfSe2mO2sglneba HfSe2mO2sglnebacini -j neb[cont] -o long (neb: copy, cont:for time limit)\
                    "
make.vas_make_incar ="\n\t    If not 'incar.key', make it, check it and modify it before run this again\
                    \n\t    based on 'incar.key', make INCAR\
                    \n\
                    "
make.vas_make_cont = "===================== VASP CONTINUE =====================================\
                    \n\t-d dir_list -j job -i incar_option -o option_poscar\
                    \n\tchange of d2d to make_cont\
                    \n\toptions:\
                    \n\t    -d old dir list\
                    \n\t    -n new dir list, (default=olddircont)\
                    \n\t    -j default='cont' \
                    \n\t    -j [sp,chg,opt,copt,vdw,noD,mag,kisti] changes only INCAR vasp_job_incar()\
                    \n\t    -j [band,dos,zpe] changes INCAR, KPOINTS,POSCAR in vasp_jobs()\
                    \n\t    -i to modify INCAR\
                    \n\t    -k to modify KPOINTS\
                    \n\t    -io INCAR option with [key, value] list: INCAR modifi\
                    \n\t**  INCAR modification by job: vdw,noD, opt,copt, mag, kisti\
                    \n\t\t    incar: no preexist values, given from cli\
                    \n\t**  INCAR: import libincar for modify\
                    \n\t    -i designate explicit INCAR\
                    \n\t    Usage:\
                    \n\t\tvas_make_cont.py -d IntHfSe2sc34AmC -n IntHfSe2sc34AmCcont -ia EDIFF 1E-5 -x 5 -N 6\
                    \n\t\tkpy vas_make_cont.py -d HfO12md -n HfO12md2 -j ini -io TEBEG 500 TEEND 2000 NSW 2000 -k g\
                    \n\tJobs Detail:\
                    \n\t    ini, cont: (def ini)\
                    \n\t\tcopy odir [INCAR, KPOINTS, POTCAR] w. given -s POSCAR\
                    \n\t\t     if not -s POSCAR, use odir/POSCAR(ini) or CONTCAR(cont)\
                    \n\t\t<eg> -d odir -j ini -s POSCAR.newmodel\
                    \n\t\t    : default ndir is newmodel\
                    \n\t    opt\
                    \n\t\tcopy odir [CONTCAR, KPOINTS, POTCAR] w. change of INCAR\
                    \n\t    ZPE:\
                    \n\t\tafter opt, calculage zpe\
                    \n\t    band:\
                    \n\t\tafter LCHARG=.TRUE., calculate band structure\
                    \n\tUsage:\
                    \n\t    pypath.sh vas_make_cont.py -d SnO2sc22FH -j band -i i\
                    \n\t    python $sbvas/vas_make_cont.py -d sc34 -nd sc34E5 -j incar -id '{\"ENCUT\": \"500\"}' -io c\
                    "
make.mod_vas        = "modules for make vasp directory\
                    \n\tcalled by 'vas_make_ini.py'\
                    \n\tfunctions:\
                    \n\t    get_vasp_repository\
                    \n\t    make_kpoints\
                    \n\t    get_atoms_4pos(pos='POSCAR')\
                    \n\t\treturns atoms_list, Natoms_list\
                    \n\t    make_mag_4pos(poscar)\
                    \n\t\tcall get_atoms_4pos\
                    \n\t\treturn MAGMOM\
                    \n\t    make_incar\
                    \n\tRun by 'python -m mod_vas -j getmag -s poscar'\
                    \n\t    to get MAGMOM\
                    "
poscar.libposcar    ="module for POSCAR modification\
                    \n\tget_atoms_poscar\
                    \n\tget_poscar(poscar)\
                    \n\t    : read input POSCAR.job and write to wdir as 'POSCAR'\
                    \n\tpos2dirname(poscar)\
                    \n\t    : input POSCAR.job gives dirname of 'job'\
                    \n\tmodify_POSCAR(poscar, atom, atoms=None)\
                    \n\t    : poscar is modified\
                    \n\t    : all the atoms except atom will be fixed for ZPE\
                    "
poscar.pos_lattice  ="[.pl] input POSCAR\
                    \n\t\treturns lattice volume, constants, angles\
                    "
poscar.pos_d2c_ase  =" Read POSCAR/CONTCAR type and write POSCAR in d/c type using ASE"
poscar.pos_sort     ="pos_sort.py POSCAR -al atom_list -z\
                    \n\tsort atoms in POSCAR\
                    \n\treturns POSCARsort\
                    \n\tOptions:\
                    \n\t    -al new atom list for sorted POSCAR\
                    \n\t    -z  True for z-sort in the atom group\
                    \n\tNB: when generate POSCAR via ASE w increasing supercell, atoms in order are replicated\
                    "
poscar.pos_modify   ="change POSCAR\
                    \n\tmain routine in libposcar.modify_POSCAR()\
                    \n\tOptions:\
                    \n\t    -j  zpe, bomb, add, addbomb\
                    \n\t        sel   sel  add    add  : how to add or select atoms\
                    \n\t        coord vel  coord vel   : only coordinate|add velocity\
                    \n\t        addbomb: append atoms after POSCAR\
                    \n\t    -a  add atoms, O12, Fe2, ...\
                    \n\t    -z  location on z-axis: top or float [float]\
                    \n\t\ttop: 4 Ang away from top of system\
                    \n\t\t z1: at z1 position\
                    \n\t\t z2: in the middle of z1 and z2\
                    \n\t    -t input Temperature for velocity distribution in K\
                    \n\t    -o  outfile for POSCAR.outname\
                    \n\t    -suf suffix to be added to POSCAR.name\
                    \n\tUsage:\
                    \n\t    pos_modify.py POSCAR.HfSe2sc34 -j bomb -a O6 -z 1.2 2.3 -t 600 -o POSCAR.HfSe2O12L2 \
                    \n\t    pos_modify.py POSCAR.HfSe2L1O36Hfsub -j sort -s 1-10 -as Se O\
                    "
poscar.libposcar    ="library for pos_modify\
                    \n\tAttributes:\
                    \n\t    def modify_POSCAR(poscar, job, mode_atoms, zpos, temp, htemp, vel_type, outf, r_crit, asort, nlevel\
                    \n\t\tcalled by pos_modify\
                    "
convert.pos2cif      ="vstsscripts/[.pl] convert vasp format(POSCAR, CONTCAR) to cif to be read in MS\
                    \n\tUsage::\
                    \n\t    pos2cif.pl inputfile [outputfile]\
                    \n\t\tinputfile: POSCAR|CONTCAR\
                    \n\t\toutputfile (default): a/CONTCAR -> aCONTCAR.cif\
                    \n\tMS: save to msi to be converted using msi2pos.pl\
                    "
convert.vas2cif     ="[.pl] link to \"pos2cif.pl\""
convert.msi2pos     ="vscripts/[.pl] convert MS(msi) format to POSCAR\
                    \n\tMany diverse forms in /vscripts\
                    \n\tUsage::\
                    \n\t    msi2pos.pl a.msi atom_list\
                    \n\t    msi2pos.pl a.msi Pb Br C N H\
                    \n\t\treturns POSCAR.a\
                    "
run.amp_env_run     ="amp_run.py in (envs) anaconda\
                    \n\t\t   when envs is not (base), detect envs and import proper module\
                    "
run.vas_qsub        = " run vasp in queue\
                    \n\t called by vas_make_cont.py\
                    "
run.vas_env         = " imported from vasp run, make scripts such as\
                    \n\tvas_make_ini.py\
                    \n\tvas_make_cont.py\
                    \n\tvas_make_d2d.py\
                    \n\tvas_make_incar.py\
                    "
clean.clean         =" "

potcar.potcar_gen   =  "potcar_gen.py -pp pseudop\
                    \n\tread POSCAR make POTCAR inside VASP dir\
                    "
potcar.genpotcar    = "linked to potcar_gen.py\
                    "
potcar.potcar_kw    = "extract k-w\
                    \n\tpotcar_kw.py POTCAR[dir] -kw KW\
                    \n\tOptions:\
                    \n\t    POTCAR or if dir, dir/POTCAR\
                    \n\t    -kw:\
                    \n\t\tENMAX -> get 130% of ENMAX\
                    \n\te.g.:\
                    \n\t    potcar_kw.py k774 [-kw ENMAX]\
                    "
incar.libincar     ="module for INCAR modification\
                    \n\tdef modify_incar(incar_in, job, dic=None, opt='ac')\
                    \n\t    incar_in is modified by job_dict\
                    "
incar.incar_change ="change incar using libincar\
                    \n\tOptions:\
                    \n\t    -j: by job has priority\
                    \n\t    -kv: by key-value pair\
                    \n\t    -o: output filename\
                    \n\t\tdefault: makes INCAR_n (new)\
                    \n\t\tif -o INCAR, makes INCAR_o (old)\
                    \n\tUsage:\
                    \n\t    incar_change.py HfSe2mO2sglneba -j cont\
                    \n\t\tcheck outfile of 'incar' by incar_diff.py\
                    "
incar.incar_diff   ="diff_incar.py INCAR1 [INCAR2] -k keys -a -s\
                    \n\tOptions:\
                    \n\t    -j: ['kw','diff'] kw-many files to check kw, diff- compare a few files\
                    \n\t    -f: files/dirs for '-j diff' or [d,f,a] for '-j kw'\
                    \n\t    -k: for keywords for '-j kw'\
                    \n\tUsage:\
                    \n\t    incar_diff.py -j kw -f f -k ivdw\
                    \n\t    incar_diff.py -j kw -f d -k ivdw\
                    \n\t    incar_diff.py FPtb2H2 FPtb2H2hb --- old style\
                    \n\t    incar_diff.py -a -k ENCUT ISTART --- old style\
                    "
dosband.dosall      =   "perl script to decompose lm-decomposded DOSCAR"
dosband.doslm       =   "extact ldos then plot\
                        \n\tOptions:\
                        \n\t    (exclusive)\
                        \n\t\t-z zmin [zmax]  to include atoms inbetween zmin ~ zmax\
                        \n\t\t    if only zmin, zmin-=dz, zmax+=dz\
                        \n\t\t-al atom list index start from 0 from ase gui\
                        \n\t\t    -1 for Tdos, such as -1, 0, 3-9 w. -ash 1 4 4\
                        \n\t\t-ash atom list shape\
                        \n\t\t    len(atomlist) == sum(atomlist_shape)\
                        \n\t    -e [f|float_value]: f for Fermi level, value for VBM shift\
                        \n\t    -p store_true for plot\
                        \n\t\tcalls myplot2D.mplot_nvector\
                        \n\tUsage:\
                        \n\t    doslm.py -z 3.69 -p\
                        \n\t    doslm.py -z 3.69 -p -eV-2.33 for VBM\
                        \n\t    doslm.py -al 8 23 9 10 21 22 -ash 2 4 -p\
                        \n\t    doslm.py -al -1 10-19 -ash 1 10 -eF\
                        \n\t    (Pt-C60-x): doslm.py -al -1 2 54 55 -ash 1 1 1 1\
                        \n\t(3) to plot ldos of slab w.r.t. VBM: obtain VBM in slab (1)\
                        \n\t    doslm.py -al 276-285 286 287 314-317 -ash 10 2 4 -e -1.169\
                        \n\t\tmakes ldos files as much as ash(atom shape)\
                        \n\t\txmgrace ldos_files -p ../dos_sc43.par\
                        "
dosband.pldos        =  "In dos directory\
                        \n\t\tread DOSCAR & plot pldos\
                        \n\t\t(slab 1) To obtain VBM of slab using bulk, run in dos-dir\
                        \n\t\t    $ pldos.py\
                        \n\t\t    To find out VBM\
                        \n\t\t\tfind CBM, which is clean, in bulk area\
                        \n\t\t\tuse bandgap in original bulk calculation\
                        \n\t\t\tCBM - bandgap(bulk) = VBM in slab\
                        \n\t\t(slab 2) To plot band structure of slab w.r.t. VBM\
                        \n\t\t    $cd ~/band-dir: modify FERMI_ENERGY.in to VBM\
                        \n\t\t\trun vaspkit again and use BAND.dat\
                        \n\t\t\t    vaspkit/21/211\
                        \n\t\t\t    xmgrace BAND.dat -p ../band_sc43.par\
                        "
dosband.pldos_nc     =   " given by Lee Y.-G. for siesta"
dosband.doslm_1l     =   "saved doslm.py for one atom list\
                        \n\t\tupgrade to multiple atomlist sets in doslm.py\
                        "
dosband.doscar_split =   "decompose DOSCAR for each atoms written by starmj"
dosband.doscar_modi  =   "remove the first line of each energy loop for removing abnormal value\
                        \n\t\tdefault: read DOSCAR\
                        \n\t\tmove original DOSCAR to DOSCAR_o\
                        \n\t\tmake a new DOSCAR with 1 line less in each block, (natom+1) less line\
                        "
analysis.vas_anal_dos=  "get pdos for every atoms"

procar.procar_kb    =   "[procar_kb.pl] renamed from 'get_Procar.pl'\
                        \n\t\tTo extract k,b sets for the related peak in dos plot\
                        \n\t\tcontrol inside script:\
                        \n\t\t    ldos or pdos\
                        \n\t\t    dos amplitude criteria for each dos and pldos\
                        \n\t\tinput: energy range <- find energy in LDOS of the atom\
                        \n\t\tUsage:\
                        \n\t\t    procar_kb.pl [d=dirname] [e=]E1[:Emax] \"atom indices\"\
                        \n\t\t\td=dirname: optional\
                        \n\t\t\tE1: center of DOS\
                        \n\t\t\tEmin:Emax for energy range\
                        \n\t\t\t\'indices': atom index which starts from 1\
                        \n\t\tE.g.:\
                        \n\t\t    procar_kb.pl e=-6:-5 '3 58'\
                        "
procar.get_procar       =   "[.pl] similar to procar_kb.pl"
procar.get_procar_kband = "[.pl] similar to get_procar.pl"
procar.procar           =   "To extract and draw band\
                        \n\t\t'to load PROCAR' can make a memory problem\
                        "
procar.procar_byline    =   "the same as procar.py\
                        \n\t\tread by line and calculate to save memory for what?\
                        "
outcar.outcar_zpe_ts    = "Read OUTCAR: calculate T*S energy\
                        \n\t\tOption:\
                        \n\t\t    d [dir with OUTCAR|fname]\
                        \n\t\t    -na number of atoms to be calculated for frequency\
                        \n\t\tUsage:\
                        \n\t\t    outcar_zpe_ts.py OUTCAR_test_1_catO2_vib -na 2\
                        "
outcar.outcar_zpe_ts_mj = "original version of zpe from mjstar"
md.oszicar              = "MD is stored in OSZICAR and analysis\
                        \n\t\tOptions: \
                        \n\t\t    -y    keys lined by '+' for Etot Efree E0pot Ekin Skin Spot\
                        \n\t\t          Etot = Efree(E0pot)+Ekin+Skin+Spot\
                        \n\t\tUsage: \
                        \n\t\t    oszicar.py mdPtbT0300t115s1-last -y Etot Skin Spot+Skin \
                        \n\t\t    oszicar.py mdPtbT0300t115s1-last -y Etot Efree+Ekin Spot+Skin -iy 2 \
                        "
analysis.vas_anal    =   "(.sh) Charge analysis of Bader\
                        \n\tjobs: bader bader2(spin) convasp dos bchg end\
                        \n\tUsage: Run just above vasp directory\
                        \n\t    vas_anal.sh bader dirname\
                        "
charge.charge_bader =   "included in vas_anal.sh\
                        \n\tRead POTCAR for ZVAL\
                        \n\tRead POSCAR for atom list\
                        \n\toutput: Bader charge in 'bader_pcharge.dat'\
                        \n\tto compare two configurations w. the same atom index\
                        \n\t    paste adir/bader_pcharge.dat bdir/bader_pcharge.dat | awk '{print $1, $6-$3}'\
                        "
charge.vas_anal     =   analysis.vas_anal

ase.ase_fconvert    =""
ase.ase_vasp        =""
ase.ase_zpe         =""
ase.ase_pos_d2c     = poscar.pos_d2c_ase
classobj_dict={'MAKE': make, 'RUN': run, 'CLEAN': clean} 

#def classify(Lclassify, work, class_name, job, fname,HL, elimit, nc, Lgraph):
def classify(Lclassify, work, class_name, job):
    
    mdir = os.path.dirname(__file__)
    print(f"List directory of {mdir} ")
    #exe, mod = dir_files(mdir)
    exe, mod, dirs, d_link = dir_all(mdir)
    sort_exe = sorted(exe)
    sort_mod = sorted(mod)
    sort_dir = sorted(dirs)

    if sort_dir:
        print("Directories:: ")
        if not Lclassify:
            for f in sort_dir:
                print(f"    {f}")
        else:
            for instance in MyClass.instances:
                for gkey in globals().keys():
                    if gkey == instance.name:
                        break
                if work != instance.name:
                    ckeys = dir_classify_n(sort_dir, instance.name, globals()[gkey], Lwrite=1) # globals()[instance.name] is not working
                else:
                    ckeys = dir_classify_n(sort_dir, instance.name, globals()[gkey], Lwrite=0)
                    for ckey in ckeys:
                        print(f"    {ckey}.py\t:: {globals()[gkey].__dict__[ckey]}")
            print("  == not classified")
            for f in sort_dir:
                print(f"    {f}")
    
    if sort_exe: 
        print("Executable:: ")
        if not Lclassify:
            for f in sort_exe:
                print(f"    {f}")
        else:
            ### confer "ini_pypbs.py", MyClass.instances is a list of string as class variables
            for instance in MyClass.instances:
                ### globals() includes MyClass() instances as keys
                for gkey in globals().keys():
                    if gkey == instance.name:
                        break
                ### work == None without -w
                if work != instance.name:
                    ckeys = dir_classify_n(sort_exe, instance.name, globals()[gkey], Lwrite=1) # globals()[instance.name] is not working
                else:
                    ckeys = dir_classify_n(sort_exe, instance.name, globals()[gkey], Lwrite=0)
                    for ckey in ckeys:
                        print(f"    {ckey}.py\t:: {globals()[gkey].__dict__[ckey]}")
            print("  == not classified")
            for f in sort_exe:
                print(f"    {f}")
    if sort_mod:
        print("Module:: ")
        if not Lclassify:
            for f in sort_mod:
                print("    {}".format(f))
        else:
            for instance in MyClass.instances:
                for gkey in globals().keys():
                    if gkey == instance.name:
                        break
                if not work or work != instance.name:
                    ckeys = dir_classify_n(sort_mod, instance.name, globals()[gkey], Lwrite=1)
                else:
                    ckeys = dir_classify_n(sort_mod, instance.name, globals()[gkey], Lwrite=0)
                    for ckey in ckeys:
                        print(f"    {ckey}.py\t:: {globals()[gkey].__dict__[ckey]}")
            print("  == not classified ")
            for f in sort_mod:
                print(f"    {f}")

    if job == 'amp':
        print("For AMP::")
        file_conversion()
        #run_amp(fname,HL, elimit, nc, Lgraph)
    ### print dictionary here
    if job in classobj_dict.keys():
        name_class = classobj_dict[job]
        for key in name_class.__dict__.keys():
            print(f" {job} :: {name_class.__dict__[key]}")

    print("\nClass Instances:: ", end='')
    for instance in MyClass.instances:
        print(f"{instance.name}", end=' ')
    print("\n\t    -w for detail")
    #print(f"#Comment: -c    for not classification")

    return 0        

def main():

    parser = argparse.ArgumentParser(description="display Usage for ~/py_ai")
    parser.add_argument('-c', '--classify', action="store_false", help="classify files ")
    parser.add_argument('-w','--work',  help="several explanation option ")
    parser.add_argument('-cn', '--cname', help="detail for each class ")
    parser.add_argument('-j','--job',  help="[val,train,test] ")
    #parser.add_argument('-f','--file',  help="input energy data file ")
    #parser.add_argument('-hl','--hidden_layer',nargs='*', default=['4','4','4'], help="list of number of Hidden Layer")
    #parser.add_argument('-el','--energy_limit',default=0.001, type=float,  help="energy_limit for training")
    #parser.add_argument('-nc','--ncore',  help="number of parallel process")
    #parser.add_argument('-g','--graph', action='store_true',  help="draw graph or not")
    args = parser.parse_args()

    #classify(args.classify, args.work, args.cname, args.job, args.file, args.hidden_layer, args.energy_limit,args.ncore, args.graph )
    classify(args.classify, args.work, args.cname, args.job)

if __name__ == "__main__":
    main()
