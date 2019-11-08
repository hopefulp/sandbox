from common import MyClass

water = MyClass()

water.order = "calcube makecube pdb2bgf makelmp_in"
water.calcube = "chem_math.py -m H2O -d 1.0 -n 64"
water.makecube = "packmol < water_n64.inp\
                    \n\t\t\tmakes a.pdb\
                    \n\t\tgopack to see input file"
water.pdb2bgf = "babel -ipdb water_n64.pdb -obgf water_n64.bgf"
water.modify_bgf = "include pbc(CRYSTX); change FF; change charge\
                    \n\t\t\tNB: pbc FF coord's bond(CONNECT) are important\
                    \n\t\t\tcheck: vmd
                    "
water.makelmp_in = "LammpsInput.pl -b water_n64.bgf -f $FF/spcew.ff -s water_n64 -t full\
                    \n\t\t\t:makes in.water_n64 data.water_n64
                    \n\t\t\tcp in.water_n64 data.water_n64 water_n64.bgf to lammps_work_dir"
water.run_lmp = "mpirun -n 4 ~/.local/bin/lmp -in in.asps -log water.log
                    "
water.vmd2poscar =         "vmd load a.bgf\
                    \n\t\t\tvmd>pbc set { }\
                    \n\t\t\tvmd>pbc box\
                    \n\t\t\tvmd load b.traj on a.bgf (stride for skip)\
                    \n\t\t\tselect one snapshot\
                    \n\t\t\tvmd> pbc wrap; move all atoms into the box\
                    \n\t\t\tvmd save to poscar"
water.vmdpos2pos = "vpos_rearrange.py n64.vmdpos -af water_n64.bgf"
water.vasp="conver vasp"
vasp.make_incar =   "make_incar -t re0 -d d3\
                    \n\t\t\makes incar.key\
                    make_incar.py -t re0 -d d3\
                    \n\t\t\hybrid runs with WAVECAR \
                    make_incar.py -t re0 -d d3 -md nvt\
                    \n\t\t\tto run MD
                    "
vasp.make_ini = "vmake_ini.py -a O H -d dirname\
                    \n\t\t\tKPOINTS=gamma, POTCAR from VaspINI by default and use 'incar.key' for INCAR"
vasp.run = "mpirun -n 4 ~/sciwares/VASP/vasp.5.4.4/bin/vasp"

sge.vasp = "qsub -N pe500 -v np=12 -v dir=pe500 $SB/pypbs/sge_vasp.csh"
pbs.vasp = "qsub mpi_vsp.sh in Nurion@KISTI\
            \n\t\t\tqsub -N dirname $SB/pypbs/pbs_vasp\
            "


