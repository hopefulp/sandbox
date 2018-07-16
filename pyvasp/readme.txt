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


