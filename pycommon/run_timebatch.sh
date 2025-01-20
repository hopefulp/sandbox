#!/usr/bin/bash

# run at crontab

dd=$(date | awk '{print $3$4}' | cut -d: -f 1)

cd /scratch/x3075a02/HfSe2
kpy vas_make_ini.py -j fake -s POSCAR.test -sj sp -al -ra -d d${dd} -n 6
