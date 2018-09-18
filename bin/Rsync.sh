#!/usr/bin/bash

home_dir=( Backup Downloads Exe Research exe-py modules scripts softwares vmware)

for dir in `echo ${home_dir[*]}` 
    do
        echo "rsync -avz $local_dir$dir /run/media/joonho/Seagate\ Backup\ Plus\ Drive/$target"
    done

rsync -avz /sharewin/* /run/media/joonho/Seagate\ Backup\ Plus\ Drive/Chi_sharewin

