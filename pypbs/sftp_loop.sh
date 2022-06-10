#!/bin/bash

stem="ampwater.tar.parta"

tmp=${stem}d
for x in {e..k}; do
        fname=${stem}${x}
        while [ -f "$tmp" ]; do
                echo waiting for $fname to be moved
                sleep 10
        done    
        echo start to download $fname
        sftp joonho@mlet2.kaist.ac.kr:/archive/joonho <<EOF
        get $fname
        exit
EOF
        echo downloaded $fname
        tmp=${fname}
done

