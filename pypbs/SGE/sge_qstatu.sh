#!/bin/bash

if [ $1 ]; then
    user=$1
else    
    user="\*"
fi    
echo $user
qstat -u \* | \
sed '/---/d' | \
awk ' BEGIN {i=0;}
{
    i++; 
    if (i==1) 
        {printf "%8s %7s %10s %10s %5s   %15s      %2s %5s %8s %5s\n",$1, $2, $3, $4, $5, $6, $7, $8, $9, $10 } 
    else 
        {printf "%8d %7.5f %10s %10s %5s %10s %8s", $1, $2, $3, $4, $5, $6, $7; 
        if ($5 == 'qw') 
            {printf  "%2ds\n", $8} 
        else
            { printf "%15s %2d\n", $8, $9}
        }
}'


