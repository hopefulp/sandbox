#!/bin/bash

#awk '{for(i=1; i<=NF; i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++) {printf "%14.6f", sum[i]/NR}; print ""; for (i=1;i<=NF;i++) {printf "%14.6f", sqrt((sumsq[i]-sum[i]^2/NR)/NR)}; print "\n"}' $1
awk '
function abs(x)
{
    return ((x < 0.0) ? -x : x)
}

{
    for (i = 1; i <= NF; i++) 
    {
        sum[i] += $i; 
        sumsq[i] += ($i)^2
    }
} 

END {
    printf "Average\t"
    for (i = 1; i <= NF; i++) 
    {
        printf "%14.6f", (sum[i] / NR)
    }; 

    printf "\nStdev\t" 
    for (i = 1; i <= NF; i++) 
    {
        value = (sumsq[i] - sum[i]^2 / NR) / NR
        if ( value >= 0.000000001 )
            printf "%14.6f", sqrt(value)
        else
            printf "%14.6f", 0.0
    }; 

    print "\n"
}' $1
