#!/bin/bash

awk -v file=$1 'BEGIN { i=0; sum=0.0}
{ if ($9)
    { print $9; sum+=$9;}}
END{printf "%s of %s: %10.6f\n", "total charge", file, sum} ' $1


