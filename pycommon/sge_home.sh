#!/bin/bash

echo "sed '1s:home:gpfs/home:' $1 > sge_$1"
echo "chmod 755 sge_$1"
