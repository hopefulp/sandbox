#!/bin/bash
### in mlet
SB=/gpfs/home/joonho/sandboxg

message="more in server"

arg1=${1:-$message}

echo "git add . -A"
git add . -A
echo "git commit -m \"$arg1\""
git commit -m "$arg1"
echo "git push"
if [ $# -eq 2 ]; then
    git push hub HEAD:master
else
    git push
fi
