#!/bin/bash

if [ $# -eq 1 ]; then
    ext=`cut -d. -f2 <<< $1`
    echo $ext
else
    echo "calculate total charge of .top or .mol2 file"
    echo "it requires one argument of [.top|.mol2]"
    exit
fi

case $ext in
    "top")
        awk '/atoms/{L=1} L{sum+=$7} /bonds/{L=0} END {printf "total charge : %10.6f\n", sum} ' $1
    ;;
    "mol2")
        awk '/ATOM/{L=1} L{sum+=$9} /BOND/{L=0} END {printf "total charge : %10.6f\n", sum} ' $1
    ;;
    *)
        echo "Error:: The extension should be one of [.top|.mol2]"
    ;;
esac

