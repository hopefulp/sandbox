#!/bin/bash

list=( 'echo' 'rm' )

if [ $# == 0 ]; then
    job='echo'
else
    job=$1
fi

str="0* ch* p* G*"
#echo rm -rf $str

case $job in
    'echo')
        echo "echo" $str
        echo input one of [ ${list[@]} ]
        ;;
    'rm')
        rm -rf $str
        ;;
    *)
        echo "input $list"
        ;;
esac

