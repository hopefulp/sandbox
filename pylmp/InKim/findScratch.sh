#!/bin/bash

display_usage() {
	printf "\n"
	printf "findScratch.sh:\tFind a specified directory in /scratch/$USER in nodes.\n"
	printf "\t\tNode lists will be taken from /home/noische/hosts\n"
	printf "\n"
	printf "usage: findScratch.sh directory\n"
	printf "\n"
}

# if only one arguments not supplied, display usage 
if [  $# -ne 1 ] 
then 
	display_usage
	exit 1
fi 
 
# check whether user had supplied -h or --help . If yes display usage 
if [[ ( $# == "--help") ||  $# == "-h" ]] 
then 
	display_usage
	exit 0
fi

scratch=/scratch/$USER

printf "$0: find $1 from /scratch/$USER in $HOST nodes\n\n"
for node in `cat /home/noische/hosts`
do
	if `echo 'test -d '"$scratch/$1"' && exit 0 || exit 1' | ssh -q "${node}" sh`; then
		printf "$node $scratch/$1 exists\n"
	else
		printf "$node not found\n"
	fi
done
