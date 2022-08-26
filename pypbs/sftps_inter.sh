#!/bin/bash

stem="ampwater.tar.parta"


for x in {g..g}; do
	fname=${stem}${x}
	tmpsize=0
	FILESIZE=$(stat -c%s $fname)
	echo $tmpsize  $FILESIZE
	while [[ $FILESIZE -gt $tmpsize  ]]; do
		sleep 10
		tmpsize=$FILESIZE
		FILESIZE=$(stat -c%s $fname)
		echo sleep 10s: $FILESIZE
	done
	echo filesize doesnot change: $tmpsize $FILESIZE
	echo file download done
	sftp joonho@platinum.kaist.ac.kr:/home/joonho <<EOF
	put $fname
	exit
EOF
	echo uploaded $fname
	echo remove $fname
	#rm $fname
done

	
