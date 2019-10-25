#!/bin/sh
# /usr/local/bin/ps2png
for i in "$@"
do
        i_new=`dirname $i`/`basename $i .ps`.png
        echo converting $i to $i_new
        #gs -dNOPAUSE -sDEVICE=png256 -sOutputFile=$i_new -q -dBATCH $i -r1600
	convert -rotate 90 $i $i_new
#        convert -resize 50% $i_new $i_new
done
#convert -transparent white test.png test.png
