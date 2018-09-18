#!/bin/bash

old_dir=/qcfs/biduri/mof-74/hydrocarbon/from-exp
new_dir=`pwd`

echo "new_dir=$new_dir"

cd $old_dir

# the file name holds the present directory "./"
find . -name "CONTCAR" > $new_dir/test.txt
	
cd $new_dir

while IFS= read -r NAME 
    do
	NAME=${NAME/.\//}
  	cp -v "$old_dir/$NAME" "$new_dir/${NAME//\//_}"
    done < "test.txt"



