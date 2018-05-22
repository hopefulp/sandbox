#!/bin/tcsh
echo "+--------------------------------------------------------------------------+"
foreach nn(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20)
echo "[$nn]";qstat -rn | grep node$nn/
end

#qstat -n | grep -P "node1/|node2/|node3/|node4/|node5/|node6/|node7/|node8/|node9/|node10/"
