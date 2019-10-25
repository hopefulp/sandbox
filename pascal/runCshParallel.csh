#!/bin/tcsh
#!/bin/csh

if ($#argv < 2) then
    echo "usage: $0 list script [numprocs=8]"
    exit(1);
endif

set runScript = $2
set list = ($1)
set nprocs = 8

if ($#list == 0) then
    echo "ERROR: No valid files found while searching $list"
    exit(1)
endif

if !((-e $runScript) && (-r $runScript)) then
    echo "ERROR: Cannot find script $runScript"
    exit(1)
endif

set nprocs = 8
if ($#argv > 2) then
    set nprocs = $3
endif

set done = 0
set counter = 0
set offset = 0;
set iter = 1;

echo Iteration $iter
while (! $done)
    if ($counter < $nprocs) then
        @ counter = $counter + 1;
        @ index = $counter + $offset
        if ($index <= $#list) then
	    echo $runScript $list[$index]
	    $runScript $list[$index] >/dev/null &
	else
	    wait;
	    set done = 1           
        endif
    else
        @ offset = $index
        set counter = 0
        if ($offset >= $#list) then
	    wait;
            set done = 1
        else
            @ iter = $iter + 1
	    wait;
            echo Iteration $iter
        endif
    endif
end

