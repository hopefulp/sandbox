
function viamp {
    dir=${1:-'.'}
    vi $dir/amp-log.txt
    }

function taila {
    dir=${1:-'.'}
    tail $dir/amp-log.txt
    }
function tailf {
    dir=${1:-'.'}
    tail -f $dir/amp-log.txt
    }
function wl {
    if [ $# -eq 4 ]; then
        dir1=$1
        shift
    else
        dir1='.'
    fi
    for dir in $@; do
        echo $dir1/$dir
        ls -al $dir1/$dir/loose | wc -l
        done
    }
#if [[ $1 =~ ^[0-9] ]]; then
#    grep "optimization un" */amp-log.txt -B 1 | awk '{if(NF==11) {print $8} else if(NF==12) {print $9} }'
#if [[ $1 =~ ^[a-zA-Z] ]]; then
        #else
        #    grep "optimization un" $1*/amp-log.txt -B 1 | awk '{if(NF==11) {print $8} else if(NF==12) {print $9} }'
function dir_trf {
    if [[ $# -eq 0 ]]; then
        grep "optimization un" */amp-log.txt -B 1 | awk '{if(NF==11) {print substr($1,1,index($1,"/")-1), $8} else if(NF==12) {print substr($1,1,index($1,"/")-1),$9} }'
    elif [[ $# -eq 1 ]]; then
        grep "optimization un" $1*/amp-log.txt -B 1 | awk '{if(NF==11) {print substr($1,1,index($1,"/")-1), $8} else if(NF==12) {print substr($1,1,index($1,"/")-1),$9} }'
    elif [ $# -eq 2 ]; then
            grep "optimization un" $1*/$2*/amp-log.txt -B 1 | awk '{if(NF==11) {print substr($1,1,index($1,"/")-1), $8} else if(NF==12) {print substr($1,1,index($1,"/")-1),$9} }'
    fi
    }
#        if [[ $1 =~ ^[0-9] ]]; then
#            grep RMSE */test_fstat_acc.txt | awk '{print $5}'
#        elif [[ $1 =~ ^[a-zA-Z] ]]; then
#        grep RMSE $1*/test_fstat_acc.txt | awk '{print $5}'
function dir_tef {
    if [[ $# -eq 0 ]]; then
        grep RMSE */test_fstat_acc.txt | awk '{print substr($1,1,index($1,"/")-1), $5}'
    elif [[ $# -eq 1 ]]; then
        grep RMSE $1*/test_fstat_acc.txt | awk '{print substr($1,1,index($1,"/")-1), $5}'
    elif [[ $# -eq 2 ]]; then
        grep RMSE $1*/$2*/test_fstat_acc.txt | awk '{print substr($1,1,index($1,"/")-1), $5}'
    fi
    }


