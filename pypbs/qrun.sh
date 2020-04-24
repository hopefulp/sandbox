#!/bin/bash

### Modify here
lsoftware=( qchem amp )
software=${lsoftware[1]}
### NData_total can be region, NData_train can be list of [ntrain,ntest] 
if [ $# -lt 1 ]; then
    exe=`basename $0`
    echo "Usage:: Do Not Skip Any INPUT !!!"
    echo "Usage:: List:                                                         Hidden_Layer    Force_Lim NData_tot(region)"
    echo "Usage:: $exe sub_Node \$ampjob \$sub_dir \$qjob     \$fin    \$np   \$mem   '\$hl'   \$el     \$fl      ntotal ntrain int '\$data'"
    echo "Usage:: $exe [qsub|di] tr|te  sub_dir  qname      f_input nproc mem[nG] HL     E_limit F_limit ND_tot ND_tr dtype dlist"
    echo "try  :: $exe [qsub|di] tr     direct   N1500      OUTCAR  16    4      '10 10' 0.001   0.01    5000   1500   int  '0 1500'"
    echo "try  :: $exe [qsub|di] tr     direct   N1500div32 OUTCAR  16    4      '10 10' 0.001   0.00    5000   1500   div  '3 2'"
    echo "        $exe     di    te      test    N1500      OUTCAR   4    4      '10 10' 0.001   0.00    5000   1500   int  '0 500'"
    echo "     :: $exe     di    te      test1   N1500      OUTCAR   4    4      '10 10' 0.0005  0.00  '4500 5000' 1500 int '0 500'"
    echo " ND_tot controls ND_train for tr ND_test for te"
    exit 1
fi    
###### Default values
### common
###
de_qjob=qname
de_dir="test"
de_inf=OUTCAR
de_np=16
### qchem
de_mem="2"     # for default
### AMP
d_hl="10 10"
d_el=0.001
d_fl=0.01
d_dtype=int
d_data="0 1500"
d_ampjob=tr
d_ndata=5000
d_ndata_tre=1500
hname=`hostname`
qsub='qsub'
if [ $hname == "chi" ]; then
    fin=${1:-$de_inf}
    np=${2:-$de_np}
#elif [ $hname == "login" ]; then
else
    if [ $1 == 'qsub' -o $1 = 'di' ]; then
        qsub=${1:-$qsub}
        shift
    fi
    ampjob=${1:-$d_ampjob}
    shift
    sub_dir=${1:-$de_dir}
    qjobname=${2:-$de_qjob}
    fin=${3:-$de_inf}
    np=${4:-$de_np}
    mem=${5:-$de_mem}
    hl=${6:-$d_hl}
    el=${7:-$d_el}
    fl=${8:-$d_fl}
    ndata=${9:-$d_ndata}
    ndatatre=${10:-$d_ndata_tre}
    dtype=${11:-$d_dtype}
    dlist=${12:-$d_data}
fi

echo $@

echo memory $mem G
echo HL $hl

if [ $ampjob == "tr" ]; then
    ampfiles=( 'amp-fingerprint-primes.ampdb' 'amp-fingerprints.ampdb' 'amp-neighborlists.ampdb' 'OUTCAR' )
elif [ $ampjob == "te" ]; then
    ampfiles=( 'amp-fingerprint-primes.ampdb' 'amp-fingerprints.ampdb' 'amp-neighborlists.ampdb' 'OUTCAR' 'amp.amp' )
fi

#if [ $hname == 'login' -a $software == 'amp' ]; then
if [[ $hname != 'chi' &&  $software == 'amp' ]]; then
    echo "mkdir $sub_dir and cd $sub_dir"
    if [ ! -e $sub_dir ]; then
        mkdir $sub_dir
    fi

    cd $sub_dir
    for ampf in ${ampfiles[@]}
        do
            if [ -e $ampf ]; then
                echo "checked $ampf"
            elif [ -e ../$ampf ]; then
                ln -s ../$ampf $ampf
                echo "made $ampf"
            else
                echo "there is not amp data here and ../"
                exit 2
            fi  
        done            
    echo "Usage:: $0 $ampjob $sub_dir $qjobname $fin   $np    $mem         '$hl'        $el         $fl     '$data'"
    echo "Usage:: $0 tr|te [pa] sub_dir qname  f_input nproc memory[nG] hidden_layer Energy_limit Force_limit Data_region"
    echo "HL = $hl"
    if [[ $qsub =~ 'qsub' ]]; then
        st="qsub_server.py amp -qj $qjobname -i $fin -nt $ndata -n $np -hl '$hl' -el $el -fl $fl -dt $dtype -dl '$dlist' -m $mem -j $ampjob"
        read -p "$st | Will you run? [enter/no]" var
        if [ -z $var ]; then
            $st
        fi
    elif [[ $qsub =~ 'di' ]]; then
        if [ $ampjob == 'tr' ]; then
            ### -nt total for training, -ntr ntrain for test
            st="amp_run.py -f $fin -j tr -nt $ndata -dt $dtype -dl $dlist -nc $np -hl $hl -el $el -fl $fl "
        else
            st="amp_run.py -f $fin -j te -nt $ndata -ntr $ndatatre -dt $dtype -dl $dlist -nc $np -hl $hl -el $el -fl $fl "
        fi
        read -p "$st | will you run? [enter/no]" var
        if [ -z $var ]; then
            $st
        fi
    fi
fi

