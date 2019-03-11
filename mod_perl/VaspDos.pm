# by J. Park
# usage &VaspDos::lm_quantum($l,$m)
# usage &VaspDos::lm_quantum($l,$m,"spin")
# without 3rd argument, to get " s     py     pz     px    dxy    dyz    dz2    dxz    dx2 "
# with    3rd argument, to get " s-u,d py-u,d ...    "

package VaspDos;

use Exporter;

@ISA=qw(Exporter);
@EXPORT=qw(lm_quantum);

### from the quantum number get column number
### input l, m, and nothing for PROCAR and something for DOSCAR
sub lm_quantum{
    my($l,$m,$nspin,$spin)=@_;
    my($ic,$fc,$im,$fm,$spin_my,@n_mag_orbitals,@n_acc_mag_orbitals,@m_order,@nl_col);

    @n_mag_orbitals=qw(1 3 5);
    @n_acc_mag_orbitals=qw(1 4 9);
    #if(defined($lspin)){	$spin=2; #print "This is applied to DOSCAR Format with spin\n";
    #}else{			$spin=1; #print "This is applied to PROCAR Format regardless of spin\n";
    #}
    ###  human ordering  (['s'],['x','y','z'],['xy','yz','xz','x2-y2','z2'])
    ###  vasp  ordering  (['s'],['y','z','x'],['xy','yz','z2','xz','x2-y2']) 
    ###  human to vasp   (  1     3   1   2     1    2    4    5    3 )
    ###  human to vasp   (  0     2   0   1     0    1    3    4    2
    @m_order=([0],[2,0,1],[0,1,3,4,2]);
    
    ### for l-angular momentum
    use Switch;
    if(defined($l)){
    	switch($l){
	    case 0	{	$ic=0;	$fc=$ic+$n_mag_orbitals[$l]*$nspin-1;
				@nl_col=&list($ic,$fc,$spin);
		}
	    else	{
		if(!defined($m)){
#			$ic=2;	$fc=$ic+$n_mag_orbitals[$l]-1;
			$ic=$n_acc_mag_orbitals[$l-1]*$nspin;	$fc=$ic+$n_mag_orbitals[$l]*$nspin-1;
			@nl_col=&list($ic,$fc,$spin);
		}else{	### $ic starter is OK but magnetic quntum number 
			$im=$n_acc_mag_orbitals[$l-1]*$nspin+$m_order[$l][$m]*$nspin;
			$fm=$im+$nspin-1;					
			@nl_col=&list($im,$fm,$spin);
   		}

    		}
    	}
    }else{
	$ic=0;
	$fc=$n_acc_mag_orbitals[2]*$nspin-1;
	@nl_col=&list($ic,$fc,$spin);
    }
#    if   (defined($lspin) and $lspin =~ /[uU]/){ $fc=$ic;}
#    elsif(defined($lspin) and $lspin =~ /[dD]/){ $ic=$fc;}
    return (@nl_col);
}

sub list{
    my($ini,$fin,$spin)=@_;
    my($i,@list);

    if(!defined($spin)){
        for($i=$ini;$i<=$fin;$i++){
            push(@list,$i);
        }
    }elsif($spin==1){
	for($i=$ini;$i<=$fin;$i+=2){
            push(@list,$i);
        }
    }elsif($spin==-1){
	for($i=$ini+1;$i<=$fin;$i+=2){
            push(@list,$i);
        }
    }
    return @list;
}

jp;
