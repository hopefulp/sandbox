#!/usr/bin/perl
# written by Joonho Park
# read a.inp and write a.out

$fname=$ARGV[0];
@fname=split(/\./,$fname);
$fname_pre=$fname[0];
$fout="$fname_pre.mol";
$fxyz="$fname_pre.xyz";
print "Overwrite mol file: $fout and $fxyz\n";
open(IN,"<$fname");
open(OUTmol,">$fout");
open(OUTxyz,">$fxyz");

if($ARGV[1] eq ""){
    $njob=1;
}else{
    $njob=$ARGV[1];
}

$job_tag="OFF";
$keyword_job="Job";

$keyword_spin="molecule";
$endword_spin="end";
$flag_spin="OFF";
$tag_spin="ON";

$keyword_opt="Optimization Cycle";
$endword_opt="OPTIMIZATION CONVERGED";
#$flag_opt="OFF";
$tag_opt="OFF";

#$tag_geometry="GEOMETRIES";	# for MOLDEN output
$keyword_ext="Step";
$flag_coord="OFF";

@atom_coord=();
$i=0;
$imol=0;
$i_atom=0;
### read three atoms's coordinates
while($line=<IN>){
    $i++;
    ### extract charge and spin
    ### if for block execution
    if($tag_spin eq "ON"){
        # for block reading
        if($flag_spin eq "OFF" ){
            if($line =~ /\$$keyword_spin/ ){ 
                print "found \$$keyword_spin\n";
                $flag_spin = "ON";
            }
            next;
        }else{
            @field=split(/\s+/,$line);
            if($field[0] eq "") {shift(@field);}
            if($imol==0){
                $charge=$field[0]; $multiplicity=$field[1];
                $imol++;
            #### just count mol
            }elsif($line !~ /$endword_spin/){
                #### for blank line in molecule section
                if($field[0] =~ /\w/) {  $imol++; }
            }else{
                $flag_spin = "OFF";
                print $charge, "***",$multiplicity,"\n";
                $tag_spin = "Done";   # got it, skip
                $natom=$imol-1;
                print "num of atom = ",$natom,"\n"; 
            }
            next;
        }
    }

    ### find n-th job here
    if($job_tag eq "OFF"){
        if($njob == 1){
            $job_tag="ON";
            next;
        }else{
            if($line =~ /$keyword_job/){
                @field=split(/\s+/,$line);
                if($field[0] eq "") {shift(@field);}
                if($field[1] == $njob){
                    print "$njob-th Job\n";
                    $job_tag="ON";
                }
            }
            next;
        }
    }
    # if job_tag eq "ON", reaches here

    ### count num of optimization step
    if($tag_opt eq "OFF"){
        if($line =~ /$keyword_opt/) {
            if($i_opt == 0) {print "Optimization job\n";}
            @field=split(/\s+/,$line);
            if($field[0] eq "") {shift(@field);}
            $i_opt=$field[2];
        }elsif($line =~ /$endword_opt/){
            $tag_opt = "Done";
            print "Optimization step: $i_opt\n";
            $n_opt=$i_opt;
        }
        next;
    }
    ### extract coordinates
    if($flag_coord eq "OFF"){
        if($line =~ /$keyword_ext/){
            @field=split(/\s+/,$line);
            if($field[0] eq "") {shift(@field);}
            if($field[1] == $n_opt){
                $flag_coord = "ON"; 
                next;
            }
        }
    }else{
        if($i_atom < $natom){
            push(@atom_coord,$line);
            $i_atom++;
        }else{
            last;
        }
    }
}
close(IN);

print @atom_coord;
for($i=0;$i<=$#atom_coord;$i++){
    if($i==0) {print OUTxyz "$natom\n"; print OUTxyz "$keyword_ext $n_opt\n";}
    print OUTxyz $atom_coord[$i];
}

for($i=0;$i<=$#atom_coord;$i++){
    if($i==0) {print OUTmol "\$molecule\n"; print OUTmol "  $charge   $multiplicity \n";}
    print OUTmol $atom_coord[$i];
    if($i==$#atom_coord) {print OUTmol "\$end\n";}
}
close(OUTmol);
close(OUTxyz);

