#!/usr/bin/perl
# written by Joonho Park
# read a.inp and write a.out

$f_full=$ARGV[0];
if($f_full =~ /\//){
    @fname=split(/\//,$f_full);
    $fname=$fname[$#fname];
}else{ $fname=$f_full; }
print "$fname\n";

@fname=split(/\./,$fname);
$fname_pre=$fname[0];
$fout="$fname_pre.mol";
$fxyz="$fname_pre.xyz";
print "Overwrite mol file: $fout and $fxyz\n";
open(IN,"<$f_full");
open(OUTmol,">$fout");
open(OUTxyz,">$fxyz");

if($ARGV[1] eq ""){
    $njob=1;
}else{
    $njob=$ARGV[1];
}

$job_tag="OFF";
$KW_job="Job";

$KW_spin="molecule";
$endword_spin="end";
$flag_spin="OFF";
$tag_spin="ON";

$KW_opt="Optimization Cycle";
$KW_atom="ATOM";
$endword_opt="OPTIMIZATION CONVERGED";
#$flag_opt="OFF";
$tag_opt="OFF";

#$tag_geometry="GEOMETRIES";	# for MOLDEN output
$KW_ext="Step";
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
            if($line =~ /\$$KW_spin/ ){ 
                print "found \$$KW_spin\n";
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
            if($line =~ /$KW_job/){
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

    ### Read between opt block
    if($tag_opt eq "OFF"){
        if($line =~ /$KW_opt/) {
            @atom_coord=();
            if($i_opt == 0) {print "Optimization job\n";}
            @field=split(/\s+/,$line);
            if($field[0] eq "") {shift(@field);}
            $i_opt=$field[2];           # extract opt step in "Optimization Cycle: num"
        ### get atom coordinates in Optimization Cycle step
        }elsif($line =~ /$KW_atom/){
            $tag_ATOM="ON";
        }elsif($tag_ATOM eq "ON"){
            @field=split(/\s+/,$line);
            if($field[0] eq "") {shift(@field);}
            $iatom=int($field[0]);
            $str = "  ".join("   ", $field[1], $field[2], $field[3], $field[4], "\n");
            push(@atom_coord, $str);
            if($iatom == $natom){ $tag_ATOM="OFF";}
        }elsif($line =~ /$endword_opt/){
            $tag_opt = "Done";
            print "Optimization step: $i_opt\n";
            $n_opt=$i_opt;
        }
        next;
    }
    ### extract coordinates: at MOLDEN FORMAT [GEOMETRIES] "Step"
    if($flag_coord eq "OFF"){
        if($line =~ /$KW_ext/){
            @field=split(/\s+/,$line);
            if($field[0] eq "") {shift(@field);}
            if($field[1] == $n_opt){
                $flag_coord = "ON"; 
                @atom_coord=();      # re-initialize if there is MOLDEN
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
    if($i==0) {print OUTxyz "$natom\n"; print OUTxyz "$KW_ext $n_opt\n";}
    print OUTxyz $atom_coord[$i];
}

for($i=0;$i<=$#atom_coord;$i++){
    if($i==0) {print OUTmol "\$molecule\n"; print OUTmol "  $charge   $multiplicity \n";}
    print OUTmol $atom_coord[$i];
    if($i==$#atom_coord) {print OUTmol "\$end\n";}
}
close(OUTmol);
close(OUTxyz);

