#!/usr/bin/perl
# written by Joonho Park
# read a.inp and write a.out

if($#ARGV < 0){
    print "Usage:: $0 qchem.out [job-id]\n";
    print "Job-id default is 1\n";
    exit(0);
}    

$fname=$ARGV[0];
@fname=split(/\./,$fname);
$fname_pre=$fname[0];
$fout="$fname_pre.mol";
$fxyz="$fname_pre.xyz";
print "Overwrite mol file: $fout and $fxyz\n";

open(IN,"<$fname");
open(OUTmol,">$fout");
open(OUTxyz,">$fxyz");

### using number of Jobs or job type
if($ARGV[1] eq ""){
    $njob=1;
    $job_type="OPT"
}else{
    $njob=$ARGV[1];
    $job_type=$ARGV[1];
}

$keyword_mol="molecule";
$endword_mol="end";
$flag_mol="OFF";
$tag_mol="ON";

$job_tag="OFF";
#$keyword_job="Job";
$keyword_job="rem";

$keyword_opt="Optimization Cycle";
$endword_opt="OPTIMIZATION CONVERGED";
$tag_opt="ON";
$flag_opt="OFF";

#$tag_geometry="GEOMETRIES";	# for MOLDEN output
$keyword_ext="Step";
$flag_coord="OFF";

$i=0;
$iatom=0;
$ijob=0;
$iopt=0;

### read three atoms's coordinates
while($line=<IN>){
    $i++;
    ### extract charge and spin
    ### if for block execution
    if($tag_mol eq "ON"){          # job for molecule with spin and charge
        # for block reading
        if($flag_mol eq "OFF" ){       # find mol section 
            if($line =~ /\$$keyword_mol/ or $line =~ /$keyword_job/){ 
                if($line =~ /\$$keyword_mol/){    
                    print "found \$$keyword_mol\n";
                    $flag_mol = "ON";      # found mol section and run in the block 
                }
                # if rem section appears first, increase job id
                if($line =~ /$keyword_job/){ 
                    $ijob++;
                }
            }
            next;
        }else{
            # line analysis in the block
            @field=split(/\s+/,$line);
            if($field[0] eq "") {shift(@field);}
            if($iatom==0){
                $charge=$field[0]; $multiplicity=$field[1];
                $iatom++;
            # arrive at the end
            }elsif($line =~ /$endword_mol/){
                $flag_mol = "OFF";
                print $charge, "***",$multiplicity,"\n";
                $tag_mol = "Done";     # finish job for molecule section
                $natom=$iatom-1;        # save num of atoms
                print "num of atom = ",$natom,"\n"; 
            #### just count num of mol if not blank
            }else{
                #### skip if blank line in mol section
                if($field[0] =~ /\w/) {  $iatom++; }
            }
            next;
        }
    }

    ### find n-th job here
    if($job_tag eq "OFF"){
        if($ijob == $njob){
            print "$njob-th Job\n";
            $job_tag="ON";          # start work on finding opt mol
        }elsif($line =~ /$keyword_job/){
            $ijob++;
            #@field=split(/\s+/,$line);
            #if($field[0] eq "") {shift(@field);}
        }
        next;
    }
    ### if job_tag is "ON" for search opt block, reaches here
    ### count num of optimization step
    if($tag_opt eq "ON"){
        ### find "Optimization Cycle" or finish keyword
        if($flag_opt eq "OFF"){
            if($line =~ /$keyword_opt/) {
                if($iopt == 0) {print "Optimization job\n";}
                @field=split(/\s+/,$line);
                if($field[0] eq "") {shift(@field);}
                $iopt=$field[2];
                $flag_opt = "ON";
                # empty coordinate lines and start saving new coord's
                @atom_coord_line=();
            }elsif($line =~ /$endword_opt/){
                $tag_opt = "Done";                      # job is done
                print "Optimization step: $iopt\n";
                $n_opt=$iopt;
                last;
            }
            next;
        ### if coord block starts,
        }else{
            @field=split(/\s+/,$line);
            if($field[0] eq "") {shift(@field);}
            if($field[0] =~ /\d/ and (1<=$field[0] and $field[0]<=$natom) ){
                $flag_coord="ON";
                push(@atom_coord_line,$line);
            }elsif($flag_coord eq "ON"){
                $flag_coord= "OFF";
                $flag_opt  = "OFF";         # saving coord block is done
            }
        }
    }
}
close(IN);

if($tag_opt eq "ON"){
    print "Optimization fails: $iopt\n";
    $n_opt=$iopt;
}    

#print @atom_coord_line;
@atom_coord=();
for($i=0;$i<=$#atom_coord_line;$i++){
    #print $atom_coord_line[$i];
    @field=split(/\s+/,$atom_coord_line[$i]);
    if($field[0] eq "") {shift(@field);}
    shift(@field);
    $newline=join("\t", @field);
    push(@atom_coord,$newline);
}    

#print join("\n", @atom_coord);
#print "\n";

### xyz file
for($i=0;$i<=$#atom_coord;$i++){
    if($i==0) {
        print OUTxyz "$natom\n"; 
        print OUTxyz "$keyword_ext $n_opt\t";
        if($tag_opt eq "ON"){   print OUTxyz "not converged \n";}
        else{                   print OUTxyz "\n";}
    }
    #print OUTxyz $atom_coord[$i], "\n";
    @field=split(/\s+/,$atom_coord[$i]);
    if($field[0] eq "") {shift(@field);}
    printf OUTxyz "%5s%10.5f%10.5f%10.5f\n", $field[0], $field[1], $field[2], $field[3];
    printf "%5s%10.5f%10.5f%10.5f\n", $field[0], $field[1], $field[2], $field[3];
}

### mol file
for($i=0;$i<=$#atom_coord;$i++){
    if($i==0) {print OUTmol "\$molecule\n"; print OUTmol "  $charge   $multiplicity \n";}
    print OUTmol $atom_coord[$i], "\n";
    if($i==$#atom_coord) {print OUTmol "\$end\n";}
}
close(OUTmol);
close(OUTxyz);

