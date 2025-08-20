#!/usr/bin/python2.7 -u
# written by Joonho Park
# ? read  DOSCAR or dir/DOSCAR if there is directory name
# cd to inside of directory and run ./dosall.pl
### usage: ./dosall.pl a=1 l=1 m=0 s=1 or -1
# read l and m quantum number
# for spin up or down s=[1|0] 1 for up and 0 for down
# for atom  a=25 or a=1:25

suffix=""
#use lib '/home/joonho/sandbox/mod_perl';

#use GetHash;

#use VaspDos;


#$Fermi_shift='T';	# T for shift, F for no-shift

#$Fermi_check=0;		# write Fermi.dat at value EF with value 1.2 * max(Tdos)

#$Max_Tdos=0;		# control of y-axis of Dos.dat


#print "Now run in the directory\n";

########## argument analysis 

if len(sys.argv) < 0:
   print Usage:: 0 [DOSCAR] for Tdos and Fermi.dat 
   print Usage:: 0 a=na[:nb] l=l m=m s=[1 for up|-1 for down]
   print "        below a are options"
#    exit(1);


fin="DOSCAR"
if ARGV[0] =~ //^d//:
   @dir = ARGV[0].split("=")
       fin="dir[1]//DOSCAR"
    shift ARGV
   print "Note:: The outfiles will be located at present directory"
#exit(1);

#for($i=0;$i <= $#ARGV; $i++){

    atom_list=""
if ARGV[i] =~ //=//:
   	(quantum,num)=&GetHash::getvalue(ARGV[i],'=')
	# usage: l=2 m=2 or just l=1 then m is summed
   if quantum == "l":
                  l=num
   elif quantum == "m":
                  m=num
   elif quantum == "s":
                  spin=num
   elif quantum == "a":
      	    atom_list=num
      if atom_list =~ //://:
         	    	    (a_start,a_end)=&GetHash::get_2values(atom_list,':')
#	    	    for($j=$a_start;$j<=$a_end;$j++){

	    	        push(a_list,$j)
   elif atom_list =~ //,//:
#                @a_list=split(', ', $atom_list)

   else:
   elif quantum == "dir":
      	    fin="num//DOSCAR"
    # if there is not non-digital (only digital), it's l or m
elif ARGV[i] !~ //\D//:
   if !defined(l):
   elif !defined(m):
    # if there is not non-word(a-Z,0-9,_) character, it's input file
   elif ARGV[i] !~ //\W//:
      if ARGV[i] == "DOSCAR":
      else:

   if !defined(fin):
      print input = fin
      sys.stdout.write("quantum number: \t")
      if defined(l):
      else:
         if defined(m):
         else:
            if defined(spin):
            else:

               print atom list: ",join(' ',a_list),"

######## ordering atoms and output file naming

#@magnetic_orbital=(['y','z','x'],['xy','yz','z2','xz','x2-y2']); VASP-DOSCAR ordering
               if defined(l):
                      suffix="_orbital[l]"
                  if defined(m):
#	$suffix.=$magnetic_orbital[$l][$m];    

                      fdos="Pdos"

               if  0 <= #a_list :
#    print "atom_list: ",$#a_list, join("  ",@a_list),"\n";
                  if #a_list == 0:
                     	suffix.="_aa_list[0]"
                  else:
                         	suffix.="_aa_list[0]-a_list[#a_list]"
                  if !defined(fdos):

                  if defined(spin):
                     if spin==1:
                        if spin==-1:
                               suffix.="_spin_l"

                        if !defined(fdos):

                           fout="fdossuffix.dat"
                           if Fermi_shift == "T":
                                  fout="fdossuffix"."F0.dat"

### TDOS is overall
### LDOS, PDOS is spin up and down divided
                           sys.stdout.write(Input: fin\t)
                           print OUT: fout

#open(IN,"< $fin");

#open(OUT,"> $fout");

                           if Fermi_shift == "F" and #a_list < 0:
#    open(EFOUT,"> Fermi.dat");


                           i=0
#$iatom=0;	# note first atom 0 is for total DOS not atom

                           iene=0
#$Nene=1000;	# a large number to skip $Nene == $iene == 0

                           tag_spin="NO"
                           delimit_atom=""
                           Tag_Atom_delimit="OFF"
###################################################  read three atoms's coordinates
                           while line=<IN>:
#    if($iene==0) {print $line;}
#    if($line eq $delimit_atom){ next LINE;}
#    chomp($line);

                              @line = line.split("\s+")
                              if line[0] == ""){shift(line:

                                 if i==0:
                                    if i==1:
                                       if i==2:
                                          if i==3:
                                             if i==4:
                                                if i==5:
                                                   	delimit_atom=line
                                                   	Emax=line[0]
                                                   	Emin=line[1]
                                                   	Nene=line[2]
                                                   	ene_Fermi=line[3]
                                                   	dnknow=line[4]
#	next LINE;

    # this part locates after $i==6  for $iene=0 and start in an atom block 
                                                if line == delimit_atom:
#   	$Tag_Atom_delimit="ON";	

#	next LINE;



### DOS LINE: 1. energy 2. tdos_@_E    3. tdos_int_E 						for non-spin
### DOS LINE: 1. energy 2. tdos_@_E_up 3. tdos_@_E_down 4. tdos_int_E_up 5. tdos_int_E_down 	for spin
    # define spin 
                                                if tag_spin == "NO":
                                                   if #line == 2:
                                                                  Lspin='F'
                                                   elif #line == 4:
                                                                  Lspin='T'
                                                   else:
                                                      	tag_spin="YES"
    # save all the components in dos line
#    print "iatom = $iatom\n";
    # $iatom == 0 is for total DOS; atoms start from 1 to # of atoms
                                                   if iatom == 0:
	push(energy,$line[0])
#    	print $iatom,"\t",$line[0],$energy[$iene],"\n";
    # save all the dos components: first TDOS then atoms
#    for($j=1;$j<=$#line;$j++){

                                                           tdos[iatom][iene][j-1]=line[j]
                                                   if iatom == 1:
                                                   iene +=1

#} continue {

                                                   i +=1
    # when you read delimiter between PDOS block, skip to next line withough anithing changed, $iene, $iatom
    # So the number of atom is same as the number of delimit
                                                   if Tag_Atom_delimit == "ON":
                                                      	Tag_Atom_delimit="OFF"
                                                      if iatom==1:
                                                         	iene=0
                                                         iatom +=1

# iatom start from 0, TDOS then, 54 atoms' PDOS, and last added +1 so it should be natom+1
                                                   if natom != iatom:
                                                      print ERROR: in counting atoms natom != $iatom
#    exit(0);

############################ END: load DOSCAR 
############################ shuffle DOS
#if(!defined($spin)){ $nspin="YES";}  # For DOSCAR, for PROCAR use "NO"
#($nc,$nd,$mc,$md)=&VaspDos::lm_quantum($l,$m,$nspin,$spin);
#### other options

                                                   Efermi=ene_Fermi

                                                   if Fermi_shift == "F":
#    $ene_Fermi = 0;		# do not shift w.r.t Fermi level

                                                          EF_test=Efermi
                                                      print "Energy is not shifted"
                                                   else:
                                                      print "Energy is shifted to Fermi Level"
                                                          EF_test=0

                                                   if Lspin == "T":
                                                          nspin=2
                                                   elif Lspin == "F":
                                                          nspin=1
                                                   else:
                                                      print "Magnetism error"
#    exit(5);


#For DOSCAR, for PROCAR use 1 , note: this is not for spinless
#column list start from 0 to 17 in the case of s(2), p(6), d(10) for spin
# For non-magnetic 
                                                   (nl_col)=&VaspDos::lm_quantum(l,m,nspin,spin)
                                                   print column list : ", join(" ",nl_col),"
#for($i=0;$i<$Nene;$i++){

                                                   if Max_Tdos < tdos[0][i][0]:
                                                      	Max_Tdos = tdos[0][i][0]
                                                print Max_Tdos = Max_Tdos
#for($i=0;$i<$Nene;$i++){

#    if($energy[$i] == 3.013) { print $tdos[0][$i][0],  $tdos[0][$i][1],  $tdos[0][$i][2],  $tdos[0][$i][3],"\n"; exit; }
                                                if i==0:
    ### This is for non-magnetic DOSCAR
                                                   if Lspin == "F":
    	#print $energy[$i]-$ene_Fermi,"\t", $tdos[0][$i][0],"\n";
    	### tdos[0] is for TDOS and others for atoms
    	### tdos[][][0] is for DOS at the energy ; tdos[][][1] is for accumulated DOS
        ### if there is no atom list, write Tdos
                                                      if #a_list < 0:
#           print OUT $energy[$i]-$ene_Fermi,"\t", $tdos[0][$i][0],"\n";

                                                         if Fermi_shift == "F" and #a_list < 0:
                                                            if Fermi_check==0 and energy[i] < EF_test:
#                   print EFOUT $energy[$i]-$ene_Fermi,"\t 0\n";

                                                            elif Fermi_check==0 and energy[i] >= EF_test:
#                   print EFOUT "$EF_test\t", $Max_Tdos*2, "\n";

                                                                                  Fermi_check=1
                                                            else:
#                   print EFOUT $energy[$i]-$ene_Fermi,"\t 0\n";

        ### For atom list, get LDOS or PDOS depending on l quantum number
                                                      else:
                                                         	        pdos=0
                                                                     pdosm=0
                                                                     ldos=0
            ### as for an energy, loop atom lists [whethe not sorted], add l columns
#            for($j=0;$j<=$#a_list;$j++){

#                for($ic=0;$ic<=$#nl_col;$ic++){

                                                                             k=nl_col[ic]
#                    $ldos += $tdos[$a_list[$j]][$i][$k];

                    #print "ldos $ldos\n";
#            print OUT   $energy[$i]-$ene_Fermi,"\t", $ldos,"\n";


    ### This is for magnetic DOSCAR
                                             else:
#	print $energy[$i]-$ene_Fermi,"\t", $tdos[0][$i][0]+$tdos[0][$i][1],"\n";
#	### if there is no atom list, write Tdos
                                                if #a_list < 0:
#	    print OUT $energy[$i]-$ene_Fermi,"\t", $tdos[0][$i][0]+$tdos[0][$i][1],"\n";

                                                   if Fermi_shift == "F" and #a_list < 0:
                                                      if Fermi_check==0 and energy[i] < EF_test:
#    	            print EFOUT $energy[$i]-$ene_Fermi,"\t 0\n";

                                                      elif Fermi_check==0 and energy[i] >= EF_test:
#	            print EFOUT "$EF_test\t", $Max_Tdos*2, "\n";

                                                         	            Fermi_check=1
                                                      else:
#	            print EFOUT $energy[$i]-$ene_Fermi,"\t 0\n";

	### For atom list, get LDOS or PDOS depending on l quantum number
                                                else:
	    # for LDOS
                                                   	    pdos=0
                                                   	    pdosm=0
                                                   	    ldos=0
	    # LDOS if atoms are selected, LDOS is spin adapted if $spin is defined
 	
	    # sum all the selected atoms	j atom index
#	    for($j=0;$j<=$#a_list;$j++){	

		#print "sum $a_list[$j] -th atom\n";
		#print join(" ",@nl_col),"\n";
		# sum all the orbitals; 	k dos colume index scans all the column
#		for($ic=0;$ic<=$#nl_col;$ic++){

                                                   		    k=nl_col[ic]
#		    $ldos += $tdos[$a_list[$j]][$i][$k];

		    #print "ldos $ldos\n";
#	    print OUT   $energy[$i]-$ene_Fermi,"\t", $ldos,"\n";

		
	 # if atoms are not selected
#	}else{
#	    # PDOS for all atoms
#	    if(defined($l)){
#	    # sum all the atoms
#	    	for($j=0;$j<=$#a_list;$j++){
#	            for($k=$nc;$k<=$nd;$k++){
#	        	$pdos += $tdos[$a_list[$j]][$i][$k];
#	        	}
#	        }
#	    }
#	print $energy[$i],"\t", $pdos,"\n";
#	print $energy[$i]-$ene_Fermi,"\t", $pdos,"\n";

#close(IN); close(OUT); 


