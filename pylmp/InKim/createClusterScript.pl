#!/usr/bin/perl -w

use FindBin qw($Bin);
use lib "$FindBin::Bin";
use strict;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use Packages::General qw(FileTester);
use Cwd;
use Cwd 'abs_path';

sub init;
sub createScriptFile;
sub showUsage;

my ($OPTS, $scriptFile, $reaxFF);

$|++;
&init;
print "Create Cluster script file $scriptFile...";
&createScriptFile($OPTS, $scriptFile);
print "Done\n";

sub createScriptFile {
    my ($opts, $fileName) = @_;
    my ($currDir) = &getcwd();
    my ($prefix, $bgffile, $USER, $forcefield, $inFile, $dataFile);
    
    ($prefix, $bgffile, $forcefield, $inFile, $dataFile) = 
	($opts->{p}, abs_path($opts->{b}), $opts->{f}, $opts->{i}, $opts->{d});
    $USER = $ENV{USER};

    open SCRIPTFILE, "> $fileName" or die "ERROR: Cannot write to $fileName: $!\n";
    print SCRIPTFILE <<DATA;
#PBS -l nodes=1:ppn=12
#PBS -l walltime=480:00:00
#!/bin/csh
#PBS -N ${prefix}

echo "Nodes:"
cat \$PBS_NODEFILE
echo "JOB: MD Simulation of ${prefix}"
echo

set prefix = $prefix
set curr_dir = $currDir
set temp_dir = /scratch/${USER}/${prefix}
set bgf_file = ${bgffile}
set force_field = "${forcefield}"
set input_file = ${inFile}
set data_file = ${dataFile}
set twoPT_file = in.\${prefix}.2pt
set atom_eng_file = in.\${prefix}.atom.eng
set grp_file = \${curr_dir}/\${prefix}.grps
set nodes = 1
set nprocs = `wc -l < \$PBS_NODEFILE`
set outputfiles = "*"
set do_2pt = 0 # change this flag to do the 2pt analysis

if(\$do_2pt) then
  set create2PT = "csh -f $Bin/createVacInp.csh"
  set do2PT = /home/noische/program/2pt/2pt_analysis
  set tstep = 0.002 #CHANGE!!! default to 2fs
  set rotational_symmetry = 2 # **IMPORTANT** change this value to match your system if doing 2PT. Defaults to 2 for water!!
  if (-e \$bgf_file) then
    if ! (-e \${grp_file}) then
      $Bin/bgf2VACgrp.pl -b \$bgf_file -f resname -c 1 -o 1 -s \$grp_file >& /dev/null
    endif
  endif
  if ! (-e \${grp_file}) then
    set grp_file = `$Bin/getDOF.pl -b \$bgf_file -v 0`
    set grp_file = "\${grp_file} \${rotational_symmetry}"
  endif
  if (-e \${curr_dir}/in.\${prefix}_singlepoint) then
    set shakeStr = `grep "fix             shakeH" \${curr_dir}/in.\${prefix}`
    cp \${curr_dir}/in.\${prefix}_singlepoint \${curr_dir}/_in.test
    sed -i 's/read_data.*\$/read_restart \${restart}/' \${curr_dir}/_in.test
    sed -i 's/variable.*//' \${curr_dir}/_in.test
    sed -i 's/group.*//' \${curr_dir}/_in.test
    sed -i 's/run.*//' \${curr_dir}/_in.test
    cat \${curr_dir}/_in.test /home/noische/scripts/dat/LAMMPS/in.lammps.2pt > \${curr_dir}/in.\${prefix}.2pt
    sed -i "s/fix             shakeH.*/\$shakeStr/" \${curr_dir}/in.\${prefix}.2pt
    cat \${curr_dir}/_in.test /home/noische/scripts/dat/LAMMPS/in.lammps.atom.eng > \${curr_dir}/in.\${prefix}.atom.eng
    sed -i "s/fix             shakeH.*/\$shakeStr/" \${curr_dir}/in.\${prefix}.atom.eng
    sed -i -r 's/(ewald.*|pppm.*)/ewald 0.001/' \${curr_dir}/in.\${prefix}.atom.eng
    rm -fr \${curr_dir}/_in.test
  endif
endif

setenv LAMMPS_PARALLEL "/opt/mpi/intel/mpich2-1.4rc2/bin/mpirun -np $nprocs -hostfile \$PBS_NODEFILE"
set executable = /opt/applic/lammps/bin/lmp_kdft-mpich2-1.4-19Jul

echo Running \$executable on \$nodes nodes \$nprocs processors launched from \$PBS_O_HOST
echo "Nodes:"
cat \$PBS_NODEFILE

cat \$PBS_NODEFILE > \$curr_dir/nodelist.\$PBS_JOBID

mkdir -p \$temp_dir/analysis \${curr_dir}/results
cd \$temp_dir
foreach i (\$input_file \$data_file \$twoPT_file \$atom_eng_file)
  if (-e \${curr_dir}/\${i}) cp \${curr_dir}/\${i} ./
end

$reaxFF

echo "LAMMPS equilibration dynamics of \${prefix}"
\$LAMMPS_PARALLEL \$executable  -in \$input_file -log \${prefix}.equil.lammps.log || goto error

if !(-e \${twoPT_file}) set do_2pt = 0
if (\$do_2pt == 1) then
  foreach i (\${prefix}.npt.*.restart)
    cd \${temp_dir}
    set sname = `basename \$i .restart`
    echo "2PT analysis using \$sname"
    \$LAMMPS_PARALLEL \$executable  -in \$twoPT_file -var restart \$i -var sname \$sname || continue
    foreach j (\${sname}.2pt.*.restart)
      \$LAMMPS_PARALLEL \$executable  -in \$atom_eng_file -var restart \$j -var sname \$sname >& /dev/null
      rm -fr \$j
    end
    cd \${temp_dir}/analysis
      \${create2PT} ../\${sname}.2pt ../\${data_file} 3 \${tstep} \${grp_file}
      \${do2PT} \${sname}.2pt.mol.grps.in >& \${sname}.2pt.mol.grps.screen.out
      cp \${sname}.2pt.mol.grps.* \${curr_dir}/results
    cd ../
    rm -fr ./\${sname}.2pt.lammps
  end
endif

exit:
echo " "
echo No errors detected
echo " "
exit(0)

error:
echo " "
echo ERROR OCCURRED: Check output file
echo " "
exit(1)

DATA

close SCRIPTFILE;
}

sub init {
    my ($i);
    getopt('pbfidsr', \%{ $OPTS });
    for $i ("p", "b", "f", "i", "d") {
	die &showUsage . "\n" if (! exists($OPTS->{$i}));
    }
    print "Initializing...";
    for $i ("b", "f", "i", "d") {
	if ($i eq "f") {
	    if ($OPTS->{$i} =~ /\s/) {
		while ($OPTS->{$i} =~ /(\S+)/g) {
		    FileTester($1);
		}
	    } else {
		FileTester($OPTS->{$i});
	    }
	} else {
	    FileTester($OPTS->{$i});
	}
    }
    print "Done\n";
    $scriptFile = $OPTS->{s};
    $reaxFF = "cp $OPTS->{r} ./ffield.reax" if (exists($OPTS->{r}));
    if (! defined($scriptFile)) {
	$scriptFile = basename($OPTS->{b});
	$scriptFile =~ s/\.\w+$//;
	$scriptFile .= "_cluster.script";
    }
    print "Done\n";
}

sub showUsage {
    return "usage: $0 -p name of job -b bgf file -f forcefield file -i lammps input file -d lammps data file\n" .
	"Options\n\tprefix: This will be the name of the job. The files will be stored in /temp1/{USER}/{name of job}\n";
}
