#!/usr/bin/perl -w

use strict;
use File::Basename qw(basename);
use Getopt::Std qw(getopt);

sub getFluidicity;
sub init;
sub doAvg;
sub getStats;
sub doStats;
sub writeData;

my ($FILES, $GROUPS, $saveName);
my ($DATA);

$|++;
($FILES, $GROUPS) = &init;
$DATA = getFluidicity($FILES, $GROUPS);
print "Collecting statistics...";
&doStats($DATA);
print "Writing stats to $saveName...";
&writeData($DATA, $saveName);
print "Done\n";

sub writeData {
  my ($data, $saveName) = @_;
  my ($i, $j, $k, $sKey, $grp, $tmp, @gList);

  open OUTFILE, "> $saveName" or die "ERROR: Cannot write to $saveName:$!\n";
  $sKey = 0;
  for $i (keys %{ $data }) {
    $sKey = length($i) if (length($i) > $sKey);
  }
  $sKey += 5;
  printf OUTFILE "%-${sKey}s", "name";
  $grp = $tmp = sprintf("%-${sKey}s","");
  for $i (keys %{ $data }) {
    for $j (sort { $a<=>$b } keys %{ $data->{$i} }) {
      printf OUTFILE "%32s", "group $j";
      $tmp .= sprintf("%16s%16s","Trans", "Rot");
      $grp .= sprintf("%8s%8s%8s%8s","avg", "+/-","avg", "+/-");
      push @gList, $j;
    }
    print OUTFILE "\n$tmp\n$grp\n";
    last;
  }
  for $i (sort { $a cmp $b } keys %{ $data }) {
    printf OUTFILE "%-${sKey}s", $i;
    for $j (@gList) {
      for $k ("TRANS", "ROT") {
	printf OUTFILE "%8.3f%8.3f",$data->{$i}{$j}{STATS}{$k}{AVG}, $data->{$i}{$j}{STATS}{$k}{STDEV};
      }
    }
    print OUTFILE "\n";
  }

  close OUTFILE;
}

sub doStats {
   my ($data) = $_[0];
   my ($i, $j, $k);

   for $i (keys %{ $data }) {
    for $j (keys %{ $data->{$i} }) {
      for $k ("TRANS", "ROT") {
	$data->{$i}{$j}{STATS}{$k} = getStats($data->{$i}{$j}{DATA}{$k});
      }
    }
  }
}

sub getFluidicity {
  my ($files, $grps) = @_;
  my ($i, $j, $tmp, $pStr, $tot, $data, $sKey, $index, $pLen, $lineD);

  $tot = scalar(@{ $files });
  for $i (0 .. $#{ $files }) {
    $pStr = "Parsing $files->[$i] (". ($i+1) . ") of $tot...";
    print "$pStr\r";
    @{ $tmp } = `grep "stlin fludicity" $files->[$i]`;
    if (@{ $tmp }) {
      $sKey = basename($files->[$i]);
      $sKey =~ s/\..*$//;
      for $j (@{ $grps }) {
	$index = $j-1;
	last if ($index > $#{ $tmp });
	chop $tmp->[$index];
	@{ $lineD } = split /\s+/, $tmp->[$index];
	push @{ $data->{$sKey}{$j}{DATA}{TRANS} }, $lineD->[2];
	push @{ $data->{$sKey}{$j}{DATA}{ROT} }, $lineD->[3];
	push @{ $data->{$sKey}{$j}{DATA}{FILE} }, $files->[$i];
      }
    }
  }
  die "ERROR: No valid data found!\n" if (! defined($data));
  $pLen = length($pStr);
  printf "%${pLen}s...Done\n",$pStr;
  
  return $data;
}

sub getStats {
    my ($data, $start) = ($_[0], $_[1]);
    my (%stats, $i, $numPoints, $sigma, $x1, $x2);

    return () if (! @{ $data });

    $start = 1 if (! $start);
    $start--;
    $numPoints = 0;
    for $i ($start .. $#{ $data }) {
        $x1 += $data->[$i];
        $x2 += ($data->[$i]*$data->[$i]);
        $numPoints++;
    }

    $x1 /= $numPoints;
    $x2 /= $numPoints;
    $sigma = ($x2 - $x1*$x1);
    if ($sigma > 0.00001)  {
        $sigma = sqrt($sigma);
    } else {
        $sigma = 0;
    }
    $stats{STDEV} = $sigma;
    $stats{AVG} = $x1;
    $stats{NUM} = $numPoints;
    $stats{TOT} = $stats{AVG} * $stats{NUM};

    return \%stats;
}

sub init {
  my (%OPTS, $grpStr, $fStr, $fList, $grpList, $i);

  getopt('fsg',\%OPTS);

  die "usage: $0 -f file(s) -g (group(s)) -s (savename)\n" if (! exists($OPTS{f}));
 
  print "Initializing...";
  ($fStr, $grpStr, $saveName) = ($OPTS{f}, $OPTS{g}, $OPTS{s});
  @{ $fList } = `ls $fStr`;
  $i = 0;
  while ($i <= $#{ $fList }) {
    chop $fList->[$i];
    if (! -e $fList->[$i] or ! -r $fList->[$i] or ! -T $fList->[$i]) {
      splice @{ $fList }, $i, 1;
      next;
    }
    $i++;
  }
  die "ERROR: No files found while searching \"$fStr\"!\n" if (! @{ $fList });
  $grpStr = "" if (! defined($grpStr));
  @{ $grpList } = ();
  while ($grpStr =~ /(\d+)/g) {
    push @{ $grpList }, $1 if ($1>0);
  }
  push @{ $grpList }, 1 if (! @{ $grpList });

  if (! defined($saveName)) {
    $saveName = basename($fList->[0]);
    $saveName =~ s/\.\w+$//;
    $saveName .= "_fluidicity_stats.dat";
  }
  print "will select fluidicity from group @{ $grpList } from $#{ $fList } files...Done\n";
  return ($fList, $grpList);
}
