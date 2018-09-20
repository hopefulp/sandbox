#!/usr/bin/perl

$f = shift;

if (defined($f)) {
  open(Q, $f) || die("File $f not found\n");

  $inside_mol = 0;
  @mol = ();
  @tot = ();
  while ($l = <Q>) {
		  $l =~ tr/a-z/A-Z/;
		  if ($l =~ /\$MOLECULE/) {
				  $inside_mol = 1;
		  } elsif ( ($inside_mol==1) && ($l =~ /\$END/) ) {
				  $inside_mol = 0;
				  last;
		  } elsif ($inside_mol == 1) {
				  push(@mol, $l);
		  }
		  push(@tot, $l);
  }
  close(Q);

} else {
  $inside_mol = 0;
  @mol = ();
  @tot = ();
  while ($l = <STDIN>) {
		  $l =~ tr/a-z/A-Z/;
		  if ($l =~ /\$MOLECULE/) {
				  $inside_mol = 1;
		  } elsif ( ($inside_mol==1) && ($l =~ /\$END/) ) {
				  $inside_mol = 0;
				  last;
		  } elsif ($inside_mol == 1) {
				  push(@mol, $l);
		  }
		  push(@tot, $l);
  }

}

#if $molecule section was not found...
if ($#mol < 0) {
#	print STDERR "Using the whole file...\n";
	@mol = ("   \n", @tot);	# don't forget about charge/mult line
}

@tmp = @mol;
@mol = ();

shift(@tmp);            # remove charge/multiplicity
while ( ($l = shift(@tmp)) && ($l ne "\n") ) {
        push(@mol, $l);
}

$nat = @mol;
while ( $l = shift(@tmp) ) {
        $l =~ s/=/ /g;
        $l =~ s/,/ /g;
        @dat = split(/\s+/, $l);
        ($dat[0] eq "") && shift(@dat);
        for ($i = 0; $i < $nat; $i++) {
                $mol[$i] =~ s/$dat[0]/$dat[1]/;
        }
}
 
$n = 1;
print("HEADER    $f\n");
while ( $l = shift(@mol) ) {
        @dat = split(/\s+/, $l);
        ($dat[0] eq "") && shift(@dat);
        $atom = shift(@dat);
        $x = shift(@dat) * 1.0;
        $y = shift(@dat) * 1.0;
        $z = shift(@dat) * 1.0;
		$res = shift(@dat);
		$atom2 = shift(@dat);
		(defined($atom2)) && ($atom = $atom2);
		$resnum = shift(@dat);
		($resnum > 0) || ($resnum = 1);
		$atom = mk_atom_name($atom);
 
        if (defined($res)) {
		  $res =~ s/!//;
		  $res = substr($res, 0, 3);
		  $type = res_check($res);
		  printf("%-6s %4d %-4s %3s  %4d     %7.3f %7.3f %7.3f  1.00  1.00\n",
                        $type, $n, $atom, $res, $resnum, $x, $y, $z);
		} else {
		  printf("HETATM %4d%-3s    2     2     %7.3f %7.3f %7.3f            \n",
                        $n, $atom, $x, $y, $z);
		}
        $n++;
}
print("END\n");

sub res_check {
  my $res = $_[0];
  my @knowns =
    ("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
	 "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
	 "ACE", "NME");
  my $i;

  for ($i=0; $i<=@knowns; $i++) {
	($res eq $knowns[$i]) && return("ATOM");
  }
  return("HETATM");
}

sub mk_atom_name {
  my $atom = $_[0];
  my $ret;
  if (substr($atom, 0, 1) =~ /\d/) {
#	  lead by numeric
	  $ret = $atom;
  } elsif (substr($atom, 0, 1) eq "H") {
#	is it 4-character wide? (true: just return, false: put a leading blank)
	($atom eq "HH31") && ($atom = "1HH3");
	($atom eq "HH32") && ($atom = "2HH3");
	($atom eq "HH33") && ($atom = "3HH3");
	($atom eq "HB1") && ($atom = "1HB");
	($atom eq "HB2") && ($atom = "2HB");
	($atom eq "HB3") && ($atom = "3HB");
	if (defined(substr($atom, 4))) {
#	  4-character wide
	  $ret = $atom;
	} elsif (substr($atom,0,1) =~ /\d/) {
#	  lead by numeric
	  $ret = $atom;
	} else {
	  $ret = sprintf(" %s", $atom);
	}
  } else {
#   non-hydrogen
#	is it 4-character wide? (true: just return, false: put a leading blank)
	$ret = (defined(substr($atom, 4))) ? $atom : sprintf(" %s", $atom);
  }
}
