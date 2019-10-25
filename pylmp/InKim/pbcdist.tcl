#!/usr/bin/tclsh
proc pbcdist {a b dim} {
	foreach i $a j $b d $dim {
		set delta [expr abs($i-$j)]
		set delta [expr {( $delta > 0.5*$d ) ?  $d - $delta : $delta}]
		lappend result [expr $delta**2]
	}

	set sum 0.0
	foreach i $result {
		set sum [expr "$sum + $i"]
	}

	return [expr sqrt($sum)]
}

#set a {1.0 0.0 0.0}
#set b {4.0 4.0 0.0}
#set dim [list 5.0 5.0 10.0]
#set dim {5.0 5.0 10.0}

#set cx 1.0
#set cy 2.0
#set x 2.0
#set y 2.0
#set dx 0.7
#set dy 0.7
#
#set a [list $cx $cy]
#set b [list $x $y]
#set dim [list $dx $dy]
#
#puts [pbcdist $a $b $dim]
