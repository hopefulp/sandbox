### global
set n_frames [molinfo top get numframes]

### initialization
set filename [lindex [molinfo top get filename] 0 1]
append filename ".water.profile"
puts "result will be saved at: $filename"
set outfile [open $filename w]

### average
proc getAvg { rList } {
	set n [llength $rList]
	if {!$n} {return 0.0}

	set avg 0.0
	foreach r $rList {
		set avg [expr $avg + $r]
	}
	set avg [expr $avg/double($n)]

	return $avg
}

### main loop
for { set fr 1 } { $fr <= $n_frames } { incr fr } {
	# load molecules
	molinfo top set frame $fr
	set cnt [atomselect top "resname CNT"]
	
	# center of mass of cnt
	set cm [measure center $cnt]
	set cx [lindex $cm 0]
	set cy [lindex $cm 1]

	# radius of cnt
	set radius 0.0
	set zmin 10000.0
	set zmax -10000.0
	set dist_sq {}
	foreach i [$cnt get {x y z}] {
		set x [lindex $i 0]
		set y [lindex $i 1]
		set z [lindex $i 2]

		# CNT radius
		lappend dist_sq [expr "($cx-$x)**2+($cy-$y)**2"]
		set avg_dist_sq [getAvg $dist_sq]
		unset dist_sq

		# CNT length
		if {$z < $zmin} {
			set zmin $z
		}
		if {$z > $zmax} {
			set zmax $z
		}
	}

	# select atoms in cnt
	set inside [atomselect top "(x-$cx)**2+(y-$cy)**2 < $avg_dist_sq and z > $zmin and z < $zmax and type OW"]
	puts $outfile "$fr [$inside num]"
	puts -nonewline "\rprocessing $fr / $n_frames"
	flush stdout
}
puts ""
close $outfile
