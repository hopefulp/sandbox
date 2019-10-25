### CNT_countWaterInside.tcl: counts water molecules near CNT walls. 20151005
### global
set n_frames [molinfo top get numframes]

### initialization
set filename [lindex [molinfo top get filename] 0 1]
append filename ".infinite_nanotube_water.profile"
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

puts $outfile "t\tinner\touter"
### main loop
for { set fr 0 } { $fr <= $n_frames } { incr fr } {
	# load molecules
	molinfo top set frame $fr
    set outer [atomselect top "(water) and (within 5 of resname CNT) and (type OW)"]
    set inner [atomselect top "(water) and (not within 5 of resname CNT) and (type OW)"]

	# select atoms in cnt
	puts $outfile "$fr [$inner num] [$outer num]"
	puts -nonewline "\rprocessing $fr / $n_frames"
	flush stdout
}
puts ""
close $outfile
