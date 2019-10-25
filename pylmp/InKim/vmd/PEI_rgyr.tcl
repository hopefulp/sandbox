set sel [atomselect top "not water and not type C_2G"]
set n_frames [molinfo top get frame]

set result {}
set filename [lindex [molinfo top get filename] 0 0]
append filename ".rg.profile"
puts "result will be saved at: $filename"

set outfile [open $filename w]

for { set i 0 } { $i <= $n_frames } { incr i } {
	#puts "frame $i"
	$sel frame $i

	# radius of gyration
	set rg [measure rgyr $sel]
	lappend result "$i\t$rg"
}

# write output to file
foreach x $result {
	puts $outfile "$x"
}

close $outfile
