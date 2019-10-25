set sel [atomselect top "(within 3 of resname PEI) and resname ION"]
set n_frames [molinfo top get frame]

set result {}
set filename [lindex [molinfo top get filename] 0 0]
append filename ".ion.profile"
puts "result will be saved at: $filename"

set outfile [open $filename w]

for { set i 0 } { $i <= $n_frames } { incr i } {
	#puts "frame $i"

	# radius of gyration
	#set rg [measure rgyr $sel]
    set sel [atomselect top "(within 3 of resname PEI) and resname ION"]
	$sel frame $i
    set rg [$sel num]
	lappend result "$i\t$rg"
}

# write output to file
foreach x $result {
	puts $outfile "$x"
}

close $outfile
