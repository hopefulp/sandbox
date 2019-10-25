### global
set sel [atomselect top "not water"]
set n_frames [molinfo top get frame]


### initialization
set result {}
set filename [lindex [molinfo top get filename] 0 0]
append filename ".atom.profile"
puts "result will be saved at: $filename"


### main loop
for { set i 0 } { $i <= $n_frames } { incr i } {
	#puts "frame $i"
	$sel frame $i

	# radius of gyration
	set rg [format %f [measure rgyr $sel]]
	#puts "Rg: $rg"

	# center of mass
	set b [measure center $sel]
	#puts "CM: $b"
	set cx [format %f [lindex $b 0]]
	set cy [format %f [lindex $b 1]]
	set cz [format %f [lindex $b 2]]

	# select inside Rg
	set inside [atomselect top "(x-$cx)**2+(y-$cy)**2+(z-$cz)**2 < $rg**2 and type O_3"]
	set outside [atomselect top "(x-$cx)**2+(y-$cy)**2+(z-$cz)**2 > $rg**2 and type O_3"]
	set all [atomselect top "type O_3"]

	# count number of atoms
	set n_inside [$inside num]
	set n_outside [$outside num]
	set n_all [$all num]
	set check [expr $n_inside + $n_outside]

	lappend result "$i\t$n_all\t$n_inside\t$n_outside\t$cx\t$cy\t$cz\t$rg"
}


### output
puts "t\tall\tinside\toutside\tCx\tCy\tCz\tRg"
foreach x $result {
	puts $x
}


### write output to file
set outfile [open $filename w]
puts $outfile "t\tall\tinside\toutside\tCx\tCy\tCz\tRg"
foreach x $result {
	puts $outfile "$x"
}
close $outfile


### END OF FUNCTION ###
