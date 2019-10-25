### global
set n_frames [molinfo top get numframes]


### initialization
set result {}
set filename [lindex [molinfo top get filename] 0 1]
append filename ".end-end_dist.profile"
puts "result will be saved at: $filename"
set outfile [open $filename w]

puts $n_frames
puts "t\tend-end_dist"

puts $outfile $n_frames
puts $outfile "t\tend-end_dist"


### input: two vectors of length three (v1 and v2)
### returns: ||v2-v1||
proc vecdistsq {v1 v2} {
	lassign $v1 x1 x2 x3
	lassign $v2 y1 y2 y3
	return [expr ($x1-$y1)**2 + ($x2-$y2)**2 + ($x3-$y3)**2]
}


### main loop
for { set fr 0 } { $fr <= $n_frames } { incr fr } {
	# load molecules
	molinfo top set frame $fr
	set sel [atomselect top "not water and not hydrogen"]
	
	set coords1 [$sel get {x y z}]
	set coords2 [$sel get {x y z}]
	set max_distsq 0

	foreach i $coords1 {
		foreach j $coords2 {
			set distsq [vecdistsq $i $j]
			if {$max_distsq < $distsq} {
				set max_distsq $distsq
			}
		}
		lvarpop j
	}

	unset coords1
	unset coords2
	unset sel
	
	### output
	set max_dist [expr sqrt($max_distsq)]
	puts "$fr $max_dist"
	puts $outfile "$fr $max_dist"
}

close $outfile

### END OF FUNCTION ###
