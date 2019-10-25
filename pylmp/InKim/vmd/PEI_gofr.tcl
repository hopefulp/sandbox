# this scripts calculates average Rg and center over half of the NPT equilibration of PEI,
# and saves gofr of O_3 and O_3F.
# version 20140211

### global
source "mathtools.tcl"

### initialization
# frame
set sel [atomselect top "not water"]
set n_frames [molinfo top get frame]
set half [expr "$n_frames/2+1"]
# variables and output filename
set result {}
set filename [lindex [molinfo top get filename] 0 0]
append filename ".gofr.O_3-O_3F.profile"
puts "** result will be saved at: $filename"
set l_rg {}; set l_cx {}; set l_cy {}; set l_cz {};


### loop over frames
for { set i $half } { $i <= $n_frames } { incr i } {
	# radius of gyration
	set rg [measure rgyr $sel]; lappend l_rg $rg

	# center of mass
	set b [measure center $sel]
	set cx [lindex $b 0]; lappend l_cx $cx
	set cy [lindex $b 1]; lappend l_cy $cy
	set cz [lindex $b 2]; lappend l_cz $cz

}
set avg_rg [getAvg $l_rg]
set avg_cx [getAvg $l_cx]
set avg_cy [getAvg $l_cy]
set avg_cz [getAvg $l_cz]


### draw gofr
# selection
set inside [atomselect top "type O_3 and (x-$avg_cx)**2+(y-$avg_cy)**2+(z-$avg_cz)<$avg_rg**2"]
set outside [atomselect top "type O_3 and (x-$avg_cx)**2+(y-$avg_cy)**2+(z-$avg_cz)>$avg_rg**2"]
set o [atomselect top "type O_3F"]
set gofr_inside [measure gofr $inside $o delta 0.1 rmax 10.0 usepbc 1 selupdate 1 first $half last -1 step 1]
set ri [lindex $gofr_inside 0]
set gri [lindex $gofr_inside 1]
set igri [lindex $gofr_inside 2]

set gofr_outside [measure gofr $outside $o delta 0.1 rmax 10.0 usepbc 1 selupdate 1 first $half last -1 step 1]
set ro [lindex $gofr_outside 0]
set gro [lindex $gofr_outside 1]
set igro [lindex $gofr_outside 2]

# write
set outfile [open $filename w]
puts $outfile "r\tgr_i\tigr_i\tgr_o\tigr_o"
foreach j $ri k $gri l $igri m $gro n $igro {
	puts $outfile [format "%.4f\t%.4f\t%.4f\t%.4f\t%.4f" $j $k $l $m $n]
}
close $outfile

### END OF FUNCTION ###
