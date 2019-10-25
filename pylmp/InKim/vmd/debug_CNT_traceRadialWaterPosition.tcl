# this script aligns principal axes of CNT to z axes
# and counts water molecules inside CNT
# version 20140211

### global
source "/qcfs/noische/scripts/vmd/mathtools.tcl"
source "/qcfs/noische/scripts/vmd/vmdtools.tcl"
source "/qcfs/noische/scripts/pbcdist.tcl"
source "/qcfs/noische/scripts/progress_bar.tcl"

### initialization
# description
puts "Debug) debug_CNT_traceRadialWaterPosition.tcl: trace radial water positions confined in BNNT."
puts "Debug) Version Feb 16 2016 by inkim"

# frame
# axis of CNT will be determined by the end ring atoms.
set cnt [atomselect top "resname BNT or resname CNT and z < 10"]
set target [atomselect top "type OW"]
set n_frames [molinfo top get numframes]

# random draw for water resid
set nwater [[atomselect top "type OW and water"] num]
set n_draw 100
set rand_resid ""
for { set n 1 } { $n <= $n_draw } { incr n } {
    lappend rand_resid [expr "int(rand()*$nwater)"]
}
puts $rand_resid

# output file
set filename [lindex [molinfo top get filename] 0 1]
append filename ".rWP.vmd.debug.$n_draw"
puts "Debug) Result will be saved at: $filename"
set outfile [open $filename w]

# pbc wrap
puts "Debug) Wrapping up the coordinates.."
#pbc wrap -all -compound resid

# progress bar
puts "Debug) *** Trace the position of target water ***"
progress_init $n_frames

puts $outfile "t $rand_resid"

### count the number of confined water molecules
for { set fr 1 } { $fr <= $n_frames } { incr fr } {
	# load molecules
	molinfo top set frame $fr
	$cnt frame $fr
    $cnt update
	$target frame $fr
	
	progress_tick $fr

	# projected center of mass of cnt
	set cm [center_of_mass $cnt]
	set cx [lindex $cm 0]
	set cy [lindex $cm 1]

	# water position
    set dist ""
    for { set n 0 } { $n < $n_draw } { incr n } {
        set i [lindex $rand_resid $n]
        set seltext "resid $i and type OW"
        set target [atomselect top $seltext]
	    $target update
        set coord [$target get {x y z}]
        set x [lindex [lindex $coord 0] 0]
        set y [lindex [lindex $coord 0] 1]
        lappend dist [expr "sqrt(($cx-$x)**2+($cy-$y)**2)"]
    }

    puts $outfile "$fr $dist"
	flush $outfile
}
puts " Done."
flush stdout
close $outfile
