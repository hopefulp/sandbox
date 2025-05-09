# this script aligns principal axes of CNT to z axes
# and counts water molecules inside CNT
# version 20140211

### global
source "/qcfs/noische/scripts/vmd/mathtools.tcl"
source "/qcfs/noische/scripts/vmd/vmdtools.tcl"
source "/qcfs/noische/scripts/pbcdist.tcl"
source "/qcfs/noische/scripts/progress_bar.tcl"
package require Orient
namespace import Orient::orient
source "/qcfs/noische/scripts/vmd/plugins/orient/orient.tcl"


### initialization
# description
puts "Info) BNNT_countWaterComplex.tcl: count the number of waters confined in BNNT."

# frame
set cnt [atomselect top "resname CNT"]
set all [atomselect top "all"]
set n_frames [molinfo top get numframes]

# output file
#set filename [lindex [molinfo top get filename] 0 1]
#append filename ".countWater.vmd.v3.profile"
#puts "Info) Result will be saved at: $filename"
#set outfile [open $filename w]

# progress bar
puts "Info) *** Sweep 1: Wrapping up the coordinates ***"
#progress_init $n_frames


### rotate the whole system by principal axes of BNNT
for { set fr 1 } { $fr <= $n_frames } { incr fr } {
	# load molecules
	molinfo top set frame $fr
	$cnt frame $fr
	$all frame $fr
	$cnt update
	$all update
	
	#progress_tick $fr

	# center of mass of cnt
	set cm [center_of_mass $cnt]
	set cx [lindex $cm 0]
	set cy [lindex $cm 1]
	set cz [lindex $cm 2]

	# get pbc
	set dim [pbc get]
	set dimx [expr [lindex [lindex $dim 0] 0]/2.0]
	set dimy [expr [lindex [lindex $dim 0] 1]/2.0]
	set dimz [expr [lindex [lindex $dim 0] 2]/2.0]


	# move com of cnt to boxcenter
	moveback $all $cm
	$all moveby "$dimx $dimy $dimz"
}

# pbc wrap
puts "Info) Wrapping up the coordinates.."
pbc wrap -all -compound resid

# progress bar
puts "Info) *** Sweep 2: Rotating the NT and count the number ***"
#progress_init $n_frames

#puts $outfile "fr\tn_wat\tr\tcx\tcy\tzmin\tzmax"


### count the number of confined water molecules
for { set fr 1 } { $fr <= $n_frames } { incr fr } {
	# load molecules
	molinfo top set frame $fr
	$cnt frame $fr
	$all frame $fr
	$cnt update
	$all update
	
	#progress_tick $fr

	# get pbc
	set dim [pbc get]
	set dimx [lindex [lindex $dim 0] 0]
	set dimy [lindex [lindex $dim 0] 1]
	set dimz [lindex [lindex $dim 0] 2]
	
	# check pbc boundaries crossing
	set temp [$cnt get {x y z}]
	set max_x [lindex [lindex [lsort -real -index 0 $temp] end] 0]
	set max_y [lindex [lindex [lsort -real -index 1 $temp] end] 1]
	set max_z [lindex [lindex [lsort -real -index 2 $temp] end] 2]
	set min_x [lindex [lindex [lsort -real -index 0 $temp] 0] 0]
	set min_y [lindex [lindex [lsort -real -index 1 $temp] 0] 1]
	set min_z [lindex [lindex [lsort -real -index 2 $temp] 0] 2]
	if { ($max_x > $dimx) || ($max_y > $dimy) || ($max_z > $dimz) || ($min_x < 0) || ($min_y < 0) || ($min_z < 0) } {
		continue
	} 
	unset temp

	# calculate MI of BNT
	set I [draw principalaxes $cnt]
    puts "$fr $I"
	set A [orient $cnt [lindex $I 2] {0 0 1}]

	# rotate system
	$all move $A

	# projected center of mass of cnt
	set cm [center_of_mass $cnt]
	set cx [lindex $cm 0]
	set cy [lindex $cm 1]
	set pcm [list $cx $cy]

	# get pbc
	set dim [pbc get]
	set dimx [expr [lindex [lindex $dim 0] 0]/2.0]
	set dimy [expr [lindex [lindex $dim 0] 1]/2.0]
	set pdim [list $dimx $dimy]

	# range of cnt
	set mm [measure minmax $cnt]

	# radius of cnt
	set dist {}
	foreach i [$cnt get {x y z}] {
		set x [lindex $i 0]
		set y [lindex $i 1]

		# BNT radius
		set dist_sq [expr "($cx-$x)**2+($cy-$y)**2"]
		lappend dist [expr sqrt($dist_sq)]
		set avg_dist [getAvg $dist]
		unset dist
	}
	set zmax [lindex [lindex $mm 1] 2]
	set zmin [lindex [lindex $mm 0] 2]
	set h [expr "$zmax-$zmin"]

	# select atoms in cnt
	set inside [atomselect top "(x-$cx)**2+(y-$cy)**2 < $avg_dist**2 and z > $zmin and z < $zmax and type OW"]
	#puts "$fr $cx $cy $h [expr sqrt($avg_dist_sq)] [$inside num]"
	set str "[format %d $fr][format %8d [$inside num]][format %8.3f $avg_dist][format %8.1f $cx][format %8.1f $cy][format %8.1f $zmin][format %8.1f $zmax]"
	#puts $outfile $str
	#flush $outfile
}
puts " Done."
flush stdout
#close $outfile
