# this script aligns principal axes of CNT to z axes
# and counts water molecules inside CNT
# version 20150204

### global
source "/qcfs/noische/scripts/vmd/mathtools.tcl"
source "/qcfs/noische/scripts/vmd/vmdtools.tcl"
source "/qcfs/noische/scripts/pbcdist.tcl"
source "/qcfs/noische/scripts/progress_bar.tcl"
package require Orient
namespace import Orient::orient


### initialization
# description
puts "Info) BNNT_getWaterShot.tcl: count the number of waters confined in BNNT."

# frame
set cnt [atomselect top "(resname BNT or resname CNT) and numbonds==2"]
set all [atomselect top "all"]
set n_frames [molinfo top get numframes]

# progress bar
puts "Info) *** Sweep 1: Wrapping up the coordinates ***"
progress_init $n_frames


### rotate the whole system by principal axes of BNNT
for { set fr 1 } { $fr <= $n_frames } { incr fr } {
	# load molecules
	molinfo top set frame $fr
	$cnt frame $fr
	$all frame $fr
	$cnt update
	$all update
	
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


### count the number of confined water molecules
for { set fr 1 } { $fr <= $n_frames } { incr fr } {
	# load molecules
	molinfo top set frame $fr
	$cnt frame $fr
	$all frame $fr
	
	# calculate MI of BNT
	set I [draw principalaxes $cnt]
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
		set xy [list $x $y]

		# BNT radius
		set dist_sq [expr "($cx-$x)**2+($cy-$y)**2"]
		lappend dist [expr sqrt($dist_sq)]
		#lappend dist [pbcdist $pcm $xy $pdim]
		set avg_dist [getAvg $dist]
		unset dist
	}
	set zmax [lindex [lindex $mm 1] 2]
	set zmin [lindex [lindex $mm 0] 2]
	set h [expr "$zmax-$zmin"]

	# select atoms in cnt
	set inside [atomselect top "(x-$cx)**2+(y-$cy)**2 < $avg_dist**2 and z > $zmin and z < $zmax and type OW"]
	puts "[$inside num]"
	puts "(x-$cx)**2+(y-$cy)**2 < $avg_dist**2 and z > $zmin and z < $zmax and type OW"
}
