# this script counts water molecules between two graphene sheets
# version 20160115
# from CNT_countWaterComplex_fromBNNT.tcl

### global
source "/qcfs/noische/scripts/vmd/mathtools.tcl"
source "/qcfs/noische/scripts/vmd/vmdtools.tcl"
source "/qcfs/noische/scripts/pbcdist.tcl"
source "/qcfs/noische/scripts/progress_bar.tcl"

proc lavg L {expr ([join $L +])/[llength $L].}

### initialization
# description
puts "Info) GRA_countWater.tcl: count the number of waters confined in two Graphene Sheet."

# frame
### Assumes two resnames: GRA for slit and GRW for reservoir walls
### 
set graphene [atomselect top "resname GRA"]
set graphene1 [atomselect top "resname GRA and resid 1"]
set graphene2 [atomselect top "resname GRA and resid 2"]
set inside_water [atomselect top "type OW and chain I"]
set wall [atomselect top "resname GRW"]
set all [atomselect top "all"]
set n_frames [molinfo top get numframes]
set vdw 1.692
set margin 10.0
set avg_dist {}
set avg_eff_dist {}
set avg_density {}
set avg_eff_density {}

# output file
set filename [lindex [molinfo top get filename] 0 1]
append filename ".GRA_density.vmd.profile"
puts "Info) Result will be saved at: $filename"
set outfile [open $filename w]

# progress bar
puts "Info) *** Sweep 1: Wrapping up the coordinates ***"
pbc wrap -all -compound resid

# progress bar
puts "Info) *** Sweep 2: Count the number ***"
progress_init $n_frames


### count the number of confined water molecules
for { set fr 1 } { $fr <= $n_frames } { incr fr } {
	# load molecules
	molinfo top set frame $fr
	$graphene frame $fr
	$all frame $fr
	$graphene update
	$all update
	
	progress_tick $fr

	# range of graphene
	#set mm [measure minmax $graphene]
    set dim [pbc get]
    set x [lindex [lindex $dim 0] 0]
    set y [lindex [lindex $dim 0] 1]

	#set xmax [lindex [lindex $mm 1] 0]
	#set xmin [lindex [lindex $mm 0] 0]
	#set ymax [lindex [lindex $mm 1] 1]
	#set ymin [lindex [lindex $mm 0] 1]
	#set zmax [lindex [lindex $mm 1] 2]
	#set zmin [lindex [lindex $mm 0] 2]

    set z1avg [lavg [$graphene1 get {z}]]
    set z2avg [lavg [$graphene2 get {z}]]

    set n_wat [$inside_water num]
    #set x [expr "$xmax-$xmin"]
    #set y [expr "$ymax-$ymin"]
    #set area [expr "($x+2*$vdw)*($y+2*$vdw)"]
    set area [expr "$x*$y"]
    set dist [expr "$z2avg-$z1avg"]
    set eff_dist [expr "$dist-2*$vdw"]
    set density [expr "(18.0154*$n_wat)/6.022/($area*$dist)*10"]
    set eff_density [expr "(18.0154*$n_wat)/6.022/($area*$eff_dist)*10"]

	# select atoms in graphene
	set str "[format %d $fr]\t[format %8.3f $x]\t[format %8.3f $y]\t[format %8.3f $area]\t[format %8.3f $dist]\t[format %8.3f $eff_dist]\t[format %12.9f $density]\t[format %12.9f $eff_density]"
	puts $outfile $str
	flush $outfile
    lappend avg_dist $dist
    lappend avg_eff_dist $eff_dist
    lappend avg_density $density
    lappend avg_eff_density $eff_density
}
puts "* average dist: [format %8.3f [lavg $avg_dist]]"
puts "* average eff_dist: [format %8.3f [lavg $avg_eff_dist]]"
puts "* average density: [format %8.3f [lavg $avg_density]]"
puts "* average eff_density: [format %8.3f [lavg $avg_eff_density]]"
puts " Done."
flush stdout
close $outfile
