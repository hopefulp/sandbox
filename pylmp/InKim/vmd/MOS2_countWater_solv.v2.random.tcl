# this script counts water molecules between two mos2 sheets
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
puts "Info) GRA_countWater.tcl: count the number of waters confined in two MoS2 Sheet."

# frame
set mos2 [atomselect top "resname MOS"]
set mos2_1 [atomselect top "resid 2 and type S_3a"]
set mos2_2 [atomselect top "resid 1 and type S_3b"]

set all [atomselect top "all"]
set n_frames [molinfo top get numframes]
set vdw [expr "3.270615945/2"]
set margin 10.0
set avg_dist {}
set avg_eff_dist {}
set avg_density {}
set avg_eff_density {}

# output file
set filename [lindex [molinfo top get filename] 0 1]
append filename ".density.solv.random.profile"
puts "Info) Result will be saved at: $filename"
set outfile [open $filename w]

# log file
append filename ".log"
puts "Info) Calculation log will be saved at: $filename"
set logfile [open $filename w]

# progress bar
puts "Info) *** Sweep 1: Wrapping up the coordinates ***"
#pbc wrap -all -compound resid -center com -sel "type OW or type HW"

# progress bar
puts "Info) *** Sweep 2: Count the number ***"
progress_init $n_frames

puts $outfile "#fr\tn_wat\tarea\tdist\tdensity\t|\tmass_margin\tarea_margin\tvolume_margin\tdensity_margin"
puts $logfile "#fr\ttrial\tx1\tx2\ty1\ty2\tdx\tdy\tmass\tarea\tdens\tnum_ow\tarea_dens"

### count the number of confined water molecules
for { set fr 1 } { $fr <= $n_frames } { incr fr } {
	# load molecules
	molinfo top set frame $fr
	$mos2 frame $fr
	$all frame $fr
	$mos2 update
	$all update
	
	progress_tick $fr

	# range of mos2 capillaries
    set mm [measure minmax $mos2]

	set xmax [lindex [lindex $mm 1] 0]
	set xmin [lindex [lindex $mm 0] 0]
	set ymax [lindex [lindex $mm 1] 1]
	set ymin [lindex [lindex $mm 0] 1]
	set zmax [lindex [lindex $mm 1] 2]
	set zmin [lindex [lindex $mm 0] 2]

    #set dim [pbc get]
    set x [expr "$xmax - $xmin"]
    set y [expr "$ymax - $ymin"]

    set z1avg [lavg [$mos2_2 get {z}]]
    set z2avg [lavg [$mos2_1 get {z}]]
    set dist [expr "$z1avg - $z2avg"]
    set eff_dist [expr "$dist - 2 * $vdw"]

    # create a random region
    set fr_density {}
    set fr_area_density {}
    set trial 0
    while { $trial <= 50 } {
        set x1 [expr "rand() * $x"]
        set x2 [expr "rand() * $x"]
        if { $x1 > $x2 } {set x1 $x2[set x2 $x1; list]}
        if { $x2 - $x1 < 30.0 } { continue }
        set dx [expr "$x2 - $x1"]

        set y1 [expr "rand() * $y"]
        set y2 [expr "rand() * $y"]
        if { $y1 > $y2 } {set y1 $y2[set y2 $y1; list]}
        if { $y2 - $y1 < 30.0 } { continue }
        set dy [expr "$y2 - $y1"]

        set area [expr "$dx * $dy"]
	    set inside [atomselect top "(x > ($x1 + $xmin) and x < ($x2 + $xmin) and y > ($y1 + $ymin) and y < ($y2 + $ymin) and z > $z2avg and z < $z1avg) and (type OW or type HW)"]
	    set inside_ow [atomselect top "(x > ($x1 + $xmin) and x < ($x2 + $xmin) and y > ($y1 + $ymin) and y < ($y2 + $ymin) and z > $z2avg and z < $z1avg) and (type OW)"]
        $inside frame $fr
        $inside update
        set num_ow [$inside_ow num]

        set mass 0
        foreach m [$inside get mass] {
           set mass [expr $mass + $m]
        }

        set density [expr "$mass / 6.022 / ($area * $eff_dist) * 10"]
        set area_density [expr "$num_ow / $area"]
        lappend fr_density $density
        lappend fr_area_density $area_density
        puts $logfile "$fr $trial [format %8.3f $x1] [format %8.3f $x2] [format %8.3f $y1] [format %8.3f $y2] [format %8.3f $dx] [format %8.3f $dy] [format %12.6f $mass] [format %12.6f $area] [format %12.6f $density] [format %6d $num_ow] [format %12.6f $area_density]"
        
        unset x1 x2 dx y1 y2 dy area num_ow mass density area_density
        set inside ""
        set inside_ow ""
        unset inside
        unset inside_ow
        incr trial
    }

    set avg_density [lavg $fr_density]
    set avg_area_density [lavg $fr_area_density]

    # select atoms in mos2
	set str "[format %5d $fr] [format %8.3f $dist] [format %8.3f $eff_dist] [format %12.6f $avg_density] [format %12.6f $avg_area_density]"
	puts $outfile $str
	flush $outfile
    flush $logfile

    set fr_density ""
    set fr_area_density ""
    unset avg_density avg_area_density fr_density fr_area_density
}
puts " Done."
flush stdout
close $outfile
close $logfile
