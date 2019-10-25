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
### Assumes two resnames: GRA and GRW
### 
set mos2 [atomselect top "resname MOS"]
set mos2_1 [atomselect top "resid 2 and type S_3a"]
set mos2_2 [atomselect top "resid 1 and type S_3b"]

set all [atomselect top "all"]
set n_frames [molinfo top get numframes]
set vdw [expr "3.267512965/2"]
set margin 5.0
set margin2 10.0

# output file
set filename [lindex [molinfo top get filename] 0 1]
append filename ".density.v2.profile"
puts "Info) Result will be saved at: $filename"
set outfile [open $filename w]

# progress bar
puts "Info) *** Sweep 1: Wrapping up the coordinates ***"
#pbc wrap -all -compound resid -sel "type OW or type HW"
pbc wrap -all -compound resid -center com

# progress bar
puts "Info) *** Sweep 2: Count the number ***"
progress_init $n_frames

puts $outfile "#fr\tn_wat\tdist\t|\tarea\tdensity\t|\tmass_margin\tarea_margin\tvolume_margin\tdensity_margin"


### count the number of confined water molecules
for { set fr 1 } { $fr <= $n_frames } { incr fr } {
	# load molecules
	molinfo top set frame $fr
	$mos2 frame $fr
	$all frame $fr
	$mos2 update
	$all update
	
	progress_tick $fr

	### range of mos2 capillaries
    # x: from cell param
    set x [lindex [lindex [pbc get] 0] 0]
    # y: from mos2
	set mm [measure minmax $mos2]
	set ymax [lindex [lindex $mm 1] 1]
	set ymin [lindex [lindex $mm 0] 1]
    set y [expr "$ymax-$ymin"]
    # z: thickness
    set z1avg [lavg [$mos2_1 get {z}]]
    set z2avg [lavg [$mos2_2 get {z}]]
    set dist [expr "$z1avg - $z2avg"]
    set eff_dist [expr "$dist - 2*$vdw"]

    set area [expr "$x * ($ymax - $ymin - (2*$vdw))"]
	set inside [atomselect top "(y > ($ymin+$vdw) and y < ($ymax-$vdw) and z > $z2avg and z < $z1avg) and (type OW or type HW)"]
	set inside_ow [atomselect top "same residue as ((y > ($ymin+$vdw) and y < ($ymax-$vdw) and z > $z2avg and z < $z1avg) and type OW)"]
    $inside frame $fr
    $inside update
    set mass 0
    foreach m [$inside get mass] {
       set mass [expr $mass + $m]
    }
    set density [expr "$mass/6.022/($area*$eff_dist)*10"]
    set area_density [expr "[$inside_ow num]/$area/3.0"]

    # margin1
	set inside_margin [atomselect top "(y > ($ymin+$margin) and y < ($ymax-$margin) and z > $z2avg and z < $z1avg) and (type OW or type HW)"]
    $inside_margin frame $fr
    $inside_margin update
    set mass_margin 0
    foreach m [$inside_margin get mass] {
       set mass_margin [expr $mass_margin + $m]
    }

    set area_margin [expr "$x * (($ymax-$margin)-($ymin+$margin))"]
    set density_margin [expr "$mass_margin/6.022/($area_margin*$eff_dist)*10"]

    # margin2
	set inside_margin2 [atomselect top "(y > ($ymin+$margin2) and y < ($ymax-$margin2) and z > $z2avg and z < $z1avg) and (type OW or type HW)"]
    $inside_margin2 frame $fr
    $inside_margin2 update
    set mass_margin2 0
    foreach m [$inside_margin2 get mass] {
       set mass_margin2 [expr $mass_margin2 + $m]
    }

    set area_margin2 [expr "$x * (($ymax-$margin2)-($ymin+$margin2))"]
    set density_margin2 [expr "$mass_margin2/6.022/($area_margin2*$eff_dist)*10"]

	# select atoms in mos2
	set str "[format %5d $fr][format %8d [$inside_ow num]][format %8.3f $dist][format %8.3f $eff_dist]\t|\t[format %12.3f $mass][format %12.3f $area][format %12.6f $density][format %12.6f $area_density]\t|\t[format %12.3f $mass_margin][format %12.3f $area_margin][format %12.6f $density_margin]\t|\t[format %12.3f $mass_margin2][format %12.3f $area_margin2][format %12.6f $density_margin2]"
	puts $outfile $str
	flush $outfile
}
puts " Done."
flush stdout
close $outfile
