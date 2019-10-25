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
puts "Info) GRA_countWater.tcl: count the number of waters confined in two MoS2 Sheet."

# frame
### Assumes two resnames: GRA and GRW
### 
set graphene [atomselect top "resname MOS"]
# top
set graphene1 [atomselect top "resid 2 and type S_3a"]
# bottom
set graphene2 [atomselect top "resid 1 and type S_3b"]

set all [atomselect top "all"]
set n_frames [molinfo top get numframes]
set vdw 1.685
set vdw [expr "3.267512965/2"]
set margin 5.0
set margin2 10.0

# output file
set filename [lindex [molinfo top get filename] 0 1]
append filename ".MOS2_countWater_solv.effective.vmd.margin10.v2.profile"
puts "Info) Result will be saved at: $filename"
set outfile [open $filename w]

# progress bar
puts "Info) *** Sweep 1: Wrapping up the coordinates ***"
pbc wrap -all -compound resid -sel "type OW or type HW"

# progress bar
puts "Info) *** Sweep 2: Count the number ***"
progress_init $n_frames

puts $outfile "#fr\tn_wat\tdist\t|\tarea\tdensity\t|\tmass_margin\tarea_margin\tvolume_margin\tdensity_margin"


### count the number of confined water molecules
for { set fr 1 } { $fr <= $n_frames } { incr fr } {
	# load molecules
	molinfo top set frame $fr
	$graphene frame $fr
	$all frame $fr
	$graphene update
	$all update
	
	progress_tick $fr

	# range of mos2 capillaries
	set mm [measure minmax $graphene]

	set xmax [lindex [lindex $mm 1] 0]
	set xmin [lindex [lindex $mm 0] 0]
	set ymax [lindex [lindex $mm 1] 1]
	set ymin [lindex [lindex $mm 0] 1]
	set zmax [lindex [lindex $mm 1] 2]
	set zmin [lindex [lindex $mm 0] 2]

    set z1avg [lavg [$graphene1 get {z}]]
    set z2avg [lavg [$graphene2 get {z}]]

    set area [expr "($ymax-$ymin+2*$vdw)*($xmax-$xmin+2*$vdw)"]
    set dist [expr "$z1avg-$z2avg"]
    set eff_dist [expr "$dist-2*$vdw"]
	#set inside [atomselect top "(x > ($xmin-$vdw) and x < ($xmax+$vdw) and y > ($ymin-$vdw) and y < ($ymax+$vdw) and z > $z2avg and z < $z1avg) and type OW"]
	set inside [atomselect top "(x > ($xmin-$vdw) and x < ($xmax+$vdw) and y > ($ymin-$vdw) and y < ($ymax+$vdw) and z > $z2avg and z < $z1avg)"]
	set inside_ow [atomselect top "(x > ($xmin-$vdw) and x < ($xmax+$vdw) and y > ($ymin-$vdw) and y < ($ymax+$vdw) and z > $z2avg and z < $z1avg) and type OW"]
    $inside frame $fr
    $inside update
    set mass 0
    foreach coord [$inside get {x y z}] m [$inside get mass] {
       # sum of the masses
       set mass [expr $mass + $m]
    }
    set density [expr "$mass/6.022/($area*$eff_dist)*10"]
    set area_density [expr "[$inside_ow num]/$area"]

    # margin1
	set inside_margin [atomselect top "x>$xmin+$margin and x<$xmax-$margin and y>$ymin+$margin and y<$ymax-$margin and z>$z2avg and z<$z1avg"]
    $inside_margin frame $fr
    $inside_margin update
    set mass_margin 0
    foreach coord [$inside_margin get {x y z}] m [$inside_margin get mass] {
       # sum of the masses
       set mass_margin [expr $mass_margin + $m]
    }

    set area_margin [expr "(($ymax-$margin)-($ymin+$margin))*(($xmax-$margin)-($xmin+$margin))"]
    set volume_margin [expr "$area_margin*($dist-2*$vdw)"]
    set density_margin [expr "$mass_margin/6.022/$volume_margin*10"]

    # margin2
	set inside_margin2 [atomselect top "x > $xmin+$margin2 and x < $xmax-$margin2 and y > $ymin+$margin2 and y < $ymax-$margin2 and z > $z2avg and z < $z1avg"]
    $inside_margin2 frame $fr
    $inside_margin2 update
    set mass_margin2 0
    foreach coord [$inside_margin2 get {x y z}] m [$inside_margin2 get mass] {
       # sum of the masses
       set mass_margin2 [expr $mass_margin2 + $m]
    }

    set area_margin2 [expr "(($ymax-$margin2)-($ymin+$margin2))*(($xmax-$margin2)-($xmin+$margin2))"]
    set volume_margin2 [expr "$area_margin2*($dist - 2*$vdw)"]
    set density_margin2 [expr "$mass_margin2/6.022/$volume_margin2*10"]

	# select atoms in graphene
	set str "[format %5d $fr][format %8d [$inside_ow num]][format %8.3f $dist][format %8.3f $eff_dist]\t|\t[format %12.3f $mass][format %12.3f $area][format %12.6f $density]\t|\t[format %12.6f $area_density]\t|\t[format %12.3f $mass_margin][format %12.3f $area_margin][format %12.6f $density_margin]\t|\t[format %12.3f $mass_margin2][format %12.3f $area_margin2][format %12.3f $volume_margin2][format %12.6f $density_margin2]"
	puts $outfile $str
	flush $outfile
}
puts " Done."
flush stdout
close $outfile
