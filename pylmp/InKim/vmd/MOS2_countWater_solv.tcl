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
set graphene1 [atomselect top "resid 1 and type S_3b"]
# bottom
set graphene2 [atomselect top "resid 2 and type S_3a"]
set all [atomselect top "all"]
set n_frames [molinfo top get numframes]
set vdw 1.685
set margin 10.0
set avg_dist {}
set avg_eff_dist {}
set avg_density {}
set avg_eff_density {}

# output file
set filename [lindex [molinfo top get filename] 0 1]
append filename ".MOS2_countWater_solv.effective.vmd.margin10.v2.profile"
puts "Info) Result will be saved at: $filename"
set outfile [open $filename w]

# progress bar
puts "Info) *** Sweep 1: Wrapping up the coordinates ***"
pbc wrap -all -compound resid

# progress bar
puts "Info) *** Sweep 2: Count the number ***"
progress_init $n_frames

puts $outfile "#fr\tn_wat\tarea\tdist\tdensity\t|\tmass_margin\tarea_margin\tvolume_margin\tdensity_margin"


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
    set mm [measure minmax $graphene]

	set xmax [lindex [lindex $mm 1] 0]
	set xmin [lindex [lindex $mm 0] 0]
	set ymax [lindex [lindex $mm 1] 1]
	set ymin [lindex [lindex $mm 0] 1]
	set zmax [lindex [lindex $mm 1] 2]
	set zmin [lindex [lindex $mm 0] 2]

    set dim [pbc get]
    set x [lindex [lindex $dim 0] 0]
    set y [lindex [lindex $dim 0] 1]

    set z1avg [lavg [$graphene1 get {z}]]
    set z2avg [lavg [$graphene2 get {z}]]

    #set area [expr "$x*$y"]
    set area [expr "($ymax-$ymin+2*$vdw)*($xmax-$xmin+2*$vdw)"]
    set dist [expr "$z1avg-$z2avg"]
	set inside [atomselect top "(x > ($xmin-$vdw) and x < ($xmax+$vdw) and y > ($ymin-$vdw) and y < ($ymax+$vdw) and z > $z2avg and z < $z1avg)"]
    $inside frame $fr
    $inside update
    set n_wat [$inside num]
    set mass 0
    foreach coord [$inside get {x y z}] m [$inside get mass] {
       set mass [expr $mass + $m]
    }
    set density [expr "$mass/6.022/($area*($dist-2*$vdw))*10"]

	set inside_margin [atomselect top "x > $xmin+$margin and x < $xmax-$margin and y > $ymin+$margin and y < $ymax-$margin and z > $z2avg and z < $z1avg"]
    $inside_margin frame $fr
    $inside_margin update
    set mass_margin 0
    foreach coord [$inside_margin get {x y z}] m [$inside_margin get mass] {
       # sum of the masses
       set mass_margin [expr $mass_margin + $m]
    }

    set area_margin [expr "(($ymax-$margin)-($ymin+$margin))*(($xmax-$margin)-($xmin+$margin))"]
    set volume_margin [expr "$area_margin*($dist - 2*$vdw)"]
    set density_margin [expr "$mass_margin/6.022/$volume_margin*10"]

	# select atoms in graphene
	set str "[format %d $fr]\t[format %8.3f $mass]\t[format %8.3f $dist]\t|\t[format %8.3f $area]\t[format %12.9f $density]\t|\t[format %8.3f $mass_margin]\t[format %12.3f $area_margin]\t[format %12.3f $volume_margin]\t[format %12.6f $density_margin]"
	puts $outfile $str
	flush $outfile
}
puts " Done."
flush stdout
close $outfile
