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
puts "Info) GRA_countWater.tcl: count the number of waters confined in two Graphene Sheet. (Kwac's structure)"

# frame
### Assumes two resnames: GRA for slit and GRW for reservoir walls
### 
set graphene [atomselect top "resname GRA"]
set graphene1 [atomselect top "resname GRA and resid 1"]
set graphene2 [atomselect top "resname GRA and resid 2"]
set wall [atomselect top "resname GRW"]
set all [atomselect top "all"]
set n_frames [molinfo top get numframes]
set margin 10.0
set vdw 3.38383824/2

# output file
set filename [lindex [molinfo top get filename] 0 1]
append filename ".GRA_countWater.kwac.profile"
puts "Info) Result will be saved at: $filename"
set outfile [open $filename w]

# progress bar
puts "Info) *** Sweep 1: Wrapping up the coordinates ***"
pbc wrap -all -compound resid

# progress bar
puts "Info) *** Sweep 2: Count the number ***"
progress_init $n_frames

puts $outfile "#fr    x   y   dist    eff_dist    volume  eff_volume  density eff_density area_density    m_volume    m_eff_volume    m_density   m_eff_density m_area_density"

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

    set z1avg [lavg [$graphene1 get {z}]]
    set z2avg [lavg [$graphene2 get {z}]]

    ### Total Density
    # volume
    set pbcx [expr "$xmax - $xmin"]
    set pbcy [expr "$ymax - $ymin"]
    set area [expr "$pbcx * $pbcy"]
    set dist [expr "$z2avg-$z1avg"]
    set eff_dist [expr "$dist - 2 * $vdw"]
    set volume [expr "$area * $dist"]
    set eff_volume [expr "$area * $eff_dist"]

    # mass
    set mass 0
    set mass_sel "x > $xmin and x < $xmax and y > $ymin and y < $ymax and z > $zmin and z < $zmax"
    set sel [atomselect top $mass_sel]
    foreach m [$sel get mass] {
        set mass [expr $mass + $m]
    }

    # density
    set density [expr "$mass/6.022/($area*$dist)*10"]
    set eff_density [expr "$mass/6.022/($area*$eff_dist)*10"]

    # area density
    set sel "x > $xmin and x < $xmax and y > $ymin and y < $ymax and z > $zmin and z < $zmax and type OW"
    set wat [atomselect top $sel]
    set n_wat [$wat num]
    if {$n_wat != 0.0} {
        set area_density [expr "$area / $n_wat"]
    } else {
        set area_density 0.0
    }

    ### Central Density
    # margin volume
    set margin_area [expr "(($ymax-$margin)-($ymin+$margin))*(($xmax-$margin)-($xmin+$margin))"]
    set margin_volume [expr "$margin_area * $dist"]
    set margin_eff_volume [expr "$margin_area * $eff_dist"]

    # margin mass
    set margin_mass 0
    set margin_mass_sel_str "x > $xmin+$margin and x < $xmax-$margin and y > $ymin+$margin and y < $ymax-$margin and z > $zmin and z < $zmax"
    set margin_mass_sel [atomselect top $margin_mass_sel_str]
    foreach coord [$margin_mass_sel get {x y z}] m [$margin_mass_sel get mass] {
        set margin_mass [expr $mass + $m]
    }

    # margin density
    set margin_density [expr "$margin_mass/6.022/($margin_area*$dist)*10"]
    set margin_eff_density [expr "$margin_mass/6.022/($margin_area*$eff_dist)*10"]

    # number of water
    set margin_sel "(x > $xmin+$margin and x < $xmax-$margin and y > $ymin+$margin and y < $ymax-$margin and z > $zmin and z < $zmax) and type OW"
    set margin_wat [atomselect top $margin_sel]
    set margin_n_wat [$margin_wat num]
    if {$margin_n_wat != 0.0} {
        set margin_area_density [expr "$margin_area / $n_wat"]
    } else {
        set margin_area_density 0.0
    }

	# write output
	set str "[format %d $fr][format %8.3f $pbcx][format %8.3f $pbcy][format %8.3f $dist][format %8.3f $eff_dist][format %12.3f $volume][format %12.3f $eff_volume][format %8.3f $density][format %8.3f $eff_density][format %8.3f $area_density][format %12.3f $margin_volume][format %12.3f $margin_eff_volume][format %8.3f $margin_density][format %8.3f $margin_eff_density][format %8.3f $margin_area_density]"
	puts $outfile $str
	flush $outfile
}
puts " Done."
flush stdout
close $outfile
