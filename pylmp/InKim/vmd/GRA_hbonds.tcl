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
puts "Info) GRA_hbonds.tcl: count the number of hbonds of water confined in two Graphene Sheet."

# frame
### Assumes two resnames: GRA and GRW
### 
set graphene [atomselect top "resname GRA"]
set graphene1 [atomselect top "resname GRA and resid 1"]
set graphene2 [atomselect top "resname GRA and resid 2"]
set wall [atomselect top "resname GRW"]
set all [atomselect top "all"]
set n_frames [molinfo top get numframes]

# output file
set filename [lindex [molinfo top get filename] 0 1]
append filename ".GRA_hbonds.vmd.v2.profile"
puts "Info) Result will be saved at: $filename"
set outfile [open $filename w]

# progress bar
puts "Info) *** Sweep 1: Wrapping up the coordinates ***"
pbc wrap -all -compound resid

# progress bar
puts "Info) *** Sweep 2: Count the number ***"
progress_init $n_frames

puts $outfile "#fr\tn_wat\tn_hbond\thbond/wat"


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

    set margin 10.0
    #set area_margin [expr "(($ymax-$margin)-($ymin+$margin))*(($xmax-$margin)-($xmin+$margin))"]
	set inside_margin [atomselect top "(x > $xmin+$margin and x < $xmax-$margin and y > $ymin+$margin and y < $ymax-$margin and z > $zmin and z < $zmax) and type OW"]
    set n_wat_margin [$inside_margin num]
    set hbonds [measure hbonds 3.0 20 $inside_margin]
    set n_hbonds [llength [lindex $hbonds 0]]
    set avrg_margin [expr "$n_hbonds.0/$n_wat_margin.0"]
    #set avrg_margin 0.0

	# select atoms in graphene
	set str "[format %d $fr]\t[format %8d [$inside_margin num]]\t[format %8d $n_hbonds]\t[format %8.3f $avrg_margin]"
	puts $outfile $str
	flush $outfile
}
puts " Done."
flush stdout
close $outfile
