# this script counts water molecules between two graphene sheets
# version 20160115
# from CNT_countWaterComplex_fromBNNT.tcl

### global
source "/qcfs/noische/scripts/vmd/mathtools.tcl"
source "/qcfs/noische/scripts/vmd/vmdtools.tcl"
source "/qcfs/noische/scripts/progress_bar.tcl"

proc lavg L {expr ([join $L +])/[llength $L].}

### initialization
# description
puts "Info) GRA_countWater.tcl: count the number of waters confined in two Graphene Sheet."

# frame
### Assumes two resnames: GRA and GRW
### 
set graphene1 [atomselect top "resname MOS and resid 1 and type S_3a"]
set graphene2 [atomselect top "resname MOS and resid 2 and type S_3b"]
set n_frames [molinfo top get numframes]

puts "Info) *** Sweep 1: Wrapping up the coordinates ***"
#pbc wrap -all -compound resid

puts "Info) *** Sweep 2: Measure Interlayer Distances ***"
progress_init $n_frames

set z1 {}
set z2 {}

### count the number of confined water molecules
for { set fr 1 } { $fr <= $n_frames } { incr fr } {
	# load molecules
	molinfo top set frame $fr
	
	progress_tick $fr

	# range of graphene
    set z1avg [lavg [$graphene1 get {z}]]
    set z2avg [lavg [$graphene2 get {z}]]

    lappend z1 $z1avg
    lappend z2 $z2avg
}

puts [expr "[lavg $z1]-[lavg $z2]"]

puts " Done."
flush stdout
