### global
source "/qcfs/noische/scripts/progress_bar.tcl"
set n_frames [molinfo top get numframes]
proc lavg L {expr ([join $L +])/[llength $L].}

### initialization
set result {}
#set filename [lindex [molinfo top get filename] 0 1]
#append filename ".sasa.1p4.profile"

progress_init $n_frames

### main loop
for { set fr 0 } { $fr <= $n_frames } { incr fr } {
	# load molecules
	molinfo top set frame $fr
	set sel [atomselect top "not water"]
    progress_tick $fr
	
    lappend result [measure sasa 1.4 $sel]
}

### output
puts "[format %8.3f [lavg $result]]"

### END OF FUNCTION ###
