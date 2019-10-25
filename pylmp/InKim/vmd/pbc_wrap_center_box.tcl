source "/qcfs/noische/scripts/progress_bar.tcl"
source "/qcfs/noische/scripts/vmd/vmdtools.tcl"
source "/qcfs/noische/scripts/vmd/mathtools.tcl"

proc wrap_center {selection} {
    set all [atomselect top "all"]
    set sel [atomselect top "$selection"]
    set n_frames [molinfo top get numframes]
    progress_init $n_frames

    for { set fr 1 } { $fr <= $n_frames } { incr fr } {
        $all frame $fr
        $all update
        $sel frame $fr
        $sel update

        progress_tick $fr
        # center of mass
        set cm [measure center $sel]

        # pbc
        set dim [pbc get]
        set x [expr [lindex [lindex $dim 0] 0]/2.0]
        set y [expr [lindex [lindex $dim 0] 1]/2.0]
        set z [expr [lindex [lindex $dim 0] 2]/2.0]

        moveback $all $cm
        $all moveby "$x $y $z"

    }

    pbc wrap -all -compound resid -center com -centersel "$selection"
}
