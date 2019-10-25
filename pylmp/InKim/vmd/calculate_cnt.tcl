# this script will calculate height and radius of CNT
# version 20151208

source "/qcfs/noische/scripts/vmd/mathtools.tcl"
source "/qcfs/noische/scripts/vmd/vmdtools.tcl"

set cnt [atomselect top "resname BNT or resname CNT or type C_2G"]

# projected center of mass of cnt
set cm [center_of_mass $cnt]
set cx [lindex $cm 0]
set cy [lindex $cm 1]

# range of cnt
set mm [measure minmax $cnt]

foreach i [$cnt get {x y z}] {
    set x [lindex $i 0]
    set y [lindex $i 1]
    set xy [list $x $y]

    # BNT radius
    set dist_sq [expr "($cx-$x)**2+($cy-$y)**2"]
    lappend dist [expr sqrt($dist_sq)]
    #lappend dist [pbcdist $pcm $xy $pdim]
    set avg_dist [getAvg $dist]
    unset dist
}
set zmax [lindex [lindex $mm 1] 2]
set zmin [lindex [lindex $mm 0] 2]
set h [expr "$zmax-$zmin"]

puts "CNT radius: $avg_dist A"
puts "CNT height: $h A"
puts "CNT selection: (x-$cx)**2+(y-$cy)**2 < $avg_dist**2 and z > $zmin and z < $zmax"
