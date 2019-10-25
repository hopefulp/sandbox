### CoM calculation
proc center_of_mass {selection} { 
	if {[$selection num] <= 0} {
		error "center_of_mass: needs a selection with atoms"
	}

	set com [veczero]
	set mass 0

	foreach coord [$selection get {x y z}] m [$selection get mass] {
		set mass [expr $mass + $m]
		set com [vecadd $com [vecscale $m $coord]]
	}

	if {$mass == 0} {
		error "center_of_mass: total mass is zero"
	}
	
	return [vecscale [expr 1.0/$mass] $com]
}

proc moveback {sel offset} {
	foreach coord [$sel get {x y z}] {
		# moveby uses vecadd
		lappend newcoords [vecsub $coord $offset]
	}
	$sel set {x y z} $newcoords 
}

proc draw_cube {minx miny minz maxx maxy maxz} {
    draw materials off
    draw color yellow

    draw line "$minx $miny $minz" "$maxx $miny $minz"
    draw line "$minx $miny $minz" "$minx $maxy $minz"
    draw line "$minx $miny $minz" "$minx $miny $maxz"

    draw line "$maxx $miny $minz" "$maxx $maxy $minz"
    draw line "$maxx $miny $minz" "$maxx $miny $maxz"

    draw line "$minx $maxy $minz" "$maxx $maxy $minz"
    draw line "$minx $maxy $minz" "$minx $maxy $maxz"

    draw line "$minx $miny $maxz" "$maxx $miny $maxz"
    draw line "$minx $miny $maxz" "$minx $maxy $maxz"

    draw line "$maxx $maxy $maxz" "$maxx $maxy $minz"
    draw line "$maxx $maxy $maxz" "$minx $maxy $maxz"
    draw line "$maxx $maxy $maxz" "$maxx $miny $maxz"
}
