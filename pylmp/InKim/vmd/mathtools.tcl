### average function
proc getAvg { rList } {
	set n [llength $rList]
	if {!$n} {return 0.0}

	set avg 0.0
	foreach r $rList {
		set avg [expr $avg + $r]
	}
	set avg [expr $avg/double($n)]

	return $avg
}

