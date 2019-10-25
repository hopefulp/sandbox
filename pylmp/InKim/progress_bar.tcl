#!/usr/bin/tclsh
proc progress_init {tot} {
   set ::progress_start     [clock seconds]
   set ::progress_last      0
   set ::progress_last_time 0
   set ::progress_tot       $tot
}

# We update if there's a 5% difference or a 5 second difference

proc progress_tick {cur} {
   set now [clock seconds]
   set tot $::progress_tot

   if {$cur > $tot} {
       set cur $tot
   }
   if {($cur >= $tot && $::progress_last < $cur) ||
       ($cur - $::progress_last) >= (0.05 * $tot) ||
       ($now - $::progress_last_time) >= 1} {
       set ::progress_last $cur
       set ::progress_last_time $now
       set percentage [expr round($cur*100/$tot)]
       set ticks [expr $percentage/2]
       if {$cur == 0} {
           set eta   ETA:[format %5s Unknown]
       } elseif {$cur >= $tot} {
           set eta   TOT:[format %5d [expr int($now - $::progress_start)]]s
       } else {
           set eta   ETA:[format %5d [expr int(($tot - $cur) * ($now - $::progress_start)/$cur)]]s
       }
       set lticks [expr 50 - $ticks]
       set str "[format %3d $percentage] % |[string repeat = $ticks]"
       append str "[string repeat . $lticks]| [format %d $cur] / [format %d $tot] | $eta\r"
       puts -nonewline stdout $str
       if {$cur >= $tot} {
           puts ""
       }
       flush stdout
   }
}

# BELOW IS THE SAMPLE
#progress_init 5000

#for {set i 0} {$i < 6200} {incr i 200} {
#   progress_tick $i
#   after 200
#}
