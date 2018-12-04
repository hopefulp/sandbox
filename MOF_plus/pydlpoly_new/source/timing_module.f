      module timing_module

c***********************************************************************
c
c      module to time execution of dl_poly
c      this is part of the pydlpoly project
c      note that some features will work only in connection with
c      pydlpoly
c
c      (C) R. Schmid 2011 , RUB
c
c***********************************************************************

      use utility_module

      implicit none

      logical   ltime
      
      integer, parameter ::  ntimers=10

      integer tcounts(ntimers)
      real(8) timers(ntimers)
      integer current_timer
      logical timeron
      real(8) last_time, accum_time
      real(8) time_init
      data timeron,ltime/.false.,.false./
      
      save ltime, tcounts, timers, current_timer, timeron, last_time
      save accum_time, time_init 

      contains

      subroutine timer_init()
c        switch timers on or reset

      implicit none

      integer i

      do i=1,ntimers
        tcounts(i) = 1
        timers(i)  = 0.0d0
        timeron= .false.
        ltime  = .true.
        current_timer = 0
      end do
       
      call timchk(0, time_init)

      end subroutine timer_init

      subroutine timer_on(t)
c        start timer t, if we are already timing switch context 

      implicit none
      integer t
      real(8) current

      call timchk(0, current)

      if (timeron) then
        timers(current_timer) = timers(current_timer) + accum_time 
     x                              + (current-last_time)
        tcounts(current_timer) = tcounts(current_timer)+1
      end if
      
      current_timer = t
      timeron= .true.
      accum_time = 0.0d0
      last_time = current

      end subroutine timer_on

      subroutine timer_stop()
c         stop current timer - it stays the current context!

      implicit none

      real(8) current

      call timchk(0, current)
      if (timeron) then
         accum_time = accum_time + (current-last_time)
         timeron= .false.
      end if

      end subroutine timer_stop

      subroutine timer_start()
c        start again (if running od nothing)

      implicit none

      real(8) current

      call timchk(0, current)
      if (.not.timeron) then
        timeron= .true.
        last_time = current
      end if
      
      end subroutine timer_start

      subroutine timer_off()
c        switch current timer off.

      implicit none

      real(8) current

      call timchk(0, current)
      if (timeron) then 
        accum_time = accum_time + (current-last_time)
        timeron=.false.
      end if
      
      timers(current_timer) = timers(current_timer) + accum_time
      tcounts(current_timer) = tcounts(current_timer) +1
      current_timer = 0

      return
      end subroutine timer_off 

      subroutine timer_since_init(tsi)

      real(8) current, tsi
      call timchk(0,current)
      tsi = current-time_init

      end subroutine timer_since_init


      end module timing_module