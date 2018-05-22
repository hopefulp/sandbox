cRS
c     revised version using a module 
c     main reason: use a specified communicator instead of the COMM_WORLD
c       
cRS
      module comm_module

      implicit none

c      this includes mpif.h!!!
      include "comms.inc"

      integer pydlpoly_comm

      integer error_code

      save pydlpoly_comm, error_code

      end module comm_module
    

      subroutine initcomms()
      
c*********************************************************************
c     
c     communication harness initialisation
c     
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c     
c     wl
c     2008/01/14 13:33:07
c     1.4
c     Exp
c     
c*********************************************************************
      use comm_module
      
      implicit none
      
c      include "comms.inc"
      
      integer ierr

CMPIU      define MPI_init MPI_init_

      call MPI_init(ierr)

      return
      end

      subroutine machine(idnode,mxnode)

c*********************************************************************
c     
c     dl_poly subroutine for obtaining charcteristics of
c     the computer on which the program is being run
c     
c     copyright daresbury laboratory 1992
c     author - w.smith july 1992
c     
c     MPI version - t.forester may 1995
c
c     wl
c     1.4
c     Exp
c*********************************************************************

      implicit none

      integer idnode,mxnode,mynode,numnodes

      idnode=mynode()
      mxnode=numnodes()

      return
      end

      integer function mynode()

c*********************************************************************
c
c     routine to determine identity of processing node 
c
c     MPI version - t.forester may 1995
c
c     wl
c     2008/01/14 13:33:07
c     1.4
c     Exp
c
c*********************************************************************
      use comm_module 

      implicit none

c      include "comms.inc"

      integer ierr

CMPIU define MPI_COMM_RANK MPI_COMM_RANK_

      call MPI_COMM_RANK(pydlpoly_comm, mynode ,ierr)

      return
      end

      integer function nodedim()

c*********************************************************************
c
c     calculate dimension of hypercube
c
c     MPI version - t.forester may 1995
c
c     wl
c     2008/01/14 13:33:07
c     1.4
c     Exp
c
c*********************************************************************
      use comm_module 

      implicit none

c      include "comms.inc"

      integer i,n,ierr,mxnode

CMPIU      define MPI_COMM_SIZE MPI_COMM_SIZE_

      call MPI_COMM_SIZE(pydlpoly_comm, mxnode ,ierr)
      n=1
      nodedim = -1
      do i=0,16

         if(n.eq.mxnode)nodedim=i
         n=2*n

      enddo

      return
      end

      integer function numnodes()

c*********************************************************************
c
c     calculate number of nodes
c
c     MPI version - t.forester may 1995
c
c     wl
c     2008/01/14 13:33:07
c     1.4
c     Exp
c
c*********************************************************************
      use comm_module 

      implicit none

c      include "comms.inc"

      integer ierr

CMPIU      define MPI_COMM_SIZE MPI_COMM_SIZE_

      call MPI_COMM_SIZE(pydlpoly_comm, numnodes, ierr)

      return
      end

      subroutine csend(tagmsg,buf,length,pe,idum)

c*********************************************************************
c
c     Intel-like  csend (double precision)
c
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c
c     wl
c     2008/01/14 13:33:07
c     1.4
c     Exp
c
c*********************************************************************
      use comm_module 

      implicit none

c      include "comms.inc"

      integer tagmsg,length,pe,idum

      integer ierr
      real(8) buf(*)

CMPIU      define MPI_send MPI_send_

      call MPI_send(buf,length,MPI_DOUBLE_PRECISION,pe,tagmsg,
     x     pydlpoly_comm,ierr)

      return
      end

      subroutine crecv(tagmsg,buf,length)

c*********************************************************************
c
c     Intel-like  crecv (double precision)
c
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c
c     wl
c     2008/01/14 13:33:07
c     1.4
c     Exp
c
c*********************************************************************
      use comm_module 

      implicit none

c      include "comms.inc"

      integer tagmsg,length

      integer ierr
      integer status(MPI_STATUS_SIZE)
      real(8) buf(*)

CMPIU      define MPI_RECV MPI_RECV_

      call MPI_RECV(buf,length,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,
     x     tagmsg,pydlpoly_comm,status,ierr)

      return 
      end

      subroutine gisum(aaa,nnn,bbb)

c***********************************************************************
c     
c     dl_poly global summation subroutine for hypercube - MPI version
c     integer version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith march 1992.
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c     
c     wl
c     2008/01/14 13:33:07
c     1.4
c     Exp
c
c***********************************************************************
      
      use setup_module
      use comm_module 

      implicit none

      integer nnn,i,ierror,iii,kk,k,k0,k1,k2,msg1,msg2
      integer aaa(nnn),bbb(nnn)

c      include "comms.inc"

      integer status(MPI_STATUS_SIZE)

CMPIU      define MPI_allreduce MPI_allreduce_

      call MPI_allreduce(aaa,bbb,nnn,MPI_INTEGER,
     x  MPI_SUM,pydlpoly_comm,ierror)

      do i = 1,nnn
        aaa(i) = bbb(i)
      enddo

      return
      end

      subroutine gdsum(aaa,nnn,bbb)

c***********************************************************************
c     
c     dl_poly global summation subroutine for MPI - hypercube assumed
c     double precision version
c     
c     copyright - daresbury laboratory 1995
c     author    - w. smith march 1992.
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c     
c     wl
c     2008/01/14 13:33:07
c     1.4
c     Exp
c
c***********************************************************************
      use comm_module 

      implicit none

      integer nnn,i,iii,kk,k1,k2,ierror
      real(8) aaa(nnn),bbb(nnn)

c      include "comms.inc"

      integer status(MPI_STATUS_SIZE)

CMPIU      define MPI_allreduce MPI_allreduce_

      call MPI_allreduce(aaa,bbb,nnn,MPI_DOUBLE_PRECISION,
     x  MPI_SUM,pydlpoly_comm,ierror)

        do i = 1,nnn
          aaa(i) = bbb(i)
        enddo

      return
      end

      subroutine gimax(aaa,nnn,bbb)

c***********************************************************************
c     
c     dl_poly global maximum subroutine for hypercube - MPI version
c     integer version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith march 1992.
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c     
c     wl
c     2008/01/14 13:33:07
c     1.4
c     Exp
c
c***********************************************************************
      
      use setup_module
      use comm_module 

      implicit none

      integer nnn,i,iii,kk,k1,k2,k,k0msg1,msg2,ierror
      integer aaa(nnn),bbb(nnn)

c      include "comms.inc"

      integer status(MPI_STATUS_SIZE)
CMPIU      define MPI_allreduce MPI_allreduce_
      
      call MPI_allreduce(aaa,bbb,nnn,MPI_INTEGER,
     x   MPI_MAX,pydlpoly_comm,ierror)
      
      do i = 1,nnn
        aaa(i) = bbb(i)
      enddo

      return
      end

      subroutine gstate(check)

c***********************************************************************
c     
c     dl_poly global status subroutine : gisum version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith       march 1992
c     MPI version -  t. forester may 1995
c     
c     wl
c     1.4
c     Exp
c***********************************************************************


      implicit none

      logical check
      integer i,idum

      i = 0
      if(.not.check) i = 1

      call gisum(i,1,idum)
      
      check = (i.eq.0)

      return
      end

      subroutine gsync()

c*********************************************************************
c     
c     barrier / synchronization routine
c
c     MPI version - t.forester may 1995
c     CPP version - w.smith
c
c     wl
c     2008/01/14 13:33:07
c     1.4
c     Exp
c
c*********************************************************************
      use comm_module 

      implicit none

      integer ierr

c      include "comms.inc"


CMPIU      define MPI_BARRIER MPI_BARRIER_

      call  MPI_BARRIER(pydlpoly_comm,ierr)

      return
      end

c*********************************************************************
c     
c     broadcast double array
c
c     C. Spickermann/R. Schmid
c
c
c*********************************************************************


      subroutine bcast_d(d_array,n)
      use comm_module 
      implicit none
c      include "comms.inc"

Cf2py intent(in) n
      integer               :: n
Cf2py intent(inout) d_array
Cf2py depend(n) d_array
      real(8),dimension(n)  :: d_array
      integer               :: ierror
      
      call MPI_bcast(d_array,n,MPI_DOUBLE_PRECISION,
     x  0,pydlpoly_comm,ierror)
      
      end subroutine
ccs
ccs
cRS
      subroutine bcast_i(i_array,n)
      use comm_module 
      implicit none
c      include "comms.inc"

Cf2py intent(in) n
      integer               :: n
Cf2py intent(inout) i_array
Cf2py depend(n) i_array
      integer,dimension(n)  :: i_array
      integer               :: ierror
      
      call MPI_bcast(i_array,n,MPI_INTEGER,
     x  0,pydlpoly_comm,ierror)
      
      end subroutine

      subroutine bcast_c(c_array,n)
      use comm_module 
      implicit none
c      include "comms.inc"

Cf2py intent(inout) c_array
      character*(*)         :: c_array
Cf2py intent(in) n
      integer               :: n
      integer               :: ierror
      
      call MPI_bcast(c_array,n,MPI_CHARACTER,
     x  0,pydlpoly_comm,ierror)
      
      end subroutine


      subroutine exitcomms(ierr)

c*********************************************************************
c
c     exitcomms: exit from communication harness
c
c     MPI version - t.forester may 1995
c     CPP version - w.smith may 1995
c
c     wl
c     2008/01/14 13:33:07
c     1.4
c     Exp
c
c*********************************************************************

      use comm_module 
      implicit none

c      include "comms.inc"

      integer ierr
      external pyerror
CMPIU      define MPI_FINALIZE MPI_FINALIZE_
c       write (*,*) "calling MPI_Finalize ... Shutting down"

c      call MPI_FINALIZE(ierr)
cRS avoid this for a clear return to python in pydlpoly
c      call exit(0)

c     call python callback with the ierr that was passed to see what we can do now
      error_code = ierr
      call pyerror()

      return
      end
