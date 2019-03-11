      module config_module

c***********************************************************************
c     
c     dl_poly module for defining simulation configuration data
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     
c     wl
c     2007/07/31 10:02:56
c     1.3
c     Exp
c     
c***********************************************************************

      use setup_module
      use error_module
      
      implicit none

      character*1 cfgname(80)
      character*1 sysname(80)
      real(8) cell(9),rcell(9),celprp(10)
      real(8) eta(9),stress(9),stresl(9),strcns(9),strbod(9)
      
      character*8, allocatable :: atmnam(:)
      real(8), allocatable :: xxx(:),yyy(:),zzz(:)
      real(8), allocatable :: vxx(:),vyy(:),vzz(:)
      real(8), allocatable :: fxx(:),fyy(:),fzz(:)
      real(8), allocatable :: flx(:),fly(:),flz(:)
ccs   sigma array for width of Gaussians added
      real(8), allocatable :: chge(:),weight(:),rmass(:),qsigma(:)
      real(8), allocatable :: qpsigma(:)
      integer, allocatable :: spindex(:)
c     potential array for storing the dE/dq_i potential terms
      real(8), allocatable   :: pot(:)
c     Hessian matrix for QEq Newton-Raphson optimization
      real(8), allocatable :: q_hess(:,:)
      real(8)                :: enc_self
      real(8)                :: tval
c     logical indicating that self interaction and dipole correction
c     have to be updated
      logical :: recharge,time
ccs
crs
c     array for core charges
      real(8), allocatable  ::  core_chge(:)
      real(8),allocatable   ::  core_pot(:)
      real(8)               ::  enc_core, enc_coreval, enc_val
      logical               ::  use_core_charge, get_core_pot      
crs

      integer, allocatable :: ltype(:),lstfrz(:)
      integer, allocatable :: neulst(:),lstneu(:)
      integer, allocatable :: lentry(:),list(:,:)
      integer, allocatable :: lstout(:),link(:)
      integer, allocatable :: lct(:),lst(:)

      real(8), allocatable :: buffer(:)

crs    new arrays to track the image of an atom
      logical  :: ltrackimg
      integer, allocatable  ::  imgx(:), imgy(:), imgz(:)
      integer, allocatable  ::  ibuffer(:)
      integer               ::  mxibuff

crs    ignore atoms in polarization
      integer, allocatable  :: pol_ignore(:) 


      save atmnam,neulst,lstneu,cfgname,sysname
      save cell,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
ccs
      save buffer,weight,chge,qsigma,q_hess,pot,ltype,lstfrz,flx,fly,flz
      save enc_self
ccs
      save lentry,list,lstout,link,lct,lst,celprp,rmass
      save eta,stress,stresl,strcns,rcell
ccs
      data time,recharge/.false.,.false./
      data enc_self/0.0d0/
ccs
crs
      save core_chge, core_pot, use_core_charge, get_core_pot
      save enc_core,enc_coreval,enc_val
      save qpsigma, spindex
      data use_core_charge,get_core_pot/.false.,.false./
crs

crs   new arrays to track the image of an atom
      data ltrackimg/.false./
      save ltrackimg
      save imgx, imgy, imgz
      save ibuffer
      save mxibuff

      save pol_ignore
crs

      contains
      
      subroutine alloc_config_arrays(idnode)

c***********************************************************************
c     
c     dl_poly subroutine for defining simulation configuration arrays
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     
c     wl
c     2007/07/31 10:02:56
c     1.3
c     Exp
c     
c***********************************************************************
      
      integer, parameter :: nnn=27

      integer i,fail,idnode
      dimension fail(nnn)

      do i=1,nnn
        fail(i)=0
      enddo

      allocate (xxx(mxatms),stat=fail(1))
      allocate (yyy(mxatms),stat=fail(2))
      allocate (zzz(mxatms),stat=fail(3))
      allocate (vxx(mxatms),stat=fail(4))
      allocate (vyy(mxatms),stat=fail(5))
      allocate (vzz(mxatms),stat=fail(6))
      allocate (fxx(mxatms),stat=fail(7))
      allocate (fyy(mxatms),stat=fail(8))
      allocate (fzz(mxatms),stat=fail(9))
      allocate (weight(mxatms),stat=fail(11))
      allocate (chge(mxatms),stat=fail(12))
ccs
      if(lgauss)then
        allocate (qsigma(mxatms),stat=fail(12))
        !if(qeq_iter)then
        allocate (pot(mxatms),stat=fail(12))
        !if(calc_hess)
        allocate (q_hess(mxatms,mxatms),stat=fail(12))
        allocate (core_chge(mxatms), stat=fail(12))
        allocate (core_pot(mxatms), stat=fail(12))
      endif
ccs
      allocate (ltype(mxatms),stat=fail(13))
      allocate (lstfrz(mxatms),stat=fail(14))
      allocate (flx(mxatms),stat=fail(15))
      allocate (fly(mxatms),stat=fail(16))
      allocate (flz(mxatms),stat=fail(17))
      allocate (atmnam(mxatms),stat=fail(18))
      allocate (neulst(mxneut),stat=fail(19))
      allocate (lstneu(mxatms),stat=fail(20))
      allocate (lstout(mxatms),stat=fail(21))
      allocate (lentry(msatms),stat=fail(22))
      allocate (list(msatms,mxlist),stat=fail(23))
      allocate (link(mxatms),stat=fail(24))
      allocate (lct(mxcell),stat=fail(25))
      allocate (lst(mxcell),stat=fail(26))
      allocate (rmass(mxatms),stat=fail(27))
      allocate (buffer(mxbuff),stat=fail(10))

      do i=1,nnn
        if(fail(i).gt.0)then
          if(idnode.eq.0)write(nrite,'(10i5)')fail
          call error(idnode,1000)
        endif
      enddo

cRS
c      write (*,*) "DEBUG: this is node  ", idnode
c      write (*,*) "DEBUG: sizes mxatms=",mxatms
c      write (*,*) "DEBUG: sizes mxatms=",msatms
c      write (*,*) "DEBUG: sizes mxlist=",mxlist

      end subroutine alloc_config_arrays


      subroutine alloc_trackimg_arrays(idnode,mxnode)

        integer :: idnode, mxnode
        integer :: fail

        ltrackimg = .true.

        mxibuff = 6*(mxatms+mxnode-1)/mxnode

        fail = 0
        allocate (imgx(mxatms), stat=fail)
        allocate (imgy(mxatms), stat=fail)
        allocate (imgz(mxatms), stat=fail)
        allocate (ibuffer(mxibuff), stat=fail)
        allocate (pol_ignore(mxatms), stat=fail)
        imgx(:) = 0
        imgy(:) = 0
        imgz(:) = 0
        pol_ignore(:) = 0
        if (fail.gt.0)then
          if(idnode.eq.0)write(nrite,'(10i5)')fail
          call error(idnode, 1000)
        endif
      end subroutine alloc_trackimg_arrays

      end module config_module
