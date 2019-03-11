      module angles_module

c***********************************************************************
c     
c     dl_poly module for defining valence angle potentials
c     copyright - daresbury laboratory
c     
c     author    - w. smith    sep 2003
c     adapted for solvation, free energy and excitation
c               - p.-a. cazade oct 2007
c
c     wl
c     2008/12/23 10:29:11
c     1.6
c     Exp
c     
c***********************************************************************

      use config_module
      use error_module
      use parse_module
      use property_module
      use setup_module
      use site_module
      use solvation_module
      use utility_module

      implicit none

      real(8), allocatable :: prmang(:,:)
      integer, allocatable :: listang(:,:)
      integer, allocatable :: numang(:),keyang(:),lstang(:,:)

      save prmang,listang,numang,keyang,lstang

      contains
      
      subroutine alloc_ang_arrays(idnode)
      
c***********************************************************************
c     
c     dl_poly subroutine for defining valence angle potential arrays
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     
c     wl
c     2008/12/23 10:29:11
c     1.6
c     Exp
c     
c***********************************************************************
      
      implicit none

      integer i,fail,idnode
      dimension fail(5)

      do i=1,5
        fail(i)=0
      enddo

      allocate (prmang(mxtang,mxpang),stat=fail(1))
      allocate (numang(mxtmls),stat=fail(2))
      allocate (keyang(mxtang),stat=fail(3))
      allocate (lstang(mxtang,3),stat=fail(4))
      allocate (listang(mxangl,4),stat=fail(5))

      do i=1,5
        if(fail(i).gt.0)call error(idnode,1010)
      enddo

      do i=1,mxtmls
         numang(i)=0
      enddo

      end subroutine alloc_ang_arrays

      subroutine define_angles
     x  (safe,idnode,itmols,nangle,nsite,engunit)

c***********************************************************************
c     
c     dl_poly subroutine for defining bond angles
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     
c     wl
c     2008/12/23 10:29:11
c     1.6
c     Exp
c     
c***********************************************************************

      implicit none

      logical safe
      character*8 keyword
      character*1 message(80)
      integer idnode,itmols,nangle,nsite,ntmp,i,iang,iang1
      integer idum,iatm1,iatm2,iatm3,isite1,isite2,isite3,ia,ja
      real(8) engunit

      ntmp=intstr(record,lenrec,idum)
      numang(itmols)=numang(itmols)+ntmp
      if(idnode.eq.0)then
        write(nrite,"(/,1x,'number of bond angles',
     x    10x,i10)")ntmp
        write(nrite,"(/,/,1x,'bond angle details:',
     x    /,/,21x,7x,'key',5x,'index',5x,'index',5x,
     x    'index',5x,'f-const',7x,'angle',/)")
      endif
      
      iang1=numang(itmols)
      do iang=1,iang1

c     read bond angle potential parameters
        
        call getrec(safe,idnode,nfield)
        if(.not.safe)return

        call copystring(record,message,80)
        call lowcase(record,4)
        call getword(keyword,record,4,lenrec)
        iatm1=intstr(record,lenrec,idum)
        iatm2=intstr(record,lenrec,idum)
        iatm3=intstr(record,lenrec,idum)

c     test for frozen atom pairs

        isite1=nsite-numsit(itmols)+iatm1
        isite2=nsite-numsit(itmols)+iatm2
        isite3=nsite-numsit(itmols)+iatm3

        if(lfzsit(isite1)*lfzsit(isite2)*
     x    lfzsit(isite3).ne.0)then
          
          numang(itmols)=numang(itmols)-1
          if(idnode.eq.0)write(nrite,'(14x,a16,40a1)')
     x      '*** frozen *** ',(message(i),i=1,40)

        else

          nangle=nangle+1
          
          if(nangle.gt.mxtang)call error(idnode,50)
          
          if(keyword(1:4).eq.'harm')then
            keyang(nangle)=1
          elseif(keyword(1:4).eq.'-hrm')then
            keyang(nangle)=-1
          elseif(keyword(1:4).eq.'quar')then
            keyang(nangle)=2
          elseif(keyword(1:4).eq.'-qur')then
            keyang(nangle)=-2
          elseif(keyword(1:4).eq.'thrm')then
            keyang(nangle)=3
          elseif(keyword(1:4).eq.'-thm')then
            keyang(nangle)=-3
          elseif(keyword(1:4).eq.'shrm')then
            keyang(nangle)=4
          elseif(keyword(1:4).eq.'-shm')then
            keyang(nangle)=-4
          elseif(keyword(1:4).eq.'bvs1')then
            keyang(nangle)=5
          elseif(keyword(1:4).eq.'-bv1')then
            keyang(nangle)=-5
          elseif(keyword(1:4).eq.'bvs2')then
            keyang(nangle)=6
          elseif(keyword(1:4).eq.'-bv2')then
            keyang(nangle)=-6
          elseif(keyword(1:4).eq.'hcos')then
            keyang(nangle)=7
          elseif(keyword(1:4).eq.'-hcs')then
            keyang(nangle)=-7
          elseif(keyword(1:4).eq.'cos ')then
            keyang(nangle)=8
          elseif(keyword(1:4).eq.'-cos')then
            keyang(nangle)=-8
          elseif(keyword(1:4).eq.'mmsb')then
            keyang(nangle)=9
          elseif(keyword(1:4).eq.'-msb')then
            keyang(nangle)=-9
          elseif(keyword(1:4).eq.'stst') then
            keyang(nangle)=10
          elseif(keyword(1:4).eq.'-sts') then
            keyang(nangle)=-10
          elseif(keyword(1:4).eq.'stbe') then
            keyang(nangle)=11
          elseif(keyword(1:4).eq.'-stb') then
            keyang(nangle)=-11
          elseif(keyword(1:4).eq.'cmps') then
            keyang(nangle)=12
          elseif(keyword(1:4).eq.'-cmp') then
            keyang(nangle)=-12
cRS MM3 sextic angle with predefined higher order terms
          elseif(keyword(1:4).eq.'mm3a') then
            keyang(nangle)=13
          elseif(keyword(1:4).eq.'-mm3') then
            keyang(nangle)=-13
cRS  Vedani and Huhta double cosine (square planar and octahedral L-Metal-L)
          elseif(keyword(1:4).eq.'cos2') then
            keyang(nangle)=14
          elseif(keyword(1:4).eq.'-cs2') then
            keyang(nangle)=-14
cJPD coupling potential between cos2 and stretch like cmps
          elseif(keyword(1:4).eq.'cost') then
            keyang(nangle)=15
          elseif(keyword(1:4).eq.'-cot') then
            keyang(nangle)=-15
          elseif(keyword(1:4).eq.'aqua') then
            keyang(nangle)=16
cRS
          else
            if(idnode.eq.0)write(nrite,*)message
            call error(idnode,440)
          endif
          
          lstang(nangle,1)=iatm1
          lstang(nangle,2)=iatm2
          lstang(nangle,3)=iatm3
          prmang(nangle,1)=dblstr(record,lenrec,idum)
          prmang(nangle,2)=dblstr(record,lenrec,idum)
          prmang(nangle,3)=dblstr(record,lenrec,idum)
          prmang(nangle,4)=dblstr(record,lenrec,idum)
          prmang(nangle,5)=dblstr(record,lenrec,idum)
          prmang(nangle,6)=dblstr(record,lenrec,idum)
          prmang(nangle,7)=dblstr(record,lenrec,idum)
          
          if(idnode.eq.0)
     x      write(nrite,"(27x,a4,3i10,1p,e12.4,0p,9f12.6)")
     x      keyword(1:4),(lstang(nangle,ia),ia=1,3),
     x      (prmang(nangle,ja),ja=1,mxpang)
          
c     convert energies to internal units
          
          prmang(nangle,1)=prmang(nangle,1)*engunit
          if(abs(keyang(nangle)).eq.2)then
            prmang(nangle,3)=prmang(nangle,3)*engunit
            prmang(nangle,4)=prmang(nangle,4)*engunit
          elseif(abs(keyang(nangle)).eq.12)then
            prmang(nangle,2)=prmang(nangle,2)*engunit            
            prmang(nangle,3)=prmang(nangle,3)*engunit
cJPD coupling potential between cos2 and stretch like cmps
          elseif(abs(keyang(nangle)).eq.15)then
            prmang(nangle,2)=prmang(nangle,2)*engunit
            prmang(nangle,3)=prmang(nangle,3)*engunit
          endif
          
c     convert angles to radians
          
          if(abs(keyang(nangle)).eq.12)then
            prmang(nangle,4)=prmang(nangle,4)*(pi/180.d0)
          elseif(abs(keyang(nangle)).eq.15)then
            prmang(nangle,4)=prmang(nangle,4)*(pi/180.d0)
          elseif(abs(keyang(nangle)).eq.16)then
            prmang(nangle,2)=prmang(nangle,2)*(pi/180.d0) 
            prmang(nangle,3)=prmang(nangle,3)/(pi/180.d0) 
            prmang(nangle,4)=prmang(nangle,4)/(pi/180.d0)
          elseif(abs(keyang(nangle)).ne.10)then
            prmang(nangle,2)=prmang(nangle,2)*(pi/180.d0) 
          endif
          
        endif

      enddo
      
      return
      end subroutine define_angles

      subroutine angfrc
     x  (lsolva,lfree,lexcite,idnode,imcon,mxnode,ntangl,engang,virang)

c***********************************************************************
c     
c     dl_poly subroutine for calculating bond angle energy and 
c     force terms in molecular dynamics.
c     
c     replicated data - blocked version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith         may 1992
c     modified  - t. forester      feb 1993
c     modified  - t.forester       nov 1994 : block data
c     modified  - t.forester       may 1995 : stress tensor 
c     modified  - p.-a.cazade      oct 2007 : solvation etc.
c     modified  - r. schmid        sep 2009 : mm3 angle
c     modified  - r. schmid        apr 2010 : cos2 potential
c     
c     wl
c     2008/12/23 10:29:11
c     1.6
c     Exp
c     
c***********************************************************************
      
      implicit none
      
      logical safe,lsolva,lfree,lexcite,lselect
      integer idnode,mxnode,imcon,ntangl,fail1,fail2
      integer ii,iang1,iang2,i,ia,ib,ic,kk,keya
      real(8) engang,virang,theta,fxc,fyc,fzc,rab,xab
cRS
      real(8) dtheta, dtheta2
      real(8) cos2c, cos2c2, cos2s, cos2s2
cRS
      real(8) strs(6)
      real(8) yab,zab,rbc,xbc,ybc,zbc,sint,cost,pterm,vterm
      real(8) gamma,gamsa,gamsc,rrbc,rrab,fxa,fya,fza
      real(8), allocatable :: xdab(:),ydab(:),zdab(:)
      real(8), allocatable :: xdbc(:),ydbc(:),zdbc(:)

cRS    mm3 sextic potential with fixed higher order terms
cRS     cubic:    -0.014
cRS     quartic:   0.000056
cRS     pentic:   -0.0000007
cRS     sextic:    0.000000022        

       real(8), parameter :: mm3cubic = -0.014d0/(pi/180.0d0)
       real(8), parameter :: mm3quart = 0.000056d0/(pi/180.0d0)**2
       real(8), parameter :: mm3pent  = -0.0000007d0/(pi/180.0d0)**3
       real(8), parameter :: mm3sext  = 0.000000022d0/(pi/180.0d0)**4

c     define angular potential function and derivative
c     using the parameters in array prmang
CVAM      
CVAM      call VTBEGIN(22, ierr)
CVAM
      allocate (xdab(msbad),ydab(msbad),zdab(msbad),stat=fail1)
      allocate (xdbc(msbad),ydbc(msbad),zdbc(msbad),stat=fail2)
      if(fail1.ne.fail2)call error(idnode,1020)

c     flag for undefined potentials

      safe=.true.

c     check size of work arrays

      if((ntangl-mxnode+1)/mxnode.gt.msbad)call error(idnode,419)

c     block indices
      
      iang1=(idnode*ntangl)/mxnode+1
      iang2=((idnode+1)*ntangl)/mxnode

c     zero accumulators
      
      engang=0.d0
      virang=0.d0
      ang_fre=0.d0

      do i=1,6
        strs(i)=0.d0
      enddo
      
      if(lsolva)then
        
        lcomp(2)=.true.
        ang_sol(:)=0.d0
        if(lexcite)ang_exc(:)=0.d0
        
      endif
      
c     calculate atom separation vectors
      
      ii=0
      do i=iang1,iang2
        
        ii=ii+1
        
c     indices of bonded atoms
        
        ia=listang(ii,2)
        ib=listang(ii,3)
        ic=listang(ii,4)

c     components of first bond vector
        
        xdab(ii)=xxx(ia)-xxx(ib)
        ydab(ii)=yyy(ia)-yyy(ib)
        zdab(ii)=zzz(ia)-zzz(ib)

c     components of second bond vector
        
        xdbc(ii)=xxx(ic)-xxx(ib)
        ydbc(ii)=yyy(ic)-yyy(ib)
        zdbc(ii)=zzz(ic)-zzz(ib)
        
      enddo
      
c     periodic boundary condition
      
      call images(imcon,0,1,ii,cell,xdab,ydab,zdab)
      call images(imcon,0,1,ii,cell,xdbc,ydbc,zdbc)
      
c     loop over all specified angle potentials
      
      ii=0
      do i=iang1,iang2
        
        ii=ii+1

c     define components of first bond vector
        
        rab=sqrt(xdab(ii)**2+ydab(ii)**2+zdab(ii)**2)
        rrab=1.d0/rab
        
        xab=xdab(ii)*rrab
        yab=ydab(ii)*rrab
        zab=zdab(ii)*rrab

c     define components of second bond vector
        
        rbc=sqrt(xdbc(ii)**2+ydbc(ii)**2+zdbc(ii)**2)
        rrbc=1.d0/rbc
        
        xbc=xdbc(ii)*rrbc
        ybc=ydbc(ii)*rrbc
        zbc=zdbc(ii)*rrbc

c     index of potential function parameters
        
        kk=listang(ii,1)

c     determine bond angle and calculate potential energy
        
        cost=(xab*xbc+yab*ybc+zab*zbc)
cRS  make sure that cost is never below -1.0 to avoid NaN
        cost=max(-1.0d0, cost)
        theta=acos(cost)
        sint=max(1.d-8,sqrt(1.d0-cost**2))
        
        keya=abs(keyang(kk))
        
        if(keya.eq.1)then
          
c     Harmonic potential
          
          pterm=0.5d0*prmang(kk,1)*(theta-prmang(kk,2))**2
          gamma=prmang(kk,1)*(theta-prmang(kk,2))/sint
          vterm=0.d0
          gamsa=0.d0
          gamsc=0.d0
          
        elseif(keya.eq.2)then

c     Quartic potential
          
          pterm=0.5d0*prmang(kk,1)*(theta-prmang(kk,2))**2+
     x      1.d0/3.d0*prmang(kk,3)*(theta-prmang(kk,2))**3+
     x      0.25d0*prmang(kk,4)*(theta-prmang(kk,2))**4
          gamma=(prmang(kk,1)*(theta-prmang(kk,2))+
     x      prmang(kk,3)*(theta-prmang(kk,2))**2+
     x      prmang(kk,4)*(theta-prmang(kk,2))**3)/sint
          vterm=0.d0
          gamsa=0.d0
          gamsc=0.d0
          
        elseif(keya.eq.3)then

c     truncated Harmonic potential
      
          pterm=0.5d0*prmang(kk,1)*(theta-prmang(kk,2))**2*
     x      exp(-(rab**8+rbc**8)/prmang(kk,3)**8)
          gamma=prmang(kk,1)*(theta-prmang(kk,2))*
     x      exp(-(rab**8+rbc**8)/prmang(kk,3)**8)/sint
          vterm=-8.d0*pterm*(rab**8+rbc**8)/prmang(kk,3)**8
          gamsa=(8.d0*pterm/prmang(kk,3)**8)*rab**7
          gamsc=(8.d0*pterm/prmang(kk,3)**8)*rbc**7
          
        elseif(keya.eq.4)then

c     screened Harmonic potential
          
          pterm=0.5d0*prmang(kk,1)*(theta-prmang(kk,2))**2*
     x      exp(-(rab/prmang(kk,3)+rbc/prmang(kk,4)))
          gamma=prmang(kk,1)*(theta-prmang(kk,2))*
     x      exp(-(rab/prmang(kk,3)+rbc/prmang(kk,4)))/sint
          vterm=-pterm*(rab/prmang(kk,3)+rbc/prmang(kk,4))
          gamsa=(pterm/prmang(kk,3))
          gamsc=(pterm/prmang(kk,4))
          
        elseif(keya.eq.5)then

c     screened vessal potential (type 1)
      
          pterm=(prmang(kk,1)/(8.d0*(prmang(kk,2)-pi)**2)*
     x      (((prmang(kk,2)-pi)**2-(theta-pi)**2)**2))*
     x      exp(-(rab/prmang(kk,3)+rbc/prmang(kk,4)))
          gamma=(prmang(kk,1)/(2.d0*(prmang(kk,2)-pi)**2)*
     x      ((prmang(kk,2)-pi)**2-(theta-pi)**2)*(theta-pi))*
     x      exp(-(rab/prmang(kk,3)+rbc/prmang(kk,4)))/sint
          vterm=-pterm*(rab/prmang(kk,3)+rbc/prmang(kk,4))
          gamsa=(pterm/prmang(kk,3))
          gamsc=(pterm/prmang(kk,4))
          
        elseif(keya.eq.6)then

c     Truncated Vessal potential (type 2)
          
          pterm=prmang(kk,1)*(theta**prmang(kk,3)*(theta-prmang(kk,2))
     x      **2*(theta+prmang(kk,2)-2.d0*pi)**2-0.5d0*prmang(kk,3)*pi
     x      **(prmang(kk,3)-1.d0)*(theta-prmang(kk,2))**2*(pi-prmang(kk,
     x      2))**3)*exp(-(rab**8+rbc**8)/prmang(kk,4)**8)
          gamma=prmang(kk,1)*(theta**(prmang(kk,3)-1.d0)*(theta-prmang
     x      (kk,2))*(theta+prmang(kk,2)-2.d0*pi)*((prmang(kk,3)+4.d0)*
     x      theta**2-2.d0*pi*(prmang(kk,3)+2.d0)*theta+prmang(kk,3)*
     x      prmang(kk,2)*(2.d0*pi-prmang(kk,2)))-prmang(kk,3)*pi**
     x      (prmang(kk,3)-1.d0)*(theta-prmang(kk,2))*(pi-prmang(kk,2))
     x      **3)*exp(-(rab**8+rbc**8)/prmang(kk,4)**8)/sint
          vterm=-8.d0*pterm*(rab**8+rbc**8)/prmang(kk,4)**8
          gamsa=(8.d0*pterm/prmang(kk,4)**8)*rab**7
          gamsc=(8.d0*pterm/prmang(kk,4)**8)*rbc**7

        elseif(keya.eq.7)then
          
c     harmonic cosine potential
          
          pterm=0.5d0*prmang(kk,1)*(cos(theta)-cos(prmang(kk,2)))**2
          gamma=-prmang(kk,1)*(cos(theta)-cos(prmang(kk,2)))
          vterm=0.d0
          gamsa=0.d0
          gamsc=0.d0
          
        elseif(keya.eq.8)then
          
c     ordinary cosine potential
          
          pterm=prmang(kk,1)*(1+cos(prmang(kk,3)*theta-prmang(kk,2)))
          gamma=-prmang(kk,1)*prmang(kk,3)*sin(prmang(kk,3)*theta-
     x      prmang(kk,2))/sint
          vterm=0.d0
          gamsa=0.d0
          gamsc=0.d0
          
        elseif(keya.eq.9)then

c     MM3 stretch-bend potential

          pterm=prmang(kk,1)*(theta-prmang(kk,2))*
     x      (rab-prmang(kk,3))*(rbc-prmang(kk,4))
          gamma=prmang(kk,1)*(rab-prmang(kk,3))*(rbc-
     x      prmang(kk,4))/sint
          gamsa=-prmang(kk,1)*(theta-prmang(kk,2))*(rbc-prmang(kk,4))
          gamsc=-prmang(kk,1)*(theta-prmang(kk,2))*(rab-prmang(kk,3))
          vterm=-(gamsa*rab+gamsc*rbc)

        elseif(keya.eq.10)then

c     compass stretch-stretch potential

          pterm=prmang(kk,1)*(rab-prmang(kk,2))*(rbc-prmang(kk,3))
          gamma=0.d0
          gamsa=-prmang(kk,1)*(rbc-prmang(kk,3))
          gamsc=-prmang(kk,1)*(rab-prmang(kk,2))
          vterm=-(gamsa*rab+gamsc*rbc)

        elseif(keya.eq.11)then

c     compass stretch-bend potential

          pterm=prmang(kk,1)*(theta-prmang(kk,2))*(rab-prmang(kk,3))
          gamma=prmang(kk,1)*(rab-prmang(kk,3))/sint
          gamsa=-prmang(kk,1)*(theta-prmang(kk,2))
          gamsc=0.d0
          vterm=-gamsa*rab

        elseif(keya.eq.12)then

c     combined compass angle potential with 3 coupling terms

          pterm=prmang(kk,1)*(rab-prmang(kk,5))*(rbc-prmang(kk,6))+
     x          (theta-prmang(kk,4))*(prmang(kk,2)*(rab-prmang(kk,5))+
     x          prmang(kk,3)*(rbc-prmang(kk,6)))
          gamma=(prmang(kk,2)*(rab-prmang(kk,5))+
     x           prmang(kk,3)*(rbc-prmang(kk,6)))/sint
          gamsa=-prmang(kk,2)*(theta-prmang(kk,4))-
     x           prmang(kk,1)*(rbc-prmang(kk,6))
          gamsc=-prmang(kk,3)*(theta-prmang(kk,4))-
     x           prmang(kk,1)*(rab-prmang(kk,5))
          vterm=-(gamsa*rab+gamsc*rbc)

        elseif(keya.eq.13)then
          
cRS    mm3 sextic potential with fixed higher order terms
cRS     cubic:    -0.014
cRS     quartic:   0.000056
cRS     pentic:   -0.0000007
cRS     sextic:    0.000000022        
          dtheta = theta-prmang(kk,2)
          dtheta2 = dtheta*dtheta
          pterm=0.5d0*prmang(kk,1)*dtheta2 * (1.0d0 + mm3cubic*dtheta
     x      + mm3quart*dtheta2 + mm3pent*dtheta*dtheta2
     x       + mm3sext*dtheta2*dtheta2)
          gamma=0.5d0*prmang(kk,1)*dtheta/sint *  (2.0d0
     x       + mm3cubic*3.0d0*dtheta
     x       + mm3quart*4.0d0*dtheta2
     x       + mm3pent *5.0d0*dtheta*dtheta2
     x       + mm3sext *6.0d0*dtheta2*dtheta2)
          vterm=0.d0
          gamsa=0.d0
          gamsc=0.d0
        
        elseif(keya.eq.14)then
          
cRS       double cosine potential according to Vedani and Huhta
cRS          E = A [1+cos(a+abs(a0))]*[1+cos(2a+a0)]
cRS          minima at 90, 180, 270 deg for octahedral systems 
cRS                             (for a0=0.0)
          cos2c = cos(theta+abs(prmang(kk,2)))
          cos2c2= cos(2.0d0*theta+prmang(kk,2))
          cos2s = sin(theta+abs(prmang(kk,2)))
          cos2s2= sin(2.0d0*theta+prmang(kk,2))
          pterm= prmang(kk,1)*(1.0d0+cos2c)*(1.0d0+cos2c2)
          gamma= -prmang(kk,1)/sint*
     x      ((1.0d0+cos2c)*2.0d0*cos2s2 + (1.0d0+cos2c2)*cos2s)
          vterm=0.d0
          gamsa=0.d0
          gamsc=0.d0
          
        elseif(keya.eq.15)then

cJPD coupling potential between cos2 and stretch like cmps
        
c          pterm=prmang(kk,1)*(rab-prmang(kk,6))*(rbc-prmang(kk,7))+
c     x          (1+cos(prmang(kk,5)*theta + prmang(kk,4)))*
c     x          (prmang(kk,2)*(rab-prmang(kk,6))+
c     x          prmang(kk,3)*(rbc-prmang(kk,7)))
c          gamma=(-sin(prmang(kk,5)*theta+prmang(kk,4))*prmang(kk,5)*
c     x           (prmang(kk,2)*(rab-prmang(kk,6))+
c     x           prmang(kk,3)*(rbc-prmang(kk,7))))/sint
c          gamsa=-(1+cos(prmang(kk,5)*theta+prmang(kk,4)))*prmang(kk,2)-
c     x           prmang(kk,1)*(rbc-prmang(kk,7))
c          gamsc=-(1+cos(prmang(kk,5)*theta+prmang(kk,4)))*prmang(kk,3)-
c     x           prmang(kk,1)*(rab-prmang(kk,6))

          pterm=prmang(kk,1)*(rab-prmang(kk,6))*(rbc-prmang(kk,7))*
     x          sin(theta)+
     x          (1+cos(prmang(kk,5)*theta + prmang(kk,4)))*
     x          (prmang(kk,2)*(rab-prmang(kk,6))+
     x          prmang(kk,3)*(rbc-prmang(kk,7)))
          gamma=(-sin(prmang(kk,5)*theta+prmang(kk,4))*prmang(kk,5)*
     x           (prmang(kk,2)*(rab-prmang(kk,6))+
     x           prmang(kk,3)*(rbc-prmang(kk,7)))
     x           +prmang(kk,1)*(rab-prmang(kk,6))*(rbc-prmang(kk,7))*
     x           cos(theta))/sint
          gamsa=-(1+cos(prmang(kk,5)*theta+prmang(kk,4)))*prmang(kk,2)-
     x           prmang(kk,1)*(rbc-prmang(kk,7))*sin(theta)
          gamsc=-(1+cos(prmang(kk,5)*theta+prmang(kk,4)))*prmang(kk,3)-
     x           prmang(kk,1)*(rab-prmang(kk,6))*sin(theta)


        elseif(keya.eq.16)then
          
cJPD variable quartic potential in the allinger MM3 formulation
          dtheta = theta-prmang(kk,2)
          dtheta2 = dtheta*dtheta
          pterm=0.5d0*prmang(kk,1)*dtheta2 * (1.0d0-prmang(kk,3)*dtheta
     x      + prmang(kk,4)*(prmang(kk,3)**2)*dtheta2)
          gamma=0.5d0*prmang(kk,1)*dtheta/sint *  (2.0d0
     x       - prmang(kk,3)*3.0d0*dtheta
     x       + prmang(kk,4)*(prmang(kk,3)**2)*4.0d0*dtheta2)
          vterm=0.d0
          gamsa=0.d0
          gamsc=0.d0
          vterm=-(gamsa*rab+gamsc*rbc)
        else

c     undefined potential

          safe=.false.
          pterm=0.d0
          vterm=0.d0
          gamma=0.d0
          gamsa=0.d0
          gamsc=0.d0
          
        endif
        
c     indices of bonded atoms
        
        ia=listang(ii,2)
        ib=listang(ii,3)
        ic=listang(ii,4)

c     set selection control
        
        lselect=.true.
        
        if(lexcite)then
          
c     selected excitation option
        
          if((atm_fre(ia).ne.1).and.(atm_fre(ib).ne.1).and.
     x      (atm_fre(ic).ne.1))then
            
c     reset selection control
            
            lselect=(atm_fre(ia)+atm_fre(ib)+atm_fre(ic).eq.0)
            
            if(lsolva)then
              ang_exc(atmolt(ia))=ang_exc(atmolt(ia))+pterm
            endif
            
          endif
          
        elseif(lfree)then
          
c     selected free energy option
          
          if((atm_fre(ia).eq.1).or.(atm_fre(ib).eq.1).or.
     x        (atm_fre(ic).eq.1))then
            
c     set hamiltonian mixing parameter

            ang_fre=ang_fre-pterm
            pterm=lambda1*pterm
            vterm=lambda1*vterm
            gamma=lambda1*gamma
            gamsa=lambda1*gamsa
            gamsc=lambda1*gamsc
            
          elseif((atm_fre(ia).eq.2).or.(atm_fre(ib).eq.2).or.
     x        (atm_fre(ic).eq.2))then
            
c     set hamiltonian mixing parameter

            ang_fre=ang_fre+pterm
            pterm=lambda2*pterm
            vterm=lambda2*vterm
            gamma=lambda2*gamma
            gamsa=lambda2*gamsa
            gamsc=lambda2*gamsc
                        
          endif
          
        endif
        
c     sum potential energy and virial
        
        if(lselect)then
          
c     sum potential energy and virial
          
          engang=engang+pterm
          virang=virang+vterm
          
c     calculate solvation energy
        
          if(lsolva)then
            ang_sol(atmolt(ia))=ang_sol(atmolt(ia))+pterm
          endif
          
c     calculate atomic forces
        
          fxa=gamma*(xbc-xab*cost)*rrab+gamsa*xab
          fya=gamma*(ybc-yab*cost)*rrab+gamsa*yab
          fza=gamma*(zbc-zab*cost)*rrab+gamsa*zab
          
          fxc=gamma*(xab-xbc*cost)*rrbc+gamsc*xbc
          fyc=gamma*(yab-ybc*cost)*rrbc+gamsc*ybc
          fzc=gamma*(zab-zbc*cost)*rrbc+gamsc*zbc
          
c     sum atomic forces
          
          fxx(ia)=fxx(ia)+fxa
          fyy(ia)=fyy(ia)+fya
          fzz(ia)=fzz(ia)+fza
          
          fxx(ib)=fxx(ib)-fxa-fxc
          fyy(ib)=fyy(ib)-fya-fyc
          fzz(ib)=fzz(ib)-fza-fzc
          
          fxx(ic)=fxx(ic)+fxc
          fyy(ic)=fyy(ic)+fyc
          fzz(ic)=fzz(ic)+fzc
          
c     calculate stress tensor
          
          strs(1)=strs(1)+rab*xab*fxa+rbc*xbc*fxc
          strs(2)=strs(2)+rab*xab*fya+rbc*xbc*fyc
          strs(3)=strs(3)+rab*xab*fza+rbc*xbc*fzc
          strs(4)=strs(4)+rab*yab*fya+rbc*ybc*fyc
          strs(5)=strs(5)+rab*yab*fza+rbc*ybc*fzc
          strs(6)=strs(6)+rab*zab*fza+rbc*zbc*fzc
          
        endif
        
      enddo
      
c     check for undefined potentials

      if(mxnode.gt.1)call gstate(safe)
      if(.not.safe)call error(idnode,440)

c     complete stress tensor
        
      stress(1)=stress(1)+strs(1)
      stress(2)=stress(2)+strs(2)
      stress(3)=stress(3)+strs(3)
      stress(4)=stress(4)+strs(2)
      stress(5)=stress(5)+strs(4)
      stress(6)=stress(6)+strs(5)
      stress(7)=stress(7)+strs(3)
      stress(8)=stress(8)+strs(5)
      stress(9)=stress(9)+strs(6)
      
c     sum up contributions to potential and virial
      
      if(mxnode.gt.1)then

        buffer(1)=engang
        buffer(2)=virang
        buffer(3)=ang_fre
        call gdsum(buffer(1),3,buffer(4))
        engang=buffer(1)
        virang=buffer(2)
        ang_fre=buffer(3)

c     sum up solvation energies
        
        if(lsolva)then

          call gdsum(ang_sol(1),mxtmls,buffer(1))
          if(lexcite)call gdsum(ang_exc(1),mxtmls,buffer(1))
          
        endif
        
      endif

      
      deallocate (xdab,ydab,zdab,stat=fail1)
      deallocate (xdbc,ydbc,zdbc,stat=fail2)
CVAM
CVAM      call VTEND(22, ierr)
CVAM
      return
      end subroutine angfrc
      
      end module angles_module
