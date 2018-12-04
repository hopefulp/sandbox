      module vdw_module

c***********************************************************************
c     
c     PYDLPOLY
c
c     vdw module derived from DL_Poly Classic
c
c     new features (incomplete):
c      - switching at cutoff to be smooth also in force
c      - remove WRONG computation of gradients from lookuptable (use cubic spline instead
c                     and a gradient computed exactly for the spline)
c      - allow for lambda switching
c      - remove neutral group implementation, solvation and free energy stuff (not used in pydlpoly .. done via molecules)
c      - added longrange_correction to be recomputed every step
c
c
c     dl_poly module for defining van der waals potential arrays
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     adapted for solvation, free energy and excitation
c               - p.-a. cazade oct 2007
c     
c     wl
c     2009/01/13 11:22:05
c     1.7
c     Exp
c     
c***********************************************************************

      use config_module
      use error_module
      use pair_module
      use parse_module
      use setup_module
      use site_module
      use solvation_module
      use utility_module
cRS for switchable molecule implementation
      use molecule_module 

      implicit none

      integer, allocatable :: ltpvdw(:),lstvdw(:)
      real(8), allocatable :: vvv(:,:),ggg(:,:),prmvdw(:,:)

      save ltpvdw,lstvdw,prmvdw,vvv,ggg

cRS  for switching functions ( .. see Steinbach, Brooks JCC 1994, 15, 667-683.)
      real(8)  vdwswitch
      logical  lvdwswitch
cRS  for alternative force calculations: lvdwintder: derivative directly from 3poitn interpolation
c                                        lvdwspline: use cubic splines instead and get derivative directly
c                                                 in this case ggg contains the second deriv of vvv
      logical lvdwintder, lvdwspline

      save  vdwswitch,lvdwswitch,lvdwintder,lvdwspline

c     switch is the factor for the switching to start (sw_on = switch*rcut)
c     as a default switch is 1.0 which means switching is off
c     Note: switch will be set in the define_system_module

      real(8), allocatable :: vvv2(:,:)
      save vvv2

cRS   params for ed3 (exp+d3)
      real(8)   dispSr6, dispS8
      save      dispSr6, dispS8
cRS

cRS   control variables for short range damping
c     the idea is to replace the very short range potential by a parabola with the proper value and derivative
c     at the low r cutoff which is determined automatically from an energy condition: if the repulsive energy is
c     above a certain cutoff value (probably determined by kT somehow) it will be truncated
c
      logical lsrdamp
      real(8) srdamp_ener
      save    lsrdamp, srdamp_ener
cRS

cRS   control and data variables for the splitting of repulsion and dispersion
c     as a policy: the attractive part is kept on vvv/vvv2 and the repulsion is going on vvvr and vvvr2

      logical lspltvdw
      data lspltvdw/.false./
      real(8) vdw_ener_rep, vdw_ener_dis
      real(8), allocatable :: vvvr(:,:), vvvr2(:,:)
      save lspltvdw, vdw_ener_dis, vdw_ener_rep
      save vvvr, vvvr2
cRS

cRS   switch to bypass calculation of bonded terms
      logical skip_bonded
      save    skip_bonded
      data    skip_bonded/.false./
cRS

      contains
      
      subroutine alloc_vdw_arrays(idnode)

      implicit none

      integer, parameter :: nnn=8

      integer i,fail,idnode
      dimension fail(nnn)

      do i=1,nnn
        fail(i)=0
      enddo

      allocate (ltpvdw(mxvdw),stat=fail(1))
      allocate (lstvdw(mxvdw),stat=fail(2))
      allocate (prmvdw(mxvdw,mxpvdw),stat=fail(3))
      allocate (vvv(mxgrid,mxvdw),stat=fail(4))
      allocate (ggg(mxgrid,mxvdw),stat=fail(5))

      allocate (vvv2(mxgrid,mxvdw),stat=fail(6))

      if (lspltvdw) then
c        write (*,*) "Using a split vdw potential with repuslion on vvvr"
        allocate (vvvr(mxgrid,mxvdw),stat=fail(7))
        allocate (vvvr2(mxgrid,mxvdw),stat=fail(8))        
      endif

      do i=1,nnn
        if(fail(i).gt.0)call error(idnode,1014)
      enddo

      end subroutine alloc_vdw_arrays

      subroutine define_van_der_waals
     x  (safe,ltable,lunits,lmols,idnode,ntpvdw,
     x  ntpatm,keyfce,dlrpot,rvdw,engunit)

c***********************************************************************
c     
c     dl_poly subroutine for defining van der Waals potentials
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     
c     wl
c     2009/01/13 11:22:05
c     1.7
c     Exp
c     
c***********************************************************************

      implicit none

      logical safe,ltable,lunits,lmols
      character*1 message(80)
      character*8 atom1,atom2,keyword

      integer ntpvdw,ntpatm,keyfce,fail,idum,ivdw
      integer itpvdw,keypot,numpar,katom1,katom2,jtpatm,keyvdw,i
      integer ntab,idnode,j
      real(8) dlrpot,rvdw,engunit
      real(8), allocatable :: parpot(:)

      allocate (parpot(mxpvdw),stat=fail)

      ntpvdw=intstr(record,lenrec,idum)

      ltable=findstring('table',record,idum)
      
      if(idnode.eq.0) then
        
        write(nrite,"(/,/,1x,'number of specified pair ',
     x    'potentials',i10)") ntpvdw

        write(nrite,"(/,/,16x,'atom 1  ','atom 2  ',3x,
     x    ' key',30x,'parameters'/,/)")
        
      endif      

      if(ntpvdw.gt.mxvdw) call error(idnode,80)
      if(.not.lunits) call error(idnode,6)
      if(.not.lmols) call error(idnode,13)
      
      do ivdw=1,mxvdw
        
        lstvdw(ivdw)=0
        ltpvdw(ivdw)=-1
        
      enddo
      
      do itpvdw=1,ntpvdw
        
        do i=1,mxpvdw
          parpot(i)=0.d0
        enddo
        
        call getrec(safe,idnode,nfield)
        if(.not.safe)return

        call copystring(record,message,80)
        call getword(atom1,record,8,lenrec)
        call getword(atom2,record,8,lenrec)
        call lowcase(record,lenrec-16)
        call getword(keyword,record,4,lenrec)

        if(keyword(1:4).eq.'12-6') then
          keypot=1
          numpar=2
        elseif(keyword(1:4).eq.'lj  ') then
          keypot=2
          numpar=2
        elseif(keyword(1:4).eq.'nm  ') then
          keypot=3
          numpar=4
        elseif(keyword(1:4).eq.'buck') then
          keypot=4
          numpar=4
        elseif(keyword(1:4).eq.'ex6d') then
          keypot=10
          numpar=4
        elseif(keyword(1:4).eq.'bhm ') then
          keypot=5
          numpar=5
        elseif(keyword(1:4).eq.'hbnd') then
          keypot=6
          numpar=2
        elseif(keyword(1:4).eq.'snm ') then
          keypot=7
          numpar=5
        elseif(keyword(1:4).eq.'mors') then
          keypot=8
          numpar=3
        elseif(keyword(1:4).eq.'wca ') then
          keypot=9
          numpar=3
        elseif(keyword(1:4).eq.'ed3 ') then
cRS       use exponential repulsion plus grimme dispersion correction D3
          keypot=11
          numpar=5
        elseif(keyword(1:4).eq.'tab ') then
          keypot=0
          numpar=0
        else
          if(idnode.eq.0) write(nrite,*) message
          call error(idnode,452)
        endif

        parpot(1)=dblstr(record,lenrec,idum)
        parpot(2)=dblstr(record,lenrec,idum)
        parpot(3)=dblstr(record,lenrec,idum)
        parpot(4)=dblstr(record,lenrec,idum)
        parpot(5)=dblstr(record,lenrec,idum)
        
        if(idnode.eq.0) 
     x    write(nrite,"(16x,2a8,2x,a4,3x,1p,9e13.5)") 
     x    atom1,atom2,keyword(1:4),(parpot(j),j=1,numpar)
        
        katom1=0
        katom2=0
        
        do jtpatm=1,ntpatm

          if(atom1.eq.unqatm(jtpatm))katom1=jtpatm
          if(atom2.eq.unqatm(jtpatm))katom2=jtpatm
          
        enddo
        
        if(katom1.eq.0.or.katom2.eq.0) then
          call  error(idnode,81)
        endif
        
        keyvdw=loc2(katom1,katom2)

c     convert energies to internal unit

        if(keyvdw.gt.mxvdw) call error(idnode,82)
        
        parpot(1)=parpot(1)*engunit
        
        if(keypot.eq.1) then
          
          parpot(2)=parpot(2)*engunit
          
        else if(keypot.eq.4) then
          
          parpot(3)=parpot(3)*engunit
          
        else if(keypot.eq.5) then
          
          parpot(4)=parpot(4)*engunit
          parpot(5)=parpot(5)*engunit
          
        else if(keypot.eq.6) then
          
          parpot(2)=parpot(2)*engunit

        else if(keypot.eq.10) then
          
          parpot(3)=parpot(3)*engunit

        else if(keypot.eq.11) then

cRS        params:  A[energy] B[1/A] C6[energy] C8[energy] RAB[A] 
          parpot(3)=parpot(3)*engunit
          parpot(4)=parpot(4)*engunit
          
        endif

        ltable=(ltable.or.(keypot.eq.0))

        if(lstvdw(keyvdw).ne.0) call error(idnode,15)
        lstvdw(keyvdw)=itpvdw
        ltpvdw(itpvdw)=keypot
        
        do i=1,mxpvdw
          
          prmvdw(itpvdw,i)=parpot(i)
          
        enddo
        
      enddo

c     generate nonbonded force arrays

      if((ntpvdw.gt.0 .and. mod(keyfce,2).eq.1).or.(keyfce.eq.2))
     x  then
        
        call forgen(ltable,idnode,ntpvdw,dlrpot,rvdw)
        
        if(ltable)then
          
          call fortab
     x      (idnode,ntpvdw,ntpatm,dlrpot,rvdw,engunit)
          
        endif
        
      endif

c     check for unspecified atom-atom potentials
      
      ntab=(ntpatm*(ntpatm+1))/2
      
      if(ntpvdw.lt.ntab) then
        
        call warning(idnode,110,0.d0,0.d0,0.d0)

        if(mxvdw.le.ntpvdw) call error(idnode,82)

        do i=1,ntab
          
          if(lstvdw(i).eq.0)then
            
            lstvdw(i)=ntpvdw+1
            
          endif
          
        enddo

c     define zero potential for undefined interactions
        
        do i=1,mxgrid
          
          ggg(i,ntpvdw+1)=0.d0
          vvv(i,ntpvdw+1)=0.d0
          
        enddo
        
      endif

      deallocate (parpot,stat=fail)

      return
      end subroutine define_van_der_waals

      subroutine forgen(ltable,idnode,ntpvdw,dlrpot,rcut)

c***********************************************************************
c     
c     dl_poly subroutine for generating potential energy and 
c     force arrays for van der waals forces only
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith may 1992.
c     
c     wl
c     2009/01/13 11:22:05
c     1.7
c     Exp
c     
c***********************************************************************
      
      implicit none

      logical ltable
      integer i,ivdw,ntpvdw,idnode
      real(8) dlrpot,rcut,rrr,ann,amm,gam,bet,eps,rr0,aaa,bbb
      real(8) ccc,ddd,eee,sig,rho,rrc
      
cRS   for shifting
      real(8) c0, c1, c2, c3, c4, c5
      real(8) sw_on, sw_off, denom, sw_off2, sw_on2
      real(8) sw_value, sw_dval, vvv_temp, ggg_temp
      real(8) rr2,rr3,rr4,rr5
cRS   for damped buckingham and ed3
      real(8) rc, t, q14, disp, ddisp, expterm,dexpterm
      real(8) cc6,cc8,q6,q8,q6_14,q8_16,t6,t8,e6,e8,c6,c8,r2,r6,r8
cRS   for short range damping
      integer icut
      real(8) a, b, r_srcut, ggg_disp

      real(8), allocatable :: u(:)
      real(8) p
      integer fail

      allocate(u(mxgrid), stat=fail)
      
      
CVAM
CVAM      call VTBEGIN(140, ierr)
CVAM

c     define grid resolution for potential arrays
      
      dlrpot=rcut/dble(mxgrid-4)
      
cRS   compute switching polynomial
      if (.not. (vdwswitch .eq. 1.0d0)) then
          sw_off = rcut
          sw_on  = rcut*vdwswitch
          sw_off2 = sw_off*sw_off
          sw_on2  = sw_on*sw_on
          denom = (sw_off-sw_on)**5
          c0 = sw_off*sw_off2 * (sw_off2-5.0d0*sw_off*sw_on+
     x              10.0d0*sw_on2)/denom
          c1 = -30.0d0 * (sw_off2*sw_on2)/denom
          c2 = 30.0d0 * (sw_off2*sw_on + sw_off*sw_on2)/denom
          c3 = -10.0d0 * (sw_off2 + 4.0d0*sw_off*sw_on + sw_on2)/denom
          c4 = 15.0d0 * (sw_off+sw_on)/denom
          c5 = -6.0d0 / denom
      endif

c     construct arrays for all types of short ranged  potential
cRS   if lspltvdw is true then compute the repulsive part seperate and put it in vvvr
c       currently only implemented for lennard/jones, 12/6 and the exp6d of MOF-FF
c
c     NOTE: ggg is the gradient/r of the repulsive part in vvvr (use only with spline)
c           ggg is used to truncate the short range repulsion

      do ivdw=1,ntpvdw
        
        if(ltpvdw(ivdw).eq.1)then
          
c       12 - 6 potential
      
          aaa=prmvdw(ivdw,1)
          bbb=prmvdw(ivdw,2)
          
          if (lspltvdw) then

            do i=1,mxgrid            
              rrr=dble(i)*dlrpot
              vvv(i,ivdw) =-bbb/rrr**6
              vvvr(i,ivdw)=aaa/rrr**12
              ggg(i,ivdw)=12.d0*aaa/rrr**12
            enddo
          
          else
          
            do i=1,mxgrid            
              rrr=dble(i)*dlrpot
              vvv(i,ivdw)=(aaa/rrr**6-bbb)/rrr**6
              ggg(i,ivdw)=6.d0*(2.d0*aaa/rrr**6-bbb)/rrr**6
            enddo
          
          endif
          
        else if(ltpvdw(ivdw).eq.2)then
          
c       lennard-jones potential
      
          eps=prmvdw(ivdw,1)
          sig=prmvdw(ivdw,2)

          if (lspltvdw) then
          
            do i=1,mxgrid
              rrr=dble(i)*dlrpot
              vvv(i,ivdw) =-4.d0*eps*(sig/rrr)**6
              vvvr(i,ivdw)= 4.d0*eps*(sig/rrr)**12
              ggg(i,ivdw) =48.d0*eps*(sig/rrr)**12
            enddo
          
          else
          
            do i=1,mxgrid
              rrr=dble(i)*dlrpot
              vvv(i,ivdw)=4.d0*eps*(sig/rrr)**6*((sig/rrr)**6-1.d0)
              ggg(i,ivdw)=24.d0*eps*(sig/rrr)**6*
     x                                     (2.d0*(sig/rrr)**6-1.d0)
            enddo
          
          endif
          
        else if(ltpvdw(ivdw).eq.3)then

c       n - m potential
      
          eps=prmvdw(ivdw,1)
          ann=max(prmvdw(ivdw,2),prmvdw(ivdw,3))
          amm=min(prmvdw(ivdw,2),prmvdw(ivdw,3))
          rr0=prmvdw(ivdw,4)

          do i=1,mxgrid
            
            rrr=dble(i)*dlrpot
            vvv(i,ivdw)=eps/(ann-amm)*(amm*(rr0/rrr)**ann-
     x        ann*(rr0/rrr)**amm)
            ggg(i,ivdw)=eps*amm*ann/(ann-amm)*((rr0/rrr)**ann-
     x        (rr0/rrr)**amm)
            
          enddo
          
        else if(ltpvdw(ivdw).eq.4)then
          
c       buckingham exp - 6 potential
      
          aaa=prmvdw(ivdw,1)
          rho=prmvdw(ivdw,2)
          ccc=prmvdw(ivdw,3)
          ddd=prmvdw(ivdw,4)

          if (lspltvdw) then
          
            do i=1,mxgrid            
              rrr=dble(i)*dlrpot-ddd
              vvv(i,ivdw) =-ccc/rrr**6
              vvvr(i,ivdw)=aaa*exp(-rrr/rho)
              ggg(i,ivdw) =rrr*aaa*exp(-rrr/rho)/rho            
            enddo
        
          else
          
            do i=1,mxgrid            
              rrr=dble(i)*dlrpot-ddd
              vvv(i,ivdw)=aaa*exp(-rrr/rho)-ccc/rrr**6
              ggg(i,ivdw)=rrr*aaa*exp(-rrr/rho)/rho-6.d0*ccc/rrr**6            
            enddo

          endif
        
        else if(ltpvdw(ivdw).eq.5)then
          
c       born-huggins-meyer exp - 6 - 8 potential
      
          aaa=prmvdw(ivdw,1)
          bbb=prmvdw(ivdw,2)
          ccc=prmvdw(ivdw,3)
          ddd=prmvdw(ivdw,4)
          eee=prmvdw(ivdw,5)

          do i=1,mxgrid
            
            rrr=dble(i)*dlrpot
            vvv(i,ivdw)=aaa*exp(bbb*(ccc-rrr))-ddd/rrr**6-eee/rrr**8
            ggg(i,ivdw)=rrr*aaa*bbb*exp(bbb*(ccc-rrr))-6.d0*ddd/rrr**6
     x        -8.d0*eee/rrr**8
            
          enddo
          
        else if(ltpvdw(ivdw).eq.6) then
          
c       Hydrogen-bond 12 - 10 potential
      
          aaa=prmvdw(ivdw,1)
          bbb=prmvdw(ivdw,2)

          do i=1,mxgrid
            
            rrr=dble(i)*dlrpot
            vvv(i,ivdw)=aaa/rrr**12-bbb/rrr**10
            ggg(i,ivdw)=12.0d0*aaa/rrr**12-10.d0*bbb/rrr**10
            
          enddo
          
        else if(ltpvdw(ivdw).eq.7) then
          
c       shifted and force corrected n - m potential (w. smith)
      
          eps=prmvdw(ivdw,1)
          ann=prmvdw(ivdw,2)
          amm=prmvdw(ivdw,3)
          rr0=prmvdw(ivdw,4)
          rrc=prmvdw(ivdw,5)
          if(rrc.lt.1.d-6)rrc=rcut
 
          if(ann.le.amm) call error(idnode,470)

          gam=rrc/rr0
          if(gam.lt.1.d0) call error(idnode,468)
          bet=gam*((gam**(amm+1.d0)-1.d0)/(gam**(ann+1.d0)-1.d0))
     x      **(1.d0/(ann-amm))
          eps=-eps*(ann-amm)/(amm*(bet**ann)*(1.d0+(ann/gam-ann-1.d0)
     x      /gam**ann)-ann*(bet**amm)*(1.d0+(amm/gam-amm-1.d0)
     x      /gam**amm))

          do i=1,mxgrid
            
            rrr=dble(i)*dlrpot
            if(rrr.gt.rrc)then

              vvv(i,ivdw)=0.d0
              ggg(i,ivdw)=0.d0

            else

              vvv(i,ivdw)=eps/(ann-amm)*(amm*(bet**ann)*((rr0/rrr)**ann-
     x          (1.d0/gam)**ann)-ann*(bet**amm)*((rr0/rrr)**amm-
     x          (1.d0/gam)**amm)+ann*amm*((rrr/(gam*rr0)-1.d0)*
     x          ((bet/gam)**ann-(bet/gam)**amm)))
              ggg(i,ivdw)=eps*amm*ann/(ann-amm)*((bet**ann)*
     x          (rr0/rrr)**ann-(bet**amm)*(rr0/rrr)**amm-rrr/
     x          (gam*rr0)*((bet/gam)**ann-(bet/gam)**amm))

            endif

          enddo
          
        else if(ltpvdw(ivdw).eq.8) then
          
c       morse potential
          
          eps=prmvdw(ivdw,1)
          rr0=prmvdw(ivdw,2)
          sig=prmvdw(ivdw,3)
          
          do i=1,mxgrid
            
            rrr=dble(i)*dlrpot
            vvv(i,ivdw)=eps*((1.d0-exp(-sig*(rrr-rr0)))**2-1.d0)
            ggg(i,ivdw)=-2.d0*rrr*eps*sig*(1.d0-exp(-sig*(rrr-rr0)))*
     x        exp(-sig*(rrr-rr0))
            
          enddo
          
        else if(ltpvdw(ivdw).eq.9) then
          
c       weeks-chandler-anderson potential
          
          eps=prmvdw(ivdw,1)
          sig=prmvdw(ivdw,2)
          rr0=prmvdw(ivdw,3)
          ddd=sig*2.d0**(1.d0/6.d0)
          
          do i=1,mxgrid
            
            rrr=dble(i)*dlrpot-rr0
            if(rrr.gt.ddd)then
              
              vvv(i,ivdw)=0.d0
              ggg(i,ivdw)=0.d0

            else if(rrr.gt.0.d0)then
              
              vvv(i,ivdw)=4.d0*eps*(sig/rrr)**6*
     x          ((sig/rrr)**6-1.d0)+eps
              ggg(i,ivdw)=24.d0*eps*(1.d0+rr0/rrr)*(sig/rrr)**6*
     x          (2.d0*(sig/rrr)**6-1.d0)
            
            endif
              
          enddo
          
        else if(ltpvdw(ivdw).eq.10)then
          
cRS       buckingham exp - 6 potential with short range damping (a la DFT dispersion correction)
c         we use here the damping form according to Stefan Grimme's D3
c         an additional cutoff parameter (in D3 it is a factor times cov. radius) is needed
c         Edisp = ccc/(r^6*t)
c             with the damping factor t = 1+6*(rc/r)^14

      
          aaa=prmvdw(ivdw,1)
          rho=prmvdw(ivdw,2)
          ccc=prmvdw(ivdw,3)
          rc =prmvdw(ivdw,4)

          if (lspltvdw) then

            do i=1,mxgrid
              rrr=dble(i)*dlrpot
              q14=(rc/rrr)**14
              t=1.0d0+6.0d0*q14
              disp  =-ccc/((rrr**6)*t)
c              ddisp =-disp*(6.0d0-84.0d0*(q14/t))
              vvv(i,ivdw) =disp
              vvvr(i,ivdw)=aaa*exp(-rrr/rho)
              ggg(i,ivdw)=rrr*aaa*exp(-rrr/rho)/rho
            enddo

          else
          
            do i=1,mxgrid
              rrr=dble(i)*dlrpot
              q14=(rc/rrr)**14
              t=1.0d0+6.0d0*q14
              disp  =-ccc/((rrr**6)*t)
              ddisp =-disp*(6.0d0-84.0d0*(q14/t))
              vvv(i,ivdw)=aaa*exp(-rrr/rho)+disp
              ggg(i,ivdw)=rrr*aaa*exp(-rrr/rho)/rho+ddisp
            enddo
          
          endif

        else if(ltpvdw(ivdw).eq.11)then
        
cRS       exp+D3 potential
c         we use an exponential A*exp(-r/rho) (with B=1/rho) plus Stefan Grimme's D3 for the dispersion
c         the parameter C6, C8 and Rab are read in. global params are Sr6 and S8 which
c         depend on the functional used to fit the repulsive part

          aaa    =prmvdw(ivdw,1)
          rho    =prmvdw(ivdw,2)
          cc6    =prmvdw(ivdw,3)
          cc8    =prmvdw(ivdw,4)*dispS8
          rc     =prmvdw(ivdw,5)

crs            write(*,*) "D3:  ",aaa,rho,cc6,cc8,rc

          do i=1,mxgrid

            rrr=dble(i)*dlrpot
crs         repulsion
            expterm = aaa*exp(-rrr/rho)
            dexpterm= -expterm/rho
crs         dispersion
            q6 = rc*dispSr6/rrr
            q8 = rc/rrr
            r2 = rrr*rrr
            r6 = r2*r2*r2
            r8 = r6*r2
            q6_14 = q6**14
            q8_16 = q8**16
            t6 = 1.0d0 + 6.0d0*q6_14
            t8 = 1.0d0 + 6.0d0*q8_16
            e6 = cc6/(r6*t6)
            e8 = cc8/(r8*t8)
            disp = e6+e8

crs            write (*,*) rrr, e6, e8, t6, t8
            
            if (lspltvdw) then
              vvvr(i,ivdw)=expterm
              vvv(i,ivdw) =-disp
              ggg(i,ivdw) =rrr*dexpterm
            else
              ddisp=e6*(6.0d0-84.0d0*q6_14/t6)
     &                  +e8*(8.0d0-96.0d0*q8_16/t8)             
              vvv(i,ivdw)=expterm-disp
              ggg(i,ivdw)=rrr*dexpterm+ddisp
            endif

          enddo



        else if(ltpvdw(ivdw).lt.100) then
          
          if(.not.ltable)call error(idnode,150)
          
        endif
        
      enddo
      
cRS apply switching to all potentials if switching is used
      if (.not. (vdwswitch.eq.1.0d0))then

        if(idnode.eq.0) 
     x    write(nrite,*) "Applying switching to vdW potentials" 
     
        do i=1,mxgrid
            rrr=dble(i)*dlrpot
            if (rrr.ge.sw_on) then
                rr2 = rrr*rrr
                rr3 = rr2*rrr
                rr4 = rr2*rr2
                rr5 = rr3*rr2
                sw_value = c5*rr5+c4*rr4+c3*rr3+c2*rr2+c1*rrr+c0
                sw_dval  = 5.0d0*c5*rr4+4.0d0*c4*rr3+3.0d0*c3*rr2+
     x                     2.0d0*c2*rrr+c1
                do ivdw=1,ntpvdw
                    vvv_temp = vvv(i,ivdw)
                    ggg_temp = ggg(i,ivdw)
                    vvv(i,ivdw) = vvv_temp*sw_value
                    ggg(i,ivdw) = ggg_temp*sw_value+vvv_temp*sw_dval
                    if (lspltvdw) then
                      vvv_temp = vvvr(i,ivdw)
                      vvvr(i,ivdw) = vvv_temp*sw_value
                    endif
                enddo    
            endif
        enddo
      
      endif
cRS end of switching block

cRS short range damping
cRS  NOTE: ggg is not the gradient but the force times r (so to get the actual gradient we need to divide by r)

      if (lsrdamp) then
        if (.not.lvdwspline) then
          write (nrite, *)  
     x      "Do not use short range vdw damping without vdwspline!"
          call exitcomms(0)
        endif
        
        
        if (lspltvdw) then
c         This is a bit painful but with lsplitvdw the repulsive part is in vvvr (ggg is the gradient of the rep only)
          do ivdw=1,ntpvdw
            do i=mxgrid,1,-1
              if (vvvr(i,ivdw).gt.srdamp_ener) then
                icut = i
                exit
              endif
            enddo
c           compute parabola E(r) = a-br^2
            r_srcut = dble(icut)*dlrpot
c           ggg is the force times r and not the gradient (despite the name)
            b = -ggg(icut,ivdw)/(-2.0d0*r_srcut*r_srcut)
            a = vvvr(icut,ivdw)+b*r_srcut*r_srcut
c           now fill low r values with parabola ..
            do i=icut-1,1,-1
              rrr = dble(i)*dlrpot
              vvvr(i,ivdw) = a - b*rrr*rrr
              ggg(i,ivdw) = 2.0d0*b*rrr*rrr
            enddo
c           NEW: if we truncate the repulsion we need to truncate the dispersion as well (for LJ-type potentials)
c                because we get largely attractive energies at close distances.
c                thus we truncate in the same way at the same point icut.
c                since we do not have a gradient for vvv (ggg is the gradient of vvvr only) we take an FD gradient 
c                here for simplicity
            ggg_disp = (vvv(icut+1,ivdw)-vvv(icut,ivdw))/dlrpot
            b = ggg_disp/(-2.0d0*r_srcut)
            a = vvv(icut,ivdw)+b*r_srcut*r_srcut
c           now fill low r values with parabola ..
            do i=icut-1,1,-1
              rrr = dble(i)*dlrpot
              vvv(i,ivdw) = a - b*rrr*rrr
            enddo
          enddo        
        else
          do ivdw=1,ntpvdw
            do i=mxgrid,1,-1
              if (vvv(i,ivdw).gt.srdamp_ener) then
                icut = i
                exit
              endif
            enddo
c           compute parabola E(r) = a-br^2
            r_srcut = dble(icut)*dlrpot
c           ggg is the force times r and not the gradient (despite the name)
            b = -ggg(icut,ivdw)/(-2.0d0*r_srcut*r_srcut)
            a = vvv(icut,ivdw)+b*r_srcut*r_srcut
c           now fill low r values with parabola ..
            do i=icut-1,1,-1
              rrr = dble(i)*dlrpot
              vvv(i,ivdw) = a - b*rrr*rrr
              ggg(i,ivdw) = 2.0d0*b*rrr*rrr
            enddo
          enddo
        endif
      endif

cRS end of short range damping

cRS use cubic spline interpolation instead
c
c          if spline interp is used we store the 2nd derivative
c          in vvv2
c          => in principle we should skip computing ggg at all
c             but this is done only once during setup
c             as ususal there is always room for improvement

      if (lvdwspline) then

        if(idnode.eq.0) 
     x    write(nrite,*) "Generating second derivs for spline interp" 
c      
        do ivdw=1,ntpvdw
          vvv2(1,ivdw)= 0.0
          u(1)  = 0.0
          do i=2,mxgrid-1
            p = 0.5d0*vvv2(i-1,ivdw)+2.0
            vvv2(i,ivdw) = -0.5d0/p
            u(i) = ((vvv(i+1,ivdw)-vvv(i,ivdw))-
     x              (vvv(i,ivdw)-vvv(i-1,ivdw)))/dlrpot
            u(i) = (6.0d0*u(i)/(2.0d0*dlrpot)-0.5d0*u(i-1))/p
          enddo
          vvv2(mxgrid,ivdw) = 0.0
          do i = mxgrid-1, 1, -1
            vvv2(i,ivdw) = vvv2(i,ivdw)*vvv2(i+1,ivdw)+u(i)
          enddo
        enddo

        if (lspltvdw) then
c         do the spline table also for the repulsive part if we splitting
          do ivdw=1,ntpvdw
            vvvr2(1,ivdw)= 0.0
            u(1)  = 0.0
            do i=2,mxgrid-1
              p = 0.5d0*vvvr2(i-1,ivdw)+2.0
              vvvr2(i,ivdw) = -0.5d0/p
              u(i) = ((vvvr(i+1,ivdw)-vvvr(i,ivdw))-
     x              (vvvr(i,ivdw)-vvvr(i-1,ivdw)))/dlrpot
              u(i) = (6.0d0*u(i)/(2.0d0*dlrpot)-0.5d0*u(i-1))/p
            enddo
            vvvr2(mxgrid,ivdw) = 0.0
            do i = mxgrid-1, 1, -1
              vvvr2(i,ivdw) = vvvr2(i,ivdw)*vvvr2(i+1,ivdw)+u(i)
            enddo
          enddo        
        endif

      endif 

      deallocate(u, stat=fail)

cRS end of spline block

CVAM
CVAM      call VTEND(140, ierr)
CVAM
      return
      end subroutine forgen

      subroutine fortab
     x  (idnode,ntpvdw,ntpatm,dlrpot,rcut,engunit)

c***********************************************************************
c     
c     dl_poly subroutine for reading potential energy and 
c     force arrays for van der waals forces only
c     
c     copyright - daresbury laboratory 1994
c     author    - w. smith march 1994
c     
c     wl
c     2009/01/13 11:22:05
c     1.7
c     Exp
c     
c***********************************************************************
      
      implicit none

      logical safe
      character*8 atom1,atom2
      integer idnode,ntpvdw,ntpatm,idum,ngrid
      integer ivdw,katom1,katom2,jtpatm,l,i,k,keyvdw
      real(8) dlrpot,rcut,engunit,delpot,cutpot,rdr,rrr,ppp
      real(8) vk0,vk1,vk2,t1,t2
CVAM
CVAM      call VTBEGIN(142, ierr)
CVAM

      if(idnode.eq.0)open (ntable,file='TABLE')

c     skip header record
      
      call getrec(safe,idnode,ntable)
      if(.not.safe)call abort_table_read(idnode,ntable)

c     read mesh resolution
      
      call getrec(safe,idnode,ntable)
      if(.not.safe)call abort_table_read(idnode,ntable)
      delpot=dblstr(record,lenrec,idum)
      cutpot=dblstr(record,lenrec,idum)
      ngrid=intstr(record,lenrec,idum)
      dlrpot=rcut/dble(mxgrid-4)
      if(abs(delpot-dlrpot)/dlrpot.le.1.d-4)delpot=dlrpot
      if((delpot.gt.dlrpot).or.(ngrid-4.ne.nint(cutpot/delpot)))then

        if(idnode.eq.0) then
          write(nrite,"('expected radial increment : ',1p,e15.7,/,
     x                '   TABLE radial increment : ',1p,e15.7,/,/,    
     x                'expected number of grid points : ',0p,i10,/,
     x                'grid points in TABLE           : ',i10)")
     x      dlrpot, delpot, mxgrid, ngrid
        endif
        
        call error(idnode,22)

      endif

      if(cutpot.lt.rcut) call error(idnode,504)
      if(idnode.eq.0) then
        if(abs(1.d0-(delpot/dlrpot)).gt.1d-7) write(nrite,
     x    "(/,' TABLE arrays resized for mxgrid=',i10)") mxgrid

      endif

c     read potential arrays for all pairs
      
      do ivdw=1,ntpvdw

c     read potential arrays if potential not already defined
        
        if(ltpvdw(ivdw).eq.0)then
          
c     read pair potential labels and long range corrections
          
          call getrec(safe,idnode,ntable)
          if(.not.safe)call abort_table_read(idnode,ntable)

          call getword(atom1,record,8,lenrec)
          call getword(atom2,record,8,lenrec)
          prmvdw(ivdw,1)=dblstr(record,lenrec,idum)
          prmvdw(ivdw,2)=dblstr(record,lenrec,idum)
          
          katom1=0
          katom2=0
          
          do jtpatm=1,ntpatm
            
            if(atom1.eq.unqatm(jtpatm))katom1=jtpatm
            if(atom2.eq.unqatm(jtpatm))katom2=jtpatm
            
          enddo
          
          if(katom1.eq.0.or.katom2.eq.0) then
            if(idnode.eq.0) 
     x        write(nrite,'(a)') '****',atom1,'***',atom2,'****'
            call  error(idnode,81)
          endif
          
          keyvdw=loc2(katom1,katom2)
          
          if(lstvdw(keyvdw).ne.ivdw) call error(idnode,23)
          
c     read potential arrays
          
          if(idnode.eq.0)then

            if(mxbuff.lt.ngrid)  then
              
              write(nrite,*) 'mxbuff must be >=',ngrid,' in fortab'
              call error(idnode,48)
              
            endif

            read(ntable,'(4e15.8)',end=100)(buffer(i),i=1,ngrid)

c     reconstruct arrays using 3pt interpolation

            rdr=1.d0/delpot
            vvv(1,ivdw)=1.d0
            ggg(1,ivdw)=0.d0
            do i=2,mxgrid
              rrr=dble(i)*dlrpot
              l=int(rrr*rdr)
              ppp=rrr*rdr-dble(l)
              vk0=buffer(l)
              vk1=buffer(l+1)
              vk2=buffer(l+2)
            
              t1=vk0+(vk1-vk0)*ppp
              t2=vk1+(vk2-vk1)*(ppp-1.0d0)
              vvv(i,ivdw)=t1+(t2-t1)*ppp*0.5d0

            enddo

            read(ntable,'(4e15.8)',end=100)(buffer(i),i=1,ngrid)

c     reconstruct ggg arrays using 3pt interpolation

            do i=2,mxgrid

              rrr=dble(i)*dlrpot
              l=int(rrr*rdr)
              ppp=rrr*rdr-dble(l)
              vk0=buffer(l)
              vk1=buffer(l+1)
              vk2=buffer(l+2)
            
              t1=vk0+(vk1-vk0)*ppp
              t2=vk1+(vk2-vk1)*(ppp-1.0d0)
            
              ggg(i,ivdw)=t1+(t2-t1)*ppp*0.5d0

            enddo

            call gdsum(vvv(1,ivdw),mxgrid,buffer)
            call gdsum(ggg(1,ivdw),mxgrid,buffer)

          else
            
            if(mxbuff.lt.mxgrid) call error(idnode,48)

            do i=1,mxgrid

              vvv(i,ivdw)=0.d0
              ggg(i,ivdw)=0.d0

            enddo

            call gdsum(vvv(1,ivdw),mxgrid,buffer)
            call gdsum(ggg(1,ivdw),mxgrid,buffer)

          endif

        endif
        
      enddo

c     convert to internal units
      
      do k=1,ntpvdw
        
        if(ltpvdw(k).eq.0)then

          do i=1,mxgrid
            
            vvv(i,k)=vvv(i,k)*engunit
            ggg(i,k)=ggg(i,k)*engunit
            
          enddo
          
        endif
        
      enddo
      
      if(idnode.eq.0)close (ntable)
      
      if(idnode.eq.0)write(nrite,'(/,/,1x,a)')
     x  'potential tables read from TABLE file'
CVAM
CVAM      call VTEND(142, ierr)
CVAM
      return
      
c     end of file error exit
      
  100 call abort_table_read(idnode,ntable)

      end subroutine fortab

      subroutine abort_table_read(idnode,ntable)

c***********************************************************************
c     
c     dl_poly error exit subroutine for reading TABLE file
c     
c     copyright - daresbury laboratory 
c     author    - w. smith   sept 2005
c     

      implicit none
      integer idnode,ntable

      if(idnode.eq.0)close (ntable)
      
      call error(idnode,24)
      
      end subroutine abort_table_read

      subroutine srfrce
     x  (lsolva,lfree,lghost,iatm,ik,engsrp,virsrp,rcut,dlrpot)

c***********************************************************************
c     
c     dl_poly subroutine for calculating short range force and
c     potential energy terms using verlet neighbour list
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith       march 1992
c     adapted   - p.-a. cazade oct 2007: solvation, free energy etc
c     
c     version 3
c     author    - t. forester    june  1993
c     stress tensor added t.forester may 1994
c     
c     wl
c     2009/01/13 11:22:05
c     1.7
c     Exp
c     
c***********************************************************************
      
      implicit none

      logical lsolva,lfree,lghost,lselect,lskip
      integer iatm,ik,m,jatm,k,l,kkk
      real(8) engsrp,virsrp,rcut,dlrpot
      real(8) fi,rcsq,rdr,strs1,strs2,strs3,strs5,strs6,strs9,ai,aj
      real(8) ab,rrr,rsq,ppp,t1,t2,vk0,vk1,vk2,gk0,gk1,gk2,gamma,fx
      real(8) fy,fz,omega
cRS      
      real(8) omega_rep, gamma_rep   
      real(8) faca, facb, dfaca, dfacb

cRS
      real(8) der2, der1
      real(8) a, b, v20, v21
cRS switchable molecules
      integer moli,molj
      real(8) lambi,lambj,lambinti
      real(8) lambrepi, lambrepj

      dimension fi(3)

CDIR$ CACHE_ALIGN fi
CVAM
CVAM      call VTBEGIN(83, ierr)
CVAM
      lskip=(lfree.or.lghost)

c     set cutoff condition for pair forces

      rcsq=rcut**2

c     interpolation spacing
      
      rdr=1.d0/dlrpot

c     initialise stress tensor accumulators

      strs1=0.d0
      strs2=0.d0
      strs3=0.d0
      strs5=0.d0
      strs6=0.d0
      strs9=0.d0

c     initialise potential energy and virial
      
      engsrp=0.d0
      virsrp=0.d0

c     store forces for iatm 
      
      ai=dble(ltype(iatm))
      fi(1)=fxx(iatm)
      fi(2)=fyy(iatm)
      fi(3)=fzz(iatm)
cRS
      if (lmolecules) then
        moli = mol_which(iatm)
        lambi = mol_lamb_vdw(moli)
        lambinti = mol_lambint_vdw(moli)
        lambrepi = mol_lamb_vdwr(moli)
      endif

c     start of primary loop for forces evaluation
      
      do m=1,ik

c     atomic and potential function indices
        
        jatm=ilist(m)
        
        if(lskip)then
          if(atm_fre(iatm)*atm_fre(jatm).eq.2)cycle
        endif
        
        aj=dble(ltype(jatm))
        
        if(ai.gt.aj) then
          ab=ai*(ai-1.d0)*0.5d0+aj+0.5d0
        else
          ab=aj*(aj-1.d0)*0.5d0+ai+0.5d0
        endif
        
        k=lstvdw(int(ab))
        
        if((ltpvdw(k).lt.100).and.(abs(vvv(1,k)).gt.1.d-10))then

c     apply truncation of potential
          
          rsq=rsqdf(m)
          
          if(rcsq.gt.rsq)then
            
            rrr=sqrt(rsq)               
            l=int(rrr*rdr)
            ppp=rrr*rdr-dble(l)

cRS
            if (.not.lvdwspline) then

c     calculate interaction energy using 3-point interpolation
            
            vk0=vvv(l,k)
            vk1=vvv(l+1,k)
            vk2=vvv(l+2,k)
            t1=vk0+(vk1-vk0)*ppp
            t2=vk1+(vk2-vk1)*(ppp-1.0d0)
            omega=t1+(t2-t1)*ppp*0.5d0
            

c     calculate forces using 3-point interpolation
            
            gk0=ggg(l,k)
            gk1=ggg(l+1,k)
            gk2=ggg(l+2,k)
            t1=gk0+(gk1-gk0)*ppp
            t2=gk1+(gk2-gk1)*(ppp-1.0d0)
            gamma=(t1+(t2-t1)*ppp*0.5d0)/rsq

            else

cRS  THIS IS THE NEW DEFAULT (because the old original dl_poly code is simply bogus:
c                             the force is computed form a seperate interpolation and is inconsistent with the energy
c                             => no energy conserving dynamics
c                             now the force is computed as the derivative of the spline interpolation)
              vk0=vvv(l,k)
              vk1=vvv(l+1,k)
              v20 = vvv2(l,k)
              v21 = vvv2(l+1,k)
              a = (((l+1)*dlrpot)-rrr)/dlrpot
              b = 1.0d0-a
              faca = (a*a*a)-a
              facb = (b*b*b)-b
              dfaca = 3.0d0*a*a-1.0d0
              dfacb = 3.0d0*b*b-1.0d0
              omega = a*vk0+b*vk1+
     x        (faca*v20+facb*v21)*dlrpot*dlrpot/6.0d0
              gamma = (vk1-vk0)/dlrpot-
     x        (dfaca*v20-dfacb*v21)*dlrpot/6.0d0
              gamma = -gamma/rrr
              if (lspltvdw) then
                vk0=vvvr(l,k)
                vk1=vvvr(l+1,k)
                v20 = vvvr2(l,k)
                v21 = vvvr2(l+1,k)
                omega_rep = a*vk0+b*vk1+
     x             (faca*v20+facb*v21)*dlrpot*dlrpot/6.0d0
                gamma_rep = (vk1-vk0)/dlrpot-
     x             (dfaca*v20-dfacb*v21)*dlrpot/6.0d0
                gamma_rep = -gamma_rep/rrr
              endif
            endif

cRS DEBUG ---> not working properly and i do not see why!!
            if (lvdwintder) then
             write (*,*) "DO NOT USE VDW FORCE INTERPOALTION!"
             call exit(0)
c             der1 = 2.0d0*vk1-0.5d0*vk2-1.5d0*vk0
c             der2 = vk0+0.5d0*vk2-1.5d0*vk1
c             write (*,*) "VDW FORCE  ", rrr, omega,
c     x         vk0+der1*ppp+der2*ppp*ppp
c             write (*,*) der1, der2, ppp
            endif
            
c     set selection control
            
            lselect=.true.
            
c     set double index
            
            if(lsolva)kkk=loc2(atmolt(iatm),atmolt(jatm))
            
            if(lghost)then
              
c     selected excitation option
              
              if((atm_fre(iatm).ne.1).and.(atm_fre(jatm).ne.1))then
                
c     reset selection control
                
                lselect=(atm_fre(iatm)+atm_fre(jatm).eq.0)
                
                if(lsolva)vdw_exc(kkk)=vdw_exc(kkk)+omega
                
              endif
              
            elseif(lfree)then
              
c     selected free energy option
              
              if((atm_fre(iatm).eq.1).or.(atm_fre(jatm).eq.1))then
                
c     set hamiltonian mixing parameter
                
                vdw_fre=vdw_fre-omega
                omega=lambda1*omega
                gamma=lambda1*gamma
                
              elseif((atm_fre(iatm).eq.2).or.(atm_fre(jatm).eq.2))then
                
c     set hamiltonian mixing parameter
                
                vdw_fre=vdw_fre+omega
                omega=lambda2*omega
                gamma=lambda2*gamma
                
              endif
cRS
            elseif(lmolecules) then

                molj=mol_which(jatm)
                if (molj.ne.moli) then

c     two atoms belong to different molecules -> scale by lambda1*lambda2
c       the force on lambda1 = lambda2*energy and vice versa
 
                   lambj    =mol_lamb_vdw(molj)
                   mol_dlam_vdw(moli)=mol_dlam_vdw(moli) + lambj*omega
                   mol_dlam_vdw(molj)=mol_dlam_vdw(molj) + lambi*omega
                   omega = omega*lambi*lambj
                   gamma = gamma*lambi*lambj
                   if (lspltvdw) then
c                    in case of split lambda then the dispersive part is in omega/gamma and
c                    the repulsive in omega_rep and gamma_rep. in this case we need to scale the latter 
c                    seperate but then add to omega/gamma
                     lambrepj =mol_lamb_vdwr(molj)
                     mol_dlam_vdwr(moli)=mol_dlam_vdwr(moli)
     x                                             + lambrepj*omega_rep
                     mol_dlam_vdwr(molj)=mol_dlam_vdwr(molj)
     x                                             + lambrepi*omega_rep
                     omega_rep = omega_rep*lambrepi*lambrepj
                     gamma_rep = gamma_rep*lambrepi*lambrepj
                     omega = omega + omega_rep
                     gamma = gamma + gamma_rep
#ifdef DEBUG
                if ((rrr.lt.1.0d0).and.((lambi*lambj).gt.0.0d0)) then
                  if (iatm.gt.jatm) then
                    write (*,*) "VDW  ",jatm,iatm,rrr,omega/418.4d0,
     x                      omega_rep/418.4d0, lambj, lambi,
     x                      lambrepj, lambrepi
                  else
                    write (*,*) "VDW  ",iatm,jatm,rrr,omega/418.4d0,
     x                      omega_rep/418.4d0, lambi, lambj,
     x                      lambrepi, lambrepj
                  end if
                endif
#endif
                   endif

                else

c      two atoms belong to the same molecule -> scale by lambint
c         currently only a scaling is applied but no force computed
c         => only useful to switch off interactions but not for TI

                   if (lspltvdw) then
                     omega = omega + omega_rep
                     gamma = gamma + gamma_rep
                   endif
                   omega = omega*lambinti
                   gamma = gamma*lambinti

                endif

              
            endif
            
            if(lselect)then
              
c     calculate potential and virial
              
              engsrp=engsrp+omega
              virsrp=virsrp-gamma*rsq
              
              if(lsolva)vdw_sol(kkk)=vdw_sol(kkk)+omega
              
c     calculate forces
              
              fx=gamma*xdf(m)
              fy=gamma*ydf(m)
              fz=gamma*zdf(m)
              
              fi(1)=fi(1)+fx
              fi(2)=fi(2)+fy
              fi(3)=fi(3)+fz
              
              fxx(jatm)=fxx(jatm)-fx
              fyy(jatm)=fyy(jatm)-fy
              fzz(jatm)=fzz(jatm)-fz
              
c     calculate stress tensor
              
              strs1=strs1+xdf(m)*fx
              strs2=strs2+xdf(m)*fy
              strs3=strs3+xdf(m)*fz
              strs5=strs5+ydf(m)*fy
              strs6=strs6+ydf(m)*fz
              strs9=strs9+zdf(m)*fz
              
            endif
            
          endif
          
        endif
        
      enddo

c     load temps back to fxx(iatm) etc
      
      fxx(iatm)=fi(1)
      fyy(iatm)=fi(2)
      fzz(iatm)=fi(3)

c     complete stress tensor
        
      stress(1)=stress(1)+strs1
      stress(2)=stress(2)+strs2
      stress(3)=stress(3)+strs3
      stress(4)=stress(4)+strs2
      stress(5)=stress(5)+strs5
      stress(6)=stress(6)+strs6
      stress(7)=stress(7)+strs3
      stress(8)=stress(8)+strs6
      stress(9)=stress(9)+strs9
CVAM
CVAM      call VTEND(83, ierr)
CVAM
      return
      end subroutine srfrce

      subroutine lrcorrect
     x  (lsolva,lfree,lghost,idnode,imcon,keyfce,natms,
     x  ntpatm,ntpvdw,elrc,engunit,virlrc,rcut,volm)
      
c*************************************************************************
c     
c     DL_POLY subroutine to evaluate long-range corrections to
c     pressure and energy in a periodic system.
c     
c     copyright daresbury laboratory 1993
c     author    - t. forester may 1993
c     adapted   - p.-a. cazade oct 2007: solvation, free energy etc
c     
c     wl
c     2009/01/13 11:22:05
c     1.7
c     Exp
c     
c***************************************************************************
      
      implicit none

      integer, parameter :: nnn=10
      logical lsolva,lfree,lghost
      integer idnode,imcon,keyfce,natms,ntpatm,i,ka,ntpvdw
      real(8) natyp,nbtyp,nctyp,ndtyp,nafrz,nbfrz,ncfrz,ndfrz
      integer ivdw,j,k,it,jt,kt,fail
      real(8) elrc,engunit,virlrc,rcut,volm,twopi,eadd,padd
      real(8) denprd,aaa,bbb,ccc,ddd,eee,eps,sig,rr0,ann,amm
      real(8) denprd1,denprd2,denprd3,denprdf
      integer, allocatable :: numtyp_sol0(:,:),numfrz_sol0(:,:)
      integer, allocatable :: numtyp_sol1(:,:),numfrz_sol1(:,:)
      integer, allocatable :: numtyp_sol2(:,:),numfrz_sol2(:,:)
      integer, allocatable :: numtyp_fre(:,:),numfrz_fre(:,:)
      real(8), allocatable :: elrc_sol0(:),elrc_exc0(:)
      
      dimension fail(nnn)
CVAM
CVAM      call VTBEGIN(148, ierr)
CVAM
      twopi=2.0d0*pi
      
c     allocate working arrays
      
      do i=1,nnn
        fail(i)=0
      enddo
      
      if(lfree.or.lghost)then
        
        allocate (numtyp_fre(mxatyp,0:2),stat=fail(1))
        allocate (numfrz_fre(mxatyp,0:2),stat=fail(2))
        allocate (elrc_exc0(mxtmls_exc2),stat=fail(3))
        
      endif
      
      if(lsolva)then
        
        allocate (elrc_sol0(mxtmls_sol2),stat=fail(4))
        allocate (numtyp_sol0(mxatyp,mxtmls),stat=fail(5))
        allocate (numfrz_sol0(mxatyp,mxtmls),stat=fail(6))
        
        if(lghost)then
          
          allocate (numtyp_sol1(mxatyp,mxtmls),stat=fail(7))
          allocate (numfrz_sol1(mxatyp,mxtmls),stat=fail(8))
          allocate (numtyp_sol2(mxatyp,mxtmls),stat=fail(9))
          allocate (numfrz_sol2(mxatyp,mxtmls),stat=fail(10))
          
        endif
        
      endif
      
      do i=1,nnn
        if(fail(i).ne.0)call error(idnode,1015)
      enddo
      
c     initalise counter arrays
      
      do i=1,ntpatm
        
        numtyp(i)=0
        numfrz(i)=0
        
      enddo
      
      if(lfree.or.lghost)then
        
        numtyp_fre(:,:)=0
        numfrz_fre(:,:)=0
        
      endif
      
      if(lsolva)then
        
        numtyp_sol0(:,:)=0
        numfrz_sol0(:,:)=0
        
        if(lghost)then
          
          numtyp_sol1(:,:)=0
          numfrz_sol1(:,:)=0
          numtyp_sol2(:,:)=0
          numfrz_sol2(:,:)=0
          
        endif
        
      endif
      
c     evaluate number density in system
      
      do i=1,natms
        
        ka=ltype(i)
        numtyp(ka)=numtyp(ka)+1
        if(lstfrz(i).ne.0)numfrz(ka)=numfrz(ka)+1
        
      enddo
      
      if(lfree.or.lghost)then
         
        do i=1,natms
          
          ka=ltype(i)
          numtyp_fre(ka,atm_fre(i))=numtyp_fre(ka,atm_fre(i))+1
          if(lstfrz(i).ne.0)
     x      numfrz_fre(ka,atm_fre(i))=numfrz_fre(ka,atm_fre(i))+1          
          
        enddo
        
      endif
      
      if(lsolva)then
        
        if(lghost)then
          
          do i=1,natms
            
            ka=ltype(i)
            
            if(atm_fre(i).eq.0)then
              
              numtyp_sol0(ka,atmolt(i))=numtyp_sol0(ka,atmolt(i))+1
              if(lstfrz(i).ne.0)
     x          numfrz_sol0(ka,atmolt(i))=numfrz_sol0(ka,atmolt(i))+1
              
            elseif(atm_fre(i).eq.1)then
              
              numtyp_sol1(ka,atmolt(i))=numtyp_sol1(ka,atmolt(i))+1
              if(lstfrz(i).ne.0)
     x          numfrz_sol1(ka,atmolt(i))=numfrz_sol1(ka,atmolt(i))+1
              
            elseif(atm_fre(i).eq.2)then
              
              numtyp_sol2(ka,atmolt(i))=numtyp_sol2(ka,atmolt(i))+1
              if(lstfrz(i).ne.0)
     x          numfrz_sol2(ka,atmolt(i))=numfrz_sol2(ka,atmolt(i))+1
              
            endif
            
          enddo
          
        else
          
          do i=1,natms
            
            ka=ltype(i)
            numtyp_sol0(ka,atmolt(i))=numtyp_sol0(ka,atmolt(i))+1
            if(lstfrz(i).ne.0)
     x        numfrz_sol0(ka,atmolt(i))=numfrz_sol0(ka,atmolt(i))+1
            
          enddo
          
        endif
        
      endif
      
c     number densities
      
      do i=1,ntpatm
        dens(i)=dble(numtyp(i))/volm
      enddo
      
c     long range corrections to energy and pressure
      
      elrc=0.d0
      elrc2=0.d0
      virlrc=0.d0
      virlrc2=0.d0
      denprdf=0.d0
      elrc_fre=0.d0
      
      if(imcon.ne.0.and.imcon.ne.6.and.ntpvdw.gt.0) then 
         
        if(mod(keyfce,2).eq.1) then
          
          ivdw=0
          
          do i=1,ntpatm
            
            do j=1,i
               
              eadd=0.d0
              padd=0.d0
              
              ivdw=ivdw+1
              k=lstvdw(ivdw)
              
              if(ltpvdw(k).eq.0) then
                
c     tabulated potential
                
                eadd=prmvdw(k,1)
                padd=-prmvdw(k,2)
                
              else if(ltpvdw(k).eq.1) then
                
c     12-6 potential
                
                aaa=prmvdw(k,1)
                bbb=prmvdw(k,2)
                
                eadd=aaa/(9.d0*rcut**9)-bbb/(3.d0*rcut**3)
                padd=12.d0*aaa/(9.d0*rcut**9)-6.d0*bbb/(3.d0*rcut**3)
                
              else if(ltpvdw(k).eq.2) then
                
c     Lennard Jones potential
      
                eps=prmvdw(k,1)
                sig=prmvdw(k,2)
                
                eadd=4.d0*eps*(sig**12/(9.d0*rcut**9)-
     x            sig**6/(3.d0*rcut**3))
                padd=4.d0*eps*(12.d0*sig**12/(9.d0*rcut**9)-
     x            2.d0*sig**6/(rcut**3))
                
              else if(ltpvdw(k).eq.3) then
                
c     n - m potential
                
                eps=prmvdw(k,1)
                ann=prmvdw(k,2)
                amm=prmvdw(k,3)
                rr0=prmvdw(k,4)
                
                eadd=eps/(ann-amm)*(amm*rr0**ann/((ann-3.d0)*
     x            rcut**(ann-3.d0))-ann*rr0**amm/((amm-3.0d0)*
     x            rcut**(amm-3.d0)))
                padd=eps/(ann-amm)*ann*amm*(rr0**ann/((ann-3.d0)*
     x            rcut**(ann-3.d0))-rr0**amm/((amm-3.0d0)*
     x            rcut**(amm-3.d0)))
                
              else if((ltpvdw(k).eq.4).or.(ltpvdw(k).eq.10)) then
                
c     buckingham exp - 6 potential
cRS or damped exp -6 potential (same long range correction as buckingham)
                
                ccc=prmvdw(k,3)
                
                eadd=-ccc/(3.d0*rcut**3)
                padd=-2.d0*ccc/(rcut**3)
                
              else if(ltpvdw(k).eq.5) then
                
c     born huggins meyer exp -6 - 8  potential
                
                ddd=prmvdw(k,4)
                eee=prmvdw(k,5)
                
                eadd=-ddd/(3.d0*rcut**3)-eee/(5.d0*rcut**5)
                padd=-2.d0*ddd/(rcut**3)-8.d0*eee/(5.d0*rcut**5)
                
              else if(ltpvdw(k).eq.6) then
                
c     hydrogen bond  12 - 10 potential
                
                aaa=prmvdw(k,1)
                bbb=prmvdw(k,2)
                
                eadd=aaa/(9.d0*rcut**9)-bbb/(7.d0*rcut**7)
                padd=12.d0*aaa/(9.d0*rcut**9)-1.d1*bbb/(7.d0*rcut**7)
                
              endif
              
              if(i.ne.j) then
                
                eadd=eadd*2.d0
                padd=padd*2.d0
                
              endif
              
              if(.not.(lfree.or.lghost))then
                
                denprd=twopi*(dble(numtyp(i))*dble(numtyp(j))-
     x            dble(numfrz(i))*dble(numfrz(j)))/volm**2
                
              else
                 
                nafrz=dble(numfrz_fre(i,0)+numfrz_fre(i,1))
                natyp=dble(numtyp_fre(i,0)+numtyp_fre(i,1))
                nbfrz=dble(numfrz_fre(j,0)+numfrz_fre(j,1))
                nbtyp=dble(numtyp_fre(j,0)+numtyp_fre(j,1))
                ncfrz=dble(numfrz_fre(i,0)+numfrz_fre(i,2))
                nctyp=dble(numtyp_fre(i,0)+numtyp_fre(i,2))
                ndfrz=dble(numfrz_fre(j,0)+numfrz_fre(j,2))
                ndtyp=dble(numtyp_fre(j,0)+numtyp_fre(j,2))
                
                if(lghost)then
                  
                  denprd=twopi*(natyp*nbtyp-nafrz*nbfrz)/volm**2
                  denprd3=twopi*(nctyp*ndtyp-ncfrz*ndfrz)/volm**2
                  
                elseif(lfree)then
                  
                  denprd1=twopi*(natyp*nbtyp-nafrz*nbfrz)/volm**2
                  denprd2=twopi*(nctyp*ndtyp-ncfrz*ndfrz)/volm**2
                  denprd=lambda1*denprd1+lambda2*denprd2
                  denprd3=lambda2*denprd1+lambda1*denprd2
                  denprdf=denprd2-denprd1
                  
                endif
                
              endif
              
              elrc=elrc+volm*denprd*eadd
              virlrc=virlrc-denprd*padd*volm
              
              if(lfree.or.lghost)then
                
                elrc2=elrc2+volm*denprd3*eadd
                virlrc2=virlrc2-denprd3*padd*volm
                if(lfree)elrc_fre=elrc_fre+volm*denprdf*eadd
                
              endif
              
              if(lsolva)then
                
                elrc_sol0(:)=0.d0
                if(lghost)elrc_exc0(:)=0.d0
                
                do it=1,mxtmls
                  
                  do jt=1,mxtmls
                    
                    kt=loc2(it,jt)
                    
                    if(lghost)then
                       
                      natyp=dble(numtyp_sol0(i,it)+numtyp_sol1(i,it))
                      nbtyp=dble(numtyp_sol0(j,jt)+numtyp_sol1(j,jt))
                      nafrz=dble(numfrz_sol0(i,it)+numfrz_sol1(i,it))
                      nbfrz=dble(numfrz_sol0(j,jt)+numfrz_sol1(j,jt))
                      
                      elrc_sol0(kt)=elrc_sol0(kt)+twopi*(natyp*
     x                nbtyp-nafrz*nbfrz)/volm**2
                      
                      nctyp=dble(numtyp_sol0(i,it)+numtyp_sol2(i,it))
                      ndtyp=dble(numtyp_sol0(j,jt)+numtyp_sol2(j,jt))
                      ncfrz=dble(numfrz_sol0(i,it)+numfrz_sol2(i,it))
                      ndfrz=dble(numfrz_sol0(j,jt)+numfrz_sol2(j,jt))
                      
                      elrc_exc0(kt)=elrc_exc0(kt)+twopi*(nctyp*
     x                ndtyp-ncfrz*ndfrz)/volm**2
                      
                    else
                      
                      natyp=dble(numtyp_sol0(i,it))
                      nbtyp=dble(numtyp_sol0(j,jt))
                      nafrz=dble(numfrz_sol0(i,it))
                      nbfrz=dble(numfrz_sol0(j,jt))
                      
                      elrc_sol0(kt)=elrc_sol0(kt)+twopi*(natyp*
     x                nbtyp-nafrz*nbfrz)/volm**2             
                      
                    endif
                    
                  enddo
                  
                enddo
                
                if(lghost)then
                   
                  elrc_sol(:)=elrc_sol(:)+volm*eadd*elrc_sol0(:)
                  elrc_exc(:)=elrc_exc(:)+volm*eadd*elrc_exc0(:)
                  
                else
                  
                  elrc_sol(:)=elrc_sol(:)+volm*eadd*elrc_sol0(:)
                  
                endif
                
              endif
              
            enddo
            
          enddo
          
          if(lfree.or.lghost)then
             
            elrc_sav=elrc
            elrc2_sav=elrc2
            virlrc_sav=virlrc
            virlrc2_sav=virlrc2
            elrc_fre_sav=elrc_fre
            
          endif
          
          volm_sav=volm
          
          if(lghost)then
             
            elrc_sol_sav(:)=elrc_sol(:)
            elrc_exc_sav(:)=elrc_exc(:)
            
          elseif(lsolva)then
            
            elrc_sol_sav(:)=elrc_sol(:)
            
          endif
          
        endif
        
      endif
      
      if(idnode.eq.0)then
        
        write(nrite,
     x    "(/,/,'long range correction for: vdw energy  ',e15.6,/,
     x    25x,': vdw pressure',e15.6)")elrc/engunit,
     x    prsunt*virlrc/(-3.d0*volm)
      
        if(lghost)
     x    write(nrite,
     x    "(/,/,'long range correction for: vdw energy  ',e15.6,/,
     x    25x,': vdw pressure',e15.6)")elrc2/engunit,
     x    prsunt*virlrc2/(-3.d0*volm)
        
      endif
      
c     deallocate work arrays
      
      if(lfree.or.lghost)
     x  deallocate (elrc_exc0,numtyp_fre,numfrz_fre,stat=fail(1))
      
      if(lsolva)then
        
        deallocate (elrc_sol0,numtyp_sol0,numfrz_sol0,stat=fail(2))
        
        if(lghost)then
          
          deallocate (numtyp_sol1,numfrz_sol1,stat=fail(3))
          deallocate (numtyp_sol2,numfrz_sol2,stat=fail(4))
          
        endif
        
      endif
CVAM
CVAM      call VTEND(148, ierr)
CVAM
      return
      end subroutine lrcorrect


cRS special lrcorrect for lambda dynamics (GCMD) 

      subroutine lrcorrect_lamb
     x  (idnode,imcon,keyfce,natms,
     x  ntpatm,ntpvdw,elrc,virlrc,rcut,volm)
      

c*************************************************************************
c     
c     DL_POLY subroutine to evaluate long-range corrections to
c     pressure and energy in a periodic system.
c
cRS   pydlpoly variant to be called every step if lambda is changing
c     changes in lambda essentially mean a change in the number density and thus in the lrcorrect energy
c     in addition a force on lambda results from the lr correction
c 
c     NOTE: this works NOT with solvation or free energy mode of DL_Poly (currently unsupported in pyldlpoly)
cRS
c     
c     copyright daresbury laboratory 1993
c     author    - t. forester may 1993
c     adapted   - p.-a. cazade oct 2007: solvation, free energy etc
c     revised   - R. Schmid 2015 for lambda calcs
c     
c***************************************************************************
      
      implicit none

      integer idnode,imcon,keyfce,natms,ntpatm,i,ka,ntpvdw
      real(8) natyp,nbtyp,nctyp,ndtyp,nafrz,nbfrz,ncfrz,ndfrz
      integer ivdw,j,k,it,jt,kt,fail
      real(8) elrc,virlrc,rcut,volm,twopi,eadd,padd
      real(8) denprd,aaa,bbb,ccc,ddd,eee,eps,sig,rr0,ann,amm
      real(8) denprd1,denprd2,denprd3,denprdf

      real(8) lamb, lamb_bs, denprd_dlam_A, denprd_dlam_B
      real(8) denprd_dlam, factor     

      twopi=2.0d0*pi
            
c     initalise counter arrays      
      do i=1,ntpatm
        numtyp(i)=0
        mol_numtyp_bs(i)=0
        numfrz(i)=0        
      enddo      
      
c     evaluate number density in system
      lamb_bs = -1.0d0
      do i=1,natms
        lamb=mol_lamb_vdw(mol_which(i))
        ka=ltype(i)
        if (lamb.eq.1.0d0) then
            numtyp(ka)=numtyp(ka)+1
        else if (lamb.gt.0.0d0) then
            mol_numtyp_bs(ka)=mol_numtyp_bs(ka)+1
c            if (lamb_bs.ge.0.0d0) then
c                if (lamb_bs.ne.lamb) then
c                    write (*,*) "WARNING!! Different lambda are not allowed"
c                endif
c            endif
            lamb_bs = lamb
        endif
        if(lstfrz(i).ne.0)numfrz(ka)=numfrz(ka)+1        
      enddo
      
      
c     number densities      
      do i=1,ntpatm
        dens(i)=(dble(numtyp(i))+(dble(mol_numtyp_bs(i))*lamb_bs))/volm
      enddo
      
c     long range corrections to energy and pressure
      elrc=0.d0
      virlrc=0.d0
      denprdf=0.d0
      
      mol_dlam_lrcorr=0.0d0
      
      if(imcon.ne.0.and.imcon.ne.6.and.ntpvdw.gt.0) then 
        if(mod(keyfce,2).eq.1) then
          ivdw=0
          do i=1,ntpatm
            do j=1,i
              eadd=0.d0
              padd=0.d0
              ivdw=ivdw+1
              k=lstvdw(ivdw)
              if(ltpvdw(k).eq.0) then
                
c     tabulated potential
                eadd=prmvdw(k,1)
                padd=-prmvdw(k,2)
                
              else if(ltpvdw(k).eq.1) then
                
c     12-6 potential
                aaa=prmvdw(k,1)
                bbb=prmvdw(k,2)
                eadd=aaa/(9.d0*rcut**9)-bbb/(3.d0*rcut**3)
                padd=12.d0*aaa/(9.d0*rcut**9)-6.d0*bbb/(3.d0*rcut**3)
                
              else if(ltpvdw(k).eq.2) then
                
c     Lennard Jones potential
                eps=prmvdw(k,1)
                sig=prmvdw(k,2)
                eadd=4.d0*eps*(sig**12/(9.d0*rcut**9)-
     x            sig**6/(3.d0*rcut**3))
                padd=4.d0*eps*(12.d0*sig**12/(9.d0*rcut**9)-
     x            2.d0*sig**6/(rcut**3))
                
              else if(ltpvdw(k).eq.3) then
                
c     n - m potential
                eps=prmvdw(k,1)
                ann=prmvdw(k,2)
                amm=prmvdw(k,3)
                rr0=prmvdw(k,4)
                eadd=eps/(ann-amm)*(amm*rr0**ann/((ann-3.d0)*
     x            rcut**(ann-3.d0))-ann*rr0**amm/((amm-3.0d0)*
     x            rcut**(amm-3.d0)))
                padd=eps/(ann-amm)*ann*amm*(rr0**ann/((ann-3.d0)*
     x            rcut**(ann-3.d0))-rr0**amm/((amm-3.0d0)*
     x            rcut**(amm-3.d0)))
                
              else if((ltpvdw(k).eq.4).or.(ltpvdw(k).eq.10)) then
                
c     buckingham exp - 6 potential
cRS   or damped exp-6 potential in MOF-FF
                ccc=prmvdw(k,3)
                eadd=-ccc/(3.d0*rcut**3)
                padd=-2.d0*ccc/(rcut**3)
                
              else if(ltpvdw(k).eq.5) then
                
c     born huggins meyer exp -6 - 8  potential
                
                ddd=prmvdw(k,4)
                eee=prmvdw(k,5)
                eadd=-ddd/(3.d0*rcut**3)-eee/(5.d0*rcut**5)
                padd=-2.d0*ddd/(rcut**3)-8.d0*eee/(5.d0*rcut**5)
                
              else if(ltpvdw(k).eq.6) then
                
c     hydrogen bond  12 - 10 potential
                
                aaa=prmvdw(k,1)
                bbb=prmvdw(k,2)
                eadd=aaa/(9.d0*rcut**9)-bbb/(7.d0*rcut**7)
                padd=12.d0*aaa/(9.d0*rcut**9)-1.d1*bbb/(7.d0*rcut**7)
                
              endif
              
              if(i.ne.j) then
                
                eadd=eadd*2.d0
                padd=padd*2.d0
                
              endif

              factor = twopi/volm**2
                
              denprd=factor*(dble(numtyp(i))*dble(numtyp(j))-
     x            dble(numfrz(i))*dble(numfrz(j)))
              
cRS  compute contrib from black sheep and the force on it from the long range correction
              denprd_dlam_A = dble(numtyp(i)*mol_numtyp_bs(j)+
     x                             numtyp(j)*mol_numtyp_bs(i))
              denprd_dlam_B = dble(mol_numtyp_bs(i)*mol_numtyp_bs(j))*
     x                             lamb_bs
              denprd_dlam = factor*(denprd_dlam_A+2.0d0*denprd_dlam_B)
                            
              denprd=denprd+factor*(denprd_dlam_A+denprd_dlam_B)*lamb_bs
              
              elrc=elrc+volm*denprd*eadd
              virlrc=virlrc-denprd*padd*volm
              
              mol_dlam_lrcorr = mol_dlam_lrcorr+volm*denprd_dlam*eadd
                            
            enddo
          enddo
cRS end of double loop over atom types i, j

c          write (*,*) elrc, lamb_bs, mol_dlam_lrcorr
                    
          volm_sav=volm
                    
        endif
      endif
cRS end of if clause checking imcon and whether keyfce requests calc of forces
            
      return
      end subroutine lrcorrect_lamb



!       subroutine srfrceneu
!      x  (lsolva,lfree,lghost,ik,engsrp,virsrp,dlrpot,rcut)
! 
! c***********************************************************************
! c     
! c     dl_poly subroutine for calculating short range force and
! c     potential energy terms using verlet neighbour list
! c     neutral group implementation
! c     
! c     parallel replicated data version
! c     
! c     copyright - daresbury laboratory 1992
! c     author    - w. smith       march 1992
! c     
! c     neutral groups
! c     author    - t. forester    march  1994
! c     
! c     adapted   - p.-a. cazade oct 2007: solvation, free energy etc
! c     
! c     wl
! c     2009/01/13 11:22:05
! c     1.7
! c     Exp
! c     
! c***********************************************************************
!       
!       implicit none
! 
!       logical lsolva,lfree,lghost,lselect,lskip
!       integer ik,m,iatm,jatm,l,k,kkk
!       real(8) engsrp,virsrp,dlrpot,rcut,rcsq,fx,fy,fz,omega,omega_exc
!       real(8) strs1,strs2,strs3,strs5,strs6,strs9,ai,aj,ak,rsq
!       real(8) rrr,ppp,vk0,vk1,vk2,t1,t2,gk0,gk1,gk2,rdlpot,gamma
!       
! CVAM
! CVAM      call VTBEGIN(115, ierr)
! CVAM
!       lskip=(lfree.or.lghost)
! 
! c     set cutoff condition for pair forces
!       
!       rcsq=rcut**2
! 
! c     reciprocal of interpolation spacing
! 
!       rdlpot=1.d0/dlrpot
! 
! c     initialise stress tensor accumulators
! 
!       strs1=0.d0
!       strs2=0.d0
!       strs3=0.d0
!       strs5=0.d0
!       strs6=0.d0
!       strs9=0.d0
! 
! c     initialise potential energy and virial
!       
!       engsrp=0.d0
!       virsrp=0.d0
! 
! c     start of primary loop for forces evaluation
!       
!       do m=1,ik
! 
! c     atomic and potential function indices
!         
!         iatm=ilist(m)
!         jatm=jlist(m)
!         
!         if(lskip)then
!           if(atm_fre(iatm)*atm_fre(jatm).eq.2)cycle
!         endif
!         
!         aj=ltype(jatm)
!         ai=ltype(iatm)
! 
!         if(ai.gt.aj) then
!           ak=(ai*(ai-1.d0)*0.5d0+aj+0.5d0)
!         else
!           ak=(aj*(aj-1.d0)*0.5d0+ai+0.5d0)
!         endif
!         k=lstvdw(int(ak))
! 
!         if(abs(vvv(1,k)).gt.1.d-10)then
! 
!           rsq=rsqdf(m)
! 
!           if(rsq.lt.rcsq) then
!               
!             rrr=sqrt(rsq)
! 
! c     determine interpolation panel for force arrays
!             
!             l=int(rrr*rdlpot)
!             ppp=rrr*rdlpot-dble(l)
! 
! c     calculate interaction energy using 3-point interpolation
!             
!             vk0=vvv(l,k)
!             vk1=vvv(l+1,k)
!             vk2=vvv(l+2,k)
!             t1=vk0+(vk1-vk0)*ppp
!             t2=vk1+(vk2-vk1)*(ppp-1.0d0)
!             omega=t1+(t2-t1)*ppp*0.5d0
! 
! c     calculate forces using 3-point interpolation
!             
!             gk0=ggg(l,k)
!             gk1=ggg(l+1,k)
!             gk2=ggg(l+2,k)
!             t1=gk0+(gk1-gk0)*ppp
!             t2=gk1+(gk2-gk1)*(ppp-1.0d0)
!             gamma=(t1+(t2-t1)*ppp*0.5d0)/rsq
!             
! 
! c     set selection control
!               
!             lselect=.true.
!             
! c     set double index
!             
!             if(lsolva)kkk=loc2(atmolt(iatm),atmolt(jatm))
!             
!             if(lghost)then
!               
! c     selected excitation option
!               
!               if((atm_fre(iatm).ne.1).and.(atm_fre(jatm).ne.1))then
!                 
! c     reset selection control
!                 
!                 lselect=(atm_fre(iatm)+atm_fre(jatm).eq.0)
!                 
!                 if(lsolva)vdw_exc(kkk)=vdw_exc(kkk)+omega
!                 
!               endif
!               
!             elseif(lfree)then
!               
! c     selected free energy option
!               
!               if((atm_fre(iatm).eq.1).or.(atm_fre(jatm).eq.1))then
!                 
! c     set hamiltonian mixing parameter
!                 
!                 omega=lambda1*omega
!                 gamma=lambda1*gamma
!                 
!               elseif((atm_fre(iatm).eq.2).or.(atm_fre(jatm).eq.2))then
!                 
! c     set hamiltonian mixing parameter
!                 
!                 omega=lambda2*omega
!                 gamma=lambda2*gamma
!                 
!               endif
!               
!             endif
!             
!             if(lselect)then
!               
! c     calculate potential energy and virial
!             
!               engsrp=omega+engsrp
!               virsrp=virsrp-gamma*rsq
!               
!               if(lsolva)vdw_sol(kkk)=vdw_sol(kkk)+omega
!               
!               fx=gamma*xdf(m)
!               fy=gamma*ydf(m)
!               fz=gamma*zdf(m)
!               
!               fxx(iatm)=fxx(iatm)+fx
!               fyy(iatm)=fyy(iatm)+fy
!               fzz(iatm)=fzz(iatm)+fz
!               
!               fxx(jatm)=fxx(jatm)-fx
!               fyy(jatm)=fyy(jatm)-fy
!               fzz(jatm)=fzz(jatm)-fz
!               
! c     calculate stress tensor
!             
!               strs1=strs1+xdf(m)*fx
!               strs2=strs2+xdf(m)*fy
!               strs3=strs3+xdf(m)*fz
!               strs5=strs5+ydf(m)*fy
!               strs6=strs6+ydf(m)*fz
!               strs9=strs9+zdf(m)*fz
!               
!             endif
!             
!           endif
!           
!         endif
!         
!       enddo
!       
! c     complete stress tensor
!       
!       stress(1)=stress(1)+strs1
!       stress(2)=stress(2)+strs2
!       stress(3)=stress(3)+strs3
!       stress(4)=stress(4)+strs2
!       stress(5)=stress(5)+strs5
!       stress(6)=stress(6)+strs6
!       stress(7)=stress(7)+strs3
!       stress(8)=stress(8)+strs6
!       stress(9)=stress(9)+strs9
! CVAM
! CVAM      call VTEND(115, ierr)
! CVAM
!       return
!       end subroutine srfrceneu

      end module vdw_module
