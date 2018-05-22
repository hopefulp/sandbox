      module bonds_module
      
c***********************************************************************
c     
c     dl_poly module for defining bond potential arrays
c     copyright - daresbury laboratory
c     
c     author    - w. smith    sep 2003
c     adapted for solvation, free energy and excitation
c     - p.-a. cazade oct 2007
c     
c     wl
c     2008/12/23 10:29:11
c     1.5
c     Exp
c     
c***********************************************************************
      
      use config_module
      use error_module
      use parse_module
      use setup_module
      use site_module
      use solvation_module
      use utility_module
      
      implicit none
      
      real(8), allocatable :: prmbnd(:,:)
      integer, allocatable :: listbnd(:,:)
      integer, allocatable :: numbonds(:),keybnd(:),lstbnd(:,:)
      
      save prmbnd,listbnd,numbonds,keybnd,lstbnd
      
      contains
      
      subroutine alloc_bnd_arrays(idnode)
      
      implicit none
      
      integer i,fail,idnode
      dimension fail(5)
      
      do i=1,5
        fail(i)=0
      enddo
      
      allocate (prmbnd(mxtbnd,mxpbnd),stat=fail(1))    !! contains bond params per type (listbnd(i,1) indexes)
      allocate (numbonds(mxtmls),stat=fail(2))         !! only used in define_system_module
      allocate (keybnd(mxtbnd),stat=fail(3))           !! contains bond type (key) (listbnd(i,1) indexes)
      allocate (lstbnd(mxtbnd,3),stat=fail(4))         !! only used in define_system_module to set up
      allocate (listbnd(mxbond,3),stat=fail(5))        !! contains actual bonds (index, atom1, atom2) NOTE: in orig code 4 fields .. no idea why
      
      do i=1,5
        if(fail(i).gt.0)call error(idnode,1030)
      enddo
      
      do i=1,mxtmls
        numbonds(i)=0
      enddo
      
      end subroutine alloc_bnd_arrays
      
      subroutine define_bonds
     x  (safe,idnode,itmols,nbonds,nsite,engunit)
      
c***********************************************************************
c     
c     dl_poly subroutine for defining bonds
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     
c     wl
c     2008/12/23 10:29:11
c     1.5
c     Exp
c     
c***********************************************************************
      
      implicit none
      
      logical safe
      character*8 keyword
      character*1 message(80)
      integer idnode,itmols,nbonds,nsite,ntmp,ibond,ibond1
      integer iatm1,iatm2,isite1,isite2,idum,i,j
      real(8) engunit
      
      ntmp=intstr(record,lenrec,idum)
      numbonds(itmols)=numbonds(itmols)+ntmp
      if(idnode.eq.0)then
        write(nrite,"(/,1x,'number of chemical bonds',
     x    7x,i10)")ntmp
        write(nrite,"(/,/,1x,'chemical bond details:',
     x    /,/,21x,7x,'key',5x,'index',5x,'index',28x,
     x    'parameters', /)")
      endif
      
      ibond1=numbonds(itmols)
      do ibond=1,ibond1
        
        call getrec(safe,idnode,nfield)
        if(.not.safe)return
        
        call copystring(record,message,80)
        call lowcase(record,4)
        call getword(keyword,record,4,lenrec)
        iatm1=intstr(record,lenrec,idum)
        iatm2=intstr(record,lenrec,idum)
        
c     test for frozen atom pairs
        
        isite1=nsite-numsit(itmols)+iatm1
        isite2=nsite-numsit(itmols)+iatm2
        
        if(lfzsit(isite1)*lfzsit(isite2).ne.0)then
          
          numbonds(itmols)=numbonds(itmols)-1
          if(idnode.eq.0)write(nrite,'(12x,a16,40a1)')
     x      '*** frozen *** ',(message(i),i=1,40)
          
        else
          
          nbonds=nbonds+1
          if(nbonds.gt.mxtbnd)call error(idnode,30)
          
          if    (keyword(1:4).eq.'harm')then
            keybnd(nbonds)=1
          elseif(keyword(1:4).eq.'-hrm')then
            keybnd(nbonds)=-1
          elseif(keyword(1:4).eq.'mors')then
            keybnd(nbonds)=2
          elseif(keyword(1:4).eq.'-mrs')then
            keybnd(nbonds)=-2
          elseif(keyword(1:4).eq.'12-6')then
            keybnd(nbonds)=3
          elseif(keyword(1:4).eq.'-126')then
            keybnd(nbonds)=-3
          elseif(keyword(1:4).eq.'rhrm')then
            keybnd(nbonds)=4
          elseif(keyword(1:4).eq.'-rhm')then
            keybnd(nbonds)=-4
          elseif(keyword(1:4).eq.'quar')then
            keybnd(nbonds)=5
          elseif(keyword(1:4).eq.'-qur')then
            keybnd(nbonds)=-5
          elseif(keyword(1:4).eq.'buck')then
            keybnd(nbonds)=6
          elseif(keyword(1:4).eq.'-bck')then
            keybnd(nbonds)=-6
          elseif(keyword(1:4).eq.'fene')then
            keybnd(nbonds)=7
          elseif(keyword(1:4).eq.'-fen')then
            keybnd(nbonds)=-7
cRS MM3 quartic bond with predefined cubic/quartic terms
          elseif(keyword(1:4).eq.'mm3b')then
            keybnd(nbonds)=8
          elseif(keyword(1:4).eq.'-mm3')then
            keybnd(nbonds)=-8
cJPD variable quartic potential in the MM3 allinger form
          elseif(keyword(1:4).eq.'aqua')then
            keybnd(nbonds)=9
cRS
          else
            if(idnode.eq.0)write(nrite,*)message
            call error(idnode,444)
          endif
          
          lstbnd(nbonds,1)=iatm1
          lstbnd(nbonds,2)=iatm2
          prmbnd(nbonds,1)=dblstr(record,lenrec,idum)
          prmbnd(nbonds,2)=dblstr(record,lenrec,idum)
          prmbnd(nbonds,3)=dblstr(record,lenrec,idum)
          prmbnd(nbonds,4)=dblstr(record,lenrec,idum)
          
          if(idnode.eq.0)
     x      write(nrite,"(27x,a4,2i10,2x,1p,10e15.6)")
     x      keyword(1:4),lstbnd(nbonds,1),
     x      lstbnd(nbonds,2),(prmbnd(nbonds,j),j=1,mxpbnd)
c     
c     convert energy units to internal units
          
          if(abs(keybnd(nbonds)).eq.3)then
            prmbnd(nbonds,2)=prmbnd(nbonds,2)*engunit
          endif
          if(abs(keybnd(nbonds)).eq.5)then
            prmbnd(nbonds,3)=prmbnd(nbonds,3)*engunit
            prmbnd(nbonds,4)=prmbnd(nbonds,4)*engunit
          endif
          if(abs(keybnd(nbonds)).eq.6)then
            prmbnd(nbonds,3)=prmbnd(nbonds,3)*engunit
          endif
          
          prmbnd(nbonds,1)=prmbnd(nbonds,1)*engunit
          
        endif
        
      enddo
      
      return
      end subroutine define_bonds
      
      subroutine bndfrc
     x  (lsolva,lfree,lexcite,idnode,imcon,mxnode,ntbond,engbnd,virbnd)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating chemical bond energy and 
c     force terms in molecular dynamics.
c     
c     replicated data - blocked  data version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith        july 1992
c     modified  - t. forester    march 1993 
c     modified  - t. forester    march 1994 
c     modified  - t. forester    may   1994 
c     modified  - t. forester    nov   1994 
c     modified  - w. smith       nov   2006
c     modified  - p.-a. cazade   oct   2007, solvation etc.
c     modified  - r. schmid      sep   2009, mm3 bond
c     
c     wl
c     2008/12/23 10:29:11
c     1.5
c     Exp
c     
c***********************************************************************
      
      implicit none
      
      logical safe,lsolva,lfree,lexcite,lselect
      integer i,fail,ibnd1,ibnd2,idnode,mxnode,ii,ia,ib,imcon
      integer keyb,kk,ntbond
      real(8) strs(6)
      real(8) rab,rrab,omega,gamma,fx,fy,fz,engbnd,virbnd
cRS
      real(8) drab, drab2
cRS
      real(8), allocatable :: xdab(:),ydab(:),zdab(:)
CVAM
CVAM      call VTBEGIN(21, ierr)
CVAM
      allocate (xdab(msbad),ydab(msbad),zdab(msbad),stat=fail)
      if(fail.ne.0)call error(idnode,1040)
      
c     flag for undefined potential
      
      safe=.true.
      
c     check size of work arrays
      
      if((ntbond-mxnode+1)/mxnode.gt.msbad)call error(idnode,418)
      
c     block indices
      
      ibnd1=(idnode*ntbond)/mxnode+1
      ibnd2=((idnode+1)*ntbond)/mxnode
      
c     initialise accumulators
      
      engbnd=0.d0
      virbnd=0.d0
      bnd_fre=0.d0
      do i=1,6
        strs(i)=0.d0
      enddo

      if(lsolva)then
        
        lcomp(1)=.true.
        bnd_sol(:)=0.d0
        if(lexcite)bnd_exc(:)=0.d0
        
      endif
      
c     calculate atom separation vectors
      
      ii=0
      do i=ibnd1,ibnd2
        
        ii=ii+1
        
c     indices of bonded atoms
        
        ia=listbnd(ii,2)
        ib=listbnd(ii,3)
        
c     components of bond vector
        
        xdab(ii)=xxx(ia)-xxx(ib)
        ydab(ii)=yyy(ia)-yyy(ib)
        zdab(ii)=zzz(ia)-zzz(ib)
        
      enddo
      
c     periodic boundary condition
      
      call images(imcon,0,1,ii,cell,xdab,ydab,zdab)
      
c     loop over all specified chemical bond potentials
      
      ii=0
      do i=ibnd1,ibnd2
        
        ii=ii+1
        
c     define components of bond vector
        
        rrab=0.d0
        rab=sqrt(xdab(ii)**2+ydab(ii)**2+zdab(ii)**2)
        if(rab.gt.1.d-6)rrab=1.d0/rab
        
c     index of potential function parameters
        
        kk=listbnd(ii,1)
        keyb=abs(keybnd(kk))
        
c     calculate scalar constant terms
        
        if(keyb.eq.0)then
          
c     null interaction
          
          omega=0.d0
          gamma=0.d0
          
        elseif(keyb.eq.1)then
          
c     harmonic potential
          
          omega=0.5d0*prmbnd(kk,1)*(rab-prmbnd(kk,2))**2
          gamma=prmbnd(kk,1)*(rab-prmbnd(kk,2))*rrab
          
        else if(keyb.eq.2)then
          
c     morse potential
          
c          omega=prmbnd(kk,1)*((1.d0-exp(-prmbnd(kk,3)*
c     x      (rab-prmbnd(kk,2))))**2-1.d0)
c put minimum to zero (=remove constant term -prmbnd(kk,1) to make it equal to tinker
          omega=prmbnd(kk,1)*((1.d0-exp(-prmbnd(kk,3)*
     x      (rab-prmbnd(kk,2))))**2)
          gamma=2.d0*prmbnd(kk,1)*prmbnd(kk,3)*(1.d0-
     x      exp(-prmbnd(kk,3)*(rab-prmbnd(kk,2))))*
     x      exp(-prmbnd(kk,3)*(rab-prmbnd(kk,2)))*rrab

        else if(keyb.eq.3)then
          
c     12-6 potential
          
          omega=(prmbnd(kk,1)*rrab**6-prmbnd(kk,2))*rrab**6
          gamma=(6.d0*prmbnd(kk,2)-12.d0*prmbnd(kk,1)*rrab**6)*
     x      rrab**8
          
        elseif(keyb.eq.4)then
          
c     restrained harmonic
          
          rab=rab-prmbnd(kk,2)
          omega=0.5d0*prmbnd(kk,1)*(min(abs(rab),prmbnd(kk,3)))**2
     x      +prmbnd(kk,1)*prmbnd(kk,3)*max(abs(rab)-prmbnd(kk,3),0.d0)
          gamma=rrab*prmbnd(kk,1)*(sign(min(abs(rab),prmbnd(kk,3)),rab))
          
        elseif(keyb.eq.5)then
          
c     quartic potential
          
          omega=0.5d0*prmbnd(kk,1)*(rab-prmbnd(kk,2))**2+
     x      1.d0/3.d0*prmbnd(kk,3)*(rab-prmbnd(kk,2))**3+
     x      0.25d0*prmbnd(kk,4)*(rab-prmbnd(kk,2))**4
          gamma=rrab*(prmbnd(kk,1)*(rab-prmbnd(kk,2))+
     x      prmbnd(kk,3)*(rab-prmbnd(kk,2))**2+
     x      prmbnd(kk,4)*(rab-prmbnd(kk,2))**3)
          
        else if(keyb.eq.6)then
          
c     buckingham exp-6 potential
          
          omega=prmbnd(kk,1)*exp(-rab/prmbnd(kk,2))-prmbnd(kk,3)*
     x      rrab**6
          gamma=-rrab*prmbnd(kk,1)*exp(-rab/prmbnd(kk,2))/prmbnd(kk,2)+
     x      6.d0*prmbnd(kk,3)*rrab**8
          
        else if(keyb.eq.7)then
          
c     FENE bond potential
          
          omega=-0.5d0*prmbnd(kk,1)*prmbnd(kk,2)**2*log(1.d0-
     x      ((rab-prmbnd(kk,3))/prmbnd(kk,2))**2)
          gamma=rrab*prmbnd(kk,1)*(rab-prmbnd(kk,3))/
     x      (1.d0-((rab-prmbnd(kk,3))/prmbnd(kk,2))**2)
          
        elseif(keyb.eq.8)then
          
cRS     mm3 quartic potential with fixed cubic and quartic terms
cRS       cubic: -2.55
cRS       quartic: 3.793125 = (7/12)*bond-cubic^2
          drab = rab-prmbnd(kk,2)
          drab2 = drab*drab
          omega=0.5d0*prmbnd(kk,1)*drab2 * (1.0d0-2.55d0*drab+ 
     x              3.793125d0*drab2)
          gamma=prmbnd(kk,1)*drab*rrab*(1.0d0-1.5d0*2.55d0*drab+
     x              2.0d0*3.793125d0*drab2)
        
        elseif(keyb.eq.9)then
cJPD    variable quartic potential in the allinger MM3 formulation
          drab = rab-prmbnd(kk,2)
          drab2 = drab*drab
          omega=0.5d0*prmbnd(kk,1)*drab2 * (1.0d0-prmbnd(kk,3)*
     x              drab + prmbnd(kk,4)*(prmbnd(kk,3)**2)*drab2)
          gamma=prmbnd(kk,1)*drab*rrab*(1.0d0-1.5d0*prmbnd(kk,3)*
     x              drab + 2.0d0*prmbnd(kk,4)*(prmbnd(kk,3)**2)*drab2)     
        else
          
c     undefined potential
          
          omega=0.d0
          gamma=0.d0
          safe=.false.
          
        endif
        
c     indices of bonded atoms
        
        ia=listbnd(ii,2)
        ib=listbnd(ii,3)
        
c     set selection control
        
        lselect=.true.
        
        if(lexcite)then
          
c     selected excitation option
          
          if((atm_fre(ia).ne.1).and.(atm_fre(ib).ne.1))then
            
c     reset selection control
            
            lselect=(atm_fre(ia)+atm_fre(ib).eq.0)
            
            if(lsolva)then
              bnd_exc(atmolt(ia))=bnd_exc(atmolt(ia))+omega
            endif
            
          endif
          
        elseif(lfree)then
          
c     selected free energy option
          
          if((atm_fre(ia).eq.1).or.(atm_fre(ib).eq.1))then
            
c     set hamiltonian mixing parameter
            
            bnd_fre=bnd_fre-omega
            omega=lambda1*omega
            gamma=lambda1*gamma
            
          elseif((atm_fre(ia).eq.2).or.(atm_fre(ib).eq.2))then
            
c     set hamiltonian mixing parameter
            
            bnd_fre=bnd_fre+omega
            omega=lambda2*omega
            gamma=lambda2*gamma
            
          endif
          
        endif
        
        if(lselect)then
          
c     calculate bond energy and virial
        
          engbnd=engbnd+omega
          virbnd=virbnd+gamma*rab*rab
          
c     calculate solvation energy
        
          if(lsolva)then
            bnd_sol(atmolt(ia))=bnd_sol(atmolt(ia))+omega
          endif
          
c     calculate forces
          
          fx=-gamma*xdab(ii)
          fy=-gamma*ydab(ii)
          fz=-gamma*zdab(ii)
          
          fxx(ia)=fxx(ia)+fx
          fyy(ia)=fyy(ia)+fy
          fzz(ia)=fzz(ia)+fz
          
          fxx(ib)=fxx(ib)-fx
          fyy(ib)=fyy(ib)-fy
          fzz(ib)=fzz(ib)-fz
          
c     calculate stress tensor
        
          strs(1)=strs(1)+xdab(ii)*fx
          strs(2)=strs(2)+xdab(ii)*fy
          strs(3)=strs(3)+xdab(ii)*fz
          strs(4)=strs(4)+ydab(ii)*fy
          strs(5)=strs(5)+ydab(ii)*fz
          strs(6)=strs(6)+zdab(ii)*fz
          
        endif
        
      enddo
      
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
      
c     check for undefined potentials
      
      if(mxnode.gt.1)call gstate(safe)
      if(.not.safe)call error(idnode,444)
      
c     sum contributions to potential and virial
      
      if(mxnode.gt.1)then
        
        buffer(1)=engbnd
        buffer(2)=virbnd
        buffer(3)=bnd_fre
        call gdsum(buffer(1),3,buffer(4))        
        engbnd=buffer(1)
        virbnd=buffer(2)
        bnd_fre=buffer(3)
        
        if(lsolva)then
          
          call gdsum(bnd_sol(1),mxtmls,buffer(1))
          if(lexcite)call gdsum(bnd_exc(1),mxtmls,buffer(1))
          
        endif
        
      endif
      
      deallocate (xdab,ydab,zdab,stat=fail)
CVAM
CVAM      call VTEND(21, ierr)
CVAM
      return
      end subroutine bndfrc
      
      end module bonds_module
