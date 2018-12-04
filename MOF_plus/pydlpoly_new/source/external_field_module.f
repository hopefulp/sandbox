      module external_field_module
      
c***********************************************************************
c     
c     dl_poly module for defining external field potential arrays
c     copyright - daresbury laboratory
c     author    - w. smith    oct 2003
c     
c     wl
c     2008/12/23 10:29:12
c     1.4
c     Exp
c     
c     RS Oct 2013
c      added "double wall potential" which adds two "harmonic walls" acting
c      differently on the two sets of atoms (first set is 1 to split_index-1
c      and the second is from split_index to the last atom)
c      NOTE: this works only in pyldpoly becaus the settings must be made directly 
c            into the fortran variables. i can not be requested from the FIELD file
c
c***********************************************************************
      
      use config_module
      use error_module
      use parse_module
      use setup_module
      use utility_module

cRS   for molecular dipoles
      use molecule_module
      use rigid_body_module
      
      implicit none
      
      real(8), allocatable :: prmfld(:)
cRS added for double wall potential
      integer              :: split_index
cRS   total dipole moment
      real(8)              :: dip(3), dip_ref(3)
      
      logical              :: lcalc_dip
      data lcalc_dip/.false./      
      save prmfld, split_index, dip, lcalc_dip, dip_ref
      
      contains
      
      subroutine alloc_fld_arrays(idnode)
      
      implicit none
      
      integer fail,idnode
      
      data fail/0/
      
      allocate (prmfld(mxfld),stat=fail)
      if(fail.ne.0)call error(idnode,1200)
      
      end subroutine alloc_fld_arrays
      
      subroutine define_external_field
     x  (safe,lunits,idnode,keyfld,engunit)
      
c***********************************************************************
c     
c     dl_poly subroutine to define external fields
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     
c     wl
c     2008/12/23 10:29:12
c     1.4
c     Exp
c     
c***********************************************************************
      
      implicit none
      
      logical safe,lunits
      character*8 keyword
      character*1 message(80)
      integer idnode,keyfld,nfld,i,k,idum
      real(8) engunit
      
      call getrec(safe,idnode,nfield)
      if(.not.safe)return
      
      call strip(record,lenrec)
      call lowcase(record,lenrec)
      call copystring(record,message,80)
      call getword(keyword,record,4,lenrec)
      
      if(keyword(1:4).eq.'elec') then
        keyfld =1 
      elseif(keyword(1:4).eq.'oshr') then
        keyfld=2
      elseif(keyword(1:4).eq.'shrx') then
        keyfld=3
      elseif(keyword(1:4).eq.'grav') then
        keyfld=4
      elseif(keyword(1:4).eq.'magn') then
        keyfld=5
      elseif(keyword(1:4).eq.'sphr') then
        keyfld=6
      elseif(keyword(1:4).eq.'zbnd') then
        keyfld=7
      else
        if(idnode.eq.0) write(nrite,*) message
        call error(idnode,454)
      endif
      
      do i = 1,mxfld
        prmfld(i)=0.d0
      enddo
      
      nfld=intstr(record,lenrec,idum)
      if(nfld.eq.0)nfld=5
      call getrec(safe,idnode,nfield)
      if(.not.safe)return
      do k=1,nfld
        
        prmfld(k)=dblstr(record,lenrec,idum)
        if(idum.gt.lenrec.and.k.lt.nfld)then
          
          call getrec(safe,idnode,nfield)
          if(.not.safe)return
          
        endif
        
      enddo
      
      if(idnode.eq.0) then
        
        write(nrite,"(/,/,1x,'external field key ',13x,a4,
     x    /,/,30x,'external field parameters')") keyword(1:4)
        write(nrite,"(2(/,1x,1p,5e15.5))") prmfld
        
      endif      
      
c     convert to internal units
      
      if(keyfld.eq.1.or.keyfld.eq.4.or.keyfld.eq.5) then
        
        if(.not.lunits)call error(idnode,6)
        
        do i = 1,3
          prmfld(i) = prmfld(i)*engunit
        enddo
        
      elseif(keyfld.eq.2.or.keyfld.eq.6.or.keyfld.eq.7) then
        
        prmfld(1) = prmfld(1)*engunit
        
      endif
      
      return
      end subroutine define_external_field
      
      subroutine extnfld
     x  (idnode,imcon,keyfld,mxnode,natms,engfld,virfld)
      
c***********************************************************************
c     
c     dl_poly routine for application of an external field
c     
c     replicated data version / block data
c     
c     copyright daresbury laboratory 1993
c     author -    t.forester october 1993
c     amended-    t.forester dec 1994
c     
c     wl
c     2008/12/23 10:29:12
c     1.4
c     Exp
c     
c***********************************************************************
      
      implicit none
      
      integer idnode,imcon,keyfld,mxnode,natms,iatm1,iatm2,i
      real(8) engfld,virfld,rz,rrr,gamma,zdif
      
cRS   helpers for double wall potential
      real(8) z1, z2, V1, V2, z1h, z2h, z1i, s1, s2, E
cRS   variables for constantD (keyfld=10)
      real(8) epsilon0, efact, DmP(3), Ef(3), diag, efact2, D(3)

      real(8) buffer(3)
    
CVAM
CVAM      call VTBEGIN(27, ierr)
CVAM
      
c     energy and virial accumulators 
      
      engfld = 0.d0
      virfld = 0.d0
    
      
c     block indices
      
      iatm1 = (idnode*natms)/mxnode + 1
      iatm2 = ((idnode+1)*natms)/mxnode
      
      if(keyfld.eq.1) then
        
c     electric field: prmfld(1-3) are field components
c     constant electric filed: E_fld = -V_cell*E*P (we use the cell dipole computed in comp_polarization)
c     here (means the ltrackimg and lcalc_pol must be set properly for this to work!)
c     NOTE: the values in dip are NOT divided by the volume
c     E is the fixed electric field in units giving DL_poly energy units (0.1 kJ/mol) when multiplied with
c     a dipole in unitcharge*Angstrom
c     the corresponding force is E*q_i
                
        do i = iatm1,iatm2
          if (pol_ignore(i).eq.0) then         
            fxx(i) = fxx(i)+ chge(i)*prmfld(1)
            fyy(i) = fyy(i)+ chge(i)*prmfld(2)
            fzz(i) = fzz(i)+ chge(i)*prmfld(3)
          endif                    
        enddo   

        engfld = engfld - prmfld(1)*dip(1)
        engfld = engfld - prmfld(2)*dip(2)
        engfld = engfld - prmfld(3)*dip(3)
c       dip is already parallel communicated .. we need to divide by the number of nodes
c       becasue it is later summed over the nodes
        engfld = engfld/mxnode
                
      elseif(keyfld.eq.10) then
        
c     constant electric displacement: prmfld(1-3) are the D components
c     E_fld = V_cell/(2.0*epsilon_0)*(D-P)^2 (P is the cell dipole computed dip in comp_polarization divided by V_cell)
c     (means the ltrackimg and lcalc_pol etc must be set properly for this to work!)
c     D is the fixed electric displacement in units matching P => [unitcharge/A^2]
c     the corresponding force is (D-P)/epsilon_0*q_i
c     the cell volume is take from celprp(10) in the config module (hoping this is properly upadeted at this point)

cJPD  we have to adjust first D to the current lattice, ie compute D out
c     of d and the the cellparameters
 
        D(1) = (prmfld(1)*cell(1)+prmfld(2)*cell(2)+prmfld(3)*cell(3))
        D(2) = (prmfld(1)*cell(4)+prmfld(2)*cell(5)+prmfld(3)*cell(6))
        D(3) = (prmfld(1)*cell(7)+prmfld(2)*cell(8)+prmfld(3)*cell(9))
        D(1) = D(1)/celprp(10)
        D(2) = D(2)/celprp(10)
        D(3) = D(3)/celprp(10)

c        write(*,*) "D", D

        DmP(1) = D(1)-(dip(1)/celprp(10))
        DmP(2) = D(2)-(dip(2)/celprp(10))
        DmP(3) = D(3)-(dip(3)/celprp(10))

c        DmP(1) = prmfld(1)-(dip(1)/celprp(10))
c        DmP(2) = prmfld(2)-(dip(2)/celprp(10))
c        DmP(3) = prmfld(3)-(dip(3)/celprp(10))

c        write (*,*) "DmP", DmP(1), DmP(2), DmP(3)
c        write(*,*) "pol",dip/celprp(10)

c           eps0 in units of unitcharge/(V*Angstrom)
        epsilon0 = 5.526348d-3
c           efact is faraday const *0.1/eps0
        efact = 9648.53082d0/epsilon0

c        write(*,*)   celprp(10)/2.0d0*efact*(DmP(1)*DmP(1)) 
c        write(*,*)   celprp(10)/2.0d0*efact*(Dmp(2)*DmP(2))
c        write(*,*)   celprp(10)/2.0d0*efact*(DmP(3)*DmP(3))
        engfld = engfld + celprp(10)/2.0d0*efact*(DmP(1)*DmP(1)) 
        engfld = engfld + celprp(10)/2.0d0*efact*(Dmp(2)*DmP(2))
        engfld = engfld + celprp(10)/2.0d0*efact*(DmP(3)*DmP(3))
c       dip is already parallel communicated .. we need to divide by the number of nodes
c       becasue it is later summed over the nodes
        engfld = engfld/mxnode

        do i = iatm1,iatm2
          if (pol_ignore(i).eq.0) then
            fxx(i) = fxx(i)+ efact*chge(i)*DmP(1)
            fyy(i) = fyy(i)+ efact*chge(i)*DmP(2)
            fzz(i) = fzz(i)+ efact*chge(i)*DmP(3)
          endif                    
        enddo

cJPD  compute maxwell stress tensor

cJPD  diag is the diagonal elements substracted on the diagonals of the
cJPD  tensor, since the DL-POLY stress tensor has units of energy and
cJPD  not of energy/V we have to multiply the whole thing by the
cJPD  current Volume

      efact2 = 1.0d0/efact 

      Ef(1) = DmP(1)*efact
      Ef(2) = DmP(2)*efact
      Ef(3) = DmP(3)*efact

c      write(*,*) "Ef",Ef
c      write(*,*) "DmP",DmP
c      write(*,*) "dip",dip
c      write(*,*) "vol",celprp(10)
c      write(*,*) "pol",dip/celprp(10)

      diag = 0.5d0 * efact2 * (Ef(1)**2+Ef(2)**2+Ef(3)**2)
      
      stress(1) = stress(1)-celprp(10)*(efact2*Ef(1)*Ef(1)-diag)
      stress(2) = stress(2)-celprp(10)*efact2*Ef(1)*Ef(2)
      stress(3) = stress(3)-celprp(10)*efact2*Ef(1)*Ef(3)
      stress(4) = stress(4)-celprp(10)*efact2*Ef(2)*Ef(1)
      stress(5) = stress(5)-celprp(10)*(efact2*Ef(2)*Ef(2)-diag)
      stress(6) = stress(6)-celprp(10)*efact2*Ef(2)*Ef(3)
      stress(7) = stress(7)-celprp(10)*efact2*Ef(3)*Ef(1)
      stress(8) = stress(8)-celprp(10)*efact2*Ef(3)*Ef(2)
      stress(9) = stress(9)-celprp(10)*(efact2*Ef(3)*Ef(3)-diag)

      
      elseif(keyfld.eq.2) then
        
c     oscillating shear: orthorhombic box:  Fx = a*cos(b.2.pi.z/L)
        
        rz = 2.d0*pi/cell(9)
        
        do i = iatm1,iatm2
          
          fxx(i) = fxx(i) + prmfld(1)*cos(prmfld(2)*zzz(i)*rz)
          
        enddo
        
      elseif(keyfld.eq.3.and.imcon.eq.6) then
        
c     continuous shear of walls : 2D periodic box (imcon=6)
c     shear rate = prmfld(1) angstrom per ps for atoms at
c     abs(z) > prmfld(2)
        
        do i=iatm1,iatm2
          
          if(abs(zzz(i)).gt.prmfld(2)) then
            
            vxx(i) = 0.5d0*sign(prmfld(1),zzz(i))
            
          endif
          
        enddo
        
      elseif(keyfld.eq.4) then
        
c     gravitational field: field components given by prmfld(1-3)
        
        do i =iatm1,iatm2
          
          fxx(i) = fxx(i) + prmfld(1)*weight(i)
          fyy(i) = fyy(i) + prmfld(2)*weight(i)
          fzz(i) = fzz(i) + prmfld(3)*weight(i)
          
        enddo
        
      elseif(keyfld.eq.5) then
        
c     magnetic field: field components given by prmfld(1-3)
        
        do i = iatm1,iatm2
          
          fxx(i)=fxx(i)+(vyy(i)*prmfld(3)-vzz(i)*prmfld(2))
     x      *chge(i)
          fyy(i)=fyy(i)+(vzz(i)*prmfld(1)-vxx(i)*prmfld(3))
     x      *chge(i)
          fzz(i)=fzz(i)+(vxx(i)*prmfld(2)-vyy(i)*prmfld(1))
     x      *chge(i)
          
        enddo
        
      elseif(keyfld.eq.6) then
        
c     containing sphere : r^(-n) potential
        
        do i = iatm1,iatm2
          
          rrr = sqrt(xxx(i)**2+yyy(i)**2+zzz(i)**2)
          if(rrr.gt.prmfld(4)) then
            rrr = prmfld(2) - rrr
            if(rrr.lt.0.d0) rrr = 0.1d0
            
            gamma  = prmfld(1)*rrr**(-prmfld(3))
            engfld = engfld + gamma
            
            gamma = -prmfld(3)*gamma/((prmfld(2)-rrr)*rrr)
            
            fxx(i)=fxx(i)+ gamma*xxx(i)
            fyy(i)=fyy(i)+ gamma*yyy(i)
            fzz(i)=fzz(i)+ gamma*zzz(i)
            
          endif
          
        enddo
        
      elseif(keyfld.eq.7) then
        
c     repulsive wall (harmonic) starting at z0
        
        do i = iatm1,iatm2
          
          if(prmfld(3)*zzz(i).gt.prmfld(3)*prmfld(2)) then
            
            zdif = zzz(i) - prmfld(2)
            gamma = -prmfld(1)*zdif
            
            fzz(i) = fzz(i) + gamma
            engfld = engfld - gamma*zdif/2.
            
          endif
          
        enddo
        
      elseif(keyfld.eq.8) then
        
cRS   double wall potential acting differntly on the two sets of atoms
c
c     wall 1: at z1 with force-factor V1
c     wall 2: at z2 with force-factor V2
c
c     prmfld:  z1, V1, z2, V2
        
        z1 = prmfld(1)
        V1 = prmfld(2)
        z2 = prmfld(3)
        V2 = prmfld(4)
        
        z1i = z1+cell(9)
        z1h = z1+(z2-z1)/2.0d0
        z2h = z2+(z1i-z2)/2.0d0
        if (z2h.ge.cell(9)) z2h = z2h-cell(9)
        
        do i = iatm1,iatm2
          
          if (i.lt.split_index) then
      
cRS         first set of atoms .. no force between z1 and z2
          
            if (zzz(i).ge.z2) then
              
              if ((zzz(i).ge.z2h).and.(z2h.gt.z2)) then

                zdif = z1i-zzz(i)
                gamma = V1*zdif
                fzz(i) = fzz(i) + gamma
                engfld = engfld + gamma*zdif/2.0d0

c                write (*,*) "set 1  + c1 ",i,atmnam(i),zzz(i),zdif
                
              else
              
                zdif = zzz(i)-z2
                gamma = V2*zdif
                fzz(i) = fzz(i) - gamma
                engfld = engfld + gamma*zdif/2.0d0

c                write (*,*) "set 1  - c2 ",i,atmnam(i),zzz(i),zdif
                                            
              endif
              
            elseif (zzz(i).lt.z1) then
            
              if ((zzz(i).lt.z2h).and.(z2h.lt.z2)) then

                zdif = zzz(i)-(z2-cell(9))
                gamma = V2*zdif
                fzz(i) = fzz(i) - gamma
                engfld = engfld + gamma*zdif/2.0d0

c                write (*,*) "set 1  - c3 ",i,atmnam(i),zzz(i),zdif

              else
              
                zdif = z1-zzz(i)
                gamma = V1*zdif
                fzz(i) = fzz(i) + gamma
                engfld = engfld + gamma*zdif/2.0d0

c                write (*,*) "set 1  + c4 ",i,atmnam(i),zzz(i),zdif
                
              endif
                          
            endif
          
          else
          
cRS         second set of atoms .. no force between z2 and z1

            if ((zzz(i).ge.z1).and.(zzz(i).lt.z2)) then

               if (zzz(i).lt.z1h) then
               
                zdif = zzz(i)-z1
                gamma = V1*zdif
                fzz(i) = fzz(i) - gamma
                engfld = engfld + gamma*zdif/2.0d0

c                write (*,*) "set 2  - c1  ",i,atmnam(i),zzz(i),zdif
                              
               else
               
                zdif = z2-zzz(i)
                gamma = V2*zdif
                fzz(i) = fzz(i) + gamma
                engfld = engfld + gamma*zdif/2.0d0
                              
c                write (*,*) "set 2  + c2 ",i,atmnam(i),zzz(i),zdif

               endif
                        
            endif
                    
          endif
          
        enddo

cRS end of double wall code        
        
      elseif(keyfld.eq.9) then
        
cRS   gaussian walls - Two gaussians potentials at z1 and z2 can be defined
c           with E1 = V1*exp(-((z-z1)/s1))**2)
c           so V1 defines the height at z1 and s1 the width (unit angstrom)
c           NOTE: both walls act equaly on all atoms
c
c     wall 1: at z1 with height V1 and width s1
c     wall 2: at z2 with height V2 and width s2
c
c     prmfld:  z1, V1, s1, z2, V2, s2
        
        z1 = prmfld(1)
        V1 = prmfld(2)
        s1 = prmfld(3)
        z2 = prmfld(4)
        V2 = prmfld(5)
        s2 = prmfld(6)
        
        do i = iatm1,iatm2
        
c         do first wall (consider boundary conditions) 
          zdif = zzz(i)-z1
          if (zdif.gt.(cell(9)/2.0d0))   zdif = zdif - cell(9)
          if (zdif.lt.(-cell(9)/2.0d0))  zdif = zdif + cell(9)
          E     = V1*exp(-(zdif/s1)**2)
          gamma = E*(-2.0d0*(zdif/s1))
          engfld = engfld + E
          fzz(i) = fzz(i) + sign(gamma, zdif)
          
c         do second wall (consider boundary conditions) 
          zdif = zzz(i)-z2
          if (zdif.gt.(cell(9)/2.0d0))   zdif = zdif - cell(9)
          if (zdif.lt.(-cell(9)/2.0d0))  zdif = zdif + cell(9)
          E     = V2*exp(-(zdif/s2)**2)
          gamma = E*(-2.0d0*(zdif/s2))
          engfld = engfld + E
          fzz(i) = fzz(i) + sign(gamma, zdif)

        enddo

      else
        
c     unidentified field potential error exit
        
        call error(idnode,454)
        
      endif
      
cRS
c     I think this was a bug ... no communication of engfld in parallel runs!!!
      if(mxnode.gt.1)then
         buffer(1) = engfld
         call gdsum(buffer(1), 1, buffer(2))
         engfld = buffer(1)
      endif

cRS
      
CVAM
CVAM      call VTEND(27, ierr)
CVAM
      return
      end subroutine extnfld

cRS ################ new stuff for polarization ###############################
cRS
cRS   needs track_image to work! so make sure this is the case on python level
cRS   NOTE: for rigid groups we do not use tracking (needs to be doen on COM level)
cRS         so for those atoms we skip the addition of img

      subroutine comp_polarization
     x         (idnode,mxnode,imcon,natms,ngrp,ntfree)

      implicit none

      integer idnode, mxnode, imcon, natms,ngrp,ntfree

      integer iatm1, iatm2, i, m
      integer ifre, ifre1, ifre2
      real(8) xyz(3), buffer(3), dcell(3), adip(3),rcell(9)
      real(8) ssx, ssy, ssz, det 

      iatm1 = (idnode*natms)/mxnode + 1
      iatm2 = ((idnode+1)*natms)/mxnode

      ifre1=(idnode*ntfree)/mxnode+1
      ifre2=((idnode+1)*ntfree)/mxnode

      dip(:) = 0.0d0
      if (lcalc_moldip) mol_dip(:,:) = 0.0d0

      if (imcon>0) then
         dcell(1) = cell(1)
         dcell(2) = cell(5)
         dcell(3) = cell(9)
      endif

      if (ngrp.eq.0) then
c this is the more direct variant if all atoms are free
        do i =iatm1,iatm2
          xyz(1) = xxx(i)
          xyz(2) = yyy(i) 
          xyz(3) = zzz(i)
          if (imcon>0) then
            if (imcon<3) then 
c             for cubic(1) and orthoromb(2) we can use the cell diagonal dcell
              xyz(1) = xyz(1) + imgx(i)*dcell(1)
              xyz(2) = xyz(2) + imgy(i)*dcell(2)
              xyz(3) = xyz(3) + imgz(i)*dcell(3)
            elseif (imcon==3) then
cJPD          transform to fractional coordinates

              call invert(cell, rcell, det)

              ssx=(rcell(1)*xyz(1)+rcell(4)*xyz(2)+rcell(7)*xyz(3))
              ssy=(rcell(2)*xyz(1)+rcell(5)*xyz(2)+rcell(8)*xyz(3))
              ssz=(rcell(3)*xyz(1)+rcell(6)*xyz(2)+rcell(9)*xyz(3))

cJPD          apply images

              ssx=ssx+imgx(i)
              ssy=ssy+imgy(i)
              ssz=ssz+imgz(i)

cJPD          transform back to real coordinates

              xyz(1)=(cell(1)*ssx+cell(4)*ssy+cell(7)*ssz)
              xyz(2)=(cell(2)*ssx+cell(5)*ssy+cell(8)*ssz)
              xyz(3)=(cell(3)*ssx+cell(6)*ssy+cell(9)*ssz)
            else
              write (*,*) "No polarization with this imcon: ", imcon
            endif
          endif
          adip(:) = xyz(:)*chge(i)
          if (lcalc_moldip) then
            m = mol_which(i)
            mol_dip(:,m) = mol_dip(:,m) + adip(:)
          endif
          if (pol_ignore(i).eq.0) then
            dip(:) = dip(:) + adip(:)
          endif
        end do
      else 
c for rigid atoms we take the rigid body atoms as such as they are "collected by image" anyway
c       first pass without correction
        do i =iatm1,iatm2
          xyz(1) = xxx(i)
          xyz(2) = yyy(i) 
          xyz(3) = zzz(i)
          adip(:) = xyz(:)*chge(i)
          if (lcalc_moldip) then
            m = mol_which(i)
            mol_dip(:,m) = mol_dip(:,m) + adip(:)
          endif
          if (pol_ignore(i).eq.0) then
            dip(:) = dip(:) + adip(:)
          endif
        end do
c       second pass with corrections only for free atoms
        do ifre=ifre1,ifre2
          i=lstfre(ifre)
          if (imcon>0) then
            if (imcon<3) then 
c             for cubic(1) and orthoromb(2) we can use the cell diagonal dcell
              xyz(1) = imgx(i)*dcell(1)
              xyz(2) = imgy(i)*dcell(2)
              xyz(3) = imgz(i)*dcell(3)
            elseif (imcon==3) then
              write (*,*) "TO BE IMPLEMNTED!!!"
            else
              write (*,*) "No polarization with this imcon: ", imcon
            endif
          endif
          adip(:) = xyz(:)*chge(i)
          if (lcalc_moldip) then
            m = mol_which(i)
            mol_dip(:,m) = mol_dip(:,m) + adip(:)
          endif
          if (pol_ignore(i).eq.0) then
            dip(:) = dip(:) + adip(:)
          endif
        end do
      endif

      call gdsum (dip, 3, buffer)
cRS   TBI: communicate mol_dip (need a buffer for it)

c     IMPORTANT: in the end subtract the refrence dipole moment (registered at the start)
      dip(:) = dip(:) - dip_ref(:)

      end subroutine comp_polarization

      
      end module external_field_module
