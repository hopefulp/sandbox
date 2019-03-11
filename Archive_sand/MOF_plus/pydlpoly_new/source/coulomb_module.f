      module coulomb_module

c***********************************************************************
c     
c     dl_poly module for defining coulomb terms
c     copyright - daresbury laboratory
c     
c     author    - w. smith    sep 2003
c     adapted for solvation, free energy and excitation
c               - p.-a. cazade oct 2007
c     
c     Added functionality for fluctuating charges (QEq) for pydlpoly project
c        including Gaussian charge distribution interactions
c        by Ch. Spickermann and R. Schmid (2010/2011)     
c     Added options to store Coulomb interactions Jij in a compact form (using neighborlists)
c        by R. Schmid (2011)
c     Comment:
c     For QEq we always need all interactions including exclusion list.
c     we store the interactions distributed (loop 1 = idnode+1 -> natoms, stride:maxnode)
c     in the same way the neighborlist is handeled. 
c
c     wl
c     2008/12/23 10:29:11
c     1.4
c     Exp
c     
c***********************************************************************

      use config_module
      use ewald_module
      use pair_module
      use property_module
      use setup_module
      use solvation_module
ccs
      use dihedral_module, only: gaupot,gauforce
ccs
cRS for switchable molecule implementation
      use molecule_module 

cRS
      use nlist_builders_module
      use exclude_module
 
      implicit none

cRS special data structures to keep Jij in a distributed compact form for QEq et al.

      real(8), allocatable ::  Jij(:)
      integer, allocatable ::  Jij_mapper(:)
      integer              ::  Jij_first
      integer              ::  Jij_size
      integer, parameter   ::  Jij_increment_size = 500
      logical              ::  store_Jij 
      logical              ::  valid_Jij

      data  store_Jij,valid_Jij/.false.,.false./
      data  Jij_size/500/

      save Jij, Jij_mapper, Jij_first, Jij_size, store_Jij, valid_Jij
cRS

      contains


cRS special function to allocate a local Jij
      subroutine init_Jij(idnode)
     
      implicit none

      integer fail, idnode 

      allocate (Jij(Jij_size), stat=fail)
c      write (*,*) "DEBUG Jij : initial allocation  ", Jij_size
      if (fail.ne.0) then
         if(idnode.eq.0)write(nrite,'(10i5)')fail
         call error(idnode,1000)
      end if
      allocate (Jij_mapper(mxatms), stat=fail)  
      if (fail.ne.0) then
         if(idnode.eq.0)write(nrite,'(10i5)')fail
         call error(idnode,1000)
      end if

      end subroutine init_Jij


      subroutine allocate_Jij(totsize,idnode)

      implicit none

      integer fail, idnode
      integer totsize

c      write (*,*) "DEBUG Jij: node , current size , requested size  ", 
c     &                              idnode, Jij_size, totsize
      if (totsize > Jij_size) then
cRS      yes we need to increase ... increase at least by the increment size
        Jij_size = Jij_size + Jij_increment_size
        if (totsize > Jij_size) then
            Jij_size = totsize
        end if
        deallocate (Jij, stat=fail)
        if (fail.ne.0) then
            if(idnode.eq.0)write(nrite,'(10i5)')fail
            call error(idnode,1000)
        end if
        allocate (Jij(Jij_size), stat=fail)
        write (*,*) "DEBUG Jij : reallocating node , new size  ", 
     &              idnode, Jij_size
        if (fail.ne.0) then
            if(idnode.eq.0)write(nrite,'(10i5)')fail
            call error(idnode,1000)
        end if
      end if
      return
      end subroutine allocate_Jij     

      subroutine create_Jij_mapper(idnode,mxnode,natms,loglnk,mapper)

      implicit none

      integer idnode,mxnode,natms
      logical loglnk
      integer i,ii,k,j,mapper 
      mapper=1 
      ii=0
      do i=idnode+1,natms,mxnode
         ii=ii+1
         Jij_mapper(i)=mapper
c        now compute how many neighbors we have and increment mapper accordingly
         mapper = mapper+lentry(ii)
c        in case of linked cells we need to check it 
         if (loglnk) then
           do k=1,nexatm(ii)
             j = lexatm(ii,k)
             if (j.gt.i) then
                mapper = mapper+1
             end if
           end do
         else
           mapper = mapper+nexatm(ii)
         end if
      enddo

      end subroutine create_Jij_mapper

      subroutine zero_Jij()
      
      implicit none
      
      integer i
      
      do i=1,Jij_size
         Jij(i)=0.0d0
      end do
      
      end subroutine zero_Jij

      subroutine apply_Jij(idnode,mxnode,natms,loglnk,lcon,lconmol)

      implicit none

      integer idnode, mxnode, natms
      logical loglnk, lcon, lconmol
      integer debug

      integer i,j,ii,jj,m
      integer Jij_index
      real(8) term
      real(8) molpot(mx_molecules)

c      real(8) t1, t0
      
c      call timchk(0, t0)
      
c     initialize pot array (where we sum up the potential on the charges)
      do i=1,natms
        pot(i) = 0.0d0
      end do

      if (lcon.and.lconmol) then
        do m=1,mx_molecules
          molpot(m) = 0.0d0
        end do
      end if
      
c      call timchk(0, t1)
c      write (*,*) "zero pot  ", idnode, t1-t0
c      t0=t1
      
c     this mimicks the double loop (outer: natoms, inner: neigborlist)
c     of the forces_module and the coul routines
c             
      ii=0
      do i=idnode+1,natms,mxnode
        ii=ii+1
        Jij_index = Jij_mapper(i)
c       first all atoms from neigborlist
        do jj=1,lentry(ii)
           j = list(ii,jj)
           term = Jij(Jij_index)
           pot(i) = pot(i) + term*chge(j)
           pot(j) = pot(j) + term*chge(i)
           Jij_index = Jij_index + 1
        end do     
c       then all atoms from exclude list     
        if (loglnk) then
           do jj=1,nexatm(ii)
             j = lexatm(ii,jj)
             if (j.gt.i) then
                term = Jij(Jij_index)
                pot(i) = pot(i) + term*chge(j)
                pot(j) = pot(j) + term*chge(i)
                Jij_index = Jij_index + 1
             end if
           end do
        else
           do jj=1,nexatm(ii)
             j = lexatm(ii,jj)
                term = Jij(Jij_index)
                pot(i) = pot(i) + term*chge(j)
                pot(j) = pot(j) + term*chge(i)
                Jij_index = Jij_index + 1
           end do
        end if
      end do      

c      call timchk(0, t1)
c      write (*,*) "compute pot  ", idnode, t1-t0, ii
c      t0=t1
      
c      now do global sum for parallel runs
      if (mxnode.gt.1) then
        call gdsum(pot, natms, buffer)
      end if

      if (lcon) then
        if (lconmol) then
c         collect total potential per atom
          do i=1,natms
            m = mol_which(i)
            molpot(m) = molpot(m) + pot(i)
          end do
          
c         parallel communicate (because molecules could be astray over nodes)
c          if (mxnode.gt.1) call gdsum(molpot, mx_molecules, buffer)

c         compute average potential for each mol
          do m=1,mx_molecules
            molpot(m) = molpot(m)/mol_natoms(m)
          end do
          
c         now correct qpot ... 
          do i=1,natms
             pot(i) = pot(i) - molpot(mol_which(i))
          end do
      
        else
c         constrain total system        
        end if
      end if
      
c      call timchk(0, t1)
c      write (*,*) "par pot  ", idnode, t1-t0
c      t0=t1

      end subroutine apply_Jij

cRS  specific helper routines for QEq (keep them here in the coulomb module which makes sense)

      subroutine constrain_qpot_permol(idnode,mxnode,qpot, natoms, 
     x                                 molpot, nmols, natpermol)
      
        implicit none
        
        integer,intent(in)    :: idnode
        integer,intent(in)    :: mxnode
        real(8),intent(inout) :: qpot(natoms)
        integer,intent(in)    :: natoms
        real(8),intent(inout) :: molpot(nmols)
        integer,intent(in)    :: nmols
        integer,intent(in)    :: natpermol(nmols)
        
        integer                 :: i, m

        do m=1,nmols
          molpot(m) = 0.0d0
        end do
        
c collect total potential per atom
        do i=idnode+1,natoms,mxnode
          molpot(mol_which(i)) = molpot(mol_which(i)) + qpot(i)
        end do
c parallel communicate (because molecules could be astray over nodes)
        if (mxnode.gt.1) call gdsum(molpot, nmols, buffer)
c compute average potential for each mol
        do m=1,nmols
          molpot(m) = molpot(m)/natpermol(m)
        end do
c now correct qpot ... probably it is faster to do this on all nodes instead of communicating a correction term
        do i=1,natoms
          qpot(i) = qpot(i) - molpot(mol_which(i))
        end do
        
      end subroutine constrain_qpot_permol
      
      
cRS

      subroutine coul0
     x  (lsolva,lfree,lghost,iatm,ik,engcpe,vircpe,rcut,epsq)

c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic force.
c     1/r potential, no truncation or damping
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1993
c     author    - t. forester february 1993
c     stress tensor - t.forester may 1994
c     adapted   - p.-a. cazade oct 2007: solvation, free energy etc
c     
c     wl
c     2008/12/23 10:29:11
c     1.4
c     Exp
c     
c***********************************************************************
      
      implicit none

      logical lsolva,lfree,lghost,lselect,lskip
      integer iatm,ik,m,jatm,kkk
      real(8) engcpe,vircpe,rcut,epsq
      real(8) rcsq,strs1,strs2,strs3,strs5,strs6,strs9,chgea,rsq
      real(8) chgprd,rrr,coul,fcoul,fi,fx,fy,fz
ccs/rs
      real(8) alpha_ij, alr, erfar_or, expar
      real(8) J, dJ, coul_prefac
      real(8) act_pot_iatm,act_pot_jatm,term
      real(8), parameter   :: twooversqrtpi = 2.0d0/sqrpi
crs for core charges
      real(8) alri, alrj, Jci, dJci, Jcj, dJcj, Jcc, dJcc
      real(8) chgprd_i, chgprd_j, chgprd_c     
  
ccs/rs

      dimension fi(3)

CDIR$ CACHE_ALIGN fi
CVAM
CVAM      call VTBEGIN(84, ierr)
CVAM
      lskip=(lfree.or.lghost)
      
c     set cutoff condition for pair forces
      
      rcsq=rcut**2

c     initialise stress tensor accumulators

      strs1=0.d0
      strs2=0.d0
      strs3=0.d0
      strs5=0.d0
      strs6=0.d0
      strs9=0.d0

c     initialise potential energy and virial
      
      engcpe=0.d0
      vircpe=0.d0     
      
c     start of primary loop for forces evaluation

      coul_prefac = r4pie0/epsq  
c      chgea=chge(iatm)*coul_prefac

ccs   Every atoms has to be considered for our purpose even if its charge
ccs   is zero...
c     if(abs(chgea).gt.1.d-10)then
ccs
        fi(1)=fxx(iatm)
        fi(2)=fyy(iatm)
        fi(3)=fzz(iatm)
        
        do m=1,ik

c     atomic index and charge product
          
          jatm=ilist(m)
          
          if(lskip)then
            if(atm_fre(iatm)*atm_fre(jatm).eq.2)cycle
          endif

c          chgprd=chgea*chge(jatm)

ccs   ...because the "potential" (dE/dq(iatm)) could enforce a charge transfer
ccs    to this site during charge iteration
c         if(abs(chgprd).gt.1.d-10)then
ccs

c     calculate interatomic distance
            
            rsq=rsqdf(m)

c     apply truncation of potential
            
            if(rcsq.gt.rsq)then
              
              rrr=sqrt(rsq)
              
c     coulomb potential and force
              
ccs   calculate erf contributions if gaussian charges are added
crs     revised version not using the gaupot function
crs     and revised definition of sigma (in angstrom)
              chgprd = coul_prefac*chge(iatm)*chge(jatm)
              if(lgauss)then
                  alpha_ij=1.0d0/sqrt((qsigma(iatm)*qsigma(iatm))+
     X                                 (qsigma(jatm)*qsigma(jatm)))
                  alr =alpha_ij*rrr
                  erfar_or   = derf(alr)/rrr
                  expar      = twooversqrtpi*alpha_ij
     x                            *exp(-(alr*alr))
                  J=erfar_or
                  dJ=(-erfar_or+expar)/rrr
crs    additional terms if core point charges are included
                  if (use_core_charge) then
                      alri = rrr/qsigma(iatm)
                      alrj = rrr/qsigma(jatm)
                      Jci = derf(alri)/rrr
                      Jcj = derf(alrj)/rrr
                      dJci = (-Jci+(twooversqrtpi/qsigma(iatm)*
     &                      exp(-(alri*alri))))/rrr
                      dJcj = (-Jcj+(twooversqrtpi/qsigma(jatm)*
     &                      exp(-(alrj*alrj))))/rrr
                      Jcc  = 1.0d0/rrr
                      dJcc = -Jcc/rrr                  
                  endif
                  if (qeq_iter) then
                      term = coul_prefac*J
                      act_pot_iatm = chge(jatm)*term
                      act_pot_jatm = chge(iatm)*term
                      if(calc_hess)then
                          q_hess(iatm,jatm) = q_hess(iatm,jatm) + term
                          q_hess(jatm,iatm) = q_hess(jatm,iatm) + term
                      endif
                  endif
                  coul   = chgprd*J 
                  fcoul  =-chgprd*dJ/rrr
                  if (use_core_charge) then
cRS   add core charge terms to coul and fcoul (note: fcoul is the negative grad divided by r)
                      chgprd_i=coul_prefac*chge(iatm)*core_chge(jatm)
                      chgprd_j=coul_prefac*chge(jatm)*core_chge(iatm)
                      chgprd_c=coul_prefac*core_chge(iatm)*
     &                                     core_chge(jatm)
                      coul = coul + chgprd_i*Jci + chgprd_j*Jcj +
     &                              chgprd_c*Jcc
                      fcoul = fcoul - (chgprd_i*dJci+chgprd_j*dJcj+
     &                                 chgprd_c*dJcc)/rrr
                      if (get_core_pot) then
cRS  this should be done once (core_pot does not change during qeq iterations)
                          core_pot(iatm)=core_pot(iatm)+
     &                       coul_prefac*Jci*core_chge(jatm)
                          core_pot(jatm)=core_pot(jatm)+
     &                       coul_prefac*Jcj*core_chge(iatm)
                          enc_coreval=enc_coreval+
     &                        chgprd_i*Jci+chgprd_j*Jcj 
                          enc_core=enc_core+chgprd_c*Jcc
                          enc_val =enc_val +chgprd*J 
cRS DEBUG
c          write (*,*) "DEBUG coulomb0:  ",iatm,jatm,chgprd,J
cRS DEBUG
                      endif                       
                  end if
              else
                  coul=chgprd/rrr
                  fcoul=coul/rsq
              endif
ccs

              
c     set selection control
              
              lselect=.true.
              
c     set double index
              
              if(lsolva)kkk=loc2(atmolt(iatm),atmolt(jatm))
              
              if(lghost)then
                
c     selected excitation option
                
                if((atm_fre(iatm).ne.1).and.(atm_fre(jatm).ne.1))then
                  
c     reset selection control
                  
                  lselect=(atm_fre(iatm)+atm_fre(jatm).eq.0)
                  
c     calculate solvation energy
                  
                  if(lsolva)cou_exc(kkk)=cou_exc(kkk)+coul
                  
                endif
                
              elseif(lfree)then
                
c     selected free energy option
                
                if((atm_fre(iatm).eq.1).or.(atm_fre(jatm).eq.1))then
                  
c     set hamiltonian mixing parameter
                  
                  cou_fre=cou_fre-coul
                  coul=lambda1*coul
                  fcoul=lambda1*fcoul
                  
                elseif((atm_fre(iatm).eq.2).or.(atm_fre(jatm).eq.2))then
                  
c     set hamiltonian mixing parameter

                  cou_fre=cou_fre+coul
                  coul=lambda2*coul
                  fcoul=lambda2*fcoul
                  
                endif
                
              endif
              
              if(lselect)then
                
c     calculate potential energy and virial
              
                engcpe=engcpe+coul
                vircpe=vircpe-coul
ccs
                if(qeq_iter)then
                  pot(iatm)=pot(iatm) + act_pot_iatm
                  pot(jatm)=pot(jatm) + act_pot_jatm
                endif
ccs
c     calculate solvation energy
              
                if(lsolva)cou_sol(kkk)=cou_sol(kkk)+coul
                
c     calculate forces
                
                fx=fcoul*xdf(m)
                fy=fcoul*ydf(m)
                fz=fcoul*zdf(m)
                
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
ccs
c          endif
ccs
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
ccs
c      endif
ccs
CVAM
CVAM      call VTEND(84, ierr)
CVAM
      return
      end subroutine coul0

      subroutine coul1
     x  (lsolva,lfree,lghost,iatm,ik,engcpe,vircpe,rcut,epsq)

c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic forces 
c     assuming a standard coulomb potential truncated at rcut
c     and shifted to zero at rcut.
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith december 1992.
c     adapted   - p.-a. cazade oct 2007: solvation, free energy etc
c     
c     stress tensor t.forester may 1994
c     
c     wl
c     2008/12/23 10:29:11
c     1.4
c     Exp
c     
c***********************************************************************
      
      implicit none

      logical lsolva,lfree,lghost,lselect,lskip
      integer iatm,ik,m,jatm,kkk
      real(8) engcpe,vircpe,rcut,epsq
      real(8) rcsq,strs1,strs2,strs3,strs5,strs6,strs9,chgea,rsq
      real(8) fi,chgprd,omega,egamma,fx,fy,fz,rrr

      dimension fi(3)

CDIR$ CACHE_ALIGN fi
CVAM
CVAM      call VTBEGIN(85, ierr)
CVAM
      lskip=(lfree.or.lghost)

c     set cutoff condition for pair forces
      
      rcsq=rcut**2

c     initialise stress tensor accumulators

      strs1=0.d0
      strs2=0.d0
      strs3=0.d0
      strs5=0.d0
      strs6=0.d0
      strs9=0.d0

c     initialise potential energy and virial
      
      engcpe=0.d0
      vircpe=0.d0
      
      chgea=chge(iatm)*r4pie0/epsq

      if(abs(chgea).gt.1.d-10)then

c     start of primary loop for forces evaluation
        
        fi(1)=fxx(iatm)
        fi(2)=fyy(iatm)
        fi(3)=fzz(iatm)
        
        do m=1,ik
          
c     atomic index and charge product
          
          jatm=ilist(m)
          
          if(lskip)then
            if(atm_fre(iatm)*atm_fre(jatm).eq.2)cycle
          endif
          
          chgprd=chgea*chge(jatm)
          if(abs(chgprd).gt.1.d-10) then

c     calculate interatomic distance
            
            rsq=rsqdf(m)

c     apply truncation of potential
            
            if(rcsq.gt.rsq)then
              
              rrr=sqrt(rsq)

c     calculate potential energy and virial

              omega=chgprd*(rcut-rrr)/(rrr*rcut)
              egamma=chgprd/(rrr*rsq)
              
c     set selection control
              
              lselect=.true.
              
c     set double index
              
              if(lsolva)kkk=loc2(atmolt(iatm),atmolt(jatm))
              
              if(lghost)then
                
c     selected excitation option
                
                if((atm_fre(iatm).ne.1).and.(atm_fre(jatm).ne.1))then
                  
c     reset selection control
                  
                  lselect=(atm_fre(iatm)+atm_fre(jatm).eq.0)
                  
c     calculate solvation energy
                  
                  if(lsolva)cou_exc(kkk)=cou_exc(kkk)+omega
                  
                endif
                
              elseif(lfree)then
                
c     selected free energy option
                
                if((atm_fre(iatm).eq.1).or.(atm_fre(jatm).eq.1))then
                  
c     set hamiltonian mixing parameter
                  
                  cou_fre=cou_fre-omega
                  omega=lambda1*omega
                  egamma=lambda1*egamma
                  
                elseif((atm_fre(iatm).eq.2).or.(atm_fre(jatm).eq.2))then
                  
c     set hamiltonian mixing parameter

                  cou_fre=cou_fre+omega
                  omega=lambda2*omega
                  egamma=lambda2*egamma
                  
                endif
                
              endif
              
              if(lselect)then

c     calculate potential energy and virial
              
                engcpe=engcpe+omega
                vircpe=vircpe-egamma*rsq
                
c     calculate solvation energy
              
                if(lsolva)cou_sol(kkk)=cou_sol(kkk)+omega
              
c     calculate forces
                
                fx=egamma*xdf(m)
                fy=egamma*ydf(m)
                fz=egamma*zdf(m)
                
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

      endif
CVAM
CVAM      call VTEND(85, ierr)
CVAM
      return
      end subroutine coul1

      subroutine coul2
     x  (lsolva,lfree,lghost,iatm,ik,engcpe,vircpe,rcut,epsq)

c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic forces 
c     assuming a distance dependant dielectric `constant'.
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1993
c     author    - t. forester    april 1993
c     stress tensor added - t.forester may 1994
c     adapted   - p.-a. cazade oct 2007: solvation, free energy etc
c     
c     wl
c     2008/12/23 10:29:11
c     1.4
c     Exp
c     
c***********************************************************************
      
      implicit none
      
      logical lsolva,lfree,lghost,lselect,lskip
      integer iatm,ik,m,jatm,kkk
      real(8) engcpe,vircpe,rcut,epsq
      real(8) fi,rcsq,strs1,strs2,strs3,strs5,strs6,strs9,chgea
      real(8) chgprd,rsq,rrsq,coul,egamma,fx,fy,fz
ccs
      real(8) act_gq,act_fgq,act_pot_iatm,act_pot_jatm,term
      integer squared_r
ccs
      
      dimension fi(3)
      
CDIR$ CACHE_ALIGN fi
CVAM
CVAM      call VTBEGIN(86, ierr)
CVAM
      lskip=(lfree.or.lghost)
      
c     set cutoff condition for pair forces
      
      rcsq=rcut**2
      
c     initialise stress tensor accumulators
      
      strs1=0.d0
      strs2=0.d0
      strs3=0.d0
      strs5=0.d0
      strs6=0.d0
      strs9=0.d0
      
c     initialise potential energy and virial
      
      engcpe=0.d0
      vircpe=0.d0

      chgea=chge(iatm)/epsq*r4pie0
ccs   Cancel switch for same reasons as in coul0
c     if(abs(chgea).gt.1.d-10)then
ccs
        fi(1)=fxx(iatm)
        fi(2)=fyy(iatm)
        fi(3)=fzz(iatm)
        
c     start of primary loop for forces evaluation
        
        do m=1,ik
          
c     atomic index and charge product
          
          jatm=ilist(m)
          
          if(lskip)then
            if(atm_fre(iatm)*atm_fre(jatm).eq.2)cycle
          endif
          
          chgprd=chgea*chge(jatm)
ccs   Cancel switch for same reasons as in coul0
c         if(abs(chgprd).gt.1.d-10)then
ccs
            
c     calculate interatomic distance
            
            rsq=rsqdf(m)
            
c     apply truncation of potential
            
            if(rcsq.gt.rsq)then
              
c     calculate potential energy and Virial
              
              rrsq=1.d0/rsq
              coul=chgprd*rrsq
              egamma=2.d0*coul*rrsq
ccs
              if(lgauss)then
                squared_r=1
                call gaupot
c    x          (squared_r,qsigma(iatm),qsigma(jatm),coul,rsq,act_gq)
     x          (squared_r,qsigma(iatm),qsigma(jatm),rsq,act_gq)
                act_pot_iatm=act_gq*chge(jatm)/epsq*r4pie0
                act_pot_jatm=act_gq*chgea
                if(calc_hess)then
                  term = act_gq/epsq*r4pie0
                  q_hess(iatm,jatm) = q_hess(iatm,jatm) + term
                  q_hess(jatm,iatm) = q_hess(jatm,iatm) + term
                endif
                call gauforce
     x          (squared_r,qsigma(iatm),qsigma(jatm),coul,egamma,
     x           rsq,act_fgq)
                coul=act_pot_iatm*chge(iatm)
                egamma=act_fgq
              endif
ccs
              
c     set selection control
              
              lselect=.true.
              
c     set double index
              
              if(lsolva)kkk=loc2(atmolt(iatm),atmolt(jatm))
              if(lghost)then
                
c     selected excitation option
                
                if((atm_fre(iatm).ne.1).and.(atm_fre(jatm).ne.1))then
                  
c     reset selection control
                  
                  lselect=(atm_fre(iatm)+atm_fre(jatm).eq.0)
                  
c     calculate solvation energy
                  
                  if(lsolva)cou_exc(kkk)=cou_exc(kkk)+coul
                  
                endif
                
              elseif(lfree)then
                
c     selected free energy option
                
                if((atm_fre(iatm).eq.1).or.(atm_fre(jatm).eq.1))then
                  
c     set hamiltonian mixing parameter
                  
                  cou_fre=cou_fre-coul
                  coul=lambda1*coul
                  egamma=lambda1*egamma
                  
                elseif((atm_fre(iatm).eq.2).or.(atm_fre(jatm).eq.2))then
                  
c     set hamiltonian mixing parameter
                  
                  cou_fre=cou_fre+coul
                  coul=lambda2*coul
                  egamma=lambda2*egamma
                  
                endif
                
              endif
              
              if(lselect)then

c     calculate potential energy and Virial
                
                engcpe=engcpe+coul
                vircpe=vircpe-2.d0*coul
ccs
                if(qeq_iter)then
                  pot(iatm)=pot(iatm)+act_pot_iatm
                  pot(jatm)=pot(jatm)+act_pot_jatm
                endif
ccs

c     calculate solvation energy
              
                if(lsolva)cou_sol(kkk)=cou_sol(kkk)+coul
              
c     calculate forces
                
                fx=egamma*xdf(m)
                fy=egamma*ydf(m)
                fz=egamma*zdf(m)
                
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
ccs
c          endif
ccs
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
ccs
c      endif
ccs
CVAM
CVAM      call VTEND(86, ierr)
CVAM
      return
      end subroutine coul2

      subroutine coul3
     x  (lsolva,lfree,lghost,iatm,ik,engcpe,vircpe,rcut,epsq)

c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic force.
c     reaction field  potential
c     Ref: M Neumann, J Chem Phys, 82, 5633, (1985)
c     adapted for fennell-gezelter coulombic model
c     by w.smith june 2007
c     Ref: CJ Fennell and JD Gezelter, J Chem Phys, 
c     124, 234104, (2006)
c     adapted   - p.-a. cazade oct 2007: solvation, free energy etc
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1995
c     author    - t. forester february 1995
c     stress tensor - t.forester   feb 1995
c     
c     wl
c     2008/12/23 10:29:11
c     1.4
c     Exp
c     
c***********************************************************************
      
      implicit none

      logical lsolva,lfree,lghost,lselect,lskip
      integer iatm,ik,m,jatm,l,kkk
      real(8) engcpe,vircpe,rcut,epsq,vcon,fcon,rdr,ppp
      real(8) fi,rcsq,b0,rfld0,rfld1,rfld2,strs1,strs2,strs3
      real(8) strs5,strs6,strs9,chgea,chgprd,rsq,coul,omega
      real(8) fx,fy,fz,fcoul,rrr,vk0,vk1,vk2,gk0,gk1,gk2,t1,t2
      dimension fi(3)
CDIR$ CACHE_ALIGN fi
CVAM
CVAM      call VTBEGIN(87, ierr)
CVAM
      lskip=(lfree.or.lghost)

        
c     reaction field terms
      
      b0=2.d0*(epsq-1.d0)/(2.d0*epsq+1.d0)
      rfld0=b0/rcut**3
      rfld1=(1.d0+b0*0.5d0)/rcut
      rfld2=rfld0*0.5d0
      
c     screened coulomb terms
        
      vcon=erc(mxegrd-4)+rfld2*rcut**2-rfld1
      fcon=rcut*fer(mxegrd-4)-rfld0*rcut
      
c     set cutoff condition for pair forces
      
      rcsq=rcut**2
      rdr=dble(mxegrd-4)/rcut

c     initialise stress tensor accumulators

      strs1=0.d0
      strs2=0.d0
      strs3=0.d0
      strs5=0.d0
      strs6=0.d0
      strs9=0.d0

c     initialise potential energy and virial
      
      engcpe=0.d0
      vircpe=0.d0
      
c     start of primary loop for forces evaluation
      
      chgea=chge(iatm)*r4pie0
      
      if(abs(chgea).gt.1.d-10)then
        
        fi(1)=fxx(iatm)
        fi(2)=fyy(iatm)
        fi(3)=fzz(iatm)
        
        do m=1,ik

c     atomic index and charge product
          
          jatm=ilist(m)
          
          if(lskip)then
            if(atm_fre(iatm)*atm_fre(jatm).eq.2)cycle
          endif
          
          chgprd=chgea*chge(jatm)
          if(abs(chgprd).gt.1.d-10)then

c     calculate interatomic distance
            
            rsq=rsqdf(m)

c     apply truncation of potential
            
            if(rcsq.gt.rsq)then
              
              rrr=sqrt(rsq)
              l=int(rrr*rdr)
              ppp=rrr*rdr-dble(l)
              
c     calculate potential energy using 3-point interpolation
              
              vk0=erc(l)
              vk1=erc(l+1)
              vk2=erc(l+2)
              t1=vk0+(vk1-vk0)*ppp
              t2=vk1+(vk2-vk1)*(ppp-1.0d0)
              omega=t1+(t2-t1)*ppp*0.5d0-vcon+fcon*(rrr-rcut)
              coul=chgprd*(omega+rfld2*rsq-rfld1)

c     calculate forces using 3-point interpolation
              
              gk0=fer(l)
              gk1=fer(l+1)
              gk2=fer(l+2)
              t1=gk0+(gk1-gk0)*ppp
              t2=gk1+(gk2-gk1)*(ppp-1.0d0)
              fcoul=chgprd*((t1+(t2-t1)*ppp*0.5d0)-fcon/rrr-rfld0)
              
c     set selection control
              
              lselect=.true.
              
c     set double index
              
              if(lsolva)kkk=loc2(atmolt(iatm),atmolt(jatm))
              
              if(lghost)then
                
c     selected excitation option
                
                if((atm_fre(iatm).ne.1).and.(atm_fre(jatm).ne.1))then
                  
c     reset selection control
                  
                  lselect=(atm_fre(iatm)+atm_fre(jatm).eq.0)
                  
c     calculate solvation energy
                  
                  if(lsolva)cou_exc(kkk)=cou_exc(kkk)+coul
                  
                endif
                
              elseif(lfree)then
                
c     selected free energy option
                
                if((atm_fre(iatm).eq.1).or.(atm_fre(jatm).eq.1))then
                  
c     set hamiltonian mixing parameter
                  
                  cou_fre=cou_fre-coul
                  coul=lambda1*coul
                  fcoul=lambda1*fcoul
                  
                elseif((atm_fre(iatm).eq.2).or.(atm_fre(jatm).eq.2))then
                  
c     set hamiltonian mixing parameter

                  cou_fre=cou_fre+coul
                  coul=lambda2*coul
                  fcoul=lambda2*fcoul
                  
                endif
                
              endif
              
              if(lselect)then

c     calculate coulombic energy and virial
              
                engcpe=engcpe+coul
                vircpe=vircpe-fcoul*rsq
              
c     calculate solvation energy
              
                if(lsolva)cou_sol(kkk)=cou_sol(kkk)+coul
                
c     calculate coulombic force
                
                fx=fcoul*xdf(m)
                fy=fcoul*ydf(m)
                fz=fcoul*zdf(m)
                
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

      endif
CVAM
CVAM      call VTEND(87, ierr)
CVAM
      return
      end subroutine coul3

      subroutine coul4
     X  (lsolva,lfree,lghost,iatm,ik,engcpe,vircpe,rcut,epsq,offs)

c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic forces 
c     assuming a force shifted coulomb potential.
c     adapted for fennell-gezelter coulombic model
c     by w.smith may 2007
c     Ref: CJ Fennell and JD Gezelter, J Chem Phys, 
c     124, 234104, (2006)
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1995
c     author    -  t.forester october  1995
c     adapted   - p.-a. cazade oct 2007: solvation, free energy etc
c     added support for gaussian charge distribs (essentially differnent
c               alpha for each pair) -  R. Schmid, C. Spickermann, 2011, RUB
c     added molecule switching option  - R. Schmid 2011 RUB
c     added temp. storage of Jij integrals - R. Schmid 2011 RUB         
c     
c     wl
c     2008/12/23 10:29:11
c     1.4
c     Exp
c     
c***********************************************************************
      
      implicit none
      
      logical lsolva,lfree,lghost,lselect,lskip
      integer iatm,ik,m,jatm,l,kkk,offs
      real(8) engcpe,vircpe,rcut,epsq,vcon,fcon,rdr,ppp
      real(8) fi,rcsq,strs1,strs2,strs3,strs5,strs6,coul
      real(8) strs9,chgea,chgprd,rsq,rrr,omega,fcoul,fx,fy,fz
      real(8) vk0,vk1,vk2,gk0,gk1,gk2,t1,t2

cRS
      real(8) coul_prefac,alpha_ij,alr,alrc,J,dJ
      real(8) erfar_or,erfarc_orc,exparc_orc,expar_or
      real(8) term, act_pot_iatm,act_pot_jatm
c      real(8), parameter   :: ang2bohr      = 1.0d0/0.529177249d0
      real(8), parameter   :: twooversqrtpi = 2.0d0/sqrpi
cRS switchable molecules
      integer moli,molj
      real(8) lambi,lambj,lambinti
cRS Jij store
      integer Jij_offset

      dimension fi(3)

CDIR$ CACHE_ALIGN fi
CVAM
CVAM      call VTBEGIN(88, ierr)
CVAM
      lskip=(lfree.or.lghost)

      if (store_Jij) then 
          Jij_offset = Jij_mapper(iatm)+offs-1
      end if

c     screened coulomb terms
        
      vcon=erc(mxegrd-4)
      fcon=rcut*fer(mxegrd-4)
      rdr=dble(mxegrd-4)/rcut
      
c     set cutoff condition for pair forces
      
      rcsq=rcut**2

c     initialise stress tensor accumulators

      strs1=0.d0
      strs2=0.d0
      strs3=0.d0
      strs5=0.d0
      strs6=0.d0
      strs9=0.d0

c     initialise potential energy and virial
      
      engcpe=0.d0
      vircpe=0.d0
      
      coul_prefac = r4pie0/epsq
      chgea=chge(iatm)*coul_prefac

cRS skip as in coul0 -> potential even if q=0
c      if(abs(chgea).gt.1.d-10)then

c     start of primary loop for forces evaluation
        
        fi(1)=fxx(iatm)
        fi(2)=fyy(iatm)
        fi(3)=fzz(iatm)
cRS
        if (lmolecules) then
           moli = mol_which(iatm)
           lambi = mol_lamb_coul(moli)
           lambinti = mol_lambint_coul(moli)
        endif

        do m=1,ik
          
c     atomic index and charge product
          
          jatm=ilist(m)
          
          if(lskip)then
            if(atm_fre(iatm)*atm_fre(jatm).eq.2)cycle
          endif
          
          chgprd=chgea*chge(jatm)

cRS skip as in coul0 -> potential even if q=0
c          if(abs(chgprd).gt.1.d-10)then

c     calculate interatomic distance
            
            rsq=rsqdf(m)

c     apply truncation of potential
            
            if(rcsq.gt.rsq)then
              
              rrr=sqrt(rsq)
              
              if (lgauss) then
cRS
c      new version: here we really want to model gaussians! however, each
c                 atom has its own sigma
c
                  alpha_ij=1.0d0/sqrt((qsigma(iatm)*qsigma(iatm))+
     X                                 (qsigma(jatm)*qsigma(jatm)))
              
                  alr =alpha_ij*rrr
                  alrc=alpha_ij*rcut
                  erfar_or   = derf(alr)/rrr
                  erfarc_orc = derf(alrc)/rcut
                  exparc_orc  = twooversqrtpi*alpha_ij
     x                            *exp(-(alrc*alrc))/rcut
                  expar_or    = twooversqrtpi*alpha_ij
     x                            *exp(-(alr*alr))/rrr
                  J=erfar_or - erfarc_orc +
     x                  ((erfarc_orc/rcut) - exparc_orc)*(rrr-rcut)
                  dJ=(-erfar_or/rrr) + expar_or +
     x                   (erfarc_orc/rcut) - exparc_orc
c only potential shifting
c                  J=erfar_or - erfarc_orc 
c                  dJ=(-erfar_or/rrr) + expar_or
c no shifting
c                  J=erfar_or
c                  dJ=(-erfar_or/rrr) + expar_or
cRS Qeq specific stuff
                  if(qeq_iter)then
                      term = coul_prefac*J
                      act_pot_iatm = chge(jatm)*term
                      act_pot_jatm = chge(iatm)*term
                      if(store_Jij) then
                        Jij(Jij_offset+m)=term
                      end if
                      if(calc_hess)then
                        q_hess(iatm,jatm) = q_hess(iatm,jatm) + term
                        q_hess(jatm,iatm) = q_hess(jatm,iatm) + term
                      endif
                  endif
cRS now compute the regular dl_poly coul energy and force/rrr
                  coul  = chgprd*J
                  fcoul = -chgprd*dJ/rrr                   
              else
              
cRS
c      original code using "global" alpha (mimicing point charges)

                  l=int(rrr*rdr)
                  ppp=rrr*rdr-dble(l)
              
c     calculate potential energy using 3-point interpolation
              
                  vk0=erc(l)
                  vk1=erc(l+1)
                  vk2=erc(l+2)
                  t1=vk0+(vk1-vk0)*ppp
                  t2=vk1+(vk2-vk1)*(ppp-1.0d0)
                  omega=t1+(t2-t1)*ppp*0.5d0
                  coul=chgprd*(omega-vcon+fcon*(rrr-rcut))
              
c     calculate forces using 3-point interpolation
              
                  gk0=fer(l)
                  gk1=fer(l+1)
                  gk2=fer(l+2)
                  t1=gk0+(gk1-gk0)*ppp
                  t2=gk1+(gk2-gk1)*(ppp-1.0d0)
                  fcoul=chgprd*((t1+(t2-t1)*ppp*0.5d0)-fcon/rrr)
              
              endif
cRS DEBUG
c          write (*,*) "DEBUG coulomb4:  ",iatm,jatm,rrr,coul/418.4d0
cRS DEBUG
              
c     set selection control
              
              lselect=.true.
              
c     set double index
              
              if(lsolva)kkk=loc2(atmolt(iatm),atmolt(jatm))
              
              if(lghost)then
                
c     selected excitation option
                
                if((atm_fre(iatm).ne.1).and.(atm_fre(jatm).ne.1))then
                  
c     reset selection control
                  
                  lselect=(atm_fre(iatm)+atm_fre(jatm).eq.0)
                  
c     calculate solvation energy
                  
                  if(lsolva)cou_exc(kkk)=cou_exc(kkk)+coul
                  
                endif
                
              elseif(lfree)then
                
c     selected free energy option
                
                if((atm_fre(iatm).eq.1).or.(atm_fre(jatm).eq.1))then
                  
c     set hamiltonian mixing parameter
                  
                  cou_fre=cou_fre-coul
                  coul=lambda1*coul
                  fcoul=lambda1*fcoul
                  
                elseif((atm_fre(iatm).eq.2).or.(atm_fre(jatm).eq.2))then
                  
c     set hamiltonian mixing parameter
                  
                  cou_fre=cou_fre+coul
                  coul=lambda2*coul
                  fcoul=lambda2*fcoul
                  
                endif
cRS
              elseif(lmolecules) then

                  molj=mol_which(jatm)
                  if (molj.ne.moli) then

c     two atoms belong to different molecules -> scale by lambda1*lambda2
c       the force on lambda1 = lambda2*energy and vice versa
 
                   lambj=mol_lamb_coul(molj)
                   mol_dlam_coul(moli)=
     x                        mol_dlam_coul(moli) + lambj*coul
                   mol_dlam_coul(molj)=
     x                        mol_dlam_coul(molj) + lambi*coul
                   coul  = coul*lambi*lambj
                   fcoul = fcoul*lambi*lambj 

                else

c      two atoms belong to the same molecule -> use lambint
c        -> only for switching, no force (not usefull for TD)

                   coul  = coul*lambinti
                   fcoul = fcoul*lambinti

                endif
                
              endif
              
              if(lselect)then

c     calculate the coulombic energy and virial
              
                engcpe=engcpe+coul
                vircpe=vircpe-fcoul*rsq
ccs
                if(qeq_iter)then
                  pot(iatm)=pot(iatm) + act_pot_iatm
                  pot(jatm)=pot(jatm) + act_pot_jatm
                endif
ccs
                
c     calculate solvation energy
                
                if(lsolva)cou_sol(kkk)=cou_sol(kkk)+coul
                
c     calculate coulombic forces
                
                fx=fcoul*xdf(m)
                fy=fcoul*ydf(m)
                fz=fcoul*zdf(m)
                
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
            
c          endif

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

c      endif
CVAM
CVAM      call VTEND(88, ierr)
CVAM
      return
      end subroutine coul4
      
      subroutine coul_nsq
     x  (lsolva,lfree,lghost,idnode,mxnode,natms,imcon,epsq,rcut,
     x  engcpe,vircpe)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic potential and forces
c     for the all-pairs algorithm beyond the range of the normal cutoff
c     i.e. the 'tertiary' forces.  frozen atom option included
c     
c     to be used with multiple_nsq
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory
c     author    - w.smith  august 2008
c     
c     wl
c     2008/12/23 10:29:11
c     1.4
c     Exp
c     
c***********************************************************************
      
      implicit none

      logical lsolva,lfree,lghost,lselect,lskip
      integer natms,idnode,mxnode,imcon,ibig,i,last,mpm2
      integer npm2,m,ii,j,idum,kkk
      real(8) engcpe,epsq,rcut,vircpe,rsq,rrr,chgprd,fcoul,coul,rct2
ccs
      real(8) act_gq,act_fgq,act_pot_i,act_pot_j
      integer squared_r
ccs

CVAM
CVAM      call VTBEGIN(x, ierr)
CVAM
      lskip=(lfree.or.lghost)
      
c     zero energy and virial
      
      engcpe=0.d0
      vircpe=0.d0
      
c     zero force arrays
      
      do i=1,natms
        
        flx(i)=0.d0
        fly(i)=0.d0
        flz(i)=0.d0
        
      enddo
      
c     zero stress tensor
      
      do i=1,9
        stresl(i)=0.d0
      enddo
      
c     zero solvation and  excitation accumulators
      
      if(lsolva)then
        
        cou_sol_lng(:)=0.d0
        
        if(lghost)then
          
          cou_exc_lng(:)=0.d0
          
        endif
        
      endif
      
c     set control variables
      
      last=natms
      mpm2=natms/2
      npm2=(natms-1)/2
      
c     set cutoff radius
      
      rct2=rcut**2
      
c     outer loop over atoms
      
      do m=1,mpm2
        
        if(m.gt.npm2)last=mpm2
        
c     inner loop over atoms
        
        ii=0
        do i=idnode+1,last,mxnode
          
c     calculate atom indices
          
          j=i+m
          if(j.gt.natms)j=j-natms
          
          if(lskip)then
            if(atm_fre(i)*atm_fre(j).eq.2)cycle
          endif
          
          if(lstfrz(i)*lstfrz(j).eq.0)then
            
            ii=ii+1
            
c     calculate interatomic displacements
            
            xdf(ii)=xxx(i)-xxx(j)
            ydf(ii)=yyy(i)-yyy(j)
            zdf(ii)=zzz(i)-zzz(j)
            
          endif
          
        enddo
        
c     apply minimum image convention
        
        call images(imcon,0,1,ii,cell,xdf,ydf,zdf)

c     calculate coulomb terms
        
        ii=0
        
        do i=idnode+1,last,mxnode
          
c     calculate atom indices
          
          j=i+m
          if(j.gt.natms)j=j-natms
          
          if(lskip)then
            if(atm_fre(i)*atm_fre(j).eq.2)cycle
          endif
          
          ii=ii+1
          if(lstfrz(i).eq.0.or.lstfrz(j).eq.0)then
            
c     reject frozen atoms and calculate interatomic distance
            
            rsq=xdf(ii)**2+ydf(ii)**2+zdf(ii)**2
            
c     running check of neighbour list array capacity
            
            if(rsq.ge.rct2)then
              
              chgprd=chge(i)*chge(j)*r4pie0/epsq
              rrr=sqrt(rsq)
              
c     calculate potential energy and force
              
              coul=chgprd/rrr
              fcoul=coul/rsq
ccs   calculate erf contributions if gaussian charges are added
              if(lgauss)then
                squared_r=0
                call gaupot
c    x          (squared_r,qsigma(i),qsigma(j),coul,rrr,act_gq)
     x          (squared_r,qsigma(i),qsigma(j),rrr,act_gq)
c     Beware of the indices here (i is fixed atom)
                act_pot_i=act_gq*chge(j)/epsq*r4pie0
                act_pot_j=act_gq*chge(i)/epsq*r4pie0
                call gauforce
     x           (squared_r,qsigma(i),qsigma(j),coul,fcoul,rrr,act_fgq)
                coul=act_pot_i*chge(i)
                fcoul=act_fgq
              endif
ccs
              
c     set selection control
              
              lselect=.true.
              
c     set double index
              
              if(lsolva)kkk=loc2(atmolt(i),atmolt(j))
              
              if(lghost)then
                
c     selected excitation option
                
                if((atm_fre(i).ne.1).and.(atm_fre(j).ne.1))then
                  
c     reset selection control
                  
                  lselect=(atm_fre(i)+atm_fre(j).eq.0)
                  
c     calculate solvation energy
                  
                  if(lsolva)cou_exc_lng(kkk)=cou_exc_lng(kkk)+coul
                  
                endif
                
              elseif(lfree)then
                
c     selected free energy option
                
                if((atm_fre(i).eq.1).or.(atm_fre(j).eq.1))then
                  
c     set hamiltonian mixing parameter
                  
                  cou_fre=cou_fre-coul
                  coul=lambda1*coul
                  fcoul=lambda1*fcoul
                  
                elseif((atm_fre(i).eq.2).or.(atm_fre(j).eq.2))then
                  
c     set hamiltonian mixing parameter

                  cou_fre=cou_fre+coul
                  coul=lambda2*coul
                  fcoul=lambda2*fcoul
                  
                endif
                
              endif
              
              if(lselect)then
                
c     calculate potential energy and virial
              
                engcpe=engcpe+coul
                vircpe=vircpe-coul

                if(qeq_iter)then
                  pot(i)=pot(i) + act_pot_i
                  pot(j)=pot(j) + act_pot_j
                endif
ccs

c     calculate solvation energy
              
                if(lsolva)cou_sol_lng(kkk)=cou_sol_lng(kkk)+coul
                
c     calculate forces
                
                flx(i)=flx(i)+fcoul*xdf(ii)
                fly(i)=fly(i)+fcoul*ydf(ii)
                flz(i)=flz(i)+fcoul*zdf(ii)             
                
                flx(j)=flx(j)-fcoul*xdf(ii)
                fly(j)=fly(j)-fcoul*ydf(ii)
                flz(j)=flz(j)-fcoul*zdf(ii)  
                
c     stress tensor
                
                stresl(1)=stresl(1)+xdf(ii)*fcoul*xdf(ii)
                stresl(2)=stresl(2)+xdf(ii)*fcoul*ydf(ii)
                stresl(3)=stresl(3)+xdf(ii)*fcoul*zdf(ii)
                stresl(5)=stresl(5)+ydf(ii)*fcoul*ydf(ii)
                stresl(6)=stresl(6)+ydf(ii)*fcoul*zdf(ii)
                stresl(9)=stresl(9)+zdf(ii)*fcoul*zdf(ii)
                
              endif
              
            endif
            
          endif
          
        enddo
          
      enddo

c     complete stress tensor

      stresl(4)=stresl(2)
      stresl(7)=stresl(3)
      stresl(8)=stresl(6)
      
CVAM
CVAM      call VTEND(x,ierr)
CVAM
      return
      end subroutine coul_nsq
      
      end module coulomb_module
