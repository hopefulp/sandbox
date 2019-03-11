      module exclude_module

c***********************************************************************
c     
c     dl_poly module for defining excluded atom arrays
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     
c     wl
c     2008/12/23 10:29:12
c     1.8
c     Exp
c     
c***********************************************************************

      use angles_module
      use bonds_module
      use config_module
      use core_shell_module
      use dihedral_module
      use error_module
      use inversion_module
      use rigid_body_module
      use setup_module
      use shake_module
      use site_module

      implicit none

      integer, allocatable :: lexatm(:,:)
      integer, allocatable :: nexatm(:),noxatm(:)
ccs
      logical, dimension(4) :: mapper
      logical                 :: chklist
      data mapper(:)/.true., .true., .true., .true./,chklist/.true./
c      data mapper(:)/.false.,.false.,.false.,.false./,chklist/.true./
ccs
      save lexatm,nexatm,noxatm
      contains
      
      subroutine alloc_exc_arrays(idnode)

c***********************************************************************
c     
c     dl_poly subroutine for allocating excluded atom arrays
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     
c     wl
c     2008/12/23 10:29:12
c     1.8
c     Exp
c     
c***********************************************************************

      implicit none

      integer, parameter :: nnn=3

      integer i,fail,idnode
      dimension fail(nnn)

      do i=1,nnn
        fail(i)=0
      enddo
      allocate (lexatm(msatms,mxexcl),stat=fail(1))
      allocate (nexatm(msatms),stat=fail(2))
      allocate (noxatm(msatms),stat=fail(3))
      do i=1,nnn
        if(fail(i).gt.0)call error(idnode,1012)
      enddo

      end subroutine alloc_exc_arrays

      subroutine exclude(idnode,mxnode,natms,ntpmls)

c***********************************************************************
c     
c     dl_poly subroutine for constructing the excluded pair
c     interaction list of the system to be simulated
c     
c     keybnd < 0 distance restraint so not excluded
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith        june 1992
c     
c     rigid body exclusions added : t.forester nov 1993
c     check on 1..4 scale factors : t.forester feb 1994
c     inversion terms added       : w.smith    jul 1996
c     
c     wl
c     2008/12/23 10:29:12
c     1.8
c     Exp
c     
c***********************************************************************
      
      implicit none

      logical newjob,check,safe,lchk
      integer idnode,mxnode,natms,ntpmls,nsatms,ibig,i,ibonds
      integer iangle,iconst,idihdr,invers,igrp,isite,ishels,ia,ib,k
      integer ic,id,jz,j,jj,ia1,ib1,kk,jk,ntpsit,nlast,lexsav,itmols
      real(8) a1,a2,a3
ccs
      integer,allocatable::  pairs(:,:)
      logical skip
      real(8) diff
      integer  l,m,n,o,counter,no_pairs
ccs
      save newjob
      
      data newjob/.true./
      
CVAM
CVAM      call VTBEGIN(77, ierr)
CVAM
      
      if(newjob)then

c     check on array allocations
        
        nsatms=(natms+mxnode-1)/mxnode
        if(nsatms.gt.msatms)then
          
          if(idnode.eq.0)write(nrite,*)'make msatms >= ',nsatms
          call error(idnode,100)
          
        endif
        
        newjob=.false.
        
      endif
      
c     variables for array bound checking

      ibig=0
      safe=.true.

c     initialise excluded atom arrays
      
      do i=1,mxsite
        
        nexsit(i)=0
        
      enddo
      
      do i=1,msatms
        
        nexatm(i)=0
        
      enddo

      do j=1,mxexcl
        
        do i=1,mxsite
          
          lexsit(i,j)=0
          
        enddo
        
        do i=1,msatms
          
          lexatm(i,j)=0
          
        enddo
        
      enddo
      

c     loop over molecules in system
      
      ibonds=0
      iangle=0
      iconst=0
      idihdr=0
      invers=0
      igrp  =0
      isite =0
      ishels=0
      
      do itmols=1,ntpmls
        
c     exclude sites on basis of chemical bonds
ccs
        if(mapper(1))then
ccs
        do i=1,numbonds(itmols)
          
          ibonds=ibonds+1
          
          if(keybnd(ibonds).gt.0)then
            
            ia=lstbnd(ibonds,1)+isite
            ib=lstbnd(ibonds,2)+isite

c     check interaction not already included

            lchk=.true.
            do jz=1,min(nexsit(ia),mxexcl)
              if(lexsit(ia,jz).eq.ib-isite)lchk=.false.
            enddo
            if(lchk)then
              nexsit(ia)=nexsit(ia)+1
              nexsit(ib)=nexsit(ib)+1
              if(max(nexsit(ia),nexsit(ib)).gt.mxexcl)then
                ibig=max(ibig,nexsit(ia),nexsit(ib))
                safe=.false.
              else
                lexsit(ia,nexsit(ia))=ib-isite
                lexsit(ib,nexsit(ib))=ia-isite
              endif
            endif
            
          endif
          
        enddo
ccs
        endif
ccs
c     exclude sites on basis of bond constraints
        
        do i=1,numcon(itmols)
          
          iconst=iconst+1
          ia=lstcon(iconst,1)+isite
          ib=lstcon(iconst,2)+isite

c     check interaction not already included

          lchk=.true.
          do jz=1,min(nexsit(ia),mxexcl)
            if(lexsit(ia,jz).eq.ib-isite)lchk=.false.
          enddo
          if(lchk)then
            nexsit(ia)=nexsit(ia)+1
            nexsit(ib)=nexsit(ib)+1
            if(max(nexsit(ia),nexsit(ib)).gt.mxexcl)then
              ibig=max(ibig,nexsit(ia),nexsit(ib))
              safe=.false.
            else
              lexsit(ia,nexsit(ia))=ib-isite
              lexsit(ib,nexsit(ib))=ia-isite
            endif
          endif
          
        enddo
        
c     exclude sites on basis of bond angles
ccs
        if(mapper(2))then
ccs
        do i=1,numang(itmols)
          
          iangle=iangle+1
          if(keyang(iangle).gt.0)then
            ia=lstang(iangle,1)+isite
            ib=lstang(iangle,2)+isite
            ic=lstang(iangle,3)+isite

c     check if already added to lists ..
c     ia - ib interaction

            lchk=.true.
            do jz=1,min(nexsit(ia),mxexcl)
              if(lexsit(ia,jz).eq.ib-isite)lchk=.false.
            enddo
            if(lchk)then
              nexsit(ia)=nexsit(ia)+1
              nexsit(ib)=nexsit(ib)+1
              if(max(nexsit(ia),nexsit(ib)).gt.mxexcl)then
                ibig=max(ibig,nexsit(ia),nexsit(ib))
                safe=.false.
              else
                lexsit(ia,nexsit(ia))=ib-isite
                lexsit(ib,nexsit(ib))=ia-isite
              endif
            endif

c     ib - ic interaction

            lchk=.true.
            do jz=1,min(nexsit(ib),mxexcl)
              if(lexsit(ib,jz).eq.ic-isite)lchk=.false.
            enddo
            if(lchk)then
              nexsit(ib)=nexsit(ib)+1
              nexsit(ic)=nexsit(ic)+1
              if(max(nexsit(ib),nexsit(ic)).gt.mxexcl)then
                ibig=max(ibig,nexsit(ib),nexsit(ic))
                safe=.false.
              else
                lexsit(ib,nexsit(ib))=ic-isite
                lexsit(ic,nexsit(ic))=ib-isite
              endif
            endif

c     ia - ic interaction

            lchk=.true.
            do jz=1,min(nexsit(ia),mxexcl)
              if(lexsit(ia,jz).eq.ic-isite)lchk=.false.
            enddo
            if(lchk)then
              nexsit(ia)=nexsit(ia)+1
              nexsit(ic)=nexsit(ic)+1
              if(max(nexsit(ia),nexsit(ic)).gt.mxexcl)then
                ibig=max(ibig,nexsit(ia),nexsit(ic))
                safe=.false.
              else
                lexsit(ia,nexsit(ia))=ic-isite
                lexsit(ic,nexsit(ic))=ia-isite
              endif
            endif

          endif
          
        enddo
ccs
        endif
ccs
c     exclude on basis of rigid groups
        
        do i=1,numgrp(itmols)
          
          igrp=igrp+1
          id=listyp(igrp)
          
          do jj=1,numgsit(id)-1
            
            ia=lstgst(igrp,jj)+isite
            
            do jk=jj+1,numgsit(id)
              
              ib=lstgst(igrp,jk)+isite
              
c     check interaction not already included

              lchk=.true.
              do jz=1,min(nexsit(ia),mxexcl)
                if(lexsit(ia,jz).eq.ib-isite)lchk=.false.
              enddo
              
              if(lchk)then
                nexsit(ia)=nexsit(ia)+1
                nexsit(ib)=nexsit(ib)+1
                if(max(nexsit(ia),nexsit(ib)).gt.mxexcl)then
                  ibig=max(ibig,nexsit(ia),nexsit(ib))
                  safe=.false.
                else
                  lexsit(ia,nexsit(ia))=ib-isite
                  lexsit(ib,nexsit(ib))=ia-isite
                endif
              endif
              
            enddo
            
          enddo
          
        enddo
        
c     exclude sites on basis of 1-4 dihedral angles
ccs
        if(mapper(3))then
ccs
        do i=1,numdih(itmols)
          
          idihdr=idihdr+1
          ia=lstdih(idihdr,1)+isite
          ib=lstdih(idihdr,2)+isite
          ic=lstdih(idihdr,3)+isite
          id=lstdih(idihdr,4)+isite
          
c     check if already added to lists ..
c     ia - ib interaction
          
          lchk=.true.
          do jz=1,min(nexsit(ia),mxexcl)
            if(lexsit(ia,jz).eq.ib-isite)lchk=.false.
          enddo
          if(lchk)then
            nexsit(ia)=nexsit(ia)+1
            nexsit(ib)=nexsit(ib)+1
            if(max(nexsit(ia),nexsit(ib)).gt.mxexcl)then
              ibig=max(ibig,nexsit(ia),nexsit(ib))
              safe=.false.
            else
              lexsit(ia,nexsit(ia))=ib-isite
              lexsit(ib,nexsit(ib))=ia-isite
            endif
          endif
          
c     ib - ic interaction
          
          lchk=.true.
          do jz=1,min(nexsit(ib),mxexcl)
            if(lexsit(ib,jz).eq.ic-isite)lchk=.false.
          enddo
          if(lchk)then
            nexsit(ib)=nexsit(ib)+1
            nexsit(ic)=nexsit(ic)+1
            if(max(nexsit(ib),nexsit(ic)).gt.mxexcl)then
              ibig=max(ibig,nexsit(ib),nexsit(ic))
              safe=.false.
            else
              lexsit(ib,nexsit(ib))=ic-isite
              lexsit(ic,nexsit(ic))=ib-isite
            endif
          endif
          
c     ia - ic interaction
          
          lchk=.true.
          do jz=1,min(nexsit(ia),mxexcl)
            if(lexsit(ia,jz).eq.ic-isite)lchk=.false.
          enddo
          if(lchk)then
            nexsit(ia)=nexsit(ia)+1
            nexsit(ic)=nexsit(ic)+1
            if(max(nexsit(ia),nexsit(ic)).gt.mxexcl)then
              ibig=max(ibig,nexsit(ia),nexsit(ic))
              safe=.false.
            else
              lexsit(ia,nexsit(ia))=ic-isite
              lexsit(ic,nexsit(ic))=ia-isite
            endif
          endif
          
c     id - ib interaction
          
            lchk=.true.
            do jz=1,min(nexsit(id),mxexcl)
              if(lexsit(id,jz).eq.ib-isite)lchk=.false.
            enddo
            if(lchk)then
              nexsit(id)=nexsit(id)+1
              nexsit(ib)=nexsit(ib)+1
              if(max(nexsit(id),nexsit(ib)).gt.mxexcl)then
                ibig=max(ibig,nexsit(id),nexsit(ib))
                safe=.false.
              else
                lexsit(id,nexsit(id))=ib-isite
                lexsit(ib,nexsit(ib))=id-isite
              endif
            endif
            
c     id - ic interaction
            
            lchk=.true.
            do jz=1,min(nexsit(id),mxexcl)
              if(lexsit(id,jz).eq.ic-isite)lchk=.false.
            enddo
            if(lchk)then
              nexsit(id)=nexsit(id)+1
              nexsit(ic)=nexsit(ic)+1
              if(max(nexsit(id),nexsit(ic)).gt.mxexcl)then
                ibig=max(ibig,nexsit(id),nexsit(ic))
                safe=.false.
              else
                lexsit(id,nexsit(id))=ic-isite
                lexsit(ic,nexsit(ic))=id-isite
              endif
            endif
            
c     ia - id interaction: may need to reset vdw and elec scale factors
            
            lchk=.true.
            do jz=1,min(nexsit(ia),mxexcl)
              if(lexsit(ia,jz).eq.id-isite)lchk=.false.
            enddo
            if(lchk)then
              
              nexsit(ia)=nexsit(ia)+1
              nexsit(id)=nexsit(id)+1
              if(max(nexsit(ia),nexsit(id)).gt.mxexcl)then
                ibig=max(ibig,nexsit(ia),nexsit(id))
                safe=.false.
              else
                lexsit(ia,nexsit(ia))=id-isite
                lexsit(id,nexsit(id))=ia-isite
              endif
              
            elseif(keydih(idihdr).ne.7)then
              
c     if already excluded reset 1..4 vdw and coulombic scale factors
cRS    BUGFIX changed scale params from 4,5 to 5,6 (to allow for cos4)
              
              check=((abs(prmdih(idihdr,5)).gt.1.d-10).or.
     x          (abs(prmdih(idihdr,6)).gt.1.d-10))
            
              if(check)then
                
                a1=dble(itmols)
                a2=dble(ia)
                a3=dble(id)
                call warning(idnode,20,a1,a2,a3)
                
                prmdih(idihdr,5)=0.d0
                prmdih(idihdr,6)=0.d0
                
              endif
              
            endif
            
        enddo
ccs
        endif
ccs
          
c     exclude sites on basis of inversion potentials
ccs
        if(mapper(4))then
ccs
        do i=1,numinv(itmols)
          
          invers=invers+1
          ia=lstinv(invers,1)+isite
          ib=lstinv(invers,2)+isite
          ic=lstinv(invers,3)+isite
          id=lstinv(invers,4)+isite

c     check if already added to lists ..
c     ia - ib interaction

          lchk=.true.
          do jz=1,min(nexsit(ia),mxexcl)
            if(lexsit(ia,jz).eq.ib-isite)lchk=.false.
          enddo
          if(lchk)then
            nexsit(ia)=nexsit(ia)+1
            nexsit(ib)=nexsit(ib)+1
            if(max(nexsit(ia),nexsit(ib)).gt.mxexcl)then
              ibig=max(ibig,nexsit(ia),nexsit(ib))
              safe=.false.
            else
              lexsit(ia,nexsit(ia))=ib-isite
              lexsit(ib,nexsit(ib))=ia-isite
            endif
          endif

c     ib - ic interaction

          lchk=.true.
          do jz=1,min(nexsit(ib),mxexcl)
            if(lexsit(ib,jz).eq.ic-isite)lchk=.false.
          enddo
          if(lchk)then
            nexsit(ib)=nexsit(ib)+1
            nexsit(ic)=nexsit(ic)+1
            if(max(nexsit(ib),nexsit(ic)).gt.mxexcl)then
              ibig=max(ibig,nexsit(ib),nexsit(ic))
              safe=.false.
            else
              lexsit(ib,nexsit(ib))=ic-isite
              lexsit(ic,nexsit(ic))=ib-isite
            endif
          endif

c     ia - ic interaction

          lchk=.true.
          do jz=1,min(nexsit(ia),mxexcl)
            if(lexsit(ia,jz).eq.ic-isite)lchk=.false.
          enddo
          if(lchk)then
            nexsit(ia)=nexsit(ia)+1
            nexsit(ic)=nexsit(ic)+1
            if(max(nexsit(ia),nexsit(ic)).gt.mxexcl)then
              ibig=max(ibig,nexsit(ia),nexsit(ic))
              safe=.false.
            else
              lexsit(ia,nexsit(ia))=ic-isite
              lexsit(ic,nexsit(ic))=ia-isite
            endif
          endif

c     id - ib interaction

            lchk=.true.
            do jz=1,min(nexsit(id),mxexcl)
              if(lexsit(id,jz).eq.ib-isite)lchk=.false.
            enddo
            if(lchk)then
              nexsit(id)=nexsit(id)+1
              nexsit(ib)=nexsit(ib)+1
              if(max(nexsit(id),nexsit(ib)).gt.mxexcl)then
                ibig=max(ibig,nexsit(id),nexsit(ib))
                safe=.false.
              else
                lexsit(id,nexsit(id))=ib-isite
                lexsit(ib,nexsit(ib))=id-isite
              endif
            endif

c     id - ic interaction

            lchk=.true.
            do jz=1,min(nexsit(id),mxexcl)
              if(lexsit(id,jz).eq.ic-isite)lchk=.false.
            enddo
            if(lchk)then
              nexsit(id)=nexsit(id)+1
              nexsit(ic)=nexsit(ic)+1
              if(max(nexsit(id),nexsit(ic)).gt.mxexcl)then
                ibig=max(ibig,nexsit(id),nexsit(ic))
                safe=.false.
              else
                lexsit(id,nexsit(id))=ic-isite
                lexsit(ic,nexsit(ic))=id-isite
              endif
            endif

c     ia - id interaction

            lchk=.true.
            do jz=1,min(nexsit(ia),mxexcl)
              if(lexsit(ia,jz).eq.id-isite)lchk=.false.
            enddo
            if(lchk)then
              nexsit(ia)=nexsit(ia)+1
              nexsit(id)=nexsit(id)+1
              if(max(nexsit(ia),nexsit(id)).gt.mxexcl)then
                ibig=max(ibig,nexsit(ia),nexsit(id))
                safe=.false.
              else
                lexsit(ia,nexsit(ia))=id-isite
                lexsit(id,nexsit(id))=ia-isite
              endif

          endif
          
        enddo
ccs
        endif
ccs
        
c     exclude sites on basis of core-shell units
        
        do i=1,numshl(itmols)
          
          ishels=ishels+1
          
          ia=lstshl(ishels,1)+isite
          ib=lstshl(ishels,2)+isite

c     check interaction not already included
          
          lchk=.true.
          do jz=1,min(nexsit(ia),mxexcl)
            if(lexsit(ia,jz).eq.ib-isite)lchk=.false.
          enddo
          if(lchk)then
            nexsit(ia)=nexsit(ia)+1
            nexsit(ib)=nexsit(ib)+1
            if(max(nexsit(ia),nexsit(ib)).gt.mxexcl)then
              ibig=max(ibig,nexsit(ia),nexsit(ib))
              safe=.false.
            else
              lexsit(ia,nexsit(ia))=ib-isite
              lexsit(ib,nexsit(ib))=ia-isite
            endif
          endif
          
c     exclude sites on basis of bonds to core-shell units

          ibonds=ibonds-numbonds(itmols)
          do kk=1,numbonds(itmols)
            
            ibonds=ibonds+1
            
            if(keybnd(ibonds).gt.0)then
              
              ia1=lstbnd(ibonds,1)+isite
              ib1=lstbnd(ibonds,2)+isite

              if(ia.eq.ia1)then

c     check interaction not already included
          
                lchk=.true.
                do jz=1,min(nexsit(ib1),mxexcl)
                  if(lexsit(ib1,jz).eq.ib-isite)lchk=.false.
                enddo
                if(lchk)then
                  nexsit(ib1)=nexsit(ib1)+1
                  nexsit(ib)=nexsit(ib)+1
                  if(max(nexsit(ib1),nexsit(ib)).gt.mxexcl)then
                    ibig=max(ibig,nexsit(ib1),nexsit(ib))
                    safe=.false.
                  else
                    lexsit(ib1,nexsit(ib1))=ib-isite
                    lexsit(ib,nexsit(ib))=ib1-isite
                  endif
                endif

              endif

              if(ia.eq.ib1)then

c     check interaction not already included
          
                lchk=.true.
                do jz=1,min(nexsit(ia1),mxexcl)
                  if(lexsit(ia1,jz).eq.ib-isite)lchk=.false.
                enddo
                if(lchk)then
                  nexsit(ia1)=nexsit(ia1)+1
                  nexsit(ib)=nexsit(ib)+1
                  if(max(nexsit(ia1),nexsit(ib)).gt.mxexcl)then
                    ibig=max(ibig,nexsit(ia1),nexsit(ib))
                    safe=.false.
                  else
                    lexsit(ia1,nexsit(ia1))=ib-isite
                    lexsit(ib,nexsit(ib))=ia1-isite
                  endif
                endif

              endif

              if(ib.eq.ia1)then

c     check interaction not already included
          
                lchk=.true.
                do jz=1,min(nexsit(ia),mxexcl)
                  if(lexsit(ia,jz).eq.ib1-isite)lchk=.false.
                enddo
                if(lchk)then
                  nexsit(ia)=nexsit(ia)+1
                  nexsit(ib1)=nexsit(ib1)+1
                  if(max(nexsit(ia),nexsit(ib1)).gt.mxexcl)then
                    ibig=max(ibig,nexsit(ia),nexsit(ib1))
                    safe=.false.
                  else
                    lexsit(ia,nexsit(ia))=ib1-isite
                    lexsit(ib1,nexsit(ib1))=ia-isite
                  endif
                endif

              endif
              if(ib.eq.ib1)then

c     check interaction not already included
          
                lchk=.true.
                do jz=1,min(nexsit(ia),mxexcl)
                  if(lexsit(ia,jz).eq.ia1-isite)lchk=.false.
                enddo
                if(lchk)then
                  nexsit(ia)=nexsit(ia)+1
                  nexsit(ia1)=nexsit(ia1)+1
                  if(max(nexsit(ia),nexsit(ia1)).gt.mxexcl)then
                    ibig=max(ibig,nexsit(ia),nexsit(ia1))
                    safe=.false.
                  else
                    lexsit(ia,nexsit(ia))=ia1-isite
                    lexsit(ia1,nexsit(ia1))=ia-isite
                  endif
                endif

              endif

            endif

          enddo
          
c     exclude sites on basis of constraint bonds to core-shell units
          
          iconst=iconst-numcon(itmols)
          do kk=1,numcon(itmols)
            
            iconst=iconst+1
            
            ia1=lstcon(iconst,1)+isite
            ib1=lstcon(iconst,2)+isite

            if(ia.eq.ia1)then

c     check interaction not already included
          
              lchk=.true.
              do jz=1,min(nexsit(ib1),mxexcl)
                if(lexsit(ib1,jz).eq.ib-isite)lchk=.false.
              enddo
              if(lchk)then
                nexsit(ib1)=nexsit(ib1)+1
                nexsit(ib)=nexsit(ib)+1
                if(max(nexsit(ib1),nexsit(ib)).gt.mxexcl)then
                  ibig=max(ibig,nexsit(ib1),nexsit(ib))
                  safe=.false.
                else
                  lexsit(ib1,nexsit(ib1))=ib-isite
                  lexsit(ib,nexsit(ib))=ib1-isite
                endif
              endif

            endif

            if(ia.eq.ib1)then

c     check interaction not already included
          
              lchk=.true.
              do jz=1,min(nexsit(ia1),mxexcl)
                if(lexsit(ia1,jz).eq.ib-isite)lchk=.false.
              enddo
              if(lchk)then
                nexsit(ia1)=nexsit(ia1)+1
                nexsit(ib)=nexsit(ib)+1
                if(max(nexsit(ia1),nexsit(ib)).gt.mxexcl)then
                  ibig=max(ibig,nexsit(ia1),nexsit(ib))
                  safe=.false.
                else
                  lexsit(ia1,nexsit(ia1))=ib-isite
                  lexsit(ib,nexsit(ib))=ia1-isite
                endif
              endif

            endif

            if(ib.eq.ia1)then

c     check interaction not already included
          
              lchk=.true.
              do jz=1,min(nexsit(ia),mxexcl)
                if(lexsit(ia,jz).eq.ib1-isite)lchk=.false.
              enddo
              if(lchk)then
                nexsit(ia)=nexsit(ia)+1
                nexsit(ib1)=nexsit(ib1)+1
                if(max(nexsit(ia),nexsit(ib1)).gt.mxexcl)then
                  ibig=max(ibig,nexsit(ia),nexsit(ib1))
                  safe=.false.
                else
                  lexsit(ia,nexsit(ia))=ib1-isite
                  lexsit(ib1,nexsit(ib1))=ia-isite
                endif
              endif

            endif
            if(ib.eq.ib1)then

c     check interaction not already included
          
              lchk=.true.
              do jz=1,min(nexsit(ia),mxexcl)
                if(lexsit(ia,jz).eq.ia1-isite)lchk=.false.
              enddo
              if(lchk)then
                nexsit(ia)=nexsit(ia)+1
                nexsit(ia1)=nexsit(ia1)+1
                if(max(nexsit(ia),nexsit(ia1)).gt.mxexcl)then
                  ibig=max(ibig,nexsit(ia),nexsit(ia1))
                  safe=.false.
                else
                  lexsit(ia,nexsit(ia))=ia1-isite
                  lexsit(ia1,nexsit(ia1))=ia-isite
                endif
              endif

            endif

          enddo
          
c     exclude sites on basis of rigid units involving  core or shell
          
          igrp=igrp-numgrp(itmols)
          do kk=1,numgrp(itmols)
          
            igrp=igrp+1
          
            id=listyp(igrp)
          
            do jj=1,numgsit(id)
            
              ia1=lstgst(igrp,jj)+isite
              if(ia1.eq.ia)then

                do jk=1,numgsit(id)
            
                  if(jk.ne.jj)then
                    ib1=lstgst(igrp,jk)+isite

c     check interaction not already included
          
                    lchk=.true.
                    do jz=1,min(nexsit(ib1),mxexcl)
                      if(lexsit(ib1,jz).eq.ib-isite)lchk=.false.
                    enddo
                    if(lchk)then
                      nexsit(ib1)=nexsit(ib1)+1
                      nexsit(ib)=nexsit(ib)+1
                      if(max(nexsit(ib1),nexsit(ib)).gt.mxexcl)then
                        ibig=max(ibig,nexsit(ib1),nexsit(ib))
                        safe=.false.
                      else
                        lexsit(ib1,nexsit(ib1))=ib-isite
                        lexsit(ib,nexsit(ib))=ib1-isite
                      endif
                    endif

                  endif

                enddo
                
              endif

              if(ia1.eq.ib)then

                do jk=1,numgsit(id)
            
                  if(jk.ne.jj)then
                    ib1=lstgst(igrp,jk)+isite

c     check interaction not already included
          
                    lchk=.true.
                    do jz=1,min(nexsit(ia),mxexcl)
                      if(lexsit(ia,jz).eq.ib1-isite)lchk=.false.
                    enddo
                    if(lchk)then
                      nexsit(ia)=nexsit(ia)+1
                      nexsit(ib1)=nexsit(ib1)+1
                      if(max(nexsit(ia),nexsit(ib1)).gt.mxexcl)then
                        ibig=max(ibig,nexsit(ia),nexsit(ib1))
                        safe=.false.
                      else
                        lexsit(ia,nexsit(ia))=ib1-isite
                        lexsit(ib1,nexsit(ib1))=ia-isite
                      endif

                    endif

                  endif

                enddo
                
              endif

            enddo

          enddo

        enddo

        isite=isite+numsit(itmols)
        
      enddo
      
      ntpsit=isite

c     check for exceeded array bounds
      
      if(mxnode.gt.1)call gstate(safe)
      if(.not.safe)then
        if(mxnode.gt.1)call gimax(ibig,1,jj)
        if(idnode.eq.0)write(nrite,*)'mxexcl must be at least ',ibig
        if(idnode.eq.0)write(nrite,*)'mxexcl is currently     ',mxexcl
        call error(idnode,65)
      endif

c     remove redundant entries from exclusion list
c     (there shouldn't be any!)
      
      do i=1,ntpsit
        
        nlast=nexsit(i)
        do j=1,nexsit(i)-1
          
          if(j.lt.nlast)then
            
            kk=j
            do k=j+1,nexsit(i)
              
              if(lexsit(i,j).eq.lexsit(i,k))then
                
                nlast=nlast-1
                lexsit(i,k)=0
                
              else if(lexsit(i,k).gt.0)then
                
                kk=kk+1
                lexsav=lexsit(i,k)
                lexsit(i,k)=0
                lexsit(i,kk)=lexsav
                
              endif
              
            enddo
            
          endif
          
        enddo
        
        nexsit(i)=nlast
        
      enddo
ccs   *****************Code for printing the excluded atom list********
c      no_pairs=mxatms*(mxatms-1)/2
c      allocate(pairs(no_pairs,2))
c      pairs(:,:)=-1
c
c      counter=0
c      do l=1,mxsite
c        do m=1,mxexcl
c          skip=.false.
c          ia=l
c          ib=lexsit(l,m)
c          do n=1,no_pairs
c            if(ia==pairs(n,1))then
c              if(ib==pairs(n,2)) skip=.true.
c            elseif(ib==pairs(n,1))then
c              if(ia==pairs(n,2)) skip=.true.
c            elseif(ib==0)then
c              skip=.true.
c            endif
c          enddo
c          if(.not.skip)then
c            counter=counter+1
c            pairs(counter,1)=ia
c            pairs(counter,2)=ib
c            diff=(xxx(ia)-xxx(ib))**2+(yyy(ia)-yyy(ib))**2
c     x       +(zzz(ia)-zzz(ib))**2
c            diff=sqrt(diff)
c            write(*,'(a24,i4,a10,i4,a12,f8.2)')
c     x       'Excluded list: index 1: ',ia,
c     x       ' index 2: ',ib,'  distance: ',diff
c          endif
c        enddo
c      enddo
ccs   *****************************************************************

CVAM
CVAM      call VTEND(77, ierr)
CVAM
      return
      end subroutine exclude

      subroutine excludeneu(idnode,mxnode,nneut)

c***********************************************************************
c     
c     dl_poly subroutine for constructing the excluded pair
c     interaction list of the system to be simulated
c     part 2 - neutral group implementation
c     
c     copyright - daresbury                1994
c     author    - t. forester        march 1994
c     
c     wl
c     2008/12/23 10:29:12
c     1.8
c     Exp
c     
c***********************************************************************
      
      implicit none

      logical lchk
      integer idnode,mxnode,nneut,ibig,iatom,jatom,last,mpm2
      integer npm2,m,ii,im,itmols,inoff,isoff,isit,iolsit,jm,jtmols
      integer jnoff,jsoff,jsit,jolsit,jn1,jno1,jsite,jsite0,in1,ino1
      integer jj0,isite,ij,idum,it
      
CVAM
CVAM      call VTBEGIN(80, ierr)
CVAM

c     construct excluded pair list for verlet neighbour correction

      ibig=0
      iatom=0
      jatom=0

c     generate all atomic pairs and check for exclusions 
c     with Brode Ahlrichs ordering of groups
      
      last=nneut
      lchk=.true.
      mpm2=nneut/2+1
      npm2=(nneut-1)/2+1
      
c     outer loop over groups
      
      do m=1,mpm2
        
        if(m.gt.npm2)last=mpm2-1
        
c     inner loop over groups - include intragroup interactions
        
        ii=0
        
        do im=idnode+1,last,mxnode
          
          ii=ii+1

c     first site in neutral group
          
          itmols=1
          inoff=0
          isoff=0
          isit=numsit(itmols)*nummols(itmols)
          iolsit=numsit(itmols)
          
c     calculate j group indices
          
          jm=im+m-1
          if(jm.gt.nneut)jm=jm-nneut
          
c     inner loop over neutral groups
          
          jtmols=1
          jnoff=0
          jsoff=0
          jsit=numsit(jtmols)*nummols(jtmols)
          jolsit=numsit(jtmols)
          
c     test first sites in neutral group
          
          jatom=neulst(jm)        

c     establish pointer to sets
          
          do while(jatom.gt.jsit)
            
            jtmols=jtmols+1
            jnoff=jsit
            jsoff=jsoff+jolsit
            jsit=jsit+nummols(jtmols)*numsit(jtmols)
            jolsit=numsit(jtmols)
            
          enddo
          
          jn1=jatom-jnoff
          jno1=(jn1/jolsit)*jolsit
          jsite=jn1-jno1
          if(jsite.eq.0)then 
            jsite=jolsit
            jno1=jno1-jolsit
          endif
          jsite=jsite+jsoff
          jsite0=jsite-1
          
          do iatom=neulst(im),neulst(im+1)-1
            
c     establish pointer to sets
            
            do while(iatom.gt.isit)
              
              itmols=itmols+1
              inoff=isit
              isoff=isoff+iolsit
              isit=isit+nummols(itmols)*numsit(itmols)
              iolsit=numsit(itmols)
              
            enddo
            
            in1=iatom-inoff
            ino1=(in1/iolsit)*iolsit
            isite=in1-ino1
            if(isite.eq.0)then 
              isite=iolsit
              ino1=ino1-iolsit
            endif
            isite=isite+isoff
            
c     test im and jm are neutral groups on same molecule
            
            if((jnoff.eq.inoff).and.(ino1.eq.jno1))then
            if(abs(im-jm).lt.iolsit)then
              
              jj0=neulst(jm)
              jsite=jsite0
              
c     special case for im=jm (ie. same group)
              
              if(im.eq.jm)then 
                
                jj0=iatom+1
                jsite=isite
                
              endif

c     test for excluded interaction
              
              do jatom=jj0,neulst(jm+1)-1
                
                jsite=jsite+1
                
                do ij=1,nexsit(isite)
                  
                  if(lexsit(isite,ij).eq.jsite-jsoff)then
                    
                    it=nexatm(ii)
                    
                    if(it+2.gt.mxexcl)then
                      
                      ibig=max(it+2,ibig)
                      nexatm(ii)=it+2
                      lchk=.false.
                      
                    else
                      
                      lexatm(ii,it+1)=iatom
                      lexatm(ii,it+2)=jatom
                      nexatm(ii)=nexatm(ii)+2
                      
                    endif

                  endif
                  
                enddo
                
              enddo
              
            endif
            endif
            
          enddo
          
        enddo
        
      enddo

c     global check
      
      call gstate(lchk)
      if(.not.lchk)then
        
        if(mxnode.gt.1)call gimax(ibig,1,idum)
        if(idnode.eq.0)write(nrite,*)'mxexcl must be at least ',ibig
        if(idnode.eq.0)write(nrite,*)'mxexcl is currently     ',mxexcl
        call error(idnode,260)

      endif
CVAM
CVAM      call VTEND(80, ierr)
CVAM
      return
      end subroutine excludeneu

      subroutine exclude_link(idnode,mxnode,ntpmls)

c***********************************************************************
c     
c     dl_poly subroutine for constructing the excluded pair
c     interaction list of the system to be simulated
c     
c     part 2 - link cell implementation
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith        june 1992
c     
c     wl
c     2008/12/23 10:29:12
c     1.8
c     Exp
c     
c***********************************************************************
      
      implicit none

      integer idnode,mxnode,ntpmls,iatom,jatom,lsite,ksite
      integer itmols,imols,isite,kk,newatm,k

c     construct excluded pair list for verlet neighbour correction
      
CVAM
CVAM      call VTBEGIN(79, ierr)
CVAM

      iatom=0
      jatom=0
      lsite=0
      ksite=0
      
      do itmols=1,ntpmls
        
        do imols=1,nummols(itmols)
          
          do isite=1,numsit(itmols)
            
            iatom=iatom+1
            
            if(mod(iatom-1,mxnode).eq.idnode)then
              
              kk=0
              jatom=jatom+1
              
              do k=1,nexsit(ksite+isite)
                
                newatm=lexsit(ksite+isite,k)+lsite
                
                kk=kk+1
                lexatm(jatom,kk)=newatm
                
              enddo
              
              nexatm(jatom)=kk
              
            endif
            
          enddo
          
          lsite=lsite+numsit(itmols)
          
        enddo
        
        ksite=ksite+numsit(itmols)
        
      enddo
CVAM
CVAM      call VTEND(79, ierr)
CVAM
      return
      end subroutine exclude_link

      subroutine exclude_atom(idnode,mxnode,natms,ntpmls)

c***********************************************************************
c     
c     dl_poly subroutine for constructing the excluded pair
c     interaction list of the system to be simulated
c     part 2 
c
c     copyright - daresbury laboratory 1992
c     author    - w. smith        june 1992
c     
c     wl
c     2008/12/23 10:29:12
c     1.8
c     Exp
c     
c***********************************************************************
            
      implicit none

      integer idnode,mxnode,natms,ntpmls,iatom,jatom,lsite
      integer ksite,itmols,isite,imols,k,kk,newatm,j,latom,i,ii

CVAM
CVAM      call VTBEGIN(78, ierr)
CVAM

c     construct excluded pair list for verlet neighbour correction
      
      iatom=0
      jatom=0
      lsite=0
      ksite=0
      
      do itmols=1,ntpmls
        
        do imols=1,nummols(itmols)
          
          do isite=1,numsit(itmols)
            
            iatom=iatom+1
            
            if(mod(iatom-1,mxnode).eq.idnode)then
              
              kk=0
              jatom=jatom+1
              
              do k=1,nexsit(ksite+isite)
                
                newatm=lexsit(ksite+isite,k)+lsite

c     keep only brode-ahlrichs combinations of indices
                
                if(((newatm.gt.iatom).and.
     x            (newatm-iatom.le.natms/2)).or.
     x            ((newatm.lt.iatom).and.
     x            (newatm+natms-iatom.le.(natms-1)/2)))then
                  
                  kk=kk+1
                  lexatm(jatom,kk)=newatm
                  
                  if(kk.gt.1)then
                    
c     sort the excluded atom list in ascending indices

                    do j=kk,2,-1
                      
                      if(lexatm(jatom,j).lt.lexatm(jatom,j-1))
     x                  then
                        latom=lexatm(jatom,j)
                        lexatm(jatom,j)=lexatm(jatom,j-1)
                        lexatm(jatom,j-1)=latom
                      endif
                      
                    enddo
                    
                  endif
                  
                endif
                
              enddo
              
              nexatm(jatom)=kk
              
            endif
            
          enddo
          
          lsite=lsite+numsit(itmols)
          
        enddo
        
        ksite=ksite+numsit(itmols)
        
      enddo
      
c     final sort into brode-ahlrichs ordering
      
      ii=0
      do i=1+idnode,natms,mxnode
        
        ii=ii+1
        do j=1,nexatm(ii)
          
          if(lexatm(ii,1).lt.i)then
            
            latom=lexatm(ii,1)
            
            do k=1,nexatm(ii)-1
              
              lexatm(ii,k)=lexatm(ii,k+1)
              
            enddo
            
            lexatm(ii,nexatm(ii))=latom
            
          endif
          
        enddo
        
      enddo
CVAM
CVAM      call VTEND(78, ierr)
CVAM
      return
      end subroutine exclude_atom
      
      end module exclude_module

