      module site_module

c***********************************************************************
c     
c     dl_poly module for defining atomic/site arrays
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     
c     wl
c     2008/12/23 10:29:12
c     1.4
c     Exp
c     
c***********************************************************************

      use error_module
      use parse_module
      use setup_module

crs
      use config_module, only : use_core_charge

cRS   NOTE: we should try to deallocate all the site arrays because in pydlpoly atoms are directly 
c           read in 8no repeating of atoms in FIELD
c     so with sites and actual atoms all is duplicated!!
cRS
      
      implicit none

      character*1, allocatable :: molnam(:,:)
      character*8, allocatable :: sitnam(:),unqatm(:)
ccs   sgmsit added in analogy to chgsit
      real(8), allocatable :: dens(:),chgsit(:),wgtsit(:),sgmsit(:)
ccs
crs
      real(8), allocatable :: coresit(:)

c     arrays for sp-basis (we do not need to read in pcharges because they make no sense for rotating molecules anyway
      real(8), allocatable :: psgmsit(:)
      
crs    
      integer, allocatable :: nexsit(:),lfzsit(:),numsit(:),ltpsit(:)
      integer, allocatable :: nugrp(:),lexsit(:,:),numgrp(:)
      integer, allocatable :: numtyp(:),numfrz(:),nummols(:)

      save numtyp,numfrz,dens,chgsit,wgtsit,sitnam,unqatm,nexsit
      save lfzsit,numsit,ltpsit,nugrp,lexsit,numgrp,molnam,nummols

      save sgmsit
      save coresit
      save psgmsit
      
      contains

      subroutine alloc_site_arrays(idnode)

      implicit none

      integer, parameter :: nnn=16

      integer i,fail,idnode
      dimension fail(nnn)

      do i=1,nnn
         fail(i)=0
      enddo
      allocate (chgsit(mxsite),stat=fail(1))
      allocate (wgtsit(mxsite),stat=fail(2))
      allocate (nexsit(mxsite),stat=fail(3))
      allocate (lfzsit(mxsite),stat=fail(4))
      allocate (nugrp(mxsite) ,stat=fail(5))
      allocate (ltpsit(mxsite),stat=fail(6))
      allocate (numsit(mxtmls),stat=fail(7))
      allocate (lexsit(mxsite,mxexcl),stat=fail(8))
      allocate (sitnam(mxsite),stat=fail(9))
      allocate (unqatm(mxsite),stat=fail(10))
      allocate (numgrp(mxtmls),stat=fail(11))
      allocate (numtyp(mxatyp),stat=fail(12))
      allocate (numfrz(mxatyp),stat=fail(13))
      allocate (dens(mxatyp),stat=fail(14))
      allocate (nummols(mxtmls),stat=fail(15))
      allocate (molnam(40,mxtmls),stat=fail(16))
ccs   allocate array for sigma of sites
      if(lgauss)then
      allocate (sgmsit(mxsite),stat=fail(1))
c      write (*,*) "DEBUG: sgmsit allocated"
        if(use_core_charge)then
          allocate (coresit(mxsite),stat=fail(1))
c         write (*,*) "DEBUG DEBUG allocating core charge", mxsite
        endif
        if(lspbasis)then
          allocate (psgmsit(mxsite),stat=fail(1))
c          write (*,*) "DEBUG: psgmsit allocated"
        endif
      endif
ccs

      do i=1,nnn
         if(fail(i).ne.0)call error(idnode,1090)
      enddo

      do i=1,mxtmls
         numsit(i)=0
      enddo

      end subroutine alloc_site_arrays

      subroutine define_atoms
     x  (safe,lneut,idnode,itmols,nsite,ksite,ntpatm)

c***********************************************************************
c     
c     dl_poly subroutine for  defining atom types in system
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     
c     wl
c     2008/12/23 10:29:12
c     1.4
c     Exp
c     
cRS NOTE: in pydlpoly OUTPUT is obselete and only used for DEBUG output
cRS       therefore we skipt any complex formatting of the site charge/sigma output
c
c***********************************************************************

      implicit none

      character*8 atom1
      character*1 message(80)
ccs   lgauss added
      logical lneut,safe,atmchk
ccs
      integer idnode,itmols,nsite,ksite,ntpatm,isite,nrept
      integer ifrz,neugp,irept,jsite,idum
      real(8) weight,charge,isigma,corecg
      real(8) psigma

      numsit(itmols)=intstr(record,lenrec,idum)
      if(idnode.eq.0) then
        write(nrite,"(/,1x,'number of atoms/sites',
     x    10x,i10)") numsit(itmols)
     
c        write (*,*) "DEBUG DEBUG"
c        write (*,*) lgauss, lspbasis,use_core_charge

ccs   provide appropriate header in OUTPUT

       if(.not.lneut) then
        if(.not.lgauss) write(nrite,"(/,/,1x,'atomic characteristics:',
     x    /,/,21x,' site',5x,'name',10x,'mass',8x,
     x    'charge',4x,'repeat',4x,'freeze'/)")
        if(lgauss) then
          if (use_core_charge) then
          write(nrite,"(/,/,1x,'atomic characteristics:',
     x    /,/,21x,' site',5x,'name',10x,'mass',8x,
     x    'charge',4x,' sigma',4x,'corecg',4x,'repeat',4x,'freeze'/)")
          else
            if(lspbasis)then
            write(nrite,"(/,/,1x,'atomic characteristics:',
     x    /,/,21x,' site',5x,'name',10x,'mass',8x,
     x    'charge',4x,' sigma',4x,'repeat',4x,'freeze'/)")          
            else
            write(nrite,"(/,/,1x,'atomic characteristics:',
     x    /,/,21x,' site',5x,'name',10x,'mass',8x,
     x    'charge',4x,' sigma',4x,'psigma',4x,'repeat',4x,
     x    'freeze'/)")
            end if
          end if 
        end if

       elseif(lneut) then

        if(lgauss) then
          write(nrite,"(/,/,1x,'QEq cannot be combined with neutral ',
     x    'groups atm...QEq switched off'/)")

          lgauss=.false.

          write(nrite,"(/,/,1x,'atomic characteristics:',
     x    /,/,21x,' site',5x,'name',10x,'mass',8x,'charge',
     x    4x,'repeat',4x,'freeze',3x,'chg grp')")
        endif
       endif
ccs
        
      endif
      
      do isite=1,numsit(itmols)
        
        if(ksite.lt.numsit(itmols))then

c     read atom name, site number, mass, charge, freeze option
          
          call getrec(safe,idnode,nfield)
          if(.not.safe)return

ccs   parse sigma values

          if(lgauss)then
            call copystring(record,message,80)
            call getword(atom1,record,8,lenrec)
            weight=dblstr(record,lenrec,idum)
            charge=dblstr(record,lenrec,idum)
            isigma=dblstr(record,lenrec,idum)
            if (use_core_charge) corecg=dblstr(record,lenrec,idum)
            if(lspbasis)then
              psigma=dblstr(record,lenrec,idum)
            endif
            nrept=intstr(record,lenrec,idum)
            ifrz =intstr(record,lenrec,idum)
            neugp=intstr(record,lenrec,idum)
            if(nrept.eq.0)nrept=1
            ksite=ksite+nrept
          else
            call copystring(record,message,80)
            call getword(atom1,record,8,lenrec)
            weight=dblstr(record,lenrec,idum)
            charge=dblstr(record,lenrec,idum)
            nrept=intstr(record,lenrec,idum)
            ifrz =intstr(record,lenrec,idum)
            neugp=intstr(record,lenrec,idum)
            if(nrept.eq.0)nrept=1
            ksite=ksite+nrept
          endif
ccs
          if(idnode.eq.0) then
            
            if(.not.lneut) then
ccs
              if(.not.lgauss)then
                write(nrite,
     x            "(21x,i5,5x,a8,2f12.5,2i10)")
     x            nsite+1,atom1,weight,charge,nrept,
     x            ifrz
              else
                if (use_core_charge) then
                write(nrite,
     x            "(21x,i5,5x,a8,4f12.5,2i10)")
     x            nsite+1,atom1,weight,charge,isigma,corecg,nrept,
     x            ifrz
                else
                  if (lspbasis) then
                    write(nrite,
     x              "(21x,i5,5x,a8,4f12.5,2i10)")
     x              nsite+1,atom1,weight,charge,isigma,
     x              psigma,
     x              nrept,ifrz                
                  else
                    write(nrite,
     x              "(21x,i5,5x,a8,3f12.5,2i10)")
     x              nsite+1,atom1,weight,charge,isigma,nrept,
     x              ifrz
                  end if
                end if
              endif
ccs

            else

              write(nrite,
     x          "(21x,i5,5x,a8,2f12.5,3i10)")
     x          nsite+1,atom1,weight,charge,nrept,
     x          ifrz,neugp

            endif

          endif
          
          do irept=1,nrept
            
            nsite=nsite+1
            if(nsite.gt.mxsite) call error(idnode,20)
            
            sitnam(nsite)=atom1
            wgtsit(nsite)=weight
            chgsit(nsite)=charge
ccs   copy sigma to site array
            sgmsit(nsite)=isigma
            if (use_core_charge) coresit(nsite)=corecg
ccs
crs   copy p-basis to site array
            if (lspbasis) then
               psgmsit(nsite)   = psigma
            endif
            lfzsit(nsite)=ifrz
            nugrp(nsite)=neugp
            
          enddo

c     establish list of unique atom types
          
          atmchk=.true.
          
          do jsite=1,ntpatm
            
            if(atom1.eq.unqatm(jsite)) then
              
              atmchk=.false.
              do irept=nsite,nsite-nrept+1,-1
                
                ltpsit(irept)=jsite

              enddo
              
            endif
            
          enddo

          if(atmchk)then
            
            ntpatm=ntpatm+1

            if(ntpatm.gt.mxatyp)call error(idnode,14)
            unqatm(ntpatm)=atom1
            
            do irept=nsite,nsite-nrept+1,-1
              
              ltpsit(irept)=ntpatm

            enddo
            
          endif
          
        endif
        
      enddo
      
      return
      end subroutine define_atoms

      subroutine check_syschg(idnode,ntpmls,sumchg)

c***********************************************************************
c     
c     dl_poly subroutine for checking the system charge
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     
c     wl
c     2008/12/23 10:29:12
c     1.4
c     Exp
c     
cRS  NOTE: in case of an SP-basis charge in the P does not count!
c***********************************************************************

      implicit none

      integer idnode,ntpmls,jsite,itmols,lsite
      real(8) sumchg

      jsite=0
      do itmols=1,ntpmls
        
        do lsite=1,numsit(itmols)
          
          jsite=jsite+1
          sumchg=sumchg+dble(nummols(itmols))*chgsit(jsite)
          if(use_core_charge) then
            sumchg=sumchg+dble(nummols(itmols))*coresit(jsite)
          end if
          
        enddo
        
      enddo
      
      if(abs(sumchg).gt.1.0d-6) then
        
        call warning(idnode,60,sumchg,0.d0,0.d0)
        
      endif
      
      return
      end subroutine check_syschg
      
      end module site_module
