      module molecule_module

c***********************************************************************
c
c      module for handling molecules or subsystems
c      it allows to switch on and off individual molecules
c      this is part of the pydlpoly project
c      
c      This module only works with pydlpoly, that means instead of  
c      controling via CONTROL file flags it needs to be initalized directly
c      via calls wrapped to python
c
c      The current version of molecules works only together with the 
c      damped dhifted force coulomb implementation "coul4"
c
c      (C) R. Schmid 2011 , RUB
c
c         - Aug 2012 : added an internal lambda (switch on and off interactions
c                      within a molecule (no force!) -> use for QM/MM
c         - May 2013 : added arrays/functions to handle kinetic energy and thermostats
c                      on a per-molecule basis
c                      NOTE: add in the future degrees of freedom per mol and a chit for a thermostat per molecule
c         - Jun 2015 : added spltting vdw lambdas (one lambda for dispersion and one for repulsion)
c
c***********************************************************************

      use setup_module
      use error_module
      use config_module

      implicit none

      logical   lmolecules
      logical   lmol_ekin
      logical   lmol_lrcorr
      integer   mx_molecules

      integer, allocatable :: mol_which(:)
      integer, allocatable :: mol_natoms(:)
      real(8), allocatable :: mol_lamb_coul(:)
      real(8), allocatable :: mol_dlam_coul(:)
      real(8), allocatable :: mol_lamb_vdw(:)
      real(8), allocatable :: mol_dlam_vdw(:)
      real(8), allocatable :: mol_lamb_vdwr(:)
      real(8), allocatable :: mol_dlam_vdwr(:)
      real(8), allocatable :: mol_lamb_ekin(:)
      real(8), allocatable :: mol_lambint_coul(:)
      real(8), allocatable :: mol_lambint_vdw(:)

cRS   per molecule dipole moments 
      real(8), allocatable :: mol_dip(:,:)    !! Note: the shape is 3,  nmolecules ... to work with python
      logical lcalc_moldip
      
      real(8), allocatable :: mol_ekin(:)
      
      real(8) :: mol_dlam_lrcorr            !! this is the force on the sinlge non-zero lambda
                                          !! note that the lrcorection force only works for a sinlge black-sheep with 0<lambda<1
      integer, allocatable :: mol_numtyp_bs(:)
      
      save lmolecules, lmol_ekin, lmol_lrcorr, mx_molecules
      save mol_lamb_coul, mol_dlam_coul, mol_lamb_vdw, mol_dlam_vdw
      save mol_lamb_vdwr, mol_dlam_vdwr
      save mol_lamb_ekin
      save mol_lambint_coul, mol_lambint_vdw
      save mol_which
      save mol_ekin
      save mol_dlam_lrcorr
      save mol_numtyp_bs

      data lcalc_moldip/.false./
      save lcalc_moldip, mol_dip
      
      contains
      
      subroutine init_molecules(nmolecules)
c        
c     initialize molecule module:
c      - set the number of molecules and allocate
c      - initalize all arrays (all atoms belong to mol 1 and lambdas all 1.)
c     in order to set to anything reasonable use the following methods
            
        implicit none
      
        integer nmolecules
        integer i, fail, allfail
        
        lmolecules = .true.
        lmol_ekin  = .false.
        lmol_lrcorr= .false.
        mx_molecules = nmolecules
        
        allfail = 0
        
        allocate (mol_which(mxatms),stat=fail)
        allfail=allfail+fail
        allocate (mol_natoms(mx_molecules),stat=fail)
        allfail=allfail+fail

        allocate (mol_lamb_coul(mx_molecules),stat=fail)
        allfail=allfail+fail
        allocate (mol_dlam_coul(mx_molecules),stat=fail)
        allfail=allfail+fail
        allocate (mol_lamb_vdw(mx_molecules),stat=fail)
        allfail=allfail+fail
        allocate (mol_dlam_vdw(mx_molecules),stat=fail)
        allfail=allfail+fail
        allocate (mol_lamb_vdwr(mx_molecules),stat=fail)
        allfail=allfail+fail
        allocate (mol_dlam_vdwr(mx_molecules),stat=fail)
        allfail=allfail+fail

cRS     internal lambdas
        allocate (mol_lambint_coul(mx_molecules), stat= fail)
        allfail= allfail+fail
        allocate (mol_lambint_vdw(mx_molecules), stat= fail)
        allfail= allfail+fail
        
cRS     molecular kinetic energy
        allocate (mol_lamb_ekin(mx_molecules),stat=fail)
        allfail=allfail+fail
        allocate (mol_ekin(mx_molecules), stat=fail)
        allfail = allfail+fail
        
cRS     number of atomtypes in BS for long range correction
        allocate (mol_numtyp_bs(mxatyp), stat=fail)
        allfail = allfail+fail
        
        if (allfail.gt.0) then
          write (*,*) "ERROR during allocation in molecule_module"
          lmolecules = .false.
          mx_molecules = 0
          return
        end if
        
c       init all atoms to belong to molecule 1 and all lambdas to 1.0
        do i=1,mxatms
          mol_which(i)=1
        end do        
        do i=1,mx_molecules
          mol_lamb_coul(i) = 1.0d0
          mol_dlam_coul(i) = 0.0d0
          mol_lamb_vdw(i)  = 1.0d0
          mol_dlam_vdw(i)  = 0.0d0
          mol_lamb_vdwr(i)  = 1.0d0
          mol_dlam_vdwr(i)  = 0.0d0
          mol_lamb_ekin(i) = 1.0d0
          mol_lambint_coul(i) = 1.0d0
          mol_lambint_vdw(i)  = 1.0d0
          
          mol_ekin(i) = 0.0d0
        end do
        
        return    
        
      end subroutine init_molecules
      
      subroutine molecules_scale_ekin(scale)
      
         implicit none
         
         real*8  ::  scale, scale2
         integer ::  i
         
         scale2 = scale**2
         do i=0,mx_molecules
            mol_ekin(i)=mol_ekin(i)*scale2
         enddo
     
      end subroutine molecules_scale_ekin   
            
      subroutine molecules_gsum_dlam()
      
        implicit none
        
        call gdsum(mol_dlam_coul,mx_molecules,buffer)
        call gdsum(mol_dlam_vdw,mx_molecules,buffer)
        call gdsum(mol_dlam_vdwr,mx_molecules,buffer)
        
        return
        
      end subroutine molecules_gsum_dlam
      
      subroutine molecules_calc_moldip

        implicit none

        integer fail

        lcalc_moldip = .true.
        allocate(mol_dip(3,mx_molecules), stat=fail)

        if (fail.gt.0) then
          write (*,*) "ERROR durin allocation of moldip array"
        endif 

      end subroutine molecules_calc_moldip

      end module 
