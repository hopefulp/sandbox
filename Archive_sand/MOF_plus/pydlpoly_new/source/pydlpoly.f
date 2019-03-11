           
c***********************************************************************
c
c     PYDLPOLY
c
c     pydlpoly is a Python wrapped and extended fork of DL_Poly
c     (Version 2 now referred to as DL_Poly Classic)
c
c     Authors of the PYDLPOLY extensions:
c     Rochus Schmid, Chrisitan Spickermann, Sareeya Bureekaew, Mohammad Alaghemandi
c     Copyright Ruhr-Universität Bochum 2010-2014
c
c     Original Copyrigh for the DL_Poly basis: 
c     ##################################################################
c     dl_poly classic is an stfc/ccp5 program package for the 
c     dynamical simulation of molecular systems.
c     
c     dl_poly is the copyright of the stfc daresbury laboratory, 
c     daresbury, warrington wa4 4ad. 
c     
c     neither the stfc, daresbury laboratory, ccp5 nor the authors
c     of this package claim that it is free from errors and do not
c     accept liability for any loss or damage that may arise from
c     its use. it is the users responsibility to verify that the 
c     package dl_poly is fit for the purpose the user intends for
c     it.
c     
c     users of this package are recommended to consult the dl_poly
c     user manual for the full description of its use and purpose.
c     
c     authors: w.smith and t.r.forester 1995
c     copyright daresbury laboratory 1995
c     ####################################################################
c
c***********************************************************************
      
c     declare required modules

      module pydlpoly
      
      use angles_module
      use bonds_module
      use config_module
      use core_shell_module
      use define_system_module
      use dihedral_module
      use driver_module
      use ewald_module
      use exclude_module
      use external_field_module
      use forces_module
      use four_body_module
      use hkewald_module
      use hyper_dynamics_module
      use integrator_module
      use inversion_module
      use metal_module
      use nlist_builders_module
      use pair_module
      use pmf_module
      use property_module
      use rigid_body_module
      use setup_module
      use shake_module
      use site_module
      use solvation_module
      use spme_module
      use temp_scalers_module
      use tersoff_module
      use tether_module
      use three_body_module
      use utility_module
      use vdw_module
      
      implicit none
      
      logical ltscal,lzeql,loptim,ltraj,lgofr,lpgr,lfcap,recycle
      logical newlst,lneut,loglnk,lnsq,lzden,lshmov,lcnb,ltad,lneb
      logical stropt,lzero,nolink,newgau,lminim,lminnow,lhit,lbpd
      logical prechk,tadall,lexcite,lsolva,lfree,lfrmas,lswitch
      logical lghost,llswitch
ccs
      logical qeq
ccs

      save ltscal,lzeql,loptim,ltraj,lgofr,lpgr,lfcap,recycle
      save newlst,lneut,loglnk,lnsq,lzden,lshmov,lcnb,ltad,lneb
      save stropt,lzero,nolink,newgau,lminim,lminnow,lhit,lbpd
      save prechk,tadall,lexcite,lsolva,lfree,lfrmas,lswitch
      save lghost,llswitch
      
      integer npage,lines,idnode,mxnode,memr,intsta,istraj,nsbzdn
      integer keyens,keyfce,keyres,keytrj,kmax1,kmax2,kmax3,multt
      integer nstack,nstbgr,nstbpo,nhko,nlatt,nstbts,nsteql,nstraj
      integer nstrun,nospl,keyfld,natms,ngrp,ntpatm,ntpmls,ntpvdw
      integer ntptbp,ntpmet,ntpfbp,nshels,imcon,levcfg,nneut,minstp
      integer ntangl,ntbond,ntcons,ntdihd,ntinv,ntpmf,nspmf,ntfree
      integer ntteth,ntshl,nstep,numacc,numrdf,nzden,nscons,i,k
      integer ntpter,keyshl,isw,keyver,keystr,keytol,numgau,khit
      integer nhit,keybpd,ntrack,nblock,blkout,numneb,nturn,mode
      integer natms2,ntghost,nsolva,isolva

      save npage,lines,idnode,mxnode,memr,intsta,istraj,nsbzdn
      save keyens,keyfce,keyres,keytrj,kmax1,kmax2,kmax3,multt
      save nstack,nstbgr,nstbpo,nhko,nlatt,nstbts,nsteql,nstraj
      save nstrun,nospl,keyfld,natms,ngrp,ntpatm,ntpmls,ntpvdw
      save ntptbp,ntpmet,ntpfbp,nshels,imcon,levcfg,nneut,minstp
      save ntangl,ntbond,ntcons,ntdihd,ntinv,ntpmf,nspmf,ntfree
      save ntteth,ntshl,nstep,numacc,numrdf,nzden,nscons,i,k
      save ntpter,keyshl,isw,keyver,keystr,keytol,numgau,khit
      save nhit,keybpd,ntrack,nblock,blkout,numneb,nturn,mode
      save natms2,ntghost,nsolva,isolva

      
      real(8) alpha,delr,epsq,fmax,press,quattol,rcut,rprim,rvdw,taup
      real(8) taut,temp,timcls,timjob,tolnce,tstep,tzero,dlrpot,drewd
      real(8) engunit,rcuttb,rctter,rcutfb,degfre,degrot,chit,conint
      real(8) elrc,virlrc,engbnd,volm,degshl,chip,virbnd,engang,virang
      real(8) engdih,virdih,enginv,virinv,engtbp,virtbp,engter,virter
      real(8) engfbp,virfbp,engsrp,virsrp,engcpe,vircpe,vircon,vircom
      real(8) engfld,virfld,engshl,virshl,shlke,engtet,virtet,virpmf
      real(8) consv,engke,engrot,sigma,virtot,engcfg
      real(8) stpeng,stpeth,stpprs,stptmp,stpvir,stpvol,width,zlen
      real(8) timelp,engmet,virmet,pass0,pass1,pass2,rlxtol,opttol
      real(8) catchrad,sprneb,deltad,tlow,engtke,ehit,xhit,yhit,zhit
      real(8) ebias,vmin,boost,heinc,tboost,hyp_units,estar
cRS added variables only defined in property_module (obsolete)
      real(8) stpcns,stprot,stpshl

      save alpha,delr,epsq,fmax,press,quattol,rcut,rprim,rvdw,taup
      save taut,temp,timcls,timjob,tolnce,tstep,tzero,dlrpot,drewd
      save engunit,rcuttb,rctter,rcutfb,degfre,degrot,chit,conint
      save elrc,virlrc,engbnd,volm,degshl,chip,virbnd,engang,virang
      save engdih,virdih,enginv,virinv,engtbp,virtbp,engter,virter
      save engfbp,virfbp,engsrp,virsrp,engcpe,vircpe,vircon,vircom
      save engfld,virfld,engshl,virshl,shlke,engtet,virtet,virpmf
      save consv,engke,engrot,sigma,virtot,engcfg
      save stpeng,stpeth,stpprs,stptmp,stpvir,stpvol,width,zlen
      save timelp,engmet,virmet,pass0,pass1,pass2,rlxtol,opttol
      save catchrad,sprneb,deltad,tlow,engtke,ehit,xhit,yhit,zhit
      save ebias,vmin,boost,heinc,tboost,hyp_units,estar
      save stpcns,stprot,stpshl

      
      real(8), save, allocatable :: tbuffer(:)
ccs
cf2py intent(callback) run_qeq
c      external run_qeq
ccs
      
      data timelp/0.d0/,lminnow/.false./,ntrack/10/
      data npage,lines/8,0/,recycle/.true./,boost/1.d0/
      data pass0/0.d0/,pass1/0.d0/,pass2/0.d0/
      data delr,epsq,press,quattol,rprim,rvdw/6*0.d0/
      data temp,timcls,timjob,tolnce,rlxtol/5*0.d0/
ccs
      data qeq/.false./
ccs
CVAM
CVAM      include 'VT.inc'
CVAM

      contains

c ########################################################
      subroutine init_dlpoly()

      implicit none

c     set up the communications
      
c        call initcomms()
        call gsync()
      
c     set up VAMPIR
      
CVAM      call VTSETUP()
CVAM      call VTTRACEON(ierr)
CVAM      call VTBEGIN(99, ierr)
      
c     determine processor identities
      
        call machine(idnode,mxnode)
      
      return
      end subroutine init_dlpoly

c     activate for limited-life executable
      
CBOMB      call bomb(idnode,2008,6,30)
      
c ########################################################
      subroutine allocate_dlpoly()

      implicit none

cRS moved this opening of OUTPUT to allocate since in init we are not yet 
cRS in the right directory!!!

c     open main printing file
      
        nrite = 60
      
        if(idnode.eq.0)open(nrite,file='OUTPUT')
ccs
        if(idnode.eq.0.and.time)open(tout,file='Timings.out')
ccs
        if(idnode.eq.0) write (nrite,
     x  "(/,20x,'DL_POLY Version 2.20',
     x  /,/,30x,'Running on ',i4,' nodes',/,/)") mxnode


        allocate (tbuffer(10),stat=memr)
      
        call parset(idnode,mxnode,tbuffer)
      
c     allocate arrays for each function
      
        call alloc_ang_arrays(idnode)
        call alloc_bnd_arrays(idnode)
        call alloc_config_arrays(idnode)
        call alloc_csh_arrays(idnode)
        call alloc_dih_arrays(idnode)
        call alloc_ewald_arrays(idnode)
        call alloc_exc_arrays(idnode)
        call alloc_exi_arrays(idnode)
        call alloc_fbp_arrays(idnode)
        call alloc_fld_arrays(idnode)
        call alloc_free_arrays(idnode)
        call alloc_hke_arrays(idnode)
        call alloc_hyper_arrays(idnode)
        call alloc_inv_arrays(idnode)
        call alloc_met_arrays(idnode)
        call alloc_pair_arrays(idnode)
        call alloc_pmf_arrays(idnode)
        call alloc_prp_arrays(idnode)
        call alloc_rgbdy_arrays(idnode)
        call alloc_shake_arrays(idnode)
        call alloc_site_arrays(idnode)
        call alloc_sol_arrays(idnode)
        call alloc_spme_arrays(idnode)
        call alloc_tbp_arrays(idnode)
        call alloc_ter_arrays(idnode)
        call alloc_tet_arrays(idnode)
        call alloc_vdw_arrays(idnode)
      return
      end subroutine allocate_dlpoly

c ########################################################
      subroutine setup_dlpoly()

      implicit none

c     start clock
      
      call timchk(0,tzero)
      
c     input the control parameters defining the simulation
      
      call simdef
     x  (lfcap,lgofr,lnsq,loptim,lzero,lminim,lpgr,ltraj,ltscal,lzeql,
     x  lzden,nolink,newgau,lhit,lbpd,ltad,lneb,prechk,tadall,lsolva,
     x  lfree,lfrmas,lexcite,lswitch,lghost,idnode,minstp,intsta,istraj,
     x  keybpd,keyens,keyfce,keyres,keyver,keytrj,kmax1,kmax2,kmax3,
     x  multt,nstack,nstbgr,nsbzdn,nstbpo,nhko,nlatt,nstbts,nsteql,
     x  nstraj,nstrun,nospl,keytol,numgau,khit,nhit,nblock,ntrack,
     x  blkout,numneb,mode,nsolva,isolva,alpha,delr,epsq,fmax,press,
     x  quattol,rcut,rprim,rvdw,taup,taut,temp,timcls,timjob,tolnce,
     x  tstep,rlxtol,opttol,zlen,ehit,xhit,yhit,zhit,ebias,vmin,heinc,
     x  catchrad,sprneb,deltad,tlow,hyp_units)
      
c     input the system force field
      
      call sysdef
     x  (lneut,lnsq,lsolva,lfree,lexcite,lswitch,lghost,idnode,keyfce,
     x  keyfld,natms,ngrp,ntpatm,ntpmls,ntpvdw,ntptbp,ntpmet,ntpfbp,
     x  ntpter,nshels,keyshl,ntghost,dlrpot,engunit,rvdw,rcuttb,rctter,
     x  rcutfb)
      
      if(ntpmet.gt.0.and.multt.gt.1)call error(idnode,153)
      
c     construct initial configuration of system
      
      call sysgen
     x  (loglnk,lneut,nolink,lfree,lfrmas,idnode,imcon,keyens,
     x  keyfce,keyres,levcfg,multt,mxnode,ntpmls,delr,rcut,volm)
      
c     construct initial bookkeeping arrays
      
      call sysbook
     x  (loglnk,lneut,lshmov,lcnb,lsolva,lghost,idnode,imcon,
     x  mxnode,natms,nneut,ngrp,nscons,ntangl,ntbond,ntcons,
     x  ntdihd,ntinv,ntpmls,ntpmf,nspmf,ntfree,ntteth,ntshl,
     x  ntghost,degfre,degrot)
      
c     reset atom numbers for excitation simulation
      
      if(lghost)then
        natms2=natms-ntghost
      else
        natms2=natms
      endif
ccs
        if(time)then
          call timchk(1,tval)
          write(tout,'(a32,f15.3)')'Job time after setup_dlpoly: ',
     x    tval
        endif
ccs
      end subroutine setup_dlpoly

      subroutine reload_field()

      implicit none

c     input the system force field
      
      call sysdef
     x  (lneut,lnsq,lsolva,lfree,lexcite,lswitch,lghost,idnode,keyfce,
     x  keyfld,natms,ngrp,ntpatm,ntpmls,ntpvdw,ntptbp,ntpmet,ntpfbp,
     x  ntpter,nshels,keyshl,ntghost,dlrpot,engunit,rvdw,rcuttb,rctter,
     x  rcutfb)
c     construct initial bookkeeping arrays
      
      call sysbook
     x  (loglnk,lneut,lshmov,lcnb,lsolva,lghost,idnode,imcon,
     x  mxnode,natms,nneut,ngrp,nscons,ntangl,ntbond,ntcons,
     x  ntdihd,ntinv,ntpmls,ntpmf,nspmf,ntfree,ntteth,ntshl,
     x  ntghost,degfre,degrot)
      
      end subroutine reload_field

c#########################################################

      subroutine init_dyn(skip_REVOLD)
      
      implicit none

      logical, intent(in) :: skip_REVOLD

c     set initial system temperature
      
      call systemp
     x  (idnode,imcon,keyres,mxnode,natms2,ngrp,nscons,ntcons,
     x  ntfree,ntshl,levcfg,keyshl,degfre,degshl,degrot,temp,
     x  tolnce)
      
c     read thermodynamic and structural data from restart file
      
      call sysinit
     x  (lgofr,lzden,lsolva,lfree,lghost,idnode,imcon,keyfce,
     x  keyres,mxnode,natms,ntshl,nstep,numacc,numrdf,ntpatm,
     x  ntpmet,ntpvdw,nzden,chip,chit,conint,elrc,engunit,virlrc,
     x  rvdw,volm,virtot,vircom,tboost, skip_REVOLD)
      
c     synchronise LRC, SIC and system charge terms for switching
      
      llswitch=.false.
      if(lswitch)then
        
        if(nstep.ge.nswitch)then
          
          if(mod((nstep-nswitch)/niswitch,2).eq.0)then
            
            call switch_atm(lfrmas)
            call switch(elrc,virlrc)
            llswitch=.true.
            
          endif
          
        endif
        
      endif
      
c     zero long range component of stress
      
      do i=1,9
        stresl(i)=0.d0
      enddo
      
c     zero contraint terms
      
      vircon=0.d0
      virpmf=0.d0
      if(lminim.or.loptim.or.ntcons.eq.0)then

        do i=1,9
          strcns(i)=0.d0
        enddo

      endif

      end subroutine init_dyn
c##############################################################
c    
c    original dl_poly dynamics routine (put back in on July 2012 to test problems with Thermostat/barostat
c    of the new "splitted" and python driven MD loop
c
c    since qeq is meanwhile full in python and this is only for testing the qeq callback by Chris is removed here


      subroutine dynamics()

      implicit none
c     time check

      call timchk(1,tzero)

c     control variable for structure optimizer
      
      keystr=0
      stropt=.false.
      
c     define target kinetic energy
cRS   put this here. it looks like this is already done in systemp
cRS   where the initial velocities are randomized
cRS   However, with this setup it should be possible to change the target
cRS   temp for each call of dynamics
      
      sigma=temp*boltz*degfre*0.5d0

ccs
      if(idnode.eq.0) open(unit=999,file='energy.out',action='write')
ccs

c     calculate initial conditions for velocity verlet
      
      if(keyver.eq.1.and.nstep.eq.0)then
        
c     kinetic stress tensor at start
        
        call dcell(cell,celprp)
        width=min(celprp(7),celprp(8),celprp(9))
        call kinstress(natms,idnode,mxnode,stress)
        engke=0.5d0*(stress(1)+stress(5)+stress(9))
        do i=1,9
          stress(i)=stress(i)/dble(mxnode)
        enddo
        
c     calculate initial forces
ccs
crs        if(qeq)then
crs          call run_qeq()
crs        endif
ccs
        call molecular_dynamics
     x    (lfcap,lgofr,lneut,lnsq,loglnk,loptim,lzeql,lzero,
     x    newlst,stropt,recycle,ltad,lsolva,lfree,lghost,
     x    idnode,imcon,keyfce,keyfld,keyshl,keystr,keytol,kmax1,
     x    kmax2,kmax3,multt,mxnode,natms,ngrp,nhko,nlatt,nneut,
     x    nospl,nscons,nstbgr,nstep,nsteql,ntangl,ntbond,ntdihd,
     x    ntfree,ntinv,ntpfbp,ntpmet,ntptbp,ntpter,ntpvdw,ntshl,
     x    ntpatm,
     x    ntteth,ntcons,numrdf,nsolva,isolva,alpha,delr,dlrpot,
     x    drewd,elrc,engang,engbnd,engcpe,engdih,engfbp,engfld,
     x    enginv,engshl,engsrp,engtbp,engter,engtet,epsq,fmax,
     x    opttol,rctter,rcut,rcutfb,rcuttb,rprim,rvdw,shlke,
     x    engcfg,temp,tstep,virang,virbnd,vircpe,virdih,
     x    virfbp,virfld,virinv,virlrc,virmet,virshl,virsrp,
     x    virtbp,virter,virtet,volm,engmet,virtot)
        
c     bias potential dynamics option - reset forces
        
        if(lbpd)call bpd_forces
     x    (natms,vmin,ebias,temp,engcfg,boost)
        
      endif
      
      if(ltad.or.(lbpd.and.keybpd.eq.2))then
        
c     construct the first reference state

        call hyper_start
     x    (lbpd,lfcap,lneut,lnsq,loglnk,lzeql,newlst,nblock,idnode,
     x    imcon,keyfce,keyfld,keyshl,keytol,kmax1,kmax2,kmax3,multt,
     x    mxnode,natms,ngrp,nhko,nlatt,nneut,nospl,nscons,nstbgr,
     x    nstep,nsteql,ntangl,ntbond,ntdihd,ntfree,ntinv,ntpfbp,
     x    ntpmet,ntptbp,ntpter,ntpvdw,ntshl,ntteth,ntcons,alpha,
     x    delr,dlrpot,drewd,elrc,virlrc,epsq,fmax,opttol,rctter,
     x    rcut,rcutfb,rcuttb,rprim,rvdw,temp,tstep,volm,sigma,
     x    tboost,hyp_units)
        
      endif
      
cc     perform selected NEB calculation
c        
c      if(lneb)then
c        
c        do i=1,numneb
c
c          call neb_driver
c     x      (lfcap,lneut,lnsq,loglnk,lzeql,newlst,lneb,bsn_1(i),
c     x      bsn_2(i),idnode,mxnode,natms,imcon,nstep,nstbgr,nsteql,
c     x      keytol,keyfce,kmax1,kmax2,kmax3,nhko,nlatt,ntpvdw,ngrp,
c     x      ntcons,ntpmet,nospl,multt,nneut,ntptbp,ntpfbp,ntpter,
c     x      keyshl,ntfree,keyfld,ntbond,ntangl,ntdihd,ntinv,ntteth,
c     x      ntshl,nscons,delr,dlrpot,engcpe,engsrp,epsq,rcut,
c     x      rprim,rvdw,vircpe,virsrp,alpha,drewd,volm,
c     x      engmet,virmet,elrc,virlrc,rcuttb,engtbp,virtbp,rcutfb,
c     x      engfbp,virfbp,rctter,engter,virter,engbnd,virbnd,
c     x      engang,virang,engdih,virdih,enginv,virinv,engtet,
c     x      virtet,engshl,shlke,virshl,engfld,virfld,engcfg,fmax,
c     x      temp,tstep,opttol,sprneb,hyp_units)
c          
c          call scan_profile(nturn,estar)
c          
c          if(idnode.eq.0)then
c            
c            write(nrite,"(1x,120('-'))")
c            write(nrite,'(1x,"TRA",3i6,1p,4e14.5)')
c     x        bsn_1(i),bsn_2(i),nturn,estar/hyp_units
c            write(nrite,"(1x,120('-'))")
c            
c          endif
c          
c        enddo
c        
cc     bypass the MD cycle for this option
c        
c        recycle=.false.
c        
c      endif
        
c***********************************************************************
c     start of molecular dynamics calculations
c***********************************************************************
      
      do while(recycle)
        
c     increase step counter
        
        nstep=nstep+1
        recycle=(nstep.lt.nstrun)

c     store velocities for free energy or solvation simulation
        
        if(keyver.eq.0)then
          
          if(lsolva)then
            
            vxo_sol(:)=vxx(:)
            vyo_sol(:)=vyy(:)
            vzo_sol(:)=vzz(:)
            
          elseif(lfree)then
            
            vxo_fre(:)=vxx(:)
            vyo_fre(:)=vyy(:)
            vzo_fre(:)=vzz(:)

          endif
          
        endif
        
c     molecular switching option for excitation
        
        if(lswitch)then
          
          if(nstep.ge.nswitch)then
            
            if(mod(nstep-nswitch,niswitch).eq.0)then
              
              call switch_atm(lfrmas)
              call switch(elrc,virlrc)
              llswitch=.not.llswitch
              
            endif
            
          endif
          
        endif
        
        
c     conserved quantity (other than K + U)
        
        consv=0.d0
        
c     energy accumulators
        
        engke=0.d0
        engrot=0.d0
        
        
c     calculate volume of simulation cell
        
        if(imcon.ne.0.and.imcon.ne.6)then
          
          call dcell(cell,celprp)
          volm=celprp(10)
          if(imcon.eq.4)then
            
            volm=0.5d0*celprp(10)
            
          elseif(imcon.eq.5)then
            
            volm=0.5d0*celprp(10)
            
          elseif(imcon.eq.7)then
            
            volm=0.5d0*celprp(10)
            
          endif
          
        else
          
          volm=0.d0
          
        endif
        
c     reset sutton chen long range corrections (constant pressure only)
        
        if(ntpmet.gt.0.and.keyens.ge.4.and.keyens.le.7) then
          
          call lrcmetal
     x      (idnode,imcon,natms,ntpatm,engunit,rvdw,volm)
          
        endif
        
c     activate the impact option at designated time step
        
        if(lhit.and.nstep.eq.nhit)call impact
     x    (khit,natms,idnode,mxnode,ehit,xhit,yhit,zhit)

        
c     integrate equations of motion stage 1 of velocity verlet
        
        if(keyver.gt.0)then
            
            if(llswitch)call copy_force(idnode,mxnode)
            
            isw=1
            call vv_integrate
     x      (lcnb,lshmov,isw,idnode,mxnode,imcon,natms2,ngrp,keyens,
     x      nscons,ntcons,ntpatm,ntfree,nspmf,ntpmf,mode,tstep,engke,
     x      engrot,tolnce,vircon,vircom,virtot,temp,press,volm,sigma,
     x      taut,taup,chit,chip,consv,conint,elrc,virlrc,virpmf)
            
            if(lghost)call update_ghost(idnode,mxnode)
            
            if(lfree.or.lghost)
     x        call lrcorrect_fre(lfree,volm,elrc,virlrc)
            if(lsolva)call lrcorrect_sol(lghost,volm)

c     scale t=0 tether reference positions (constant pressure only)
          
            if(keyens.ge.4.and.keyens.le.7) then
            
              call xscale(idnode,mxnode,natms,keyens,imcon,tstep)
            
            endif
          
        endif
ccs
crs        if(qeq)then
crs          call run_qeq()
crs        endif
ccs
        if(keyshl.ne.2)then

cRS added ntpatm for recomputing lrcorrection when lambda is changing
          call molecular_dynamics
     x      (lfcap,lgofr,lneut,lnsq,loglnk,loptim,lzeql,lzero,
     x      newlst,stropt,recycle,ltad,lsolva,lfree,lghost,
     x      idnode,imcon,keyfce,keyfld,keyshl,keystr,keytol,kmax1,
     x      kmax2,kmax3,multt,mxnode,natms,ngrp,nhko,nlatt,nneut,
     x      nospl,nscons,nstbgr,nstep,nsteql,ntangl,ntbond,ntdihd,
     x      ntfree,ntinv,ntpfbp,ntpmet,ntptbp,ntpter,ntpvdw,ntshl,
     x      ntpatm,
     x      ntteth,ntcons,numrdf,nsolva,isolva,alpha,delr,dlrpot,
     x      drewd,elrc,engang,engbnd,engcpe,engdih,engfbp,engfld,
     x      enginv,engshl,engsrp,engtbp,engter,engtet,epsq,fmax,
     x      opttol,rctter,rcut,rcutfb,rcuttb,rprim,rvdw,shlke,
     x      engcfg,temp,tstep,virang,virbnd,vircpe,virdih,
     x      virfbp,virfld,virinv,virlrc,virmet,virshl,virsrp,
     x      virtbp,virter,virtet,volm,engmet,virtot)
          
        else


          call shell_relaxation
     x      (lfcap,lgofr,lneut,lnsq,loglnk,lzeql,newlst,ltad,lsolva,
     x      lfree,lghost,idnode,imcon,keyfce,keyfld,keyshl,
     x      kmax1,kmax2,kmax3,multt,mxnode,natms,nhko,nlatt,nneut,
     x      nospl,nstbgr,nstep,nsteql,ntangl,ntbond,ntdihd,ntinv,
     x      ntpfbp,ntpmet,ntptbp,ntpter,ntpvdw,ntshl,
     x      ntpatm,ntteth,numrdf,
     x      ntpmls,nsolva,isolva,alpha,delr,dlrpot,drewd,elrc,engang,
     x      engbnd,engcpe,engdih,engfbp,engfld,enginv,engshl,engsrp,
     x      engtbp,engter,engtet,epsq,fmax,rctter,rcut,rcutfb,rcuttb,
     x      rprim,rvdw,shlke,engcfg,temp,tstep,virang,virbnd,vircpe,
     x      virdih,virfbp,virfld,virinv,virlrc,virmet,virshl,virsrp,
     x      virtbp,virter,virtet,volm,engmet,virtot,rlxtol,pass0,
     x      pass1,pass2)
          
        endif
        
c     bias potential dynamics option - reset forces
        
        if(lbpd)call bpd_forces
     x    (natms,vmin,ebias,temp,engcfg,boost)
        
c     switching option for excitation simulation
        
        if(llswitch)call copy_force(idnode,mxnode)
        
c     integrate equations of motion
        
        if(keyver.eq.0)then
          
c     integrate equations of motion by leapfrog verlet
          
          call lf_integrate
     x      (lcnb,lshmov,idnode,mxnode,imcon,natms2,ngrp,keyens,
     x      nscons,ntcons,ntpatm,ntfree,nspmf,ntpmf,mode,tstep,engke,
     x      engrot,tolnce,quattol,vircon,vircom,virtot,temp,press,
     x      volm,sigma,taut,taup,chit,chip,consv,conint,elrc,
     x      virlrc,virpmf)
          
        else if(keyver.gt.0)then
          
c     integrate equations of motion by velocity verlet (stage 2)
          
          isw=2
          call vv_integrate
     x      (lcnb,lshmov,isw,idnode,mxnode,imcon,natms2,ngrp,keyens,
     x      nscons,ntcons,ntpatm,ntfree,nspmf,ntpmf,mode,tstep,engke,
     x      engrot,tolnce,vircon,vircom,virtot,temp,press,
     x      volm,sigma,taut,taup,chit,chip,consv,conint,elrc,
     x      virlrc,virpmf)
          
        endif
        
c     update the atomic positions for the ghost molecule
        
        if(lghost)call update_ghost(idnode,mxnode)
        
c     long range correction adjustment for free energy and solvation
        
        if(lsolva)call lrcorrect_sol(lghost,volm)
        
        if(lfree.or.lghost)
     x    call lrcorrect_fre(lfree,volm,elrc,virlrc)
        if(lsolva)call lrcorrect_sol(lghost,volm)
        
c     application of transition analysis procedures
        
        if(ltad.or.(lbpd.and.keybpd.eq.2))then
          
          engtke=engke+engrot

          call hyper_driver
     x      (ltad,lbpd,recycle,lfcap,lneut,lnsq,loglnk,lzeql,newlst,
     x      prechk,tadall,nblock,ntrack,idnode,imcon,keyfce,keyfld,
     x      keyshl,keytol,kmax1,kmax2,kmax3,multt,mxnode,natms,ngrp,
     x      ntcons,nhko,nlatt,nneut,nospl,nscons,nstbgr,nstep,
     x      nsteql,ntangl,ntbond,ntdihd,ntfree,ntinv,ntpfbp,
     x      ntpmet,ntptbp,ntpter,ntpvdw,ntshl,ntteth,blkout,
     x      alpha,delr,dlrpot,drewd,elrc,virlrc,epsq,fmax,
     x      opttol,rctter,rcut,rcutfb,rcuttb,rprim,rvdw,temp,
     x      tstep,volm,engcfg,catchrad,sprneb,deltad,tlow,engtke,
     x      tolnce,tboost,hyp_units)
          
        endif
        
c     calculate shell kinetic energy
        
        if(keyshl.eq.1)then

          call corshl(idnode,mxnode,ntshl,shlke)
          
        endif
        
c     scale t=0 tether reference positions (constant pressure only)
        
        if(keyver.eq.0.and.keyens.ge.4.and.keyens.le.7) then
          
          call xscale(idnode,mxnode,natms,keyens,imcon,tstep)
          
        endif
        
c     apply temperature scaling
        
        if((ltscal.and.nstep.le.nsteql).and.
     x    mod(nstep-nsteql,nstbts).eq.0)then
          
          chit=0.d0
          chip=0.d0
          do i=1,9
            eta(i)=0.d0
          enddo
          
          if(keyshl.eq.1) then
            
            do k=1,4
              
              call vscaleg(idnode,mxnode,imcon,natms2,ngrp,sigma)
              call shlqnch(idnode,mxnode,ntshl,temp)
              
            enddo
            
          else
            
            call vscaleg(idnode,mxnode,imcon,natms2,ngrp,sigma)
            
          endif
          
        endif
        
c     reset atom velocities at intervals if required
        
        if(newgau)then
          
          if(mod(nstep,numgau).eq.0)call regauss
     x      (idnode,imcon,mxnode,natms2,ngrp,nscons,ntcons,
     x      ntshl,keyshl,sigma,temp,tolnce)
          
        endif

c     calculate physical quantities


        if(nstep.gt.0)call static
     x    (lbpd,lzeql,idnode,intsta,imcon,keyens,natms,nstack,
     x    nstep,nsteql,ntpatm,numacc,mxnode,nblock,consv,degfre,
     x    degrot,engang,engbnd,engcpe,engdih,enginv,engke,engrot,
     x    engsrp,engunit,engcfg,stpeng,stpeth,stpprs,stptmp,stpvir,
     x    stpvol,tstep,virbnd,engfbp,vircom,vircon,vircpe,virsrp,
     x    engfld,virfld,engtbp,virtbp,virpmf,virshl,engshl,engtet,
     x    virtet,degshl,shlke,virang,width,engmet,virmet,engter,
     x    virter,boost,tboost,ebias,heinc)
        
c     z density calculation
        
        if(lzden.and.((.not.lzeql).or.(nstep.gt.nsteql))) then
          
          call zden0(idnode,natms,mxnode,nzden,zlen)
          
        endif
        
c     terminate program if boundary conditions violated
        
        if(imcon.gt.0.and.rcut.gt.width)then
          
          levcfg=2
          call revive
     x      (lgofr,lzden,idnode,imcon,mxnode,natms,levcfg,nstep,nzden,
     x      numacc,numrdf,chip,chit,conint,tstep,engcfg,virtot,vircom,
     x      tboost)
          call error(idnode,95)
          
        endif
        
c     line-printer output every nstbpo steps
        
CVAM
CVAM          call VTBEGIN(68, ierr)
CVAM
        
        if(nstep.eq.1.or.(nstep.gt.1.and.mod(nstep,nstbpo).eq.0))then
           
          call timchk(0,timelp)
          if(idnode.eq.0)then
            
            if(mod(lines,npage).eq.0)
     x        write(nrite,"(1x,120('-'),
     x        /,/,1x,'    step',8x,'eng_tot',7x,'temp_tot',8x,
     x        'eng_cfg',8x,'eng_vdw',8x,'eng_cou',8x,'eng_bnd',
     x        8x,'eng_ang',8x,'eng_dih',8x,'eng_tet',/,1x,
     x        'time(ps)',8x,' eng_pv',7x,'temp_rot',8x,'vir_cfg',
     x        8x,'vir_vdw',8x,'vir_cou',8x,'vir_bnd',8x,'vir_ang',
     x        8x,'vir_con',8x,'vir_tet',/,1x,'cpu  (s)',9x,
     x        'volume',7x,'temp_shl',8x,'eng_shl',8x,'vir_shl',
     x        10x,'alpha',11x,'beta',10x,'gamma',8x,'vir_pmf',
     x        10x,'press',/,/,
     x        1x,120('-'))")          
            write(nrite,"(1x,i8,1p,9e15.6,/,1x,0p,f8.1,1p,9e15.6,
     x        /,1x,0p,f8.1,1p,9e15.6)")
     x        nstep,(stpval(i),i=1,9),
     x        dble(nstep)*tstep,(stpval(i),i=10,18),
     x        timelp,(stpval(i),i=19,27)
            write(nrite,"(/,1x,' rolling',1p,9e15.6,/,1x,'averages',
     x        1p,9e15.6,/,9x,9e15.6)") (ravval(i),i=1,27)
            write(nrite,"(1x,120('-'))")

          endif
          
          lines=lines+1
          
        endif
ccs
        if(idnode.eq.0) write(999,'(i8,f14.4)') nstep,stpval(1)
ccs
CVAM
CVAM          call VTEND(68, ierr)
CVAM
c     report end of equilibration period
        
        if(nstep.ge.nsteql)then
          
          if((ltscal.and.idnode.eq.0).and.(nstep.eq.nsteql))
     x      write(nrite,"(/,/,1x,'switching off temperature ',
     x      'scaling at step ',i6,/,/,/,1x,120('-'))") nstep
          ltscal=.false.
          
        endif
        
c     write trajectory data
        
        if(ltraj.and.nstep.ge.nstraj) then
          if(idnode.eq.0.and.mod(nstep-nstraj,istraj).eq.0)then
            
            call traject
     x        (ltraj,idnode,imcon,istraj,keytrj,natms,
     x        nstraj,nstep,tstep)
            
          endif
          
        endif
        
c     write solvation energy file
        
        if(lsolva.and.nstep.ge.nsolva)then
          
          if(mod(nstep-nsolva,isolva).eq.0)then
            
            call solva_temp(idnode,mxnode,natms2,keyver)
            call solvation_write(lexcite,lswitch,idnode,natms,
     x        nstep,nsolva,isolva,tstep,engunit,elrc)
            
          endif
          
        endif
        
c     write free energy file
        
        if(lfree.and.nstep.ge.nfrn)then
          
          if(mod(nstep-nfrn,ifrn).eq.0)then
            
            call free_kinetic(lfrmas,idnode,mxnode,keyver)
            call free_energy_write(idnode,nstep,engunit)
          
          endif
          
        endif
        
c     save restart data in event of system crash
        
        if(mod(nstep,ndump).eq.0.and.nstep.ne.nstrun)then
          
          levcfg=2
          call revive
     x      (lgofr,lzden,idnode,imcon,mxnode,natms,levcfg,nstep,nzden,
     x      numacc,numrdf,chip,chit,conint,tstep,engcfg,virtot,vircom,
     x      tboost)
          
          if(ltad.or.lbpd)call hyper_close(idnode,mxnode,natms,nsteql)
          
        endif
        
c     cycle time check
        
        call timchk(0,timelp)
        recycle=(recycle.and.timjob-timelp.gt.timcls)
        
      enddo
      
c***********************************************************************
c     end of molecular dynamics calculations
c***********************************************************************
      
c     last time check
        
      call timchk(0,timelp)
        
      if(idnode.eq.0)write(nrite,
     x  "(/,/,1x,'run terminating. elapsed cpu time = ',f13.3,
     x  ', job time = ',f13.3,', close time = ',f13.3,/)")
     x  timelp,timjob,timcls
      
c     shell relaxation convergence statistics
      
      if(keyshl.eq.2)then
        
        if(idnode.eq.0)write(nrite,
     x    "(/,/,1x,'shell relaxation statistics : average cycles = ',
     x    f8.3,' maximum cycles = ',f8.3)")pass1,pass2
      
      endif

c     produce summary of simulation
      
      levcfg=2
      call result
     x  (lbpd,lgofr,lpgr,lzden,idnode,imcon,keyens,mxnode,natms,
     x  levcfg,nzden,nstep,ntpatm,numacc,numrdf,chip,chit,conint,
     x  rcut,tstep,engcfg,volm,virtot,vircom,zlen,tboost)
      
c     write hyperdynamics restart file
      
      if(ltad.or.lbpd)call hyper_close(idnode,mxnode,natms,nsteql)
ccs
      if(idnode.eq.0) close(unit=999)
ccs
      return
      end subroutine dynamics

c ########################################################
      subroutine finish_dlpoly()
      
c     close output channels

      
      if(idnode.eq.0) then
        write(nrite,"(/,/,1x,'ending dlpoly  ... ')")
        close (nrite)
        close (nstats)
        close (nhist)
        close (nevnt)
ccs
        if(time) close(unit=10)
ccs
        
      endif
      
c     terminate job
CVAM
CVAM        call VTEND(99, ierr)
CVAM
      call exitcomms(0)


      return
      
      end subroutine finish_dlpoly

c ################ EXTRA ROUTINES #####################################

      subroutine calc_energy_force()

c     time check

c      call timchk(1,tzero)

      call molecular_dynamics
     x    (lfcap,lgofr,lneut,lnsq,loglnk,loptim,lzeql,lzero,
     x    newlst,stropt,recycle,ltad,lsolva,lfree,lghost,
     x    idnode,imcon,keyfce,keyfld,keyshl,keystr,keytol,kmax1,
     x    kmax2,kmax3,multt,mxnode,natms,ngrp,nhko,nlatt,nneut,
     x    nospl,nscons,nstbgr,nstep,nsteql,ntangl,ntbond,ntdihd,
     x    ntfree,ntinv,ntpfbp,ntpmet,ntptbp,ntpter,ntpvdw,ntshl,
     x    ntpatm,
     x    ntteth,ntcons,numrdf,nsolva,isolva,alpha,delr,dlrpot,
     x    drewd,elrc,engang,engbnd,engcpe,engdih,engfbp,engfld,
     x    enginv,engshl,engsrp,engtbp,engter,engtet,epsq,fmax,
     x    opttol,rctter,rcut,rcutfb,rcuttb,rprim,rvdw,shlke,
     x    engcfg,temp,tstep,virang,virbnd,vircpe,virdih,
     x    virfbp,virfld,virinv,virlrc,virmet,virshl,virsrp,
     x    virtbp,virter,virtet,volm,engmet,virtot)
        
        
      return
      end subroutine calc_energy_force

      subroutine write_config(fname,level)
      character*(*) :: fname
      integer     :: level
      if(idnode.eq.0)call config_write(fname,level,imcon,natms,engcfg)
      end subroutine write_config

      subroutine update_celprp
      call dcell(cell,celprp)
      return
      end subroutine update_celprp

      subroutine get_atmnam(name, i)
      character*8  :: name
      integer      :: i
      name = atmnam(i)
      return
      end subroutine get_atmnam
      
      subroutine init_dlpoly_random()
        implicit none
        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed

        call random_seed(size = n)
        allocate(seed(n))

        call system_clock(count=clock)

        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        call random_seed(put = seed)

        deallocate(seed)

c        write (*,*) "Intialized Random with System Time"
        return
      end subroutine init_dlpoly_random

cRS helper function ot trigger recalc of vdw lookuptables if params have changed
      subroutine update_vdw_tables()
      implicit none
      
c we silently assume that there are no numeric vdw potentials (ltable = False)
      call forgen(.false., idnode, ntpvdw, dlrpot, rvdw)
      
      if (idnode.eq.0) then
        write (*,*) "VDW lookuptables have been recomputed!"
      endif
      
      end subroutine update_vdw_tables

cRS ####################################################################
cRS
c      in the following we split the loop into individual routines
c      one to start up the dynamics, one to be called before the 
c      energy call, one after and one to finish up (write out averages)
c      
c      if averages are doen on the fly we actually do not need this 
c      part. it looks like with nstep=0 (and resetting numacc=0)
c      we can flush the averages.
c      note also that now nstrun and recycle have no meaning any more!
c      
c      REMARK: recycle is passed to molecular_dynamics as "cycle" because
c          in case of optimiser or things like that we need to say if it is 
c          optimized (this should go evetnually as soon as we get the legacy
c          optimizers out)
c       
c *********
      subroutine dyns_start()

      implicit none
      
cRS is only doing something if nstep=0 and we use verlet!!!
ccs
cf2py intent(callback) run_qeq
c      external run_qeq
ccs
      keystr=0
      stropt=.false.
      sigma=temp*boltz*degfre*0.5d0

cRS  ---> only in the first step
c     calculate initial conditions for velocity verlet

      if(nstep.eq.0)then
        
c     kinetic stress tensor at start
        
          call dcell(cell,celprp)
          width=min(celprp(7),celprp(8),celprp(9))
          call kinstress(natms,idnode,mxnode,stress)
          engke=0.5d0*(stress(1)+stress(5)+stress(9))
          do i=1,9
             stress(i)=stress(i)/dble(mxnode)
          enddo
        
cRS NOTE: this should not be necessary. the point is to compute the force
c         in the first step before the first propagation step.
c         however, pydlpoly we should have done this after startup anyway
c         here, no extra force routines are called!!! (like qeq etc.)
c         -> this should go away in the future

c          call molecular_dynamics
c     x    (lfcap,lgofr,lneut,lnsq,loglnk,loptim,lzeql,lzero,
c     x    newlst,stropt,recycle,ltad,lsolva,lfree,lghost,
c     x    idnode,imcon,keyfce,keyfld,keyshl,keystr,keytol,kmax1,
c     x    kmax2,kmax3,multt,mxnode,natms,ngrp,nhko,nlatt,nneut,
c     x    nospl,nscons,nstbgr,nstep,nsteql,ntangl,ntbond,ntdihd,
c     x    ntfree,ntinv,ntpfbp,ntpmet,ntptbp,ntpter,ntpvdw,ntshl,
c     x    ntteth,ntcons,numrdf,nsolva,isolva,alpha,delr,dlrpot,
c     x    drewd,elrc,engang,engbnd,engcpe,engdih,engfbp,engfld,
c     x    enginv,engshl,engsrp,engtbp,engter,engtet,epsq,fmax,
c     x    opttol,rctter,rcut,rcutfb,rcuttb,rprim,rvdw,shlke,
c     x    engcfg,temp,tstep,virang,virbnd,vircpe,virdih,
c     x    virfbp,virfld,virinv,virlrc,virmet,virshl,virsrp,
c     x    virtbp,virter,virtet,volm,engmet,virtot)
        
      endif
cRS ---> end of things in the first step (if verlet)
      end subroutine dyns_start
        
c *********        
        
      subroutine dyns_step_fh() 
      
      implicit none 
        
cRS   leave this condition here for the time being ...
cRS      but since the while loop is gone it should have no effect        
c     increase step counter
        
        nstep=nstep+1
        recycle=(nstep.lt.nstrun)

c     conserved quantity (other than K + U)        
        consv=0.d0
        
c     energy accumulators        
        engke=0.d0
        engrot=0.d0
        
        
c     calculate volume of simulation cell        
        if(imcon.ne.0.and.imcon.ne.6)then          
          call dcell(cell,celprp)
          volm=celprp(10)
          if(imcon.eq.4)then            
            volm=0.5d0*celprp(10)
          elseif(imcon.eq.5)then            
            volm=0.5d0*celprp(10)            
          elseif(imcon.eq.7)then
            volm=0.5d0*celprp(10)
          endif
        else          
          volm=0.d0          
        endif
                
c     integrate equations of motion stage 1 of velocity verlet
        
        isw=1
        call vv_integrate
     x      (lcnb,lshmov,isw,idnode,mxnode,imcon,natms2,ngrp,keyens,
     x      nscons,ntcons,ntpatm,ntfree,nspmf,ntpmf,mode,tstep,engke,
     x      engrot,tolnce,vircon,vircom,virtot,temp,press,volm,sigma,
     x      taut,taup,chit,chip,consv,conint,elrc,virlrc,virpmf)
            

c     scale t=0 tether reference positions (constant pressure only)
        if(keyens.ge.4.and.keyens.le.7) then            
            call xscale(idnode,mxnode,natms,keyens,imcon,tstep)
        endif
          
        return

      end subroutine dyns_step_fh

      subroutine dyns_step_sh()          
                
      implicit none          
                
c     integrate equations of motion by velocity verlet (stage 2)
        isw=2
        call vv_integrate
     x      (lcnb,lshmov,isw,idnode,mxnode,imcon,natms2,ngrp,keyens,
     x      nscons,ntcons,ntpatm,ntfree,nspmf,ntpmf,mode,tstep,engke,
     x      engrot,tolnce,vircon,vircom,virtot,temp,press,
     x      volm,sigma,taut,taup,chit,chip,consv,conint,elrc,
     x      virlrc,virpmf)
      
c     apply temperature scaling
        if((ltscal.and.nstep.le.nsteql).and.
     x    mod(nstep-nsteql,nstbts).eq.0)then          
          chit=0.d0
          chip=0.d0
          do i=1,9
            eta(i)=0.d0
          enddo
          call vscaleg(idnode,mxnode,imcon,natms2,ngrp,sigma)
        endif
        
c     reset atom velocities at intervals if required
        if(newgau)then          
          if(mod(nstep,numgau).eq.0)call regauss
     x      (idnode,imcon,mxnode,natms2,ngrp,nscons,ntcons,
     x      ntshl,keyshl,sigma,temp,tolnce)          
        endif

! c     calculate physical quantities
!         if(nstep.gt.0)call static
!      x    (lbpd,lzeql,idnode,intsta,imcon,keyens,natms,nstack,
!      x    nstep,nsteql,ntpatm,numacc,mxnode,nblock,consv,degfre,
!      x    degrot,engang,engbnd,engcpe,engdih,enginv,engke,engrot,
!      x    engsrp,engunit,engcfg,stpeng,stpeth,stpprs,stptmp,stpvir,
!      x    stpvol,tstep,virbnd,engfbp,vircom,vircon,vircpe,virsrp,
!      x    engfld,virfld,engtbp,virtbp,virpmf,virshl,engshl,engtet,
!      x    virtet,degshl,shlke,virang,width,engmet,virmet,engter,
!      x    virter,boost,tboost,ebias,heinc)

cRS add here what is needed from static to get ekin and pressure etc
c     calculate cell volume and minimum cell half-width
        if(imcon.eq.0)then
          width=0.d0
          stpvol=0.d0
          do i=1,10
            celprp(i)=0.d0
          enddo
        else
          call dcell(cell,celprp)
          stpvol=celprp(10)
          width=min(celprp(7),celprp(8),celprp(9))/2.d0
          if(imcon.eq.4)then
            stpvol=0.5d0*celprp(10)
            width=sqrt(3.d0)*cell(1)/4.d0
          elseif(imcon.eq.5)then
            stpvol=0.5d0*celprp(10)
            width=cell(1)/2.d0
          elseif(imcon.eq.6)then
            width=min(celprp(7),celprp(8))/2.d0
          elseif(imcon.eq.7)then
            stpvol=0.5d0*celprp(10)
          endif
        endif
c       energetic properties of system
        stpvir=virsrp+vircpe+virbnd+vircon+vircom+virtbp+virang
     x    +virshl+virtet+virter+virmet+virfld
        stpeng=engcfg+engke+engrot
        stprot=2.d0*engrot/(boltz*max(1.d0,degrot))
        stpshl=2.d0*shlke/(boltz*max(1.d0,degshl))
        stptmp=2.d0*(engke+engrot)/(boltz*degfre)
        stpprs=0.d0
        if(imcon.gt.0)stpprs=(2.d0*engke-stpvir)/(3.d0*stpvol)
        stpeth=stpeng+stpprs*stpvol
        stpcns=stpeng+consv
c       convert pressure to units of katm
        stpprs=stpprs*prsunt

        
c     z density calculation
c        if(lzden.and.((.not.lzeql).or.(nstep.gt.nsteql))) then
c          call zden0(idnode,natms,mxnode,nzden,zlen)          
c        endif
        
c     terminate program if boundary conditions violated
        if(imcon.gt.0.and.rcut.gt.width)then
c          levcfg=2
c          call revive
c     x      (lgofr,lzden,idnode,imcon,mxnode,natms,levcfg,nstep,nzden,
c     x      numacc,numrdf,chip,chit,conint,tstep,engcfg,virtot,vircom,
c     x      tboost)
          call error(idnode,95)
        endif
        
! c write out intermediate results
!         if(nstep.eq.1.or.(nstep.gt.1.and.mod(nstep,nstbpo).eq.0))then
!           call timchk(0,timelp)
!           if(idnode.eq.0)then
!             if(mod(lines,npage).eq.0)
!      x        write(nrite,"(1x,120('-'),
!      x        /,/,1x,'    step',8x,'eng_tot',7x,'temp_tot',8x,
!      x        'eng_cfg',8x,'eng_vdw',8x,'eng_cou',8x,'eng_bnd',
!      x        8x,'eng_ang',8x,'eng_dih',8x,'eng_tet',/,1x,
!      x        'time(ps)',8x,' eng_pv',7x,'temp_rot',8x,'vir_cfg',
!      x        8x,'vir_vdw',8x,'vir_cou',8x,'vir_bnd',8x,'vir_ang',
!      x        8x,'vir_con',8x,'vir_tet',/,1x,'cpu  (s)',9x,
!      x        'volume',7x,'temp_shl',8x,'eng_shl',8x,'vir_shl',
!      x        10x,'alpha',11x,'beta',10x,'gamma',8x,'vir_pmf',
!      x        10x,'press',/,/,
!      x        1x,120('-'))")          
!             write(nrite,"(1x,i8,1p,9e15.6,/,1x,0p,f8.1,1p,9e15.6,
!      x        /,1x,0p,f8.1,1p,9e15.6)")
!      x        nstep,(stpval(i),i=1,9),
!      x        dble(nstep)*tstep,(stpval(i),i=10,18),
!      x        timelp,(stpval(i),i=19,27)
!             write(nrite,"(/,1x,' rolling',1p,9e15.6,/,1x,'averages',
!      x        1p,9e15.6,/,9x,9e15.6)") (ravval(i),i=1,27)
!             write(nrite,"(1x,120('-'))")
!           endif
!           lines=lines+1
!         endif

c     report end of equilibration period
c        if(nstep.ge.nsteql)then          
c          if((ltscal.and.idnode.eq.0).and.(nstep.eq.nsteql))
c     x      write(nrite,"(/,/,1x,'switching off temperature ',
c     x      'scaling at step ',i6,/,/,/,1x,120('-'))") nstep
c          ltscal=.false.
c        endif
        
c     write trajectory data
c        if(ltraj.and.nstep.ge.nstraj) then
c          if(idnode.eq.0.and.mod(nstep-nstraj,istraj).eq.0)then
c            call traject
c     x        (ltraj,idnode,imcon,istraj,keytrj,natms,
c     x        nstraj,nstep,tstep)
c          endif
c        endif
                
c     save restart data in event of system crash
c        if(mod(nstep,ndump).eq.0.and.nstep.ne.nstrun)then
c          levcfg=2
c          call revive
c     x      (lgofr,lzden,idnode,imcon,mxnode,natms,levcfg,nstep,nzden,
c     x      numacc,numrdf,chip,chit,conint,tstep,engcfg,virtot,vircom,
c     x      tboost)
c        endif
        
      end subroutine dyns_step_sh
      

c *********        

c      subroutine dyns_finish()

c      implicit none

c     produce summary of simulation
c      levcfg=2
c      call result
c     x  (lbpd,lgofr,lpgr,lzden,idnode,imcon,keyens,mxnode,natms,
c     x  levcfg,nzden,nstep,ntpatm,numacc,numrdf,chip,chit,conint,
c     x  rcut,tstep,engcfg,volm,virtot,vircom,zlen,tboost)
      
c      return
      
c      end subroutine dyns_finish



cRS ####################################################################


      end module pydlpoly



