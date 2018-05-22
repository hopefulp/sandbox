c***********************************************************************
c
c     Module pypikaia (R. Schmid, RUB, 2014)
c
c     Changes:
c     - use random_number instead of the ran0 functions
c     - convert to a module to be used from python
c     - make array allocatable (instead of fixed sized arrays)
c       -> call an init function ...
c     - main pikaia routine does only one generation (loop done in python)
c     - add options for restart files and setup with an exisiting paramter set
c
c
c***********************************************************************


      module pypikaia
      
      implicit none

cRS Meaning of control variables
c           np       - number of individuals in a population (default
c                      is 100)
c           ngen     - number of generations over which solution is
c                      to evolve (default is 500)
c           nd       - number of significant digits (i.e., number of
c                      genes) retained in chromosomal encoding (default
c                      is 6)  (Note: This number is limited by the
c                      machine floating point precision.  Most 32-bit
c                      floating point representations have only 6 full
c                      digits of precision.  To achieve greater preci-
c                      sion this routine could be converted to double
c                      precision, but note that this would also require
c                      a double precision random number generator, which
c                      likely would not have more than 9 digits of
c                      precision if it used 4-byte integers internally.)
c           pcross   - crossover probability; must be  <= 1.0 (default
c                      is 0.85). If crossover takes place, either one
c                      or two splicing points are used, with equal
c                      probabilities
c           imut     - mutation mode; 1/2/3/4/5 (default is 2)
c                      1=one-point mutation, fixed rate
c                      2=one-point, adjustable rate based on fitness 
c                      3=one-point, adjustable rate based on distance
c                      4=one-point+creep, fixed rate
c                      5=one-point+creep, adjustable rate based on fitness
c                      6=one-point+creep, adjustable rate based on distance
c           pmut     - initial mutation rate; should be small (default
c                      is 0.005) (Note: the mutation rate is the proba-
c                      bility that any one gene locus will mutate in
c                      any one generation.)
c           pmutmn   - minimum mutation rate; must be >= 0.0 (default
c                      is 0.0005)
c           pmutmx   - maximum mutation rate; must be <= 1.0 (default
c                      is 0.25)
c           fdif     - relative fitness differential; range from 0
c                      (none) to 1 (maximum).  (default is 1.)
c           irep     - reproduction plan; 1/2/3=Full generational
c                      replacement/Steady-state-replace-random/Steady-
c                      state-replace-worst (default is 3)
c           ielite   - elitism flag; 0/1=off/on (default is 0)
c                      (Applies only to reproduction plans 1 and 2)
c           ivrb     - printed output 0/1/2=None/Minimal/Verbose
c                      (default is 0)

      real(8), save, allocatable :: oldph(:,:)
      real(8), save, allocatable :: newph(:,:)
      real(8), save, allocatable :: fitns(:)
      integer, save, allocatable :: ifit(:), jfit(:)
      
      integer, save              :: imut, irep, ielite, ivrb
      real(8), save              :: pcross,pmut,pmutmn,pmutmx,fdif
      
      logical, save              :: ftns_calculated = .false.
      logical, save              :: ftns_sorted     = .false.
            
      contains

c***********************************************************************
      subroutine setup(n, np, nd, seed, status)
c=======================================================================
      
      implicit none
      
      integer n, np, nd
      integer fail
      integer status
      integer seed
      
      integer, dimension(:), allocatable :: seeds
      integer nseed, iseed
      
      allocate(oldph(n,np), stat=fail)
      if (fail.gt.0) then
        write (*,*) "error during memory allocation  ", fail
        status = 100
      end if
      allocate(newph(n,np), stat=fail)
      if (fail.gt.0) then
        write (*,*) "error during memory allocation  ", fail
        status = 101
      end if
      allocate(fitns(np), stat=fail)
      if (fail.gt.0) then
        write (*,*) "error during memory allocation  ", fail
        status = 102
        end if
      allocate(ifit(np), stat=fail)
      if (fail.gt.0) then
        write (*,*) "error during memory allocation  ", fail
        status = 103
      end if
      allocate(jfit(np), stat=fail)
      if (fail.gt.0) then
        write (*,*) "error during memory allocation  ", fail
        status = 104
      end if
 
         write(*,2) np,n,nd,pcross,pmut,pmutmn,pmutmx,fdif
    2    format(/1x,60('*'),/,
     +      ' *',13x,'PIKAIA Genetic Algorithm Report ',13x,'*',/,
     +           1x,60('*'),//,
     +      '       Individuals per generation: ',i4,/,
     +      '    Number of Chromosome segments: ',i4,/,
     +      '    Length of Chromosome segments: ',i4,/,
     +      '            Crossover probability: ',f9.4,/,
     +      '            Initial mutation rate: ',f9.4,/,
     +      '            Minimum mutation rate: ',f9.4,/,
     +      '            Maximum mutation rate: ',f9.4,/,
     +      '    Relative fitness differential: ',f9.4)
         if (imut.eq.1) write(*,3) 'Uniform, Constant Rate'
         if (imut.eq.2) write(*,3) 'Uniform, Variable Rate (F)'
         if (imut.eq.3) write(*,3) 'Uniform, Variable Rate (D)'
         if (imut.eq.4) write(*,3) 'Uniform+Creep, Constant Rate'
         if (imut.eq.5) write(*,3) 'Uniform+Creep, Variable Rate (F)'
         if (imut.eq.6) write(*,3) 'Uniform+Creep, Variable Rate (D)'
    3    format(
     +      '                    Mutation Mode: ',A)
         if (irep.eq.1) write(*,4) 'Full generational replacement'
         if (irep.eq.2) write(*,4) 'Steady-state-replace-random'
         if (irep.eq.3) write(*,4) 'Steady-state-replace-worst'
    4    format(
     +      '                Reproduction Plan: ',A)
 
c     Check some control values
      if (imut.ne.1 .and. imut.ne.2 .and. imut.ne.3 .and. imut.ne.4
     +    .and. imut.ne.5 .and. imut.ne.6) then
         write(*,10)
         status = 5
      endif
   10 format(' ERROR: illegal value for imut (ctrl(5))')
 
      if (fdif.gt.1.) then
         write(*,11)
         status = 9
      endif
   11 format(' ERROR: illegal value for fdif (ctrl(9))')
 
      if (irep.ne.1 .and. irep.ne.2 .and. irep.ne.3) then
         write(*,12)
         status = 10
      endif
   12 format(' ERROR: illegal value for irep (ctrl(10))')
 
      if (pcross.gt.1.0 .or. pcross.lt.0.) then
         write(*,13)
         status = 4
      endif
   13 format(' ERROR: illegal value for pcross (ctrl(4))')
 
      if (ielite.ne.0 .and. ielite.ne.1) then
         write(*,14)
         status = 11
      endif
   14 format(' ERROR: illegal value for ielite (ctrl(11))')
 
      if (irep.eq.1 .and. imut.eq.1 .and. pmut.gt.0.5 .and.
     +    ielite.eq.0) then
         write(*,15)
      endif
   15 format(' WARNING: dangerously high value for pmut (ctrl(6));',
     +      /' (Should enforce elitism with ctrl(11)=1.)')
 
      if (irep.eq.1 .and. imut.eq.2 .and. pmutmx.gt.0.5 .and.
     +    ielite.eq.0) then
         write(*,16)
      endif
   16 format(' WARNING: dangerously high value for pmutmx (ctrl(8));',
     +      /' (Should enforce elitism with ctrl(11)=1.)')
 
      if (fdif.lt.0.33 .and. irep.ne.3) then
         write(*,17)
      endif
   17 format(' WARNING: dangerously low value of fdif (ctrl(9))')
 
      if (mod(np,2).gt.0) then
         np=np-1
         write(*,18) np
      endif
   18 format(' WARNING: decreasing population size (ctrl(1)) to np=',i4)

cRS   use seed to initalize the random_number generator if seed is larger then 0
      if (seed.gt.0) then
        write (*,19) seed
   19 format(' Initalizing random seed with ',i8)
        call random_seed(size = nseed)
        allocate(seeds(nseed))
        seeds = seed + 37 * (/ (iseed - 1, iseed = 1, nseed) /)
        call random_seed(put = seeds)
        deallocate(seeds)
          
      end if
      
      end subroutine setup

      subroutine initialize(ff, n, np)

      implicit none
 
c     Input:
      integer   n, np
      real(8)   ff   
      external  ff
      
      integer   ip
      
cRS   we assume that the oldph is somehow initalized and we compute the 
cRS   initial fitness only if needed
      if (.not.ftns_calculated) then
        do ip=1,np
          fitns(ip) = ff(n,oldph(1,ip))
        end do
        ftns_calculated = .true.
      end if
c     Rank initial population by fitness order if needed
      if (.not.ftns_sorted) then
        call rnkpop(np,fitns,ifit,jfit)
        ftns_sorted = .true.
      end if
      
      end subroutine initialize
      
c***********************************************************************
      subroutine pikaia(ff,n,np,nd,ngen)
c=======================================================================
c     Optimization (maximization) of user-supplied "fitness" function
c     ff  over n-dimensional parameter space  x  using a basic genetic
c     algorithm method.
c
c     Paul Charbonneau & Barry Knapp
c     High Altitude Observatory
c     National Center for Atmospheric Research
c     Boulder CO 80307-3000
c     <paulchar@hao.ucar.edu>
c     <knapp@hao.ucar.edu>
c
c     Version 1.2   [ 2002 April 3 ]
c
c     Genetic algorithms are heuristic search techniques that
c     incorporate in a computational setting, the biological notion
c     of evolution by means of natural selection.  This subroutine
c     implements the three basic operations of selection, crossover,
c     and mutation, operating on "genotypes" encoded as strings.
c
c     Version 1.2 differs from version 1.0 (December 1995) in that
c     it includes (1) two-point crossover, (2) creep mutation, and
c     (3) dynamical adjustment of the mutation rate based on metric
c     distance in parameter space.
c
c     References:
c
c        Charbonneau, Paul. "An introduction to gemetic algorithms for
c           numerical optimization", NCAR Technical Note TN-450+IA
c           (April 2002)
c
c        Charbonneau, Paul. "Release Notes for PIKAIA 1.2",
c           NCAR Technical Note TN-451+STR (April 2002)
c
c        Charbonneau, Paul, and Knapp, Barry. "A User's Guide
c           to PIKAIA 1.0" NCAR Technical Note TN-418+IA
c           (December 1995)
c
c        Goldberg, David E.  Genetic Algorithms in Search, Optimization,
c           & Machine Learning.  Addison-Wesley, 1989.
c
c        Davis, Lawrence, ed.  Handbook of Genetic Algorithms.
c           Van Nostrand Reinhold, 1991.
c
c=======================================================================
c     USES: ff, urand, setctl, report, rnkpop, select, encode, decode,
c           cross, mutate, genrep, stdrep, newpop, adjmut
      implicit none
 
c     Input:
      integer   n
      real(8)      ff   
      external  ff
c
c      o Integer  n  is the parameter space dimension, i.e., the number
c        of adjustable parameters. 
c
c      o Function  ff  is a user-supplied scalar function of n vari-
c        ables, which must have the calling sequence f = ff(n,x), where
c        x is a real parameter array of length n.  This function must
c        be written so as to bound all parameters to the interval [0,1];
c        that is, the user must determine a priori bounds for the para-
c        meter space, and ff must use these bounds to perform the appro-
c        priate scalings to recover true parameter values in the
c        a priori ranges.
c
c        By convention, ff should return higher values for more optimal
c        parameter values (i.e., individuals which are more "fit").
c        For example, in fitting a function through data points, ff
c        could return the inverse of chi**2.
c
c
c     Local variables
      integer        np, nd, ngen, k, ip, ig,
     +               ip1, ip2, new, newtot
c
      real(8)        ph(n,2)
c
      integer        gn1(n*nd), gn2(n*nd)
c
c
c     if necessary inital fitness needs to be computed and ranked
      call initialize(ff, n, np)
 
c     Main Population Loop 
      newtot=0
      do ip=1,np/2
 
c           1. pick two parents
            call select(np,jfit,fdif,ip1)
   21       call select(np,jfit,fdif,ip2)
            if (ip1.eq.ip2) goto 21
 
c           2. encode parent phenotypes
            call encode(n,nd,oldph(1,ip1),gn1)
            call encode(n,nd,oldph(1,ip2),gn2)
 
c           3. breed
            call cross(n,nd,pcross,gn1,gn2)
            call mutate(n,nd,pmut,gn1,imut)
            call mutate(n,nd,pmut,gn2,imut)
 
c           4. decode offspring genotypes
            call decode(n,nd,gn1,ph(1,1))
            call decode(n,nd,gn2,ph(1,2))
 
c           5. insert into population
            if (irep.eq.1) then
               call genrep(n,np,ip,ph)
            else
               call stdrep(ff,n,np,irep,ielite,
     +                     ph,ifit,jfit,new)
               newtot = newtot+new
            endif
 
c        End of Main Population Loop
      end do
 
c     if running full generational replacement: swap populations
      if (irep.eq.1)
     +      call newpop(ff,ielite,n,np, ifit,jfit,newtot)
 
c        adjust mutation rate?
         if (imut.eq.2 .or. imut.eq.3 .or. imut.eq.5 .or. imut.eq.6)
     +      call adjmut(n,np,ifit,pmutmn,pmutmx,
     +                  pmut,imut)
c
         if (ivrb.gt.0) call report
     +      (ivrb,n,np,nd,ifit,pmut,ig,newtot)
c
      end subroutine

c********************************************************************
      subroutine report
     +   (ivrb,n,np,nd,ifit,pmut,ig,nnew)
c
c     Write generation report to standard output
c
      implicit none
 
c     Input:
      integer np,ifit(np),ivrb,n,nd,ig,nnew
      real(8) pmut
c
c     Output: none
c
c     Local
      real(8) bestft,pmutpv
      save bestft,pmutpv
      integer ndpwr,k
      logical rpt
      data bestft,pmutpv /0,0/
c
      rpt=.false.
 
      if (pmut.ne.pmutpv) then
         pmutpv=pmut
         rpt=.true.
      endif
 
      if (fitns(ifit(np)).ne.bestft) then
         bestft=fitns(ifit(np))
         rpt=.true.
      endif
 
      if (rpt .or. ivrb.ge.2) then
 
c        Power of 10 to make integer genotypes for display
         ndpwr = nint(10.**nd)
 
         write(*,'(/i6,i6,f10.6,4f10.6)') ig,nnew,pmut,
     +      fitns(ifit(np)), fitns(ifit(np-1)), fitns(ifit(np/2))
         do 15 k=1,n
            write(*,'(22x,3i10)')
     +         nint(ndpwr*oldph(k,ifit(np  ))),
     +         nint(ndpwr*oldph(k,ifit(np-1))),
     +         nint(ndpwr*oldph(k,ifit(np/2)))
   15    continue
 
      endif
      end subroutine

c**********************************************************************
c                         GENETICS MODULE
c**********************************************************************
c
c     ENCODE:    encodes phenotype into genotype
c                called by: PIKAIA
c
c     DECODE:    decodes genotype into phenotype
c                called by: PIKAIA
c
c     CROSS:     Breeds two offspring from two parents
c                called by: PIKAIA
c
c     MUTATE:    Introduces random mutation in a genotype
c                called by: PIKAIA
c
c     ADJMUT:    Implements variable mutation rate
c                called by: PIKAIA
c
c**********************************************************************
      subroutine encode(n,nd,ph,gn)
c======================================================================
c     encode phenotype parameters into integer genotype
c     ph(k) are x,y coordinates [ 0 < x,y < 1 ]
c======================================================================
c
      implicit none
c
c     Inputs:
      integer   n, nd
      real(8)      ph(n)
c
c     Output:
      integer   gn(n*nd)
c
c     Local:
      integer   ip, i, j, ii
      real(8)      z
c
      z=10.**nd
      ii=0
      do 1 i=1,n
         ip=int(ph(i)*z)
         do 2 j=nd,1,-1
            gn(ii+j)=mod(ip,10)
            ip=ip/10
    2   continue
        ii=ii+nd
    1 continue
 
      return
      end subroutine
 
c**********************************************************************
      subroutine decode(n,nd,gn,ph)
c======================================================================
c     decode genotype into phenotype parameters
c     ph(k) are x,y coordinates [ 0 < x,y < 1 ]
c======================================================================
c
      implicit none
c
c     Inputs:
      integer   n, nd, gn(n*nd)
c
c     Output:
      real(8)      ph(n)
c
c     Local:
      integer   ip, i, j, ii
      real(8)      z
c
      z=10.**(-nd)
      ii=0
      do 1 i=1,n
         ip=0
         do 2 j=1,nd
            ip=10*ip+gn(ii+j)
    2    continue
         ph(i)=ip*z
         ii=ii+nd
    1 continue
 
      return
      end subroutine
 
c**********************************************************************
      subroutine cross(n,nd,pcross,gn1,gn2)
c======================================================================
c     breeds two parent chromosomes into two offspring chromosomes
c     breeding occurs through crossover. If the crossover probability
c     test yields true (crossover taking place), either one-point or
c     two-point crossover is used, with equal probabilities.
c
c     Compatibility with version 1.0: To enforce 100% use of one-point
c     crossover, un-comment appropriate line in source code below
c======================================================================
c
      implicit none
c
c     Inputs:
      integer        n, nd
      real(8)           pcross
c
c     Input/Output:
      integer        gn1(n*nd), gn2(n*nd)
c
c     Local:
      integer        i, ispl, ispl2, itmp, t
c
 
 
c     Use crossover probability to decide whether a crossover occurs
      if (urand().lt.pcross) then

c        Compute first crossover point
         ispl=int(urand()*n*nd)+1
 
c        Now choose between one-point and two-point crossover 
         if (urand().lt.0.5) then
            ispl2=n*nd
         else
            ispl2=int(urand()*n*nd)+1
c           Un-comment following line to enforce one-point crossover
c           ispl2=n*nd
            if (ispl2.lt.ispl) then
               itmp=ispl2
               ispl2=ispl
               ispl=itmp
            endif
         endif 
 
c        Swap genes from ispl to ispl2
         do 10 i=ispl,ispl2
            t=gn2(i)
            gn2(i)=gn1(i)
            gn1(i)=t
   10    continue
      endif
 
      return
      end subroutine
 
c**********************************************************************
      subroutine mutate(n,nd,pmut,gn,imut)
c======================================================================
c     Mutations occur at rate pmut at all gene loci
c        imut=1    Uniform mutation, constant rate
c        imut=2    Uniform mutation, variable rate based on fitness
c        imut=3    Uniform mutation, variable rate based on distance
c        imut=4    Uniform or creep mutation, constant rate
c        imut=5    Uniform or creep mutation, variable rate based on
c                  fitness
c        imut=6    Uniform or creep mutation, variable rate based on
c                  distance
c======================================================================
c
      implicit none
c
c     Input:
      integer        n, nd, imut
      real(8)           pmut
c
c     Input/Output:
      integer        gn(n*nd)
c
c     Local:
      integer        i,j,k,l,ist,inc,loc,kk

c
c
c     Decide which type of mutation is to occur
      if(imut.ge.4.and.urand().le.0.5)then

c     CREEP MUTATION OPERATOR
c     Subject each locus to random +/- 1 increment at the rate pmut
         do 1 i=1,n
            do 2 j=1,nd
               if (urand().lt.pmut) then
c     Construct integer
                  loc=(i-1)*nd+j
                  inc=nint ( urand() )*2-1
                  ist=(i-1)*nd+1
                  gn(loc)=gn(loc)+inc
c                 write(*,*) ist,loc,inc
c     This is where we carry over the one (up to two digits)
c     first take care of decrement below 0 case
                  if(inc.lt.0 .and. gn(loc).lt.0)then
                     if(j.eq.1)then
                        gn(loc)=0
                     else
                        do 3 k=loc,ist+1,-1
                           gn(k)=9
                           gn(k-1)=gn(k-1)-1
                        if( gn(k-1).ge.0 )goto 4
    3                   continue
c    we popped under 0.00000 lower bound; fix it up
                        if( gn(ist).lt.0.)then
                           do 5 l=ist,loc
                              gn(l)=0
    5                      continue
                        endif
    4                   continue
                     endif
                  endif
                  if(inc.gt.0 .and. gn(loc).gt.9)then
                     if(j.eq.1)then
                        gn(loc)=9
                     else
                        do 6 k=loc,ist+1,-1
                           gn(k)=0
                           gn(k-1)=gn(k-1)+1
                           if( gn(k-1).le.9 )goto 7
    6                   continue
c                       we popped over 9.99999 upper bound; fix it up
                       if( gn(ist).gt.9 )then
                           do 8 l=ist,loc
                              gn(l)=9
    8                      continue
                        endif
    7                   continue
                     endif
                  endif
               endif
    2       continue
    1    continue

      else

c     UNIFORM MUTATION OPERATOR
c     Subject each locus to random mutation at the rate pmut
         do 10 i=1,n*nd
            if (urand().lt.pmut) then
               gn(i)=int(urand()*10.)
            endif
   10    continue
      endif

      return
      end subroutine
 
c**********************************************************************
      subroutine adjmut(n,np,ifit,pmutmn,pmutmx,
     +                  pmut,imut)
c======================================================================
c     dynamical adjustment of mutation rate;
c        imut=2 or imut=5 : adjustment based on fitness differential
c                           between best and median individuals
c        imut=3 or imut=6 : adjustment based on metric distance
c                           between best and median individuals
c======================================================================
c
      implicit none
c
c     Input:
      integer        n, np, ifit(np), imut
      real(8)        pmutmn, pmutmx
c
c     Input/Output:
      real(8)           pmut
c
c     Local:
      integer        i
      real(8)           rdif, rdiflo, rdifhi, delta
      parameter      (rdiflo=0.05, rdifhi=0.25, delta=1.5)

      if(imut.eq.2.or.imut.eq.5)then
c     Adjustment based on fitness differential 
         rdif=abs(fitns(ifit(np))-fitns(ifit(np/2)))/
     +           (fitns(ifit(np))+fitns(ifit(np/2)))
      else if(imut.eq.3.or.imut.eq.6)then
c     Adjustment based on normalized metric distance
         rdif=0.
         do 1 i=1,n
            rdif=rdif+( oldph(i,ifit(np))-oldph(i,ifit(np/2)) )**2
    1    continue
         rdif=sqrt( rdif ) / float(n)
      endif

      if(rdif.le.rdiflo)then
         pmut=min(pmutmx,pmut*delta)
      else if(rdif.ge.rdifhi)then
         pmut=max(pmutmn,pmut/delta)
      endif
 
      return
      end subroutine


c**********************************************************************
c                       REPRODUCTION MODULE
c**********************************************************************
c
c     SELECT:   Parent selection by roulette wheel algorithm
c               called by: PIKAIA
c
c     RNKPOP:   Ranks initial population
c               called by: PIKAIA, NEWPOP
c
c     GENREP:   Inserts offspring into population, for full
c               generational replacement
c               called by: PIKAIA
c
c     STDREP:   Inserts offspring into population, for steady-state
c               reproduction
c               called by: PIKAIA
c               calls:     FF
c
c     NEWPOP:   Replaces old generation with new generation
c               called by: PIKAIA
c               calls:     FF, RNKPOP
c
c**********************************************************************
      subroutine select(np,jfit,fdif,idad)
c======================================================================
c     Selects a parent from the population, using roulette wheel
c     algorithm with the relative fitnesses of the phenotypes as
c     the "hit" probabilities [see Davis 1991, chap. 1].
c======================================================================
c     USES: urand
      implicit none
c
c     Input:
      integer        np, jfit(np)
      real(8)           fdif
c
c     Output:
      integer        idad
c
c     Local:
      integer        np1, i
      real(8)           dice, rtfit
c
c
c
      np1 = np+1
      dice = urand()*np*np1
      rtfit = 0.
      do 1 i=1,np
         rtfit = rtfit+np1+fdif*(np1-2*jfit(i))
         if (rtfit.ge.dice) then
            idad=i
            goto 2
         endif
    1 continue
c     Assert: loop will never exit by falling through
 
    2 return
      end subroutine
 
c**********************************************************************
      subroutine rnkpop(n,arrin,indx,rank)
c======================================================================
c     Calls external sort routine to produce key index and rank order
c     of input array arrin (which is not altered).
c======================================================================
c     USES: rqsort
      implicit none
c
c     Input
      integer    n
      real(8)       arrin(n)
c
c     Output
      integer    indx(n),rank(n)
c
c     Local
      integer    i
c
c
c
c     Compute the key index
      call rqsort(n,arrin,indx)
c
c     ...and the rank order
      do 1 i=1,n
         rank(indx(i)) = n-i+1
    1 continue
      return
      end subroutine
 
c***********************************************************************
      subroutine genrep(n,np,ip,ph)
c=======================================================================
c     full generational replacement: accumulate offspring into new
c     population array
c=======================================================================
c
      implicit none
 
c     Input:
      integer        n, np, ip
      real(8)           ph(n,2)
c
c
c     Local:
      integer        i1, i2, k
c
c
c     Insert one offspring pair into new population
      i1=2*ip-1
      i2=i1+1
      do 1 k=1,n
         newph(k,i1)=ph(k,1)
         newph(k,i2)=ph(k,2)
    1 continue
 
      return
      end subroutine
 
c**********************************************************************
      subroutine stdrep
     +   (ff,n,np,irep,ielite,ph,ifit,jfit,nnew)
c======================================================================
c     steady-state reproduction: insert offspring pair into population
c     only if they are fit enough (replace-random if irep=2 or
c     replace-worst if irep=3).
c======================================================================
c     USES: ff, urand
      implicit none
c
c     Input:
      integer        n, np, irep, ielite
      real(8)        ff, ph(n,2)
      external       ff
c
c     Input/Output:
      integer        ifit(np), jfit(np)
c
c     Output:
      integer        nnew
 
c     Local:
      integer        i, j, k, i1, if1
      real(8)           fit
c
c
cf2py real(8) f, x(n)
cf2py f = ff(n, x)
c
      nnew = 0
      do 1 j=1,2
 
c        1. compute offspring fitness (with caller's fitness function)
         fit=ff(n,ph(1,j))
 
c        2. if fit enough, insert in population
         do 20 i=np,1,-1
            if (fit.gt.fitns(ifit(i))) then
 
c              make sure the phenotype is not already in the population
               if (i.lt.np) then
                  do 5 k=1,n
                     if (oldph(k,ifit(i+1)).ne.ph(k,j)) goto 6
    5             continue
                  goto 1
    6             continue
               endif
 
c              offspring is fit enough for insertion, and is unique
 
c              (i) insert phenotype at appropriate place in population
               if (irep.eq.3) then
                  i1=1
               else if (ielite.eq.0 .or. i.eq.np) then
                  i1=int(urand()*np)+1
               else
                  i1=int(urand()*(np-1))+1
               endif
               if1 = ifit(i1)
               fitns(if1)=fit
               do 21 k=1,n
                  oldph(k,if1)=ph(k,j)
   21          continue
 
c              (ii) shift and update ranking arrays
               if (i.lt.i1) then
 
c                 shift up
                  jfit(if1)=np-i
                  do 22 k=i1-1,i+1,-1
                     jfit(ifit(k))=jfit(ifit(k))-1
                     ifit(k+1)=ifit(k)
   22             continue
                  ifit(i+1)=if1
               else
 
c                 shift down
                  jfit(if1)=np-i+1
                  do 23 k=i1+1,i
                     jfit(ifit(k))=jfit(ifit(k))+1
                     ifit(k-1)=ifit(k)
   23             continue
                  ifit(i)=if1
               endif
               nnew = nnew+1
               goto 1
            endif
   20    continue
 
    1 continue
 
      return
      end subroutine
 
c**********************************************************************
      subroutine newpop
     +   (ff,ielite,n,np,ifit,jfit,nnew)
c======================================================================
c     replaces old population by new; recomputes fitnesses & ranks
c======================================================================
c     USES: ff, rnkpop
      implicit none
c
c     Input:
      integer        np, n, ielite
      real(8)        ff
      external       ff
c
c     Output:
      integer        ifit(np), jfit(np), nnew
c
c     Local:
      integer        i, k
c
cf2py real(8) f, x(n)
cf2py f = ff(n, x)
c
      nnew = np
 
c     if using elitism, introduce in new population fittest of old
c     population (if greater than fitness of the individual it is
c     to replace)
      if (ielite.eq.1 .and. ff(n,newph(1,1)).lt.fitns(ifit(np))) then
         do 1 k=1,n
            newph(k,1)=oldph(k,ifit(np))
    1    continue
         nnew = nnew-1
      endif
 
c     replace population
      do 2 i=1,np
         do 3 k=1,n
            oldph(k,i)=newph(k,i)
    3    continue
 
c        get fitness using caller's fitness function
         fitns(i)=ff(n,oldph(1,i))
    2 continue
 
c     compute new population fitness rank order
      call rnkpop(np,fitns,ifit,jfit)
 
      return
      end subroutine

c**********************************************************************
c     Utility routines
c
c**********************************************************************
      subroutine rqsort(n,a,p)
c======================================================================
c     Return integer array p which indexes array a in increasing order.
c     Array a is not disturbed.  The Quicksort algorithm is used.
c
c     B. G. Knapp, 86/12/23
c
c     Reference: N. Wirth, Algorithms and Data Structures,
c     Prentice-Hall, 1986
c======================================================================
      implicit none

c     Input:
      integer   n
      real(8)   a(n)

c     Output:
      integer   p(n)

c     Constants
      integer   LGN, Q
      parameter (LGN=32, Q=11)
c        (LGN = log base 2 of maximum n;
c         Q = smallest subfile to use quicksort on)

c     Local:
      real(8)   x
      integer   stackl(LGN),stackr(LGN),s,t,l,m,r,i,j
c
cf2py intent(in) a
cf2py intent(out) p
c
c     Initialize the stack
      stackl(1)=1
      stackr(1)=n
      s=1

c     Initialize the pointer array
      do 1 i=1,n
         p(i)=i
    1 continue

    2 if (s.gt.0) then
         l=stackl(s)
         r=stackr(s)
         s=s-1

    3    if ((r-l).lt.Q) then

c           Use straight insertion
            do 6 i=l+1,r
               t = p(i)
               x = a(t)
               do 4 j=i-1,l,-1
                  if (a(p(j)).le.x) goto 5
                  p(j+1) = p(j)
    4          continue
               j=l-1
    5          p(j+1) = t
    6       continue
         else

c           Use quicksort, with pivot as median of a(l), a(m), a(r)
            m=(l+r)/2
            t=p(m)
            if (a(t).lt.a(p(l))) then
               p(m)=p(l)
               p(l)=t
               t=p(m)
            endif
            if (a(t).gt.a(p(r))) then
               p(m)=p(r)
               p(r)=t
               t=p(m)
               if (a(t).lt.a(p(l))) then
                  p(m)=p(l)
                  p(l)=t
                  t=p(m)
               endif
            endif

c           Partition
            x=a(t)
            i=l+1
            j=r-1
    7       if (i.le.j) then
    8          if (a(p(i)).lt.x) then
                  i=i+1
                  goto 8
               endif
    9          if (x.lt.a(p(j))) then
                  j=j-1
                  goto 9
               endif
               if (i.le.j) then
                  t=p(i)
                  p(i)=p(j)
                  p(j)=t
                  i=i+1
                  j=j-1
               endif
               goto 7
            endif

c           Stack the larger subfile
            s=s+1
            if ((j-l).gt.(r-i)) then
               stackl(s)=l
               stackr(s)=j
               l=i
            else
               stackl(s)=i
               stackr(s)=r
               r=j
            endif
            goto 3
         endif
         goto 2
      endif
      return
      end subroutine
c**********************************************************************
      function urand()
c=====================================================================
c     Return the next pseudo-random deviate from a sequence which is
c     uniformly distributed in the interval [0,1]
c
c     Uses the function ran0, the "minimal standard" random number
c     generator of Park and Miller (Comm. ACM 31, 1192-1201, Oct 1988;
c     Comm. ACM 36 No. 7, 105-110, July 1993).
c=====================================================================
      implicit none
c
c     Input - none
c
c     Output
      real(8)     urand
c
c     Local
c      integer  iseed
c      real     ran0
c      external ran0
c
c     Common block to make iseed visible to rninit (and to save
c     it between calls)
c      common /rnseed/ iseed
c
c      urand = ran0( iseed )

      call random_number(urand) 

      return
      end function
c*********************************************************************

      end module pypikaia 
