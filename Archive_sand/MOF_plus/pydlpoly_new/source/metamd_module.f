      module metamd_module

c***********************************************************************
c     
c     pydlpoly module for metamd of self diffusion
c     copyright - RUB
c     author    - Michael Jacob    2012/2013
c                 Rochus Schmid    2013 (revisions)
c     
c     This code has initially been written by Michael Jacob
c       starting with a version using analytic gaussians. ater he converted to
c       using a numeric grid (first with scipy lin interpolation for energy and force -> not energy conserving)
c       Now the code is designed to use tricubic as a 3D interpolator (isotropic!!) => no force array needed any more 
c     
c***********************************************************************

      
      implicit none
c     here we could have some module variables      

      contains

      
      subroutine drop_energy(grid,pos,W,start,resolution,celllength,
     &                    gwidthfact,ncalc,startcalc,energylist)
        
        implicit none
        integer, intent(in) :: grid
        real*8, dimension(3), intent(in)                 :: pos
        real*8, intent(in)                               :: W
        real*8, dimension(3), intent(in)                 :: start
        real*8, dimension(3), intent(in)                 :: resolution
        real*8, dimension(3), intent(in)                 :: celllength
        real*8, intent(in)                               :: gwidthfact
        integer, intent(in)                              :: ncalc
        integer, intent(in)                              :: startcalc
        real*8, dimension(grid,grid,grid), intent(inout) :: energylist
        
        real*8, dimension(3) :: thispos
        real*8, dimension(3) :: d
        real*8  :: d2
        integer :: x, y, z
        integer :: j
        integer :: myx
        
        do z=1,grid
            do y=1,grid
                do x=1,grid
                    energylist(z,y,x)=0.0d0
                end do
            end do
        end do
        
        do z=1,grid
            do y=1,grid
                do x=1,ncalc
                    myx=x+startcalc
                    if (myx .LE. grid) then
                        thispos(1)=(myx-1)*resolution(1)
                        thispos(2)=(y-1)*resolution(2)
                        thispos(3)=(z-1)*resolution(3)
                        d2 = 0.0d0
                        do j=1,3
                            d(j)=thispos(j)-pos(j)
                            d(j)=d(j)-dnint(d(j)/celllength(j))*
     &                                             celllength(j)
                            d2 = d2 +(d(j)*d(j))
                        end do
                        energylist(z,y,myx)=exp(-d2/gwidthfact)*W
                    end if
                end do
            end do
        end do
      end subroutine drop_energy
      
      

      end module metamd_module