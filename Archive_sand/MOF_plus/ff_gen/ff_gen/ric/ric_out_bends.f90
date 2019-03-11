module ric_out_bends

  use ric_common
  use ric_math

  implicit none
  private

  public :: ric_out_bends_bmat

contains

  integer function ric_out_bends_bmat(hmat,coords,defs,ibrs,bmat,vals)
    real(dp), intent(in)    :: hmat(3,3)
    real(dp), intent(in)    :: coords(:,:)
    integer , intent(in)    :: defs(:,:)
    integer , intent(in)    :: ibrs(:)
    real(dp), intent(inout) :: bmat(:,:)
    real(dp), intent(out)   :: vals(:)

    integer :: nout_bend, i

    ! Get the number of out-of-plane bends
    nout_bend = size(defs,dim=2)

    do i = 1, nout_bend
      call ric_out_bend(hmat,coords,defs(:,i),bmat(:,ibrs(i)),vals(i))
    end do

    ric_out_bends_bmat = 0

  end function

!  TODO: DEBUG!!!
!  subroutine ric_out_bend_new(hmat,coords,def,bmat,val)
!    real(dp), intent(in)    :: hmat(3,3)
!    real(dp), intent(in)    :: coords(:,:)
!    integer , intent(in)    :: def(4)
!    real(dp), intent(inout) :: bmat(:)
!    real(dp), intent(out)   :: val
!
!    integer  :: ia1(3), ia2(3), ia3(3), ia4(3), i
!    real(dp) :: a1(3), a2(3), a3(3), a4(3),      &
!                v12(3), v13(3), v14(3), d12, d13, d14, &
!                c(3), vc(3,4), dc, dc1, dc2, dc3, dc4, &
!                n(3), d, ax1(3), ax2(3), am(3), lm(3)
!
!    ! Get atom positions
!    a1 = coords(:,def(1)) ! Central atom
!    a2 = coords(:,def(2))
!    a3 = coords(:,def(3))
!    a4 = coords(:,def(4))
!
!    ! Compute vectors to the central atom
!    call rel_vec(hmat,a1,a2,norm_vec=v12,dist=d12) ! 1 --> 2
!    call rel_vec(hmat,a1,a3,norm_vec=v13,dist=d13) ! 1 --> 3
!    call rel_vec(hmat,a1,a4,norm_vec=v14,dist=d14) ! 1 --> 4
!
!    ! TODO: check planar angles between vectors
!
!    ! Compute plane normal vector
!    n = normalize(cross_prod(v13-v12,v14-v13))
!
!    ! Compute out-of-plane angle value
!    ! Note: all terminal vectors {v12,v12,v14} are normalized,
!    !       so the angles with the normal vector, n, are the same.
!    val = dot_product(n,v12)
!    val = min(max(val,-1._dp),1._dp)
!    val = acos(val) - acos(0._dp) ! acos(0._dp) = pi
!
!    ! Compute B-matrix indices
!    ia1 = (/( 3*(def(1)-1)+i, i=1,3 )/)
!    ia2 = (/( 3*(def(2)-1)+i, i=1,3 )/)
!    ia3 = (/( 3*(def(3)-1)+i, i=1,3 )/)
!    ia4 = (/( 3*(def(4)-1)+i, i=1,3 )/)
!
!    ! Compute B-matrix elements (with angular momenta)
!    bmat(ia2) = normalize(cross_prod(cross_prod(n,v12),v12))/(3._dp*d12)
!    bmat(ia3) = normalize(cross_prod(cross_prod(n,v13),v13))/(3._dp*d13)
!    bmat(ia4) = normalize(cross_prod(cross_prod(n,v14),v14))/(3._dp*d14)
!    bmat(ia1) = -bmat(ia2) - bmat(ia3) - bmat(ia4)

!    ! Compute vectors to the centre
!    c = 0.25_dp*(a1 + a2 + a3 + a4)
!    call rel_vec(hmat,c,a1,vec=vc(:,1),dist=dc1) ! c --> 1
!    call rel_vec(hmat,c,a2,vec=vc(:,2),dist=dc2) ! c --> 2
!    call rel_vec(hmat,c,a3,vec=vc(:,3),dist=dc3) ! c --> 3
!    call rel_vec(hmat,c,a4,vec=vc(:,4),dist=dc4) ! c --> 4
!
!    ! Compute angular momenta
!    lm = bmat(ia1) +  bmat(ia2) + bmat(ia3) + bmat(ia4)
!    am = cross_prod(bmat(ia1),vc(:,1)) + cross_prod(bmat(ia2),vc(:,2)) &
!       + cross_prod(bmat(ia3),vc(:,3)) + cross_prod(bmat(ia4),vc(:,4))
!    write(*,*) 'lm', lm
!    write(*,*) 'am', am
!
!    !dc = dc1 + dc2 + dc3 + dc4
!    !bmat(ia1) = bmat(ia1) + cross_prod(am,vc(:,1))*dc1/dc
!    !bmat(ia2) = bmat(ia2) + cross_prod(am,vc(:,2))*dc2/dc
!    !bmat(ia3) = bmat(ia3) + cross_prod(am,vc(:,3))*dc3/dc
!    !bmat(ia4) = bmat(ia4) + cross_prod(am,vc(:,4))*dc4/dc
!
!    !bmat(ia1) = bmat(ia1) - .25_dp*lm
!    !bmat(ia2) = bmat(ia2) - .25_dp*lm
!    !bmat(ia3) = bmat(ia3) - .25_dp*lm
!    !bmat(ia4) = bmat(ia4) - .25_dp*lm
!
!    lm = bmat(ia1) +  bmat(ia2) + bmat(ia3) + bmat(ia4)
!    am = cross_prod(bmat(ia1),vc(:,1)) + cross_prod(bmat(ia2),vc(:,2)) &
!       + cross_prod(bmat(ia3),vc(:,3)) + cross_prod(bmat(ia4),vc(:,4))
!    write(*,*) 'lm', lm
!    write(*,*) 'am', am
!
!  end subroutine

  ! This is from ff_gen_internals.f (wagin4)
  subroutine ric_out_bend(hmat,coords,def,bmat,val)
    real(dp), intent(in)    :: hmat(3,3)
    real(dp), intent(in)    :: coords(:,:)
    integer , intent(in)    :: def(4)
    real(dp), intent(inout) :: bmat(:)
    real(dp), intent(out)   :: val


    integer  :: i, ia1(3), ia2(3), ia3(3), ia4(3)
    real(dp) :: a1(3) , a2(3) , a3(3) , a4(3) , cc(3), &
                a1_(3), a2_(3), a3_(3), &
                e41(3) , e42(3) , e43(3) , os23(3), &
                e41_(3), e42_(3), &
                pn, r41 , r42, r43, r41_, r1, ss, c, s

    ! Get atom positions
    a1 = coords(:,def(2))
    a2 = coords(:,def(3))
    a3 = coords(:,def(4))
    a4 = coords(:,def(1)) ! Central atom

    call rel_vec(hmat,a4,a1,norm_vec=e41,dist=r41)
    call rel_vec(hmat,a4,a2,norm_vec=e42,dist=r42)
    call rel_vec(hmat,a4,a3,norm_vec=e43,dist=r43)

    pn   = sqrt(max(1._dp - cdc(hmat,a3,a4,a2)**2,0._dp))
    os23 = cross_prod(e42,e43)
    ss   = sum(e41*os23)
    val  = acos(ss/pn)
    cc   = e41 - os23*ss/pn**2 + a4

    call rel_vec(hmat,a4,cc,norm_vec=e41_,dist=r41_)

    a1_ = cross_prod(e41_,e42 )
    a2_ = cross_prod(e43 ,e41_)
    a3_ = a4 + os23

    call rel_vec(hmat,a4,a1 ,norm_vec=e41_,dist=r1)
    call rel_vec(hmat,a4,a3_,norm_vec=e42_)
    c = cdc(hmat,a1,a4,a3_)
    s = sqrt(max(1._dp - c**2,0._dp))

    ! Compute B-matrix indices
    ia1 = (/( 3*(def(2)-1)+i, i=1,3 )/)
    ia2 = (/( 3*(def(3)-1)+i, i=1,3 )/)
    ia3 = (/( 3*(def(4)-1)+i, i=1,3 )/)
    ia4 = (/( 3*(def(1)-1)+i, i=1,3 )/)

    ! Compute B-matrix elements
    bmat(ia1) = (c*e41_ - e42_)/(r1*s)
    bmat(ia2) = -a2_/(r42*pn)
    bmat(ia3) = -a1_/(r43*pn)
    bmat(ia4) = -bmat(ia1) - bmat(ia2) - bmat(ia3)

  end subroutine

end module

