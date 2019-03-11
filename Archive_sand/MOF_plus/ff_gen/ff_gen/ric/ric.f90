module red_int_coords

  !use ric_common
  use ric_stretches
  use ric_in_bends
  use ric_out_bends
  use ric_lin_bends
  use ric_torsions
  use ric_eckart

  implicit none

  ! Type of floating-point numbers
  ! Note: this cannot be used from ric_common as pyf don't support module constants
  integer, parameter :: dp = selected_real_kind(15,307)

  ! System description in the Cartesian coordinates
  real(dp)              :: cart_hmat(3,3)    ! h-matrix (lattice vectors in columns)
  real(dp), allocatable :: cart_coords(:,:)  ! atomic coordinates (3, # of atoms)
  real(dp), allocatable :: cart_hessian(:,:) ! Hessian matrix (3 x # of atoms, 3 x # of atoms)

  real(dp), allocatable :: atomic_masses(:)  ! (# of atoms)

  ! Definitions of the redundant internal coordinates (ric = redundant internal coordinates)
  ! This replace
  !  common /config/ ibb(noi)
  !  common /atomtors/ ijtors
  ! as it is hard-follow and bug-prone.
  integer , allocatable :: ric_def_stretches(:,:)  ! definitions of bond stretches (2, # of stretches)
  integer , allocatable :: ric_def_in_bends(:,:)   ! definitions of in-plane bends (3, # of in-plane bends)
  integer , allocatable :: ric_def_out_bends(:,:)  ! definitions of out-of-plane bends (4, # of out-of-plane bends)
  integer , allocatable :: ric_def_lin_bends(:,:)  ! definitions of linear bends (3, # of linear bends)
  integer , allocatable :: ric_def_torsions(:,:)   ! definitions of torsion (12, # of torsions)
  integer , allocatable :: ric_def_eckart_trans(:) ! enable Eckart translations coords
  integer , allocatable :: ric_def_eckart_rots(:)  ! enable Eckart rotation coords

  ! Linear bend reference atom indices and axes
  integer , allocatable :: ric_lin_bend_inds(:)   ! (# of linear bends)
  real(dp), allocatable :: ric_lin_bend_axes(:,:) ! (3, # of linear bends)

  ! Torsion dihedral angle definition for value calculations
  integer, allocatable :: ric_torsion_ivals(:,:) ! (2, # of torsions)

  ! B-matrix indices of redundant internal coordinates (ibr = indices of B-matrix rows)
  ! This will allow to keep the ordering of the coordinates as in the old code for comparison, but it completely
  ! decouple evaluation of B-matrix elements allowing OpenMP parallelization.
  integer , allocatable :: ric_ibr_stretches(:)    ! bonds strech indices (# of stretches)
  integer , allocatable :: ric_ibr_in_bends(:)     ! in-plane bend indices (# of in-plane bends)
  integer , allocatable :: ric_ibr_out_bends(:)    ! out-of-plane indices (# of out-of-plane bends)
  integer , allocatable :: ric_ibr_lin_bends(:)    ! linear bend indices (# of linear bends)
  integer , allocatable :: ric_ibr_torsions(:)     ! torsions indices (# of torsions)
  integer , allocatable :: ric_ibr_eckart_trans(:) ! Eckart rotation indices
  integer , allocatable :: ric_ibr_eckart_rots(:)  ! Eckart translation indices

  ! Relation between the Cartesian and redundant internal coordinates
  real(dp), allocatable :: bmat(:,:)     ! B-matrix (cart => ric) (# of ric coords, 3 x # of atoms)
  real(dp), allocatable :: bmat_inv(:,:) ! Inverse B-matrix (ric => cart) (3 x # of atoms, # of ric coords)
  integer               :: rank          ! Rank of B-matrix

  ! System description in the redundant internal coordinates
  real(dp), allocatable :: ric_val_stretches(:)
  real(dp), allocatable :: ric_val_in_bends(:)
  real(dp), allocatable :: ric_val_out_bends(:)
  real(dp), allocatable :: ric_val_lin_bends(:)
  real(dp), allocatable :: ric_val_torsions(:)
  real(dp), allocatable :: ric_val_eckarts(:)
  real(dp), allocatable :: ric_hessian(:,:) ! (# of ric coords, # of ric coords)

contains

  integer function bmat_construct()

    integer :: stat

    ! stuff from ff_gen_internal.f, ff_gen_vector.f, etc., split into comprehensible size functions

    if (allocated(ric_def_stretches)) &
      stat = ric_stretches_bmat(cart_hmat,cart_coords,ric_def_stretches, &
                                ric_ibr_stretches,bmat,ric_val_stretches)

    if (allocated(ric_def_in_bends)) &
      stat = ric_in_bends_bmat(cart_hmat,cart_coords,ric_def_in_bends, &
                               ric_ibr_in_bends,bmat,ric_val_in_bends)

    if (allocated(ric_def_out_bends)) &
      stat = ric_out_bends_bmat(cart_hmat,cart_coords,ric_def_out_bends, &
                                ric_ibr_out_bends,bmat,ric_val_out_bends)

    if (allocated(ric_def_lin_bends)) &
      stat = ric_lin_bends_bmat(cart_hmat,cart_coords,ric_def_lin_bends,  &
                                ric_lin_bend_inds, ric_lin_bend_axes,     &
                                ric_ibr_lin_bends,bmat,ric_val_lin_bends)

    if (allocated(ric_def_torsions)) &
      stat = ric_torsions_bmat(cart_hmat,cart_coords,ric_def_torsions,  &
                               ric_torsion_ivals,ric_ibr_torsions,bmat, &
                               ric_val_torsions)

    ! Evaluate B-matrix elements for Eckart translation coordinates
    if (allocated(ric_def_eckart_trans)) &
      stat = ric_eckart_trans_bmat(cart_coords,atomic_masses,ric_def_eckart_trans, &
                                   ric_ibr_eckart_trans,bmat,ric_val_eckarts)

    ! Evaluate B-matrix elements for Eckart rotation coordinates
    if (allocated(ric_def_eckart_rots)) &
      stat = ric_eckart_rot_bmat(cart_coords,atomic_masses,ric_def_eckart_rots, &
                                 ric_ibr_eckart_rots,bmat,ric_val_eckarts)

    bmat_construct = stat

  end function

  integer function bmat_invert()

    integer :: nric, ncart, i, nwork, info
    real(dp), allocatable :: bmat_(:,:), rhs(:,:), sing(:), work(:)

    ! stuff from ff_gen_b_invert.f

!     Make a copy of B-matrix, since dgelss overwrites it.
!     A = b

    ! Get numbers
    ncart = size(bmat,dim=1) ! # of redundant internat coordinates
    nric  = size(bmat,dim=2) ! # of Cartesian coordinates

    ! Make a copy of B-matrix, since dgelss overwrites it
    allocate(bmat_(nric,ncart))
    bmat_ = transpose(bmat) ! HACK

!     Construct the left-hand-side matrix as a diagonal matrix
!      C = 0.0d0
!      do i = 1,ir
!        C(i,i) = 1.0d0
!      end do

    ! Construct right hand side matrix
    allocate(rhs(nric,nric))
    rhs = 0._dp
    forall(i=1:nric) rhs(i,i) = 1._dp

    ! Allocate singular value and work arrays
    allocate(sing(ncart))
    nwork = 3*ncart + max(2*ncart,nric)
    allocate(work(nwork))

!     Invert B-matrix.
!      call dgelss(ir,k,ir,A,ms,C,ms,S,1.0e-10,rank,work,ms*ms,
!      *           info)

!   Compute pseudo inverse of B-matrix
    call dgelss(nric  , & ! M
                ncart , & ! N
                nric  , & ! NRHS = M
                bmat_ , & ! A(LDA,N) = A(M,N)
                nric  , & ! LDA = M
                rhs   , & ! B(LDB,NRHS) = B(M,M)
                nric  , & ! LDB = M
                sing  , & ! S(N)
                -1._dp, & ! RCOND
                rank  , & ! RANK
                work  , & ! WORK(NWORK)
                nwork , & ! NWORK
                info    ) ! INFO
    bmat_inv = rhs(1:ncart,:)

!     Inverted and transposed B-matrix is in the vv array
!      vv(1:ir,1:k) = transpose(C(1:k,1:ir))

!      if (info .eq. 0) then
!        write(iout,'(a,/)') 'done!'
!      else
!        write(iout,'(a,i5,/)') 'dgelss error: ', info
!        idm = 0
!        return
!      end if

!      Check if B-matrix is rank-defincient
!       if (rank .lt. ir) idm = 0

    bmat_invert = info

  end function

  integer function hessian_project()

    integer :: nric, ncart, stat
    real(dp), allocatable :: tmp(:,:)

    ! stuff from ff_gen_trans_hess.f

    ! Get numbers
    ncart = size(bmat,dim=1) ! # of redundant internat coordinates
    nric  = size(bmat,dim=2) ! # of Cartesian coordinates

    ! Allocate a temporary array
    allocate(tmp(ncart,nric),stat=stat)


!     Multiply the Hessian matrix and the inverted B-matrix
!     taking into accout that the former is symmetrix.
!      call dsymm('L','U',ir,k,1.0d0,v,ms,vv,ms,0.0d0,tmp,ms)

!     Multilpy the transposed inverted B-matrix and the matrix
!     from previous operations.
!      call dgemm('T','N',k,k,ir,1.0d0,vv,ms,tmp,ms,0.0d0,v,ms)

!   Multiply Cartesian Hessian matrix and the B-matrix
!   taking into accout that the former is symmetrix
    call dsymm('L'         , & ! SIDE
               'U'         , & ! UPLO
               ncart       , & ! M
               nric        , & ! N
               1._dp       , & ! ALPHA
               cart_hessian, & ! A(LDA,ka) = A(M,M)
               ncart       , & ! LDA = M
               bmat_inv    , & ! B(LDB,N) = B(M,N)
               ncart       , & ! LDB = M
               0._dp       , & ! BETA
               tmp         , & ! C(LDC,N) = C(M,N)
               ncart         ) ! LDC = M

!   Multilpy the transposed inverted B-matrix and the matrix
!   from previous operations
    call dgemm('T'        , & ! TRANSA
               'N'        , & ! TRANSB
               nric       , & ! M
               nric       , & ! N
               ncart      , & ! K
               1._dp      , & ! ALPHA
               bmat_inv   , & ! A(LDA,ka) = A(K,M)
               ncart      , & ! LDA = K
               tmp        , & ! B(LDB,kb) = B(K,N)
               ncart      , & ! LDB = K
               0._dp      , & ! BETA
               ric_hessian, & ! C(LDC,N) = C(M,N)
               nric         ) ! LDC = M

    hessian_project = stat

  end function

end module

