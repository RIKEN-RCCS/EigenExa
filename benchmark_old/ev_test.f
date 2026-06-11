      subroutine ev_test(N, NVEC, A, lda, w, z, ldz, msolver, mode)
      use MPI
!     $    use OMP_LIB
      use eigen_libs_mod
      use eigen_blacs_mod
      use, intrinsic :: ieee_arithmetic
      implicit none

      integer,intent(in) :: N, NVEC, msolver
      real(8),intent(in) :: w(*)
      real(8),intent(inout) :: z(*)
      real(8),intent(inout) :: A(*)
      integer,intent(in) :: lda,ldz
      character(*),intent(in) :: mode

      real(8),parameter :: one = 1.0D0
      real(8),parameter :: zero = 0.0D0

      integer :: NPROW,NPCOL,larray
      integer :: ictxt
      integer :: iam1,nprocs
      integer :: myrow1,mycol1
      integer :: numr_row,numr_col
      integer :: ierr,inod
      integer :: i,j
      real(8) :: s_time,e_time
      integer :: nm, nx, NB

      integer :: numroc,iceil,indxl2g
      external :: numroc,iceil,indxl2g
      real(8) :: pdlange
      external :: pdlange

      real(8) :: err1, err2, rr, m_epsilon
      logical :: err_flag(2)
      integer :: gcol
      real(8) :: DUMWORK
      real(8), allocatable :: Y(:)
      real(8) :: buf(6)

      integer, parameter :: DESC_DIM = 9
      integer, dimension(DESC_DIM) :: a0desc, z0desc
      integer, dimension(DESC_DIM) :: a1desc, z1desc

      integer :: i_1,i_2,i_3
      integer :: j_1,j_2,j_3
      integer :: iam, myrow, mycol      
      integer(8) :: idx

      m_epsilon = get_constant_eps()
      if (msolver == 2) then
         m_epsilon = get_constant_eps_f32()
      endif

      call eigen_get_id   ( iam, myrow, mycol )            
      call eigen_get_procs( nprocs, NPROW, NPCOL )            
      inod=iam-1; myrow1=myrow-1; mycol1=mycol-1

      call eigen_get_matdims( n, nm, nx )

      larray = nm * nx
*     
      NB = eigen_NB
      nm = (((n-1)/NB)/NPROW+1)*NB
      nx = (((NVEC-1)/NB)/NPCOL+1)*NB

      larray = max( larray, nm * nx )

      ictxt = eigen_get_blacs_context()
      call descinit( a0desc, n, n, 1, 1, 0, 0, ictxt, lda, ierr )
      call descinit( z0desc, n, n, 1, 1, 0, 0, ictxt, ldz, ierr )
      call descinit( a1desc, n, n, NB, NB, 0, 0, ictxt, nm, ierr )
      call descinit( z1desc, n, n, NB, NB, 0, 0, ictxt, nm, ierr )
*
* local buffers are elimibnated by sharing other arrays, which
* are allocated largely according to the get_matdims mode='O'.
*
      allocate(Y(larray))
!
! a ==> Y
! z ==> A
!
      call pdgemr2d( n, n, a, 1, 1, a0desc,
     &     Y, 1, 1, a1desc, ictxt)
      call pdgemr2d( n, n, z, 1, 1, z0desc,
     &     A, 1, 1, z1desc, ictxt)

! 0. fill zero in out of index
!
!$OMP PARALLEL DO PRIVATE(i,i_1,j,j_1,idx)
      do i_1=1, nx
         i = (MOD(i_1-1,NB)+1)+(((i_1-1)/NB+1-1)*NPCOL+MYCOL-1)*NB
         if ( i <= nvec ) then
            do j_1=1,nm
               j = (MOD(j_1-1,NB)+1)+(((j_1-1)/NB+1-1)*NPROW+MYROW-1)*NB
               idx = j_1+INT(i_1-1,8)*nm
               if ( j > n ) then
                  Y(idx) = ZERO
                  A(idx) = ZERO
               else
                  if ( ieee_is_nan( Y(idx) ) ) then
                     print*,"???? NAN Y",i_1,i,j_1,j
                  end if
                  if ( ieee_is_nan( A(idx) ) ) then
                     print*,"???? NAN A",i_1,i,j_1,j
                  end if
               end if
            end do
         else
           idx = INT(i_1-1,8)*nm
           Y(1+idx:nm+idx) = ZERO
           A(1+idx:nm+idx) = ZERO
         end if
           Z(1+idx:nm+idx) = ZERO
      end do
!$OMP END PARALLEL DO

! 1. |AZ-ZW|
!    Z = YA = az
      call pdgemm('N', 'N', n, NVEC, n,
     &     one,
     &     Y, 1, 1, a1desc,
     &     A, 1, 1, z1desc,
     &     zero,
     &     Z, 1, 1, z1desc)

!    rr = \| Y \|_{F} = \| a \|_{F}
      rr = Fnorm2_local( n, n, Y(1:), nm, NB )
      deallocate( Y )

!    Z = Z - WA = az - zw
!$OMP PARALLEL DO PRIVATE(i,j_1,j_2,idx)
      do i_1=1, nx
         j_1 = MOD(i_1-1,NB)+1
         j_2 = (i_1-1)/NB*NPCOL+(MYCOL-1)
         i   = j_1 + j_2*NB
         idx = 1 + INT(i_1-1,8)*nm
         call DAXPY(nm, -w(i), A(idx), 1, Z(idx), 1)
      end do
!$OMP END PARALLEL DO
!    err1 = \| Z \|_{F} = \| az-zw \|_{F}
      err1 = Fnorm2_local( n, NVEC, Z, nm, NB )

! 2. |ZZ-I|
!    Z = AA = zz
      call pdgemm('T', 'N', NVEC, NVEC, n,
     &     one,
     &     A, 1, 1, z1desc,
     &     A, 1, 1, z1desc,
     &     zero,
     &     Z, 1, 1, z1desc)

!    Z = Z - eye = zz - eye
!$OMP PARALLEL DO PRIVATE(i,j,j_1,j_2,idx)
      do i_1=1,nx
         j_1 = MOD(i_1-1,NB)+1
         j_2 = (i_1-1)/NB*NPCOL+(MYCOL-1)
         j   = MOD(j_2,NPROW)+1
         if (j == MYROW) then
            idx = (j_1 + (j_2/NPROW)*NB) + INT(i_1-1,8)*nm
            Z(idx) = Z(idx) - one
         end if
      end do
!$OMP END PARALLEL DO
!    err2 = \| Z \|_{F} = \| zz-eye \|_{F}
      err2 = Fnorm2_local( NVEC, NVEC, Z, nm, NB )
!      err2 = pdlange( 'F', NVEC, NVEC, Z, 1, 1, z1desc, DUMWORK )

      buf(1) = rr
      buf(2) = err1
      buf(3) = err2
      call MPI_Allreduce(buf(1), buf(4), 3, MPI_DOUBLE_PRECISION,
     &       MPI_SUM, MPI_COMM_WORLD, ierr)
      rr   = sqrt(buf(4))
      err1 = sqrt(buf(5))
      err2 = sqrt(buf(6))

      if ( inod == 0 ) then

         err_flag(1:2) = .false.

         if ( mode == 'X' .or. mode == 'A' ) then
         rr = err1/(dble(n)*rr*m_epsilon)
         print*,'|AZ-ZW|_{F}/Ne|A|_{F}=', rr
         if ( rr < 768D0 ) then ! '768D0' is a magic number
          print*,"*** Residual Error Test ***   : PASSED"
         else
          print*,"***=============================******"
          print*,"*** Residual Error Test ***   : FAILED"
          print*,"***=============================******"
          if ( rr > 1024 ) err_flag(1) = .true.
         end if
         end if

         if ( mode == 'X' .or. mode == 'A' .or.
     &        mode == 'S' .or. mode == 'T' .or. mode == 'R' ) then
         rr = err2/(dble(n)*m_epsilon)
         print*,'|ZZ-I|_{F}/Ne=        ', rr
         if ( rr < 8D0 ) then ! '8D0' is a magic number
          print*,"*** Orthogonality  Test ***   : PASSED"
         else
          print*,"***=============================******"
          print*,"*** Orthogonality  Test ***   : FAILED"
          print*,"***=============================******"
          if ( rr > 32 ) err_flag(2) = .true.
         end if
         end if

         call flush(6)

         do i=2,2
          err_flag(1) = err_flag(1) .or. err_flag(i)
         end do
         if ( err_flag(1) ) then
            call MPI_Abort( MPI_COMM_WORLD, 1, ierr )
         end if

      end if

      return
      contains

      real(8) function Fnorm2_local( M, N, A, LDA, NB )
      use MPI
!$    use OMP_LIB
      integer, intent(IN) :: M, N, LDA, NB
      real(8), intent(IN) :: A(1:LDA,1:*)
!
      integer :: lx, ly, i, j, nth, i_1, j_1
      real(8) :: r
      real(8), pointer :: rr(:)
!
      ly = (((N-1)/NB+1-1)/NPCOL+1)*NB
      lx = (((M-1)/NB+1-1)/NPROW+1)*NB

!$OMP PARALLEL PRIVATE(i,j,i_1,j_1,r)
!$OMP MASTER
!$    allocate( rr(1:omp_get_num_threads()) )
!$OMP END MASTER
!$OMP BARRIER

      r = zero
!$OMP DO COLLAPSE(2)
      do i_1=1,ly
        do j_1=1,lx
          i = (MOD(i_1-1,NB)+1)+(((i_1-1)/NB+1-1)*NPCOL+MYCOL-1)*NB
          j = (MOD(j_1-1,NB)+1)+(((j_1-1)/NB+1-1)*NPROW+MYROW-1)*NB
          if ( i <= n .and. j <= m ) then
            r = r + A(j_1,i_1)**2
          end if
        end do
      end do
!$OMP END DO
!$    rr(omp_get_thread_num()+1) = r
!$OMP BARRIER

!$OMP MASTER
!$    do i=2, omp_get_num_threads()
!$      r = r + rr(i)
!$    end do
!$    deallocate( rr )
      Fnorm2_local = r
!$OMP END MASTER
!$OMP END PARALLEL

!=      print*,MYROW,MYCOL, "Fnorm=",Fnorm2_local

      return
      end function Fnorm2_local

      real(4) function get_constant_eps_f32()
     &     result(r)
!DIR$ ATTRIBUTES FORCEINLINE :: get_constant_eps_f32
 
      integer(4), target :: const_eps_i4
      real(4), pointer   :: const_eps_r4
      data const_eps_i4 / z'34000000' /
 
      call c_f_pointer( c_loc(const_eps_i4), const_eps_r4 )
      r = const_eps_r4
 
      return
 
      end function get_constant_eps_f32

      end subroutine ev_test

