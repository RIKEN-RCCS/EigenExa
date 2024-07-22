!--------*---------*---------*---------*---------*---------*---------*-*
!     
! File: mat_set.f
!
! Purpose
! =======
!
! < purpose of this module ... >
!
!
! Copyright(C) 2012-2021 RIKEN.
! Copyright(C) 2011-2012 Toshiyuki Imamura
!                        Graduate School of Informatics and Engineering,
!                        The University of Electro-Communications.
! Copyright (C) 2011- 2015 Japan Atomic Energy Agency.
! 
! Redistribution  and  use  in  source and binary forms, with or without
! modification,  are  permitted  provided  that the following conditions
! are met:
! 
! * Redistributions  of  source  code  must  retain  the above copyright
!   notice,  this  list  of  conditions  and  the  following  disclaimer.
! * Redistributions  in  binary  form must reproduce the above copyright
!   notice,  this list of conditions and the following disclaimer in the
!   documentation  and/or other materials provided with the distribution.
! 
! THIS  SOFTWARE  IS  PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! ``AS IS''  AND  ANY  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A  PARTICULAR  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL,  EXEMPLARY,  OR  CONSEQUENTIAL  DAMAGES  (INCLUDING,  BUT NOT
! LIMITED  TO,  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA,  OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY  OF  LIABILITY,  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF  THIS  SOFTWARE,  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine mat_set(n, a, nm, mtype)

      use eigen_libs_mod
      use eigen_blacs_mod
      use mpi
!$    use omp_lib

      implicit none

      integer, intent(in)    :: n, nm, mtype
      real(8), intent(out)   :: a(1:nm,*)

!     Parameters BLACS array descritor(the position of entry tags), etc
      integer, parameter     :: BLOCK_CYCLIC_2D = 1
      integer, parameter     :: DLEN_  = 9
      integer, parameter     :: DTYPE_ = 1
      integer, parameter     :: CTXT_  = 2
      integer, parameter     :: M_     = 3
      integer, parameter     :: N_     = 4
      integer, parameter     :: MB_    = 5
      integer, parameter     :: NB_    = 6
      integer, parameter     :: RSRC_  = 7
      integer, parameter     :: CSRC_  = 8
      integer, parameter     :: LLD_   = 9

      real(8), parameter     :: ZERO= 0.0D0
      real(8), parameter     :: ONE = 1.0D0

      integer                :: DESCA(DLEN_)
      integer                :: NPROW, NPCOL

      real(8), pointer       :: as(:,:), w(:)
      integer, pointer       :: irow(:), icol(:)
      real(8)                :: t, s, hi_, hj_

      integer                :: COMM, x_COMM, y_COMM
      integer                :: nnod, x_nnod, y_nnod
      integer                :: inod, x_inod, y_inod
      integer                :: iloop_sta, iloop_end
      integer                :: jloop_sta, jloop_end
      integer                :: i, i_1
      integer                :: j, j_1
      integer                :: k, k0
      integer                :: nx, info, ICTXT, ierr
      integer, pointer       :: iseed(:)

      external               :: DESCINIT, PDTRAN

      real(8)                :: PAI, EPS, EPS2, EPS4, theta

      integer :: n1, n2, ne, ierror, ns(1:3)

      character(len=1024) :: file_name
      character(len=1024) :: buff


      if (mtype < -2 .or. mtype > 10) then
        print*,"Illeagal Matrix Type is specified."
        ierror = -1
        call MPI_Abort(COMM, MPI_ERR_OTHER, ierror)
      end if

      call eigen_get_comm (COMM, x_COMM, y_COMM)
      call eigen_get_procs(nnod, x_nnod, y_nnod)
      call eigen_get_id   (inod, x_inod, y_inod)

      jloop_sta = eigen_loop_start(1, 'X')
      jloop_end = eigen_loop_end  (n, 'X')
      iloop_sta = eigen_loop_start(1, 'Y')
      iloop_end = eigen_loop_end  (n, 'Y')

      PAI  = get_constant_pai()
      EPS  = get_constant_eps()
      EPS2 = sqrt(EPS)
      EPS4 = sqrt(EPS2)

      if (mtype == 0) then

        if (iloop_sta <= iloop_end) then
!$OMP PARALLEL DO
!$OMP+         PRIVATE(i,j,i_1,j_1)
          do i_1 = iloop_sta, iloop_end
            i = eigen_translate_l2g(i_1, 'Y')
            do j_1 = jloop_sta, jloop_end
              j = eigen_translate_l2g(j_1, 'X')
              a(j_1, i_1) = DBLE(min(i,j))
            end do
          end do
!$OMP END PARALLEL DO
        end if

      end if

      if (mtype == 1) then

        if (iloop_sta <= iloop_end) then
!$OMP PARALLEL DO
!$OMP+         PRIVATE(i,j,i_1,j_1)
          do i_1 = iloop_sta, iloop_end
            i = eigen_translate_l2g(i_1, 'Y')
            do j_1 = jloop_sta, jloop_end
              j = eigen_translate_l2g(j_1, 'X')
              if ( i == j ) then
                a(j_1, i_1) = -7.2D+00
              else
                a(j_1, i_1) = -3.0D+00/(i-j)**2
              end if
            end do
          end do
!$OMP END PARALLEL DO
        end if

      end if

      if (mtype == 2) then

        NPROW = x_nnod
        NPCOL = y_nnod

        ICTXT = eigen_get_blacs_context( )

        call DESCINIT(DESCA, n, n, 1, 1, 0, 0, ICTXT, nm, INFO)

        nx = (n-1)/NPCOL+1
        allocate(as(1:nm,nx))

        call random_seed(size = i)
        allocate(iseed(i))
        iseed(1:i) = inod
        call random_seed(put = iseed)
        deallocate(iseed)

        do i_1 = iloop_sta, iloop_end
          do j_1 = jloop_sta, jloop_end
            call random_number(t)
            as(j_1, i_1) = t
            a (j_1, i_1) = t
          end do
        end do

        call PDTRAN(n, n,
     &       ONE, as, 1, 1, DESCA, ONE, a, 1, 1, DESCA)

        deallocate(as)

      end if

      if (mtype == 3) then

         if (iloop_sta <= iloop_end) then
!$OMP PARALLEL DO
!$OMP+         PRIVATE(i,j,i_1,j_1)
           do i_1 = iloop_sta, iloop_end
             i = eigen_translate_l2g(i_1, 'Y')
             do j_1 = jloop_sta, jloop_end
               j = eigen_translate_l2g(j_1, 'X')
               a(j_1, i_1) = DBLE(n+1-max(i,j))
             end do
           end do
!$OMP END PARALLEL DO
         end if

      end if

      if (4 <= mtype .and. mtype <= 10) then

        allocate (w(1:n))
        if (mtype /= 10 .or. inod == 1) then
          call w_set(n, w, mtype)
        end if
        if (mtype == 10) then
          call MPI_Bcast(w, n, MPI_DOUBLE_PRECISION, 0, COMM, ierr)
        end if
        call helmert_trans(n, a, nm, w, 
     &       iloop_sta, iloop_end, jloop_sta, jloop_end)
        deallocate (w)

      end if

      if (mtype < 0) then

        if (mtype == -1) then
          file_name = "A.mtx"
        else
          file_name = "B.mtx"
        end if

        if (iloop_sta <= iloop_end) then
!$OMP PARALLEL DO
!$OMP+         PRIVATE(i_1,j_1)
          do i_1 = iloop_sta, iloop_end
            do j_1 = jloop_sta, jloop_end
              a(j_1, i_1) = ZERO
            end do
          end do
!$OMP END PARALLEL DO
        end if

        if (inod == 1) then
          open(unit=11, file=file_name, status='old', recl=1024,
     &         iostat=ierror)
          if (ierror /= 0) then
            print*,"Can not open input file [",trim(file_name),"]"
            call MPI_Abort( MPI_COMM_WORLD, 1, ierror )
          end if
          do
            read(UNIT=11, FMT='(a1024)',
     &           ERR=10000, END=10000 ) buff
            if ( buff(1:1) /= '%') exit
          end do
          read(buff, *) n1, n2, ne
        endif

        ns(1) = n1
        ns(2) = n2
        ns(3) = ne
        call MPI_Bcast(ns, 3, MPI_INTEGER, 0, COMM, ierr)
        n1 = ns(1)
        n2 = ns(2)
        ne = ns(3)

        if (n1 /= n .or. n2 /= n) then
          print*,"Matrix size inconsistency has been found."
          print*,"Illeagal MTX data is specified."
          ierror = -2
          call  MPI_Abort(COMM, MPI_ERR_OTHER, ierror)
        end if

        allocate(w(4096), irow(4096), icol(4096))

        do k0 = 1, ne, 4096

          if (inod == 1) then
            do k = 1, min(4096,ne-k0+1)
!     read( UNIT=11, FMT='(a1024)',
!     &                      ERR=10000, END=10000 ) buff
!     read( buff, * ) n1, n2, t
              read(unit=11, fmt=*,
     &             err=10000, end=10000) n1, n2, t
              irow(k) = n1
              icol(k) = n2
              w   (k) = t
            end do
          end if

          call MPI_Bcast(w,    4096, MPI_DOUBLE_PRECISION, 0, COMM,
     &         ierr)
          call MPI_Bcast(irow, 4096, MPI_INTEGER,          0, COMM,
     &         ierr)
          call MPI_Bcast(icol, 4096, MPI_INTEGER,          0, COMM,
     &         ierr)

!$OMP PARALLEL DO
!$OMP+         PRIVATE(i,j,k,n1,n2,t,i_1,j_1)
          do k = 1, min(4096,ne-k0+1)

            n1 = irow(k)
            n2 = icol(k)
            t  = w   (k)

            i = eigen_owner_node(n1, 'Y')
            j = eigen_owner_node(n2, 'X')
            if (i == y_inod .and. j == x_inod) then
              i_1 = eigen_translate_g2l(n1, 'Y')
              j_1 = eigen_translate_g2l(n2, 'X')
              a(j_1, i_1) = t
            end if

            if (n1 /= n2) then
              i = eigen_owner_node(n2, 'Y')
              j = eigen_owner_node(n1, 'X')
              if (i == y_inod .and. j == x_inod) then
                i_1 = eigen_translate_g2l(n2, 'Y')
                j_1 = eigen_translate_g2l(n1, 'X')
                a(j_1, i_1) = t
              end if
            end if

          end do
!$OMP END PARALLEL DO

        end do

        deallocate(w, irow, icol)

10000   continue
        if (inod == 1) then
          close(11)
        end if

      end if

      return

      contains

      subroutine helmert_trans(n, a, lda, w,
     &     iloop_sta, iloop_end, jloop_sta, jloop_end)

      integer, intent(in)    :: n, lda
      integer, intent(in)    :: iloop_sta, iloop_end
      integer, intent(in)    :: jloop_sta, jloop_end
      real(8), intent(out)   :: a(lda, *)
      real(8), intent(in)    :: w(*)

      integer                :: i_1, j_1, i, j
      real(8)                :: hi_, hj_, t, s, scale
      real(8), pointer       :: w_(:), ht(:)
      real(8), pointer       :: hi(:), hj(:)
      integer, pointer       :: iseed(:)


      allocate (w_(1:n), ht(1:n))

      scale = ZERO
!$OMP PARALLEL DO
!$OMP+         REDUCTION(max: scale)
      do k=1,n
        scale = max(scale, ABS(w(k)))
      end do
!$OMP END PARALLEL DO
      if (scale < ONE) then
        scale = ONE
      end if
!$OMP PARALLEL DO
      do k=1,n
        w_(k) = w(k) / scale
      end do
!$OMP END PARALLEL DO

      call random_seed(size = i)
      allocate(iseed(i))
      iseed(1:i) = 0
      call random_seed(put = iseed)
      deallocate(iseed)
      do k=1,n
        call random_number(t)
        i = int(t*n)+1
        call random_number(t)
        j = int(t*n)+1
        t = w_(i); w_(i) = w_(j); w_(j) = t
      end do

!$OMP PARALLEL
!$OMP+         PRIVATE(i,j,k, i_1,j_1, s,t, hi,hj, hi_,hj_)

      allocate(hi(1:n), hj(1:n))

!$OMP DO
      do i=1,n
        ht(i) = sqrt(dble(i))
      end do
!$OMP END DO

!$OMP DO
      do i_1 = iloop_sta, iloop_end
        i = eigen_translate_l2g(i_1, 'Y')

        if ( i==1 ) then
          hi_ = ht(n)           ! SQRT(DBLE(n))
          hi(1:n)   = ONE / hi_
        else
          s   = dble(i-1)
!         hi_ = s * sqrt(ONE+ONE/s)
          hi_ = ht(i-1) * ht(i) ! SQRT(s) * SQRT(s+ONE)
          hi(1:i-1) =  ONE / hi_
          hi(i)     = -s / hi_
          hi(i+1:n) =  ZERO
        end if

        do j_1 = jloop_sta, jloop_end
          j = eigen_translate_l2g(j_1, 'X')

          if (j==1) then
            hj_ = ht(n)         ! SQRT(DBLE(n))
            hj(1:n)   = ONE / hj_
          else
            s   = dble(j-1)
!           hj_ = s * sqrt(ONE+ONE/s)
            hj_ = ht(j-1)*ht(j) ! SQRT(s) * SQRT(s+ONE)
            hj(1:j-1) =  ONE / hj_
            hj(j)     = -s / hj_
            hj(j+1:n) =  ZERO
          end if

          t = ZERO
!$OMP SIMD
          do k=1,n
            hj(k) = hj(k) * hi(k)
            t = t + w_(k) * hj(k)
          end do
          a(j_1, i_1) = t       ! / (hi_ * hj_)

        end do
      end do
!$OMP END DO

!$OMP DO
      do i_1 = iloop_sta, iloop_end
        do j_1 = jloop_sta, jloop_end
          a(j_1, i_1) = a(j_1, i_1) * scale
        end do
      end do
!$OMP END DO

      deallocate(hi, hj)

!$OMP END PARALLEL

      deallocate(w_, ht)

      return

      end subroutine  helmert_trans

      end subroutine  mat_set

!----------------------------------------------------------------------
!----------------------------------------------------------------------

      subroutine mat_dim_get(mtype, n)

      use eigen_libs_mod
      use mpi
!$    use omp_lib

      implicit none

      integer, intent(in)    :: mtype
      integer, intent(out)   :: n

      integer                :: n1, n2, ne, ierror, ierr, ns(3)
      integer                :: COMM, x_COMM, y_COMM
      integer                :: nnod, x_nnod, y_nnod
      integer                :: inod, x_inod, y_inod

      character(len=1024)    :: file_name
      character(len=1024)    :: buff


      if (mtype == -1 .or. mtype == -2) then

        if (mtype == -1) then
          file_name = "A.mtx"
        else
          file_name = "B.mtx"
        end if

        call eigen_get_comm (COMM, x_COMM, y_COMM)
        call eigen_get_procs(nnod, x_nnod, y_nnod)
        call eigen_get_id   (inod, x_inod, y_inod)

        if (inod == 1) then
          open(unit=11, file=file_name, status='old', recl=1024,
     &         iostat=ierror)
          if (ierror /= 0) then
            print*,"Can not open input file [",trim(file_name),"]"
            call MPI_Abort(MPI_COMM_WORLD, 1, ierror)
          end if
          do
            read(unit=11, fmt='(a1024)',
     &           err=10000, end=10000) buff
            if (buff(1:1) /= '%') exit
          end do
          read(buff, *) n1, n2, ne
        end if

        ns(1) = n1
        ns(2) = n2
        ns(3) = ne
        call MPI_Bcast(ns, 3, MPI_INTEGER, 0, COMM, ierr)
        n1 = ns(1) 
        n2 = ns(2)
        ne = ns(3) 

        if (n1 /= n2) then
          print*,"Matrix size inconsistency has been found."
          print*,"Non-square MTX data is specified."
          n = -1
        else
          n = n1
        end if

      end if

10000 continue
      if (inod == 1) then
        close(11)
      end if

      return

      end subroutine mat_dim_get

!----------------------------------------------------------------------
!----------------------------------------------------------------------

      subroutine print_mat_name_list()

      implicit none

      integer :: mtype


      do mtype = 0, 10
        call  print_mat_name(mtype)
      end do
      call  print_mat_name(-1)
      call  print_mat_name(-2)

      return

      end subroutine print_mat_name_list

!----------------------------------------------------------------------
!----------------------------------------------------------------------

      subroutine print_mat_name(mtype)

      implicit none

      integer, intent(in)    :: mtype

      character*(100)        :: message


      message = " (unkonwn or invalid)"
      select case (mtype)
      case (0)
        message = " (Frank matrix)"
      case (1) 
        message = " (Toeplitz matrix)"
      case (2)
        message = " (Random matrix)"
      case (3)
        message = " (Frank matrix 2)"
      case (4)
        message = " (W: 0, 1, ..., n-1)"
      case (5)
        message = " (W: sin(PAI*5*i/(n-1)+EPS^1/4)^3)"
      case (6)
        message = " (W: MOD(i,5)+MOD(i,2))"
      case (7)
        message = " (W: same as Frank matrix)"
      case (8)
        message = " (W: Uniform Distribution, [0,1))"
      case (9)
        message = " (W: Gauss Distribution, m=0,s=1)"
      case (10)
        message = " (W: Read from the data file 'W.dat')"
      case (-1)
        message = " (Read from the data file 'A.mtx')"
      case (-2)
        message = " (Read from the data file 'B.mtx')"
      end select

      print*, "Matrix type = ", mtype, trim(message)

      return

      end subroutine print_mat_name

!----------------------------------------------------------------------
!----------------------------------------------------------------------

      subroutine w_set(n, w, mtype)

      use eigen_libs_mod
      use mpi

      implicit none

      integer, intent(in)    :: n, mtype
      real(8), intent(out)   :: w(1:n)

      real(8), parameter     :: ZERO = 0.0D0
      real(8), parameter     :: ONE  = 1.0D0

      real(8)                :: PAI, EPS, EPS2, EPS4
      real(8)                :: ax, bx, x, y, z, theta, s, t
      real(8), pointer       :: ww(:)
      integer, pointer       :: iseed(:)

      integer                :: i, j, k, kk, l, ierror


      if (mtype == 1 .or. mtype == 2) then

        return

      end if

      PAI  = get_constant_pai()
      EPS  = get_constant_eps()
      EPS2 = SQRT(EPS)
      EPS4 = SQRT(EPS2)

      if (mtype == 0 .or. mtype == 3 .or. mtype == 7) then

!$OMP PARALLEL DO
!$OMP+         PRIVATE(j,theta)
        do i = 1, n
          j = n-i
          theta = PAI*(2*j+1)/(2*n+1)
          w(i) = 5D-1/(ONE-cos(theta))
        end do
!$OMP END PARALLEL DO

      end if

      if (mtype == 4) then

!$OMP PARALLEL DO
        do i = 1, n
          w(i) = dble(i-1)
        end do
!$OMP END PARALLEL DO

      end if

      if (mtype == 5) then

!$OMP PARALLEL DO
!$OMP+         PRIVATE(theta)
        do i = 1, n
          theta = PAI*5*i/(n-1)+EPS4
          w(i) = sin(theta)**3
        end do
!$OMP END PARALLEL DO

      end if

      if (mtype == 6) then

!$OMP PARALLEL DO
        do i = 1, n
          w(i) = dble(mod(i, 5)+mod(i, 2))
        end do
!$OMP END PARALLEL DO

      end if

      if (mtype == 8) then

        call random_seed(size = i)
        allocate(iseed(i))
        iseed(1:i) = 0
        call random_seed(put = iseed)
        deallocate(iseed)

        do i = 1, n
          call random_number(t)
          w(i) = t
        end do

      end if

      if ( mtype == 9 ) then

        call random_seed(size = i)
        allocate(iseed(i))
        iseed(1:i) = 0
        call random_seed(put = iseed)
        deallocate(iseed)

        do i = 1, n
          call random_number(t)
          call random_number(s)
          w(i) = sqrt(-2*log(s)) * sin(2*PAI*s)
        end do

      end if

      if (mtype == 10) then

        open(unit=11, file='W.dat', status='old', recl=1024,
     &       iostat=ierror)
        if (ierror /= 0) then
          print*,"Can not open input file [",'W.dat',"]"
          call MPI_Abort(MPI_COMM_WORLD, 1, ierror)
        end if
        read(11, *) w(1:n)
        close(11)

      end if

      return

      end subroutine w_set
