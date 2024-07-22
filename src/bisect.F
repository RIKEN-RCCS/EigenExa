!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: bisect_mod
!
! Purpose
! =======
!
! bisect_mod manages the modules for computing eigenvalues of a
! symmetric tridiagonal matrix by the bisection method.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module bisect_mod

      use eigen_libs_mod
      use comm_mod
      use eigen_devel_mod
      use mpi
!$    use omp_lib

      implicit none
      private


      public  :: eigen_bisect

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_bisect
!
! Purpose
! =======
!
! eigen_bisect invokes the main body of the bisection solver,
! eigen_bisect_body, to solve the eigenvalues for the symmetric
! tridiagonal matrix.
!
!
! Arguments
! =========
!
! d       (input) real(8) array, dimension(n)
!         d contains the diagonal elements of the symmetric tridiagonal
!         matrix.
!
! e       (input) real(8) array, dimension(n)
!         e contains the off-diagonal elements of the symmetric
!         tridiagonal matrix.
!
! w_out   (input/output) real(8) array, dimension(n)
!         On entry, when mode is nozero, w_out contains
!         appropriximation of the eigenvalues which are treated as the
!         initial intervals of the bisection method.
!         On exit, w_out contains eigenvalues of the input matrix.
!
! n       (input) integer
!         The dimension of the symmetric tridiagonal matrix. N >= 0.
!
! mode    (output) integer
!         = 0: without approximated interval.
!         !=0: refine w_out passed as the initial intervals 
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_bisect(d, e, w_out, n, mode)

      real(8), intent(in)    ::  d(1:n)
      real(8), intent(in)    ::  e(1:n)
      real(8), intent(inout) ::  w_out(1:n)
      integer, intent(in)    ::  n
      integer, intent(in)    ::  mode

      real(8), allocatable   ::  wk(:)
      real(8), allocatable   ::  w(:)
      integer                ::  nv
      integer                ::  ierr

#if TIMER_PRINT>1
      real(8)                ::  timer
#endif

#if TIMER_PRINT>1
      timer = eigen_get_wtime()
#endif

      nv = (n-1)/TRD_nnod+1

      allocate(wk(1:nv*TRD_nnod+2),w(1:nv), stat=ierr)
      if (ierr /= 0) then
        call eigen_abort( "Memory allocation error [BISEC].", 1 )
      end if

      call eigen_bisect_body(d, e, w_out, n, wk, w, mode)

      deallocate(wk)
      deallocate(w)

#if TIMER_PRINT>1
      timer = eigen_get_wtime() - timer
      if (TRD_inod == 1) print*, "BISECT ", timer
#endif
      return

      contains

!------------------------------------------------------------------
!------------------------------------------------------------------
      subroutine  eigen_bisect_body(d, e, w_out, n, wk, w, mode)

#define	DEBUG		0

      real(8), intent(in)    ::  d(1:n)
      real(8), intent(in)    ::  e(1:n)
      real(8), intent(inout) ::  w_out(1:n)
      integer, intent(in)    ::  n
      real(8), intent(inout) ::  wk(*)
      real(8), intent(inout) ::  w(*)
      integer, intent(in)    ::  mode

      integer, parameter     ::  ITR_MAX = 128

      real(8)                ::  u, v, x, t
      real(8)                ::  x_lower, x_upper
      real(8)                ::  m_epsilon, epsilon
#if DEBUG
      real(8)                ::  PAI, theta
#endif
      real(8)                ::  lb, ub, lb_, ub_
      integer                ::  i, j, k, s
      integer                ::  i_1,i_2,i_3, nv
      integer                ::  sturm_val
      integer                ::  itr


      if (n <= 0) return
      if (n == 1) then
        w_out(1) = d(1)
        return
      end if
     
      m_epsilon = get_constant_eps()
#if DEBUG
      PAI = get_constant_pai()
#endif

      x_lower = d(1)
      x_upper = d(1)
      u = ZERO
!$OMP PARALLEL DO
!$OMP+         PRIVATE(u, v, t)
!$OMP+         REDUCTION(MAX:x_upper)
!$OMP+         REDUCTION(MIN:x_lower)
      do i = 1, n
        u = abs(e(i))
        t = u
        if (i == n) then
          u = ZERO
        else
          u = abs(e(i+1))
        end if
        t = t + u
        x_lower = min(x_lower, d(i)-t)
        x_upper = max(x_upper, d(i)+t)
      end do
!$OMP END PARALLEL DO
!     
      t = ZERO
!$OMP PARALLEL DO
!$OMP+         REDUCTION(MAX:t)
      do i = 1, n
        t = max(t, abs(e(i)))
      end do                    ! j
!$OMP END PARALLEL DO
      epsilon = m_epsilon * t

      x   = (abs(x_lower)+abs(x_upper))*m_epsilon
      lb_ = (x_lower - x) - epsilon
      ub_ = (x_upper + x) + epsilon

      nv = (n-1)/TRD_nnod+1
      i_2 = nv*(TRD_inod-1)+1
      i_3 = min(n, i_2+nv-1)

      if (mode == 0) then
!$OMP PARALLEL DO
        do i=1,n
          w_out(i) = (lb_*(n+1-i)+ub_*(i+1-1))/(n+1)
        end do
!$OMP END PARALLEL DO
      end if

      if (i_2 <= i_3) then
!$OMP PARALLEL DO
!$OMP+         PRIVATE(s)
        do i=i_2,i_3
          s = sturm(w_out(i), d, e, n, m_epsilon)
          w(i-i_2+1) = dble(s)
        end do
!$OMP END PARALLEL DO
      end if
      call allgather_dbl(w(1), wk(2), nv, 1, TRD_COMM_WORLD)
      s = sturm(lb_, d, e, n, m_epsilon)
      k = sturm(ub_, d, e, n, m_epsilon)
      wk(1  ) = dble(s)
      wk(n+2) = dble(k)

      if ( i_2 <= i_3 ) then
#if !DEBUG
!$OMP PARALLEL DO
!$OMP+         PRIVATE(i, j, k, s, sturm_val, t, x, lb, ub)
#endif
      do i_1 = i_2, i_3
        j = i_1

        sturm_val = j
!     
!     Find i such that "wk(i) < sturm_val <= wk(i+1)" holds
        s = int(wk((j+1)))
        if (s < sturm_val) then
          k = n+1
          do i=j+1,n
            if (int(wk((i+1))) >= sturm_val) then
              k = i
              exit
            end if
          end do
          if (k == n+1) then
            lb = w_out(k-1)
            ub = ub_
          else
            lb = w_out(k-1)
            ub = w_out(k  )
          end if
        else
!     sturm_val <= s
          k = 0
          do i=j-1,1,-1
            if (int(wk(i+1)) < sturm_val) then
              k = i
              exit
            end if
          end do
!     w_out(k) < val <= w_out(k+1)
          if (k == 0) then
            lb = lb_
            ub = w_out(1)
          else
            lb = w_out(k)
            ub = w_out(k+1)
          end if
        end if

        if (mode == 1) then
          x = max(lb, ub - epsilon)
          if (x > lb) then
            s = sturm(x, d, e, n, m_epsilon)
            if (s < sturm_val) then
              lb = x
            else
              ub = x
            end if
          end if
          x = min(ub, lb + epsilon)
          if (x < ub) then
            s = sturm(x, d, e, n, m_epsilon)
            if (s < sturm_val) then
              lb = x
            else
              ub = x
            end if
          end if
        end if

        x = lb
        do itr=1,ITR_MAX
          t = x
          x = (lb+ub)/2
          if (x == t) exit
          s = sturm(x, d, e, n, m_epsilon)
          if (s < sturm_val) then
            lb = x
          else
            ub = x
          end if
        end do

#if DEBUG
        theta = PAI*(2*(n-j)+1)/(2*n+1)
        u = 5D-1/(1d0-cos(theta))

        print*,j,":",w_out(j),"->",x,"/",u,abs(x-u)/abs(x)
#endif

        w(j-i_2+1) = x

      end do                    ! i_1
#if !DEBUG
!$OMP END PARALLEL DO
#endif
      end if

      call allgather_dbl(w(1), wk(1), nv, 2, TRD_COMM_WORLD)

!$OMP PARALLEL DO
      do j=1,n
        w_out(j) = wk(j)
      end do
!$OMP END PARALLEL DO
!     
#if 1
      call lazy_qsort(n, w_out(1))
#endif
!     
      return

      end subroutine eigen_bisect_body

!------------------------------------------------------------------
!------------------------------------------------------------------
      integer function sturm(x, d, e, n, epsilon)
     &     result(ret)

      integer, intent(in)    ::  n
      real(8), intent(in)    ::  x
      real(8), intent(in)    ::  d(1:n)
      real(8), intent(in)    ::  e(1:n)
      real(8), intent(in)    ::  epsilon

      integer                ::  i, j
      real(8)                ::  t

      real(8)                ::  H2,H3,H4, d0,e0

      H2 = ONE; H3 = ONE

      j = 0
!OCL UXSIMD
      do i = 1, n
        d0 = (d(i)-x)
        e0 = e(i-1+1)
        H4 = (d0*H3 - (e0**2)*H2)
        if (abs(H4) > epsilon * (ABS(d0)+ABS(e0))) then
          t = ONE / H4
          H2 = H3 * t
          H3 = ONE
          if (H2 < 0) j = j + 1
        else
          H2 = H3
          H3 = epsilon
        end if
      end do
      ret = j

      return

      end function sturm

!------------------------------------------------------------------
!------------------------------------------------------------------
#if 1
      recursive subroutine lazy_qsort(n, w)

      integer, intent(in)    :: n
      real(8), intent(inout) :: w(n)

      real(8) :: r, t
      integer :: i, m, mid


      if (n < 2) return

      mid = (1+n)/2

      m = mid; r = w(m)
      do i=1,mid-1
        if (w(i) > r) then
          w(m) = w(i); w(i) = r; m = i
        end if
      end do
      do i=mid+1,n
        if (w(i) < r) then
          w(m) = w(i); w(i) = r; m = i
        end if
      end do

      call lazy_qsort(m-1, w(1)  )
      call lazy_qsort(n-m, w(m+1))

      return 

      end subroutine lazy_qsort
#endif
!------------------------------------------------------------------
!------------------------------------------------------------------
      end subroutine eigen_bisect

      end module bisect_mod
