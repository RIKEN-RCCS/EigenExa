!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: bisect2_mod
!
! Purpose
! =======
!
! bisect2_mod manages the modules for computing eigenvalues of a
! symmetric penta-diagonal matrix by the bisection method.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module bisect2_mod

      use eigen_libs_mod
      use comm_mod
      use eigen_devel_mod
      use mpi
!$    use omp_lib

      implicit none
      private


      public  :: eigen_bisect2

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_bisect2
!
! Purpose
! =======
!
! eigen_bisect2 invokes the main body of the bisection solver,
! eigen_bisect2_body, to solve the eigenvalues for the symmetric
! penta-diagonal matrix.
!
!
! Arguments
! =========
!
! d       (input) real(8) array, dimension(n)
!         d contains the diagonal elements of the symmetric
!         penta-diagonal matrix.
!
! e       (input) real(8) array, dimension(n)
!         e contains the first off-diagonal elements of the symmetric
!         penta-diagonal matrix.

! f       (input/output) real(8) array, dimension(n)
!         f contains the second off-diagonal elements of the symmetric
!         penta-diagonal matrix.
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

      subroutine eigen_bisect2(d, e, f, w_out, n, mode)

      real(8), intent(in)    ::  d(1:n)
      real(8), intent(in)    ::  e(1:n)
      real(8), intent(inout) ::  f(1:n)
      real(8), intent(inout) ::  w_out(1:n)
      integer, intent(in)    ::  n
      integer, intent(in)    ::  mode

      real(8), allocatable   ::  wk(:)
      real(8), allocatable   ::  w(:)
      integer                ::  nv
      integer                ::  ierr

#if TIMER_PRINT
      real(8)                ::  timer
#endif

#if TIMER_PRINT
      timer = eigen_get_wtime()
#endif
      nv = (n-1)/TRD_nnod+1

      allocate(wk(1:nv*TRD_nnod+2),w(1:nv), stat=ierr)
      if (ierr /= 0) then
        call eigen_abort( "Memory allocation error [BISEC].", 1 )
      end if

      call eigen_bisect2_body(d, e, f, w_out, n, wk, w, mode)

      deallocate(wk)
      deallocate(w)

#if TIMER_PRINT
      timer = eigen_get_wtime() - timer
      if (TRD_inod == 1) print*,"BISECT2 ", timer
#endif
      return

      contains

!------------------------------------------------------------------
!------------------------------------------------------------------

      subroutine eigen_bisect2_body(d, e, f, w_out, n, wk, w, mode)

#define	DEBUG		0

      real(8), intent(in)    ::  d(1:n)
      real(8), intent(in)    ::  e(1:n)
      real(8), intent(inout) ::  f(1:n)
      real(8), intent(inout) ::  w_out(1:n)
      integer, intent(in)    ::  n
      real(8), intent(inout) ::  wk(*)
      real(8), intent(inout) ::  w(*)
      integer, intent(in)    ::  mode

      integer, parameter     ::  ITRMAX = 128

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
      v = ZERO
!$OMP PARALLEL DO
!$OMP+         PRIVATE(u, v, t)
!$OMP+         REDUCTION(MAX:x_upper)
!$OMP+         REDUCTION(MIN:x_lower)
      do i = 1, n
        u = abs(e(i))
        v = abs(f(i))
        t = u + v
        if (i == n) then
          u = ZERO
          v = ZERO
        else if (i == n-1) then
          u = abs(e(i+1))
          v = ZERO
        else
          u = abs(e(i+1))
          v = abs(f(i+2))
        end if
        t = t + u + v
        x_lower = min(x_lower, d(i)-t)
        x_upper = max(x_upper, d(i)+t)
      end do
!$OMP END PARALLEL DO

      t = ZERO
!$OMP PARALLEL DO
!$OMP+         REDUCTION(MAX:t)
      do i = 1, n
        t = max(t, abs(e(i))+abs(f(i)))
      end do                    ! j
!$OMP END PARALLEL DO
      epsilon = m_epsilon * t

      x  = (abs(x_lower)+abs(x_upper))*m_epsilon
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
          s = sturm2(w_out(i), d, e, f, n, epsilon)
          w(i-i_2+1) = dble(s)
        end do
!$OMP END PARALLEL DO
      end if
      call allgather_dbl(w(1), wk(2), nv, 1, TRD_COMM_WORLD)
      s = sturm2(lb_, d, e, f, n, epsilon)
      k = sturm2(ub_, d, e, f, n, epsilon)
      wk(1  ) = dble(s)
      wk(n+2) = dble(k)
#if DEBUG
      print*,wk(1:n+2)
      print*,"LB_",lb_
      print*,"UB_",ub_
#endif

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
            lb = min(w_out(k-1), w_out(j))
            ub = ub_
#if DEBUG
            print*,"T1",lb,ub
#endif
          else
            lb = min(w_out(k-1), w_out(j))
            ub = max(w_out(k), w_out(j))
#if DEBUG
            print*,"T2",lb,ub
#endif
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
            ub = max(w_out(k+1), w_out(j))
#if DEBUG
            print*,"T3",lb,ub
#endif
          else
            lb = min(w_out(k), w_out(j))
            ub = max(w_out(k+1), w_out(j))
#if DEBUG
            print*,"T4",lb,ub
#endif
          end if
        end if
#if DEBUG
        print*,j,"++++  :",lb,".",ub,",", w_out(j)
#endif

        if (mode == 1) then
          x = max(lb, ub - epsilon)
          if (x > lb) then
            s = sturm2(x, d, e, f, n, epsilon)
            if (s < sturm_val) then
              lb = x
#if DEBUG
              print*,j,"LB:",lb,".",s,sturm_val
#endif
            else
              ub = x
#if DEBUG
              print*,j,"UB:",ub,".",s,sturm_val
#endif
            end if
          end if
          x = min(ub, lb + epsilon)
          if (x < ub) then
            s = sturm2(x, d, e, f, n, epsilon)
            if (s < sturm_val) then
              lb = x
#if DEBUG
              print*,j,"LB:",lb,".",s,sturm_val
#endif
            else
              ub = x
#if DEBUG
              print*,j,"UB:",ub,".",s,sturm_val
#endif
            end if
          end if
        end if

        x = lb
        do itr=1,ITRMAX
          t = x
          x = (lb+ub)/2
          if (x == t) exit
          s = sturm2(x, d, e, f, n, epsilon)
          if (s < sturm_val) then
            lb = x
          else
            ub = x
          end if
#if DEBUG
          print*,j,":",lb,".",x,".",ub,s,sturm_val
#endif
        end do

#if DEBUG
!     theta = PAI*(2*(n-j)+1)/(2*n+1)
!     u = 5D-1/(1d0-cos(theta))
!     *-
        print*,j,":",w_out(j),"->",x
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

#if 1
      call lazy_qsort(n, w_out(1))
#endif

      return

      end subroutine eigen_bisect2_body

!------------------------------------------------------------------
!------------------------------------------------------------------
      integer function sturm2(x, d, e, f, n, epsilon)
     &     result(ret)

      real(8), intent(in)    ::  x
      real(8), intent(in)    ::  d(1:n)
      real(8), intent(in)    ::  e(1:n)
      real(8), intent(in)    ::  f(1:n)
      integer, intent(in)    ::  n
      real(8), intent(in)    ::  epsilon

!     Assumption
!     e(1) = f(1) = f(2) = 0

      ret = sturm2_LDLT(x, d, e, f, n,  epsilon)

      return

      end function sturm2

!------------------------------------------------------------------
!------------------------------------------------------------------

      integer function sturm2_LDLT(x, d, e, f, n, epsilon)
     &     result(ret)

      real(8), intent(in)    ::  x
      real(8), intent(in)    ::  d(1:n)
      real(8), intent(in)    ::  e(1:n)
      real(8), intent(in)    ::  f(1:n)
      integer, intent(in)    ::  n
      real(8), intent(in)    ::  epsilon

      integer                ::  i, j, k
      real(8)                ::  d0, e0, e1, f0, f1, g0, g1, t

#define	p(i,j)		p_ i _ j

      real(8)                ::  p_1_1
      real(8)                ::  p_2_1, p_2_2
      real(8)                ::  p_3_1, p_3_2, p_3_3
      real(8)                ::  p_4_1, p_4_2, p_4_3, p_4_4

      logical                ::  trail
!     
!     
!     
!     e(1) = f(1) = f(2) = 0
!     
!===============================================================
!     
!     We take a diagonal-neibor pivoting strategy
!     in order to stabilize the numerical process to count 
!     the number of negative signs for the diagonal matrix
!     resulted in the LDL decomposition
!     
      j = 0
!     
      p(1,1) = ONE
      p(2,1) = ZERO; p(2,2) = ONE
      p(3,1) = ZERO; p(3,2) = ZERO; p(3,3) = ONE
      p(4,1) = ZERO; p(4,2) = ZERO; p(4,3) = ZERO; p(4,4) = ONE

      trail = .false.
!OCL UXSIMD
      do i = 1, n

        if (trail) then

          e0 = p(2,1)
          if (abs(e0) < epsilon) e0 = epsilon
          f0 = p(3,1)
          g0 = p(4,1)
          f1 = p(3,2)
          g1 = p(4,2)

          p(1,1) = p(3,3) - 2*f0*f1/e0
          p(2,1) = p(3,2) - (g0*f1+g1*f0)/e0
          p(2,2) = p(3,3) - 2*g0*g1/e0
          p(3,1) = f(i-1)
          p(3,2) = e(i-1)
          p(3,3) = d(i-1)-x
          p(4,1) = ZERO
          p(4,2) = f(i)
          p(4,3) = e(i)
          p(4,4) = d(i)-x

          if (abs(e0) > epsilon) then
            j = j + 1
          else
            j = j + 2
          end if
          trail = .false.

        else

          if (abs(p(1,1)) < abs(p(2,2))) then
            t = p(1,1); p(1,1) = p(2,2); p(2,2) = t
            t = p(3,1); p(3,1) = p(3,2); p(3,2) = t
            p(4,1) = p(4,2); p(4,2) = ZERO ! since p(4,1) was set to 
!     0
          end if

          if (p(1,1) /= ZERO) then
!     
!     If p(1,1)==0, then
!     we can do block LDL decomposition for a 2x2 submatrix.
!     
            d0 = p(1,1)
            e0 = p(2,1)
            f0 = p(3,1)
            g0 = p(4,1)

            p(1,1) = p(2,2) - e0*e0/d0
            p(2,1) = p(3,2) - e0*f0/d0
            p(2,2) = p(3,3) - f0*f0/d0
            p(3,1) = p(4,2) - e0*g0/d0
            p(3,2) = p(4,3) - f0*g0/d0
            p(3,3) = p(4,4) - g0*g0/d0
            p(4,1) = ZERO
            p(4,2) = f(i)
            p(4,3) = e(i)
            p(4,4) = d(i)-x

            if (d0 < ZERO) j = j + 1
            trail = .false.

          else                  ! p(1,1) ==0 && p(2,2) == 0

            trail = .true.

          end if
        end if

      end do

      if (trail) then

        e0 = p(2,1)
!     
!     When p(2,1) is zero, we need another vipoting.
!     The present version, however, treat the p(2,1) element
!     as a small number (machine epsilon) in order to reduce
!     more complex pivoting cost.
!     
        if (abs(e0) < epsilon) e0 = epsilon
        f0 = p(3,1)
        g0 = p(4,1)
        f1 = p(3,2)
        g1 = p(4,2)

        p(1,1) = p(3,3) - 2*f0*f1/e0
        p(2,1) = p(3,2) - (g0*f1+g1*f0)/e0
        p(2,2) = p(3,3) - 2*g0*g1/e0
        p(3,1) = f(n)
        p(3,2) = e(n)
        p(3,3) = d(n)-x

        if (abs(e0) > epsilon) then
          j = j + 1
        else
          j = j + 2
        end if
        trail = .false.

      else

        if (abs(p(1,1)) < abs(p(2,2))) then
          t = p(1,1); p(1,1) = p(2,2); p(2,2) = t
          t = p(3,1); p(3,1) = p(3,2); p(3,2) = t
          p(4,1) = p(4,2); p(4,2) = ZERO ! since p(4,1) was set to 0
        end if

        if (p(1,1) /= ZERO) then

          d0 = p(1,1)
          e0 = p(2,1)
          f0 = p(3,1)
          g0 = p(4,1)

          p(1,1) = p(2,2) - e0*e0/d0
          p(2,1) = p(3,2) - e0*f0/d0
          p(2,2) = p(3,3) - f0*f0/d0
          p(3,1) = p(4,2) - e0*g0/d0
          p(3,2) = p(4,3) - f0*g0/d0
          p(3,3) = p(4,4) - g0*g0/d0

          if (d0 < ZERO) j = j + 1
          trail = .false.

        else                    ! p(1,1) ==0 && p(2,2) == 0

          trail = .true.

        end if
      end if

      if (trail) then

        e0 = p(2,1)
        if (abs(e0) < epsilon) e0 = epsilon
        f0 = p(3,1)
        g0 = p(4,1)
        f1 = p(3,2)
        g1 = p(4,2)

        p(1,1) = p(3,3) - 2*f0*f1/e0
        p(2,1) = p(3,2) - (g0*f1+g1*f0)/e0
        p(2,2) = p(3,3) - 2*g0*g1/e0

        if (abs(e0) > epsilon) then
          j = j + 1
        else
          j = j + 2
        end if
        trail = .false.

      else

        if (abs(p(1,1)) < abs(p(2,2))) then
          t = p(1,1); p(1,1) = p(2,2); p(2,2) = t
          t = p(3,1); p(3,1) = p(3,2); p(3,2) = t
        end if

        if (p(1,1) /= ZERO) then

          d0 = p(1,1)
          e0 = p(2,1)
          f0 = p(3,1)

          p(1,1) = p(2,2) - e0*e0/d0
          p(2,1) = p(3,2) - e0*f0/d0
          p(2,2) = p(3,3) - f0*f0/d0

          if (d0 < ZERO) j = j + 1
          trail = .false.

        else                    ! p(1,1) ==0 && p(2,2) == 0

          trail = .true.

        end if
      end if

      if (trail) then

        e0 = p(2,1)
        if (abs(e0) < epsilon) e0 = epsilon
        f0 = p(3,1)
        f1 = p(3,2)

        p(1,1) = p(3,3) - 2*f0*f1/e0

        if (abs(e0) > epsilon) then
          j = j + 1
        else
          j = j + 2
        end if
        trail = .false.

      else

        if (abs(p(1,1)) < abs(p(2,2))) then
          t = p(1,1); p(1,1) = p(2,2); p(2,2) = t
        end if

        if (p(1,1) /= ZERO) then

          d0 = p(1,1)
          e0 = p(2,1)

          p(1,1) = p(2,2) - e0*e0/d0

          if (d0 < ZERO) j = j + 1
          trail = .false.

        else                    ! p(1,1) ==0 && p(2,2) == 0

          trail = .true.

        end if
      end if

      if (trail) then

        e0 = p(2,1)
        if (abs(e0) < epsilon) e0 = epsilon

        if (abs(e0) > epsilon) then
          j = j + 1
        else
          j = j + 2
        end if

      else

        d0 = p(1,1)

        if (d0 < ZERO) j = j + 1

      end if

#undef  p

      ret =  j

      return

      end function sturm2_LDLT

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

      end subroutine eigen_bisect2

      end module bisect2_mod
