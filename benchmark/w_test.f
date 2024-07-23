!--------*---------*---------*---------*---------*---------*---------*-*
!     
! File: w_test.f
!
! Purpose
! =======
!
! < purpose of this module ... >
!
!
! Copyright(C) 2012-2024 RIKEN.
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

      subroutine w_test(n, w, mtype)

      use eigen_libs_mod
      use mpi
      use, intrinsic :: ieee_arithmetic

      implicit none

      integer, intent(in)    :: n, mtype
      real(8)                :: w(1:n)

      real(8)                :: PAI, EPS, EPS2, EPS4
      real(8)                :: ax, bx, x, y, z, theta, s, t
      real(8)                :: awmin, awmax, kappa
      real(8), pointer       :: ww(:)
      integer, pointer       :: iseed(:)

      real(8), parameter     :: ONE  = 1D0
      real(8), parameter     :: ZERO = 0D0

      integer                :: i, j, k, kk, l, ierr


      if (mtype < -2 .or. mtype > 10) then
        print*,"Invalid matrix type ",mtype
        return
      end if

      if (mtype == -1 .or. mtype == -2 .or.
     &    mtype ==  1 .or. mtype ==  2 ) then
        print*,"-----------------------------------------------"
        print*,"*** Eigenvalue Error Test *** : SKIP"
        print*,"Because the eigenvalues are not given by user"
        print*,"or it is hard to solve them analytically."
        print*,"-----------------------------------------------"
        return
      end if

      PAI  = get_constant_pai()
      EPS  = get_constant_eps()
      EPS2 = sqrt(EPS)
      EPS4 = sqrt(EPS2)

      do i = 1, n
        if ( ieee_is_nan(w(i)) ) then
           print*,"NAN detected :: EV(",i,") = ",w(i)
           print*,"NAN detected."
           call MPI_Abort( MPI_COMM_WORLD, 1, ierr )
        end if
      enddo

      allocate (ww(1:n))
      call w_set(n, ww, mtype)

      awmin = abs(ww(n))
      awmax = abs(ww(n))
!
! Bucket sort
!
      do i = 1, n-1
        awmin = min(abs(ww(i)), awmin)
        awmax = max(abs(ww(i)), awmax)
        do j = i+1, n
          if (ww(i) > ww(j)) then
            x = ww(i); ww(i) = ww(j); ww(j) = x
          end if
        end do
      end do

      ax=ZERO; k  = 1
      bx=ZERO; kk = 1
      do i = 1, n
        x = ww(i)
        y = abs(w(i)-x)
        if (x == ZERO) then
          z = ZERO ! we do not check the relative error in this case
        else
          z = y / abs(x)
        end if
        if (ax < z) then
          ax = max(ax, z)
          k = i
        end if
        if (bx < y) then
          bx = max(bx, y)
          kk = i
        end if
      end do


      print*,"-----------------------------------------------"

      print*, "cond(A)=|w_max|/|w_min|=", awmax,"/",awmin
      if (awmin > ZERO) then
        kappa = awmax/awmin
        print*, "       =", kappa
      else
        kappa = ZERO
      end if

      print*, "max|w(i)-w(i).true|/|w.true|=", ax, w(k)
      if (ax < EPS2) then
        print*, "*** Eigenvalue Relative Error *** : PASSED"
      else if (ax < EPS4) then
        print*, "*** Eigenvalue Relative Error *** : CAUTION"
      else
        print*, "*** Eigenvalue Relative Error *** : FAILED"
        if (abs(w(k)) < EPS) then
          print*, " |w| is too small, so it is not severe."
        end if
      end if

      print*, "max|w(i)-w(i).true|         =", bx, w(kk)
      if (bx < EPS2) then
        print*, "*** Eigenvalue Absolute Error *** : PASSED"
      else if (bx < EPS4) then
        print*, "*** Eigenvalue Absolute Error *** : CAUTION"
      else
        print*, "*** Eigenvalue Absolute Error *** : FAILED"
        if (ax < EPS4) then
          if (2* kappa * EPS2 >= ONE) then
            print*, " Do not mind it. Condition number is too large."
          else
            print*, " Do not mind it. Relative error is small enough"
     &           //  "."
          end if
        else
          print*, " Both absolute test and relative test are failed."
          print*, " Since it is severe, and some bugs might be hidden"
     &         //  "."
          print*, " Please inform this result to the developper soon."
        end if
      end if

      print*,"-----------------------------------------------"


      deallocate(ww)

      return

      end subroutine w_test
