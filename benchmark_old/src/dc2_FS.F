!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: dc2_mod
!
! Purpose
! =======
!
! dc2_mod manages the modules for divide and conquer algrithm for
! solving a symmetric tridiagonal matrix.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module dc2_FS_mod

      use, intrinsic :: iso_c_binding

      implicit none
      private

      public :: dc2_FS
      public :: dc2_FS_f32
      interface 
         subroutine dc2_FS(n, nvec, d, e, z, ldz, info, ret) 
     &      bind(c, name="dc2_FS_f64")

            use, intrinsic :: iso_c_binding
            integer(c_int), intent(in),value :: n, nvec, ldz
            real(c_double), intent(inout) :: d(1:n)
            real(c_double), intent(inout) :: e(1:n-1)
            real(c_double), intent(out)   :: z(1:ldz,*)
            integer(c_int), intent(out) :: info
            real(c_double), intent(out) :: ret
         end subroutine

         subroutine dc2_FS_f32(n, nvec, d, e, z, ldz, info, ret) 
     &      bind(c, name="dc2_FS_f32")
       
            use, intrinsic :: iso_c_binding
            integer(c_int), intent(in),value :: n, nvec, ldz
            real(c_float), intent(inout) :: d(1:n)
            real(c_float), intent(inout) :: e(1:n-1)
            real(c_float), intent(out)   :: z(1:ldz,*)
            integer(c_int), intent(out) :: info
            real(c_float), intent(out) :: ret
         end subroutine
      end interface

      contains
      end module dc2_FS_mod
