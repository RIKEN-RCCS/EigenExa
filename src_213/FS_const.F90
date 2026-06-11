!>
!> @file   FS_MERGE_D.F90
!> @brief  subroutine FS_MERGE_D
!!>
!
!
!
!>module FS_const_mod
!>
!> @brief @n
!> Purpose @n
!> ======= @n
!> definitions of const parameters
!> 
module FS_const_mod

  implicit none

  public

  real(kind(0.0d0)), parameter :: ZERO   =  0.0D+00
  real(kind(0.0d0)), parameter :: HALF   =  0.5D+00
  real(kind(0.0d0)), parameter :: ONE    =  1.0D+00
  real(kind(0.0d0)), parameter :: TWO    =  2.0D+00
  real(kind(0.0d0)), parameter :: THREE  =  3.0D+00
  real(kind(0.0d0)), parameter :: FOUR   =  4.0D+00
  real(kind(0.0d0)), parameter :: FIVE   =  5.0D+00
  real(kind(0.0d0)), parameter :: SIX    =  6.0D+00
  real(kind(0.0d0)), parameter :: SEVEN  =  7.0D+00
  real(kind(0.0d0)), parameter :: EIGHT  =  8.0D+00
  real(kind(0.0d0)), parameter :: NINE   =  9.0D+00
  real(kind(0.0d0)), parameter :: TEN    =  1.0D+01

  real(kind(0.0d0)), parameter :: MHALF  =  -0.5D+00
  real(kind(0.0d0)), parameter :: MONE   =  -1.0D+00
  real(kind(0.0d0)), parameter :: MTWO   =  -2.0D+00


end module FS_const_mod

