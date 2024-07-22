      module eigen_HRBAK_mod
!$    use OMP_LIB
      implicit NONE
!=
         real(8) :: TRBK_TIME_BCAST
         real(8) :: TRBK_TIME_REDUC
         real(8) :: TRBK_TIME_FR
         real(8) :: TRBK_TIME_TRBK1
         real(8) :: TRBK_TIME_TRBK1x
         real(8) :: TRBK_TIME_TRBK1y
         real(8) :: TRBK_TIME_TRBK2

!$       integer(kind=omp_lock_kind) :: TRBK_LOCK
!$       integer(8)                  :: TRBK_MASK(1:2)

         integer :: DO_OVERLAP_BCAST_LEVEL
         integer :: TRBK_TIME_COUNTER
         integer :: TRBK_TIME_INTERVAL
         integer :: TRBK_TIME_NEXT
         integer :: TRBK_SWITCHED
!=
!     This must be attributed with pointer for Fortran 200x late
         complex(8), pointer :: HRBK_buf(:)

!=
      end module eigen_HRBAK_mod
