!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: CSTAB_mod
!
! Purpose
! =======
!
! CSTAB_mod mnages the modules for Cache stabilization routines.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module CSTAB_mod


      implicit none
      private

      public :: CSTAB_get_optdim
      public :: CSTAB_adjust_base
      public :: CSTAB_adjust_page
      public :: get_delta
      public :: CSTAB_round_offset
      public :: CSTAB_round_offset_h

      interface get_delta
        module procedure get_delta_r
        module procedure get_delta_h
      end interface get_delta
      interface CSTAB_adjust_base
        module procedure CSTAB_adjust_base_r
        module procedure CSTAB_adjust_base_h
      end interface CSTAB_adjust_base
      interface CSTAB_adjust_page
        module procedure CSTAB_adjust_page_r
        module procedure CSTAB_adjust_page_h
      end interface CSTAB_adjust_page

!--------*---------*---------*---------*---------*---------*---------*-*

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: CSTAB_get_optdim
!
! Purpose
! =======
!
! Return the recommendational dimension of the array to minimize the
! risk of Cache thrashing
!
!
! Arguments
! =========
!
! n_min   (input) integer
!         minimum dimension of the array
!
! n_unroll (input) integer
!          possible unroll factor
!
! delta_L1 (input) integer
!          window size of L1
!
! delta_L2 (input) integer
!          window size of L2
!
! n_opt   (output) integer
!         recommendational dimension of the array
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine CSTAB_get_optdim(n_min, n_unroll,
     &     delta_L1, delta_L2, n_opt)

      integer, intent(in)    :: n_min
      integer, intent(in)    :: n_unroll
      integer, intent(in)    :: delta_L1
      integer, intent(in)    :: delta_L2
      integer, intent(out)   :: n_opt

      include 'CSTAB.h'

      real(8), parameter     :: ONE = 1d0

      integer                :: n_delta
      integer                :: n_opt2
      integer                :: n_delta2
      integer                :: i,j,k


      n_opt=n_min

      do

        n_opt   = (n_opt-1)/L1_WINDOW+1
        n_opt   = (n_opt/2)*2+1
        n_opt   = n_opt*L1_WINDOW

        n_delta = 0

        do i=1,int((n_unroll*1.2d0-ONE)/L1_WAY+1)

          k=mod(i*n_opt+L1_LSIZE/2,L1_LSIZE)-L1_LSIZE/2
          if (abs(k)<=delta_L1/2) then
            n_delta=(delta_L1/2-k-1)/i+1
            goto 10000
          end if

        end do

        do i=1,int((n_unroll*1.2d0-ONE)/L2_WAY+1)

          k=mod(i*n_opt+L2_LSIZE/2,L2_LSIZE)-L2_LSIZE/2
          if (abs(k)<=delta_L2/2) then
            n_delta=(delta_L2/2-k-1)/i+1
            goto 10000
          end if

        end do

10000   continue

        if (n_delta==0) exit
        n_opt = n_opt + n_delta

      end do

      return

      end subroutine CSTAB_get_optdim
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine get_delta_r(a,b,offset)

      real(8), intent(in)    :: a(*)
      real(8), intent(in)    :: b(*)
      integer, intent(out)   :: offset
      integer, intrinsic     :: loc

      offset = INT(loc(a) - loc(b), 8)
      
      end subroutine get_delta_r

      subroutine get_delta_h(a,b,offset)

      complex(8), intent(in) :: a(*)
      complex(8), intent(in) :: b(*)
      integer, intent(out)   :: offset
      integer, intrinsic     :: loc
      
      offset = INT(loc(a) - loc(b), 8)
      
      end subroutine get_delta_h
!--------*---------*---------*---------*---------*---------*---------*-*

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: CSTAB_adjust_base
!
! Purpose
! =======
!
! Adjust the cache line position
!
!
! Arguments
! =========
!
! a       (input/output) real array, dimension(*)
!         top position of the array a
!
! b       (input/output) real array, dimension(*)
!         top position of the array b
!
! offset  (input/output) integer
!         offset to adjust arrays a and b
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine CSTAB_adjust_base_r(a,b,offset)

      real(8), intent(in)    :: a(*)
      real(8), intent(in)    :: b(*)
      integer, intent(out)   :: offset

      include 'CSTAB.h'


      call get_delta_r(a(1),b(1),offset)
      offset=(offset/8)

      if (offset>0) then
        offset=mod(L2_LSIZE-mod(+offset,L2_LSIZE),L2_LSIZE)
      else
        offset=mod(L2_LSIZE+mod(-offset,L2_LSIZE),L2_LSIZE)
      end if

      return

      end subroutine CSTAB_adjust_base_r

      subroutine CSTAB_adjust_base_h(a,b,offset)

      complex(8), intent(in) :: a(*)
      complex(8), intent(in) :: b(*)
      integer, intent(out)   :: offset
      integer                :: L2_LSIZE_h

      include 'CSTAB.h'


      call get_delta_h(a(1),b(1),offset)
      offset=(offset/16)
      L2_LSIZE_h = L2_LSIZE/2

      if(offset>0)then
        offset=MOD(L2_LSIZE_h-MOD(+offset,L2_LSIZE_h),L2_LSIZE_h)
      else
        offset=MOD(L2_LSIZE_h+MOD(-offset,L2_LSIZE_h),L2_LSIZE_h)
      endif

      return

      end subroutine CSTAB_adjust_base_h

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: CSTAB_adjust_page
!
! Purpose
! =======
!
! Adjust the page position
!
!
! Arguments
! =========
!
! a       (input/output) real array, dimension(*)
!         top position of the array a
!
! b       (input/output) real array, dimension(*)
!         top position of the array b
!
! offset  (input/output) integer
!         offset to adjust arrays a and b
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine CSTAB_adjust_page_r(a,b,offset)

      real(8)                :: a(*)
      real(8)                :: b(*)
      integer                :: offset

      include 'CSTAB.h'


      call get_delta_r(a(1),b(1),offset)
      offset=(offset/8)
      if (offset>0) then
        offset = mod(PAGE_LSIZE-mod(+offset,PAGE_LSIZE),PAGE_LSIZE)
      else
        offset = mod(PAGE_LSIZE+mod(-offset,PAGE_LSIZE),PAGE_LSIZE)
      end if

      return

      end subroutine CSTAB_adjust_page_r

      subroutine CSTAB_adjust_page_h(a,b,offset)

      complex(8)             :: a(*)
      complex(8)             :: b(*)
      integer                :: offset
      integer                :: PAGE_LSIZE_h

      include 'CSTAB.h'


      call get_delta_h(a(1),b(1),offset)
      offset=(offset/16)
      PAGE_LSIZE_h = PAGE_LSIZE/2
      if (offset>0) then
        offset= mod(PAGE_LSIZE_h-mod(+offset,PAGE_LSIZE_h),PAGE_LSIZE_h)
      else
        offset= mod(PAGE_LSIZE_h+mod(-offset,PAGE_LSIZE_h),PAGE_LSIZE_h)
      end if

      return

      end subroutine CSTAB_adjust_page_h

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: CSTAB_round_offset
!
! Purpose
! =======
!
! Round up the offset by pagesize
!
!
! Arguments
! =========
!
! offset  (input/output) integer
!         offset
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine CSTAB_round_offset(offset)

      integer, intent(inout) :: offset

      include 'CSTAB.h'


      offset = mod(offset,L2_LSIZE)

      return

      end subroutine CSTAB_round_offset

      subroutine CSTAB_round_offset_h(offset)

      integer, intent(inout) :: offset
      integer                :: L2_LSIZE_h

      include 'CSTAB.h'


      L2_LSIZE_h = L2_LSIZE/2
      offset = mod(offset,L2_LSIZE_h)

      return

      end subroutine CSTAB_round_offset_h

      end module CSTAB_mod
