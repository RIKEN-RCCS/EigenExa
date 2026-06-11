module blacs_c_wrapper
  use iso_c_binding
  implicit none
contains

  subroutine blacs_get_c(context, what, ictxt) bind(C, name="blacs_get_c")
    integer(c_int), value :: context, what
    integer(c_int) :: ictxt
    call BLACS_GET(context, what, ictxt)
  end subroutine blacs_get_c

  subroutine blacs_gridmap_c(ictxt, usermap, ldup, nprow, npcol) bind(C, name="blacs_gridmap_c")
    integer(c_int), value :: ictxt
    integer(c_int), dimension(*) :: usermap
    integer(c_int), value :: ldup, nprow, npcol
    call BLACS_GRIDMAP(ictxt, usermap, ldup, nprow, npcol)
  end subroutine blacs_gridmap_c

  subroutine blacs_gridexit_c(ictxt) bind(C, name="blacs_gridexit_c")
    integer(c_int), value :: ictxt
    call BLACS_GRIDEXIT(ictxt)
  end subroutine blacs_gridexit_c

end module blacs_c_wrapper

