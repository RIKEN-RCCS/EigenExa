subroutine eigen_libs_eigen_init()
  use eigen_libs_mod
  implicit none
  call eigen_init()
end subroutine eigen_libs_eigen_init

subroutine eigen_libs_eigen_free()
  use eigen_libs_mod
  implicit none
  call eigen_free()
end subroutine eigen_libs_eigen_free

subroutine eigen_libs0_eigen_show_version()
  use eigen_libs0_mod
  implicit none
  call eigen_show_version()
end subroutine eigen_libs0_eigen_show_version

integer function eigen_blacs_eigen_get_blacs_context()
  use eigen_blacs_mod
  implicit none
  eigen_blacs_eigen_get_blacs_context = eigen_get_blacs_context()
end function eigen_blacs_eigen_get_blacs_context

subroutine eigen_libs_eigen_s(n,nvec,a,lda,w,z,ldz,m_forward,m_backward,mode)
  use eigen_libs_mod
  implicit none
  integer, intent(in)           :: n
  integer, intent(in), optional :: nvec
  real(8), intent(inout)        :: a(lda,*)
  integer, intent(in)           :: lda
  real(8), intent(out)          :: w(*)
  real(8), intent(out)          :: z(ldz,*)
  integer, intent(in)           :: ldz
  integer, intent(in), optional :: m_forward
  integer, intent(in), optional :: m_backward
  character, intent(in), optional :: mode
  call eigen_s(n,nvec,a,lda,w,z,ldz,m_forward,m_backward,mode)
end subroutine eigen_libs_eigen_s

subroutine eigen_libs_eigen_sx(n,nvec,a,lda,w,z,ldz,m_forward,m_backward,mode)
  use eigen_libs_mod
  implicit none
  integer, intent(in)           :: n
  integer, intent(in), optional :: nvec
  real(8), intent(inout)        :: a(lda,*)
  integer, intent(in)           :: lda
  real(8), intent(out)          :: w(*)
  real(8), intent(out)          :: z(ldz,*)
  integer, intent(in)           :: ldz
  integer, intent(in), optional :: m_forward
  integer, intent(in), optional :: m_backward
  character, intent(in), optional :: mode
  call eigen_sx(n,nvec,a,lda,w,z,ldz,m_forward,m_backward,mode)
end subroutine eigen_libs_eigen_sx

subroutine eigen_libs_eigen_h(n,nvec,a,lda,w,z,ldz,m_forward,m_backward,mode)
  use eigen_libs_mod
  implicit none
  integer, intent(in)           :: n
  integer, intent(in), optional :: nvec
  complex(8), intent(inout)     :: a(lda,*)
  integer, intent(in)           :: lda
  real(8), intent(out)          :: w(*)
  complex(8), intent(out)       :: z(ldz,*)
  integer, intent(in)           :: ldz
  integer, intent(in), optional :: m_forward
  integer, intent(in), optional :: m_backward
  character, intent(in), optional :: mode
  call eigen_h(n,nvec,a,lda,w,z,ldz,m_forward,m_backward,mode)
end subroutine eigen_libs_eigen_h

subroutine eigen_libs0_eigen_get_version(version,date,vcode)
  use eigen_libs0_mod
  implicit none
  integer, intent(out) :: version
  character(*), intent(out), optional :: date
  character(*), intent(out), optional :: vcode
  call eigen_get_version(version,date,vcode)
end subroutine eigen_libs0_eigen_get_version

!subroutine eigen_libs_eigen_get_matdims(n,nx,ny,m_forward,m_backward,mode)
!  use eigen_libs_mod
!  implicit none
!  integer, intent(in)           :: n
!  integer, intent(out)          :: nx, ny
!  integer, intent(in), optional :: m_forward
!  integer, intent(in), optional :: m_backward
!  character, intent(in), optional :: mode
!  call eigen_get_matdims(n,nx,ny,m_forward,m_backward,mode)
!end subroutine eigen_libs_eigen_get_matdims

subroutine eigen_libs_eigen_get_matdims(n,nx,ny,m_forward,m_backward,mode)
  use eigen_libs_mod
  implicit none
  integer, intent(in)           :: n
  integer, intent(out)          :: nx, ny
  integer, intent(in), optional :: m_forward, m_backward
  character, intent(in), optional :: mode
  call eigen_get_matdims(n,nx,ny,m_forward,m_backward,mode)
end subroutine eigen_libs_eigen_get_matdims

integer function eigen_libs0_eigen_memory_internal(n,lda,ldz,m1_opt,m0_opt)
  use eigen_libs0_mod
  implicit none
  integer, intent(in)           :: n, lda, ldz
  integer, intent(in), optional :: m1_opt, m0_opt
  eigen_libs0_eigen_memory_internal = eigen_memory_internal(n,lda,ldz,m1_opt,m0_opt)
end function eigen_libs0_eigen_memory_internal

subroutine eigen_libs0_eigen_get_comm(comm,x_comm,y_comm)
  use eigen_libs0_mod
  implicit none
  integer, intent(out)   ::  comm
  integer, intent(out)   ::  x_comm
  integer, intent(out)   ::  y_comm
  call eigen_get_comm(comm,x_comm,y_comm)
end subroutine eigen_libs0_eigen_get_comm

subroutine eigen_libs0_eigen_get_procs(procs,x_procs,y_procs)
  use eigen_libs0_mod
  implicit none
  integer, intent(out)   ::  procs
  integer, intent(out)   ::  x_procs
  integer, intent(out)   ::  y_procs
  call eigen_get_procs(procs,x_procs,y_procs)
end subroutine eigen_libs0_eigen_get_procs

subroutine eigen_libs0_eigen_get_id(id,x_id,y_id)
  use eigen_libs0_mod
  implicit none
  integer, intent(out)   ::  id
  integer, intent(out)   ::  x_id
  integer, intent(out)   ::  y_id
  call eigen_get_id(id,x_id,y_id)
end subroutine eigen_libs0_eigen_get_id

integer function eigen_libs0_eigen_loop_start(istart,nnod,inod)
  use eigen_libs0_mod
  implicit none
  integer, intent(in)    ::  istart
  integer, intent(in)    ::  nnod
  integer, intent(in)    ::  inod
  eigen_libs0_eigen_loop_start = eigen_loop_start(istart,nnod,inod)
end function eigen_libs0_eigen_loop_start

integer function eigen_libs0_eigen_loop_end(iend,nnod,inod)
  use eigen_libs0_mod
  implicit none
  integer, intent(in)    ::  iend
  integer, intent(in)    ::  nnod
  integer, intent(in)    ::  inod
  eigen_libs0_eigen_loop_end = eigen_loop_end(iend,nnod,inod)
end function eigen_libs0_eigen_loop_end

subroutine eigen_libs0_eigen_loop_info(istart,iend,lstart,lend,nnod,inod)
  use eigen_libs0_mod
  implicit none
  integer, intent(in)    ::  istart
  integer, intent(in)    ::  iend
  integer, intent(out)   ::  lstart
  integer, intent(out)   ::  lend
  integer, intent(in)    ::  nnod
  integer, intent(in)    ::  inod
  call eigen_loop_info(istart,iend,lstart,lend,nnod,inod)
end subroutine eigen_libs0_eigen_loop_info

integer function eigen_libs0_eigen_translate_l2g(ictr,nnod,inod)
  use eigen_libs0_mod
  implicit none
  integer, intent(in)    ::  ictr
  integer, intent(in)    ::  nnod
  integer, intent(in)    ::  inod
  eigen_libs0_eigen_translate_l2g = eigen_translate_l2g(ictr,nnod,inod)
end function eigen_libs0_eigen_translate_l2g

integer function eigen_libs0_eigen_translate_g2l(ictr,nnod,inod)
  use eigen_libs0_mod
  implicit none
  integer, intent(in)    ::  ictr
  integer, intent(in)    ::  nnod
  integer, intent(in)    ::  inod
  eigen_libs0_eigen_translate_g2l = eigen_translate_g2l(ictr,nnod,inod)
end function

integer function eigen_libs0_eigen_owner_node(ictr,nnod,inod)
  use eigen_libs0_mod
  implicit none
  integer, intent(in)    ::  ictr
  integer, intent(in)    ::  nnod
  integer, intent(in)    ::  inod
  eigen_libs0_eigen_owner_node = eigen_owner_node(ictr,nnod,inod)
end function eigen_libs0_eigen_owner_node

integer function eigen_libs0_eigen_owner_index(ictr,nnod,inod)
  use eigen_libs0_mod
  implicit none
  integer, intent(in)    ::  ictr
  integer, intent(in)    ::  nnod
  integer, intent(in)    ::  inod
  eigen_libs0_eigen_owner_index = eigen_owner_index(ictr,nnod,inod)
end function eigen_libs0_eigen_owner_index

integer function eigen_libs0_eigen_convert_id_xy2w(xinod,yinod)
  use eigen_libs0_mod
  implicit none
  integer, intent(in)    ::  xinod
  integer, intent(in)    ::  yinod
  eigen_libs0_eigen_convert_id_xy2w = eigen_convert_id_xy2w(xinod,yinod)
end function eigen_libs0_eigen_convert_id_xy2w

subroutine eigen_libs0_eigen_convert_id_w2xy(inod,xinod,yinod)
  use eigen_libs0_mod
  implicit none
  integer, intent(in)    ::  inod
  integer, intent(out)   ::  xinod
  integer, intent(out)   ::  yinod
  call eigen_convert_id_w2xy(inod,xinod,yinod)
end subroutine eigen_libs0_eigen_convert_id_w2xy

subroutine eigen_libs0_eigen_get_errinfo(info)
  use eigen_libs0_mod
  implicit none
  integer, intent(inout) ::  info
  call eigen_get_errinfo(info)
end subroutine eigen_libs0_eigen_get_errinfo
