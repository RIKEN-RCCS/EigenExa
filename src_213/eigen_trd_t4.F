!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_trd_t4_mod
!
! Purpose
! =======
!
! eigen_trd_t4_mod manages the modules for the generation of a seed
! Householder reflector.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module eigen_trd_t4_mod

      use eigen_libs_mod
      use comm_mod
      use eigen_devel_mod
!$    use omp_lib

      implicit none
      private


      public :: eigen_trd_compute_u

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_trd_compute_u
!
! Purpose
! =======
!
! Generator of a Householder reflector
!
!
! Arguments
! =========
!
! d_out  (output) real(8)
!        d_out contains diagonal elements of the output tri-diagonal
!        matrix.
!
! u_x    (input/output) real(8) array, dimension(nv)
!        On entry, u_x contains the specfic row-vector from a to be
!        broadcasted.
!        On exit, u_x returns the reflector vector broadcasted.
!
! u_y    (input/output) real(8) array, dimension(nv)
!        On entry, u_y contains the specfic row-vector from a to be
!        broadcasted.
!        On exit, u_y returns the reflector vector broadcasted.
!
! nv     (input) integer
!        The leading dimension of the working arrays, u_x, u_y, ...
!
! u_t    (input/output) real(8) array, dimension(*)
!        working buffer
!
! v_t    (input/output) real(8) array, dimension(*)
!        working buffer
!
! i      (input) integer
!        The global loop counter
!
! flag_concat (input) logical
!        Flag for concatenation of multiple broadcasts
!
! flag_sample (input) logical
!        Flag for sampleing of cost of broadcasts
!
! v_x    (input) real(8) array, dimension(*)
!        v_x contains a counterpart vector generated in last iteration
!
! v_y    (output) real(8) array, dimension(*)
!        v_y contains a counterpart vector generated in last iteration
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_trd_compute_u(
     &     d_out,
     &     u_x, u_y, nv,
     &     u_t, v_t, i, flag_concat, flag_sample, v_x, v_y)

      real(8), intent(out)   ::  d_out
      real(8), intent(inout) ::  u_x(1:nv)
      real(8), intent(inout) ::  u_y(1:nv)
      integer, intent(in)    ::  nv
      real(8), intent(inout) ::  u_t(*)
      real(8), intent(inout) ::  v_t(*)
      integer, intent(in)    ::  i
      logical, intent(in)    ::  flag_concat
      logical, intent(in)    ::  flag_sample
      real(8), intent(inout) ::  v_x(*)
      real(8), intent(out)   ::  v_y(*)

      integer                ::  x_pos, y_pos
      integer                ::  j_1, j_2, j_3
      integer                ::  L


      L = i-1

!
! Pick up the diagonal element from the target vector in panel (u_x)
! into d_out, then nullfy the element on the panel.
!
      j_1 = eigen_owner_index(i, 'X')
      if (j_1 > 0) then
        d_out = u_x(j_1)
        u_x(j_1) = ZERO
      else
        d_out = ZERO
      end if


      if (flag_overlap .or. flag_oncache) then

        x_pos = eigen_translate_g2l(L, 'X')
        call datacast_dbl(u_y(1), u_x(1), v_t(1), u_t(1), L, 1)

        call eigen_vector_zeropad_x(u_x(1), L)
        call eigen_vector_zeropad_y(u_y(1), L)

      else

        if (flag_concat) then

          x_pos = eigen_translate_g2l(i, 'X')

          if (.not. flag_oncache) then
            if (flag_sample) then
              timer_t1 = eigen_get_wtime()
            end if
          end if
          call datacast_dbl2(u_y(1), v_y(1), u_x(1), v_x(1),
     &                       v_t(1), u_t(1), i, 1)
          if (.not. flag_oncache) then
          if (flag_sample) then
            timer_t1 = eigen_get_wtime() - timer_t1 - Bcast_Overhead_x
            timer_t1 = timer_t1/2 + Bcast_cont_Overhead_x
          end if
          end if

        else

          x_pos = eigen_translate_g2l(L, 'X')
          call datacast_dbl(u_y(1), u_x(1),
     &                      v_t(1), u_t(1), L, 1)

        end if

        call eigen_vector_zeropad_x(u_x(1), L)
        call eigen_vector_zeropad_y(u_y(1), L)
        if (flag_concat) then
          call eigen_vector_zeropad_y(v_y(1), i)
        end if

      end if

!$OMP BARRIER

      if (.not. flag_oncache) then
      if (flag_overlap) then

        if (flag_concat) then

          x_pos = eigen_translate_g2l(i, 'X')

          if (flag_sample) then
            timer_t1 = eigen_get_wtime()
          end if
          call datacast_dbl(v_y(1), v_x(1),
     &                      v_y(1-nv), u_t(1), i, 1)
          if (flag_sample) then
            timer_t1 = eigen_get_wtime() - timer_t1
          end if

          call eigen_vector_zeropad_y(v_y(1), i)

        end if

      end if
      end if

      return

      end subroutine eigen_trd_compute_u

      end module eigen_trd_t4_mod
