       subroutine  eigen_hrd_u(
     &               a, nm,
     &               u_x, u_y, nv,
     &               u_t, v_t, i, beta, e)
!$     use OMP_LIB
       use eigen_devel_mod
       use eigen_libs_mod
       use comm_mod
       implicit NONE

       integer, intent(in)    ::  nm, nv, i
       complex(8), intent(inout) ::  a(1:nm)
       complex(8), intent(inout) ::  u_x(1:nv), u_y(1:nv)
       complex(8), intent(inout) ::  u_t(*), v_t(*)
       complex(8), intent(out)   ::  beta
       real(8), intent(out)   ::  e(*)

       complex(8)             ::  a_n, u_n
       real(8)                ::  anorm2, g_n

       complex(8)             ::  t, s
       real(8)                ::  tr(4), sr(4)

       integer                ::  x_owner_nod, x_pos
       integer                ::  y_owner_nod, y_pos
       integer                ::  j_1, j_2, j_3
       integer                ::  jj_1, jj_2, jj_3
       integer                ::  L


       L = i-1

       x_owner_nod = eigen_owner_node (L, 'X')
       y_owner_nod = eigen_owner_node (i, 'Y')

       x_pos       = eigen_translate_g2l(L, 'X')
       y_pos       = eigen_translate_g2l(i, 'Y')

       j_2         = eigen_loop_start(1, 'X')
       j_3         = eigen_loop_end  (L, 'X')
!
! u=...
!
       if ( y_owner_nod == y_inod ) then
          anorm2 = ZERO
          do  j_1=j_2,j_3
             t = a(j_1)
             anorm2   =  anorm2 + DBLE(t * dconjg(t))
             u_x(j_1) =  t
          end do ! j_1
          u_x(x_pos+1) =  DCMPLX(anorm2, ZERO)
       end if

       call bcast_cmplx(u_x(1), x_pos+1, y_owner_nod, y_COMM_WORLD)

       anorm2 = DBLE(u_x(x_pos+1))
       u_x(x_pos+1) = ZEROZ

       if ( x_owner_nod == x_inod ) then
          a_n =  u_x(x_pos)
       else
          a_n =  ZEROZ
       end if

       tr(1) =  anorm2
       tr(2) =  DBLE(a_n)
       tr(3) =  DIMAG(a_n)
       call reduce_dbl(tr(1), sr(1), 3, 1, x_COMM_WORLD)
       anorm2 =  tr(1)
       a_n    =  DCMPLX(tr(2), tr(3))

       if ( anorm2 /= ZERO ) then
!
! Similar to NAG and LAPACK approaches as \gamma:=-1, but
! w:=x+\nu e_1 and \sigma:=1/beta=1/(\xi_1+\nu)\nu are utilized,
! in short, just differs from the normalization fator from them.
!
! Theoretical points are from:
! [0] The Wilksinson Book
! [1] D. D. Moriison, J. of the ACM, 7(2) 1960, though it describes
! the triagularization method
! [2] A. R. Gourlay, the Computer Journal, 13(4), 1970
! [3] R. Lehoucq, ACM TOMS, 22(4) 1996
!
          g_n   = -SIGN(SQRT(anorm2), DBLE(a_n))
          u_n   =  a_n - g_n
          beta  = -u_n * g_n
          e (i) = g_n

          if ( x_owner_nod == x_inod ) then
             u_x(x_pos) =  u_n
             if ( y_owner_nod == y_inod ) then
                a(x_pos) =  u_x(x_pos)
             end if
          else
             if ( j_3 < x_pos ) then
                u_x(j_3+1:x_pos) = ZEROZ
             end if
          end if

          call datacast_cmplx(u_y(1), u_x(1), u_t(1), v_t(1), x_pos)

       else

          beta  =  ONEZ
          e (i) =  ZERO

       end if

       return
       end subroutine  eigen_hrd_u

