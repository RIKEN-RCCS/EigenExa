       subroutine eigen_hrd_vo(
     &               u_x, v_x, v_y,
     &               ux, vx, nv,
     &               u_t, v_t, beta,
     &               i, i_base, m0)
!$     use OMP_LIB
       use eigen_devel_mod
       use eigen_libs_mod
       implicit none
!
       integer, intent(in)    ::  nv, i, i_base, m0
       complex(8), intent(in)    ::  u_x(*)
       complex(8), intent(out)   ::  v_x(*), v_y(*)
       complex(8), intent(in)    ::  ux(1:nv,*)
       complex(8), intent(in)    ::  vx(1:nv,*)
       complex(8), intent(inout) ::  u_t(*)
       complex(8), intent(inout) ::  v_t(*)
       complex(8), intent(in)    ::  beta
!
       integer                ::  j, k, L, n, LL
       integer                ::  k_1, k_2, k_3
!
       integer                ::  i_1, i_2, i_3, i_4
       integer                ::  j_1, j_2, j_3, j_4
       integer                ::  l_1, l_2, l_3, l_4
       integer                ::  jj_1, jj_2, jj_3, jj_4
       integer                ::  LX
!
       include 'CSTAB.h'
!
       complex(8)                ::  w0
       complex(8)                ::  u0, v0
       complex(8)                ::  u1, v1
       complex(8)                ::  u2, v2
       complex(8)                ::  ux0, vx0
       complex(8)                ::  ux1, vx1
       complex(8)                ::  ux2, vx2
!
       complex(8)                ::  alpha, prod_uv
       integer                ::  x_pos, x_owner_nod
       integer                ::  y_pos, y_owner_nod



          k_1 = i - i_base
          k_2 = m0

          L = i-1
          n  = eigen_translate_g2l(L, 'X')

! FOR attention to unexpected overflow or NAN
          j_3 = eigen_loop_end    (L, 'X')
          if ( j_3 < n ) then
             v_x(j_3+1:n) = ZEROZ ! in case
          end if

          if ( k_2 <= k_1 ) then

!$OMP BARRIER

          else
!
! v=v-(UV+VU)u
!
             l_2 = k_2-k_1
             do j=1,l_2*2+0
                u_t(j) = ZEROZ
             end do

             l_4 = MOD(k_2-k_1, 3)+k_1+1
             LX  = L1_LSIZE*L1_WAY/16

             LL = (n-1)/y_nnod+1
             LL = ((LL-1)/2+1)*2

             jj_2 = 1
             jj_3 = n

             k_3 = eigen_loop_end(L, 'X')

             do jj_1=jj_2,jj_3,LX
                j_2 = jj_1; j_3 = MIN(jj_1+LX-1, jj_3)

                if(l_4-1==k_1+1)then
                   l_1 = k_1+1                           ! 0

                   j = l_1-k_1
                   u0 = u_t(2*(j+0)-1+0)
                   v0 = u_t(2*(j+0)-0+0)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
!IBM* INDEPENDENT
!IBM* ASSERT(NODEPS)
                   do j_1=j_2,j_3
                      w0 = u_x(j_1)
                      u0 = u0+dconjg(vx(j_1,l_1+0))*w0
                      v0 = v0+dconjg(ux(j_1,l_1+0))*w0
                   end do! j_1
                   u_t(2*(j+0)-1+0) = u0
                   u_t(2*(j+0)-0+0) = v0
                end if
                if(l_4-2==k_1+1)then
                   l_1 = k_1+1                           ! 1

                   j = l_1-k_1
                   u0 = u_t(2*(j+0)-1+0)
                   v0 = u_t(2*(j+0)-0+0)
                   u1 = u_t(2*(j+1)-1+0)
                   v1 = u_t(2*(j+1)-0+0)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
!IBM* INDEPENDENT
!IBM* ASSERT(NODEPS)
                   do j_1=j_2,j_3
                      w0 = u_x(j_1)
                      u0 = u0+dconjg(vx(j_1,l_1+0))*w0
                      v0 = v0+dconjg(ux(j_1,l_1+0))*w0
                      u1 = u1+dconjg(vx(j_1,l_1+1))*w0
                      v1 = v1+dconjg(ux(j_1,l_1+1))*w0
                   end do! j_1
                   u_t(2*(j+0)-1+0) = u0
                   u_t(2*(j+0)-0+0) = v0
                   u_t(2*(j+1)-1+0) = u1
                   u_t(2*(j+1)-0+0) = v1
                end if
                do l_1=l_4,k_2,3                  ! 2

                   j = l_1-k_1
                   u0 = u_t(2*(j+0)-1+0)
                   v0 = u_t(2*(j+0)-0+0)
                   u1 = u_t(2*(j+1)-1+0)
                   v1 = u_t(2*(j+1)-0+0)
                   u2 = u_t(2*(j+2)-1+0)
                   v2 = u_t(2*(j+2)-0+0)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
!IBM* INDEPENDENT
!IBM* ASSERT(NODEPS)
                   do j_1=j_2,j_3
                      w0 = u_x(j_1)
                      u0 = u0+dconjg(vx(j_1,l_1+0))*w0
                      v0 = v0+dconjg(ux(j_1,l_1+0))*w0
                      u1 = u1+dconjg(vx(j_1,l_1+1))*w0
                      v1 = v1+dconjg(ux(j_1,l_1+1))*w0
                      u2 = u2+dconjg(vx(j_1,l_1+2))*w0
                      v2 = v2+dconjg(ux(j_1,l_1+2))*w0
                   end do! j_1
                   u_t(2*(j+0)-1+0) = u0
                   u_t(2*(j+0)-0+0) = v0
                   u_t(2*(j+1)-1+0) = u1
                   u_t(2*(j+1)-0+0) = v1
                   u_t(2*(j+2)-1+0) = u2
                   u_t(2*(j+2)-0+0) = v2
                end do! l_1

             end do! jj_1

!$OMP BARRIER

             j_3 = (k_2-k_1)*2+0
             call reduce_cmplx(u_t, v_t, j_3, 1, x_COMM_WORLD)

             jj_2 = 1
             jj_3 = n

             do jj_1=jj_2,jj_3,LX
                j_2 = jj_1; j_3 = MIN(jj_1+LX-1, jj_3)

                if(l_4-1==k_1+1)then
                   l_1 = k_1+1                           ! 0

                   j = l_1-k_1

                   u0 = u_t(2*(j+0)-1+0)
                   v0 = u_t(2*(j+0)-0+0)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
!IBM* INDEPENDENT
!IBM* ASSERT(NODEPS)
                   do j_1=j_2,j_3
                      w0 = v_x(j_1)
                      ux0 = ux(j_1,l_1+0)
                      vx0 = vx(j_1,l_1+0)
                      w0 = w0
     &                   -ux0*u0
     &                   -vx0*v0
                      v_x(j_1) = w0
                   end do! j_1
                end if
                if(l_4-2==k_1+1)then
                   l_1 = k_1+1                           ! 1

                   j = l_1-k_1

                   u0 = u_t(2*(j+0)-1+0)
                   v0 = u_t(2*(j+0)-0+0)
                   u1 = u_t(2*(j+1)-1+0)
                   v1 = u_t(2*(j+1)-0+0)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
!IBM* INDEPENDENT
!IBM* ASSERT(NODEPS)
                   do j_1=j_2,j_3
                      w0 = v_x(j_1)
                      ux0 = ux(j_1,l_1+0)
                      vx0 = vx(j_1,l_1+0)
                      w0 = w0
     &                   -ux0*u0
     &                   -vx0*v0
                      ux1 = ux(j_1,l_1+1)
                      vx1 = vx(j_1,l_1+1)
                      w0 = w0
     &                   -ux1*u1
     &                   -vx1*v1
                      v_x(j_1) = w0
                   end do! j_1
                end if
                do l_1=l_4,k_2,3                  ! 2

                   j = l_1-k_1

                   u0 = u_t(2*(j+0)-1+0)
                   v0 = u_t(2*(j+0)-0+0)
                   u1 = u_t(2*(j+1)-1+0)
                   v1 = u_t(2*(j+1)-0+0)
                   u2 = u_t(2*(j+2)-1+0)
                   v2 = u_t(2*(j+2)-0+0)
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
!IBM* INDEPENDENT
!IBM* ASSERT(NODEPS)
                   do j_1=j_2,j_3
                      w0 = v_x(j_1)
                      ux0 = ux(j_1,l_1+0)
                      vx0 = vx(j_1,l_1+0)
                      w0 = w0
     &                   -ux0*u0
     &                   -vx0*v0
                      ux1 = ux(j_1,l_1+1)
                      vx1 = vx(j_1,l_1+1)
                      ux2 = ux(j_1,l_1+2)
                      vx2 = vx(j_1,l_1+2)
                      w0 = w0
     &                   -ux1*u1
     &                   -vx1*v1
     &                   -ux2*u2
     &                   -vx2*v2
                      v_x(j_1) = w0
                   end do! j_1
                end do! l_1

             end do! jj_1

          end if

             j_3 = eigen_loop_end(L, 'X')
             prod_uv = ZEROZ
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
!IBM* INDEPENDENT
!IBM* ASSERT(NODEPS)
             do j_1=1,j_3
                prod_uv = prod_uv + v_x(j_1) * dconjg(u_x(j_1))
             end do! j_1

             call reduce_cmplx(prod_uv, v_t, 1, 1, x_COMM_WORLD)

          if ( beta /= ZEROZ ) then

             alpha = prod_uv/(2*beta)
!DIR$ IVDEP
!IBM* INDEPENDENT
!IBM* ASSERT(NODEPS)
             do j_1=1,j_3
                v_x(j_1) = (v_x(j_1)-alpha*u_x(j_1))/dconjg(beta)
             end do! j_1

          end if

             x_pos = eigen_translate_g2l(L, 'X')

             call datacast_cmplx(v_y(1), v_x(1), u_t(1), v_t(1), x_pos)

! FOR attention to unexpected overflow or NAN
             j_3 = eigen_loop_end     (L, 'X')
             n   = eigen_translate_g2l(L, 'X')
             if ( j_3 < n ) then
                v_x(j_3+1:n) = ZEROZ ! in case
             end if

             j_3 = eigen_loop_end     (L, 'Y')
             n   = eigen_translate_g2l(L, 'Y')
             if ( j_3 < n ) then
                v_y(j_3+1:n) = ZEROZ ! in case
             end if


      return
      end subroutine  eigen_hrd_vo

