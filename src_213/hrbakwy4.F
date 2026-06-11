       subroutine eigen_hrbakwyx(n, nvec, a, nma0, z, nmz0, beta, m0,
     &            iblk)
!$     use OMP_LIB
       use MPI
       use eigen_devel_mod
       use eigen_libs_mod
       use eigen_hrbak_mod
       use CSTAB_mod
       implicit NONE

       integer, intent(in)    :: n, nvec, nma0, nmz0, m0, iblk
       complex(8), intent(in)    :: a(*)
       complex(8), intent(out)   :: z(*)
       complex(8), intent(inout) :: beta(*)

       complex(8), allocatable   :: v1(:), ss1(:)
       complex(8), allocatable   :: v2(:), ss2(:)
       complex(8), allocatable   :: v3(:), ss3(:)
       complex(8), allocatable   :: tt(:)

       integer                :: nma, nmz, m
       integer                :: na, nm, kx
       integer                :: i_v1, i_v2, i_v3
       integer                :: i_s1, i_s2, i_s3
       integer                :: i_t
       integer                :: i_stat, ierr

       real(8) :: d1, d2, d3

       include 'CSTAB.h'
*
       nma  = nma0
       nmz  = nmz0
       m    = MIN(nsm, m0)
       if ( m < 1 ) m = 1
*
       call eigen_timer_reset( )
*
       na   = (n-1)/y_nnod+1
       na   = na  +MOD(na  -1,2)
       call CSTAB_get_optdim(nma,9,16*4,16*6,nm)
*
#if BOOST_BY_CACHE_ALIGNMENT
       allocate(
     &          v1(1:MAX(nm*m,n)+n_columns),
     &          ss1(1:na*m+ns0+n_columns),
     &          v2(1:MAX(nm*m,n)+n_columns),
     &          ss2(1:na*m+ns0+n_columns),
     &          v3(1:MAX(nm*m,n)+n_columns),
     &          ss3(1:na*m+ns0+n_columns),
     &          tt(1:na*m+ns0+n_columns),
     &          stat=i_stat)
#else
      allocate(
     &     v1(1:max(nm*m,n)),
     &     ss1(1:na*m+ns0),
     &     v2(1:max(nm*m,n)),
     &     ss2(1:na*m+ns0),
     &     v3(1:max(nm*m,n)),
     &     ss3(1:na*m+ns0),
     &     tt(1:na*m+ns0),
     &     stat=i_stat)
#endif
       if ( i_stat /= 0 ) then
         call eigen_abort( "Memory allocation error [HBK].", 1 )
       end if
*
       v1 = ZERO
       ss1 = ZERO
       v2 = ZERO
       ss2 = ZERO
       v3 = ZERO
       ss3 = ZERO
       tt = ZERO
*
#if BOOST_BY_CACHE_ALIGNMENT
       call CSTAB_adjust_base(v1,  z, i_v1)
       call CSTAB_adjust_base(ss1, z, i_s1)
       call CSTAB_adjust_base(v2,  z, i_v2)
       call CSTAB_adjust_base(ss2, z, i_s2)
       call CSTAB_adjust_base(v3,  z, i_v3)
       call CSTAB_adjust_base(ss3, z, i_s3)
       call CSTAB_adjust_base(tt,  z, i_t)
!          kx = (L1_WINDOW/8)
          kx = (L1_WINDOW/16)
!     &           +(L1_WINDOW)
!     &           +(L1_LSIZE/8)
     &           +(L1_LSIZE)
!     &           +(L2_LSIZE/8)
     &           +(L2_LSIZE/16)
          kx = kx / 2
       i_v1 = i_v1 + kx*3
       i_s1 = i_s1 + kx*1
       i_v2 = i_v2 + kx*5
       i_s2 = i_s2 + kx*1
       i_v3 = i_v3 + kx*7
       i_s3 = i_s3 + kx*1
       i_t = i_t + kx*1
       call CSTAB_round_offset_h(i_v1)
       call CSTAB_round_offset_h(i_s1)
       call CSTAB_round_offset_h(i_v2)
       call CSTAB_round_offset_h(i_s2)
       call CSTAB_round_offset_h(i_v3)
       call CSTAB_round_offset_h(i_s3)
       call CSTAB_round_offset_h(i_t)
#else
      i_v1=0
      i_s1=0
      i_v2=0
      i_s2=0
      i_v3=0
      i_s3=0
      i_t=0
#endif
*
       allocate( HRBK_buf(1:m*nm), stat=i_stat )
       if (i_stat /= 0) then
         call eigen_abort( "Memory allocation error [HBK].", 1 )
       end if
*
       call MPI_Barrier(TRD_COMM_WORLD,ierr)
       d1 = eigen_get_wtime()
!$OMP PARALLEL
       call hrbakwy_body( n, nvec,
     $      a, nma,
     $      z, nmz,
     $      beta,
     $      v1(1+i_v1:), v2(1+i_v2:), v3(1+i_v3:), nm, m,
     $      ss1(1+i_s1:), ss2(1+i_s2:), ss3(1+i_s3:),
     $      tt(1+i_t:), iblk, na )
!$OMP END PARALLEL
       call MPI_Barrier(TRD_COMM_WORLD,ierr)
       d2 = eigen_get_wtime()
*
       deallocate( HRBK_buf )
*
       deallocate( v1 )
       deallocate( ss1 )
       deallocate( v2 )
       deallocate( ss2 )
       deallocate( v3 )
       deallocate( ss3 )
       deallocate( tt )
*
       d3 = eigen_timer_print(
     &      'EigenExa (Back-Hermete-transformation)' )

#if TIMER_PRINT
       IF(TRD_inod==1)print*,"HRBAKWY",d2-d1
#endif
*
       return
       end
*
*
*
       subroutine hrbakwy_body(n, nvec, a, nma, z, nmz, beta,
     &             v1, v2, v3, nm, m, ss1, ss2, ss3, tt, iblk, nss)
!$     use OMP_LIB
       use MPI
       use eigen_devel_mod
       use eigen_libs_mod
       use comm_mod
       use eigen_HRBAK_mod
       implicit NONE

       integer, intent(in)    :: n, nvec, nma, nmz
       integer, intent(in)    :: nm, m, iblk, nss
       complex(8), intent(in)    :: a(1:nma,*)
       complex(8), intent(inout) :: z(1:nmz,*)
       complex(8), intent(inout) :: beta(1:n)

       complex(8), intent(out)   :: v1(1:nm,*), ss1(*)
       complex(8), intent(out)   :: v2(1:nm,*), ss2(*)
       complex(8), intent(out)   :: v3(1:nm,*), ss3(*)
       complex(8), intent(out)   :: tt(*)
*
       integer                :: nodes(0:nsm-1)

       complex(8), allocatable   :: wk(:)
       integer                :: nx, lwk, ierr

       integer                :: i, L, mode
       integer                :: x_pos, x_root
       integer                :: i_1, i_2, i_3, i_4
       integer                :: j_1, j_2, j_3, j_4
       complex(8)             :: s0, s1, s2, s3
       real(8)                :: d0, d1, d2, ds, de, dcom, dx, dy, dz

*
       dx = 0D0; dy = 0D0; dz = 0D0

!$OMP MASTER
       call MPI_Barrier(TRD_COMM_WORLD,ierr)
!$OMP END MASTER
!$OMP BARRIER
!$OMP MASTER
#if TIMER_PRINT
       d1 = eigen_get_wtime()
#endif
       dcom=0.0D0+00
*-
#if AT_BCAST_OVERLAP
!$     DO_OVERLAP_BCAST_LEVEL  = 2
!$     TRBK_TIME_COUNTER  = 0
!$     TRBK_TIME_INTERVAL = 1
!$     TRBK_TIME_NEXT     = 2
!$     TRBK_SWITCHED      = 0
!$     TRBK_TIME_TRBK1    = 0D0
!$     TRBK_TIME_TRBK1x   = 0D0
!$     TRBK_TIME_BCAST    = 0D0
!$     TRBK_TIME_REDUC    = 1D0
!$     TRBK_TIME_FR       = 0D0
!===
!$     call sync_other_than_master_init( TRBK_LOCK, TRBK_MASK )
!===
#else
!$     DO_OVERLAP_BCAST_LEVEL  = 0
#endif
*-
       lwk=MAX(((m-1)/y_nnod+1)*((n-1)/x_nnod+1),n)
       allocate(wk(lwk))
!$OMP END MASTER
!$OMP BARRIER
*
       i_2 = eigen_loop_start(1+iblk, 'Y')
       i_3 = eigen_loop_end  (n,      'Y')
*
!$OMP DO
       do i=1,n
          l = i-iblk
          if ( i >= 1+iblk .and.
     &       mod(i-1, y_nnod)+1 == y_inod .and.
     &       mod(l-1, x_nnod)+1 == x_inod ) then
            i_1 = (i-y_inod)/y_nnod+1
            j_1 = (l-x_inod)/x_nnod+1
            beta(i) = a(j_1,i_1)*beta(i)
          else
            beta(i) = ZEROZ
          end if
       end do! i_1
!$OMP END DO
!$OMP MASTER
       call reduce_cmplx( beta, wk, n, 1, TRD_COMM_WORLD )
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO
       do i=1,n
          if ( beta(i) == ZEROZ ) then
             beta(i) = ONEZ
          endif
       end do! i_1
!$OMP END DO
*
       i_2 = eigen_loop_start(1,    'Y')
       i_3 = eigen_loop_end  (nvec, 'Y')
*
*----------------------------------------------------------------
*
       nx = MIN(MOD(n-(1+iblk)+1,m)+(1+iblk)-1,n)
*
       do i=(1+iblk),nx

!$OMP BARRIER
*
          if ( beta(i) == ZEROZ ) cycle
*
          j_2 = eigen_loop_start(1,      'X')
          j_3 = eigen_loop_end  (i-iblk, 'X')
          i_4=MOD(i_3-i_2+1,4)+i_2
*
!$OMP MASTER
          ds = eigen_get_wtime()
          nodes(0) = eigen_owner_node(i, 'Y')
          if ( nodes(0) == y_inod ) then
             i_1 = eigen_translate_g2l(i, 'Y')
             do j_1=j_2,j_3
                v1(j_1,1) = a(j_1,i_1)
             end do! j_1
          end if
          call bcast_cmplx(v1(j_2,1), j_3-j_2+1, nodes(0), y_COMM_WORLD)

          de = eigen_get_wtime()
          dcom = dcom + (de - ds)
!$OMP END MASTER
*
!$OMP BARRIER
*
!$OMP MASTER
          call barrier( y_COMM_WORLD )
          call barrier( x_COMM_WORLD )
!$OMP END MASTER
!$    if ( omp_get_num_threads() == 1 .OR.
!$   +     omp_get_thread_num() == 1 ) then
          if ( i_4 == i_2 + 1 ) then
             i_1 = i_2
             s0=0.0D+00
             do j_1=j_2,j_3
                s0=s0+dconjg(v1(j_1,1))*z(j_1,i_1+0)
             end do! j_1
             ss1(i_1+0)=s0
          end if
          if ( i_4 == i_2 + 2 ) then
             i_1 = i_2
             s0=0.0D+00
             s1=0.0D+00
             do j_1=j_2,j_3
                s0=s0+dconjg(v1(j_1,1))*z(j_1,i_1+0)
                s1=s1+dconjg(v1(j_1,1))*z(j_1,i_1+1)
             end do! j_1
             ss1(i_1+0)=s0
             ss1(i_1+1)=s1
          end if
          if ( i_4 == i_2 + 3 ) then
             i_1 = i_2
             s0=0.0D+00
             s1=0.0D+00
             s2=0.0D+00
             do j_1=j_2,j_3
                s0=s0+dconjg(v1(j_1,1))*z(j_1,i_1+0)
                s1=s1+dconjg(v1(j_1,1))*z(j_1,i_1+1)
                s2=s2+dconjg(v1(j_1,1))*z(j_1,i_1+2)
             end do! j_1
             ss1(i_1+0)=s0
             ss1(i_1+1)=s1
             ss1(i_1+2)=s2
          end if
!$        end if
!$OMP DO SCHEDULE(STATIC)
          do i_1=i_4,i_3,4
             s0=0.0D+00
             s1=0.0D+00
             s2=0.0D+00
             s3=0.0D+00
             do j_1=j_2,j_3
                s0=s0+dconjg(v1(j_1,1))*z(j_1,i_1+0)
                s1=s1+dconjg(v1(j_1,1))*z(j_1,i_1+1)
                s2=s2+dconjg(v1(j_1,1))*z(j_1,i_1+2)
                s3=s3+dconjg(v1(j_1,1))*z(j_1,i_1+3)
             end do! j_1
             ss1(i_1+0)=s0
             ss1(i_1+1)=s1
             ss1(i_1+2)=s2
             ss1(i_1+3)=s3
          end do! i_1
!$OMP END DO
*
!$OMP BARRIER
*
!$OMP MASTER
          ds = eigen_get_wtime()
          call reduce_cmplx(ss1(i_2),tt, i_3-i_2+1, 1, x_COMM_WORLD)
          de = eigen_get_wtime()
          dcom = dcom + (de - ds)

          do i_1=i_2,i_3
             ss1(i_1) = ss1(i_1) / dconjg(beta(i))
          end do! i_1
!$OMP END MASTER

!$OMP BARRIER

          j_2 = eigen_loop_start(1,      'X')
          j_3 = eigen_loop_end  (i-iblk, 'X')
          i_4=MOD(i_3-i_2+1,4)+i_2

!$OMP MASTER
          if ( i_4 == i_2 + 1 ) then
             i_1 = i_2
             s0 = ss1(i_1+0)
             do j_1=j_2,j_3
                z(j_1,i_1+0) = z(j_1,i_1+0) + s0 * v1(j_1,1)
             end do! j_1
          end if
          if ( i_4 == i_2 + 2 ) then
             i_1 = i_2
             s0 = ss1(i_1+0)
             s1 = ss1(i_1+1)
             do j_1=j_2,j_3
                z(j_1,i_1+0) = z(j_1,i_1+0) + s0 * v1(j_1,1)
                z(j_1,i_1+1) = z(j_1,i_1+1) + s1 * v1(j_1,1)
             end do! j_1
          end if
          if ( i_4 == i_2 + 3 ) then
             i_1 = i_2
             s0 = ss1(i_1+0)
             s1 = ss1(i_1+1)
             s2 = ss1(i_1+2)
             do j_1=j_2,j_3
                z(j_1,i_1+0) = z(j_1,i_1+0) + s0 * v1(j_1,1)
                z(j_1,i_1+1) = z(j_1,i_1+1) + s1 * v1(j_1,1)
                z(j_1,i_1+2) = z(j_1,i_1+2) + s2 * v1(j_1,1)
             end do! j_1
          end if
!$OMP END MASTER
!$OMP DO
          do i_1=i_4,i_3,4
             s0 = ss1(i_1+0)
             s1 = ss1(i_1+1)
             s2 = ss1(i_1+2)
             s3 = ss1(i_1+3)
             do j_1=j_2,j_3
                z(j_1,i_1+0) = z(j_1,i_1+0) + s0 * v1(j_1,1)
                z(j_1,i_1+1) = z(j_1,i_1+1) + s1 * v1(j_1,1)
                z(j_1,i_1+2) = z(j_1,i_1+2) + s2 * v1(j_1,1)
                z(j_1,i_1+3) = z(j_1,i_1+3) + s3 * v1(j_1,1)
             end do! j_1
          end do! i_1
!$OMP ENDDO
*
!$OMP BARRIER
*
       end do
*
*----------------------------------------------------------------
*
       d2 = eigen_get_wtime()
       d0 = d2 -d1

#define	USE_BCASTW	1

       mode = 0; i = nx+1
*
!$OMP MASTER
          if ( nx+1 <= n ) then
             ds = eigen_get_wtime()
             call hrbakwy_datacast(i_3, m, i, a(1,1), nma,
     $                          v1(1,1), nm, ss1(1),
#if !USE_BCASTW
     $                          wk(1),
#endif
     $                          iblk)
             de = eigen_get_wtime()
             dcom = dcom + (de - ds)
          end if
!$OMP END MASTER
!$OMP BARRIER
!$OMP MASTER
          if ( nx+1+m <= n ) then
             ds = eigen_get_wtime()
             call hrbakwy_datacast(i_3, m, i+m, a(1,1), nma,
     $                          v2(1,1), nm, ss2(1),
#if !USE_BCASTW
     $                          wk(1),
#endif
     $                          iblk)
             de = eigen_get_wtime()
             dcom = dcom + (de - ds)
          end if
!$OMP END MASTER
*
       do i=nx+1, n, m
*
          if ( mode == 0 ) then
             call hrbakwy_body2(i_3, z(1,1), nmz,
     $            beta(1), v1(1,1), nm, m,
     $            i, ss1(1), tt, nss, iblk,
     $            dcom, dx, dy, dz)
          endif
          if ( mode == 1 ) then
             call hrbakwy_body2(i_3, z(1,1), nmz,
     $            beta(1), v2(1,1), nm, m,
     $            i, ss2(1), tt, nss, iblk,
     $            dcom, dx, dy, dz)
          endif
          if ( mode == 2 ) then
             call hrbakwy_body2(i_3, z(1,1), nmz,
     $            beta(1), v3(1,1), nm, m,
     $            i, ss3(1), tt, nss, iblk,
     $            dcom, dx, dy, dz)
          endif
*
!$OMP MASTER
          if ( i+2*m <= n ) then
          ds = eigen_get_wtime()
          if ( mode == 0 ) then
             call hrbakwy_datacast(i_3, m, i+2*m, a(1,1), nma,
     $                             v3(1,1), nm, ss3(1),
#if !USE_BCASTW
     $                             wk(1),
#endif
     $                             iblk)
          end if
          if ( mode == 1 ) then
             call hrbakwy_datacast(i_3, m, i+2*m, a(1,1), nma,
     $                             v1(1,1), nm, ss1(1),
#if !USE_BCASTW
     $                             wk(1),
#endif
     $                             iblk)
          end if
          if ( mode == 2 ) then
             call hrbakwy_datacast(i_3, m, i+2*m, a(1,1), nma,
     $                             v2(1,1), nm, ss2(1),
#if !USE_BCASTW
     $                             wk(1),
#endif
     $                             iblk)
          end if
          de = eigen_get_wtime()
          dcom = dcom + (de - ds)
          end if
!$OMP END MASTER
*
          mode = MOD( mode + 1, 3 )
*
       end do
*
*----------------------------------------------------------------
*
!$     if ( omp_get_num_threads() > 1 ) then
!$OMP BARRIER
!$        if ( omp_get_thread_num() == 1 ) then
!$           ss1(1) = dx; ss1(2) = dy
!$        endif
!$OMP BARRIER
!$        if ( omp_get_thread_num() == 0 ) then
!$           dx = ss1(1); dy = ss1(2)
!$        endif
!$     endif
!$OMP BARRIER
*
#if AT_BCAST_OVERLAP
!===
!$OMP MASTER
!$     call omp_destroy_lock( TRBK_LOCK )
!$OMP END MASTER
!===
#endif
!$OMP MASTER
       call MPI_Barrier(TRD_COMM_WORLD,ierr)
!$OMP END MASTER
!$OMP BARRIER
!$OMP MASTER
       deallocate(wk)

#if TIMER_PRINT
       d2 = eigen_get_wtime()

       if ( 1 == TRD_inod ) then
          print*,"HRBAK=",(d2-d1)
          print*,"COMM=",dcom
          print*,"   ",d2-d1,(2d0*nvec*n*n)/(d2-d1)*1d-9,"GFLOPS"
          print*,"   ",dx,(1d0*nvec*n*n)/(dx)*1d-9,"GFLOPS"
          print*,"   ",dy,(1d0*nvec*n*n)/(dy)*1d-9,"GFLOPS"
          print*,"   ",d0,dz
       end if
#endif
!$OMP END MASTER
*
*
       return

       end subroutine

       subroutine hrbakwy_datacast(i_3, m, i, a, nma, v, nm,
     $                             ss,
#if !USE_BCASTW
     $                             wk,
#endif
     $                             iblk)
!$     use OMP_LIB
       use MPI
       use eigen_devel_mod
       use eigen_libs_mod
       use eigen_HRBAK_mod
       implicit NONE

       integer, intent(in)   :: i_3, m, i, nma, nm, iblk
       complex(8), intent(in)   :: a(nma,*)
       complex(8), intent(out)  :: ss(*)
#if !USE_BCASTW
       complex(8), intent(out)  :: wk(*)
#endif
       complex(8), intent(out)  :: v(nm,*)

       integer :: j, iy, k0
       integer :: i_1
       integer :: j_1, j_2, j_3, j_4, j_5
       integer :: nodes(0:nsm-1)

!$     real(8) :: ds, de
       real(8) :: du1,du2
       integer :: req(2), IERR

*
!$      ds = omp_get_wtime()

        ss(1:i_3*m+ns0) = 0.0D+0

#if USE_BCASTW
        j_3 = MIN(eigen_loop_end (i+m-1-iblk, 'X'), nm)
!-
      if ( MOD( m, y_nnod ) == 0 ) then
!-
        do j=0,m-1
           if ( y_inod == eigen_owner_node(i+j, 'Y') ) then
              iy  = eigen_owner_node(i+j, 'Y')
              i_1 = eigen_translate_g2l(i+j, y_nnod, iy)
              k0  = (j/y_nnod)*j_3
              v(1+k0:j_3+k0, 1) = a(1:j_3, i_1)
           end if
        end do

        k0  = (m/y_nnod)*j_3
        call allgather_cmplx(v, HRBK_buf, k0, y_COMM_WORLD)

        do j=0,m-1
           iy = eigen_owner_node(i+j, 'Y')
           k0 = ((iy-1)*(m/y_nnod)+(j/y_nnod))*j_3
           v(1:j_3, j+1) = HRBK_buf(1+k0:j_3+k0)
           j_2 = eigen_loop_start(i+j,  'X')
           if ( j_2 <= j_3 ) v(j_2:j_3, j+1) = ZERO
        end do
!-
      else
!-
        do j=0,m-1
           if ( y_inod == eigen_owner_node(i+j, 'Y') ) then
              iy  = eigen_owner_node(i+j, 'Y')
              i_1 = eigen_translate_g2l(i+j, y_nnod, iy)
              v(1:j_3, j+1) = a(1:j_3, i_1)
              j_2 = eigen_loop_start(i+j,  'X')
              if ( j_2 <= j_3 ) v(j_2:j_3, j+1) = ZERO
           end if
        end do
        do j=0,m-1,y_nnod
           iy  = eigen_owner_node(i+j, 'Y')
           k0  = MIN(m-j, y_nnod)
           call bcastw_cmplx(v(1,j+1), j_3, iy,
     &                     nm, k0, HRBK_buf, y_COMM_WORLD)
        end do
!-
      end if
#else
!-
        IF ( m > y_nnod .AND. y_nnod > 1 ) THEN

          do j=0,m-1
             nodes(j) = eigen_owner_node(i+j, 'Y')
          enddo

          do iy=1,y_nnod

             j_2 = eigen_loop_start(1,         'X')
             j_3 = eigen_loop_end  (i+m-1-iblk,'X')

             k0=0
             do j=0,m-1
             if ( nodes(j) == iy ) then
                i_1 = eigen_translate_g2l(i+j, y_nnod, iy)
                if ( y_inod == iy ) then
                   do j_1=j_2,j_3
                      wk(k0+j_1) = a(j_1, i_1)
                   end do! k
                endif
                k0=k0+(j_3-j_2+1)
             endif
             enddo

             call bcast_cmplx(wk, k0, iy, y_COMM_WORLD)

             k0=0
             do j=0,m-1
             if ( nodes(j) == iy ) then
                do j_1=j_2,j_3
                   v(j_1, j+1) = wk(k0+j_1)
                end do! k
                k0=k0+(j_3-j_2+1)
             endif
             enddo

             do j=0,m-1
             if ( nodes(j) == iy ) then
                j_4 = eigen_loop_start(i+j,       'X')
                j_5 = eigen_loop_end  (i+m-1-iblk,'X')
!OCL XFILL
                do j_1=j_4,j_5
                   v(j_1, j+1) = 0.0D+00
                end do
             endif
             enddo

          enddo

        ELSE

          do j=0,m-1
             nodes(j) = eigen_owner_node(i+j, 'Y')
             if ( nodes(j) == y_inod ) then
                i_1 = eigen_translate_g2l(i+j, 'Y')
                j_2 = eigen_loop_start(1,         'X')
                j_3 = eigen_loop_end  (i+m-1-iblk,'X')
                do j_1=j_2,j_3
                   v(j_1, j+1) = a(j_1, i_1)
                end do! k
                j_2 = eigen_loop_start(i+j,       'X')
                j_3 = eigen_loop_end  (i+m-1-iblk,'X')
!CDIR NOVECTOR
                do j_1=j_2,j_3
                   v(j_1, j+1) = 0.0D+00
                end do
             end if
          end do
          j_2 = eigen_loop_start(1,          'X')
          j_3 = eigen_loop_end  (i+m-1-iblk, 'X')
          do j=0,m-1
             call bcast_cmplx(v(1,j+1), j_3-j_2+1,
     &                      nodes(j), y_COMM_WORLD)
          enddo

        ENDIF
!-
#endif

!$      de = omp_get_wtime()
!$      TRBK_TIME_BCAST = de - ds
*
       return
       end subroutine
