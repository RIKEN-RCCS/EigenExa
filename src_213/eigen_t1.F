!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_t1_mod
!
! Purpose
! =======
!
! eigen_t1_mod manages the modules for the Rank-2k update routine.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module eigen_t1_mod

      use eigen_libs_mod
      use eigen_devel_mod
!$    use omp_lib

      implicit none
      private

      public :: eigen_common_2update

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_common_2update
!
! Purpose
! =======
!
! parallel Rank-2k update, A := A - Vy*U - Uy*V
!
!
! Arguments
! =========
!
! ar      (input/output) real(8) array, dimension(nma,*)
!         a contains the element of upper triangluar part of the
!         symmetric matrix.
!
! nma     (input) integer
!         The leading dimenstion of the array a.
!
! ur      (input) real(8) array, dimension(nmv,*)
!         ur contains m_size vectors of U.
!
! uyr     (input) real(8) array, dimension(nmv,*)
!         uyr contains m_size vectors of U^t.
!
! vr      (input) real(8) array, dimension(nmv,*)
!         vr contains m_size vectors of V.
!
! vyr     (input) real(8) array, dimension(nmv,*)
!         vyr contains m_size vectors of V^t.
!
! nmv     (input) integer
!         The leading dimenstion of the array U and V.
!
! m_size  (input) integer
!         The number of vectors bunched.
!
! i_base  (input) integer
!         The dimenion of the sub-matrix of A
!     
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine  eigen_common_2update(
     &     ar, nma,
     &     ur, uyr, vr, vyr, nmv,
     &     m_size, i_base)

      real(8), intent(inout) :: ar (nma, *)
      integer, intent(in)    :: nma
      real(8), intent(in)    :: ur (nmv, *)
      real(8), intent(in)    :: uyr(nmv, *)
      real(8), intent(in)    :: vr (nmv, *)
      real(8), intent(in)    :: vyr(nmv, *)
      integer, intent(in)    :: nmv
      integer, intent(in)    :: m_size
      integer, intent(in)    :: i_base

#if __FUJITSU
      integer                :: BLAS_CHUNK_N
      integer                :: BLAS_CHUNK_M
#else
      integer, parameter     :: BLAS_CHUNK_N = 8*4
      integer, parameter     :: BLAS_CHUNK_M = 96
#endif

      integer                :: k1, k2, m, n
      integer                :: i_1,i_2,i_3,i_4
      integer                :: j_1,j_2,j_3,j_4
      integer                :: i_loop_sta, i_loop_end
      integer                :: j_loop_sta, j_loop_end
      integer                :: blk_sizem, blk_sizen, blk_sizek
      integer                :: ii_step
      integer                :: local_rank, local_size

      intrinsic              :: min, max
      external               :: dgemm


      if (i_base <= 0) return

      local_rank = 0
      local_size = 1
!$    local_rank = omp_get_thread_num()
!$    local_size = omp_get_num_threads()

      m  = m_size
      n  = eigen_translate_g2l(i_base, 'Y')

!     12 + 8 + 4 = 24 = 8*3
#if __FUJITSU
! best parameters confrmed via preliminary benchmark
      BLAS_CHUNK_M = 96
      BLAS_CHUNK_N = 40
#endif

      ii_step = 0

!***************************************************
!*** LOOP TYPE 1
!***************************************************
!
! do j1= 1, i_base
!   do i1= j1, i_base
!
! do j_1= loop_start(1,X), loop_end(i_base,X)
!   g_1 = translate_l2g(j_1,X)
!   do i_1= loop_start(g_1_,Y), loop_end(i_base,Y)
!
! do j_1= loop_start(1,X), loop_end(i_base,X), CHUNK_M
!   g_1 = translate_l2g(j_1,X)
!   do i_1= loop_start(g_1_,Y), loop_end(i_base,Y), CHUNK_N
!     j_2=j_1
!     j_3=min(j_1+CHUNK_M-1,loop_end(i_base,X))
!     i_2=i_1
!     i_3=min(i_1+CHUNK_N-1,loop_end(i_base,Y))
!
!   /* basically matrix to be updated is inside (j_2:j_3,i_2:i_3) */
!      --[][][][]
!          [][][]
!            [][]
!              [] ==> (j_2:min(j_3,loop_end(l2g(i_3,Y),X)),
!              --      max(i_2,loop_start(l2g(j_2,X),Y)):i_3) */
!     j2 = j_2
!     j3 = min(j_3,loop_end(l2g(i_3,Y),X))
!     i2 = max(i_2,loop_start(l2g(j_2,X),Y)) =i_2
!     i3 = i_3
!
!     blk_sizem = j3-j2+1
!     blk_sizen = i3-i2+1
!     blk_sizek = m
!
!   enddo
! enddo

#if 0
      j_loop_sta = 1
      j_loop_end = eigen_loop_end  (i_base, 'X')

      do j_1 = j_loop_sta, j_loop_end, BLAS_CHUNK_M

        k1 = eigen_translate_l2g(j_1, 'X')

        i_loop_sta = eigen_loop_start(k1,     'Y')
        i_loop_end = eigen_loop_end  (i_base, 'Y')

        do i_1 = i_loop_sta, i_loop_end, BLAS_CHUNK_N

          j_4 = min(j_1+BLAS_CHUNK_M-1, j_loop_end)
          i_4 = min(i_1+BLAS_CHUNK_N-1, i_loop_end)

          k2  = eigen_translate_l2g(i_4, 'Y')

          j_2 = j_1
          j_3 = min(j_4, eigen_loop_end  (k2, 'X'))

          i_2 = i_1
          i_3 = i_4

          blk_sizem = j_3-j_2+1
          blk_sizen = i_3-i_2+1
          blk_sizek = m

          if (blk_sizem > 0 .and. blk_sizen > 0) then

            if (mod(ii_step, local_size) == local_rank) then

#if __FUJITSU
#define dgemm   dl_gemm
#endif

              call dgemm('N','T',
     &             blk_sizem, blk_sizen, blk_sizek,
     &             MONE, ur (j_2, 1),   nmv,
     &             vyr(i_2, 1),   nmv,
     &             ONE,  ar (j_2, i_2), nma)

              call dgemm('N','T',
     &             blk_sizem, blk_sizen, blk_sizek,
     &             MONE, vr (j_2, 1),   nmv,
     &             uyr(i_2, 1),   nmv,
     &             ONE,  ar (j_2, i_2), nma)

            end if

            ii_step = ii_step + 1

          end if

        end do                  ! i_1

      end do                    ! j_1
#endif

!***************************************************
!*** LOOP TYPE 2
!***************************************************
!
! do i1= 1, i_base
!   do j1= 1, i1
!
! do i_1= loop_start(1,Y), loop_end(i_base,Y)
!   g_1 = translate_l2g(i_1,Y)
!   do j_1= loop_start(1_,X), loop_end(g_1,X)
!
! do i_1= loop_start(1,Y), loop_end(i_base,Y), CHUNK_N
!   g_1 = translate_l2g(i_1+CHUNK_N-1,Y)
!   do j_1= loop_start(1_,X), loop_end(g_1,X), CHUNK_M
!     j_2=j_1
!     j_3=min(j_1+CHUNK_M-1,loop_end(g_1,X))
!     i_2=i_1
!     i_3=min(i_1+CHUNK_N-1,loop_end(i_base,Y))
!
! /*  (j_2:min(j_3,loop_end(l2g(i_3,Y),X)),
!     max(i_2,loop_start(l2g(j_2,X),Y)):i_3) */
!
!     j2 = j_2
!     j3 = min(j_3,loop_end(l2g(i_2,Y),X)) =j_3
!     i2 = max(i_2,loop_start(l2g(j_2,X),Y))
!     i3 = i_3
!
!     blk_sizem = j3-j2+1
!     blk_sizen = i3-i2+1
!     blk_sizek = m

#if 1
      i_loop_sta = 1
      i_loop_end = eigen_loop_end  (i_base, 'Y')

      do i_1 = i_loop_sta, i_loop_end, BLAS_CHUNK_N

        i_4 = min(i_1+BLAS_CHUNK_N-1, i_loop_end)
        k1  = eigen_translate_l2g(i_4, 'Y')

        j_loop_sta = 1
        j_loop_end = eigen_loop_end  (k1, 'X')

        do j_1 = j_loop_sta, j_loop_end, BLAS_CHUNK_M

          j_4 = min(j_1+BLAS_CHUNK_M-1, j_loop_end)
          k2  = eigen_translate_l2g(j_1, 'X')

          j_2 = j_1
          j_3 = j_4

          i_2 = max(i_1, eigen_loop_start  (k2, 'Y'))
          i_3 = i_4

          blk_sizem = j_3-j_2+1
          blk_sizen = i_3-i_2+1
          blk_sizek = m

          if (blk_sizem > 0 .and. blk_sizen > 0) then

            if (mod(ii_step, local_size) == local_rank) then

#if __FUJITSU
#define dgemm   dl_gemm
#endif

              call dgemm('N','T',
     &             blk_sizem, blk_sizen, blk_sizek,
     &             MONE, ur (j_2, 1),   nmv,
     &             vyr(i_2, 1),   nmv,
     &             ONE,  ar (j_2, i_2), nma)

              call dgemm('N','T',
     &             blk_sizem, blk_sizen, blk_sizek,
     &             MONE, vr (j_2, 1),   nmv,
     &             uyr(i_2, 1),   nmv,
     &             ONE,  ar (j_2, i_2), nma)

            end if

            ii_step = ii_step + 1

          end if

        end do                  ! j_1

      end do                    ! i_1
#endif

      return
      end subroutine  eigen_common_2update

      end module eigen_t1_mod
