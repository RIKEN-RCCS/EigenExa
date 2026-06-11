
!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: dc_redist2
!
! Purpose
! =======
!
! dc_redist2 is a utility function, which redistributes the data-format
! from 2D-cyclic(NB) to 2D-cyclic(1).
!
!
! Arguments
! =========
!
! n       (input) integer
!         Dimension of the vectors
!
! nvec    (input) integer
!         The number of vectors to be transformed
!
! NB      (input) integer
!         Block factor
!
! a       (input) real(8) array, dimension(lda,*)
!         send buffer
!
! lda     (input) integer
!         The leading dimenstion of the array a
!
! b       (output) real(8) array, dimension(ldb,*)
!         recv buffer
!
! ldb     (input) integer
!         The leading dimenstion of the array b
!
! wk      (input/output) real(8) array, dimension(lwk)
!         working buffer
!
! lwk     (input) integer
!         The length of the working buffer
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine dc_redist2(n, nvec, NB, a, lda, b, ldb, wk, lwk)

      use eigen_devel_mod
      use mpi

      implicit none

      integer, intent(in)    :: n
      integer, intent(in)    :: nvec
      integer, intent(in)    :: NB
      real(8), intent(in)    :: a(lda,*)
      integer, intent(in)    :: lda
      real(8), intent(out)   :: b(ldb,*)
      integer, intent(in)    :: ldb
      real(8), intent(inout) :: wk(1:lwk)
      integer, intent(in)    :: lwk

      real(8), allocatable   :: wk1(:), wk2(:)
      integer, allocatable   :: ir_sz(:,:)
      integer, allocatable   :: l_(:)

      integer                :: i,i0,j,k,k0,l,lx,m,m0, IERR
      integer                :: iblk_, jblk_
      integer                :: iblk,  jblk
      integer                :: idist, ir_size, is_size
      integer                :: his_rank, her_rank

      integer                :: PACK, MPACK, temp(4)

      real(8)                :: d1, d2, aa(10)


      aa = 0

      iblk_ = (nvec-1)/y_nnod+1
      iblk  = (iblk_-1)/NB+1
      jblk_ = (n-1)/x_nnod+1
      jblk  = (jblk_-1)/NB+1


      d1 = eigen_get_wtime()
      allocate(ir_sz(x_nnod,jblk))
      d2 = eigen_get_wtime()
      aa(1) = aa(1) + (d2-d1)

      d1=eigen_get_wtime()

      PACK  = 0
      MPACK = 0
      lx    = 0

      if (TRD_inod == 1) then
        do m = 4,1,-1
          lx = 0
          PACK = m
          do i0=1,jblk,PACK
            do idist = 1, x_nnod-1

              her_rank = mod(x_inod-1+idist+x_nnod,x_nnod)+1
              his_rank = mod(x_inod-1-idist+x_nnod,x_nnod)+1

              l = 0
              do i=i0,min(jblk,i0+PACK-1)
                do j=1,NB
                  K=j+((x_inod-1)+(i-1)*x_nnod)*NB
                  if (K <= n .and. mod(K-1,x_nnod)+1 == her_rank) then
                    l = l + 1
                  end if
                end do
              end do
              lx = max(l, lx)

              l = 0
              do i=i0,min(jblk,i0+PACK-1)
                do j=1,NB
                  K=j+((his_rank-1)+(i-1)*x_nnod)*NB
                  if (K <= n .and. mod(K-1,x_nnod)+1 == x_inod) then
                    l = l + 1
                  end if
                end do
              end do
              lx = max(l, lx)

            end do
          end do
          lx = max(1,lx)

          MPACK = max(1,lwk/(2*lx))
          MPACK = min(MPACK, iblk_)
          if (2*lx*MPACK <= lwk) exit
        end do
      end if

      temp(1) = PACK
      temp(2) = MPACK
      temp(3) = lx
      call MPI_Bcast(temp, 3, MPI_INTEGER, 0,
     &     TRD_COMM_WORLD, IERR)
      PACK  = temp(1)
      MPACK = temp(2)
      lx    = temp(3)


      do i0=1,jblk,PACK
        do idist = 1, x_nnod-1
          his_rank = mod(x_inod-1-idist+x_nnod,x_nnod)+1
          l = 0
          do i=i0,min(jblk,i0+PACK-1)
            do j=1,NB
              K=j+((his_rank-1)+(i-1)*x_nnod)*NB
              if (K <= n .and. mod(K-1,x_nnod)+1 == x_inod) then
                l = l + 1
              end if
            end do
          end do
          ir_sz(idist,i0) = l
        end do
      end do


      d2 = eigen_get_wtime()
      aa(3) = aa(3) + (d2-d1)


      MPACK = max(1,lwk/(2*lx))
      MPACK = min(MPACK, iblk_)
      if (2*lx*MPACK <= lwk) then
        call dc_redist2_sub(n, nvec, NB, a, lda, b, ldb,
     &       wk(1), wk(1+lx*MPACK),
     &       ir_sz, x_nnod, PACK, MPACK, aa)
      else
        MPACK = max(1,lwk/(1*lx))
        MPACK = min(MPACK, iblk_)
        if (lx*MPACK <= lwk) then
          d1 = eigen_get_wtime()
          allocate(wk1(1:lx*MPACK))
          allocate(wk2(1:lx*MPACK))
          d2 = eigen_get_wtime()
          aa(1) = aa(1) + (d2-d1)
          call dc_redist2_sub(n, nvec, NB, a, lda, b, ldb,
     &         wk, wk2,
     &         ir_sz, x_nnod, PACK, MPACK, aa)
          d1 = eigen_get_wtime()
          deallocate(wk1)
          deallocate(wk2)
          d2 = eigen_get_wtime()
          aa(1) = aa(1) + (d2-d1)
        else
          MPACK = max(1,lwk/(2*lx))
          MPACK = min(MPACK, iblk_)
          d1 = eigen_get_wtime()
          allocate(wk1(1:lx*MPACK))
          allocate(wk2(1:lx*MPACK))
          d2 = eigen_get_wtime()
          aa(1) = aa(1) + (d2-d1)
          call dc_redist2_sub(n, nvec, NB, a, lda, b, ldb,
     &         wk1, wk2,
     &         ir_sz, x_nnod, PACK, MPACK, aa)
          d1 = eigen_get_wtime()
          deallocate(wk1)
          deallocate(wk2)
          d2 = eigen_get_wtime()
          aa(1) = aa(1) + (d2-d1)
        end if
      end if
      
      d1 = eigen_get_wtime()
      deallocate(ir_sz)
      d2 = eigen_get_wtime()
      aa(1) = aa(1) + (d2-d1)

!     if (TRD_inod == 1) print *, PACK, MPACK, aa(1:6)

      return
      end subroutine dc_redist2

!------------------------------------------------------------------
!------------------------------------------------------------------

      subroutine dc_redist2_sub(n, nvec, NB, a, lda, b, ldb,
     &     wk1, wk2, ir_sz, ldr, PACK, MPACK, aa)

      use comm_mod
      use eigen_devel_mod
      use mpi
!$    use omp_lib

      implicit none

      integer, intent(in)    :: n
      integer, intent(in)    :: nvec
      integer, intent(in)    :: NB
      real(8), intent(in)    :: a(lda,*)
      integer, intent(in)    :: lda
      real(8), intent(out)   :: b(ldb,*)
      integer, intent(in)    :: ldb
      real(8), intent(inout) :: wk1(*)
      real(8), intent(inout) :: wk2(*)
      integer, intent(in)    :: ir_sz(ldr,*)
      integer, intent(in)    :: ldr
      integer, intent(in)    :: PACK
      integer, intent(in)    :: MPACK
      real(8), intent(inout) :: aa(*)

      integer, allocatable   :: l_(:)

      integer                :: i, i0, i1, i2, j, j0, j1, j2
      integer                :: k, k0, k1, k2, l, L0, lx, m, m0
      integer                :: iblk_, jblk_
      integer                :: iblk,  jblk
      integer                :: idist, ir_size, is_size
      integer                :: iq(2), his_rank, her_rank

      real(8)                :: d1, d2


      iblk_ = (nvec-1)/y_nnod+1
      iblk  = (iblk_-1)/NB+1
      jblk_ = (n-1)/x_nnod+1
      jblk  = (jblk_-1)/NB+1


      allocate(l_(1:iblk_+1))

      do i0=1,jblk,PACK

        K1 = mod((x_inod-1)*NB,x_nnod)
        j2 = mod(x_inod-K1+x_nnod-1,x_nnod)+1
        i1 = (n-j2-(x_inod-1)*NB)/(x_nnod*NB)+1

        d1 = eigen_get_wtime()
!$OMP PARALLEL DO
!$OMP+         PRIVATE(i,j,j0,K,K0,L0,m)
        do m=1,iblk_
#if 0
          do i=i0,min(jblk,i0+PACK-1)
            do j=1,NB
              K=j+((x_inod-1)+(i-1)*x_nnod)*NB
              if (K <= n .and. mod(K-1,x_nnod)+1 == x_inod) then
                K0=(K-1)/x_nnod+1
                b(K0,m) = a(j+NB*(i-1),m)
              end if
            end do
          end do
#else
#if 0
          do i=i0,min(jblk,i0+PACK-1)
            K=j2+((x_inod-1)+(i-1)*x_nnod)*NB
            if (K <= n) then
              K0=(K-1)/x_nnod+1-1
              do j=j2,NB,x_nnod
                K  = j+((x_inod-1)+(i-1)*x_nnod)*NB
                if (K > n) exit
                K0 = K0 + 1
                b(K0,m) = a(j+NB*(i-1),m)
              end do
            end if
          end do
#else
          j0 = (NB-j2+x_nnod)/x_nnod
          do i=i0,min(i1,jblk,i0+PACK-1)
            K  = j2+((x_inod-1)+(i-1)*x_nnod)*NB
            K0 = (n-K+x_nnod)/x_nnod
            L0 = min(j0,K0)
            K0 = (K-1)/x_nnod+1-1
            do j=1,L0
              b(j+K0,m) = a((j-1)*x_nnod+j2+NB*(i-1),m)
            end do
          end do
#endif
#endif
        end do
!$OMP END PARALLEL DO
        d2 = eigen_get_wtime()
        aa(2) = aa(2) + (d2-d1)

        do idist = 1, x_nnod-1

          her_rank = mod(x_inod-1+idist+x_nnod,x_nnod)+1
          his_rank = mod(x_inod-1-idist+x_nnod,x_nnod)+1

          K1 = mod((x_inod-1)*NB,x_nnod)
          j1 = mod(her_rank-K1+x_nnod-1,x_nnod)+1
          i1 = (n-j1-(x_inod-1)*NB)/(x_nnod*NB)+1

          K2 = mod((his_rank-1)*NB,x_nnod)
          j2 = mod(x_inod-K2+x_nnod-1,x_nnod)+1
          i2 = (n-j2-(his_rank-1)*NB)/(x_nnod*NB)+1

          do m0=1,iblk_,MPACK

            ir_size=ir_sz(idist,i0)*min(iblk_-m0+1,MPACK)
            if (ir_size > 0) then

              d1 = eigen_get_wtime()
              call irecv_dbl(wk2, ir_size, his_rank, iq(2),
     &             x_COMM_WORLD)
              d2 = eigen_get_wtime()
              aa(4) = aa(4) + (d2-d1)

            end if


            l_(1)=0
!$OMP PARALLEL DO
!$OMP+         PRIVATE(i,j,j0,K,k0,l,m)
            do m=m0,min(iblk_,m0+MPACK-1)
              l = 0
              do i=i0,min(i1,jblk,i0+PACK-1)
                K  = j1+((x_inod-1)+(i-1)*x_nnod)*NB
                j0 = (NB-j1+x_nnod)/x_nnod
                k0 = (n-K+x_nnod)/x_nnod
                l  = l + min(j0,k0)
              end do
              l_(m-m0+1+1)=l
            end do
!$OMP END PARALLEL DO
            l=0
            do m=m0+1,min(iblk_,m0+MPACK-1)+1
              l=l_(m-m0+1)+l_(m-m0+1-1)
              l_(m-m0+1)=l
            end do


            is_size=l
            if (is_size > 0) then

              d1 = eigen_get_wtime()
!$OMP PARALLEL DO
!$OMP+         PRIVATE(i,j,j0,K,K0,l,L0,m)
              do m=m0,min(iblk_,m0+MPACK-1)
                l=l_(m-m0+1)
#if 0
                do i=i0,min(jblk,i0+PACK-1)
                  j0=NB+1
                  do j=1,NB
                    K=j+((x_inod-1)+(i-1)*x_nnod)*NB
                    if (K <= n .and. mod(K-1,x_nnod)+1 == her_rank)then
                      j0 = j
                      exit
                    end if
                  end do
                  do j=j0,NB,x_nnod
                    K=j+((x_inod-1)+(i-1)*x_nnod)*NB
                    if (K <= n) then
                      l = l + 1
                      wk1(l) = a(j+NB*(i-1),m)
                    end if
                  end do
                end do
#else
#if 0
                do i=i0,min(jblk,i0+PACK-1)
                  K  = j1+((x_inod-1)+(i-1)*x_nnod)*NB
                  if (K <= n) then
                    do j=j1,NB,x_nnod
                      l = l + 1
                      wk1(l) = a(j+NB*(i-1),m)
                      if (K > n) exit
                      K  = K + x_nnod
                    end do
                  end if
                end do
#else
                j0 = (NB-j1+x_nnod)/x_nnod
                do i=i0,min(i1,jblk,i0+PACK-1)
                  K  = j1+((x_inod-1)+(i-1)*x_nnod)*NB
                  K0 = (n-K+x_nnod)/x_nnod
                  L0 = min(j0,K0)
                  do j=1,L0
                    wk1(j+l) = a((j-1)*x_nnod+j1+NB*(i-1),m)
                  end do
                  l = l + L0
                end do
#endif
#endif
              end do
!$OMP END PARALLEL DO
              d2 = eigen_get_wtime()
              aa(5) = aa(5) + (d2-d1)

              d1 = eigen_get_wtime()
              call isend_dbl(wk1, is_size, her_rank, iq(1),
     &             x_COMM_WORLD)
              d2 = eigen_get_wtime()
              aa(4) = aa(4) + (d2-d1)

            end if


            ir_size=ir_sz(idist,i0)*min(iblk_-m0+1,MPACK)
            if (ir_size > 0) then

              l_(1)=0
!$OMP PARALLEL DO
!$OMP+         PRIVATE(i,j,j0,K,K0,l,m)
              do m=m0,min(iblk_,m0+MPACK-1)
                l = 0
                do i=i0,min(i2,jblk,i0+PACK-1)
                  K  = j2+((his_rank-1)+(i-1)*x_nnod)*NB
                  j0 = (NB-j2+x_nnod)/x_nnod
                  K0 = (n-K+x_nnod)/x_nnod
                  l  = l + min(j0,K0)
                end do
                l_(m-m0+1+1)=l
              end do
!$OMP END PARALLEL DO
              do m=m0+1,min(iblk_,m0+MPACK-1)+1
                l=l_(m-m0+1)+l_(m-m0+1-1)
                l_(m-m0+1)=l
              end do

              d1 = eigen_get_wtime()
              call wait_dbl(iq(2))
              d2 = eigen_get_wtime()
              aa(4) = aa(4) + (d2-d1)

              d1 = eigen_get_wtime()
!$OMP PARALLEL DO
!$OMP+         PRIVATE(i,j,j0,K,K0,l,L0,m)
              do m=m0,min(iblk_,m0+MPACK-1)
                l=l_(m-m0+1)
#if 0
                do i=i0,min(jblk,i0+PACK-1)
                  j0=NB+1
                  do j=1,NB
                    K=j+((his_rank-1)+(i-1)*x_nnod)*NB
                    if (K <= n .and. mod(K-1,x_nnod)+1 == x_inod) then
                      j0 = j
                      exit
                    end if
                  end do
                  K=j0+((his_rank-1)+(i-1)*x_nnod)*NB
                  K0=(K-1)/x_nnod+1-1
                  do j=j0,NB,x_nnod
                    K=j+((his_rank-1)+(i-1)*x_nnod)*NB
                    if (K <= n) then
                      l = l + 1
                      K0 = K0 + 1
                      b(K0,m) = wk2(l)
                    end if
                  end do
                end do

#else
#if 0
                do i=i0,min(jblk,i0+PACK-1)
                  K  = j2+((his_rank-1)+(i-1)*x_nnod)*NB
                  if (K <= n) then
                    K0=(K-1)/x_nnod+1-1
                    do j=j2,NB,x_nnod
                      K  = j+((his_rank-1)+(i-1)*x_nnod)*NB
                      if (K > n) exit
                      l = l + 1
                      K0 = K0 + 1
                      b(K0,m) = wk2(l)
                    end do
                  end if
                end do
#else
                j0 = (NB-j2+x_nnod)/x_nnod
                do i=i0,min(i2,jblk,i0+PACK-1)
                  K  = j2+((his_rank-1)+(i-1)*x_nnod)*NB
                  K0 = (n-K+x_nnod)/x_nnod
                  L0 = min(j0,K0)
                  K0 = (K-1)/x_nnod+1-1
                  do j=1,L0
                    b(j+K0,m) = wk2(j+l)
                  end do
                  l = l + L0
                end do
#endif
#endif
              end do
!$OMP END PARALLEL DO
              d2 = eigen_get_wtime()
              aa(6) = aa(6) + (d2-d1)

            end if

            if (is_size > 0) then

              d1 = eigen_get_wtime()
              call wait_dbl(iq(1))
              d2 = eigen_get_wtime()
              aa(4) = aa(4) + (d2-d1)

            end if

          end do

        end do

      end do

      d1 = eigen_get_wtime()
      deallocate(l_)
      d2 = eigen_get_wtime()
      aa(1) = aa(1) + (d2-d1)

      return
      end subroutine dc_redist2_sub
