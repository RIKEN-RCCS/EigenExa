
!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: dc_redist1_r
!
! Purpose
! =======
!
! dc_redist1 is a utility function, which redistributes the data-format
! from 2D-cyclic(1) to 2D-cyclic(NB).
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
! a       (output) real(8) array, dimension(ldm,*)
!         recv buffer
!
! b       (input) real(8) array, dimension(ldm,*)
!         send buffer
!
! ldm     (input) integer
!         leading dimensions
!
! wk      (input/output) real(8) array, dimension(lwk)
!         working buffer
!
! lwk     (input) integer
!         The length of the working buffer
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine dc_redist1_r(n, NB, a, b, ldm, wk, lwk)

      use eigen_devel_mod
      use mpi
!$    use omp_lib

      implicit none

      integer, intent(in)    :: n
      integer, intent(in)    :: NB
      real(8), intent(out)   :: a(ldm,*)
      real(8), intent(in)    :: b(ldm,*)
      integer, intent(in)    :: ldm
      real(8), intent(inout) :: wk(1:lwk)
      integer, intent(in)    :: lwk

      real(8), allocatable   :: wk1(:, :), wk2(:, :)

      integer                :: i,j,j0,j1,k,k0,l,lx,IERR
      integer                :: iNQ
      integer                :: iblk_, jblk_
      integer                :: iblk, jblk
      integer                :: idist, ir_size, is_size
      integer                :: his_rank, her_rank
      integer                :: temp(2)

      integer                :: local_size, local_rank, th_counter
      integer                :: NQ, NBQ
      integer                :: j_end, j_begin


      iblk_ = (n-1)/y_nnod+1
      iblk  = (iblk_-1)/NB+1
      jblk_ = (n-1)/x_nnod+1
      jblk  = (jblk_-1)/NB+1

!     
!     Calculate the size of required buffer for re-distribution;
!     divide the Block width by NQ and get the buffer shape
!     as (1:ldm, lx) where lx is the max number of vectors to be
!     transfered such that ldm * lx * 2 <= lwk (lwk is the size of
!     working vector passed from the callee procedure.

      NQ = 0
      lx = 0

      if (y_inod == 1) then
        do iNQ=1,NB

          NQ = iNQ
          NBQ = (NB-1)/NQ+1

!     Emulate the communication and evaluate data size
          lx = 0
!=#if !defined(__INTEL_COMPILER)
!$OMP PARALLEL DO
!$OMP+         PRIVATE(i,j0,j1,idist,his_rank,her_rank,j_end,j_begin,l)
!$OMP+         REDUCTION(MAX:lx)
!=#endif
          do i=1,iblk
            do j0=1,NB,NBQ; j1 = min(NB,j0+NBQ-1)
            do idist = 1, y_nnod-1

!     stride communication with the upper and the lower processess

              his_rank = mod(y_inod-1-idist+y_nnod,y_nnod)+1
              her_rank = mod(y_inod-1+idist+y_nnod,y_nnod)+1

!     the number of vectors which should be sent to her_rank
              call dc1_compute_loop_iterations_x(n, i, j0, j1, NB,
     &             her_rank, y_inod, y_nnod, j_end, j_begin, l)
              lx = max(l, lx)

!     the number of vectors which should be recieved from his_rank
              call dc1_compute_loop_iterations_x(n, i, j0, j1, NB,
     &                 y_inod, his_rank, y_nnod, j_end, j_begin, l)
              lx = max(l, lx)

            end do
            end do
          end do

          if (ldm*lx*2 <= lwk) then
            exit
          end if

        end do
      end if

      temp(1) = NQ
      temp(2) = lx
      call MPI_Bcast(temp(1), 2, MPI_INTEGER, 0,
     &     y_COMM_WORLD, IERR)
      NQ  = temp(1)
      NBQ = (NB-1)/NQ+1
      lx  = temp(2)


      if (ldm*lx*2 <= lwk) then
!     The most eco case
        call dc_redist1_sub_x(n, NB, a, b, ldm,
     &        wk(1), wk(1+lx*ldm),
     &        NBQ)
      else if (ldm*lx <= lwk) then
!     One buffer is allocated
        allocate(wk2(1:ldm, 1:lx))
        call dc_redist1_sub_x(n, NB, a, b, ldm,
     &        wk(1), wk2,
     &        NBQ)
        deallocate(wk2)
      else
!     Both buffers are allocated
        allocate(wk1(1:ldm, 1:lx))
        allocate(wk2(1:ldm, 1:lx))
        call dc_redist1_sub_x(n, NB, a, b, ldm,
     &       wk1, wk2,
     &       NBQ)
        deallocate(wk1)
        deallocate(wk2)
      end if

      return
      contains

!------------------------------------------------------------------
!------------------------------------------------------------------

      subroutine dc_redist1_sub_x(n, NB, a, b, ldm, wk1, wk2, NBQ)

      use comm_mod
      use eigen_devel_mod
      use mpi
!$    use omp_lib

      implicit none

      integer, intent(in)    :: n
      integer, intent(in)    :: NB
      real(8), intent(out)   :: a(ldm,*)
      real(8), intent(in)    :: b(ldm,*)
      integer, intent(in)    :: ldm
      real(8), intent(inout) :: wk1(ldm,*)
      real(8), intent(inout) :: wk2(ldm,*)
      integer, intent(in)    :: NBQ

      integer                :: i,j,j0,j1,k,k0,k1,l,lx
      integer                :: iblk_, jblk_
      integer                :: iblk,  jblk
      integer                :: idist, ir_size, is_size
      integer                :: iq_r, iq_s, his_rank, her_rank
      integer                :: j_begin1, j_end1, lx1
      integer                :: j_begin2, j_end2, lx2


      iblk_ = (n-1)/y_nnod+1
      iblk  = (iblk_-1)/NB+1
      jblk_ = (n-1)/x_nnod+1
      jblk  = (jblk_-1)/NB+1

      do i=1,iblk
        do j0=1,NB,NBQ; j1 = min(NB,j0+NBQ-1)

          call dc1_compute_loop_iterations_x(n, i, j0, j1, NB,
     &         y_inod, y_inod, y_nnod, j_end1, j_begin1, lx)

          K=j_begin1+((y_inod-1)+(i-1)*y_nnod)*NB
          K0=(K-1)/y_nnod+1-1
          K1=j_begin1+NB*(i-1)
!$OMP PARALLEL DO
!$OMP+         PRIVATE(j,L)
          do j=1,lx
            L=K1+(j-1)*y_nnod
            a(1:ldm, L) = b(1:ldm, K0+j) 
          end do
!$OMP END PARALLEL DO

          do idist = 1, y_nnod-1

            her_rank = mod(y_inod-1+idist+y_nnod,y_nnod)+1
            his_rank = mod(y_inod-1-idist+y_nnod,y_nnod)+1

            call dc1_compute_loop_iterations_x(n, i, j0, j1, NB,
     &           y_inod, his_rank, y_nnod, j_end1, j_begin1, lx1)
            ir_size=ldm*lx1

            call dc1_compute_loop_iterations_x(n, i, j0, j1, NB,
     &           her_rank, y_inod, y_nnod, j_end2, j_begin2, lx2)
            is_size=ldm*lx2


            if (is_size > 0) then
              call irecv_dbl(wk1, is_size, her_rank, iq_s,
     &             y_COMM_WORLD)
            end if

            if (ir_size > 0) then
              K=j_begin1+((his_rank-1)+(i-1)*y_nnod)*NB
              K0=(K-1)/y_nnod+1-1
!$OMP PARALLEL DO
!$OMP+         PRIVATE(j)
              do j=1,lx1
                wk2(1:ldm, j) = b(1:ldm, K0+j) 
              end do
!$OMP END PARALLEL DO

              call isend_dbl(wk2, ir_size, his_rank, iq_r,
     &             y_COMM_WORLD)
            end if


            if (is_size > 0) then
              call wait_dbl(iq_s)

              K1=j_begin2+NB*(i-1)
!$OMP PARALLEL DO
!$OMP+         PRIVATE(j,L)
              do j=1,lx2
                L=K1+(j-1)*y_nnod
                a(1:ldm, L) = wk1(1:ldm, j) 
              end do
!$OMP END PARALLEL DO
            end if

            if (ir_size > 0) then
              call wait_dbl(iq_r)
            end if

          end do

        end do

      end do

      return
      end subroutine  dc_redist1_sub_x

!------------------------------------------------------------------
!------------------------------------------------------------------

      subroutine dc1_compute_loop_iterations_x(n, i, j0, j1, NB,
     &     id_dest, id_from, y_nnod, j_end, j_begin, lx)

      implicit none

      integer, intent(in)  :: n
      integer, intent(in)  :: i
      integer, intent(in)  :: j0
      integer, intent(in)  :: j1
      integer, intent(in)  :: NB
      integer, intent(in)  :: id_dest
      integer, intent(in)  :: id_from
      integer, intent(in)  :: y_nnod
      integer, intent(out) :: j_end
      integer, intent(out) :: j_begin
      integer, intent(out) :: lx

      integer              :: kc
      integer              :: kn, km, kx

!     
!  This routine returns the pair of loop_begin and loop_end,
!  also the loop count, when the follwing first loop structure
!  is translated to the second loop-structure
!        do j=j0,j1
!           K=j+((id_from-1)+(i-1)*y_nnod)*NB
!           if (K <= n .and. mod(K-1,y_nnod)+1 == id_dest) then
!              ......
!           end if
!        end do
! ==>
!        do j=j_begin,j_end,y_nnod
!           K=j+((id_from-1)+(i-1)*y_nnod)*NB
!           ......
!        end do
!     

      kx = (id_from-1)*NB
      kn = kx+(i-1)*y_nnod*NB
      km = mod(j0+kx-1, y_nnod)

      j_begin = j0 + mod((id_dest-1)+y_nnod-km, y_nnod)
      j_end   = j1 ! min(n - kn, j1)

!     
! lx returns the number of iterations
!     
      kc = 0
      if (j_begin <= j_end) then
        kc = (j_end - j_begin + y_nnod)/y_nnod
      end if
      lx = kc

      return
      end subroutine dc1_compute_loop_iterations_x

      end subroutine  dc_redist1_r
