!>
!> @file   FS_PDLAED1.F90
!> @brief  subroutine FS_PDLAED1
!!>
!
!
!
!> subroutine FS_PDLAED1
!>
!> @brief @n
!>  Purpose @n
!>  ======= @n
!>  FS_PDLAED1 computes the updated eigensystem of a diagonal           @n
!>  matrix after modification by a rank-one symmetric matrix,           @n
!>  in parallel.                                                        @n
!>                                                                      @n
!>  T = Q(in) ( D(in) + RHO * Z*Z' ) Q'(in) = Q(out) * D(out) * Q'(out) @n
!>                                                                      @n
!>  where Z = Q'u, u is a vector of length N with ones in the           @n
!>  N1 and N1 + 1 th elements and zeros elsewhere.                      @n
!>                                                                      @n
!>  The eigenvectors of the original matrix are stored in Q, and the    @n
!>  eigenvalues are in D.  The algorithm consists of three stages:      @n
!>                                                                      @n
!>  The first stage consists of deflating the size of the problem       @n
!>  when there are multiple eigenvalues or if there is a zero in        @n
!>  the Z vector.  For each such occurence the dimension of the         @n
!>  secular equation problem is reduced by one.  This stage is          @n
!>  performed by the routine FS_PDLAED2.                                @n
!>                                                                      @n
!>  The second stage consists of calculating the updated                @n
!>  eigenvalues. This is done by finding the roots of the secular       @n
!>  equation via the routine DLAED4 (as called by FS_PDLAED3).          @n
!>  This routine also calculates the eigenvectors of the current        @n
!>  problem.                                                            @n
!>                                                                      @n
!>  The final stage consists of computing the updated eigenvectors      @n
!>  directly using the updated eigenvalues.  The eigenvectors for       @n
!>  the current problem are multiplied with the eigenvectors from       @n
!>  the overall problem by the routine FS_PDLAED3.
!
!  Arguments
!  =========
!>
!> @param[in]     N        (global input) INTEGER @n
!>                         The order of the tridiagonal matrix T.  N >= 0.
!>
!> @param[in]     N1       (input) INTEGER @n
!>                         The location of the last eigenvalue in the leading sub-matrix. @n
!>                         min(1,N) <= N1 <= N.
!>
!> @param[in,out] D        (global input/output) DOUBLE PRECISION array, dimension (N) @n
!>                         On entry, the diagonal elements of the tridiagonal matrix.  @n
!>                         On exit, if INFO = 0, the eigenvalues in descending order.
!>
!> @param[in,out] Q        (local output) DOUBLE PRECISION array,                    @n
!>                         global dimension (N, N),                                  @n
!>                         local dimension (LDQ, NQ)                                 @n
!>                         Q  contains the orthonormal eigenvectors of the symmetric @n
!>                         tridiagonal matrix.
!>
!> @param[in]     LDQ      (local input) INTEGER @n
!>                         The leading dimension of the array Q.  LDQ >= max(1,NP).
!>
!> @param[in]     SUBTREE  (input) type(bt_node) @n
!>                         sub-tree information of merge block.
!>
!> @param[in,out] RHO      (input) DOUBLE PRECISION @n
!>                         The subdiagonal entry used to create the rank-1 modification.
!>
!> @param         WORK     (local workspace/output) DOUBLE PRECISION array, dimension (*)
!>
!> @param         IWORK    (local workspace/output) INTEGER array, dimension (*)
!>
!> @param[out]    INFO     (global output) INTEGER @n
!>                         = 0: successful exit   @n
!>                         /=0: error exit
!>
!> @param[out]    prof     (global output) type(FS_prof) @n
!>                         profiling information of each subroutines.
!>
!> @note This routine is modified from ScaLAPACK PDLAED1.f
!>
subroutine FS_PDLAED1(N, N1, D, Q, LDQ, SUBTREE, RHO, WORK, IWORK, INFO, prof)
      use FS_const_mod
      use FS_libs_mod
      use FS_prof_mod
      use FS_dividing_mod
      implicit none
!
!     .. Scalar Arguments ..
      integer          , intent(in)    :: N, N1, LDQ
      type(bt_node)    , intent(in)    :: SUBTREE
      real(kind(0.0d0)), intent(inout) :: RHO
      integer          , intent(out)   :: INFO
      type(FS_prof)    , intent(inout) :: prof
!
!     .. Array Arguments ..
      real(kind(0.0d0)), intent(inout) :: D(N)
      real(kind(0.0d0)), intent(inout) :: Q(LDQ,*)
      ! work
      real(kind(0.0d0)) :: WORK(*)
      integer           :: IWORK(*)
!     ..
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
      integer :: NBLK, NB, NP, NQ
      integer :: NPROW, NPCOL, MYROW, MYCOL
      integer :: LDQ2, LDU, LSENDQ2
      integer :: IZ, IDLAMDA, IW, IPQ2, IPU, IBUF, ISENDQ2, IRECVQ2
      integer :: ICTOT, IPSM, INDX, INDXC, INDXP, INDCOL, COLTYP, &
                 INDROW, INDXR, INDXCB
      integer :: IZWORK, IDWORK, IEND, IIEND
      integer :: I, K
      integer :: SC
      integer :: ierr
!     ..
!     .. Local Arrays ..
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
      INFO = 0
!
#ifdef _DEBUGLOG
!      CALL MPI_Barrier(SUBTREE%MERGE_COMM,ierr)
      if( FS_MYRANK.eq.0 ) then
        write(*,'(a)') "FS_PDLAED1 start."
      endif
#endif
#if TIMER_PRINT
      CALL FS_prof_start(prof, 30)
#endif
!
!     Quick return if possible
      IF( N.EQ.0 ) THEN
        GO TO 20
      END IF
!
!     The following values are  integer pointers which indicate
!     the portion of the workspace used by a particular array
!     in FS_PDLAED2 and FS_PDLAED3.
!
      CALL FS_GRIDINFO(SUBTREE, NPROW, NPCOL, MYROW, MYCOL)
      NBLK    = FS_get_NBLK(SUBTREE)
      NB      = FS_get_NB(SUBTREE)
      NP      = (NBLK / NPROW) * NB
      NQ      = (NBLK / NPCOL) * NB
      LDQ2    = MAX(NP,1)
      LDU     = NQ
      LSENDQ2 = LDQ2*NQ
!
      IZ      = 1
      IDLAMDA = IZ + N
      IW      = IDLAMDA + N
      IPQ2    = IW + N
      IPU     = IPQ2 + LDQ2*NQ
      IBUF    = IPU + LDU*NQ
      ISENDQ2 = IBUF + 4*N
      IRECVQ2 = ISENDQ2 + LSENDQ2
      IEND    = IRECVQ2 + LSENDQ2 - 1
!
      ICTOT  = 1
      IPSM   = ICTOT + NPCOL*4
      INDX   = IPSM + NPCOL*4
      INDXC  = INDX + N
      INDXP  = INDXC + N
      INDCOL = INDXP + N
      COLTYP = INDCOL + N
      INDROW = COLTYP + N
      INDXR  = INDROW + N
      INDXCB = INDXR + N
      IIEND  = INDXCB + N - 1
!
!     for FS_MERGE_D, FS_PDLAEDZ, FS_REDUCE_ZD
      IZWORK = IZ + N*2
      IDWORK = IZWORK + N
!
!     merge D
!
#if TIMER_PRINT
      CALL FS_prof_start(prof, 31)
#endif
      call FS_MERGE_D(N, D, SUBTREE, WORK(IDWORK))
#if TIMER_PRINT
      CALL FS_prof_end(prof, 31)
#endif
!
!     Form the z-vector which consists of the last row of Q_1 and the
!     first row of Q_2.
!
      CALL FS_PDLAEDZ( N, N1, Q, LDQ, SUBTREE, WORK(IZWORK), prof )
!
!     MPI_ALLREDUCE D and Z
!
      CALL FS_REDUCE_ZD( N, SUBTREE, WORK(IZWORK), WORK(IZ), D, prof )
!     
!     Deflate eigenvalues.
!     
      CALL FS_PDLAED2( K, N, N1, D, Q, LDQ, SUBTREE, RHO,  &
                       WORK(IZ), WORK(IW), WORK(IDLAMDA), WORK(IPQ2), &
                       LDQ2, IWORK(INDX), IWORK(ICTOT), &
                       WORK(IBUF), IWORK(COLTYP), IWORK(INDCOL), &
                       IWORK(INDXC), IWORK(INDXP), IWORK(IPSM), &
                       INFO, prof )
!     
!     Solve Secular Equation.
!     
      IF( K.NE.0 ) THEN

        SC = NB !NBでブロッキング
!       SC = NQ !ブロッキングなし
        CALL FS_PDLAED3( K, N, N1, D, RHO, WORK(IDLAMDA), WORK(IW), Q, &
                         LDQ, SUBTREE, WORK(IPQ2), LDQ2, WORK(IPU), &
                         LDU, SC, IWORK(INDX), IWORK(ICTOT), &
                         WORK(ISENDQ2), WORK(IRECVQ2), LSENDQ2, &
                         WORK(IZ), WORK(IBUF), IWORK(INDROW), &
                         IWORK(INDCOL), IWORK(INDXC), IWORK(INDXR), &
                         IWORK(INDXCB), &
                         INFO, prof )

      END IF
!
 20   CONTINUE
!
#if TIMER_PRINT
      CALL FS_prof_end(prof, 30)
#endif
#ifdef _DEBUGLOG
!      CALL MPI_Barrier(SUBTREE%MERGE_COMM,ierr)
      if( FS_MYRANK.eq.0 ) then
        write(*,'(a,i0)') "FS_PDLAED1 end. INFO=", INFO
      endif
#endif
      RETURN
!
!     End of FS_PDLAED1
!
end subroutine FS_PDLAED1

