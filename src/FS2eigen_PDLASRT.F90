!>
!> @file   FS2eigen_PDLASRT.F90
!> @brief  subroutine FS2eigen_PDLASRT
!!>
!
!
!
!> subroutine FS2eigen_PDLASRT
!>
!> @brief @n
!>  Purpose @n
!>  ======= @n
!>  FS2eigen_PDLASRT Sort the numbers in D in increasing order and the @n
!>  corresponding vectors in Q.
!
!  Arguments
!  =========
!>
!> @param[in]     N        (global input) INTEGER @n
!>                         The number of columns to be operated on i.e the number of @n
!>                         columns of the distributed submatrix sub( Q ). N >= 0.
!>
!> @param[in,out] D        (global input/output) DOUBLE PRECISION array, dimension (N) @n
!>                         On exit, the number in D are sorted in increasing order.
!>
!> @param[in,out] Q        (input/output) DOUBLE PRECISION pointer into the local memory @n
!>                         to an array of dimension (LDQ, NQ). This array contains the   @n
!>                         local pieces of the distributed matrix sub( A ) to be copied  @n
!>                         from.
!>
!> @param[in]     LDQ      (local input) INTEGER @n
!>                         The leading dimension of the array Q.  LDQ >= max(1,NP).
!>
!> @param[in]     SUBTREE  (input) type(bt_node) @n
!>                         sub-tree information of merge block.
!>
!> @param         IBUF     (workspace) INTEGR array, dimension (FS_NBROW*FS_NBCOL)
!>
!> @param         RBUF     (workspace) DOUBLE PRECISION array, dimension (N)
!>
!> @param         TBUF     (workspace) TYPE(GpositionValue) array, dimension (FS_NBROW*FS_NBCOL)
!>
!> @param         INDX     (workspace) INTEGER array, dimension (N)
!>
!> @param[out]    INFO     (global output) INTEGER @n
!>                         = 0: successful exit    @n
!>                         /=0: error exit
!>
!> @param[out]    prof     (global output) type(FS_prof) @n
!>                         profiling information of each subroutines.
!>
!> @note This routine is modified from ScaLAPACK PDLASRT.f
!>

module FS2eigen
  use mpi
  implicit none
  integer,parameter :: GpositionValue_BYTE_SIZE=4+4+8 ! GpositionValueの大きさ
!  integer,parameter :: GpositionValue_BYTE_SIZE=4+4+4 ! GpositionValueの大きさ
  type GpositionValue
     integer GRow
     integer GCol
     real(kind(0.0d0)) MatrixValue
!     real MatrixValue
  end type GpositionValue
  type CommBuf
     integer Rank   ! eigen準拠なのでランク番号の開始は1から 
     integer Ndata
     integer ireq
     integer ista(MPI_STATUS_SIZE)
     logical iflag
     type(GpositionValue), pointer :: bufp(:)
!     type(GpositionValue), allocatable ::  buf(:)
!     type(GpositionValue), allocatable ::  bufp(:)
  end type CommBuf

  type RANKLIST
     integer index
     integer, pointer :: LID(:)
!     integer, allocatable :: LID(:)
!     integer, allocatable :: LCOL(:)
!     integer, allocatable :: LROW(:)
  end type RANKLIST

contains
  subroutine FS2eigen_init_send(NP,COMM_INFO,COMM_BUF)
    implicit none
    integer NP ! プロセス数
    integer COMM_INFO(NP) 
    type(CommBuf) COMM_BUF(:)
    integer I
    integer NRANK
    NRANK = 0
    do I = 1,NP
       if(COMM_INFO(I).ne.0)then
          NRANK=NRANK+1
          COMM_BUF(NRANK)%iflag = .true.
          COMM_BUF(NRANK)%Rank  = I
          !         COMM_BUF(NRANK)%Ndata = 0
          COMM_BUF(NRANK)%Ndata = COMM_INFO(I)
          COMM_BUF(NRANK)%bufp => null() 

          !          write(*,*)"COMM_INFO(I)",COMM_INFO(I)   
          ! allocate( COMM_BUF(NRANK)%buf( COMM_INFO(I) ) )
!          allocate( COMM_BUF(NRANK)%buf( 1 ) )
       endif
    enddo

  end subroutine FS2eigen_init_send

  subroutine FS2eigen_init_recv(NP,COMM_INFO,COMM_BUF)
    implicit none
    integer NP ! プロセス数
    integer COMM_INFO(NP) 
    type(CommBuf) COMM_BUF(:)
    integer I
    integer NRANK
    NRANK = 0
    do I = 1,NP
       if(COMM_INFO(I).ne.0)then
          NRANK=NRANK+1
          COMM_BUF(NRANK)%iflag = .true.
          COMM_BUF(NRANK)%Rank  = I
          !         COMM_BUF(NRANK)%Ndata = 0
          COMM_BUF(NRANK)%Ndata = COMM_INFO(I)
          !          write(*,*)"COMM_INFO(I)",COMM_INFO(I)   
!          allocate( COMM_BUF(NRANK)%buf( COMM_INFO(I) ) )
!          allocate( COMM_BUF(NRANK)%buf( 1 ) )
          COMM_BUF(NRANK)%bufp => null() 

       endif
    enddo

  end subroutine FS2eigen_init_recv

  ! subroutine FS2eigen_free_comm(NRANK,COMM_BUF)
  !   implicit none
  !   integer NRANK
  !   type(CommBuf) COMM_BUF(:)
  !   integer I

  !   do I=1,NRANK
  !      deallocate( COMM_BUF(I)%buf ) 
  !   enddo

  ! end subroutine FS2eigen_free_comm

!  subroutine FS2eigen_isend(COMM_SEND_DATA,comm)
  subroutine FS2eigen_isend(COMM_SEND_DATA,SENDBUF,comm)
    use mpi
    implicit none
    type(CommBuf) COMM_SEND_DATA
    type(GpositionValue) SENDBUF(:)
    integer comm
    integer ncount
    integer srank
    integer ierr

    ncount = GpositionValue_BYTE_SIZE * COMM_SEND_DATA%Ndata
    srank  = COMM_SEND_DATA%Rank-1
    COMM_SEND_DATA%iflag=.false.
    ! write(*,*)"Send",COMM_SEND_DATA%Rank, COMM_SEND_DATA%Ndata 
    ! return
    call MPI_Isend(&
!         COMM_SEND_DATA%buf, &
         SENDBUF, &
         ncount, &
         MPI_BYTE, &
         srank, &
         1, &  ! COMM_SEND_DATA%Rank 
         comm, &
         COMM_SEND_DATA%ireq, &
         ierr)

  end subroutine FS2eigen_isend

!  subroutine FS2eigen_send(COMM_SEND_DATA,comm)
  subroutine FS2eigen_send(COMM_SEND_DATA,SENDBUF,comm)
    use mpi
    implicit none
    type(CommBuf) COMM_SEND_DATA
    type(GpositionValue) SENDBUF(:)
    integer comm
    integer ncount
    integer srank
    integer ierr

    ncount = GpositionValue_BYTE_SIZE * COMM_SEND_DATA%Ndata
    srank  = COMM_SEND_DATA%Rank-1
    COMM_SEND_DATA%iflag=.false.
    ! write(*,*)"Send",COMM_SEND_DATA%Rank, COMM_SEND_DATA%Ndata 
    ! return
    call MPI_send(&
!         COMM_SEND_DATA%buf, &
         SENDBUF, &
         ncount, &
         MPI_BYTE, &
         srank, &
         1, &  ! COMM_SEND_DATA%Rank 
         comm, &
         ierr)

  end subroutine FS2eigen_send


  subroutine FS2eigen_irecv(COMM_RECV_DATA,RECVBUF,comm)
    use mpi
    implicit none
    type(CommBuf) COMM_RECV_DATA
    type(GpositionValue) RECVBUF(:)
    integer comm
    integer ncount
    integer rrank
    integer ierr

    ncount = GpositionValue_BYTE_SIZE * COMM_RECV_DATA%Ndata
    rrank  = COMM_RECV_DATA%Rank-1
    COMM_RECV_DATA%iflag=.false.
    ! write(*,*)"Recv",COMM_RECV_DATA%Rank, COMM_RECV_DATA%Ndata 
    ! return

    call MPI_Irecv(&
!         COMM_RECV_DATA%buf, &
         RECVBUF, &
         ncount, &
         MPI_BYTE, &
         rrank, &
         1, &  ! COMM_RECV_DATA%Rank 
         comm, &
         COMM_RECV_DATA%ireq, &
         ierr)

  end subroutine FS2eigen_irecv

end module FS2eigen

subroutine FS2eigen_PDLASRT( N, D, Q, LDQ, SUBTREE, &
      IBUF, RBUF, TBUF, &
      INDX, INFO, prof )
  use FS_libs_mod
  use FS_dividing_mod

  use eigen_libs_mod
  use eigen_devel_mod,only : eigen_get_grid_major 
  use FS2eigen
  !$ use omp_lib
  use mpi
  implicit none
  !
  !     .. Scalar Arguments ..
  integer          , intent(in)    :: N, LDQ
  type(bt_node)    , intent(in)    :: SUBTREE
  integer          , intent(out)   :: INFO
  type(FS_prof)    , intent(inout) :: prof
  !
  !     .. Array Arguments ..
  real(kind(0.0d0)), intent(inout) :: Q(LDQ,*)
  real(kind(0.0d0)), intent(inout) :: D(N)

  ! work
  integer,target :: IBUF(*)
  real(kind(0.0d0)),target :: RBUF(*)
  type(GpositionValue),target :: TBUF(*)
  integer :: INDX(*)

  !     .. Local Scalars ..
  integer :: I, J, K
  integer :: LROW,LCOL,GROW,GCOL
  integer :: FS_NP
  integer :: FS_NB,FS_NBLK
  integer :: FS_NBROW,FS_NBCOL
  integer :: FS_NPROW, FS_NPCOL
  integer :: FS_MYROW, FS_MYCOL 

  integer :: eigen_NP,eigen_MYRANK
  integer :: eigen_NPROW, eigen_NPCOL
  integer :: eigen_MYROW, eigen_MYCOL 
  character*1 :: eigen_GRID_major
  integer eigen_comm, eigen_x_comm, eigen_y_comm

  integer :: I0
  integer :: K0,K1
  integer :: PROW,PCOL,PN
  integer ierr 

  integer :: eigen_rank_XY2COMM ! function
  integer Ndata,Ndata0
  integer SEND_NRANK ! 送信相手の総数
  integer RECV_NRANK ! 受信相手の総数
  integer FS_NBROW_MAX,FS_NBCOL_MAX
  integer,allocatable :: COMM_SEND_INFO(:) 
  integer,allocatable :: COMM_RECV_INFO(:) 
  type(CommBuf),allocatable :: COMM_SEND_DATA(:)
  type(CommBuf),allocatable :: COMM_RECV_DATA(:)
  integer,allocatable :: LCOL2GCOL_INDEX(:) 
  integer,allocatable :: LROW2GROW_INDEX(:) 
  logical :: ICOMM_FLAG,ICOMM_FLAG0
  integer :: SEND_MYRANK
  integer :: SEND_MAXSIZE
  integer :: RECV_MAXSIZE
  type(GpositionValue),pointer ::  SENDBUF(:)
  type(RANKLIST),allocatable :: SENDRANK_LIST(:)
  ! -----
  real(kind(0.0d0)) :: PROF_TIME(40)
  real(kind(0.0d0)) :: STIME,ETIME
  integer ipointer

  SENDBUF => null() 
  INFO = -1

  PROF_TIME=0.0d0
  ! -----

#ifdef _DEBUGLOG
  if( FS_MYRANK.eq.0 ) then
     write(*,'(a)') "FS2eigen_PDLASRT start."
  endif
#endif
#if TIMER_PRINT
  CALL FS_prof_start(prof, 70)
#endif

  ! ----
  ! eigen側の情報           
  call eigen_get_procs(eigen_NP    ,eigen_NPROW,eigen_NPCOL) ! call eigen_get_procs(nnod, x_nnod, y_nnod)
  call eigen_get_id   (eigen_MYRANK,eigen_MYROW,eigen_MYCOL) ! call eigen_get_id   (inod, x_inod, y_inod)
  call eigen_get_grid_major(eigen_GRID_major)
  call eigen_get_comm(eigen_comm, eigen_x_comm, eigen_y_comm)
  allocate(COMM_SEND_INFO(eigen_NP))
  allocate(COMM_RECV_INFO(eigen_NP))

  !$OMP PARALLEL DO
  do I = 1,eigen_NP
     COMM_SEND_INFO(I) = 0
     COMM_RECV_INFO(I) = -1
  enddo

  ! -----
  ! 固有値を全プロセスに送信する
  call MPI_Bcast(D,N,MPI_DOUBLE_PRECISION,0,eigen_comm,ierr)

  ! -----
  !
  !     Sort the eigenvalues in D
  !
  CALL DLAPST( 'I', N, D, INDX, INFO )
  !$OMP PARALLEL PRIVATE(I)
  !$OMP DO
  DO I = 1, N
     RBUF( I ) = D( INDX(I) )
  END DO
  !$OMP END DO

  !$OMP DO
  DO I = 1, N
     D( I ) = RBUF( I )
  END DO
  !$OMP END DO
  !$OMP END PARALLEL

  ! call MPI_Barrier(eigen_comm,ierr)
  STIME=MPI_Wtime()

  if(FS_COMM_MEMBER)then
     CALL FS_GRIDINFO(SUBTREE, FS_NPROW, FS_NPCOL, FS_MYROW, FS_MYCOL)
     FS_NP=FS_NPROW*FS_NPCOL
     ! 行列内のブロック数
     FS_NBLK  = FS_get_NBLK(SUBTREE)
     ! ブロックのサイズ
     FS_NB    = FS_get_NB(SUBTREE)
     ! プロセスあたりの列数
     FS_NBROW    = (FS_NBLK / FS_NPROW) * FS_NB
     ! プロセスあたりの行数
     FS_NBCOL    = (FS_NBLK / FS_NPCOL) * FS_NB
  else
     FS_NP=0; FS_NBLK=0; FS_NB=0
     FS_MYROW=0; FS_NBROW=0; FS_NBROW_MAX=0
     FS_MYCOL=0; FS_NBCOL=0; FS_NBCOL_MAX=0
  endif

  ! call MPI_Barrier(eigen_comm,ierr)
  ETIME=MPI_Wtime()
  PROF_TIME(1)=ETIME-STIME 
  STIME=ETIME

  if(FS_COMM_MEMBER)then
     allocate( LCOL2GCOL_INDEX(FS_NBCOL) )
     !$OMP PARALLEL DO  PRIVATE(I,J)  
     do J = 1,FS_NBCOL
        LCOL2GCOL_INDEX(J)=-1
     enddo

     FS_NBCOL_MAX = FS_NBCOL
     do J = 1,FS_NBCOL
        LCOL = J
        GCOL = FS_INDXL2G('C',LCOL,FS_MYCOL,SUBTREE) 
        if(N<GCOL)then
           FS_NBCOL_MAX = J -1
           exit
        endif
     enddo
     ! -----
     allocate( LROW2GROW_INDEX(FS_NBROW) )
     !$OMP PARALLEL DO  PRIVATE(I,J)  
     do I = 1,FS_NBROW
        LROW2GROW_INDEX(I)=-1
     enddo

     ! プロセス数で割り切れない場合に拡張した領域を除いた計算領域の総数
     FS_NBROW_MAX=FS_NBROW
     do I = 1,FS_NBROW
        LROW = I
        GROW = FS_INDXL2G('R',LROW,FS_MYROW,SUBTREE) 
        if(N<GROW)then
           FS_NBROW_MAX = I -1 
           exit
        endif
     enddo

     !$OMP PARALLEL DO  PRIVATE(I,LROW,GROW  )
     do I = 1,FS_NBROW_MAX
        LROW = I
        GROW = FS_INDXL2G('R',LROW,FS_MYROW,SUBTREE) 
        LROW2GROW_INDEX(I) = GROW
     enddo

     ! -----
     ! ROW <--> X
     ! COL <--> Y
     ! GROW,GCOLはプロセス数で行列次数を割り切れない場合に対応するために拡張した次数の行列番号を返すので,
     ! 拡張した範囲を省くために行列次数を超えたらexitする

     !$OMP PARALLEL DO PRIVATE(I,J,K,LROW,GROW,LCOL,GCOL,PROW,PCOL,PN) REDUCTION(+:COMM_SEND_INFO) schedule(dynamic,1)
     do J = 1,FS_NBCOL_MAX

        LCOL = J
        ! 固有値を並び替える前の列番号を取得
        GCOL = FS_INDXL2G('C',LCOL,FS_MYCOL,SUBTREE) 

        ! 固有値を並び変えた後の列番号に変換
        do K=1,N
           if(INDX(K).eq.GCOL)then
              GCOL = K
              LCOL2GCOL_INDEX(J) = K
              exit
           endif
        enddo

        ! グローバル情報から再分散後に担当するノードを求める
        PCOL = eigen_owner_node(GCOL,eigen_NPCOL,eigen_MYCOL) ! 引数にeigen_MYCOLがあるが参照されない
        do I = 1,FS_NBROW_MAX
           LROW = I
           GROW = LROW2GROW_INDEX(I)
           PROW = eigen_owner_node(GROW,eigen_NPROW,eigen_MYROW) ! 引数にeigen_MYROWがあるが参照されない
           PN   = eigen_rank_XY2COMM(eigen_GRID_major,PROW,PCOL)
           COMM_SEND_INFO(PN)=COMM_SEND_INFO(PN)+1
        enddo
     enddo
  endif

  ! call MPI_Barrier(eigen_comm,ierr)
  ETIME=MPI_Wtime()
  PROF_TIME(2)=ETIME-STIME 
  STIME=ETIME

  ! 各ランクごとの送信リストをALLTOALLすることで受信リストを作成する.
  call MPI_Alltoall( &
       COMM_SEND_INFO, 1, MPI_INTEGER, &
       COMM_RECV_INFO, 1, MPI_INTEGER, eigen_comm, ierr)
  ETIME=MPI_Wtime()
  PROF_TIME(3)=ETIME-STIME 
  STIME=ETIME

  SEND_NRANK=0
  SEND_MAXSIZE=1
  !$OMP PARALLEL DO PRIVATE(I) REDUCTION(+:SEND_NRANK) REDUCTION(max:SEND_MAXSIZE)
  do I = 1,eigen_NP
     if(COMM_SEND_INFO(I).ne.0)then
        SEND_NRANK=SEND_NRANK+1
        if(SEND_MAXSIZE<COMM_SEND_INFO(I))then
           SEND_MAXSIZE=COMM_SEND_INFO(I)
        endif
     endif
  enddo

  ! call MPI_Barrier(eigen_comm,ierr)
  ETIME=MPI_Wtime()
  PROF_TIME(4)=ETIME-STIME 
  STIME=ETIME

  ipointer=1
  if(SEND_NRANK.ne.0)then
     allocate( COMM_SEND_DATA(SEND_NRANK) )
     call FS2eigen_init_send(eigen_NP,COMM_SEND_INFO,COMM_SEND_DATA)

     allocate( SENDRANK_LIST(SEND_NRANK) )
     do K = 1,SEND_NRANK
#if 1
        SENDRANK_LIST(K)%LID=>IBUF(ipointer:ipointer+COMM_SEND_DATA(K)%Ndata-1)
#else
        allocate( SENDRANK_LIST(K)%LID( COMM_SEND_DATA(K)%Ndata ) )
#endif
        ipointer = ipointer + COMM_SEND_DATA(K)%Ndata
        do J = 1,COMM_SEND_DATA(K)%Ndata
           SENDRANK_LIST(K)%LID(J)=-1 
        enddo
     enddo
     ipointer = (ipointer + 4) / 4 + 1

#if 1
     SENDBUF => TBUF(ipointer:ipointer+SEND_MAXSIZE-1)
#else
     allocate(SENDBUF(SEND_MAXSIZE))
#endif
     ipointer = ipointer + SEND_MAXSIZE
     !$OMP PARALLEL DO PRIVATE(I) 
     do I = 1,SEND_MAXSIZE
        SENDBUF(I)%GRow=-1
        SENDBUF(I)%GCol=-1
        SENDBUF(I)%MatrixValue=0.0d0
     enddo

  endif

  ! call MPI_Barrier(eigen_comm,ierr)
  ETIME=MPI_Wtime()
  PROF_TIME(5)=ETIME-STIME 
  STIME=ETIME
  RECV_NRANK=0
  RECV_MAXSIZE=1
  !$OMP PARALLEL DO PRIVATE(I) REDUCTION(+:RECV_NRANK) REDUCTION(max:RECV_MAXSIZE)
  do I = 1,eigen_NP
     if(COMM_RECV_INFO(I).ne.0)then
        RECV_NRANK=RECV_NRANK+1
        if(RECV_MAXSIZE<COMM_RECV_INFO(I))then
           RECV_MAXSIZE=COMM_RECV_INFO(I)
        endif
      endif
  enddo

  ETIME=MPI_Wtime()
  PROF_TIME(6)=ETIME-STIME 
  STIME=ETIME

  ! ----
  ! 受信バッファの設定
  if(.not. FS_COMM_MEMBER)then
     ipointer = 1
  endif
  if(RECV_NRANK.ne.0)then
     allocate(COMM_RECV_DATA(RECV_NRANK))
     call FS2eigen_init_recv(eigen_NP,COMM_RECV_INFO,COMM_RECV_DATA)
     do K = 1,RECV_NRANK
        if(COMM_RECV_DATA(K)%Rank.ne.eigen_MYRANK)then
#if 1
           COMM_RECV_DATA(K)%bufp => TBUF(ipointer:ipointer+COMM_RECV_DATA(K)%Ndata-1)
#else
           allocate(COMM_RECV_DATA(K)%bufp(COMM_RECV_DATA(K)%Ndata) ,stat=ierr) 
#endif
           ipointer=ipointer+COMM_RECV_DATA(K)%Ndata
           call FS2eigen_irecv(COMM_RECV_DATA(K),COMM_RECV_DATA(K)%bufp,eigen_comm)
        endif
     enddo
  endif

  ! ----
  ! 送信データのパック
  if(FS_COMM_MEMBER)then

     ! -----
     ! COMM_SEND_INFOを送信先のプロセスランクからローカルな送信先インデックスに変換するリストを作成
     do K = 1,SEND_NRANK
        PN = COMM_SEND_DATA(K)%Rank
        COMM_SEND_INFO(PN) = K
        SENDRANK_LIST(K)%index = 0
     end do

     ! -----
     ! 各スレッドが担当する送信先プロセスごとのデータ数を計算
     do J = 1,FS_NBCOL_MAX
        do I = 1,FS_NBROW_MAX
           LCOL = J
           ! 固有値を並び変えた後の列番号に変換
           GCOL = LCOL2GCOL_INDEX(J) 
           PCOL = eigen_owner_node(GCOL,eigen_NPCOL,eigen_MYCOL) ! 引数にeigen_MYCOLがあるが参照されない

           LROW = I
           GROW = LROW2GROW_INDEX(I)

           PROW = eigen_owner_node(GROW,eigen_NPROW,eigen_MYROW) ! 引数にeigen_MYROWがあるが参照されない
           ! 送信先のランク
           PN   = eigen_rank_XY2COMM(eigen_GRID_major,PROW,PCOL) 

           SENDRANK_LIST( COMM_SEND_INFO(PN) )%index =  SENDRANK_LIST( COMM_SEND_INFO(PN) )%index + 1 
           ! 送信先ランクのグローバル行,列番号 
           SENDRANK_LIST( COMM_SEND_INFO(PN) )%LID( SENDRANK_LIST( COMM_SEND_INFO(PN) )%index ) = ( (I-1) + (J-1)*FS_NBROW_MAX ) + 1

        enddo
     enddo

     ETIME=MPI_Wtime()
     PROF_TIME(12)=ETIME-STIME 
     STIME=ETIME

     !------
     ! 送信先が被らないように初回の通信相手を選択する
     I0=1
     K0=1
     if(0<SEND_NRANK)then
        I0=-1
        K0=-1
        ! ------
        do K = 1,eigen_NP
           K1 = mod(K+eigen_MYRANK-1,eigen_NP) + 1
           do I = 1,SEND_NRANK
              if(K1.eq.COMM_SEND_DATA(I)%Rank)then
                 I0=I ! 最初の通信相手
                 exit
              endif
           enddo
           if(I0.ne.-1)then
              exit
           endif
        enddo
     endif

     do K0 = I0,SEND_NRANK+I0-1
        K=mod(K0,SEND_NRANK) + 1
        if(COMM_SEND_DATA(K)%Rank.ne.eigen_MYRANK)then
           do J = 1,COMM_SEND_DATA(K)%Ndata
              LCOL = (SENDRANK_LIST( K )%LID( J ) - 1) / FS_NBROW_MAX + 1
              GCOL = LCOL2GCOL_INDEX(LCOL) 
              LROW = mod( (SENDRANK_LIST( K )%LID( J ) - 1) ,FS_NBROW_MAX) + 1
              GROW = LROW2GROW_INDEX(LROW)

              SENDBUF(J)%GRow = GROW
              SENDBUF(J)%GCol = GCOL
              SENDBUF(J)%MatrixValue = Q(LROW, LCOL)
           enddo
           call FS2eigen_send(COMM_SEND_DATA(K),SENDBUF,eigen_comm)
        endif
     enddo

     ! Recv側のWait
     do K0 = 1,RECV_NRANK
        if(COMM_RECV_DATA(K0)%Rank.ne.eigen_MYRANK)then
           call MPI_Wait(COMM_RECV_DATA(K0)%ireq, MPI_STATUS_IGNORE, ierr)
        endif
     enddo

     ! ----
     ! 受信データのアンパック
     ! 送信先が自分自身であればデータコピー

     do K0=1,SEND_NRANK
        if(eigen_MYRANK.eq.COMM_SEND_DATA(K0)%Rank)then
           do I=1,COMM_SEND_DATA(K0)%Ndata
              LCOL = (SENDRANK_LIST( K0 )%LID( I ) - 1) / FS_NBROW_MAX + 1
              GCOL = LCOL2GCOL_INDEX( LCOL )
              LROW = mod( (SENDRANK_LIST( K0 )%LID( I ) - 1) ,FS_NBROW_MAX) + 1
              GROW = LROW2GROW_INDEX( LROW )
              
              SENDBUF(I)%GRow = GROW
              SENDBUF(I)%GCol = GCOL
              SENDBUF(I)%MatrixValue = Q(LROW, LCOL)
           enddo
           do I=1,COMM_SEND_DATA(K0)%Ndata
              GCOL = SENDBUF(I)%GCol
              GROW = SENDBUF(I)%GRow 
              LCOL = eigen_translate_g2l(GCOL,eigen_NPCOL,eigen_MYCOL)
              LROW = eigen_translate_g2l(GROW,eigen_NPROW,eigen_MYROW)
              Q(LROW, LCOL) = SENDBUF(I)%MatrixValue 
           enddo
           exit
        endif
     enddo

     !$OMP PARALLEL DO PRIVATE(K,K0,I,Ndata,GCOL,GROW,LCOL,LROW),schedule(dynamic,1)
     do K = 1,RECV_NRANK
        if(COMM_RECV_DATA(K)%Rank.ne.eigen_MYRANK)then
           Ndata=COMM_RECV_DATA(K)%Ndata
           do I=1,Ndata
              GCOL = COMM_RECV_DATA(K)%bufp(I)%GCol
              GROW = COMM_RECV_DATA(K)%bufp(I)%GRow              
              LCOL = eigen_translate_g2l(GCOL,eigen_NPCOL,eigen_MYCOL)
              LROW = eigen_translate_g2l(GROW,eigen_NPROW,eigen_MYROW)
              Q(LROW, LCOL) = COMM_RECV_DATA(K)%bufp(I)%MatrixValue
           enddo
        endif
     enddo

  else
     if(RECV_NRANK.ne.0)then
        do K = 1 , RECV_NRANK
           call MPI_Wait(COMM_RECV_DATA(K)%ireq, MPI_STATUS_IGNORE, ierr)              
           Ndata=COMM_RECV_DATA(K)%Ndata
           !$OMP PARALLEL DO PRIVATE(I,GCOL,GROW,LCOL,LROW)
           do I=1,Ndata
              GCOL = COMM_RECV_DATA(K)%bufp(I)%GCol
              GROW = COMM_RECV_DATA(K)%bufp(I)%GRow                 
              LCOL = eigen_translate_g2l(GCOL,eigen_NPCOL,eigen_MYCOL)
              LROW = eigen_translate_g2l(GROW,eigen_NPROW,eigen_MYROW)
              Q(LROW, LCOL) = COMM_RECV_DATA(K)%bufp(I)%MatrixValue
           enddo
        enddo
     endif
  endif
      
  ! call MPI_Barrier(eigen_comm,ierr)
  ETIME=MPI_Wtime()
  PROF_TIME(7)=ETIME-STIME 
  STIME=ETIME

  ! call MPI_Barrier(eigen_comm,ierr)
  ETIME=MPI_Wtime()
  PROF_TIME(9)=ETIME-STIME 
  STIME=ETIME

#if TIMER_PRINT
  CALL FS_prof_end(prof, 70)
#endif

#ifdef _DEBUGLOG
  if( FS_MYRANK.eq.0 ) then
     write(*,'(a)') "FS_PDLASRT end."
  endif
#endif

  if(SEND_NRANK.ne.0)then
     deallocate(COMM_SEND_DATA)
     do K = 1,SEND_NRANK
!        SENDRANK_LIST(K)%LID => null() 
        nullify( SENDRANK_LIST(K)%LID )
     enddo
#if 1
     !     SENDBUF  => null()
     nullify( SENDBUF )
#else
     deallocate( SENDBUF )
#endif

     deallocate(SENDRANK_LIST)
  end if
  if(RECV_NRANK.ne.0)then
     do K = 1,RECV_NRANK
!        COMM_RECV_DATA(K)%bufp => null() 
        nullify( COMM_RECV_DATA(K)%bufp )
     enddo
     deallocate(COMM_RECV_DATA)
   end if

  if(FS_COMM_MEMBER)then
     deallocate(LCOL2GCOL_INDEX)
     deallocate(LROW2GROW_INDEX)     
  endif

  ! call MPI_Barrier(eigen_comm,ierr)
  ETIME=MPI_Wtime()
  PROF_TIME(10)=ETIME-STIME 
  STIME=ETIME
  INFO = 0
  return
end subroutine FS2eigen_PDLASRT

! eigenの2次元グリッドランク(x_inod,y_inod)に対応する1次元ランク(inod)を返す
integer function eigen_rank_XY2COMM(GRID_major,x_inod,y_inod) & 
     result(inod)
  use eigen_libs_mod,only : eigen_get_procs
  implicit none

  character*1 :: GRID_major
  ! eigen準拠なのでランク番号の開始は1から 
  integer x_inod ! input  : 2次元グリッドのx方向プロセスランク
  integer y_inod ! input  : 2次元グリッドのy方向プロセスランク
  integer nnod, x_nnod, y_nnod

  call eigen_get_procs(nnod, x_nnod, y_nnod)
  if (GRID_major == 'R') then
     inod = y_inod + (x_inod-1) * y_nnod
  else
     inod = x_inod + (y_inod-1) * x_nnod
  end if
  return
end function eigen_rank_XY2COMM
