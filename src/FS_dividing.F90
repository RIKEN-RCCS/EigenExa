!>
!> @file   FS_dividing.F90
!> @brief  module FS_dividing_mod
!!>
!
!
!
!
!> module FS_dividing_mod
!> @brief  @n
!> Purpose @n
!> ======= @n
!> FS_dividing_mod divide process grid recursive and computes rank-1 modification
module FS_dividing_mod
      use mpi
      use FS_libs_mod
      use FS_prof_mod
      implicit none

      !> type bt_node
      !> @brief structure of tree node
      type, public :: bt_node
        integer                :: layer                = 1        !< layer number (>=1)
        logical                :: direction_horizontal = .true.   !< divided direction (.true.:horizontal)
        integer                :: nstart               = 1        !< start index of N at own layer
        integer                :: nend                 = 1        !< end   index of N at own layer(for extended N)
        integer                :: nend_active          = 1        !< end   index of N at own layer
        integer                :: proc_istart          = 1        !< process start number of direction i
        integer                :: proc_iend            = 1        !< process end   number of direction i
        integer                :: proc_jstart          = 1        !< process start number of direction j
        integer                :: proc_jend            = 1        !< process end   number of direction j
        integer                :: block_start          = 1        !< merge block start number(refer to procs_i/j)
        integer                :: block_end            = 1        !< merge block end   number(refer to procs_i/j)
        type(bt_node),pointer  :: sub_bt_node(:)                  !< sub tree node
        type(bt_node),pointer  :: parent_node                     !< parent node
        integer,pointer        :: procs_i(:)                      !< process No. list of row
        integer,pointer        :: procs_j(:)                      !< process No. list of column
        integer                :: nnod   = 0                      !< nprocs of communicator
        integer                :: x_nnod = 0                      !< nprocs of X direction communicator
        integer                :: y_nnod = 0                      !< nprocs of Y direction communicator
        integer                :: inod   = 0                      !< inod in MERGE_COMM(1～)
        integer                :: x_inod = 0                      !< x_inod in MERGE_COMM_X(1～)
        integer                :: y_inod = 0                      !< y_inod in MERGE_COMM_Y(1～)
        integer                :: div_bit  = -1                   !< bit stream of divided direction
        integer                :: div_nbit = 0                    !< number of dights of div_bit

        integer                :: MERGE_GROUP   = MPI_UNDEFINED   !< MERGE_COMM group
        integer                :: MERGE_GROUP_X = MPI_UNDEFINED   !< MERGE_COMM_X group
        integer                :: MERGE_GROUP_Y = MPI_UNDEFINED   !< MERGE_COMM_Y group
        integer,allocatable    :: group_processranklist(:)        !< list to convert from group rank to communicator rank  
        integer,allocatable    :: group_X_processranklist(:)      !< list to convert from X group rank to communicator rank  
        integer,allocatable    :: group_Y_processranklist(:)      !< list to convert from Y group rank to communicator rank  

      end type bt_node

      integer bt_id ! loop variable 
contains

  !---------------------------------------------------------------------------------
  !> subroutine FS_create_hint
  !> @brief create default hint of dividind tree
  !> @param[out] hint   (output) creted hint. dimension = number of tree layer
  subroutine FS_create_hint(hint)
        implicit none
        logical, intent(out) :: hint(*)

        integer :: nnod,x_nnod,y_nnod
        integer :: layer

        ! process grid
        call FS_get_procs(nnod,x_nnod,y_nnod)

        ! while
        layer = 0

#if _TREEDIV==1
        ! 縦分割優先
        ! 縦->横->縦->横->...
        do while( x_nnod*y_nnod.ge.1 )
          layer = layer + 1
          if( y_nnod.ge.x_nnod ) then
            hint(layer) = .false.
            y_nnod = y_nnod / 2
          else
            hint(layer) = .true.
            x_nnod = x_nnod / 2
          endif
        enddo
#elif _TREEDIV==2
        ! 横分割先行
        ! 横->横->縦->縦
        do while( x_nnod*y_nnod.ge.1 )
          layer = layer + 1
          if( x_nnod.ge.2 ) then
            hint(layer) = .true.
            x_nnod = x_nnod / 2
          else
            hint(layer) = .false.
            y_nnod = y_nnod / 2
          endif
        enddo
#elif _TREEDIV==3
        ! 縦分割先行
        ! 縦->縦->横->横
        do while( x_nnod*y_nnod.ge.1 )
          layer = layer + 1
          if( y_nnod.ge.2 ) then
            hint(layer) = .false.
            y_nnod = y_nnod / 2
          else
            hint(layer) = .true.
            x_nnod = x_nnod / 2
          endif
        enddo
#else
        ! 横分割優先(デフォルト)
        ! 横->縦->横->縦->...
        do while( x_nnod*y_nnod.ge.1 )
          layer = layer + 1
          if( x_nnod.ge.y_nnod ) then
            hint(layer) = .true.
            x_nnod = x_nnod / 2
          else
            hint(layer) = .false.
            y_nnod = y_nnod / 2
          endif
        enddo
#endif

#ifdef _DEBUGLOG
        if( FS_myrank.eq.0 ) then
          write(*,'(a,1000(1x,l))') "procdiv=",hint(1:layer)
        endif
#endif

        return
  end subroutine FS_create_hint

  !---------------------------------------------------------------------------------
  !> subroutine FS_dividing
  !> @brief main routine of dividing tree
  !>
  !> @param[in]     n      (global input) INTEGER @n
  !>                       The order of the tridiagonal matrix T.  N >= 0.
  !>
  !> @param[in,out] d      (global input/output) DOUBLE PRECISION array, dimension (N) @n
  !>                       On entry, the diagonal elements of the tridiagonal matrix.  @n
  !>                       On exit, rank-1 modification.
  !>
  !> @param[in]     e      (global input) DOUBLE PRECISION array, dimension (N-1) @n
  !>                       the subdiagonal elements of the tridiagonal matrix.
  !>
  !> @param[out]    tree   (output) type(bt_node) @n
  !>                       tree information
  !> @param[in]     hint   (input) LOGICAL array, dimension = number of tree layer @n
  !>                       tree divide pattern
  !>
  !> @param[out]    info   (global output) INTEGER @n
  !>                       = 0: successful exit   @n
  !>                       /=0: error exit
  !>
  !> @param[out]    prof   (global output) type(FS_prof) @n
  !>                       profiling information of each subroutines.
  !>
  subroutine FS_dividing(n,d,e,tree,hint,info,prof)
        implicit none

        integer          , intent(in)        :: n
        real(kind(0.0d0)), intent(inout)     :: d(1:n)
        real(kind(0.0d0)), intent(in)        :: e(1:n-1)
        type(bt_node)    , intent(out)       :: tree
        logical          , intent(in)        :: hint(*)
        integer          , intent(inout)     :: info
        type(FS_prof)    , intent(inout)     :: prof

        integer :: nnod,x_nnod,y_nnod
        integer :: Next, i
        integer :: ierr

#ifdef _DEBUGLOG
        if( FS_MYRANK.eq.0 ) then
          write(*,'(a)') "FS_dividing start."
        endif
#endif
#if TIMER_PRINT
        call FS_prof_start(prof, 21)
#endif

        ! process grid
        call FS_get_procs(nnod,x_nnod,y_nnod)

        ! check N
        Next = n
        if( mod(n,nnod).ne.0 ) then
          Next = (n / nnod + 1) * nnod
        endif

        ! root node setting
        tree%layer = 1
        tree%direction_horizontal = hint(tree%layer)
        tree%nstart = 1
        tree%nend   = Next
        tree%nend_active = n
        tree%proc_istart = 1
        tree%proc_iend   = x_nnod
        tree%proc_jstart = 1
        tree%proc_jend   = y_nnod
        tree%block_start = 1
        tree%block_end   = nnod
        allocate( tree%procs_i(nnod) )
        allocate( tree%procs_j(nnod) )
        tree%procs_i = -1
        tree%procs_j = -1
        tree%parent_node  => null()

        ! divide recursive      
        bt_id = 1
        call FS_dividing_recursive(n,d,e,tree,hint,info,prof)

        ! create communicator for merge block
#if TIMER_PRINT>1
        call FS_prof_start(prof, 22)
#endif
        call FS_create_merge_comm(tree)
#if TIMER_PRINT>1
        call FS_prof_end(prof, 22)
#endif

#if TIMER_PRINT>1
        call FS_prof_start(prof, 23)
        call MPI_Barrier(FS_COMM_WORLD,ierr)
        call FS_prof_end(prof, 23)
#endif

#ifdef _PRINTTREE
        if( FS_MYRANK.eq.0 ) then
          call print_tree(6, tree)
        endif
#endif

#if TIMER_PRINT
        call FS_prof_end(prof, 21)
#endif

#ifdef _DEBUGLOG
        if( FS_MYRANK.eq.0 ) then
          write(*,'(a,i0)') "FS_dividing end. INFO=", info
        endif
#endif

        return
  end subroutine FS_dividing


  !---------------------------------------------------------------------------------
  !> subroutine FS_dividing_recursive
  !> @brief create sub-tree recursive
  !>
  !> @param[in]     n      (global input) INTEGER @n
  !>                       The order of the tridiagonal matrix T.  N >= 0.
  !>
  !> @param[in,out] d      (global input/output) DOUBLE PRECISION array, dimension (N) @n
  !>                       On entry, the diagonal elements of the tridiagonal matrix.  @n
  !>                       On exit, rank-1 modification.
  !>
  !> @param[in]     e      (global input) DOUBLE PRECISION array, dimension (N-1) @n
  !>                       the subdiagonal elements of the tridiagonal matrix.
  !>
  !> @param[in,out] tree   (input/output) type(bt_node) @n
  !>                       tree information
  !>
  !> @param[in]     hint   (input) LOGICAL array, dimension = number of tree layer @n
  !>                       tree divide pattern
  !>
  !> @param[out]    info   (global output) INTEGER @n
  !>                       = 0: successful exit   @n
  !>                       /=0: error exit
  !>
  !> @param[out]    prof   (global output) type(FS_prof) @n
  !>                       profiling information of each subroutines.
  !>
  recursive subroutine FS_dividing_recursive(n,d,e,tree,hint,info,prof)
        implicit none

        integer          , intent(in)            :: n
        real(kind(0.0d0)), intent(inout)         :: d(1:n)
        real(kind(0.0d0)), intent(in)            :: e(1:n-1)
        type(bt_node)    , intent(inout), target :: tree
        logical          , intent(in)            :: hint(*)
        integer          , intent(inout)         :: info
        type(FS_prof)    , intent(inout)         :: prof

        integer :: i
        integer :: lnod,x_lnod,y_lnod !ローカル
        integer :: nnod,x_nnod,y_nnod !全体
        integer :: inod,x_inod,y_inod
        type(bt_node), pointer :: subptr

        info = 0
        tree%sub_bt_node  => null()

        call FS_get_procs(nnod,x_nnod,y_nnod)
        call FS_get_id(inod,x_inod,y_inod)

        x_lnod = tree%proc_iend - tree%proc_istart + 1
        y_lnod = tree%proc_jend - tree%proc_jstart + 1
        lnod = x_lnod * y_lnod
        if( lnod.eq.0 ) then
          ! error
          info = -1
          goto 10
        endif

        if( lnod.eq.1 ) then
          ! reach to leaf

          ! grid process No. of row & column
          ! merge block position
!          do i=1,nnod
           do i=bt_id,nnod
            if( tree%procs_i(i) < 0 ) then

              ! grid process No. of row & column
              tree%procs_i(i) = tree%proc_istart
              tree%procs_j(i) = tree%proc_jstart

              ! merge block position
              tree%block_start = i
              tree%block_end   = i

              ! rank-1 modification
              if( tree%nend.lt.n ) then
                d(tree%nend  ) = d(tree%nend  ) - abs( e(tree%nend) )
                d(tree%nend+1) = d(tree%nend+1) - abs( e(tree%nend) )
              endif

              exit
            endif
          enddo
          bt_id = i
          goto 10
        endif

        if( (      hint(tree%layer) .and. x_lnod.eq.1 ) .or. &
            (.not. hint(tree%layer) .and. y_lnod.eq.1 ) ) then
          ! error
          info = -2
          goto 10
        endif

        ! divide
        allocate( tree%sub_bt_node(2) )
        do i=1,2
          subptr => tree%sub_bt_node(i)
          subptr%layer = tree%layer+1
          subptr%direction_horizontal = hint(subptr%layer)
          subptr%nstart = tree%nstart + (i-1) * (tree%nend - tree%nstart + 1) / 2
          subptr%nend   = tree%nend   - (2-i) * (tree%nend - tree%nstart + 1) / 2
          subptr%nend_active = max(min(subptr%nend, n), subptr%nstart-1)
          if( hint(tree%layer) ) then
            subptr%proc_istart = tree%proc_istart + (i-1) * (tree%proc_iend - tree%proc_istart + 1) / 2
            subptr%proc_iend   = tree%proc_iend   - (2-i) * (tree%proc_iend - tree%proc_istart + 1) / 2
            subptr%proc_jstart = tree%proc_jstart
            subptr%proc_jend   = tree%proc_jend
          else
            subptr%proc_istart = tree%proc_istart
            subptr%proc_iend   = tree%proc_iend
            subptr%proc_jstart = tree%proc_jstart + (i-1) * (tree%proc_jend - tree%proc_jstart + 1) / 2
            subptr%proc_jend   = tree%proc_jend   - (2-i) * (tree%proc_jend - tree%proc_jstart + 1) / 2
          endif
          subptr%parent_node => tree
          subptr%procs_i => tree%procs_i
          subptr%procs_j => tree%procs_j

          ! divide recursive      
          call FS_dividing_recursive(n,d,e,subptr,hint,info,prof)
          if( info.ne.0 ) then
            goto 10
          endif
        enddo

        ! merge block position
        tree%block_start = min( tree%sub_bt_node(1)%block_start, tree%sub_bt_node(2)%block_start )
        tree%block_end   = max( tree%sub_bt_node(1)%block_end  , tree%sub_bt_node(2)%block_end   )

        ! set bit stream of divided direction
        call FS_dividing_setBitStream(tree)

  10    continue
        return
  end subroutine FS_dividing_recursive


  !---------------------------------------------------------------------------------
  !> subroutine FS_dividing_setBitStream
  !> @brief set bit stream of tree dividing direction of all child node to leaf
  !>
  !> @param[in,out] tree   (input/output) type(bt_node) @n
  !>                       tree information
  !>
  subroutine FS_dividing_setBitStream(tree)
        implicit none
        type(bt_node), intent(inout), target  :: tree
        type(bt_node),pointer :: node

        if( .not. associated( tree%sub_bt_node ) ) then
          return
        endif

        tree%div_bit = 0
        if( .not. tree%direction_horizontal ) then
          tree%div_bit = IBSET(tree%div_bit,0)
        endif
        tree%div_nbit = 1

        node => tree%sub_bt_node(1)
        do while( associated( node%sub_bt_node ) )
          tree%div_bit = ISHFT(tree%div_bit,1)
          if( .not. node%direction_horizontal ) then
            tree%div_bit = IBSET(tree%div_bit,0)
          endif
          tree%div_nbit = tree%div_nbit + 1
          node=>node%sub_bt_node(1)
        enddo

        return
  end subroutine FS_dividing_setBitStream


  !---------------------------------------------------------------------------------
  !> function FS_node_included
  !> @brief 自プロセスがノードに含まれるかチェックする
  !> @param[in] node   (input) tree node
  !> @retval    true    含まれる
  !> @retval    false   含まれない
  logical function FS_node_included(node)
    implicit none
    type(bt_node), intent(in), target  :: node
    integer :: inod,x_inod,y_inod

    FS_node_included = .false.
    call FS_get_id(inod,x_inod,y_inod)

    if( x_inod .lt. node%proc_istart .or. x_inod .gt. node%proc_iend ) then
      return
    endif

    if( y_inod .lt. node%proc_jstart .or. y_inod .gt. node%proc_jend ) then
      return
    endif

    FS_node_included = .true.

    return
  end function FS_node_included


  !---------------------------------------------------------------------------------
  !> subroutine FS_dividing_getleaf
  !> @brief search leaf node of own process recursive.
  !>
  !> @param[in]     tree   (input) type(bt_node) @n
  !>                       tree information
  !>
  !> @param[out]    leaf   (output) type(bt_node) @n
  !>                       leaf node pointer
  !>
  !> @param[out]    info   (global output) INTEGER @n
  !>                       = 0: successful exit   @n
  !>                       /=0: error exit
  !>
  recursive subroutine FS_dividing_getleaf(tree,leaf,info)
        implicit none

        type(bt_node), intent(in) , target  :: tree
        type(bt_node), intent(out), pointer :: leaf
!        integer, intent(out) :: info
        integer, intent(inout) :: info

        integer :: nnod,x_nnod,y_nnod 
        integer :: inod,x_inod,y_inod
        integer :: i

        if( associated(leaf) ) then
          return
        endif

        call FS_get_procs(nnod,x_nnod,y_nnod)
        call FS_get_id(inod,x_inod,y_inod)

        if( x_inod.eq.tree%proc_istart .and. x_inod.eq.tree%proc_iend .and. &
            y_inod.eq.tree%proc_jstart .and. y_inod.eq.tree%proc_jend ) then
          leaf => tree
          return
        else
          if( associated( tree%sub_bt_node ) ) then
            do i=1,2
              call FS_dividing_getleaf(tree%sub_bt_node(i),leaf,info)
              if( associated(leaf) ) then
                return
              endif
            enddo
          endif
        endif

        if( tree%layer.eq.1 .and. .not. associated(leaf) ) then
          info = 9999
        endif

        return
  end subroutine FS_dividing_getleaf

  !---------------------------------------------------------------------------------
  !> subroutine FS_create_merge_comm
  !> @brief create local merge group 
  !>
  !> @param[in,out] node   (input) type(bt_node) @n
  !>                       top tree pointer
  !>
  subroutine FS_create_merge_comm(node)
        implicit none
        type(bt_node), intent(inout) :: node

        integer :: inod,x_inod,y_inod
        integer :: nnod,x_nnod,y_nnod
        integer :: ierr

        ! get process info
        call FS_get_id(inod,x_inod,y_inod)
        call FS_get_procs(nnod, x_nnod, y_nnod)

        node%MERGE_GROUP = FS_GROUP
        ! create X,Y group
        if( node%MERGE_GROUP.ne.MPI_UNDEFINED ) then
          node%inod   = inod
          node%x_inod = x_inod
          node%y_inod = y_inod
          node%nnod   = nnod
          node%x_nnod = x_nnod
          node%y_nnod = y_nnod
          call FS_create_mergeXY_group(node)
        endif

        ! subtree
        call FS_create_merge_comm_recursive(node)

        return
  end subroutine FS_create_merge_comm

  !---------------------------------------------------------------------------------
  !> subroutine FS_create_mergeXY_group
  !> @brief create local merge X,Y group 
  !>
  !> @param[in,out] PARENT_NODE   (input) type(bt_node)  node pointer
  !>
  subroutine FS_create_mergeXY_group(node)
    implicit none
    type(bt_node), intent(inout) :: node

    integer :: inod,x_inod,y_inod
    integer :: ierr
    character*1 :: order

    integer :: nnod,x_nnod,y_nnod
    integer :: nrank
    integer :: ii,jj
    integer,allocatable :: ranklist_group(:)

    call FS_get_id(inod,x_inod,y_inod)
    call FS_get_grid_major(order)
    call FS_get_procs(nnod, x_nnod, y_nnod)
    allocate(ranklist_group(nnod))

    ! --------
    ! i <--> x

    if( order == 'R' ) then
       nrank=0
       jj = y_inod
       do ii=node%proc_istart, node%proc_iend
          nrank=nrank+1
          ranklist_group(nrank) = (ii-1) * y_nnod             + (jj-1)
       enddo
    else
       nrank=0
       jj = y_inod
       do ii=node%proc_istart, node%proc_iend
          nrank=nrank+1
          ranklist_group(nrank) = (jj-1) * x_nnod             + (ii-1)
       enddo
    endif
    call MPI_Group_incl(FS_GROUP, nrank, ranklist_group, node%MERGE_GROUP_X, ierr)

    do ii=1,nrank
       ranklist_group(ii)=ii-1
    enddo
    allocate(node%group_X_processranklist(nrank))
    call MPI_GROUP_TRANSLATE_RANKS(&
         node%MERGE_GROUP_X, &
         nrank, &
         ranklist_group, &
         FS_GROUP, &
         node%group_X_processranklist, &
         ierr)

    ! j <--> y
    if( order == 'R' ) then
       nrank=0
       ii = x_inod
       do jj=node%proc_jstart, node%proc_jend
          nrank=nrank+1
          ranklist_group(nrank) = (ii-1) * y_nnod             + (jj-1)
       enddo
    else
       nrank=0
       ii = x_inod
       do jj=node%proc_jstart, node%proc_jend
          nrank=nrank+1
          ranklist_group(nrank) = (jj-1) * x_nnod             + (ii-1)
       enddo
    endif
    call MPI_Group_incl(FS_GROUP, nrank, ranklist_group, node%MERGE_GROUP_Y, ierr)
    do ii=1,nrank
       ranklist_group(ii)=ii-1
    enddo
    allocate(node%group_Y_processranklist(nrank))
    call MPI_GROUP_TRANSLATE_RANKS(&
         node%MERGE_GROUP_Y, &
         nrank, &
         ranklist_group, &
         FS_GROUP, &
         node%group_Y_processranklist, &
         ierr)

    deallocate(ranklist_group)
  end subroutine FS_create_mergeXY_group

  !---------------------------------------------------------------------------------
  !> subroutine FS_create_merge_comm_recursive
  !> @brief create local merge group (recursive)
  !>
  !> @param[in] PARENT_NODE   (input) type(bt_node)  node pointer
  !>
  recursive subroutine FS_create_merge_comm_recursive(PARENT_NODE)
        implicit none
        type(bt_node), intent(in) :: PARENT_NODE

        integer :: nnod,x_nnod,y_nnod !全体
        integer :: inod,x_inod,y_inod

        type(bt_node), pointer :: node
        integer :: group_world, group_merge
        integer :: nrank, ni, nj
        integer,allocatable :: ranklist(:)
        integer :: n, i, j, ii, jj, ierr
        character*1 :: order
        integer :: NEW_COMM

        integer,allocatable :: ranklist_group(:)

        ! get process info
        call FS_get_id(inod,x_inod,y_inod)
        call FS_get_grid_major(order)
        call FS_get_procs(nnod, x_nnod, y_nnod)

        if( .not. associated( PARENT_NODE%sub_bt_node ) ) then
          return
        endif

        ! child node loop
        do n=1,2
          node => PARENT_NODE%sub_bt_node(n)

          ! create rank list
          ni = node%proc_iend - node%proc_istart + 1
          nj = node%proc_jend - node%proc_jstart + 1
          nrank = ni * nj
          allocate( ranklist(nrank) )
          allocate( ranklist_group(nrank) )
          nrank = 0
          if( order == 'R' ) then
            do ii=node%proc_istart, node%proc_iend
            do jj=node%proc_jstart, node%proc_jend
              i = ii - PARENT_NODE%proc_istart
              j = jj - PARENT_NODE%proc_jstart
              nrank = nrank + 1
              ranklist(nrank)       = i      * PARENT_NODE%y_nnod + j
              ranklist_group(nrank) = (ii-1) * y_nnod             + (jj-1)
            enddo
            enddo
          else
            do jj=node%proc_jstart, node%proc_jend
            do ii=node%proc_istart, node%proc_iend
              i = ii - PARENT_NODE%proc_istart
              j = jj - PARENT_NODE%proc_jstart
              nrank = nrank + 1
              ranklist(nrank)       = j      * PARENT_NODE%x_nnod + i
              ranklist_group(nrank) = (jj-1) * x_nnod             + (ii-1)
            enddo
            enddo
          endif

          ! check number of process
          if( nrank.gt.1 ) then

            ! create new group 
            call MPI_Group_incl(FS_GROUP, nrank, ranklist_group, node%MERGE_GROUP, ierr)
            allocate(node%group_processranklist(nrank))
            do i=1,nrank
               ranklist_group(i)=i-1
            enddo
            call MPI_GROUP_TRANSLATE_RANKS(&
                 node%MERGE_GROUP, &
                 nrank, &
                 ranklist_group, &
                 FS_GROUP, &
                 node%group_processranklist, &
                 ierr)

            ! create X,Y group
            if( node%MERGE_GROUP.ne.MPI_UNDEFINED) then

              ! ------
              node%nnod   = ni * nj
              node%x_nnod = ni
              node%y_nnod = nj
              call MPI_Group_rank( &
                   node%MERGE_GROUP,node%inod,ierr)

              node%inod = node%inod + 1
              if (order == 'R') then
                 !     row-major
                 node%x_inod =    (node%inod-1)/node%y_nnod +1
                 node%y_inod = mod(node%inod-1, node%y_nnod)+1
              else
                 !     column-major
                 node%x_inod = mod(node%inod-1 , node%x_nnod)+1
                 node%y_inod =    (node%inod-1)/node%x_nnod +1
              end if
              call FS_create_mergeXY_group(node)

            endif

          endif

          ! deallocate
          deallocate( ranklist )
          deallocate( ranklist_group )

        end do !n

        ! subtree
        do n=1,2
          node => PARENT_NODE%sub_bt_node(n)
          if( FS_node_included(node) ) then
            call FS_create_merge_comm_recursive(node)
          endif
        enddo

        return
  end subroutine FS_create_merge_comm_recursive

  !---------------------------------------------------------------------------------
  !> subroutine FS_dividing_free
  !> @brief deallocate tree information
  !>
  !> @param[in,out] TREE   (input) type(bt_node) @n
  !>                       On entry, top node pointer
  !>                       On exit, destroyed.
  !>
  recursive subroutine FS_dividing_free(TREE)
        implicit none
        type(bt_node), intent(inout) :: TREE

        type(bt_node), pointer :: NODE
        integer :: ierr

        ! 子供を解放
        if( associated(TREE%sub_bt_node) ) then
          call FS_dividing_free(TREE%sub_bt_node(1))
          call FS_dividing_free(TREE%sub_bt_node(2))
        endif

        ! 自身の解放
        if( TREE%MERGE_GROUP_X.ne.MPI_UNDEFINED ) then
          call MPI_Group_free(TREE%MERGE_GROUP_X,ierr)
          deallocate(TREE%group_X_processranklist)
        endif
        if( TREE%MERGE_GROUP_Y.ne.MPI_UNDEFINED ) then
          call MPI_Group_free(TREE%MERGE_GROUP_Y,ierr)
          deallocate(TREE%group_Y_processranklist)
        endif
        if( (TREE%MERGE_GROUP.ne.FS_GROUP) .and. (TREE%MERGE_GROUP.ne.MPI_UNDEFINED) ) then
          call MPI_Group_free(TREE%MERGE_GROUP,ierr)
          deallocate(TREE%group_processranklist)
        endif
        if( associated(TREE%sub_bt_node) ) then
          deallocate(TREE%sub_bt_node)
          TREE%sub_bt_node => null()
        endif
        TREE%MERGE_GROUP   = MPI_UNDEFINED
        TREE%MERGE_GROUP_X = MPI_UNDEFINED
        TREE%MERGE_GROUP_Y = MPI_UNDEFINED

        TREE%parent_node => null()
        if( TREE%layer.eq.1 ) then
          if( associated( TREE%procs_i ) ) deallocate( TREE%procs_i )
          if( associated( TREE%procs_j ) ) deallocate( TREE%procs_j )
          TREE%procs_i => null()
          TREE%procs_j => null()
        endif

        return
  end subroutine FS_dividing_free

  !---------------------------------------------------------------------------------
  !> マージブロックのNを取得
  !> 全体次数Nが割り切れないとき、拡張したNの範囲で取得する
  !> function FS_get_N
  !> @brief get matrix size of merge block
  !>
  !> @param[in]     node   (input) type(bt_node) @n
  !>                       node pointer of merge block
  !>
  !> @return matrix size of merge block
  !>
  integer function FS_get_N(node) result(N)
        implicit none
        type(bt_node), intent(in) :: node

        N = node%nend - node%nstart + 1

        return
  end function FS_get_N

  !---------------------------------------------------------------------------------
  !> マージブロックのNを取得
  !> 全体次数Nが割り切れないとき、本来の次数Nの範囲で取得する
  !> function FS_get_N_active
  !> @brief get matrix size of merge block
  !>
  !> @param[in]     node   (input) type(bt_node) @n
  !>                       node pointer of merge block
  !>
  !> @return matrix size of merge block
  !>
  integer function FS_get_N_active(node) result(N)
        implicit none
        type(bt_node), intent(in) :: node

        N = node%nend_active - node%nstart + 1

        return
  end function FS_get_N_active

  !---------------------------------------------------------------------------------
  !> マージブロック内の行/列ブロック数を取得
  !> function FS_get_NBLK
  !> @brief get number of row/col block
  !>
  !> @param[in]     node   (input) type(bt_node) @n
  !>                       node pointer of merge block
  !>
  !> @return number of row/col block
  !>
  integer function FS_get_NBLK(node) result(NBLK)
        implicit none
        type(bt_node), intent(in) :: node

        NBLK = node%block_end - node%block_start + 1

        return
  end function FS_get_NBLK

  !---------------------------------------------------------------------------------
  ! マージブロックの1ブロックの行/列次数
  !> function FS_get_NB
  !> @brief get matrix size of one block
  !>
  !> @param[in]     node   (input) type(bt_node) @n
  !>                       node pointer of merge block
  !>
  !> @return matrix size of one block
  !>
  integer function FS_get_NB(node) result(NB)
        implicit none
        type(bt_node), intent(in) :: node

        NB = FS_get_N(node) / FS_get_NBLK(node)

        return
  end function FS_get_NB

  !---------------------------------------------------------------------------------
  ! マージブロックにおける自プロセス担当のQの全体先頭インデクスを取得
  !> subroutine FS_get_QTOP
  !> @brief get top index of Q in merge block
  !>
  !> @param[in]     node   (input) type(bt_node) @n
  !>                       node pointer of merge block
  !>
  !> @param[out]    IPQ    (output) INTEGER @n
  !>                       top index of 1st dimension.
  !>
  !> @param[out]    JPQ    (output) INTEGER @n
  !>                       top index of 2nd dimension.
  !>
  subroutine FS_get_QTOP(node, IPQ, JPQ)

        implicit none
        type(bt_node), intent(in), target :: node
        integer, intent(out) :: IPQ, JPQ

        integer :: NPROW, NPCOL, MYROW, MYCOL, N, NB
        integer :: I, J, II, JJ, ROW, COL
        type(bt_node), pointer :: ROOTNODE

        ! プロセス情報取得
        call FS_GRIDINFO(node, NPROW, NPCOL, MYROW, MYCOL)
        N  = FS_get_N(node)
        NB = FS_get_NB(node)

        ! 自プロセス行を検索
        ! この時点でIIはマージブロック内のグローバルインデクスになる
        DO I=1,N,NB
          call FS_INFOG1L('R',I,node,II,ROW)
          IF( ROW.EQ.MYROW ) THEN
            II = I
            EXIT
          END IF
        END DO

        ! 自プロセス列を検索
        ! この時点でJJはマージブロック内のグローバルインデクスになる
        DO J=1,N,NB
          call FS_INFOG1L('C',J,node,JJ,COL)
          IF( COL.EQ.MYCOL ) THEN
            JJ = J
            EXIT
          END IF
        END DO

        ! 全体インデクスに変換
        II = II + node%nstart - 1
        JJ = JJ + node%nstart - 1

        ! ルートを検索
        ROOTNODE => node
        DO WHILE( associated(ROOTNODE%parent_node) )
          ROOTNODE => ROOTNODE%parent_node
        END DO

        ! ルートブロックのローカルインデクスに変換
        IPQ = FS_INDXG2L('R', II, ROOTNODE)
        JPQ = FS_INDXG2L('C', JJ, ROOTNODE)

        return
  end subroutine FS_get_QTOP

  !---------------------------------------------------------------------------------
  ! マージブロック内のプロセス情報を取得
  !> subroutine FS_GRIDINFO
  !> @brief get process grid information
  !>
  !> @param[in]     node   (input) type(bt_node) @n
  !>                       node pointer of merge block
  !>
  !> @param[out]    NPROW  (output) INTEGER @n
  !>                       number of process grid row
  !>
  !> @param[out]    NPCOL  (output) INTEGER @n
  !>                       number of process grid column
  !>
  !> @param[out]    MYROW  (output) INTEGER @n
  !>                       row process index of own process (>=0)
  !>
  !> @param[out]    MYCOL  (output) INTEGER @n
  !>                       column process index of own process (>=0)
  !>
  subroutine FS_GRIDINFO(node, NPROW, NPCOL, MYROW, MYCOL)
        implicit none
        type(bt_node), intent(in) :: node
        integer, intent(out) :: NPROW, NPCOL, MYROW, MYCOL

        NPROW = node%x_nnod
        NPCOL = node%y_nnod
        MYROW = node%x_inod - 1
        MYCOL = node%y_inod - 1

        return
  end subroutine FS_GRIDINFO

  !---------------------------------------------------------------------------------
  ! 自プロセスに含まれないときでもLINDXにはROCSRCにおけるローカルインデクスが格納される
  ! ROCSRCにはCOMM_X/Yにおけるランク番号(0～)が入る
  !> subroutine FS_INFOG1L
  !> @brief convert index global to local
  !>
  !> @param[in]     COMP   (input) character @n
  !>                       set flag. 'R':row, 'C':columun
  !>
  !> @param[in]     GINDX  (input) INTEGER @n
  !>                       global index
  !>
  !> @param[in]     node   (input) type(bt_node) @n
  !>                       node pointer of merge block
  !>
  !> @param[out]    LINDX  (output) INTEGER @n
  !>                       local index
  !>
  !> @param[out]    ROCSRC (output) INTEGER @n
  !>                       row/column index of process grid include GINDX 
  !>
  subroutine FS_INFOG1L(COMP, GINDX, node, LINDX, ROCSRC)
        implicit none
        character, intent(in)  :: COMP
        integer  , intent(in)  :: GINDX
        integer  , intent(out) :: LINDX, ROCSRC
        type(bt_node), intent(in) :: node

        integer :: NB, IBLK, IBIT0, IBIT1, LBLK
        integer :: i
        logical, external :: LSAME

        ! 1ブロックの次数
        NB = FS_get_NB(node)

        ! グローバルインデクスが該当するブロック位置(0～)
        ! GINDXはマージブロック内での相対インデクス
        IBLK = (GINDX-1)/NB

        ! IBLKからdの[0]の桁だけ取り出し = IBIT0
        ! IBLKからdの[1]の桁だけ取り出し = IBIT1
        IBIT0 = 0
        IBIT1 = 0
        if( node%div_nbit .gt. 0 ) then
          do i=node%div_nbit,1,-1
            if( .not. BTEST(node%div_bit, i-1) ) then
              IBIT0 = ISHFT(IBIT0, 1)
              if( BTEST(IBLK, i-1) ) then
                IBIT0 = IBSET(IBIT0, 0)
              endif
            else
              IBIT1 = ISHFT(IBIT1, 1)
              if( BTEST(IBLK, i-1) ) then
                IBIT1 = IBSET(IBIT1, 0)
              endif
            endif
          enddo
        endif

        ! 行/列の場合分け
        if( LSAME(COMP,'R') ) then
          ! IBLKからdの[0]の桁だけ取り出し = プロセス行番号(0～)
          ! IBLKからdの[1]の桁だけ取り出し = プロセス行内のブロック位置(0～)
          ROCSRC = IBIT0
          LBLK   = IBIT1
        else
          ! IBLKからdの[1]の桁だけ取り出し = プロセス行番号(0～)
          ! IBLKからdの[0]の桁だけ取り出し = プロセス行内のブロック位置(0～)
          ROCSRC = IBIT1
          LBLK   = IBIT0
        endif

        ! ローカルインデクス
        LINDX = LBLK*NB + MOD(GINDX-1,NB) + 1

        return
  end subroutine FS_INFOG1L

  !---------------------------------------------------------------------------------
  !> subroutine FS_INFOG2L
  !> @brief convert index global to local
  !>
  !> @param[in]     GRINDX (input) INTEGER @n
  !>                       global row index
  !>
  !> @param[in]     GCINDX (input) INTEGER @n
  !>                       global column index
  !>
  !> @param[in]     node   (input) type(bt_node) @n
  !>                       node pointer of merge block
  !>
  !> @param[out]    LRINDX (output) INTEGER @n
  !>                       local row index
  !>
  !> @param[out]    LCINDX (output) INTEGER @n
  !>                       local column index
  !>
  !> @param[out]    RSRC   (output) INTEGER @n
  !>                       row index of process grid include GRINDX 
  !>
  !> @param[out]    CSRC   (output) INTEGER @n
  !>                       column index of process grid include GCINDX 
  !>
  subroutine FS_INFOG2L(GRINDX, GCINDX, node, LRINDX, LCINDX, RSRC, CSRC)
        implicit none
        integer, intent(in)  :: GRINDX, GCINDX
        integer, intent(out) :: LRINDX, LCINDX, RSRC, CSRC
        type(bt_node), intent(in) :: node

        call FS_INFOG1L('R', GRINDX, node, LRINDX, RSRC)
        call FS_INFOG1L('I', GCINDX, node, LCINDX, CSRC)

        return
  end subroutine FS_INFOG2L

  !---------------------------------------------------------------------------------
  !> function FS_INDXG2L
  !> @brief convert index global to local
  !>
  !> @param[in]     COMP   (input) character @n
  !>                       set flag. 'R':row, 'C':columun
  !>
  !> @param[in]     GINDX  (input) INTEGER @n
  !>                       global index
  !>
  !> @param[in]     node   (input) type(bt_node) @n
  !>                       node pointer of merge block
  !>
  !> @return local index
  !>
  integer function FS_INDXG2L(COMP, GINDX, node) result(LINDX)
        implicit none
        character, intent(in)  :: COMP
        integer  , intent(in)  :: GINDX
        type(bt_node), intent(in) :: node

        integer :: ROCSRC

        call FS_INFOG1L(COMP, GINDX, node, LINDX, ROCSRC)

        return
  end function FS_INDXG2L

  !---------------------------------------------------------------------------------
  !> function FS_INDXL2G
  !> @brief convert index local to global
  !>
  !> @param[in]     COMP   (input) character @n
  !>                       set flag. 'R':row, 'C':columun
  !>
  !> @param[in]     LINDX  (input) INTEGER @n
  !>                       local index
  !>
  !> @param[in]     MYROC  (input) INTEGER @n
  !>                       row/column index of process grid include LINDX
  !>
  !> @param[in]     node   (input) type(bt_node) @n
  !>                       node pointer of merge block
  !>
  !> @return global index
  !>
  integer function FS_INDXL2G(COMP, LINDX, MYROC, node) result(GINDX)
        implicit none
        character, intent(in) :: COMP
        integer  , intent(in) :: LINDX
        integer  , intent(in) :: MYROC
        type(bt_node), intent(in) :: node

        integer :: NB, LBLK, IBLK, IBIT0, IBIT1, IPNT0, IPNT1
        integer :: i
        logical, external :: LSAME

        ! 1ブロックの次数
        NB = FS_get_NB(node)

        ! ローカルインデクスが該当するブロック位置(0～)
        ! GRINDX,GCINDXはマージブロック内での相対インデクス
        LBLK = (LINDX-1)/NB

        ! 行/列の場合分け
        if( LSAME(COMP,'R') ) then
          ! IBLKからdの[0]の桁だけ取り出し = プロセス行番号(0～)
          ! IBLKからdの[1]の桁だけ取り出し = プロセス行内のブロック位置(0～)
          IBIT0 = MYROC
          IBIT1 = LBLK
        else
          ! IBLKからdの[1]の桁だけ取り出し = プロセス行番号(0～)
          ! IBLKからdの[0]の桁だけ取り出し = プロセス行内のブロック位置(0～)
          IBIT1 = MYROC
          IBIT0 = LBLK
        endif

        ! 全体でのブロック位置を取得
        IBLK = 0
        IPNT0 = 0
        IPNT1 = 0
        if( node%div_nbit .gt. 0 ) then
          do i=1,node%div_nbit
            if( .not. BTEST(node%div_bit, i-1) ) then
              if( BTEST(IBIT0, IPNT0) ) then
                IBLK = IBSET(IBLK, i-1)
              endif
              IPNT0 = IPNT0 + 1
            else
              if( BTEST(IBIT1, IPNT1) ) then
                IBLK  = IBSET(IBLK, i-1)
              endif
              IPNT1 = IPNT1 + 1
            endif
          enddo
        endif

        ! グローバルインデクス
        GINDX = IBLK * NB + mod(LINDX-1, NB) + 1

        return
  end function FS_INDXL2G

  !---------------------------------------------------------------------------------
  !> subroutine print_tree
  !> @brief output log of tree information recursive
  !>
  !> @param[in]     kout  (input) INTEGER @n
  !>                      file No.
  !>
  !> @param[in]     tree  (input) type(bt_node) @n
  !>                      tree pointer
  !> 
  recursive subroutine print_tree(kout, tree)
        implicit none
        integer, intent(in) :: kout
        type(bt_node), target, intent(in) :: tree
        integer :: nnod,x_nnod,y_nnod 
        integer :: inod,x_inod,y_inod
        integer :: i

        call FS_get_procs(nnod,x_nnod,y_nnod)
        call FS_get_id(inod,x_inod,y_inod)

        if( tree%layer.eq.1 ) then
          write(kout,*) "nnod,(x_nnod,y_nnod)=",nnod,"(",x_nnod,y_nnod,")"
          write(kout,*) "inod,(x_inod,y_inod)=",inod,"(",x_inod,y_inod,")"
        endif

        call print_node(kout, tree)

        if( associated(tree%sub_bt_node) ) then
          do i=1,size(tree%sub_bt_node)
            call print_tree(kout, tree%sub_bt_node(i))
          enddo
        endif

        return
  end subroutine print_tree

  !---------------------------------------------------------------------------------
  !> subroutine print_node
  !> @brief output log of node information
  !>
  !> @param[in]     kout  (input) INTEGER @n
  !>                      file No.
  !>
  !> @param[in]     node  (input) type(bt_node) @n
  !>                      node pointer
  !> 
  subroutine print_node(kout, node)
        implicit none
        integer, intent(in) :: kout
        type(bt_node), target, intent(in) :: node

        write(kout,'(a)') "******************************"
        write(kout,*) "layer               =",node%layer
        write(kout,*) "direction_horizontal=",node%direction_horizontal
        write(kout,*) "nstart              =",node%nstart
        write(kout,*) "nend                =",node%nend
        write(kout,*) "nend_active         =",node%nend_active
        write(kout,*) "proc_istart         =",node%proc_istart
        write(kout,*) "proc_iend           =",node%proc_iend
        write(kout,*) "proc_jstart         =",node%proc_jstart
        write(kout,*) "proc_jend           =",node%proc_jend
        write(kout,*) "block_start         =",node%block_start
        write(kout,*) "block_end           =",node%block_end
        write(kout,'(1x,a,10000i3)') "procs_i =",node%procs_i
        write(kout,'(1x,a,10000i3)') "procs_j =",node%procs_j
        write(kout,*) "parent_node =",associated(node%parent_node)
        write(kout,*) "child_node  =",associated(node%sub_bt_node)
        write(kout,*) "merge procs    = ", node%nnod
        write(kout,*) "merge procs X  = ", node%x_nnod
        write(kout,*) "merge procs Y  = ", node%y_nnod
        write(kout,*) "merge rankid   = ", node%inod
        write(kout,*) "merge rankid X = ", node%x_inod
        write(kout,*) "merge rankid Y = ", node%y_inod
        call bitprint(kout, "bit stream     = ", node%div_bit, node%div_nbit)
        write(kout,*) "#dights of bit = ", node%div_nbit

        return
  end subroutine print_node

  subroutine bitprint(kout, title, ibit, nbit)
    implicit none
    integer, intent(in)       :: kout, ibit, nbit
    character(*), intent(in)  :: title
    character(nbit) :: cbit
    integer :: i

    cbit = ""
    do i=nbit,1,-1
      if( BTEST(ibit,i-1) ) then
        cbit = trim(cbit) // "1"
      else
        cbit = trim(cbit) // "0"
      endif
    enddo

!    write(kout,'(2(1x,a),1x,i)') trim(title), trim(cbit), ibit
    write(kout,*) trim(title), trim(cbit), ibit

    return
  end subroutine bitprint

end module FS_dividing_mod
