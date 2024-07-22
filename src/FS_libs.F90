!>
!> @file   FS_libs.F90
!> @brief  module FS_libs_mod
!!>
!
!
!
!> module FS_libs_mod
!> @brief  @n
!> Purpose @n
!> ======= @n
!> FS_libs_mod include common variables and subroutines
!>
module FS_libs_mod
  use mpi
  implicit none

  ! comm world
  integer :: FS_COMM_WORLD = MPI_COMM_WORLD  !< comm world
  integer :: FS_MYRANK = 0                   !< rank No.
  logical :: FS_COMM_MEMBER = .FALSE.        !< comm member fla

  integer :: FS_GROUP = MPI_UNDEFINED !< FS_COMM_WORLD group

  type  process_grid 
     integer :: nnod,x_nnod,y_nnod
     integer :: inod,x_inod,y_inod
  end type process_grid
  type(process_grid) :: FS_node

  !> type version_t
  !> @brief version info
  type, public :: version_t
     integer       :: Major_Version !< Major version
     integer       :: Minor_Version !< Minor version
     integer       :: Patch_Level   !< Patchlevel 0=none,1=a,2=b,...
     character(32) :: date          !< Release date
     character(32) :: vcode         !< Version code name
  end type  version_t

  !> version info
  type(version_t) :: FS_Version & ! 1.0
       = version_t (            &
       1, 1, 0,                 & ! Major, Minor, Patchlevel
       'Mar 31, 2019',          & ! Release date
       'FS proto'               & ! Version code
       )

  !> grid major
  character(1) :: FS_GRID_major = 'C'

#ifdef WRITE_INPUT_VEC
  !> matrix type
  integer :: mat_type
#endif

  interface

       subroutine eigen_FS(n, nvec, a, lda, w, z, ldz, &
          m_forward, m_backward, mode)
       integer,   intent(in)           :: n
       integer,   intent(in)           :: nvec
       real(8),   intent(inout)        :: a(lda,*)
       integer,   intent(in)           :: lda
       real(8),   intent(out)          :: w(1:n)
       real(8),   intent(out)          :: z(ldz,*)
       integer,   intent(in)           :: ldz
       integer,   intent(in), optional :: m_forward
       integer,   intent(in), optional :: m_backward
       character(*), intent(in), optional :: mode
       end subroutine eigen_FS

  end interface

contains

  !--------*---------*---------*---------*---------*---------*---------*-*
  !> subroutine FS_get_version
  !> @brief get version information
  !> @param[out] version   version
  !> @param[out] date      date
  !> @param[out] vcode     vcode
  subroutine FS_get_version(version, date, vcode)

    integer,      intent(out)           :: version
    character(*), intent(out), optional :: date
    character(*), intent(out), optional :: vcode


    version = FS_Version%Major_Version * 100 &
         + FS_Version%Minor_Version * 10  &
         + FS_Version%Patch_Level

    if (present(date)) then
       date = FS_Version%date
    end if

    if (present(vcode)) then
       vcode = FS_Version%vcode
    end if

    return

  end subroutine FS_get_version

  !--------*---------*---------*---------*---------*---------*---------*-*
  !> subroutine FS_show_version
  !> @brief print version information
  subroutine FS_show_version()

    use eigen_libs0_mod,only : eigen_get_id
    character(256) :: version
    character(1  ) :: patchlevel
    integer         :: i
    integer         :: id, x_id, y_id

    call eigen_get_id(id, x_id, y_id)

    i = min(26, FS_Version%Patch_Level) + 1
    patchlevel = (" abcdefghijklmnopqrstuvwxyz*" (i:i))

    write(version, '(I1,A,I1,A)')     &
         FS_Version%Major_Version, &
         '.',FS_Version%Minor_Version, trim(patchlevel)

    if (id == 1) then
       print*, "## FS version (", trim(version), &
            ") / (", trim(FS_Version%date),         &
            ") / (", trim(FS_Version%vcode), ")"
    end if

    return

  end subroutine FS_show_version

  !--------*---------*---------*---------*---------*---------*---------*-*
  ! init
  !--------*---------*---------*---------*---------*---------*---------*-*
  !> subroutine FS_init
  !> @brief initialize routine of FS
  !> @param[in]  order  process grid order. 'R':row-major, 'C':column-major
  subroutine FS_init(comm,order)
    use eigen_libs0_mod
    implicit none

    integer, intent(in), optional :: comm
    character(*), intent(in), optional :: order

    integer :: comm0
!    integer :: x_comm, y_comm

    integer :: p,color
    integer :: eigen_comm,eigen_x_comm, eigen_y_comm
    integer :: nnod, x_nnod, y_nnod
    integer :: inod, x_inod, y_inod
    integer :: ierr

    if (present(comm)) then
       comm0=comm
    else
       comm0=MPI_COMM_WORLD
    endif

    if (present(order)) then
       FS_GRID_major = order(1:1)
    else
       FS_GRID_major = 'C'
    endif
    if (FS_GRID_major == 'R' .or. FS_GRID_major == 'r') then
       FS_GRID_major = 'R'
    else
       FS_GRID_major = 'C'
    end if

    !      call eigen_init(order=FS_GRID_major,comm=comm0)
    call eigen_init0(order=FS_GRID_major,comm=comm0)
    call eigen_get_comm(eigen_comm, eigen_x_comm, eigen_y_comm)
    ! -----
    ! FS_COMMM_WORLDの設定
    call eigen_get_procs(nnod, x_nnod, y_nnod)
    call eigen_get_id   (inod, x_inod, y_inod)

    p=INT(log(dble(nnod))/log(2.0d0))

    if(inod<=2**p)then
       color = 0
       FS_COMM_MEMBER = .TRUE.
    else
       color = 1
       FS_COMM_MEMBER = .FALSE.
    endif
    call MPI_Comm_split(eigen_comm, color, inod, FS_COMM_WORLD, ierr)

    if(FS_COMM_MEMBER)then
       call MPI_COMM_RANK(FS_COMM_WORLD, FS_MYRANK, ierr)
       call MPI_Comm_group(FS_COMM_WORLD, FS_GROUP, ierr)

       call MPI_Comm_size(FS_COMM_WORLD, nnod, ierr)
       call FS_init_cartesian(FS_GRID_major,nnod,FS_MYRANK+1)
    else
       FS_MYRANK      = -1
       FS_node%nnod   = -1
       FS_node%x_nnod = -1
       FS_node%y_nnod = -1
       FS_node%inod   = -1
       FS_node%x_inod = -1
       FS_node%y_inod = -1
    endif

    return
  end subroutine FS_init

  subroutine FS_init_cartesian(GRID_major,nnod,inod)
    implicit none 
    character(1) :: GRID_major
    integer :: nnod,inod ! 引数

    integer :: i,k
    integer :: x_nnod,y_nnod
    integer :: x_inod,y_inod

    !        if (TRD_COMM_WORLD /= MPI_COMM_NULL) then
    !---- Setup 2D process map ---
    x_nnod = int(sqrt(dble(nnod)))
    i = 1                 ! minimum factor, x_nnod must be
    !     multiple of k
    if (mod(nnod, i) == 0) then
       k = i
    else
       k = 1
    end if
    do
       if (x_nnod <= k) exit
       if (mod(x_nnod, k) == 0 .and. &
            mod(nnod, x_nnod) == 0) exit
       x_nnod = x_nnod-1
    end do                !!
    y_nnod = nnod/x_nnod

    if (GRID_major == 'R') then
       !     row-major
       x_inod =    (inod-1)/y_nnod +1
       y_inod = mod(inod-1, y_nnod)+1
    else
       !     column-major
       x_inod = mod(inod-1, x_nnod)+1
       y_inod =    (inod-1)/x_nnod +1
    end if
    !        endif

    FS_node%nnod   = nnod
    FS_node%x_nnod = x_nnod
    FS_node%y_nnod = y_nnod
    FS_node%inod   = inod
    FS_node%x_inod = x_inod
    FS_node%y_inod = y_inod

    return

  end subroutine FS_init_cartesian


  !--------*---------*---------*---------*---------*---------*---------*-*
  ! free
  !--------*---------*---------*---------*---------*---------*---------*-*
  !> subroutine FS_free
  !> @brief free work
  subroutine FS_free()
    use eigen_libs0_mod,only : eigen_free0
    implicit none
    integer ierr
    ! eigen_libs free
    call eigen_free0()
    !      call MPI_Group_free(FS_GROUP,ierr)
    call MPI_Comm_free(FS_COMM_WORLD,ierr)

    return
  end subroutine FS_free

  !--------*---------*---------*---------*---------*---------*---------*-*
  ! calculate work array size
  !--------*---------*---------*---------*---------*---------*---------*-*
  !> subroutine FS_WorkSize
  !> @brief calculate work array size for FS_PDSTEDC
  !> @param[in]  N       The order of the tridiagonal matrix T.
  !> @param[out] LWORK   work array size for real
  !> @param[out] LIWORK  work array size for integer
  subroutine FS_WorkSize(N, LWORK, LIWORK)
    use FS_const_mod
    implicit none

    integer, intent(in)  :: N
    integer(8), intent(out) :: LWORK, LIWORK

    integer :: nnod, x_nnod, y_nnod
    integer :: NP, NQ

    call FS_get_procs(nnod, x_nnod, y_nnod)
    call FS_get_matdims(N, NP, NQ)

    LWORK  = 1 + 7*N + 3*INT(NP,8)*INT(NQ,8) + INT(NQ,8)*INT(NQ,8)
    LIWORK = 1 + 8*N + 2*4*y_nnod

    return
  end subroutine FS_WorkSize

  !--------*---------*---------*---------*---------*---------*---------*-*
  ! get procs
  !--------*---------*---------*---------*---------*---------*---------*-*
  !> subroutine FS_get_procs
  !> @brief get number of process and process grid size
  !> @param[out] nnod    number of process
  !> @param[out] x_nnod  number of row of process grid
  !> @param[out] y_nnod  number of column of process grid
  subroutine FS_get_procs(nnod, x_nnod, y_nnod)
    implicit none

    integer, intent(out) :: nnod, x_nnod, y_nnod

    !      call eigen_get_procs(nnod,x_nnod,y_nnod)
    nnod   = FS_node%nnod 
    x_nnod = FS_node%x_nnod 
    y_nnod = FS_node%y_nnod 
    return
  end subroutine FS_get_procs

  !--------*---------*---------*---------*---------*---------*---------*-*
  ! get id
  !--------*---------*---------*---------*---------*---------*---------*-*
  !> subroutine FS_get_id
  !> @brief get number of process and process grid size
  !> @param[out] inod    process No. (>=1)
  !> @param[out] x_inod  row index of process grid (>=1)
  !> @param[out] y_inod  column index of process grid (>=1)
  subroutine FS_get_id(inod, x_inod, y_inod)
    implicit none

    integer, intent(out) :: inod, x_inod, y_inod

    !      call eigen_get_id(inod,x_inod,y_inod)
    inod   = FS_node%inod
    x_inod = FS_node%x_inod
    y_inod = FS_node%y_inod

    return
  end subroutine FS_get_id

  !--------*---------*---------*---------*---------*---------*---------*-*
  ! get mat dims
  !--------*---------*---------*---------*---------*---------*---------*-*
  !> subroutine FS_get_matdims
  !> @brief get matrix size
  !> @param[in]  N   The order of the tridiagonal matrix T.
  !> @param[out] nx  row size (1st dimension)
  !> @param[out] ny  column size (2nd dimension)
  subroutine FS_get_matdims(n, nx, ny)
    implicit none
    !
    integer, intent(in)    :: n
    integer, intent(out)   :: nx, ny
    !
    integer :: nnod, x_nnod, y_nnod
    integer :: n1
    !
    call FS_get_procs(nnod, x_nnod, y_nnod)
    !
    n1 = n / nnod
    if( mod(n, nnod) .ne. 0 ) then
       n1 = n1 + 1
    endif
    nx = n1 * (nnod / x_nnod)
    ny = n1 * (nnod / y_nnod)

    return
  end subroutine FS_get_matdims

  !--------*---------*---------*---------*---------*---------*---------*-*
  ! get GRID_major
  !--------*---------*---------*---------*---------*---------*---------*-*
  !> subroutine FS_get_grid_major
  !> @brief get grid major
  !> @param[out] Major   process grid order. 'R':row-major, 'C':column-major
  subroutine FS_get_grid_major(Major)
    implicit none

    character(*), intent(out) :: Major

    Major(1:1) = FS_GRID_major(1:1)
    return

  end subroutine FS_get_grid_major


end module FS_libs_mod
