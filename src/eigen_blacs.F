!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_blacs_mod
!
! Purpose
! =======
!
! eigen_blacs_mod manages the modules for BLACS.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      module eigen_blacs_mod

      use mpi
      use eigen_devel_mod

      implicit none
      private

      public :: eigen_blacs_init
      public :: eigen_blacs_exit
      public :: eigen_get_blacs_context

      integer, private :: BLACS_ICONTXT_FOR_EIGENEXA

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_blacs_init
!
! Purpose
! =======
!
! eigen_blocs_init initializes and sets up the grid mapping.
!
!
! Arguments
! =========
!
! TRD_COMM_WORLD (input) integer
!         communicator
!
! x_nnod  (input) integer
!         dimension of x-dir
!
! y_nnod  (input) integer
!         dimension of y-dir
!
! GRID_major (input) character(*)
!         'C': column major
!         'R': row major
!     
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_blacs_init(TRD_COMM_WORLD,
     &     x_nnod, y_nnod,
     &     GRID_major)

      integer,      intent(in) ::  TRD_COMM_WORLD
      integer,      intent(in) ::  x_nnod, y_nnod
      character(*), intent(in) ::  GRID_major

      integer, allocatable ::  tmpgrid(:,:)
      integer, allocatable ::  kk0(:), kk1(:)
      integer          ::  group0, group1
      integer          ::  i, j, k, ierr


      call BLACS_GET(0, 0, BLACS_ICONTXT_FOR_EIGENEXA)

      allocate(tmpgrid(1:x_nnod, 1:y_nnod),
     &         kk0(1:x_nnod), kk1(1:x_nnod), stat=ierr)
      if (ierr /= 0) then
        call eigen_abort( "Memory allocation error [init].", 1 )
        return
      end if

      call MPI_Comm_group(MPI_COMM_WORLD, group0, ierr)
      call MPI_Comm_group(TRD_COMM_WORLD, group1, ierr)

      if (Grid_major == 'R') then
        do j = 1, y_nnod
          do i = 1, x_nnod
            kk1(i) = j+(i-1)*y_nnod-1
          end do
          k = x_nnod
          call MPI_Group_translate_ranks(group1, k, kk1,
     &         group0, kk0, ierr)
          do i = 1, x_nnod
            tmpgrid(i, j) = kk0(i)
          end do
        end do
      else
        do j = 1, y_nnod
          do i = 1, x_nnod
            kk1(i) = i+(j-1)*x_nnod-1
          end do
          k = x_nnod
          call MPI_Group_translate_ranks(group1, k, kk1,
     &         group0, kk0, ierr)
          do i = 1, x_nnod
            tmpgrid(i, j) = kk0(i)
          end do
        end do
      end if

      call MPI_Barrier(TRD_COMM_WORLD, ierr)

      call BLACS_GRIDMAP(BLACS_ICONTXT_FOR_EIGENEXA,
     &     tmpgrid, x_nnod, x_nnod, y_nnod)

      call MPI_Barrier(TRD_COMM_WORLD, ierr)

      call MPI_Group_free(group0, ierr)
      call MPI_Group_free(group1, ierr)

      deallocate(tmpgrid)
      deallocate(kk0)
      deallocate(kk1)

      return

      end subroutine eigen_blacs_init

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_blacs_exit
!
! Purpose
! =======
!
! eigen_blacs_exit finalizes the utility of the BLACS library
!
!
! Arguments
! =========
!     
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_blacs_exit()

      call BLACS_GRIDEXIT(BLACS_ICONTXT_FOR_EIGENEXA)
!     commented out on ver 2.2a
!     call BLACS_EXIT(1)

      return

      end subroutine eigen_blacs_exit

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: eigen_get_blacs_context
!
! Purpose
! =======
!
! eigen_get_blacs_context returns available BLACS context.
!
!
! Arguments
! =========
!
!
!--------*---------*---------*---------*---------*---------*---------*-*

      integer function eigen_get_blacs_context()


      eigen_get_blacs_context = BLACS_ICONTXT_FOR_EIGENEXA

      return

      end function eigen_get_blacs_context

      end module eigen_blacs_mod
