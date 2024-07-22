       subroutine  mat_set_h( n, a, nm )
       use MPI
       use eigen_libs_mod
       implicit NONE

       integer, intent(IN)    :: n, nm
       complex(8), intent(OUT):: a(1:nm,*)

       integer                :: COMM, x_COMM, y_COMM
       integer                :: nnod, x_nnod, y_nnod
       integer                :: inod, x_inod, y_inod
       integer                :: i, i_1, i_2, i_3
       integer                :: j, j_1, j_2, j_3

       complex(8)             :: a_ij, a_ji
       real(8)                :: re,im
       integer                :: seedsize
       integer, allocatable   :: seed(:)
       integer                :: DESCA(9)
       complex(8), allocatable:: as(:,:)
       integer                :: ICTXT, INFO, ny

          call  eigen_get_comm ( COMM, x_COMM, y_COMM )
          call  eigen_get_procs( nnod, x_nnod, y_nnod )
          call  eigen_get_id   ( inod, x_inod, y_inod )

          ICTXT = eigen_get_blacs_context( )
          call DESCINIT(DESCA, n, n, 1, 1, 0, 0, ICTXT, nm, INFO)

          i_2 = eigen_loop_start( 1, 'Y' )
          i_3 = eigen_loop_end  ( n, 'Y' )
          j_2 = eigen_loop_start( 1, 'X' )
          j_3 = eigen_loop_end  ( n, 'X' )

          call  random_seed( size=seedsize )
          allocate( seed(seedsize) )
          seed(:) = inod
          call  random_seed( put=seed )

          ny = (n-1)/y_nnod + 1
          allocate( as(1:nm, 1:ny) )
          a (1:nm, 1:ny) = DCMPLX(0d0, 0d0)
          as(1:nm, 1:ny) = DCMPLX(0d0, 0d0)

          do i_1 = i_2, i_3
             do j_1 = j_2, j_3
                i = y_inod + y_nnod * (i_1 - 1)
                j = x_inod + x_nnod * (j_1 - 1)

                call random_number(re)
                call random_number(im)
                re = (re - 0.5d0) * 2.0d0 / 2
                im = (im - 0.5d0) * 2.0d0 / 2
                if ( i == j ) then
                   a_ji = dcmplx(re,0.0d0)
                else
                   a_ji = dcmplx(re, im)
                end if

                a (j_1, i_1) = a_ji
                as(j_1, i_1) = a_ji

             end do
          end do
          call PZTRANC(n, n,
     &       1D0, as, 1, 1, DESCA, 1D0, a, 1, 1, DESCA)

          deallocate( as )
          deallocate( seed )

       end subroutine  mat_set_h

