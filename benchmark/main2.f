!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Program: EigenExa_benchmark
!
! Purpose
! =======
!
! Benchmarking and verification of the EigenExa solver kernels
!
!
! Copyright(C) 2012-2021 RIKEN.
! Copyright(C) 2011-2012 Toshiyuki Imamura
!                        Graduate School of Informatics and Engineering,
!                        The University of Electro-Communications.
! Copyright (C) 2011- 2015 Japan Atomic Energy Agency.
! 
! Redistribution  and  use  in  source and binary forms, with or without
! modification,  are  permitted  provided  that the following conditions
! are met:
! 
! * Redistributions  of  source  code  must  retain  the above copyright
!   notice,  this  list  of  conditions  and  the  following  disclaimer.
! * Redistributions  in  binary  form must reproduce the above copyright
!   notice,  this list of conditions and the following disclaimer in the
!   documentation  and/or other materials provided with the distribution.
! 
! THIS  SOFTWARE  IS  PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! ``AS IS''  AND  ANY  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A  PARTICULAR  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL,  EXEMPLARY,  OR  CONSEQUENTIAL  DAMAGES  (INCLUDING,  BUT NOT
! LIMITED  TO,  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA,  OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY  OF  LIABILITY,  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF  THIS  SOFTWARE,  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      program EigenExa_benchmark

      use eigen_libs_mod
      use mpi
!$    use omp_lib

      implicit none

      real(8), parameter :: ONE = 1D0

      real(8), pointer   :: a(:,:),z(:,:),w(:)
      real(8)            :: PAI, EPS, EPS2, EPS4
      real(8)            :: d1, d2, flops, commu, time

      integer            :: nnod, x_nnod, y_nnod
      integer            :: inod, x_inod, y_inod
      integer            :: i, i_inod, i_nnod, itr, ierror, istat
      integer            :: n, nvec, m, mb, nall, mtype, msolver, merror
      integer            :: nm, nx, ny
      integer(8)         :: imem
      integer            :: new_comm, color, key, n_comms
      integer            :: new_dim_x, new_dim_y, dims(2), coords(2)
      logical            :: periods(2), reorder
      character*1        :: mode = ' ', grid = ' '
      character*256      :: input_file = ' '
      logical            :: check_accuracy = .true.
      character*1024     :: in_buffer
      integer            :: pbuffer(10)

      integer            :: nargc
      character*256      :: argv

      integer            :: version
      character*128      :: date, vcode
      integer iij
!-
!----------------------------------------------------------------------
!-
!      call MPI_Init_thread( MPI_THREAD_MULTIPLE, i, ierror )
!      call MPI_Init(ierror)
      call MPI_Init_thread( MPI_THREAD_SERIALIZED, i, ierror )
      call MPI_Comm_rank( MPI_COMM_WORLD, i_inod, ierror )
      call MPI_Comm_size( MPI_COMM_WORLD, i_nnod, ierror )
      if ( i_inod == 0 ) call system( "date" )
!-
!----------------------------------------------------------------------
!-
      nargc = command_argument_count()
      i = 0
      do
        i = i+1
        call get_command_argument(i, argv)
        if (len_trim(argv) == 0) exit
        argv = trim(argv)

        if (argv(1:2) == '-h') then
          if (i_inod == 0) then

            call print_help_message()

          end if
          goto 99999
        end if

        if (argv(1:2) == '-L') then
          if (i_inod == 0) then

            call print_header_message()

            print*,"The list of test matrices:"
            print*,""

            call print_mat_name_list()

            print*,""

          end if
          goto 99999
        end if

        if (argv(1:2) == '-c') then
          check_accuracy = .true.
          cycle
        end if

        if (argv(1:2) == '-n') then
          check_accuracy = .false.
          cycle
        end if

        if (argv(1:2) == '-f') then
          i = i+1
          call get_command_argument(i, argv)
          if (len_trim(argv) == 0) exit
          input_file = trim(argv)
          cycle
        end if

        if (argv(1:2) == '-g') then
          i = i+1
          call get_command_argument(i, argv)
          if (len_trim(argv) == 0) exit
          grid = argv(1:1)
          if (grid == 'x') grid = ' '
          cycle
        end if

        if (argv(1:2) == '-x') then
          i = i+1
          call get_command_argument(i, argv)
          if (len_trim(argv) == 0) exit
          read(argv, *) new_dim_x
          i = i+1
          call get_command_argument(i, argv)
          if (len_trim(argv) == 0) exit
          read(argv, *) new_dim_y
          if ( i_inod == 0 ) print*,new_dim_x, new_dim_y
          grid = 'x'
          cycle
        end if

        print*,"Command option '", argv(1:2), "' is inavailable. ",
     &       "This benchmark code skips it, and continues to proceed."

      end do
!     
!     User can specify the process grid shape by -g option
!     
      n_comms = 1
      color   = 0
      select case (grid)
      case ('R', 'r')
        call eigen_init(order='r')
        new_comm = MPI_COMM_WORLD
      case ('C', 'c')
        call eigen_init(order='c')
        new_comm = MPI_COMM_WORLD
      case ('A', 'a')
        n_comms = i_nnod
        color   = i_inod
        call eigen_init(MPI_COMM_SELF)
        new_comm = MPI_COMM_SELF
      case ('1', '2', '3', '4', '5', '6', '7', '8', '9')
        read(grid, *) n_comms
        color = mod(i_inod, n_comms)
        key   = i_inod / n_comms
        call MPI_Comm_split(MPI_COMM_WORLD,
     &       color, key,
     &       new_comm, ierror)
        if (color /= 0) new_comm = MPI_COMM_NULL
        call eigen_init(new_comm)
      case ('x')
        if (new_dim_x * new_dim_y /= i_nnod) then
          if (i_inod == 0) then
            print*,"Illegal dimensions are specified."
          end if
          call MPI_Abort(MPI_COMM_WORLD, 1, ierror)
        end if
        if (new_dim_x > new_dim_y) then
          if (i_inod == 0) then
            print*,"This process map is not supported."
            print*,"Px should be smaller than Py"
          end if
          call MPI_Abort(MPI_COMM_WORLD, 1, ierror)
        end if
        dims(1:2)    = (/ new_dim_x, new_dim_y /)
        periods(1:2) = .false.
        reorder      = .false.
        call MPI_Cart_create(MPI_COMM_WORLD,
     &       2, dims, periods, reorder,
     &       new_comm, ierror)
        call eigen_init(new_comm)
      case default
        call eigen_init()
        new_comm = MPI_COMM_WORLD
      end select

      call eigen_get_version(version, date=date ) ! , vcode=vcode)
!     
!-f option specifies the input file, default is 'IN'
!     
      if (input_file(1:1) == ' ') then
        input_file = 'IN'
      end if
!-
!----------------------------------------------------------------------
!-
      if (i_inod == 0) then
        open(unit=10, file=input_file, status='old', recl=1024,
     &       iostat=ierror)
        if (ierror /= 0) then
          print*,"Can not open input file [",trim(input_file),"]"
          call MPI_Abort(MPI_COMM_WORLD, 1, ierror)
        end if
        print*,"INPUT FILE='",trim(input_file),"'"
      end if
!-
!----------------------------------------------------------------------
!-
      call eigen_get_procs(nnod, x_nnod, y_nnod)
      call eigen_get_id   (inod, x_inod, y_inod)
!-
!----------------------------------------------------------------------
!-
      PAI  = get_constant_pai()
      EPS  = get_constant_eps()
      EPS2 = sqrt(EPS)
      EPS4 = sqrt(EPS2)
!-
!----------------------------------------------------------------------
!-
      do itr=1,1000000
        if (i_inod == 0) then
!     
!     Input file format
!     
!     N bx by mode matrix solver
!     
!     N      : matrix dimension
!     nvec   : the number of eigenvectors to be computed
!     bx     : block width for the forward transformation
!     by     : block width for the backward transformation
!     mode   : eigensolver mode { 0 : only eigenvalues }
!     { 1 : eigenvalues and corresponding eigenvectors}
!     { 2 : mode 1 + accuracy improvement for eigenvalue}
!     matrix : test matrix { 11 types, 0 ... 10 }
!     solver : { 0 : eigen_sx, new algorithm, faster on the K }
!     { 1 : eigen_s,  conventional algorithm }
!     
!     if a line starts from '!', the line is treated as a comment
!     
          do
            read(unit=10, fmt='(a1024)',
     &           err=10000, end=10000) in_buffer
            if (in_buffer(1:1) /= '!') EXIT
          end do
          backspace(10)
          read(10, fmt=*,
     &         err=10000, end=10000)
     &         n, nvec, m, mb, nall, mtype, msolver, merror
        end if
        goto 20000
10000   continue
        n = -1
20000   continue

!     n = 64
!     nvec = n
!     m = 32
!     mb = 128
!     nall = 1
!     mtype = 0
!     msolver = 0
!     merror = 0

        pbuffer(1) = n
        pbuffer(2) = nvec
        pbuffer(3) = m
        pbuffer(4) = mb
        pbuffer(5) = nall
        pbuffer(6) = mtype
        pbuffer(7) = msolver
        pbuffer(8) = merror

        call MPI_Bcast(pbuffer, 8, MPI_INTEGER,
     &       0, MPI_COMM_WORLD, ierror)

        n       = pbuffer(1)
        nvec    = pbuffer(2)
        m       = pbuffer(3)
        mb      = pbuffer(4)
        nall    = pbuffer(5)
        mtype   = pbuffer(6)
        msolver = pbuffer(7)
        merror  = pbuffer(8)

        if (n <= 0) exit
        if (n < nvec) then
#if DEBUG
          if (i_inod == 0) then
            print*,"Nvec (the number of eigenvectors to be computed) is"
            print*,"larger than N. It is adjusted equivalent to N."
          end if
#endif
          nvec = min(nvec, n)
        end if

        mode = 'A'
        select case (nall)
        case (0)
          mode = 'N'            ! only eigenvalues, no eigenvector
        case (1)
          mode = 'A'            ! all the eigenpairs
        case (2)
          mode = 'X'            ! mode 'A' + accuracy improvement
        case (3)
          mode = 'S'            ! skip DC but internal Z=identity
        case (4)
          mode = 'T'            ! run DC but skip TRBAK
        case (5)
          mode = 'C'            ! skip DC and TRBAK return X=identity
        case (6)
          mode = 'R'            ! read D,E,F do DC, and return
        case default
          mode = 'A'            ! mode 'A' + accuracy improvement
        end select
        check_accuracy = (merror == 1)

#if DEBUG
      if (i_inod == 0) then
        print*,"N := ",n
        print*,"NV:= ",nvec
        print*,"mf:= ",m
        print*,"mb:= ",mb
        print*,"nall:= ",nall
        print*,"mtyp:= ",mtype
        print*,"msol:= ",msolver
        print*,"merr:= ",merror
      end if
#endif
!-
!----------------------------------------------------------------------
!-
        if (new_comm /= MPI_COMM_NULL) then

          if (mtype < 0) then

            call mat_dim_get(mtype, n)
            if (n <= 0) then
              call flush(6)
              exit
            end if

          end if
! mode='O' is neccessary for performance
! and memory reduction at error check
          call eigen_get_matdims(n, nm, ny, mode='O')

          if (nm <= 0 .or. ny <= 0) then
            print*,"oversized problem", nm, ny
            call flush(6)
            exit
          end if

          imem = eigen_memory_internal(n, nm, nm, m, 128)
          if (imem <= 0) then
            print*,"oversized problem",imem
            call flush(6)
            exit
          end if

          allocate(
     &         a(nm, ny),
     &         z(nm, ny),
     &         w(n),
     &         stat=istat )
          if (istat /= 0) then
            print*,"Memory exhausted"
            call flush(6)
            exit
          end if

          call mat_set(n, a(1,1), nm, mtype)

        end if
!-
!----------------------------------------------------------------------
!-
        call MPI_Barrier(MPI_COMM_WORLD, ierror)
        d1 = MPI_Wtime()

        if (new_comm /= MPI_COMM_NULL) then
          if (msolver == 0) then
            call eigen_sx(n, nvec, a, nm, w, z, nm,
     &                m_forward=m, m_backward=mb, mode=mode)
          else
            call eigen_s (n, nvec, a, nm, w, z, nm,
!           call eigen_s0 (n, nvec, a, nm, w, z, nm,
     &                m_forward=m, m_backward=mb, mode=mode)
          end if

          flops = a(1, 1)
          commu = a(3, 1)

        end if

        call MPI_Barrier(MPI_COMM_WORLD, ierror)
        d2 = MPI_Wtime()
        time = d2 -d1
!-
!----------------------------------------------------------------------
!-
        call MPI_Barrier(MPI_COMM_WORLD, ierror)

        if ( i_inod == 0 ) then
          print*,"==================================================="
     &         //   "==="
          call eigen_show_version()

          if (msolver == 0) then
            print*,"Solver = eigen_sx / via penta-diagonal format"
          else
            print*,"Solver = eigen_s  / via tri-diagonal format"
          end if
          print*,"Block width = ", m, "/", mb
          print*,"NUM.OF.PROCESS=",nnod,"(",x_nnod,y_nnod,")"
!$OMP PARALLEL
!$OMP MASTER
!$          print*,"NUM.OF.THREADS=", omp_get_num_threads( )
!$OMP END MASTER
!$OMP END PARALLEL
          print*,"Matrix dimension = ",N
          call print_mat_name(mtype)
          print*,"Internally required memory = ",imem," [Byte]"
          print*,"The number of eigenvectors computed = ",nvec
          select case (nall)
          case (0)
            print*, "mode 'N' :: only eigenvalues, no eigenvector"
          case (1)
            print*, "mode 'A' :: all the eigenpairs"
          case (2)
            print*, "mode 'X' :: mode 'A' + accuracy improvement"
          case (3)
            print*, "mode 'S' :: skip DC but set Z as Identity"
          case (4)
            print*, "mode 'T' :: run DC but skip TRBAK"
          case (5)
            print*, "mode 'C' :: skip DC and TRBAK return X=identity"
          case (6)
            print*, "mode 'R' :: read D,E,F do DC, and return"
          end select
          print*,"Elapsed time = ",time," [sec]"
          print*,"FLOP         = ",abs(flops)
          print*,"Performance  = ",(abs(flops)/time)*1D-9," [GFLOPS]"
          if (flops <= 0) then
            print*,"* Since FLOPs on D&C could not be counted up"
     &           //" correctly, above performance"
            print*,"  is lower than the actual, which could be"
     &           //" 10-25 % higher :"
            print*," (", 
     &           (abs(flops)/time)*1D-9 * 1.1, "-",
     &           (abs(flops)/time)*1D-9 * 1.25,
     &           ")"
          end if
          if (commu >= 0) then
            print*,"Communication time = ",commu," [sec]"
            print*,"Ratio Communication/Elapsed time = ",commu/time
          end if

          if (check_accuracy) then
            call w_test( n, w, mtype)
          end if

        end if

        if (check_accuracy) then
          if (mode /= 'N' .and. nvec == N) then
            do i=0,n_comms-1

              call MPI_Barrier(MPI_COMM_WORLD, ierror)
              d1 = MPI_Wtime()

              if (color == i) then
                if (new_comm /= MPI_COMM_NULL) then
                  if (n_comms > 1) then
                    if (inod == 1) then
                      print*,"> Group ",i
                    end if
                  end if
                  if ( nvec > 0 ) then
                    call mat_set(n, a(1,1), nm, mtype)
                    call ev_test(n, nvec, 
     &                   a(1,1), nm, w(1), z(1,1), nm, mode)
                  end if
                end if
              end if

              call MPI_Barrier(MPI_COMM_WORLD, ierror)
              d2 = MPI_Wtime()
              if (inod == 1) then
                print*,"Time spent for error check : ",d2-d1," [sec]"
              end if

            end do
          end if
        end if

        if (i_inod == 0) then
          print*,"==================================================="
     &         //  "==="
          print*,""
        end if

        call MPI_Barrier(MPI_COMM_WORLD, ierror)
!-
!----------------------------------------------------------------------
!-
        if (new_comm /= MPI_COMM_NULL) then
          deallocate( a, z, w )
        end if

        call MPI_Barrier(MPI_COMM_WORLD, ierror)

      end do
!-
!----------------------------------------------------------------------
!-
      call MPI_Barrier(MPI_COMM_WORLD, ierror)
!-
!----------------------------------------------------------------------
!-
      if (i_inod == 0) then
        close(unit=10)
      end if
!-
!----------------------------------------------------------------------
!-
      call eigen_free()
!-
!----------------------------------------------------------------------
!-
      if (i_inod == 0) then
        print*,"Benchmark completed"
      end if
!-
!----------------------------------------------------------------------
!-
99999 continue
!-
!----------------------------------------------------------------------
!-
      call MPI_Finalize(ierror)

      contains
!-
!----------------------------------------------------------------------
!-
      subroutine print_header_message()

      print*,""
      print*," eigenexa_benchmark [options]"
      print*,""

      return

      end subroutine  print_header_message

      subroutine print_help_message()

      call print_header_message()

      print*,"options:"
      print*,""
      print*," -h              displays this help and exit"
      print*," -L              displays the list of test matrices"
      print*," -f input_file   uses input_file"
      print*,"                   default is ./IN"
      print*," -c or -n        check or non-check accuracy."
      print*,"                   default is check."
      print*," -g mode         sets the process grid as follows"
      print*,"    R, r           MPI_COMM_WORLD row-major mode"
      print*,"    C, c           MPI_COMM_WORLD column-major mode"
      print*,"    A, a           MPI_COMM_SELF (embarrasingly parallel)"
      print*,"    1, 2,... 9     splitted MPI_COMM_WORLD"
      print*,"                      with the color=MOD(RANK,{number})"
      print*," -x dimX dimY    sets the cartecian shape (dimX,dimY)"
      print*,"                      dimX <= dimY must be hold."
      print*,""

      return

      end subroutine print_help_message
!-
!----------------------------------------------------------------------
!-
      end program EigenExa_benchmark
