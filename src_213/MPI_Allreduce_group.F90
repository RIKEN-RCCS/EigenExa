module MPI_Group_property
  implicit none

  type MPI_Group_type
     integer,pointer :: sbuf_i(:)
     integer,pointer :: rbuf_i(:)
     real(kind(0.0d0)),pointer  :: sbuf_r8(:)
     real(kind(0.0d0)),pointer  :: rbuf_r8(:)

     integer icount
     integer idatatype
     integer iop
     integer icomm
     integer igroup
     integer igroup_rank
     integer igroup_size
     integer ierr

     ! -----

     integer icomm_group
     integer icomm_group_rank
     integer icomm_group_size
     ! icomm_group_ranklist --> groupランクとプロセスランクの変換を行う
     integer,allocatable :: icomm_group_ranklist(:)

  end type MPI_Group_type

end module MPI_Group_property

module Group_Allreduce_main
  use mpi
  use MPI_Group_property
  implicit none

  interface swap_ptr
     module procedure swap_ptr_integer
     module procedure swap_ptr_real8
  end interface swap_ptr
contains

  subroutine swap_ptr_real8(ptr1,ptr2)
    real(kind(0.0d0)),pointer :: tmp(:)
    real(kind(0.0d0)),pointer :: ptr1(:)
    real(kind(0.0d0)),pointer :: ptr2(:)
    tmp  => ptr1
    ptr1 => ptr2
    ptr2 => tmp
  end subroutine swap_ptr_real8

  subroutine swap_ptr_integer(ptr1,ptr2)
    integer,pointer :: tmp(:)
    integer,pointer :: ptr1(:)
    integer,pointer :: ptr2(:)
    tmp  => ptr1
    ptr1 => ptr2
    ptr2 => tmp
  end subroutine swap_ptr_integer

  subroutine comm_op_integer(&
       n,ihead, &
       sbuf,rbuf,&
       icomm,iop,&
       myrank,ipair)
    implicit none
    integer n
    integer ihead
    integer sbuf(n)
    integer rbuf(n)
    integer icomm
    integer iop
    integer myrank
    integer ipair
    ! ------
    integer ihead_s,ihead_r
    integer icount_s,icount_r
    integer ireq_s,ireq_r
    integer itag
    integer ierr
    integer i

    itag=1
    if(myrank<ipair)then
       ihead_r  = 1
       ihead_s  = n/2 + 1
       icount_r = n/2
       icount_s = n -  n/2 
    else
       ihead_r  = n/2 + 1
       ihead_s  = 1
       icount_r = n - n/2 
       icount_s = n/2
    endif
    call MPI_Isend(&
         sbuf(ihead_s),icount_s,MPI_INTEGER,ipair, &
         itag,icomm,ireq_s,ierr)

    call MPI_Irecv(&
         rbuf,icount_r,MPI_INTEGER,ipair,&
         itag,icomm,ireq_r,ierr)

    call MPI_Wait(ireq_s,MPI_STATUS_IGNORE,ierr)
    call MPI_Wait(ireq_r,MPI_STATUS_IGNORE,ierr)

    select case(iop)
    case(MPI_SUM)
       do i=1,icount_r
          sbuf(ihead_r + i -1) = sbuf(ihead_r + i - 1)+rbuf(i)
       enddo
    case(MPI_PROD)
       do i=1,icount_r
          sbuf(ihead_r + i -1) = sbuf(ihead_r + i - 1)*rbuf(i)
       enddo
    end select
    n = icount_r
    ihead = ihead + ihead_r -1
  end subroutine comm_op_integer

  subroutine comm_op_rev_integer(&
       n,ihead, &
       sbuf,rbuf,&
       icomm,iop,&
       myrank,ipair)
    implicit none
    integer n
    integer ihead
    integer sbuf(n)
    integer rbuf(n)
    integer icomm
    integer iop
    integer myrank
    integer ipair
    ! ------
    integer ihead_s,ihead_r
    integer icount_s,icount_r
    integer ireq_s,ireq_r
    integer itag
    integer ierr
    integer i

    itag=1
    if(myrank<ipair)then
       ihead_r  = n/2 + 1
       ihead_s  = 1
       icount_r = n - n/2 
       icount_s = n/2
    else
       ihead_r  = 1
       ihead_s  = n/2 + 1
       icount_r = n/2
       icount_s = n - n/2 
    endif
    call MPI_Isend(&
         sbuf(ihead_s),icount_s,MPI_INTEGER,ipair, &
         itag,icomm,ireq_s,ierr)

    call MPI_Irecv(&
         rbuf,icount_r,MPI_INTEGER,ipair,&
         itag,icomm,ireq_r,ierr)

    call MPI_Wait(ireq_s,MPI_STATUS_IGNORE,ierr)
    call MPI_Wait(ireq_r,MPI_STATUS_IGNORE,ierr)

    do i=1,icount_r
       sbuf(ihead_r + i -1) = rbuf(i)
    enddo

  end subroutine comm_op_rev_integer

  ! -----
  subroutine comm_op_real8(&
       n,ihead, &
       sbuf,rbuf,&
       icomm,iop,&
       myrank,ipair)
    implicit none
    integer n
    integer ihead
    real(kind(0.0d0)) sbuf(n)
    real(kind(0.0d0)) rbuf(n)
    integer icomm
    integer iop
    integer myrank
    integer ipair
    ! ------
    integer ihead_s,ihead_r
    integer icount_s,icount_r
    integer ireq_s,ireq_r
    integer itag
    integer ierr
    integer i

    itag=1
    if(myrank<ipair)then
       ihead_r  = 1
       ihead_s  = n/2 + 1
       icount_r = n/2
       icount_s = n -  n/2 
    else
       ihead_r  = n/2 + 1
       ihead_s  = 1
       icount_r = n - n/2 
       icount_s = n/2
    endif
    call MPI_Isend(&
         sbuf(ihead_s),icount_s,MPI_DOUBLE_PRECISION,ipair, &
         itag,icomm,ireq_s,ierr)

    call MPI_Irecv(&
         rbuf,icount_r,MPI_DOUBLE_PRECISION,ipair,&
         itag,icomm,ireq_r,ierr)

    call MPI_Wait(ireq_s,MPI_STATUS_IGNORE,ierr)
    call MPI_Wait(ireq_r,MPI_STATUS_IGNORE,ierr)

    select case(iop)
    case(MPI_SUM)
       do i=1,icount_r
          sbuf(ihead_r + i -1) = sbuf(ihead_r + i - 1)+rbuf(i)
       enddo
    case(MPI_PROD)
       do i=1,icount_r
          sbuf(ihead_r + i -1) = sbuf(ihead_r + i - 1)*rbuf(i)
       enddo
    end select
    n = icount_r
    ihead = ihead + ihead_r -1
  end subroutine comm_op_real8

  subroutine comm_op_rev_real8(&
       n,ihead, &
       sbuf,rbuf,&
       icomm,iop,&
       myrank,ipair)
    implicit none
    integer n
    integer ihead
    real(kind(0.0d0)) sbuf(n)
    real(kind(0.0d0)) rbuf(n)
    integer icomm
    integer iop
    integer myrank
    integer ipair
    ! ------
    integer ihead_s,ihead_r
    integer icount_s,icount_r
    integer ireq_s,ireq_r
    integer itag
    integer ierr
    integer i

    itag=1
    if(myrank<ipair)then
       ihead_r  = n/2 + 1
       ihead_s  = 1
       icount_r = n - n/2 
       icount_s = n/2
    else
       ihead_r  = 1
       ihead_s  = n/2 + 1
       icount_r = n/2
       icount_s = n - n/2 
    endif
    call MPI_Isend(&
         sbuf(ihead_s),icount_s,MPI_DOUBLE_PRECISION,ipair, &
         itag,icomm,ireq_s,ierr)

    call MPI_Irecv(&
         rbuf,icount_r,MPI_DOUBLE_PRECISION,ipair,&
         itag,icomm,ireq_r,ierr)

    call MPI_Wait(ireq_s,MPI_STATUS_IGNORE,ierr)
    call MPI_Wait(ireq_r,MPI_STATUS_IGNORE,ierr)

    do i=1,icount_r
       sbuf(ihead_r + i -1) = rbuf(i)
    enddo

  end subroutine comm_op_rev_real8

  subroutine Group_Allreduce(mygroup)
    implicit none
    type(MPI_Group_type) mygroup

    integer i

    integer myrank
    integer nprocess
    integer istep
    integer istep0 
    integer ipair
    integer icount
    integer ihead
    integer level
    integer,allocatable :: icount_level(:)
    integer,allocatable :: ihead_level(:)
    integer,allocatable :: ipair_level(:)

    myrank   = mygroup%igroup_rank
    nprocess = mygroup%igroup_size
    icount   = mygroup%icount

    ! ---
    level=1
    i = mygroup%igroup_size
    do while(i.ne.1)
       i=i/2
       level=level+1
       if(i/2.eq.0)then
          exit
       endif
    enddo
    allocate(icount_level(level))
    allocate(ihead_level(level))
    allocate(ipair_level(level))

    icount_level = -1
    ihead_level  = -1
    ipair_level  = -1
    ! ---

    ! 2べきのプロセスでのみ動作
    i = mygroup%igroup_size
    istep = mygroup%igroup_size/2
    if(istep <= myrank)then
       istep=-istep
    endif
    istep0 = 0
    level  = 0
    ihead  = 1
    do while(i.ne.1)
       i=i/2
       ipair = myrank + istep 
       level  = level +1
       icount_level(level) = icount
       ihead_level(level)  = ihead
       ipair_level(level)  = ipair
       select case(mygroup%idatatype)
       case(MPI_INTEGER)
          call comm_op_integer( &
               icount,ihead, &
               mygroup%sbuf_i(ihead:), &
               mygroup%rbuf_i, &
               mygroup%icomm,  &
               mygroup%iop,    &
               mygroup%icomm_group_ranklist(myrank+1), &
               mygroup%icomm_group_ranklist(ipair+1)   &
               )
       case(MPI_DOUBLE_PRECISION)
          call comm_op_real8( &
               icount,ihead, &
               mygroup%sbuf_r8(ihead:), &
               mygroup%rbuf_r8, &
               mygroup%icomm,  &
               mygroup%iop,    &
               mygroup%icomm_group_ranklist(myrank+1), &
               mygroup%icomm_group_ranklist(ipair+1)   &
               )
       end select
       if(i/2.eq.0)then
          exit
       endif
       if(istep<0)then
          istep0 = istep0 + i
       endif
       istep=abs(istep/2)
       if( (istep+istep0) <= myrank)then
          istep=-istep
       endif

    enddo
    ! ------
    ! 逆回転
    i=2
    do while(i.le.mygroup%igroup_size)
       icount = icount_level(level) 
       ihead  = ihead_level(level) 
       ipair  = ipair_level(level) 
       select case(mygroup%idatatype)
       case(MPI_INTEGER)
          call comm_op_rev_integer( &
               icount,ihead, &
               mygroup%sbuf_i(ihead:), &
               mygroup%rbuf_i, &
               mygroup%icomm,  &
               mygroup%iop,    &
               mygroup%icomm_group_ranklist(myrank+1), &
               mygroup%icomm_group_ranklist(ipair+1)   &
               )
       case(MPI_DOUBLE_PRECISION)
          call comm_op_rev_real8( &
               icount,ihead, &
               mygroup%sbuf_r8(ihead:), &
               mygroup%rbuf_r8, &
               mygroup%icomm,  &
               mygroup%iop,    &
               mygroup%icomm_group_ranklist(myrank+1), &
               mygroup%icomm_group_ranklist(ipair+1)   &
               )
       end select

       level  = level - 1
       i=i*2
    enddo
    ! ------

    select case(mygroup%idatatype)
    case(MPI_INTEGER)
       do i = 1,mygroup%icount
          mygroup%rbuf_i(i)=mygroup%sbuf_i(i)
       enddo
    case(MPI_DOUBLE_PRECISION)
       do i = 1,mygroup%icount
          mygroup%rbuf_r8(i)=mygroup%sbuf_r8(i)
       enddo
    end select

    deallocate(icount_level)
    deallocate(ihead_level)
    deallocate(ipair_level)

  end subroutine Group_Allreduce


end module Group_Allreduce_main

!module Group_Allreduce(mygroup)

module FS_MPI_Group
  use mpi
  use MPI_Group_property
  implicit none

  interface MPI_Group_Allreduce
     module procedure MPI_Group_Allreduce_integer
     module procedure MPI_Group_Allreduce_real8
  end interface MPI_Group_Allreduce

  interface set_group
     module procedure set_group_integer
     module procedure set_group_real8
  end interface set_group

  interface free_group
     module procedure free_group_integer
     module procedure free_group_real8
  end interface free_group

  public :: MPI_Group_Allreduce
  public :: set_group
  public :: free_group


  ! -------
  !   interface @
  !      module procedure @_integer
  !      module procedure @_real8
  !   end interface
  ! -------

contains
  subroutine set_group2comm_ranklist(mygroup)
    implicit none
    type(MPI_Group_type) mygroup
    integer,allocatable :: group_ranks(:)
    integer i

    allocate(group_ranks(mygroup%igroup_size))
    call MPI_Comm_group( &
         mygroup%icomm,mygroup%icomm_group,mygroup%ierr)
    call MPI_Group_size( &
         mygroup%icomm_group,mygroup%icomm_group_size,mygroup%ierr)
    call MPI_Group_rank( &
         mygroup%icomm_group,mygroup%icomm_group_rank,mygroup%ierr)
    allocate(mygroup%icomm_group_ranklist(mygroup%igroup_size))

    do i=1,mygroup%igroup_size
       group_ranks(i)=i-1
       mygroup%icomm_group_ranklist(i)=-1
    enddo
    call MPI_GROUP_TRANSLATE_RANKS(&
         mygroup%igroup, &
         mygroup%igroup_size, &
         group_ranks, &
         mygroup%icomm_group, &
         mygroup%icomm_group_ranklist, &
         mygroup%ierr)
    ! do i=1, mygroup%igroup_size
    !    write(*,*)i,group_ranks(i),mygroup%icomm_group_ranklist(i)
    ! enddo
    deallocate(group_ranks)

  end subroutine set_group2comm_ranklist

  subroutine free_group_integer(mygroup,rbuf,icount)
    implicit none
    type(MPI_Group_type) mygroup
    integer icount
    integer rbuf(icount)
    integer i

    deallocate(mygroup%icomm_group_ranklist)
    call MPI_Group_free(mygroup%icomm_group,mygroup%ierr)
    if( icount < mygroup%igroup_size )then
       do i=1,icount
          rbuf(i)=mygroup%rbuf_i(i)
       enddo
       deallocate(mygroup%sbuf_i)
       deallocate(mygroup%rbuf_i)
    endif

  end subroutine free_group_integer

  subroutine free_group_real8(mygroup,rbuf,icount)
    implicit none
    type(MPI_Group_type) mygroup
    integer icount
    real(kind(0.0d0)) rbuf(icount)
    integer i

    deallocate(mygroup%icomm_group_ranklist)
    call MPI_Group_free(mygroup%icomm_group,mygroup%ierr)
    if( icount < mygroup%igroup_size )then
       do i=1,icount
          rbuf(i)=mygroup%rbuf_r8(i)
       enddo
       deallocate(mygroup%sbuf_r8)
       deallocate(mygroup%rbuf_r8)
    endif

  end subroutine free_group_real8

  subroutine set_group_real8(&
       sbuf,rbuf,icount,&
       idatatype,iop,icomm,igroup,&
       ierr,&
       mygroup)
    implicit none
    integer,intent(IN) :: icount
    real(kind(0.0d0)),target     :: sbuf(icount)
    real(kind(0.0d0)),target     :: rbuf(icount)
    integer,intent(IN) :: idatatype
    integer,intent(IN) :: iop
    integer,intent(IN) :: icomm
    integer,intent(IN) :: igroup
    integer,intent(INOUT) :: ierr

    type(MPI_Group_type) mygroup

    integer i

    mygroup%icount = icount 
    mygroup%idatatype = idatatype
    mygroup%iop = iop
    mygroup%icomm = icomm
    mygroup%igroup = igroup
    mygroup%ierr = 0
    mygroup%ierr = ierr 

    call MPI_Group_size( &
         mygroup%igroup,mygroup%igroup_size,mygroup%ierr)
    call MPI_Group_rank( &
         mygroup%igroup,mygroup%igroup_rank,mygroup%ierr)
    call set_group2comm_ranklist(mygroup)

    if( mygroup%icount < mygroup%igroup_size )then
       allocate(mygroup%sbuf_r8(mygroup%igroup_size))
       allocate(mygroup%rbuf_r8(mygroup%igroup_size))

       do i=1,mygroup%igroup_size
          mygroup%rbuf_r8(i)=0
       enddo
       do i=1,mygroup%icount
          mygroup%sbuf_r8(i) = sbuf(i)
       enddo
       do i=mygroup%icount+1,mygroup%igroup_size
          mygroup%sbuf_r8(i) = 0
       enddo
       mygroup%icount = mygroup%igroup_size
    else
       mygroup%sbuf_r8 => sbuf
       mygroup%rbuf_r8 => rbuf
    endif
    mygroup%sbuf_i => null() 
    mygroup%rbuf_i => null()
  end subroutine set_group_real8

  subroutine set_group_integer(&
       sbuf,rbuf,icount,&
       idatatype,iop,icomm,igroup,&
       ierr,&
       mygroup)
    implicit none
    integer,intent(IN) :: icount
    integer,target     :: sbuf(icount)
    integer,target     :: rbuf(icount)
    integer,intent(IN) :: idatatype
    integer,intent(IN) :: iop
    integer,intent(IN) :: icomm
    integer,intent(IN) :: igroup
    integer,intent(INOUT) :: ierr

    type(MPI_Group_type) mygroup

    integer i

    mygroup%icount = icount 
    mygroup%idatatype = idatatype
    mygroup%iop = iop
    mygroup%icomm = icomm
    mygroup%igroup = igroup
    mygroup%ierr = 0
    mygroup%ierr = ierr 

    call MPI_Group_size( &
         mygroup%igroup,mygroup%igroup_size,mygroup%ierr)
    call MPI_Group_rank( &
         mygroup%igroup,mygroup%igroup_rank,mygroup%ierr)
    call set_group2comm_ranklist(mygroup)

    if( mygroup%icount < mygroup%igroup_size )then
       allocate(mygroup%sbuf_i(mygroup%igroup_size))
       allocate(mygroup%rbuf_i(mygroup%igroup_size))

       do i=1,mygroup%igroup_size
          mygroup%rbuf_i(i)=0
       enddo
       do i=1,mygroup%icount
          mygroup%sbuf_i(i) = sbuf(i)
       enddo
       do i=mygroup%icount+1,mygroup%igroup_size
          mygroup%sbuf_i(i) = 0
       enddo
       mygroup%icount = mygroup%igroup_size
    else
       mygroup%sbuf_i => sbuf
       mygroup%rbuf_i => rbuf
    endif
    mygroup%sbuf_r8 => null() 
    mygroup%rbuf_r8 => null() 

  end subroutine set_group_integer

  subroutine MPI_Group_Allreduce_integer(&
       sbuf,rbuf,icount,&
       idatatype,iop,icomm,igroup,&
       ierr)
    use Group_Allreduce_main
    implicit none
    integer,intent(IN) :: icount
    integer,target     :: sbuf(icount)
    integer,target     :: rbuf(icount)
    integer,intent(IN) :: idatatype
    integer,intent(IN) :: iop
    integer,intent(IN) :: icomm
    integer,intent(IN) :: igroup
    integer,intent(OUT) :: ierr

    type(MPI_Group_type) mygroup

    call set_group(&
         sbuf,rbuf,icount,&
         idatatype,iop,icomm,igroup,&
         ierr,&
         mygroup)

    call Group_Allreduce(mygroup)

    call free_group(mygroup,rbuf,icount)

  end subroutine MPI_Group_Allreduce_integer

  subroutine MPI_Group_Allreduce_real8( &
       sbuf,rbuf,icount,&
       idatatype,iop,icomm,igroup,&
       ierr)
    use Group_Allreduce_main
    implicit none
    integer,intent(IN) :: icount
    real(kind(0.0d0))      :: sbuf(icount)
    real(kind(0.0d0))      :: rbuf(icount)
    integer,intent(IN) :: idatatype
    integer,intent(IN) :: iop
    integer,intent(IN) :: icomm
    integer,intent(IN) :: igroup
    integer,intent(OUT) :: ierr

    type(MPI_Group_type) mygroup
    call set_group(&
         sbuf,rbuf,icount,&
         idatatype,iop,icomm,igroup,&
         ierr,&
         mygroup)
    call Group_Allreduce(mygroup)
    call free_group(mygroup,rbuf,icount)
  end subroutine MPI_Group_Allreduce_real8

end module FS_MPI_Group

! -------
! program main
!   use mpi
!   use MPI_Group
!   implicit none

!   integer npes
!   integer myrank
!   integer ierr

!   integer group_world
!   integer nrank
!   integer,allocatable :: ranklist(:)
!   integer new_group
!   integer new_comm

!   integer i
!   integer,parameter :: icount=5
!   integer sbuf_i(icount)
!   integer rbuf_i(icount)
!   real(kind(0.0d0)) sbuf_r(icount)
!   real(kind(0.0d0)) rbuf_r(icount)

!   call MPI_Init(ierr)
!   call MPI_Comm_size(MPI_COMM_WORLD,npes,ierr)
!   call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)

!   allocate(ranklist(npes))
!   ranklist=0

!   nrank=npes/2
!   if(mod(myrank,2)==0)then
!      do i=1,npes/2
!         ranklist(i) = (i-1)*2
!      enddo
!      if(mod(npes,2)==1)then
!         nrank=nrank+1
!         ranklist( npes/2+1 ) = ( npes/2+1 -1)*2
!      endif
!   else
!      do i=1,npes/2
!         ranklist(i) = (i-1)*2 + 1
!      enddo
!   endif

!   call MPI_Comm_group(MPI_COMM_WORLD, group_world, ierr)
!   call MPI_Group_incl(group_world, nrank, ranklist, new_group, ierr)

!   sbuf_i = 1 *myrank
!   rbuf_i = -1
!   sbuf_r = 1.0d0
!   rbuf_r = -1.0d0
!   call MPI_Group_Allreduce(&
!        sbuf_i,rbuf_i,icount,&
!        MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,new_group,&
!        ierr)
!   call MPI_Group_Allreduce(&
!        sbuf_r,rbuf_r,icount,&
!        MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,new_group,&
!        ierr)

!   if(myrank.eq.0)then
!      do i=1,icount
!         write(*,*)i,sbuf_i(i),rbuf_i(i)   
!      enddo
!      do i=1,icount
!         write(*,*)i,sbuf_r(i),rbuf_r(i)   
!      enddo
!   endif

!   call MPI_Group_free(new_group,ierr)
!   call MPI_Finalize(ierr)
!   deallocate(ranklist)

! end program main
