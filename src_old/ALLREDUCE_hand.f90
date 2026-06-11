! Floor of binary logarithm
! If argument n is less than 1 then assign dummy(-9999) to n.
! If n is 2**i and i is integer(kind=intType) then assign exactly i to n.
! Otherwise, assign floor(log_2 (n)) to n.
subroutine floor_log2(n)
  implicit none
  integer,intent(inout) :: n
  integer :: power_of_2, power

  if(n.le.0)then
     ! If n<=0, no log2(n) exists. So, assign n to dummy(-9999).
     n=-9999
     return
  endif

  power_of_2=2
  do power=0,29
     if(n.lt.power_of_2)then
        n=power
        return
     endif
     power_of_2=power_of_2*2
  enddo
  ! 4 Byte integer type variable is definitely less than 2**31
  n=30

end subroutine floor_log2

subroutine ALLREDUCE_bin(comm1,s,R)
  use mpi
  implicit none
  integer s
  integer comm1
  integer i,j,k
  integer i_,j_

  integer s2
  integer ierr
  integer irank,nprocess
  integer stage
  integer stat_recv(mpi_status_size),stat_send(mpi_status_size)

  integer req_send,req_recv
  integer tag,nrk
  real(8) rbuff(s),sbuff(s)
  !    real(8) A(2*s,s)

  real(8) tmp
  integer itmp
  integer istage,irank_stage,irank0,irank0_master
  real(8) R(s)

  !     real(8) tau(2*s*s)
  !     real(8) work(2*s*s)
  !     integer lwork
  integer ss_int4
  !    lwork=2*s*s

  call MPI_Comm_size(comm1,nprocess, ierr)
  call MPI_Comm_rank(comm1,irank, ierr)

  ! rbuff=-99.0d0
  ! sbuff=irank
  !     A=0.0d0
  do k=1,s
     sbuff(k)=R(k)
  enddo
  !     do i=1,s
  !        do j=1,s
  !           A(j,i)=R(j+s*(i-1))
  !        enddo
  !     enddo

  stage=0
  irank0=0
  irank0_master=0

  istage=nprocess
  call floor_log2(istage)

  irank_stage=irank
  do i=istage,0,-1

     if(irank_stage<2**i)then
        ! ----
        ! irank0_master --> 集約された三角行列を持つランク
        ! irank0_master - irank0 --> 集約した三角行列の送信先ランク
        ! -----
        j=0
        stage=1
        do while(stage<2**i)
           j=j+1
           stage=2**j
           !           if(mod(irank,stage).eq.0)then
           if(mod(irank,stage).eq.0)then
              if(irank.eq.(nprocess-1))then
                 nrk=irank0_master
                 stage=2**i
                 !! 送信先が自プロセスでない場合
                 if(nrk.ne.irank)then
                    !write(*,*)"Send stage = ",stage,irank,nrk
                    ss_int4=s
                    call MPI_Isend(sbuff,ss_int4, MPI_DOUBLE, nrk, nrk,  comm1, req_send,ierr)
                    call MPI_Wait(req_send, stat_send,ierr)
                 endif

                 exit
              else
                 nrk=irank+stage/2
                 !! 受信先がある場合
                 if(nrk<(nprocess))then
                    !!write(*,*)"Recv stage = ",stage,irank,nrk
                    ss_int4=s
                    call MPI_Irecv(rbuff,ss_int4,MPI_DOUBLE,nrk,irank,comm1,req_recv,ierr)
                    call MPI_Wait(req_recv, stat_recv,ierr)
                    !write(*,*)"Recv stage = ",stage,irank,nrk,rbuff

                    do k=1,s
                       R(k) = R(k) + rbuff(k)
                       sbuff(k) = R(k)
                    enddo

                 endif
              endif
           else
              nrk=irank-stage/2
              !write(*,*)"Send stage = ",stage,irank,nrk   
              ss_int4=s
              call MPI_Isend(sbuff, ss_int4, MPI_DOUBLE, nrk, nrk,  comm1, req_send,ierr)
              call MPI_Wait(req_send, stat_send,ierr)
              exit
           endif

           !              write(*,*)stage,nrk,irank
        end do
        ! -----
        ! ここで1ステージ上のツリーに集約を行う

        if((irank0_master.eq.irank).and.(irank.ne.0))then
           if((irank+2**i)<nprocess)then
              nrk=irank+2**i
              ss_int4=s
              call MPI_Irecv(rbuff,ss_int4, MPI_DOUBLE, nrk, irank,   comm1, req_recv,ierr)
              call MPI_Wait(req_recv, stat_recv,ierr)             
              !write(*,*)"Recv stage = ",stage,irank,nrk,rbuff

              do k=1,s
                 R(k) = R(k) + rbuff(k)
                 sbuff(k) = R(k)
              enddo

           endif

           nrk=irank0_master - irank0
           !write(*,*)"Send stage = ",stage,irank,nrk
           ss_int4=s
           call MPI_Isend(sbuff, ss_int4, MPI_DOUBLE, nrk, nrk,  comm1, req_send,ierr)
           call MPI_Wait(req_send, stat_send,ierr)
        endif
        if((irank.eq.0).and.(2**(i)<nprocess))then
           nrk=irank+2**(i)
           ss_int4=s
           !!write(*,*)"Recv stage = ",stage,irank,nrk
           call MPI_Irecv(rbuff, ss_int4, MPI_DOUBLE, nrk, irank,   comm1, req_recv,ierr)
           call MPI_Wait(req_recv, stat_recv,ierr)             

           do k=1,s
              R(k) = R(k) + rbuff(k)
              sbuff(k) = R(k)
           enddo

           !write(*,*)"Recv stage = ",stage,irank,nrk,rbuff
        endif
        exit
     else
        irank0=2**(i)
        irank0_master=irank0_master+irank0
        irank_stage=mod(irank_stage,2**i)
     endif
  enddo
  ss_int4=s
  call MPI_Bcast(R,ss_int4,MPI_DOUBLE,0,comm1, ierr)

end subroutine ALLREDUCE_BIN

subroutine ALLREDUCE_prod(comm,n,buff,buff0)
  use mpi
  implicit none
  integer comm
  integer n
  real(8) buff(n)
  real(8) buff0(n)

  real(8) buff_c(n)
  real(8) yy,tt
  integer ierr 
  integer nprocess,myrank
  integer irank
  integer i,j
  integer req_send

  real(8),allocatable :: rbuff(:)
  integer,allocatable :: req_recv(:)
  integer stat_recv(mpi_status_size),stat_send(mpi_status_size)

  call MPI_Comm_size(comm,nprocess, ierr)
  call MPI_Comm_rank(comm,myrank, ierr)

  if(myrank.eq.0)then
     allocate(rbuff(n*nprocess))
     allocate(req_recv(nprocess))
     do i=1,n
        rbuff(i) = buff(i)
     enddo
     do i=1+1,nprocess
        irank = i-1
        call MPI_Irecv(rbuff(1+irank*n),n,MPI_DOUBLE,irank,irank,comm,req_recv(irank),ierr)
     enddo

     do i=1+1,nprocess
        irank = i-1
        call MPI_Wait(req_recv(irank), stat_recv,ierr)
     enddo
     do i=1,n
        buff_c(i)=0.0d0
     enddo
     do j=1+1,nprocess
        irank = j-1
        do i=1,n
           buff(i)=buff(i)*rbuff(i+irank*n)
        enddo
     enddo
  else
     irank = 0
     call MPI_Isend(buff,n, MPI_DOUBLE, 0, myrank,  comm, req_send,ierr)
     call MPI_Wait(req_send, stat_send,ierr)
  endif
  call MPI_Bcast(buff,n,MPI_DOUBLE,0,comm, ierr)

  if(myrank.eq.0)then
     deallocate(rbuff)
     deallocate(req_recv)
  endif
  buff0=buff
end subroutine ALLREDUCE_prod

subroutine ALLREDUCE_sum(comm,n,buff,buff0)
  use mpi
  implicit none
  integer comm
  integer n
  real(8) buff(n)
  real(8) buff0(n)

  real(8) buff_c(n)
  real(8) yy,tt
  integer ierr 
  integer nprocess,myrank
  integer irank
  integer i,j
  integer req_send

  real(8),allocatable :: rbuff(:)
  integer,allocatable :: req_recv(:)
  integer stat_recv(mpi_status_size),stat_send(mpi_status_size)

  call MPI_Comm_size(comm,nprocess, ierr)
  call MPI_Comm_rank(comm,myrank, ierr)

  if(myrank.eq.0)then
     allocate(rbuff(n*nprocess))
     allocate(req_recv(nprocess))
     do i=1,n
        rbuff(i) = buff(i)
     enddo
     do i=1+1,nprocess
        irank = i-1
        call MPI_Irecv(rbuff(1+irank*n),n,MPI_DOUBLE,irank,irank,comm,req_recv(irank),ierr)
     enddo

     do i=1+1,nprocess
        irank = i-1
        call MPI_Wait(req_recv(irank), stat_recv,ierr)
     enddo
     do i=1,n
        buff_c(i)=0.0d0
     enddo
     do j=1+1,nprocess
        irank = j-1
#if 0
        do i=1,n
           buff(i)=buff(i)+rbuff(i+irank*n)
        enddo
#else
        do i=1,n
           yy = rbuff(i+irank*n) -  buff_c(i)
           tt = buff(i) + yy
           buff_c(i) = (tt - buff(i)) - yy
           buff(i) = tt
        enddo
#endif
     enddo
  else
     irank = 0
     call MPI_Isend(buff,n, MPI_DOUBLE, 0, myrank,  comm, req_send,ierr)
     call MPI_Wait(req_send, stat_send,ierr)
  endif
  call MPI_Bcast(buff,n,MPI_DOUBLE,0,comm, ierr)

  if(myrank.eq.0)then
     deallocate(rbuff)
     deallocate(req_recv)
  endif

  buff0=buff

end subroutine ALLREDUCE_sum

program main
  use mpi
  implicit none
  integer ierr
  real(8) buff0(100),buff1(100),buff2(100)
  integer n 
  integer myrank
  integer i

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,myrank, ierr)

  n =100
  buff0=0.0d0
  buff1=0.0d0
  buff2=0.0d0
  call random_number(buff1)
  buff2=buff1
  buff0=0.0d0

  call MPI_Allreduce(buff1,buff0,n,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
  !  call ALLREDUCE_sum(MPI_COMM_WORLD,n,buff2,buff1)
  !  buff2=buff1
  call ALLREDUCE_bin(MPI_COMM_WORLD,n,buff2)
  if(myrank.eq.0)then
     do i=1,n
        write(*,*)i,buff0(i),buff2(i),( buff0(i) - buff2(i) )/buff0(i)
     enddo
  endif
  call MPI_Finalize(ierr)

end program main

