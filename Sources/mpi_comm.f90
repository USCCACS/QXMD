



module parallel
!-----------------------------------------------------------------------
! type declaration of variables in mpi_comm.f
!-----------------------------------------------------------------------
implicit none

integer :: nodes, iam
integer :: world_
save

end module




!-----------------------------------------------------------------------
subroutine set_world( world )
!-----------------------------------------------------------------------
use parallel
implicit none
integer :: world

world_ = world

return
end




!-----------------------------------------------------------------------
subroutine get_world( world )
!-----------------------------------------------------------------------
use parallel
implicit none
integer :: world

world = world_

return
end




!-----------------------------------------------------------------------
subroutine start_parallel( nprocs, myid, world )
!-----------------------------------------------------------------------
use parallel
implicit none
include 'mpif.h'
integer :: nprocs, myid, world
integer :: ierr, nbuf, ierror
character(72)  :: exefile
character(256) :: buf

call MPI_INIT(ierr)
if(ierr /= MPI_SUCCESS) then
  write(*,*) 'MPI_INIT: ierr =',ierr
  call MPI_ERROR_STRING(ierr,buf,nbuf,ierror)
  write(*,*) buf
endif
world_ = MPI_COMM_WORLD

call comm_rank( iam, MPI_COMM_WORLD )
!      call MPI_COMM_RANK(MPI_COMM_WORLD,iam,ierr)
!      if(ierr.ne.MPI_SUCCESS) then
!        write(*,*) 'MPI_COMM_RANK: ierr =',ierr
!        call MPI_ERROR_STRING(ierr,buf,nbuf,ierror)
!        write(*,*) buf
!      endif

call comm_size( nodes, MPI_COMM_WORLD )
!      call MPI_COMM_SIZE(MPI_COMM_WORLD,nodes,ierr)
!      if(ierr.ne.MPI_SUCCESS) then
!        write(*,*) 'MPI_COMM_SIZE: ierr =',ierr
!        call MPI_ERROR_STRING(ierr,buf,nbuf,ierror)
!        write(*,*) buf
!      endif

!      write(*,*) ' iam =',iam,'  of  ',nodes
nprocs= nodes
myid  = iam
world = MPI_COMM_WORLD


return
end


!-----------------------------------------------------------------------
subroutine end_parallel(info)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: info, ierr

!-----Finalize the MPI environment
call MPI_FINALIZE(ierr)
info = ierr

return
end


!-----------------------------------------------------------------------
subroutine abort_parallel(info)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: info, ierr, world

!-----Finalize the MPI environment
call get_world( world )

call MPI_Abort(world, ierr)
info = ierr

return
end


!-----------------------------------------------------------------------
subroutine comm_rank( myid, world )
!-----------------------------------------------------------------------
use parallel
implicit none
include 'mpif.h'
integer :: myid, world
integer :: ierr, nbuf, ierror
character(256) :: buf

if( world /= MPI_COMM_NULL ) then
   call MPI_COMM_RANK( world, myid, ierr )
   if( ierr /= MPI_SUCCESS ) then
      write(*,*) 'MPI_COMM_RANK: ierr =',ierr
      call MPI_ERROR_STRING(ierr,buf,nbuf,ierror)
      write(*,*) buf
   end if
  else
   myid = MPI_UNDEFINED
end if

return
end


!-----------------------------------------------------------------------
subroutine comm_size( nprocs, world )
!-----------------------------------------------------------------------
use parallel
implicit none
include 'mpif.h'
integer :: nprocs, world
integer :: ierr, nbuf, ierror
character(256) :: buf

if( world /= MPI_COMM_NULL ) then
   call MPI_COMM_SIZE( world, nprocs, ierr )
   if( ierr /= MPI_SUCCESS ) then
      write(*,*) 'MPI_COMM_SIZE: ierr =',ierr
      call MPI_ERROR_STRING(ierr,buf,nbuf,ierror)
      write(*,*) buf
   end if
  else
   nprocs = MPI_UNDEFINED
end if

return
end


!-----------------------------------------------------------------------
subroutine comm_group( world, group )
!-----------------------------------------------------------------------
use parallel
implicit none
include 'mpif.h'
integer :: world, group
integer :: ierr

call MPI_COMM_GROUP( world, group, ierr )

return
end


!-----------------------------------------------------------------------
subroutine comm_split( oldcomm, color, key, newcomm )
!-----------------------------------------------------------------------
use parallel
implicit none
include 'mpif.h'
integer :: oldcomm, color, key, newcomm
integer :: ierr, nbuf, ierror
character(256) :: buf

if( color < 0 ) color = MPI_UNDEFINED
call MPI_COMM_SPLIT( oldcomm, color, key, newcomm, ierr )
if( ierr /= MPI_SUCCESS ) then
    write(*,*) 'comm_split: ierr =',ierr
    call MPI_ERROR_STRING(ierr,buf,nbuf,ierror)
    write(*,*) buf
    stop
end if

return
end


!-----------------------------------------------------------------------
subroutine group_excl( group, extract, extract_ids, subgroup )
!-----------------------------------------------------------------------
use parallel
implicit none
include 'mpif.h'
integer :: group, extract, subgroup
integer, dimension(extract) :: extract_ids
integer :: ierr, nbuf, ierror
character(256) :: buf

call MPI_GROUP_EXCL( group, extract, extract_ids, subgroup, ierr )
if( ierr /= MPI_SUCCESS ) then
    write(*,*) 'group_excl: ierr =',ierr
    call MPI_ERROR_STRING(ierr,buf,nbuf,ierror)
    write(*,*) buf
    stop
end if

return
end


!-----------------------------------------------------------------------
subroutine comm_create( world, group_sub, world_sub )
!-----------------------------------------------------------------------
use parallel
implicit none
include 'mpif.h'
integer :: world, group_sub, world_sub
integer :: ierr, nbuf, ierror
character(256) :: buf

call MPI_COMM_CREATE( world, group_sub, world_sub, ierr )
if( ierr /= MPI_SUCCESS ) then
    write(*,*) 'comm_create: ierr =',ierr
    call MPI_ERROR_STRING(ierr,buf,nbuf,ierror)
    write(*,*) buf
    stop
end if

return
end


!-----------------------------------------------------------------------
subroutine group_free( group )
!-----------------------------------------------------------------------
use parallel
implicit none
include 'mpif.h'
integer :: group
integer :: ierr

call MPI_GROUP_FREE( group, ierr )

return
end


!-----------------------------------------------------------------------
subroutine comm_free( world )
!-----------------------------------------------------------------------
use parallel
implicit none
include 'mpif.h'
integer :: world
integer :: ierr

if( world /= MPI_COMM_NULL ) then
    call MPI_COMM_FREE( world, ierr )
end if

return
end


!-----------------------------------------------------------------------
subroutine comm_dup( world, world_sub )
!-----------------------------------------------------------------------
use parallel
implicit none
include 'mpif.h'
integer :: world, world_sub
integer :: ierr, nbuf, ierror
character(256) :: buf

call MPI_COMM_DUP( world, world_sub, ierr )
if( ierr /= MPI_SUCCESS ) then
    write(*,*) 'comm_dup: ierr =',ierr
    call MPI_ERROR_STRING(ierr,buf,nbuf,ierror)
    write(*,*) buf
    stop
end if

return
end




!-----------------------------------------------------------------------
subroutine csend(mesg_tag,mesg_buff,size,inode,pid)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: mesg_tag,size,pid,inode,tid,ierr
character :: mesg_buff(*)

call MPI_SEND(mesg_buff,size,MPI_BYTE,inode,mesg_tag,  &
&             world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine cdsend(mesg_tag,mesg_buff,size,inode,pid)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: mesg_tag,size,pid,inode,tid,ierr
real*8  :: mesg_buff(*)

call MPI_SEND(mesg_buff,size,MPI_DOUBLE_PRECISION,inode,mesg_tag,  &
&             world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine crsend(mesg_tag,mesg_buff,size,inode,pid)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: mesg_tag,size,pid,inode,tid,ierr
real    :: mesg_buff(*)

call MPI_SEND(mesg_buff,size,MPI_REAL,inode,mesg_tag,  &
&             world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine cisend(mesg_tag,mesg_buff,size,inode,pid)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: mesg_tag,size,pid,inode,tid,ierr
integer :: mesg_buff(*)

call MPI_SEND(mesg_buff,size,MPI_INTEGER,inode,mesg_tag,  &
&             world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine clsend(mesg_tag,mesg_buff,size,inode,pid)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: mesg_tag,size,pid,inode,tid,ierr
logical :: mesg_buff(*)

call MPI_SEND(mesg_buff,size,MPI_LOGICAL,inode,mesg_tag,  &
&             world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine crecv(mesg_tag,mesg_buff,size,pid)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: mesg_tag,size,bufid,pid
character :: mesg_buff(*)
integer :: stat(MPI_STATUS_SIZE) 

call MPI_RECV(mesg_buff,size,MPI_BYTE,MPI_ANY_SOURCE,  &
&             mesg_tag,world_,stat,ierr)

return
end


!-----------------------------------------------------------------------
subroutine cdrecv(mesg_tag,mesg_buff,size,pid)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: mesg_tag,size,bufid,pid
real*8  :: mesg_buff(*)
integer :: stat(MPI_STATUS_SIZE) 

call MPI_RECV(mesg_buff,size,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,  &
&             mesg_tag,world_,stat,ierr)

return
end


!-----------------------------------------------------------------------
subroutine crrecv(mesg_tag,mesg_buff,size,pid)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: mesg_tag,size,bufid,pid
real    :: mesg_buff(*)
integer :: stat(MPI_STATUS_SIZE) 

call MPI_RECV(mesg_buff,size,MPI_REAL,MPI_ANY_SOURCE,  &
&             mesg_tag,world_,stat,ierr)

return
end


!-----------------------------------------------------------------------
subroutine cirecv(mesg_tag,mesg_buff,size,pid)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: mesg_tag,size,bufid,pid
integer :: mesg_buff(*)
integer :: stat(MPI_STATUS_SIZE) 

call MPI_RECV(mesg_buff,size,MPI_INTEGER,MPI_ANY_SOURCE,  &
&             mesg_tag,world_,stat,ierr)

return
end


!-----------------------------------------------------------------------
subroutine clrecv(mesg_tag,mesg_buff,size,pid)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: mesg_tag,size,bufid,pid
logical :: mesg_buff(*)
integer :: stat(MPI_STATUS_SIZE) 

call MPI_RECV(mesg_buff,size,MPI_LOGICAL,MPI_ANY_SOURCE,  &
&             mesg_tag,world_,stat,ierr)

return
end


!-----------------------------------------------------------------------
subroutine cdrecvs(mesg_tag,mesg_buff,size,source,pid)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: mesg_tag,size,source,pid
real*8  :: mesg_buff(*)
integer :: stat(MPI_STATUS_SIZE) 

call MPI_RECV(mesg_buff,size,MPI_DOUBLE_PRECISION,source,  &
&             mesg_tag,world_,stat,ierr)

return
end


!-----------------------------------------------------------------------
subroutine cirecvs(mesg_tag,mesg_buff,size,source,pid)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: mesg_tag,size,source,pid
integer :: mesg_buff(*)
integer :: stat(MPI_STATUS_SIZE) 

call MPI_RECV(mesg_buff,size,MPI_INTEGER,source,  &
&             mesg_tag,world_,stat,ierr)

return
end


!-----------------------------------------------------------------------
subroutine clrecvs(mesg_tag,mesg_buff,size,source,pid)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: mesg_tag,size,source,pid
logical :: mesg_buff(*)
integer :: stat(MPI_STATUS_SIZE) 

call MPI_RECV(mesg_buff,size,MPI_LOGICAL,source,  &
&             mesg_tag,world_,stat,ierr)

return
end


!-----------------------------------------------------------------------
subroutine gdsum(mesg_buff,n,temp)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
real*8  :: mesg_buff(*),temp(*)

if(nodes.eq.1) return
do i=1,n
   temp(i)=mesg_buff(i)
end do
call MPI_ALLREDUCE(temp,mesg_buff,n,MPI_DOUBLE_PRECISION,  &
&    MPI_SUM,world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine dsum(mesg_buff,n,temp,root)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
real*8  :: mesg_buff(*),temp(*)
integer :: root

if(nodes.eq.1) return
do i=1,n
   temp(i)=mesg_buff(i)
end do
call MPI_REDUCE(temp,mesg_buff,n,MPI_DOUBLE_PRECISION,  &
&    MPI_SUM,root,world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine gisum(mesg_buff,n,temp)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: mesg_buff(*),temp(*)

if(nodes.eq.1) return
do i=1,n
   temp(i)=mesg_buff(i)
end do
call MPI_ALLREDUCE(temp,mesg_buff,n,MPI_INTEGER,  &
&    MPI_SUM,world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine gdmax(mesg_buff)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
real*8  :: mesg_buff,temp

if(nodes.eq.1) return
temp = mesg_buff
call MPI_ALLREDUCE(temp,mesg_buff,1,MPI_DOUBLE_PRECISION,MPI_MAX,  &
&     world_,ierr)

!	if(iam.eq.0) then
!	  do imesg=1,nodes-1
!	    call crecv(1000,temp,8,0)
!	    if(mesg_buff.lt.temp) mesg_buff=temp
!	  end do
!	  do imesg=1,nodes-1
!	     call csend(2000,mesg_buff,8,imesg,0)
!	  end do
!	else
!	  call csend(1000,mesg_buff,8,0,0)
!	  call crecv(2000,mesg_buff,8,0)
!	end if

return
end


!-----------------------------------------------------------------------
subroutine gdmaxpl(mesg_buff,n,temp)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: n
real*8, dimension(n) :: mesg_buff,temp

if(nodes.eq.1) return
do i=1,n
   temp(i)=mesg_buff(i)
end do
call MPI_ALLREDUCE(temp,mesg_buff,n,MPI_DOUBLE_PRECISION,MPI_MAX,  &
&     world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine gimax(mesg_buff)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: mesg_buff,temp

if(nodes.eq.1) return
temp = mesg_buff
call MPI_ALLREDUCE(temp,mesg_buff,1,MPI_INTEGER,MPI_MAX,  &
&    world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine gimaxpl(mesg_buff,n,temp)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: n
integer, dimension(n) ::  mesg_buff,temp

if(nodes.eq.1) return
do i=1,n
   temp(i)=mesg_buff(i)
end do
call MPI_ALLREDUCE(temp,mesg_buff,n,MPI_INTEGER,MPI_MAX,  &
&    world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine gdmin(mesg_buff)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
real*8  :: mesg_buff,temp

if(nodes.eq.1) return
temp = mesg_buff
call MPI_ALLREDUCE(temp,mesg_buff,1,MPI_DOUBLE_PRECISION,MPI_MIN,  &
&    world_,ierr)
!      if(iam.eq.0) then
!         do imesg=1,nodes-1
!            call crecv(1000,temp,8,0)
!            if(mesg_buff.gt.temp) mesg_buff=temp
!          end do
!          do imesg=1,nodes-1
!             call csend(2000,mesg_buff,8,imesg,0)
!          end do
!        else
!          call csend(1000,mesg_buff,8,0,0)
!          call crecv(2000,mesg_buff,8,0)
!        end if

return
end


!-----------------------------------------------------------------------
subroutine gdminpl(mesg_buff,n,temp)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: n
real*8, dimension(n) :: mesg_buff,temp

if(nodes.eq.1) return
temp(:) = mesg_buff(:)
call MPI_ALLREDUCE(temp,mesg_buff,n,MPI_DOUBLE_PRECISION,MPI_MIN,  &
&    world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine gimin(mesg_buff)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: mesg_buff,temp

if(nodes.eq.1) return
temp = mesg_buff
call MPI_ALLREDUCE(temp,mesg_buff,1,MPI_INTEGER,MPI_MIN,  &
&    world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine giminpl(mesg_buff,n,temp)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: n
real*8, dimension(n) :: mesg_buff,temp

if(nodes.eq.1) return
temp(:) = mesg_buff(:)
call MPI_ALLREDUCE(temp,mesg_buff,n,MPI_INTEGER,MPI_MIN,  &
&    world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine dbcast(mesg_buff,n,root)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
real*8  :: mesg_buff(*)
integer :: root

if(nodes.eq.1) return
call MPI_BCAST(mesg_buff,n,MPI_DOUBLE_PRECISION,  &
&              root,world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine rbcast(mesg_buff,n,root)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
real    :: mesg_buff(*)
integer :: root

if(nodes.eq.1) return
call MPI_BCAST(mesg_buff,n,MPI_REAL,  &
&              root,world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine ibcast(mesg_buff,n,root)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: mesg_buff(*)
integer :: root

if(nodes.eq.1) return
call MPI_BCAST(mesg_buff,n,MPI_INTEGER,  &
&              root,world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine lbcast(mesg_buff,n,root)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
logical :: mesg_buff(*)
integer :: root

if(nodes.eq.1) return
call MPI_BCAST(mesg_buff,n,MPI_LOGICAL,  &
&              root,world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine mespasd(mesg_tag,buffs,buffr,nsd,nrc,inode,pid)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: mesg_tag,rreq,sreq,pid,inode,tid,ierr
integer :: rstat(MPI_STATUS_SIZE),sstat(MPI_STATUS_SIZE)
real*8  ::  buffs(*), buffr(*)

if( nrc.gt.0 )  &
& call MPI_IRECV(buffr,nrc,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,  &
&             mesg_tag,world_,rreq,ierr)
if( nsd.gt.0 )  &
& call MPI_ISEND(buffs,nsd,MPI_DOUBLE_PRECISION,inode,mesg_tag,  &
&             world_,sreq,ierr)
!-------Wait for recv & send to complete
if( nrc.gt.0 )  &
&   call MPI_WAIT(rreq,rstat,ierr)
if( nsd.gt.0 )  &
&   call MPI_WAIT(sreq,sstat,ierr)

return
end


!-----------------------------------------------------------------------
subroutine mespasr(mesg_tag,buffs,buffr,nsd,nrc,inode,pid)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: mesg_tag,rreq,sreq,pid,inode,tid,ierr
integer :: rstat(MPI_STATUS_SIZE),sstat(MPI_STATUS_SIZE)
real    :: buffs(*), buffr(*)

if( nrc.gt.0 )  &
& call MPI_IRECV(buffr,nrc,MPI_REAL,MPI_ANY_SOURCE,  &
&             mesg_tag,world_,rreq,ierr)
if( nsd.gt.0 )  &
& call MPI_ISEND(buffs,nsd,MPI_REAL,inode,mesg_tag,  &
&             world_,sreq,ierr)
!-------Wait for recv & send to complete
if( nrc.gt.0 )  &
&   call MPI_WAIT(rreq,rstat,ierr)
if( nsd.gt.0 )  &
&   call MPI_WAIT(sreq,sstat,ierr)

return
end


!-----------------------------------------------------------------------
subroutine mespasi(mesg_tag,buffs,buffr,nsd,nrc,inode,pid)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: mesg_tag,rreq,sreq,pid,inode,tid,ierr
integer :: rstat(MPI_STATUS_SIZE),sstat(MPI_STATUS_SIZE)
integer :: buffs(*), buffr(*)

if( nrc.gt.0 )  &
& call MPI_IRECV(buffr,nrc,MPI_INTEGER,MPI_ANY_SOURCE,  &
&             mesg_tag,world_,rreq,ierr)
if( nsd.gt.0 )  &
& call MPI_ISEND(buffs,nsd,MPI_INTEGER,inode,mesg_tag,  &
&             world_,sreq,ierr)
!-------Wait for recv & send to complete
if( nrc.gt.0 )  &
&   call MPI_WAIT(rreq,rstat,ierr)
if( nsd.gt.0 )  &
&   call MPI_WAIT(sreq,sstat,ierr)

return
end


!-----------------------------------------------------------------------
subroutine mespasl(mesg_tag,buffs,buffr,nsd,nrc,inode,pid)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: mesg_tag,rreq,sreq,pid,inode,tid,ierr
integer :: rstat(MPI_STATUS_SIZE),sstat(MPI_STATUS_SIZE)
logical :: buffs(*), buffr(*)

if( nrc.gt.0 )  &
& call MPI_IRECV(buffr,nrc,MPI_LOGICAL,MPI_ANY_SOURCE,  &
&             mesg_tag,world_,rreq,ierr)
if( nsd.gt.0 )  &
& call MPI_ISEND(buffs,nsd,MPI_LOGICAL,inode,mesg_tag,  &
&             world_,sreq,ierr)
!-------Wait for recv & send to complete
if( nrc.gt.0 )  &
&   call MPI_WAIT(rreq,rstat,ierr)
if( nsd.gt.0 )  &
&   call MPI_WAIT(sreq,sstat,ierr)

return
end


!-----------------------------------------------------------------------
subroutine gcol(x,msize,x_tot)
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: msize, ierr
character :: x(*), x_tot(*)

call MPI_ALLGATHER(x,msize,MPI_BYTE,  &
&     x_tot,msize,MPI_BYTE,world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine alldgatherv( sendbuf,sendcount,recvbuf,recvcounts,      &
&                       rdispls )
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
real*8  :: sendbuf(*), recvbuf(*)
integer :: sendcount, recvcounts(*), rdispls(*)

if( nodes.eq.1) then
    do i = 1, sendcount
       recvbuf(i) = sendbuf(i)
    enddo
    return
endif

call MPI_ALLGATHERV(sendbuf,sendcount,MPI_DOUBLE_PRECISION,  &
&          recvbuf,recvcounts,rdispls,MPI_DOUBLE_PRECISION,  &
&    world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine dgatherv( sendbuf,sendcount,recvbuf,recvcounts,         &
&                    rdispls,root )
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
real*8  :: sendbuf(*), recvbuf(*)
integer :: sendcount, recvcounts(*), rdispls(*), root

if( nodes.eq.1) then
    do i = 1, sendcount
       recvbuf(i) = sendbuf(i)
    enddo
    return
endif

call MPI_GATHERV(sendbuf,sendcount,MPI_DOUBLE_PRECISION,  &
&          recvbuf,recvcounts,rdispls,MPI_DOUBLE_PRECISION,  &
&    root,world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine allrgatherv( sendbuf,sendcount,recvbuf,recvcounts,      &
&                       rdispls )
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
real    :: sendbuf(*), recvbuf(*)
integer :: sendcount, recvcounts(*), rdispls(*)

if( nodes.eq.1) then
    do i = 1, sendcount
       recvbuf(i) = sendbuf(i)
    enddo
    return
endif

call MPI_ALLGATHERV(sendbuf,sendcount,MPI_REAL,  &
&          recvbuf,recvcounts,rdispls,MPI_REAL,  &
&    world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine alligatherv( sendbuf,sendcount,recvbuf,recvcounts,      &
&                       rdispls )
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
integer :: sendbuf(*), recvbuf(*)
integer :: sendcount, recvcounts(*), rdispls(*)

if( nodes.eq.1) then
    do i = 1, sendcount
       recvbuf(i) = sendbuf(i)
    enddo
    return
endif

call MPI_ALLGATHERV(sendbuf,sendcount,MPI_INTEGER,  &
&          recvbuf,recvcounts,rdispls,MPI_INTEGER,  &
&    world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine alllgatherv( sendbuf,sendcount,recvbuf,recvcounts,      &
&                       rdispls )
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
logical :: sendbuf(*), recvbuf(*)
integer :: sendcount, recvcounts(*), rdispls(*)

if( nodes.eq.1) then
    do i = 1, sendcount
       recvbuf(i) = sendbuf(i)
    enddo
    return
endif

call MPI_ALLGATHERV(sendbuf,sendcount,MPI_LOGICAL,  &
&          recvbuf,recvcounts,rdispls,MPI_LOGICAL,  &
&    world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine dscatterv( sendbuf,sendcounts,displs,recvbuf,           &
&                     recvcount,root )
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'
real*8  :: sendbuf(*), recvbuf(*)
integer :: sendcounts(*), recvcount, displs(*), root

if( nodes.eq.1) then
    do i = 1, recvcount
       recvbuf(i) = sendbuf(i)
    enddo
    return
endif

call MPI_SCATTERV(sendbuf,sendcounts,displs,MPI_DOUBLE_PRECISION,  &
&          recvbuf,recvcount,MPI_DOUBLE_PRECISION,  &
&    root,world_,ierr)

return
end


!-----------------------------------------------------------------------
subroutine gsync
!-----------------------------------------------------------------------
use parallel
include 'mpif.h'

call MPI_BARRIER(world_,ierr)

return
end


!-----------------------------------------------------------------------
double precision function timecnt()
!-----------------------------------------------------------------------
use parallel
implicit real*8 (a-h,o-z)
include 'mpif.h'
#ifdef CRAY_T3E
!-----------------------------------------------------------------------
!     *** for Cray T3E ***
real :: timef
timecnt = 1.d-03 * timef()
#else
#ifdef VPP
!-----------------------------------------------------------------------
!     *** for Fujitsu VPP ***
real*4 :: tarray(2)
ct=etime(tarray)
timecnt=tarray(1)+tarray(2)
#else
!-----------------------------------------------------------------------
timecnt=MPI_WTIME()
#endif
#endif

return
end
