

!-----------------------------------------------------------------------
subroutine set_world( world )
!-----------------------------------------------------------------------
implicit none
integer :: world

return
end


!-----------------------------------------------------------------------
subroutine get_world( world )
!-----------------------------------------------------------------------
implicit none
integer :: world
world = 0
return
end


!----------------------------------------------------------------------c
subroutine start_parallel( nprocs, myid, world )
!----------------------------------------------------------------------c
implicit none
integer :: nprocs, myid, world
nprocs= 1
myid  = 0
world = 0
!      write(*,*) ' iam =',iam,'  of  ',nodes
return
end


!----------------------------------------------------------------------c
subroutine end_parallel(info)
!----------------------------------------------------------------------c
implicit none
integer :: info
return
end


!-----------------------------------------------------------------------
subroutine comm_rank( myid, world )
!-----------------------------------------------------------------------
implicit none
integer :: myid, world
myid  = 0
return
end


!-----------------------------------------------------------------------
subroutine comm_size( nprocs, world )
!-----------------------------------------------------------------------
implicit none
integer :: nprocs, world
nprocs= 1
return
end


!-----------------------------------------------------------------------
subroutine comm_group( world, group )
!-----------------------------------------------------------------------
implicit none
integer :: world, group
group = 0
return
end


!-----------------------------------------------------------------------
subroutine comm_split( oldcomm, color, key, newcomm )
!-----------------------------------------------------------------------
implicit none
integer :: oldcomm, color, key, newcomm
newcomm = 0
return
end


!-----------------------------------------------------------------------
subroutine group_excl( group, extract, extract_ids, subgroup )
!-----------------------------------------------------------------------
implicit none
integer :: group, extract, subgroup
integer, dimension(extract) :: extract_ids
subgroup = 0
return
end


!-----------------------------------------------------------------------
subroutine comm_create( world, group_sub, world_sub )
!-----------------------------------------------------------------------
implicit none
integer :: world, group_sub, world_sub
world_sub = 0
return
end


!-----------------------------------------------------------------------
subroutine group_free( group )
!-----------------------------------------------------------------------
implicit none
integer :: group
return
end


!-----------------------------------------------------------------------
subroutine comm_free( world )
!-----------------------------------------------------------------------
implicit none
integer :: world
return
end


!-----------------------------------------------------------------------
subroutine comm_dup( world, world_sub )
!-----------------------------------------------------------------------
integer :: world, world_sub
world_sub = 0
return
end


!----------------------------------------------------------------------c
subroutine csend(mesg_type,mesg_buff,size,inode,pid)
!----------------------------------------------------------------------c
integer mesg_type,size,pid
character mesg_buff(*)
return
end


!----------------------------------------------------------------------c
subroutine cdsend(mesg_type,mesg_buff,size,inode,pid)
!----------------------------------------------------------------------c
integer mesg_type,size,pid
real*8  mesg_buff(*)
return
end


!----------------------------------------------------------------------c
subroutine crsend(mesg_type,mesg_buff,size,inode,pid)
!----------------------------------------------------------------------c
integer mesg_type,size,pid
real    mesg_buff(*)
return
end


!----------------------------------------------------------------------c
subroutine cisend(mesg_type,mesg_buff,size,inode,pid)
!----------------------------------------------------------------------c
integer mesg_type,size,pid
integer mesg_buff(*)
return
end


!-----------------------------------------------------------------------
subroutine clsend(mesg_tag,mesg_buff,size,inode,pid)
!-----------------------------------------------------------------------
integer :: mesg_tag,size,pid,inode,tid,ierr
logical :: mesg_buff(*)
return
end


!----------------------------------------------------------------------c
subroutine crecv(mesg_type,mesg_buff,size,pid)
!----------------------------------------------------------------------c
integer mesg_type,size,pid
character mesg_buff(*)
return
end


!----------------------------------------------------------------------c
subroutine cdrecv(mesg_type,mesg_buff,size,pid)
!----------------------------------------------------------------------c
integer mesg_type,size,pid
real*8  mesg_buff(*)
return
end


!----------------------------------------------------------------------c
subroutine crrecv(mesg_type,mesg_buff,size,pid)
!----------------------------------------------------------------------c
integer mesg_type,size,pid
real    mesg_buff(*)
return
end


!----------------------------------------------------------------------c
subroutine cirecv(mesg_type,mesg_buff,size,pid)
!----------------------------------------------------------------------c
integer mesg_type,size,pid
integer mesg_buff(*)
return
end


!-----------------------------------------------------------------------
subroutine clrecv(mesg_tag,mesg_buff,size,pid)
!-----------------------------------------------------------------------
integer :: mesg_tag,size,bufid,pid
logical :: mesg_buff(*)
return
end


!----------------------------------------------------------------------c
subroutine gdsum(mesg_buff,n,temp)
!----------------------------------------------------------------------c
real*8 mesg_buff(*),temp(*)
return
end


!----------------------------------------------------------------------c
subroutine dsum(mesg_buff,n,temp,root)
!----------------------------------------------------------------------c
real*8 mesg_buff(*),temp(*)
integer root
return
end


!----------------------------------------------------------------------c
subroutine gisum(mesg_buff,n,temp)
!----------------------------------------------------------------------c
integer mesg_buff(*),temp(*)
return
end


!----------------------------------------------------------------------c
subroutine gdmax(mesg_buff)
!----------------------------------------------------------------------c
real*8 mesg_buff
return
end


!-----------------------------------------------------------------------
subroutine gdmaxpl(mesg_buff,n,temp)
!-----------------------------------------------------------------------
integer :: n
real*8, dimension(n) :: mesg_buff,temp
return
end


!----------------------------------------------------------------------c
subroutine gimax(mesg_buff)
!----------------------------------------------------------------------c
return
end


!----------------------------------------------------------------------c
subroutine gimaxpl(mesg_buff,n,temp)
!----------------------------------------------------------------------c
integer :: n
integer, dimension(n) :: mesg_buff,temp
return
end


!----------------------------------------------------------------------c
subroutine gdmin(mesg_buff)
!----------------------------------------------------------------------c
real*8 mesg_buff(*)
return
end


!-----------------------------------------------------------------------
subroutine gdminpl(mesg_buff,n,temp)
!-----------------------------------------------------------------------
integer :: n
real*8, dimension(n) :: mesg_buff,temp
return
end


!----------------------------------------------------------------------c
subroutine gimin(mesg_buff)
!----------------------------------------------------------------------c
return
end


!-----------------------------------------------------------------------
subroutine dbcast(mesg_buff,n,root)
!-----------------------------------------------------------------------
real*8 mesg_buff(*)
integer root
return
end


!-----------------------------------------------------------------------
subroutine rbcast(mesg_buff,n,root)
!-----------------------------------------------------------------------
real    mesg_buff(*)
integer root
return
end


!-----------------------------------------------------------------------
subroutine ibcast(mesg_buff,n,root)
!-----------------------------------------------------------------------
integer mesg_buff(*)
integer root
return
end


!-----------------------------------------------------------------------
subroutine lbcast(mesg_buff,n,root)
!-----------------------------------------------------------------------
logical :: mesg_buff(*)
integer :: root
return
end


!-----------------------------------------------------------------------
subroutine mespasd(mesg_tag,buffs,buffr,nsd,nrc,inode,pid)
!-----------------------------------------------------------------------
integer mesg_tag, pid
real*8  buffs(*), buffr(*)

do i = 1, nsd
   buffr(i) = buffs(i)
enddo

return
end


!-----------------------------------------------------------------------
subroutine mespasr(mesg_tag,buffs,buffr,nsd,nrc,inode,pid)
!-----------------------------------------------------------------------
integer mesg_tag, pid
real    buffs(*), buffr(*)

do i = 1, nsd
   buffr(i) = buffs(i)
enddo

return
end


!-----------------------------------------------------------------------
subroutine mespasi(mesg_tag,buffs,buffr,nsd,nrc,inode,pid)
!-----------------------------------------------------------------------
integer mesg_tag, pid
integer buffs(*), buffr(*)

do i = 1, nsd
   buffr(i) = buffs(i)
enddo

return
end


!-----------------------------------------------------------------------
subroutine mespasl(mesg_tag,buffs,buffr,nsd,nrc,inode,pid)
!-----------------------------------------------------------------------
integer mesg_tag, pid
logical :: buffs(*), buffr(*)

do i = 1, nsd
   buffr(i) = buffs(i)
enddo

return
end


!----------------------------------------------------------------------c
subroutine gcol(x,msize,x_tot)
!----------------------------------------------------------------------c
character x(*),x_tot(*)
do i=1,msize
  x_tot(i)=x(i)
end do
return
end


!-----------------------------------------------------------------------
subroutine alldgatherv( sendbuf,sendcount,recvbuf,recvcounts,  &
&                       rdispls )
!-----------------------------------------------------------------------
real*8 sendbuf(*), recvbuf(*)
integer sendcount, recvcounts(*), rdispls(*)

do i = 1, sendcount
   recvbuf(i) = sendbuf(i)
enddo

return
end


!-----------------------------------------------------------------------
subroutine dgatherv( sendbuf,sendcount,recvbuf,recvcounts,         &
&                    rdispls,root )
!-----------------------------------------------------------------------
real*8  :: sendbuf(*), recvbuf(*)
integer :: sendcount, recvcounts(*), rdispls(*), root

do i = 1, sendcount
   recvbuf(i) = sendbuf(i)
enddo

return
end


!-----------------------------------------------------------------------
subroutine allrgatherv( sendbuf,sendcount,recvbuf,recvcounts,  &
&                       rdispls )
!-----------------------------------------------------------------------
real    sendbuf(*), recvbuf(*)
integer sendcount, recvcounts(*), rdispls(*)

do i = 1, sendcount
   recvbuf(i) = sendbuf(i)
enddo

return
end


!-----------------------------------------------------------------------
subroutine alligatherv( sendbuf,sendcount,recvbuf,recvcounts,      &
&                       rdispls )
!-----------------------------------------------------------------------
integer :: sendbuf(*), recvbuf(*)
integer :: sendcount, recvcounts(*), rdispls(*)

do i = 1, sendcount
   recvbuf(i) = sendbuf(i)
enddo

return
end


!-----------------------------------------------------------------------
subroutine alllgatherv( sendbuf,sendcount,recvbuf,recvcounts,      &
&                       rdispls )
!-----------------------------------------------------------------------
logical :: sendbuf(*), recvbuf(*)
integer :: sendcount, recvcounts(*), rdispls(*)

do i = 1, sendcount
   recvbuf(i) = sendbuf(i)
enddo

return
end


!-----------------------------------------------------------------------
subroutine dscatterv( sendbuf,sendcounts,displs,recvbuf,  &
&                     recvcount,root )
!-----------------------------------------------------------------------
real*8 sendbuf(*), recvbuf(*)
integer sendcounts(*), recvcount, displs(*), root

do i = 1, recvcount
   recvbuf(i) = sendbuf(i)
enddo

return
end


!----------------------------------------------------------------------c
subroutine gsync
!----------------------------------------------------------------------c
return
end



#if ! IBM_SP3
!----------------------------------------------------------------------c
double precision function timecnt()
!----------------------------------------------------------------------c
implicit real*8 (a-h,o-z)
real   tarray(2), etime
timecnt=etime(tarray)
!      timecnt=dclock ()
return
end
#else
!----------------------------------------------------------------------c
double precision function timecnt()
!----------------------------------------------------------------------c
implicit real*8 (a-h,o-z)
real   tarray(2), etime
timecnt = 0.d0
!      timecnt=dclock ()
return
end
#endif
