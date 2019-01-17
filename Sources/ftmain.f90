



module communicator
!-----------------------------------------------------------------------
! type declaration of variables for communicator
!-----------------------------------------------------------------------
implicit none

integer :: world            ! global, default communicator
integer :: nodes            ! the number of nodes in communicator world
integer :: myid             ! node ID in communicator world

integer :: world_qm         ! communicator
integer :: nodes_qm         ! the number of nodes in communicator world_qm
integer :: myid_qm          ! node ID in communicator world_qm
integer :: id_qm1           ! root ID in communicator world_qm
logical :: lqm_node         ! .true. = QM node

integer :: world_qm_un      ! communicator for unifying over domains
integer :: nodes_qm_un      ! the number of nodes in communicator world_qm_un
integer :: myid_qm_un       ! node ID in communicator world_qm_un

integer :: world_kd         ! communicator for k-point decomposition
integer :: nodes_kd         ! the number of nodes in communicator world_kd
integer :: myid_kd          ! node ID in communicator world_kd
integer :: nkd              ! ID for k-points decomposition
integer :: npk = 1          ! the number of k-point decomposition domains

integer :: world_un         ! communicator for unifying over k-points decomposition domain
integer :: nodes_un         ! the number of nodes in communicator world_un
integer :: myid_un          ! node ID in communicator world_un

integer :: world_pw         ! communicator for plane-wave decomposition
integer :: nodes_pw         ! the number of nodes in communicator world_pw
integer :: myid_pw          ! node ID in communicator world_pw
integer :: nppw = 1         ! the number of plane-wave decomposition domains

integer :: world_lr         ! communicator for linear-response TDDFT
integer :: nodes_lr         ! the number of nodes in communicator world_lr
integer :: myid_lr          ! node ID in communicator world_lr
integer :: nplr = 1         ! the number of linear-response-TDDFT domains

integer :: world_md         ! communicator
integer :: nodes_md         ! the number of nodes in communicator world_md
integer :: myid_md          ! node ID in communicator world_md
integer :: id_md1           ! root ID in communicator world_md
logical :: lmd_node         ! .true. = MD node

logical :: loverlap         ! .true. = MD nodes overlap with QM nodes
logical :: lpureMD          ! .true. = pure classical MD

save

end module




program ftmain
!-----------------------------------------------------------------------
!                                                        since '94/07/31
!                                                       update '09/09/14
!
!     Electronic structure calculation 
!       with the plane-wave        method for solving Khon-Sham equation,
!        and the finite difference method for solving Poisson equation.
!
!     External potential is calculated by pseudopotentials.
!
!         make pwp      for serial calculation
!         make pwpmpi   for MPI
!-----------------------------------------------------------------------
use communicator
implicit none

integer :: ierr


!-----Start parallel environment and keep my node ID
call start_parallel( nodes, myid, world )


!-----sub-main program
call submain


!-----Finalize the parallel environment
call end_parallel( ierr )


stop
end program




subroutine submain
!-----------------------------------------------------------------------
!     sub-main program
!-----------------------------------------------------------------------
use communicator
implicit none

!-----declare local variables
real*8  :: timecnt
real*8  :: ct, ct0
integer :: npx = 1,    npy = 1,    npz = 1
integer :: npx_md = 1, npy_md = 1, npz_md = 1
integer :: iogpsz = 1  ! I/O group size
integer :: nproc, nproc_lr, nproc_pw, nproc_qm, nproc_md
logical :: lreadMDatoms
integer :: i
integer :: date_time(8)
character(10) :: big_ben(3)
integer :: ierror = 0

integer, parameter :: nfmax = 80    !  No. of data files
integer, dimension(nfmax) :: nfile
! =  &
!  &          (/ 06, 10, 11, 12, 13, 14, 15, 16, 17, 18, &
!  &             19, 20, 21, 22, 23, 24, 25, 26, 27, 28, &
!  &             29, 30, 31, 32, 33, 34, 35, 36, 37, 38, &
!  &             39, 40, 41, 42, 43, 44, 45, 46, 47, 48, &
!  &             49, 50, 51, 52, 53, 54, 55, 56, 57, 58 /)
integer, parameter :: nfqm = 31     ! starting number for QM-node IO
                                    ! nfile(  01)~(nfqm-1) for MD nodes
                                    ! nfile(nfqm)~(nfmax ) for QM nodes


!---program starting time
call date_and_time( big_ben(1), big_ben(2), big_ben(3), date_time )

ct  = timecnt()
ct0 = ct

!-----obtain unit numbers
nfile(1) = 6
do i = 2, nfmax
   call allocate_unit_number( nfile(i) )
end do
!if( myid == 0 ) then
!    do i = 1, nfmax
!       write(nfile(1),'(a,i3,a,i3)') '  nfile(',i,' ) = ', nfile(i)
!       write(nfile(2),'(a,i3,a,i3)') '  nfile(',i,' ) = ', nfile(i)
!    end do
!end if


!-----Allocate I/0 files for log & scratch files
!call open_logfiles( nfile, myid, nodes )
call open_md_logfile( nfile, myid, nodes )


if( myid == 0 ) then
    write(nfile(1),'(a,i2.2,a,i2.2,a,i2.2,a,i3.3,a,i4,a,i2.2,a,i2.2)') &
& ' Program starting at ', &
& date_time(5),':',date_time(6),':',date_time(7),'.',date_time(8), &
& ' on ', date_time(1),'/',date_time(2),'/',date_time(3)
    write(nfile(2),'(a,i2.2,a,i2.2,a,i2.2,a,i3.3,a,i4,a,i2.2,a,i2.2)') &
& ' Program starting at ', &
& date_time(5),':',date_time(6),':',date_time(7),'.',date_time(8), &
& ' on ', date_time(1),'/',date_time(2),'/',date_time(3)
    write(nfile(1),*) 'Starting--nodes = ', nodes
    write(nfile(2),*) 'Starting--nodes = ', nodes
end if


!-----get parallel nodes
call get_nps( nfile, myid, nodes,  &
& npx, npy, npz, npk, nppw, nplr, npx_md, npy_md, npz_md, iogpsz, lreadMDatoms )

nproc    = npx*npy*npz
nproc_pw = nproc    * npk
nproc_lr = nproc_pw * nppw
nproc_qm = nproc_lr * nplr
nproc_md = npx_md*npy_md*npz_md
loverlap = nproc_qm + nproc_md > nodes
lpureMD  = nproc_qm <= 0


if( myid == 0 ) then
    write(nfile(1),*) 'QM / MD --nodes = ', nproc_qm, nproc_md
    write(nfile(2),*) 'QM / MD --nodes = ', nproc_qm, nproc_md
end if

!-----error trap
if( nproc_qm > nodes .or. nproc_md > nodes )  &
&   call fstop( nfile, myid, nodes, 'too few nodes' )
if( nproc_md <= 0 )  &
&   call fstop( nfile, myid, nodes, 'MD nodes must be > 0' )


!-----define and create new communicators
call get_newcomm( nfile,  &
& nproc_qm, nproc, nproc_pw, nproc_lr, nproc_md  )


if( lqm_node ) then

    !-----set communicator
    call set_world( world_qm )

    !-----id of k-points-decomposition domain
    nkd = myid_un

    !-----Allocate I/0 files for log & scratch files
    call open_qm_logfile( nfile(nfqm), myid, nodes,  &
& myid_qm, nodes_qm, myid_qm_un, nodes_qm_un, myid_kd, nodes_kd, nkd,  &
& myid_pw, nodes_pw, myid_lr, nodes_lr, ierror )

end if

!-----set communicator
call set_world( world )

!----- error trap
ierror = abs(ierror)
call gimax(ierror)
if( ierror /= 0 )  &
&   call fstop( nfile, myid, nodes, 'cannot open file qm_log' )


!-----setup data
call qmmd_setup( nfile, nfqm,  &
& npx, npy, npz, npx_md, npy_md, npz_md, iogpsz )


!-----main routine
call qmmd_main( nfile, nfqm )


!-----save data
call qmmd_save( nfile, nfqm )




!-----set communicator
call set_world( world )

!------internode syncronization
call gsync


!-----release communicator
if( lqm_node ) then
    call comm_free( world_qm )
    call comm_free( world_kd )
    call comm_free( world_un )
end if

!-----release communicator
call comm_free( world_md )


!---program ending time
call date_and_time( big_ben(1), big_ben(2), big_ben(3), date_time )

ct = timecnt()
if( myid == 0 ) then
    write(nfile(1),*) ' '
    write(nfile(1),'(a,i2.2,a,i2.2,a,i2.2,a,i3.3,a,i4,a,i2.2,a,i2.2)') &
& '  Program ending at ', &
& date_time(5),':',date_time(6),':',date_time(7),'.',date_time(8), &
& ' on ', date_time(1),'/',date_time(2),'/',date_time(3)
    write(nfile(2),*) ' '
    write(nfile(2),'(a,i2.2,a,i2.2,a,i2.2,a,i3.3,a,i4,a,i2.2,a,i2.2)') &
& '  Program ending at ', &
& date_time(5),':',date_time(6),':',date_time(7),'.',date_time(8), &
& ' on ', date_time(1),'/',date_time(2),'/',date_time(3)
    write(nfile(1),1000) ct - ct0
    write(nfile(2),1000) ct - ct0
end if
1000 format(/'  total cpu-time ( sec )  :',f14.6/)


return
end subroutine




subroutine get_newcomm( nfile,  &
& nproc_qm, nproc, nproc_pw, nproc_lr, nproc_md  )
!-----------------------------------------------------------------------
!     define and create new communicators
!-----------------------------------------------------------------------
use communicator
implicit none
integer :: nfile(*)
integer :: nproc, nproc_lr, nproc_pw, nproc_qm, nproc_md

!integer :: world            ! global, default communicator
!integer :: nodes            ! the number of nodes in communicator world
!integer :: myid             ! node ID in communicator world
!
!integer :: world_qm         ! communicator
!integer :: nodes_qm         ! the number of nodes in communicator world_qm
!integer :: myid_qm          ! node ID in communicator world_qm
!integer :: id_qm1, id_qm2   ! start/end node ID in communicator world_qm
!logical :: lqm_node         ! .true. = QM node
!
!integer :: world_kd         ! communicator for k-points decomposition
!integer :: nodes_kd         ! the number of nodes in communicator world_kd
!integer :: myid_kd          ! node ID in communicator world_kd
!
!integer :: world_un         ! communicator for unifying over k-points decomposition domain
!integer :: nodes_un         ! the number of nodes in communicator world_un
!integer :: myid_un          ! node ID in communicator world_un
!
!integer :: world_pw         ! communicator for plane-wave decomposition
!integer :: nodes_pw         ! the number of nodes in communicator world_pw
!integer :: myid_pw          ! node ID in communicator world_pw
!
!integer :: world_lr         ! communicator for linear-response TDDFT
!integer :: nodes_lr         ! the number of nodes in communicator world_lr
!integer :: myid_lr          ! node ID in communicator world_lr
!
!integer :: world_md         ! communicator
!integer :: nodes_md         ! the number of nodes in communicator world_md
!integer :: myid_md          ! node ID in communicator world_md
!integer :: id_md1, id_md2   ! start/end node ID in communicator world_md
!logical :: lmd_node         ! .true. = MD node


!-----declare local variables
integer :: id_qm2   ! end node ID in communicator world_qm
integer :: id_md2   ! end node ID in communicator world_md
!integer, allocatable, dimension(:) :: extract_ids  ! node IDs to be extracted
!integer :: group            ! group of communicator world
!integer :: iworld           ! temporal communicator
!integer :: igroup           ! temporal group
integer :: nextract, nmd
integer :: i
!integer :: npk, nkd, nppw, nplr, npw, nlr, node_s, node_e, ii
!integer :: nun, npwun
integer :: inode, nodemx
integer :: ibuf(9) = 0
integer :: color
integer :: status


!-----extract group from communicator
!call comm_group( world, group )


!!-----allocate memory for extract_ids
!allocate( extract_ids(0:nodes-1), stat=status )
!if( status /= 0 ) call fstop( nfile, myid, nodes,  &
!& 'memory allocation error for extract_ids' )
!
!do i = 0, nodes - 1
!   extract_ids(i) = i
!end do

!-----create new communicator for QM nodes
nextract = nodes - nproc_qm

id_qm1 = nextract
id_qm2 = nodes - 1
lqm_node = myid >= id_qm1 .and. myid <= id_qm2

if( nextract > 0 .and. nextract < nodes ) then

!       !-----create new group by extracting some nodes
!       call group_excl( group, nextract, extract_ids(0), igroup )
!
!       !-----create new communicator from new group
!       call comm_create( world, igroup, world_qm )
!
!       !-----release group
!       call group_free( igroup )
    if( lqm_node ) then
        color = 1
    else
        color = -1
    end if
    call comm_split( world, color, 1, world_qm )

else

   !-----duplicate communicator with all its cached information
   call comm_dup( world, world_qm )

end if


!-----create new communicator for MD nodes
nextract = nodes - nproc_md
id_md1 = 0
id_md2 = nproc_md - 1
lmd_node = myid >= id_md1 .and. myid <= id_md2


if( nextract > 0 ) then

!   !-----create new group by extracting some nodes
!   call group_excl( group, nextract, extract_ids(nproc_md),  &
!&                   igroup )
!
!   !-----create new communicator from new group
!   call comm_create( world, igroup, world_md )
!
!   !-----release group
!   call group_free( igroup )
    if( lmd_node ) then
        color = 1
    else
        color = -1
    end if
    call comm_split( world, color, 1, world_md )

 else

   !-----duplicate communicator with all its cached information
   call comm_dup( world, world_md )

end if


!-----create new communicator for unifying over space decomposition domain

if( nproc_qm > 0 ) then

    if( lqm_node ) then
        nmd = mod(myid-id_qm1,nproc_qm)
    else
        nmd = -1
    end if
    call comm_split( world, nmd, 1, world_qm_un )

else

    !-----duplicate communicator with all its cached information
    call comm_dup( world, world_qm_un )
end if


!-----create new communicator for k-points decomposition
!if( nproc > 0 ) then
!    npk  = nproc_pw/nproc
!    nppw = nproc_lr/nproc_pw
!    nplr = nproc_qm/nproc_lr
!  else
!    npk  = 0
!    nppw = 0
!    nplr = 0
!end if

if( npk*nppw*nplr > 1 ) then

!    nextract = id_qm1 + nproc*(npk-1) + nproc_pw*(nppw-1) + nproc_lr*(nplr-1)
!    do nlr = 0, nplr - 1
!    do npw = 0, nppw - 1
!    do nkd = 0, npk - 1
!
!       node_s = nproc*nkd + nproc_pw*npw + nproc_lr*nlr
!       node_e = node_s + nproc - 1 
!       ii = id_qm1 - 1
!       do i = 0, node_s - 1
!          ii = ii + 1
!          extract_ids(ii) = i + id_qm1
!       end do
!       do i = node_e + 1, nproc_qm - 1
!          ii = ii + 1
!          extract_ids(ii) = i + id_qm1
!       end do
!
!       !-----create new group by extracting some nodes
!       call group_excl( group, nextract, extract_ids(0), igroup )
!
!       !-----create new communicator from new group
!       call comm_create( world, igroup, iworld )
!
!       !-----release group
!       call group_free( igroup )
!
!       if( myid >= node_s + id_qm1 .and.  &
!&          myid <= node_e + id_qm1 ) then
!           world_kd = iworld
!       end if
!
!    end do
!    end do
!    end do
    if( lqm_node ) then
        color = (myid-id_qm1)/nproc + id_qm1
    else
        color = -1
    end if
    call comm_split( world, color, 1, world_kd )

 else

    if( lqm_node ) then
        !-----duplicate communicator with all its cached information
        call comm_dup( world_qm, world_kd )
    end if

end if


!-----create new communicator for unifying over k-points decomposition domain
if( nproc_qm > 1 ) then

!    nextract = id_qm1 + (nproc-1)*npk + nproc_pw*(nppw-1) + nproc_lr*(nplr-1)
!    do nlr = 0, nplr - 1
!    do npw = 0, nppw - 1
!    do nun = 0, nproc - 1
!
!       node_s = nproc_pw*npw + nproc_lr*nlr
!       node_e = node_s + nproc_pw - 1 
!       ii = id_qm1 - 1
!       do i = 0, nproc_qm - 1
!          if( mod(i,nproc) /= nun .or. i < node_s .or. i > node_e ) then
!              ii = ii + 1
!              extract_ids(ii) = i + id_qm1
!          end if
!       end do
!
!       !-----create new group by extracting some nodes
!       call group_excl( group, nextract, extract_ids(0), igroup )
!
!       !-----create new communicator from new group
!       call comm_create( world, igroup, iworld )
!
!       !-----release group
!       call group_free( igroup )
!
!       if( lqm_node .and. mod(myid-id_qm1,nproc) == nun &
!&    .and. myid-id_qm1 >= node_s .and. myid-id_qm1 <= node_e ) then
!           world_un = iworld
!       end if
!
!    end do
!    end do
!    end do
    if( lqm_node ) then
        color = (myid-id_qm1)/nproc_pw
        color = color*nproc + mod(myid-id_qm1,nproc) + id_qm1
    else
        color = -1
    end if
    call comm_split( world, color, 1, world_un )

  else

    if( lqm_node ) then
        !-----duplicate communicator with all its cached information
        call comm_dup( world_qm, world_un )
    end if

end if


!-----create new communicator for plane-wave decomposition domain
if( nproc_qm > 1 ) then

!    nextract = id_qm1 + (nproc_pw-1)*nppw + nproc_lr*(nplr-1)
!    do nlr = 0, nplr - 1
!    do npwun = 0, nproc_pw - 1
!
!       node_s = nproc_lr*nlr
!       node_e = node_s + nproc_lr - 1 
!       ii = id_qm1 - 1
!       do i = 0, nproc_qm - 1
!          if( mod(i,nproc_pw) /= npwun .or. i < node_s .or. i > node_e ) then
!              ii = ii + 1
!              extract_ids(ii) = i + id_qm1
!          end if
!       end do
!
!       !-----create new group by extracting some nodes
!       call group_excl( group, nextract, extract_ids(0), igroup )
!
!       !-----create new communicator from new group
!       call comm_create( world, igroup, iworld )
!
!       !-----release group
!       call group_free( igroup )
!
!       if( lqm_node .and. mod(myid-id_qm1,nproc_pw) == npwun  &
!&    .and. myid-id_qm1 >= node_s .and. myid-id_qm1 <= node_e ) then
!           world_pw = iworld
!       end if
!
!    end do
!    end do
    if( lqm_node ) then
        color = (myid-id_qm1)/nproc_lr
        color = color*nproc_pw + mod(myid-id_qm1,nproc_pw) + id_qm1
    else
        color = -1
    end if
    call comm_split( world, color, 1, world_pw )

  else

    if( lqm_node ) then
        !-----duplicate communicator with all its cached information
        call comm_dup( world_qm, world_pw )
    end if

end if


!-----create new communicator for linear-response TDDFT domain
if( nproc_qm > 1 ) then

!    nextract = id_qm1 + (nproc_lr-1)*nplr
!    do npwun = 0, nproc_lr - 1
!
!       ii = id_qm1 - 1
!       do i = 0, nproc_qm - 1
!          if( mod(i,nproc_lr) /= npwun ) then
!              ii = ii + 1
!              extract_ids(ii) = i + id_qm1
!          end if
!       end do
!
!       !-----create new group by extracting some nodes
!       call group_excl( group, nextract, extract_ids(0), igroup )
!
!       !-----create new communicator from new group
!       call comm_create( world, igroup, iworld )
!
!       !-----release group
!       call group_free( igroup )
!
!       if( lqm_node .and. mod(myid-id_qm1,nproc_lr) == npwun ) then
!           world_lr = iworld
!       end if
!
!    end do
    if( lqm_node ) then
        color = (myid-id_qm1)/nproc_qm
        color = color*nproc_lr + mod(myid-id_qm1,nproc_lr) + id_qm1
    else
        color = -1
    end if
    call comm_split( world, color, 1, world_lr )

  else

    if( lqm_node ) then
        !-----duplicate communicator with all its cached information
        call comm_dup( world_qm, world_lr )
    end if

end if


!-----release group
!call group_free( group )

!!-----deallocate memory for extract_ids
!deallocate( extract_ids, stat=status )


if( lqm_node ) then
    call comm_size( nodes_qm, world_qm )
    call comm_rank( myid_qm,  world_qm )

    call comm_size( nodes_kd, world_kd )
    call comm_rank( myid_kd,  world_kd )

    call comm_size( nodes_un, world_un )
    call comm_rank( myid_un,  world_un )

    call comm_size( nodes_pw, world_pw )
    call comm_rank( myid_pw,  world_pw )

    call comm_size( nodes_lr, world_lr )
    call comm_rank( myid_lr,  world_lr )

    call comm_size( nodes_qm_un, world_qm_un )
    call comm_rank( myid_qm_un,  world_qm_un )

!    write(*,*) 'iam, qm=', myid, myid_qm, nodes_qm, world_qm
end if

if( lmd_node ) then
    call comm_size( nodes_md, world_md )
    call comm_rank( myid_md,  world_md )

!    write(*,*) 'iam, md=', myid, myid_md, nodes_md, world_md
end if


!--- check communicator setup
if( lmd_node ) then
    ibuf(1) = 1
    ibuf(3) = myid_md
end if
if( lqm_node ) then
    ibuf(2) = 1
    ibuf(4) = myid_qm
    ibuf(5) = myid_qm_un
    ibuf(6) = myid_kd
    ibuf(7) = myid_un
    ibuf(8) = myid_pw
    ibuf(9) = myid_lr
end if

nodemx = 195
if( myid == 0 ) then

    write(nfile(1),'(1x,a)')  &
& '    myid myid_md myid_qm myid_qm_un myid_kd myid_un myid_pw myid_lr'
    write(nfile(2),'(1x,a)')  &
& '    myid myid_md myid_qm myid_qm_un myid_kd myid_un myid_pw myid_lr'

!    do inode = 0, nodes - 1
    do inode = 0, min( nodes - 1, nodemx )
       if( inode /= myid ) then
           call cirecv(10,ibuf,2,0)
           if( ibuf(1) == 1 ) call cirecv(11,ibuf(3),1,0)
           if( ibuf(2) == 1 ) call cirecv(12,ibuf(4),6,0)
       end if
       do i = 1, 2
          if( ibuf(1) == 1 .and. ibuf(2) == 1 ) then
              write(nfile(i),'(4i8,i11,6i8)') inode, ibuf(3:9)
          else if( ibuf(1) == 1 .and. ibuf(2) == 0 ) then
              write(nfile(i),'(7i8)') inode, ibuf(3)
          else if( ibuf(1) == 0 .and. ibuf(2) == 1 ) then
              write(nfile(i),'(i8,8x,2i8,i11,6i8)') inode, ibuf(4:9)
          end if
       end do
       !------internode syncronization
       call gsync
    end do
    if( nodes - 1 > nodemx ) then
        write(nfile(1),'(1x,a)') '    ......'
        write(nfile(2),'(1x,a)') '    ......'
    end if

else

!    do inode = 0, nodes - 1
    do inode = 0, min( nodes - 1, nodemx )
       if( inode == myid ) then
           call cisend(10,ibuf,2,0,0)
           if( ibuf(1) == 1 ) call cisend(11,ibuf(3),1,0,0)
           if( ibuf(2) == 1 ) call cisend(12,ibuf(4),6,0,0)
       end if
       !------internode syncronization
       call gsync
    end do

end if

!-----Finalize the parallel environment
!call end_parallel( status )
!stop


return
end subroutine




subroutine get_nps( nfile, myid, nodes,  &
& npx, npy, npz, npk, nppw, nplr, npx_md, npy_md, npz_md, iogpsz, lreadMDatoms )
!-----------------------------------------------------------------------
!     get parallel nodes from an input file
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: npx, npy, npz, npk, nppw, nplr, npx_md, npy_md, npz_md, iogpsz
logical :: lreadMDatoms

!-----declare local variables
character(10) :: word
character(50) :: fname1
character(50) :: fname = 'control/filename'
integer :: iunit
integer :: ierror, istat


call allocate_unit_number( iunit )

!--- get file name : fname1
if( myid == 0 ) then
    write(nfile(1),*)  &
&   'open file(in get_nps): ', fname(1:len_trim(fname))
    write(nfile(2),*)  &
&   'open file(in get_nps): ', fname(1:len_trim(fname))
end if
open( iunit, file=fname, status='old', action='read',  &
&     iostat=ierror )

!--- global maximum
call gimax(ierror)

blockif_a1: if( ierror == 0 ) then

   read(iunit,*,iostat=istat) fname1
   !----- error trap
   if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                   'error in file :'//fname )

else blockif_a1

   !----- error trap
   call fstop( nfile, myid, nodes, 'cannot open file :'//fname )

end if blockif_a1
close(iunit)


!--- open input file
if( myid == 0 ) then
    write(nfile(1),*)  &
&   'open file(in get_nps): ', fname1(1:len_trim(fname1))
    write(nfile(2),*)  &
&   'open file(in get_nps): ', fname1(1:len_trim(fname1))
end if
open( iunit, file=fname1, status='old', action='read',  &
&     iostat=ierror )

!--- global maximum
call gimax(ierror)

blockif_b1: if( ierror == 0 ) then

   lreadMDatoms = .false.

   !----- read input data
   readdo: do
      read(iunit,'(a10)',iostat=istat) word
      blockif_b2: if( istat == 0 ) then

         readif: if( word == '*parallel ' ) then

            call read_prll( nfile, myid, nodes, iunit,  &
& npx, npy, npz, npk, nppw, nplr, npx_md, npy_md, npz_md, iogpsz )

         else if( word == '*MD atoms ' ) then readif

            lreadMDatoms = .true.

         end if readif

      else if( istat < 0 ) then blockif_b2

         !----- exit, if at EOF
         exit

      else blockif_b2

         !----- error trap
         call fstop( nfile, myid, nodes,'error in file :'//fname1)

      end if blockif_b2

   end do readdo

else blockif_b1

   !----- error trap
   call fstop( nfile, myid, nodes, 'cannot open file :'//fname1 )

end if blockif_b1
close(iunit)

call deallocate_unit_number( iunit )


return
end subroutine




subroutine read_prll( nfile, myid, nodes, iunit,  &
& npx, npy, npz, npk, nppw, nplr, npx_md, npy_md, npz_md, iogpsz )
!-----------------------------------------------------------------------
! read variables for parallel nodes
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
integer :: npx, npy, npz, npk, nppw, nplr, npx_md, npy_md, npz_md, iogpsz
character(10) :: word
integer :: istat


readdo: do
   read(iunit,'(a10)',iostat=istat) word

   readif: if( istat == 0 ) then

      selectif: if( word == '(QM-nodes)' ) then

         read(iunit,*,iostat=istat) npx, npy, npz
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0001 in read_prll' )

      else if( word == '(k-points)' ) then selectif

         read(iunit,*,iostat=istat) npk
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0002 in read_prll' )

      else if( word == '(plane wav' ) then selectif

         read(iunit,*,iostat=istat) nppw
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0003 in read_prll' )

      else if( word == '(linear-re' ) then selectif

         read(iunit,*,iostat=istat) nplr
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0004 in read_prll' )

      else if( word == '(MD-nodes)' ) then selectif

         read(iunit,*,iostat=istat) npx_md, npy_md, npz_md
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0005 in read_prll' )

      else if( word == '(I/O group' ) then selectif

         read(iunit,*,iostat=istat) iogpsz
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,          &
&                         'error-0006 in read_prll' )

      else if( word == '*end      ' ) then selectif

         exit

      endif selectif

   else readif

      !----- error trap
      call fstop( nfile, myid, nodes, 'error-0000 in read_prll' )

   end if readif
end do readdo


return
end subroutine




subroutine qmmd_setup( nfile, nfqm,  &
& npx, npy, npz, npx_md, npy_md, npz_md, iogpsz )
!-----------------------------------------------------------------------
!     setup data
!-----------------------------------------------------------------------
use communicator
implicit none
integer :: nfile(*)
integer :: nfqm
integer :: npx, npy, npz
integer :: npx_md, npy_md, npz_md, iogpsz

!-----declare local variables
integer :: ierror
real*8  :: ct, ct0, timecnt


ct0 = timecnt()


ierror = 0
!-----initial settting for MD nodes
if( lmd_node ) then

    !-----set communicator
    call set_world( world_md )

    !------ MD calculation 
    call remd_set( nfile, myid_md, nodes_md,  &
& npx_md, npy_md, npz_md, iogpsz, lpureMD, ierror )

end if


!-----if MD nodes overlap with QM nodes, all QM nodes must wait
if( loverlap ) then
    !-----set communicator
    call set_world( world )

    !------internode syncronization
    !call gsync

    !-----error trap
    call gimax(ierror)
end if


!-----initial settting for QM nodes
if( lqm_node .and. ierror == 0 ) then

    !-----set communicator
    call set_world( world_kd )

    !------ LDA/GGA calculation 
    call pwlda_set( nfile(nfqm), myid_kd, nodes_kd,  &
& npx, npy, npz, myid_qm_un, npx_md, npy_md, npz_md, iogpsz, ierror )

    call pwlda_set2( nfile(nfqm), myid_kd, nodes_kd, ierror )

end if


!-----set communicator
call set_world( world )

!-----error trap
call gimax(ierror)
if( ierror /= 0 ) then
    if( myid == 0 ) then
        write(nfile(1),*) 'error in initial setting:', ierror
        write(nfile(2),*) 'error in initial setting:', ierror
    end if
    call fstop( nfile, myid, nodes, ' ' )
end if


!-----send coordinates for QM atoms in MD nodes to QM nodes
!-----to determine displacement vectors, d = r_in_QM_nodes - r_in_MD_nodes
if( .not.lpureMD ) then
    call send_positions( nfile, myid, nodes, nfqm, world,  &
& myid_qm, nodes_qm, lqm_node, world_qm, id_qm1,  &
& myid_md, nodes_md, lmd_node, world_md, id_md1, loverlap, .false. )
end if

!===stop here===
!call fstop( nfile, myid, nodes, 'stop in qmmd_setup' )
!===============


ct = timecnt()
if( myid.eq.0 ) then
    write(nfile(1),*) ' setup data ... done : cpu-time :', ct- ct0
    write(nfile(2),*) ' setup data ... done : cpu-time :', ct- ct0
end if


return
end subroutine




subroutine qmmd_main( nfile, nfqm )
!-----------------------------------------------------------------------
!     send coordinates for QM atoms in MD nodes to QM nodes
!-----------------------------------------------------------------------
use communicator
implicit none
integer :: nfile(*)
integer :: nfqm

!-----declare local variables
integer :: nstep, nstop, ifmd, nstep_ini
integer :: ierror
real*8  :: ct, ct0, timecnt
integer :: imd, isccnt
real*8  :: prevene, ctmd0
logical :: lvsend = .false.
logical :: ltimecnt = .true.
logical :: lmdstop, lmdstoprun
integer :: iunit_stop, idmd
character(50) :: fstopname = 'STOP'


ltimecnt = .not.lpureMD

!-----set communicator
call set_world( world )
idmd = id_md1
!-----unify basic parameters
if( myid == idmd ) then
    call set_param_md(  nfile, myid, nodes,  &
& nstep, nstop, ifmd, nstep_ini, lmdstop )
end if

call ibcast(nstep,1,idmd)
call ibcast(nstop,1,idmd)
call ibcast(ifmd,1,idmd)
call ibcast(nstep_ini,1,idmd)
call lbcast(lmdstop,1,idmd)

ierror = 0
!-----check basic parameters
if( lmd_node ) then
    call check_param_in_md( nfile, myid, nodes,  &
& nstep, nstop, ifmd, ierror )
end if
if( lqm_node .and. ierror == 0 ) then
    call check_param_in_qm( nfile, myid, nodes,  &
& nstep, nstop, ifmd, ierror )
end if


!-----error trap
call gimax(ierror)
if( ierror /= 0 ) then
    if( myid == 0 ) then
        write(nfile(1),*) 'error in setting parameters:', ierror
        write(nfile(2),*) 'error in setting parameters:', ierror
    end if
    call fstop( nfile, myid, nodes, ' ' )
end if


!-----check if velocities need to be sent
if( .not.lpureMD ) then
    if( myid == id_qm1 ) then
        call set_lvsend(  nfile, myid, nodes, lvsend )
    end if
    call lbcast(lvsend,1,id_qm1)
end if


!-----create a STOP file, if necessary
lmdstoprun = .false.
if( myid == 0 ) then
    call allocate_unit_number( iunit_stop )
    if( lmdstop ) then
        open( iunit_stop, file=fstopname, iostat=ierror )
        write(iunit_stop,*) lmdstoprun
        close( iunit_stop )
    else
        open( iunit_stop, file=fstopname, iostat=ierror )
        close(iunit_stop,status='delete')
        call deallocate_unit_number( iunit_stop )
    end if
end if


!-----start MD calculations
ct0 = timecnt()
isccnt = 0
prevene = 1.d+10
ctmd0  = ct0
if( ifmd == 0 ) nstop = 1
if( ifmd >= 1 .and. nstep == nstep_ini ) then
    nstep = nstep - 1
    nstop = nstop + 1
end if
mddo: do imd = 1, nstop
!-----------------------------------------------------------------------
   if( ifmd >= 1 ) nstep = nstep + 1
   !-----set time step
   if( lmd_node ) then
       call set_nstep_in_md( nfile, myid, nodes, nstep )
   end if
   if( lqm_node ) then
       call set_nstep_in_qm( nfile, myid, nodes, nstep )
   end if

   !-----update atomic configuration
   if( lmd_node ) then

       !-----set communicator
       call set_world( world_md )

       !------ MD calculation 
       call remd_updtconfig( nfile, myid_md, nodes_md, ierror )

   end if


   !-----error trap
   call check_error_in_MD( nfile, myid, nodes, world, ierror )

   if( ierror < 0 ) then
       ierror = 0
       return
   end if


   !-----send coordinates of QM atoms in MD nodes to QM nodes
   if( .not.lpureMD ) then
       call send_positions( nfile, myid, nodes, nfqm, world,  &
& myid_qm, nodes_qm, lqm_node, world_qm, id_qm1,  &
& myid_md, nodes_md, lmd_node, world_md, id_md1, loverlap, .true. )
   end if


   if( ltimecnt ) then
       ct = timecnt()
       if( myid == 0 ) then
           write(nfile(1),*) ' '
           write(nfile(2),*) ' '
           write(nfile(1),*) '                   setup coordinates',  &
&                            ' : cpu-time :', ct - ct0
           write(nfile(2),*) '                   setup coordinates',  &
&                            ' : cpu-time :', ct - ct0
       end if
       ct0 = ct
   end if


   !-----calculate classical forces
   if( lmd_node ) then

       !-----set communicator
       call set_world( world_md )

       !------ MD calculation 
       call remd_force( nfile, myid_md, nodes_md )

   end if


   if( ltimecnt ) then
       ct = timecnt()
       if( myid == 0 ) then
           write(nfile(1),*) '          classical forces  ... done',  &
&                            ' : cpu-time :', ct - ct0
           write(nfile(2),*) '          classical forces  ... done',  &
&                            ' : cpu-time :', ct - ct0
           write(nfile(1),*) ' MD nodes wait for QM forces.'
!           write(nfile(2),*) ' MD nodes wait for QM forces.'
       end if
       ct0 = ct
   end if


   !-----if MD nodes overlap with QM nodes, all QM nodes must wait
   if( loverlap ) then
       !-----set communicator
       call set_world( world )

       !------internode syncronization
       call gsync
   end if


   !-----calculate QM forces
   if( lqm_node ) then

       !-----set communicator
       call set_world( world_kd )

       !------ LDA/GGA calculation 
       call pwlda( nfile(nfqm), myid_kd, nodes_kd, ierror )

   end if


   !-----error trap
   call check_error_in_MD( nfile, myid, nodes, world, ierror )


   !-----send forces for QM atoms from QM nodes to MD nodes
   if( .not.lpureMD ) then
       call send_forces( nfile, myid, nodes, world,  &
& myid_qm, nodes_qm, lqm_node, world_qm, id_qm1,  &
& myid_md, nodes_md, lmd_node, world_md, id_md1, loverlap )

       !-----send velocities back to MD nodes, if necessary
       if( lvsend ) call send_velocities_to_MD( nfile, myid, nodes, world,  &
& myid_qm, nodes_qm, lqm_node, world_qm, id_qm1,  &
& myid_md, nodes_md, lmd_node, world_md, id_md1, loverlap )
   end if


   if( ltimecnt ) then
       ct = timecnt()
       if( myid == 0 ) then
           write(nfile(1),*) '                 QM forces  ... done',  &
&                            ' : cpu-time :', ct - ct0
           write(nfile(2),*) '                 QM forces  ... done',  &
&                            ' : cpu-time :', ct - ct0
       end if
       ct0 = ct
   end if


   !-----calculate forces & energy, output data to files
   if( lmd_node ) then

       !-----set communicator
       call set_world( world_md )

       !------ MD calculation 
       call remd_semifinal( nfile, myid_md, nodes_md, ierror )

   end if


   !-----send velocities of QM atoms in MD nodes to QM nodes, if necessary
   if( .not.lpureMD ) then
       if( lvsend )  &
& call send_velocities_to_QM( nfile, myid, nodes, nfqm, world,  &
& myid_qm, nodes_qm, lqm_node, world_qm, id_qm1,  &
& myid_md, nodes_md, lmd_node, world_md, id_md1, loverlap )
   end if


   !-----error trap
   call check_error_in_MD( nfile, myid, nodes, world, ierror )


   if( ifmd >= 1 ) then
       if( ltimecnt ) then
           ct = timecnt()
           if( myid.eq.0 ) then
           write(nfile(1),*) ' '
           write(nfile(2),*) ' '
           write(nfile(1),*) '                             MD step',  &
&                            ' : cpu-time :', ct - ctmd0
           write(nfile(2),*) '                             MD step',  &
&                            ' : cpu-time :', ct - ctmd0
           end if
           ct0   = ct
           ctmd0 = ct
       end if
   endif


   if( ierror < 0 ) then
       ierror = 0
       return
   end if


   !-----check the STOP file, if necessary
   if( lmdstop ) then
       !-----set communicator
       call set_world( world )
       if( myid == 0 ) then
!           rewind(iunit_stop)
           open( iunit_stop, file=fstopname, iostat=ierror )
           if( ierror == 0 ) read(iunit_stop,*,iostat=ierror) lmdstoprun
           if( ierror /= 0 ) lmdstoprun = .false.
           close( iunit_stop )
       end if
       call lbcast(lmdstoprun,1,0)

       !---if lmdstoprun = .true., stop MD run
       if( lmdstoprun ) exit
   end if
!-----------------------------------------------------------------------
end do mddo


return
end subroutine




subroutine check_error_in_MD( nfile, myid, nodes, world, ierror )
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*)
integer :: world            ! global, default communicator
integer :: nodes            ! the number of nodes in communicator world
integer :: myid             ! node ID in communicator world
integer :: ierror
integer :: ierror_max


!-----set communicator
call set_world( world )

!-----error trap
ierror_max = ierror
call gimax(ierror_max)
if( ierror_max > 0 ) then
    if( myid == 0 ) then
        write(nfile(1),*) 'error in MD routines:', ierror_max
        write(nfile(2),*) 'error in MD routines:', ierror_max
    end if
    call fstop( nfile, myid, nodes, ' ' )
end if


!-----check ierror: stop the calculation if ierror < 0
call gimin(ierror)


return
end subroutine




subroutine send_positions( nfile, myid, nodes, nfqm, world,  &
& myid_qm, nodes_qm, lqm_node, world_qm, id_qm1,  &
& myid_md, nodes_md, lmd_node, world_md, id_md1, loverlap, lupdate )
!-----------------------------------------------------------------------
!     send coordinates of QM atoms in MD nodes to QM nodes
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*)
integer :: nfqm

integer :: world            ! global, default communicator
integer :: nodes            ! the number of nodes in communicator world
integer :: myid             ! node ID in communicator world

integer :: world_qm         ! communicator
integer :: nodes_qm         ! the number of nodes in communicator world_qm
integer :: myid_qm          ! node ID in communicator world_qm
integer :: id_qm1           ! root node ID in communicator world_qm
logical :: lqm_node         ! .true. = QM node

integer :: world_md         ! communicator
integer :: nodes_md         ! the number of nodes in communicator world_md
integer :: myid_md          ! node ID in communicator world_md
integer :: id_md1           ! root node ID in communicator world_md
logical :: lmd_node         ! .true. = MD node

logical :: loverlap
logical :: lupdate

!-----declare local variables
integer :: ierror


ierror = 0
if( lmd_node ) then

    !-----set communicator
    call set_world( world_md )

    !------ MD calculation 
    call remd_setQMpositions( nfile, myid_md, nodes_md )

end if


!-----set communicator
call set_world( world )

!-----if MD nodes overlap with QM nodes, all QM nodes must wait
if( loverlap ) then
    !------internode syncronization
    call gsync
end if

!-----send coordinates to QM nodes
if( id_qm1 == id_md1 ) then

    !-----if root nodes are the same, do not send but just copy data
    if( myid == id_md1 ) then
        call copyQMpositions( nfile, myid, nodes, ierror )
    end if

  else

    if( myid == id_md1 ) then
        !-----send data to QM nodes
        call sendQMpositions( nfile, myid, nodes, id_qm1, ierror )
    else if( myid == id_qm1 ) then
        !-----recieve data from MD nodes
        call recvQMpositions( nfile, myid, nodes, ierror )
    end if

end if

!-----set communicator
call set_world( world )

!-----error trap
call gimax(ierror)
if( ierror /= 0 )  &
&   call fstop( nfile, myid, nodes, 'data transfer error in send_positions' )


!-----unify transferred coordinates in QM nodes
if( lqm_node ) then

    !-----set communicator
    call set_world( world_qm )

    !------ LDA/GGA calculation 
    call unifQMpositions( nfile(nfqm), myid_qm, nodes_qm,  &
& lupdate, ierror )

end if


!-----set communicator
call set_world( world )

!-----error trap
call gimax(ierror)
if( ierror /= 0 )  &
&   call fstop( nfile, myid, nodes, 'data transfer error(2) in send_positions' )


return
end subroutine




subroutine send_velocities_to_QM( nfile, myid, nodes, nfqm, world,  &
& myid_qm, nodes_qm, lqm_node, world_qm, id_qm1,  &
& myid_md, nodes_md, lmd_node, world_md, id_md1, loverlap )
!-----------------------------------------------------------------------
!     send velocities of QM atoms in MD nodes to QM nodes
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*)
integer :: nfqm

integer :: world            ! global, default communicator
integer :: nodes            ! the number of nodes in communicator world
integer :: myid             ! node ID in communicator world

integer :: world_qm         ! communicator
integer :: nodes_qm         ! the number of nodes in communicator world_qm
integer :: myid_qm          ! node ID in communicator world_qm
integer :: id_qm1           ! root node ID in communicator world_qm
logical :: lqm_node         ! .true. = QM node

integer :: world_md         ! communicator
integer :: nodes_md         ! the number of nodes in communicator world_md
integer :: myid_md          ! node ID in communicator world_md
integer :: id_md1           ! root node ID in communicator world_md
logical :: lmd_node         ! .true. = MD node

logical :: loverlap

!-----declare local variables
integer :: ierror


ierror = 0
if( lmd_node ) then

    !-----set communicator
    call set_world( world_md )

    !------ MD calculation 
    call remd_setQMvelocities( nfile, myid_md, nodes_md )

end if


!-----set communicator
call set_world( world )

!-----if MD nodes overlap with QM nodes, all QM nodes must wait
if( loverlap ) then
    !------internode syncronization
    call gsync
end if

!-----send velocities to QM nodes
if( id_qm1 == id_md1 ) then

    !-----if root nodes are the same, do not send but just copy data
    if( myid == id_md1 )  &
&       call copyQMvelocities( nfile, myid, nodes, ierror )

  else

    if( myid == id_md1 ) then
        !-----send data to QM nodes
        call sendQMvelocities( nfile, myid, nodes, id_qm1, ierror )
    else if( myid == id_qm1 ) then
        !-----recieve data from MD nodes
        call recvQMvelocities( nfile, myid, nodes, ierror )
    end if

end if

!-----error trap
call gimax(ierror)
if( ierror /= 0 )  &
&   call fstop( nfile, myid, nodes, 'velocity transfer error' )


!-----unify transferred velocities in QM nodes
if( lqm_node ) then

    !-----set communicator
    call set_world( world_qm )

    !------ LDA/GGA calculation 
    call unifQMvelocities( nfile(nfqm), myid_qm, nodes_qm, ierror )

end if


!-----set communicator
call set_world( world )

!-----error trap
call gimax(ierror)
if( ierror /= 0 )  &
&   call fstop( nfile, myid, nodes, 'velocity transfer error(2)' )


return
end subroutine




subroutine send_forces( nfile, myid, nodes, world,  &
& myid_qm, nodes_qm, lqm_node, world_qm, id_qm1,  &
& myid_md, nodes_md, lmd_node, world_md, id_md1, loverlap )
!-----------------------------------------------------------------------
!     send forces for QM atoms from QM nodes to MD nodes
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*)

integer :: world            ! global, default communicator
integer :: nodes            ! the number of nodes in communicator world
integer :: myid             ! node ID in communicator world

integer :: world_qm         ! communicator
integer :: nodes_qm         ! the number of nodes in communicator world_qm
integer :: myid_qm          ! node ID in communicator world_qm
integer :: id_qm1           ! root node ID in communicator world_qm
logical :: lqm_node         ! .true. = QM node

integer :: world_md         ! communicator
integer :: nodes_md         ! the number of nodes in communicator world_md
integer :: myid_md          ! node ID in communicator world_md
integer :: id_md1           ! root node ID in communicator world_md
logical :: lmd_node         ! .true. = MD node

logical :: loverlap

!-----declare local variables
integer :: ierror


ierror = 0
!-----set communicator
call set_world( world )

!-----send forces to MD nodes
if( id_qm1 == id_md1 ) then

    !-----if root nodes are the same, do not send but just copy data
    if( myid == id_qm1 )  &
&       call copyQMforces( nfile, myid, nodes, ierror )

  else

    if( myid == id_qm1 ) then
        !-----send data to MD nodes
        call sendQMforces( nfile, myid, nodes, id_md1, ierror )
    else if( myid == id_md1 ) then
        !-----recieve data from QM nodes
        call recvQMforces( nfile, myid, nodes, ierror )
    end if

end if

!-----error trap
call gimax(ierror)
if( ierror /= 0 )  &
&   call fstop( nfile, myid, nodes, 'force transfer error' )


!-----unify transferred forces in MD nodes
if( lmd_node ) then

    !-----set communicator
    call set_world( world_md )

    !------ MD calculation 
    call unifQMforces( nfile, myid_md, nodes_md, ierror )

end if


!-----set communicator
call set_world( world )

!-----error trap
call gimax(ierror)
if( ierror /= 0 )  &
&   call fstop( nfile, myid, nodes, 'force transfer error(2)' )


return
end subroutine




subroutine send_velocities_to_MD( nfile, myid, nodes, world,  &
& myid_qm, nodes_qm, lqm_node, world_qm, id_qm1,  &
& myid_md, nodes_md, lmd_node, world_md, id_md1, loverlap )
!-----------------------------------------------------------------------
!     send velocities back to MD nodes
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*)

integer :: world            ! global, default communicator
integer :: nodes            ! the number of nodes in communicator world
integer :: myid             ! node ID in communicator world

integer :: world_qm         ! communicator
integer :: nodes_qm         ! the number of nodes in communicator world_qm
integer :: myid_qm          ! node ID in communicator world_qm
integer :: id_qm1           ! root node ID in communicator world_qm
logical :: lqm_node         ! .true. = QM node

integer :: world_md         ! communicator
integer :: nodes_md         ! the number of nodes in communicator world_md
integer :: myid_md          ! node ID in communicator world_md
integer :: id_md1           ! root node ID in communicator world_md
logical :: lmd_node         ! .true. = MD node

logical :: loverlap

!-----declare local variables
integer :: ierror


ierror = 0
!-----set communicator
call set_world( world )

!-----send velocities to MD nodes
if( id_qm1 == id_md1 ) then

    !-----if root nodes are the same, do not send but just copy data
    if( myid == id_qm1 )  &
&       call copy_velocities_toMD( nfile, myid, nodes, ierror )

  else

    if( myid == id_qm1 ) then
        !-----send data to MD nodes
        call send_velocities_toMD( nfile, myid, nodes, id_md1, ierror )
    else if( myid == id_md1 ) then
        !-----recieve data from QM nodes
        call recv_velocities_fromQM( nfile, myid, nodes, ierror )
    end if

end if

!-----error trap
call gimax(ierror)
if( ierror /= 0 )  &
&   call fstop( nfile, myid, nodes, 'velocity transfer error' )


!-----unify transferred forces in MD nodes
if( lmd_node ) then

    !-----set communicator
    call set_world( world_md )

    !------ MD calculation 
    call unify_velocities_inMD( nfile, myid_md, nodes_md, ierror )

end if


!-----set communicator
call set_world( world )

!-----error trap
call gimax(ierror)
if( ierror /= 0 )  &
&   call fstop( nfile, myid, nodes, 'velocity transfer error(2)' )


return
end subroutine




subroutine qmmd_save( nfile, nfqm )
!-----------------------------------------------------------------------
!     save data
!-----------------------------------------------------------------------
use communicator
implicit none
integer :: nfile(*)
integer :: nfqm

!-----declare local variables
integer :: ierror = 0
real*8  :: ct, ct0, timecnt


ct0 = timecnt()

!-----save data in MD nodes
if( lmd_node ) then

    !-----set communicator
    call set_world( world_md )

    !------ MD calculation 
    call remd_save( nfile, myid_md, nodes_md, ierror )

end if


!-----if MD nodes overlap with QM nodes, all QM nodes must wait
if( loverlap ) then
    !-----set communicator
    call set_world( world )

    !------internode syncronization
    call gsync
end if


!-----save data in QM nodes
if( lqm_node .and. ierror == 0 ) then

    !-----set communicator
    call set_world( world_kd )

    !------ LDA/GGA calculation 
    call pwlda_save( nfile(nfqm), myid_kd, nodes_kd, ierror )

end if


!-----set communicator
call set_world( world )

!-----error trap
call gimax(ierror)
if( ierror /= 0 ) then
    if( myid == 0 ) then
        write(nfile(1),*) 'error while saving data:', ierror
        write(nfile(2),*) 'error while saving data:', ierror
    end if
    call fstop( nfile, myid, nodes, ' ' )
end if


ct = timecnt()
if( myid.eq.0 ) then
    write(nfile(1),*) ' dump data ... done  : cpu-time :', ct- ct0
    write(nfile(2),*) ' dump data ... done  : cpu-time :', ct- ct0
endif


return
end subroutine




subroutine get_gworld( myid_, nodes_ )
!-----------------------------------------------------------------------
!     set world
!-----------------------------------------------------------------------
use communicator
implicit none
integer :: myid_, nodes_


!-----set communicator
call set_world( world )

myid_  = myid
nodes_ = nodes


return
end subroutine




subroutine get_worldqm( myid_, nodes_ )
!-----------------------------------------------------------------------
!     set world_qm
!-----------------------------------------------------------------------
use communicator
implicit none
integer :: myid_, nodes_


!-----set communicator
call set_world( world_qm )

myid_  = myid_qm
nodes_ = nodes_qm


return
end subroutine




subroutine get_worldqmun( myid_, nodes_ )
!-----------------------------------------------------------------------
!     set world_qm_un
!-----------------------------------------------------------------------
use communicator
implicit none
integer :: myid_, nodes_


!-----set communicator
call set_world( world_qm_un )

myid_  = myid_qm_un
nodes_ = nodes_qm_un


return
end subroutine



subroutine get_worldkd( myid_, nodes_, nkd_ )
!-----------------------------------------------------------------------
!     set world_kd
!-----------------------------------------------------------------------
use communicator
implicit none
integer :: myid_, nodes_, nkd_


!-----set communicator
call set_world( world_kd )

myid_  = myid_kd
nodes_ = nodes_kd
nkd_   = nkd


return
end subroutine




subroutine get_worldun( myid_, nodes_ )
!-----------------------------------------------------------------------
!     set world_un
!-----------------------------------------------------------------------
use communicator
implicit none
integer :: myid_, nodes_


!-----set communicator
call set_world( world_un )

myid_  = myid_un
nodes_ = nodes_un


return
end subroutine




subroutine get_worldpw( myid_, nodes_ )
!-----------------------------------------------------------------------
!     set world_pw
!-----------------------------------------------------------------------
use communicator
implicit none
integer :: myid_, nodes_


!-----set communicator
call set_world( world_pw )

myid_  = myid_pw
nodes_ = nodes_pw


return
end subroutine




subroutine get_worldlr( myid_, nodes_ )
!-----------------------------------------------------------------------
!     set world_pw
!-----------------------------------------------------------------------
use communicator
implicit none
integer :: myid_, nodes_


!-----set communicator
call set_world( world_lr )

myid_  = myid_lr
nodes_ = nodes_lr


return
end subroutine




subroutine get_worldmd( myid_, nodes_ )
!-----------------------------------------------------------------------
!     set world_md
!-----------------------------------------------------------------------
use communicator
implicit none
integer :: myid_, nodes_


!-----set communicator
call set_world( world_md )

myid_  = myid_md
nodes_ = nodes_md


return
end subroutine




subroutine qm_gsync
!-----------------------------------------------------------------------
!     synchronization in world_qm
!-----------------------------------------------------------------------
use communicator
implicit none
integer :: current_world


!-----get current communicator
call get_world( current_world )

!-----set communicator
call set_world( world_qm )

call gsync

!-----set communicator
call set_world( current_world )


return
end subroutine
