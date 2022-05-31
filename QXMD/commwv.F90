



module commwv_variables
!-----------------------------------------------------------------------
! type declaration for planewave decomposition
!-----------------------------------------------------------------------
implicit none

logical :: lnoncollinear

save

end module




module planewave_decomp_variables
!-----------------------------------------------------------------------
! type declaration for planewave decomposition
!-----------------------------------------------------------------------
implicit none

integer :: myid_pw, nodes_pw
integer :: nplw5nod1, nplw5nod2, nplw5nod
integer, allocatable, dimension(:) :: nplw5cnt, nplw5dsp
integer, allocatable, dimension(:) :: nplw5excnt, nplw5exdsp

save

end module




subroutine set_lnoncollinear_in_commwv( nfile, myid, nodes, lnoncollinear_ )
!-----------------------------------------------------------------------
!     allocate memory for planewave-decomposition variables
!-----------------------------------------------------------------------
use commwv_variables
implicit none
integer :: nfile(*), myid, nodes
logical :: lnoncollinear_

lnoncollinear = lnoncollinear_

return
end subroutine




subroutine planewave_decomp_variables_alloc( nfile, alloc_mem )
!-----------------------------------------------------------------------
!     allocate memory for planewave-decomposition variables
!-----------------------------------------------------------------------
use planewave_decomp_variables
implicit none
integer, dimension(*) :: nfile
real*8  :: alloc_mem

!-----declare local variables
integer :: myid, nodes, nkd
integer :: status
real*8  :: the_mem
integer :: nproc


!-----set communicator
call get_worldpw( myid, nodes )

myid_pw  = myid
nodes_pw = nodes
nproc    = nodes
!------allocate memory
allocate( nplw5cnt(nproc), nplw5dsp(nproc), &
& nplw5excnt(nproc), nplw5exdsp(nproc), &
& stat=status )

the_mem = 4.d0 * nproc * 2

!-----set communicator
call get_worldkd( myid, nodes, nkd )

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'planewave_decomp_variables_alloc', .true. )


return
end subroutine




subroutine planewave_decomp( nfile, &
& nplw5ex, nplw5 )
!-----------------------------------------------------------------------
!     set indexes for planewave decompositions
!-----------------------------------------------------------------------
use planewave_decomp_variables
implicit none
integer, dimension(*) :: nfile
integer :: nplw5ex, nplw5

!------declare local variables
integer :: myid, nodes
integer :: inode, i, nkd


!-----set communicator
call get_worldpw( myid, nodes )

do inode = 0, nodes - 1
   call divnod( nplw5+1, nodes, inode, nplw5nod1, nplw5nod2, nplw5nod )
   nplw5cnt(inode+1) = nplw5nod
end do
nplw5dsp(1) = 0
do i = 2, nodes
   nplw5dsp(i) = nplw5dsp(i-1) + nplw5cnt(i-1)
end do
do i = 1, nodes
   nplw5excnt(i) = 2 * nplw5cnt(i)
   nplw5exdsp(i) = 2 * nplw5dsp(i)
end do

call divnod( nplw5+1, nodes, myid, nplw5nod1, nplw5nod2, nplw5nod )

!-----set communicator
call get_worldkd( myid, nodes, nkd )


return
end subroutine




subroutine band_decomp( nfile, myid, nodes,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, nbxxxx,  &
& iod, iodg, iods, iodsg )
!-----------------------------------------------------------------------
!     set indexes for band decompositions
!-----------------------------------------------------------------------
implicit none

integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: nband, nbnod1, nbnod2, nbnod
integer, dimension(nodes) :: nbncnt, nbndsp
integer :: nbxxxx
integer :: iod(nbnod), iodg(nband), iods(nbnod), iodsg(nband)

!------declare local variables
integer :: inode, i


!--- index for band decomposition --------------------------------------
!--- total No. of bands        : nband
!--- No. of bands in each node : nbnod
!call decompose( nfile, myid, nodes,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp )


!-----assign bands to each node
call asgnbds( nfile, myid, nodes,  &
& iod, iodg, iods, iodsg, nband, nbnod1, nbnod2, nbnod,  &
& nbncnt, nbndsp )


#if PCHOLESKY
nbxxxx = nbnod
#else
nbxxxx = nband
#endif


return
end subroutine




subroutine atom_decomp( nfile, myid, nodes,  &
& natom, natnod1, natnod2, natnod, natcnt, natdsp, iatoit,  &
& ntype, nhk1, nhk2, nhk1_nat, nhk2_nat, ioa, ioag, ioas, ioasg )
!-----------------------------------------------------------------------
!     set indexes for atom decompositions
!-----------------------------------------------------------------------
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: natom, natnod1, natnod2, natnod, natcnt(nodes), natdsp(nodes)
integer :: iatoit(natom)
integer :: ntype
integer :: nhk1(ntype), nhk2(ntype), nhk1_nat(ntype), nhk2_nat(ntype)
integer :: ioa(natnod), ioag(natom), ioas(natnod), ioasg(natom)

!------declare local variables
integer :: it, ia, iaa, i


!--- index for atom decomposition --------------------------------------
!--- total No. of atoms        : natom
!--- No. of atoms in each node : natnod
!call decompose( nfile, myid, nodes,  &
!& natom, natnod1, natnod2, natnod, natcnt, natdsp )


!-----assign bands to each node
call asgnbds( nfile, myid, nodes,  &
& ioa, ioag, ioas, ioasg, natom, natnod1, natnod2, natnod,  &
& natcnt, natdsp )

!!---for check
!do i = 1, natnod
!   ioa(i) = ioas(i)
!enddo
!do i = 1, natom
!   ioag(i) = i
!enddo

!-----set nhk1_nod, nhk2_nod
do it = 1, ntype
   nhk2_nat(it) = 0
end do
do iaa = 1, natnod
   ia = ioa(iaa)
   it = iatoit(ia)
   nhk2_nat(it) = nhk2_nat(it) + 1
end do
!do it = 1, ntype
!   do ia = nhk1(it), nhk2(it)
!      if( ia >= natnod1 .and. ia <= natnod2 )  &
!&         nhk2_nat(it) = nhk2_nat(it) + 1
!   end do
!end do
nhk1_nat(1) = 1
nhk2_nat(1) = nhk1_nat(1) + nhk2_nat(1) - 1
do it = 2, ntype
   nhk1_nat(it) = nhk2_nat(it-1) + 1
   nhk2_nat(it) = nhk1_nat(it) + nhk2_nat(it) - 1
end do


return
end subroutine




subroutine decompose( nfile, myid, nodes,  &
& nvar, nod1, nod2, nod, ncnt, ndsp )
!-----------------------------------------------------------------------
!     set indexes for decompositions
!-----------------------------------------------------------------------
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: nvar, nod1, nod2, nod
integer, dimension(nodes) :: ncnt, ndsp

!------declare local variables
integer :: inode, i


!--- index for decomposition --------------------------------------
!--- total No. of variables        : nvar
!--- No. of variables in each node : nod
do inode = 0, nodes - 1
   call divnod( nvar, nodes, inode, nod1, nod2, nod )
   ncnt(inode+1) = nod
enddo
ndsp(1) = 0
do i = 2, nodes
   ndsp(i) = ndsp(i-1) + ncnt(i-1)
enddo

call divnod( nvar, nodes, myid, nod1, nod2, nod )


return
end subroutine




subroutine index_decomp( nfile, myid, nodes,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& nplwex, nplw, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& ncgcnt, ncgdsp, nspnod, lnoncollinear )
!-----------------------------------------------------------------------
!     set indexes for G decompositions
!-----------------------------------------------------------------------
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: nband, nbnod1, nbnod2, nbnod
integer, dimension(nodes) :: nbncnt, nbndsp
integer :: nplwex, nplw, npnod1, npnod2, npnod
integer, dimension(nodes) :: nplcnt, npldsp
integer, dimension(nodes) :: ncgcnt, ncgdsp
integer :: nspnod
logical :: lnoncollinear

!------declare local variables
integer :: nplwex_, nplw_


if( .not.lnoncollinear ) then
    call index_decomp2( nfile, myid, nodes,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& nplwex, nplw, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& ncgcnt, ncgdsp, nspnod )
else
    nplwex_ = nplwex * 2
    nplw_   = nplw + nplw + 1
    call index_decomp2( nfile, myid, nodes,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& nplwex_, nplw_, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& ncgcnt, ncgdsp, nspnod )
end if


return
end subroutine




subroutine index_decomp2( nfile, myid, nodes,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& nplwex, nplw, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& ncgcnt, ncgdsp, nspnod )
!-----------------------------------------------------------------------
!     set indexes for G decompositions
!-----------------------------------------------------------------------
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: nband, nbnod1, nbnod2, nbnod
integer, dimension(nodes) :: nbncnt, nbndsp
integer :: nplwex, nplw, npnod1, npnod2, npnod
integer, dimension(nodes) :: nplcnt, npldsp
integer, dimension(nodes) :: ncgcnt, ncgdsp
integer :: nspnod

!------declare local variables
integer :: inode, i


!--- index for G decomposition ----------------------------------------
!--- total No. of components        : nplw+1
!--- No. of components in each node : npnod
do inode = 0, nodes - 1
   call divnod( nplw+1, nodes, inode, npnod1, npnod2, npnod )
   nplcnt(inode+1) = npnod
end do
npldsp(1) = 0
do i = 2, nodes
   npldsp(i) = npldsp(i-1) + nplcnt(i-1)
end do

call divnod( nplw+1, nodes, myid, npnod1, npnod2, npnod )


!--- index for unifying wavefunctions ----------------------------------
do inode = 0, nodes - 1
   ncgcnt(inode+1) = nplwex*nbncnt(inode+1)
end do
ncgdsp(1) = 0
do i = 2, nodes
   ncgdsp(i) = ncgdsp(i-1) + ncgcnt(i-1)
end do

nspnod = max( nplwex*nbnod, 2*npnod*nband )


return
end subroutine




subroutine index_decomp7( nfile, myid, nodes,  &
& nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& nplw7ex, nplw7, npnod71, npnod72, npnod7, nplcnt7, npldsp7,  &
& ncgcnt7, ncgdsp7, nspnod7, lnoncollinear )
!-----------------------------------------------------------------------
!     set indexes for G decompositions related to nplw7
!-----------------------------------------------------------------------
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: nplw
integer :: nband, nbnod1, nbnod2, nbnod
integer, dimension(nodes) :: nbncnt, nbndsp
integer :: nplw7ex, nplw7, npnod71, npnod72, npnod7
integer, dimension(nodes) :: nplcnt7, npldsp7
integer, dimension(nodes) :: ncgcnt7, ncgdsp7
integer :: nspnod7
logical :: lnoncollinear

!------declare local variables
integer :: nplw_, nplw7_, nplw7ex_, i


if( .not.lnoncollinear ) then
    call index_decomp72( nfile, myid, nodes,  &
& nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& nplw7ex, nplw7, npnod71, npnod72, npnod7, nplcnt7, npldsp7,  &
& ncgcnt7, ncgdsp7, nspnod7 )
else
!    nplwex_  = nplwex * 2
    nplw_    = nplw + nplw + 1
!    nplw7_   = nplw7 + nplw7 + 1
    nplw7_   = nplw7
    nplw7ex_ = nplw7ex * 2
    !---set index for myid < nodes/2
    call index_decomp72( nfile, myid, nodes,  &
& nplw_, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& nplw7ex_, nplw7_, npnod71, npnod72, npnod7, nplcnt7, npldsp7,  &
& ncgcnt7, ncgdsp7, nspnod7 )
    if( nodes > 1 ) then
        do i = 1, nodes/2
           nplcnt7(i+nodes/2) = nplcnt7(i)
        end do
        do i = nodes/2+1, nodes
           npldsp7(i) = npldsp7(i-1) + nplcnt7(i-1)
        end do
        if( myid >= nodes/2 ) then
            npnod7  = nplcnt7(myid+1)
            npnod71 = npldsp7(myid+1) + 1
            npnod72 = npldsp7(myid+1) + nplcnt7(myid+1)
            nspnod7 = max( nspnod7, 2*npnod7*nband )
        end if
    else
        !---Caution!! Special care is needed for nodes = 1
        nplcnt7(1) = nplcnt7(1) * 2
        npnod7  = nplcnt7(1)
        npnod71 = 1
        npnod72 = nplcnt7(1)
    end if
end if


return
end subroutine




subroutine index_decomp72( nfile, myid, nodes,  &
& nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& nplw7ex, nplw7, npnod71, npnod72, npnod7, nplcnt7, npldsp7,  &
& ncgcnt7, ncgdsp7, nspnod7 )
!-----------------------------------------------------------------------
!     set indexes for G decompositions related to nplw7
!-----------------------------------------------------------------------
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: nplw
integer :: nband, nbnod1, nbnod2, nbnod
integer, dimension(nodes) :: nbncnt, nbndsp
integer :: nplw7ex, nplw7, npnod71, npnod72, npnod7
integer, dimension(nodes) :: nplcnt7, npldsp7
integer, dimension(nodes) :: ncgcnt7, ncgdsp7
integer :: nspnod7

!------declare local variables
integer :: inode, i
integer :: npnod1, npnod2, npnod, nplsum, nplw71


nplw71 = nplw7 + 1

!--- index for G decomposition ----------------------------------------
!--- total No. of components        : nplw
!--- No. of components in each node : npnod
nplsum = 0
do inode = 0, nodes - 1
   call divnod( nplw+1, nodes, inode, npnod1, npnod2, npnod )
   nplsum = nplsum + npnod
   if( nplsum <= nplw71 ) then
       nplcnt7(inode+1) = npnod
     else
       if( nplsum - nplw71 < npnod ) then
           nplcnt7(inode+1) = npnod - ( nplsum - nplw71 )
         else
           nplcnt7(inode+1) = 0
       end if
   end if
end do
npldsp7(1) = 0
do i = 2, nodes
   npldsp7(i) = npldsp7(i-1) + nplcnt7(i-1)
end do

npnod7  = nplcnt7(myid+1)
npnod71 = npldsp7(myid+1) + 1
npnod72 = npldsp7(myid+1) + nplcnt7(myid+1)


!--- index for unifying wavefunctions ----------------------------------
do inode = 0, nodes - 1
   ncgcnt7(inode+1) = nplw7ex*nbncnt(inode+1)
end do
ncgdsp7(1) = 0
do i = 2, nodes
   ncgdsp7(i) = ncgdsp7(i-1) + ncgcnt7(i-1)
end do

nspnod7 = max( nplw7ex*nbnod, 2*npnod7*nband, 1 )


return
end subroutine




subroutine divnod( ndata, node, myid, ndnod1, ndnod2, ndnod )
!-----------------------------------------------------------------------
! ( input )
!     ndata ...... total No. of data points
!     node  ...... No. of nodes
!     myid  ...... my ID
!
! ( output )
!     ndnod1 ...... starting pointer
!     ndnod2 ...... ending   pointer
!     ndnod  ...... No. of data points
!-----------------------------------------------------------------------

!      ndnod  = ndata/node
!      ndnod1 = ndnod*myid
!      msrest = mod(ndata,node) - 1
!      if( myid.le.msrest ) then
!          ndnod1 = ndnod1 + myid
!          ndnod  = ndnod + 1
!        else
!          ndnod1 = ndnod1 + msrest + 1
!      endif
!      ndnod1 = ndnod1 + 1
!      ndnod2 = ndnod1 + ndnod - 1

ndnod2  = (myid+1)*ndata/node + 1
if( mod((myid+1)*ndata,node).eq.0 ) ndnod2 = ndnod2 - 1
if( myid.eq.0 ) then
    ndnod1 = 1
  else
    ndnod1  = myid*ndata/node + 2
    if( mod(myid*ndata,node).eq.0 ) ndnod1 = ndnod1 - 1
endif
ndnod = ndnod2 - ndnod1 + 1

return
end




subroutine fftmsh( nfile, myid, nodes,  &
& mftnod, mftdsp, mfd2ft, ntotfd, mftwrk, noddatx,  &
& mshglb, kfft1, kfft2, kfft3,  &
& mshnx, mshny, mshnz, mshnod,  &
& mshx1, mshy1, mshz1, mshx, mshy, mshz )
!-----------------------------------------------------------------------
!   set mesh for FFT
!-----------------------------------------------------------------------
use outfile
implicit real*8 ( a-h, o-z )
dimension nfile(*)
dimension mftnod(*), mftdsp(*), mfd2ft(*)
dimension mftwrk(noddatx,3), mshglb(kfft1,kfft2,kfft3)
!!!      dimension mshgnx(*), mshgny(*), mshgnz(*)
dimension mshnx(mshnod), mshny(mshnod), mshnz(mshnod)
dimension ibuf(3)


ct0 = timecnt()

kfft0 = kfft1*kfft2*kfft3
do iz = 1, kfft3
do iy = 1, kfft2
do ix = 1, kfft1
   mshglb(ix,iy,iz) = 0
enddo
enddo
enddo
if( myid.eq.0 ) then
    mftnod(1) = mshnod
    ibuf(1)   = mshx1
    ibuf(2)   = mshy1
    ibuf(3)   = mshz1
    ntotfd = 0
    do m1 = 1, mftnod(1)
       ix = mshnx(m1)
       ix = ibuf(1) + ix - 1
       iy = mshny(m1)
       iy = ibuf(2) + iy - 1
       iz = mshnz(m1)
       iz = ibuf(3) + iz - 1
!             mfd2ft(m1+ntotfd) = ix + kfft1*( (iy-1) + kfft2*(iz-1) )
       mshglb(ix,iy,iz) = m1+ntotfd
    enddo
    do node_id = 1, nodes - 1
       call cirecv(100,mftnod(node_id+1),1,0)
       call cirecv(110,ibuf,3,0)
       call cirecv(120,mftwrk(1,1),mftnod(node_id+1),0)
       call cirecv(130,mftwrk(1,2),mftnod(node_id+1),0)
       call cirecv(140,mftwrk(1,3),mftnod(node_id+1),0)
       ntotfd = ntotfd + mftnod(node_id)
       do m1 = 1, mftnod(node_id+1)
          ix = mftwrk(m1,1)
          ix = ibuf(1) + ix - 1
          iy = mftwrk(m1,2)
          iy = ibuf(2) + iy - 1
          iz = mftwrk(m1,3)
          iz = ibuf(3) + iz - 1
!                mfd2ft(m1+ntotfd) = ix + kfft1*( (iy-1) + kfft2*(iz-1) )
          mshglb(ix,iy,iz) = m1+ntotfd
       enddo
       call gsync
    enddo
  else
    do node_id = 1, nodes - 1
       if( myid.eq.node_id ) then
           call cisend(100,mshnod,1,0,0)
           ibuf(1) = mshx1
           ibuf(2) = mshy1
           ibuf(3) = mshz1
           call cisend(110,ibuf,3,0,0)
           call cisend(120,mshnx,mshnod,0,0)
           call cisend(130,mshny,mshnod,0,0)
           call cisend(140,mshnz,mshnod,0,0)
       endif
       call gsync
    enddo
endif

call ibcast(mftnod,nodes,0)

mftdsp(1) = 0
do i = 2, nodes
   mftdsp(i) = mftdsp(i-1) + mftnod(i-1)
enddo

ntotfd = 0
do i = 1, nodes
   ntotfd = ntotfd + mftnod(i)
enddo

!      call ibcast(mfd2ft,ntotfd,0)
call ibcast(mshglb,kfft0,0)
do iz = 1, kfft3
do iy = 1, kfft2
do ix = 1, kfft1
   m1 = mshglb(ix,iy,iz)
   if( m1.ge.1 ) then
       mfd2ft(m1) = ix + kfft1*( (iy-1) + kfft2*(iz-1) )
!!!             mshgnx(m1) = ix
!!!             mshgny(m1) = iy
!!!             mshgnz(m1) = iz
   endif
enddo
enddo
enddo

ct = timecnt()
do i = 1, 2
if( loutfile(i) ) then
    write(nfile(i),*) ' No. of global FD meshes :', ntotfd
    write(nfile(i),*) ' set mesh points (FD->FFT) : cpu-time :',  &
&                      ct - ct0
end if
end do
ct0 = ct


return
end




subroutine asgnbds( nfile, myid, nodes,  &
& iod, iodg, iods, iodsg, nband, nbnod1, nbnod2, nbnod,  &
& nbncnt, nbndsp )
!-----------------------------------------------------------------------
!    assign bands to each node
!
!     iod,  iodg  ...... for solving eigenvalue problem
!     iods, iodsg ...... for Unitary transformation
!
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension nfile(*)
integer :: iod(nbnod), iodg(nband), iods(nbnod), iodsg(nband)
integer :: nbncnt(nodes), nbndsp(nodes)


nbdiv = nband/nodes
do i = 1, nbdiv
   if( mod(i,2).eq.1 ) then
       iod(i) = nodes*( i - 1 ) + myid + 1
     else
       iod(i) = nodes*i - myid
   endif
enddo

ibi = nband
do inode = 0, nodes - 1
   do i = nbdiv + 1, nbncnt(inode+1)
      if( myid.eq.inode ) iod(i) = ibi
      ibi = ibi - 1
   enddo
enddo

ii = 0
ibi = nband
do inode = 0, nodes - 1
   do i = 1, nbdiv
      ii = ii + 1
      if( mod(i,2).eq.1 ) then
          iodg(ii) = nodes*( i - 1 ) + inode + 1
        else
          iodg(ii) = nodes*i - inode
      endif
   enddo
   do i = nbdiv + 1, nbncnt(inode+1)
      ii = ii + 1
      iodg(ii) = ibi
      ibi = ibi - 1
   enddo
enddo


do i = nbnod1, nbnod2
   iods(i-nbnod1+1) = i
enddo
do i = 1, nband
   iodsg(i) = i
enddo


return
end




module all2all_out
!-----------------------------------------------------------------------
! type declaration of shared variables in commwv.f
!-----------------------------------------------------------------------
implicit none

logical :: la2aout = .true.
real*8,  dimension(10) :: tc = 0.d0

save

end module




subroutine set_all2all_out( la2aout_ )
!-----------------------------------------------------------------------
use all2all_out
logical :: la2aout_

la2aout = la2aout_

return
end




subroutine unifywv( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, idstnd, ncgcnt, ncgdsp, nspnod, npnnbnx, ltimecnt )
!-----------------------------------------------------------------------
!    unify wavefunctions decomposed by components
!
! (input)
!    gdcr : wavefunctions decomposed by components
!
! (output)
!    cgjr : unified wavefunctions
!    gdcr : wavefunctions decomposed by bands
!-----------------------------------------------------------------------
use commwv_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: ct0
integer :: nplwex, nplw
integer :: nband, nbnod1, nbnod2, nbnod, nbncnt(nodes), nbndsp(nodes)
integer :: npnod1, npnod2, npnod, nplcnt(nodes), npldsp(nodes)
integer :: nspnod, npnnbnx
#if CCREAL4
real*4  :: cgjr(*)
real*4  :: gdcr(nspnod)
real*4  :: bufcr(npnnbnx)
#else
real*8  :: cgjr(*)
real*8  :: gdcr(nspnod)
real*8  :: bufcr(npnnbnx)
#endif
integer :: iod(nbnod), iodg(nband), idstnd(nodes)
integer :: ncgcnt(nodes), ncgdsp(nodes)
logical :: ltimecnt


if( .not.lnoncollinear ) then
    call unifywv2( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, idstnd, ncgcnt, ncgdsp, nspnod, npnnbnx, ltimecnt )
else
    !-----noncollinear magnetism
    call unifywv2( nfile, myid, nodes, ct0,  &
& cgjr, nplwex*2, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, idstnd, ncgcnt, ncgdsp, nspnod, npnnbnx, ltimecnt )
end if


return
end




subroutine unifywv2( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, idstnd, ncgcnt, ncgdsp, nspnod, npnnbnx, ltimecnt )
!-----------------------------------------------------------------------
!    unify wavefunctions decomposed by components
!
! (input)
!    gdcr : wavefunctions decomposed by components
!
! (output)
!    cgjr : unified wavefunctions
!    gdcr : wavefunctions decomposed by bands
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension nfile(*)
#if CCREAL4
real*4    cgjr(nplwex,nband)
real*4    gdcr(nspnod)
real*4    bufcr(npnnbnx)
#else
dimension cgjr(nplwex,nband)
dimension gdcr(nspnod)
dimension bufcr(npnnbnx)
#endif
dimension iod(nbnod), iodg(nband), nbndsp(nodes), nbncnt(nodes)
dimension nplcnt(nodes), npldsp(nodes)
dimension ncgcnt(nodes), ncgdsp(nodes)
dimension idstnd(nodes)
logical ltimecnt
save iunify
data iunify / 3 /

if( iunify.eq.1 ) then
!--- collect all data to node 0, and broadcast to all nodes
    call unif1wv( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, idstnd, nspnod, npnnbnx, ltimecnt )
else if( iunify.eq.2 ) then
!--- all to all communication 
    call unif2wv( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, idstnd, nspnod, npnnbnx, ltimecnt )
else
!--- by alldgather
    call unif3wv( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, idstnd, ncgcnt, ncgdsp, nspnod, npnnbnx, ltimecnt )
endif

return
end




subroutine unif1wv( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, idstnd, nspnod, npnnbnx, ltimecnt )
!-----------------------------------------------------------------------
!    unify wavefunctions decomposed by components  << version 1 >>
!
!        collect all data to node 0, and broadcast to all nodes
!
! (input)
!    gdcr : wavefunctions decomposed by components
!
! (output)
!    cgjr : unified wavefunctions
!-----------------------------------------------------------------------
use outfile
implicit real*8 ( a-h, o-z )
dimension nfile(*)
#if CCREAL4
real*4    cgjr(nplwex,nband)
real*4    gdcr(*)
real*4    bufcr(*)
#else
dimension cgjr(nplwex,nband)
dimension gdcr(*)
dimension bufcr(*)
#endif
dimension iod(*), iodg(*), nbndsp(*), nbncnt(*)
dimension nplcnt(*), npldsp(*)
dimension idstnd(*)
logical ltimecnt


if( myid.eq.0 ) then
    i1 = 2*npnod1 - 1
    i2 = 2*npnod2
    i3 = i2 - i1 + 1
    do n = 1, nband
       n0 = i3*( n - 1 )
       do i = i1, i2
          cgjr(i,n) = gdcr(i-i1+1+n0)
       enddo
    enddo
    do node_id = 1, nodes - 1
       i1 = 2*( npldsp(node_id+1) + 1 ) - 1
       i2 = 2*( npldsp(node_id+1) + nplcnt(node_id+1) )
       i3 = i2 - i1 + 1
       nrc = i3*nband
#if CCREAL4
       call crrecv(100,gdcr,nrc,0)
#else
       call cdrecv(100,gdcr,nrc,0)
#endif
       do n = 1, nband
          n0 = i3*( n - 1 )
          do i = i1, i2
             cgjr(i,n) = gdcr(i-i1+1+n0)
          enddo
       enddo
       call gsync
    enddo
  else
    do node_id = 1, nodes - 1
       if( myid.eq.node_id ) then
           i1 = 2*npnod1 - 1
           i2 = 2*npnod2
           i3 = i2 - i1 + 1
           nrc = i3*nband
#if CCREAL4
           call crsend(100,gdcr,nrc,0,0)
#else
           call cdsend(100,gdcr,nrc,0,0)
#endif
       endif
       call gsync
    enddo
endif

#if CCREAL4
call rbcast(cgjr,nplwex*nband,0)
#else
call dbcast(cgjr,nplwex*nband,0)
#endif

if( ltimecnt ) then
    ct = timecnt()
    do i = 1, 2
    if( loutfile(i) ) then
#if CCREAL4
        write(nfile(i),*) '                         unify w.f. ',  &
&                         ' : cpu-time :', ct-ct0
#else
        write(nfile(i),*) '                  real*4 unify w.f. ',  &
&                         ' : cpu-time :', ct-ct0
#endif
    end if
    end do
    ct0 = ct
end if


return
end




subroutine unif2wv( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, idstnd, nspnod, npnnbnx, ltimecnt )
!-----------------------------------------------------------------------
!    unify wavefunctions decomposed by components  << version 2 >>
!
!    all to all communication 
!
! (input)
!    gdcr : wavefunctions decomposed by components
!
! (output)
!    cgjr : unified wavefunctions
!-----------------------------------------------------------------------
use outfile
implicit real*8 ( a-h, o-z )
dimension nfile(*)
#if CCREAL4
real*4    cgjr(nplwex,nband)
real*4    gdcr(*)
real*4    bufcr(*)
#else
dimension cgjr(nplwex,nband)
dimension gdcr(*)
dimension bufcr(*)
#endif
dimension iod(*), iodg(*), nbndsp(*), nbncnt(*)
dimension nplcnt(*), npldsp(*)
dimension idstnd(*)
logical ltimecnt


!    --- to convert G decomposition to band decomposition ---
#if CCREAL4
call rgdtobd( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, idstnd, nspnod, npnnbnx, ltimecnt )
#else
call gdtobd2( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, idstnd, ltimecnt )
#endif


!--- within myself
ib0  = nbndsp(myid+1)
call unif21( nfile, myid, nodes,  &
& cgjr(1,ib0+1), nplwex, nbnod, gdcr )


do isnd = 1, nodes - 1
!--- destination node
   idst = idstnd(isnd)

!--- set No. of sending data (components)
   is1 = 1
   is2 = nplwex
   is3 = is2 - is1 + 1

!--- set No. of recieving data (bands)
   ibs3 = nbncnt(idst+1)
   ib0  = nbndsp(idst+1)

   call unif22( nfile, myid, nodes, isnd,  &
& is1, is2, is3, ibs3, idst, cgjr(1,ib0+1), nplwex, nbnod, gdcr )

!---------Internode synchronization
!         call gsync
enddo


if( ltimecnt ) then
    ct = timecnt()
    do i = 1, 2
    if( loutfile(i) ) then
#if CCREAL4
        write(nfile(i),*) '          unify w.f. (by all to all)',  &
&                         ' : cpu-time :', ct-ct0
#else
        write(nfile(i),*) '   real*4 unify w.f. (by all to all)',  &
&                         ' : cpu-time :', ct-ct0
#endif
    end if
    end do
    ct0 = ct
end if


return
end




subroutine unif21( nfile, myid, nodes, cgjr, nplwex, nbnod, gdcr )
!-----------------------------------------------------------------------
!    unify wavefunctions decomposed by components  << version 2 >>
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension nfile(*)
#if CCREAL4
real*4    cgjr(nplwex,*)
real*4    gdcr(nplwex,*)
#else
dimension cgjr(nplwex,*)
dimension gdcr(nplwex,*)
#endif

do i = 1, nbnod
   do ig = 1, nplwex
      cgjr(ig,i) = gdcr(ig,i)
   enddo
enddo

return
end




subroutine unif22( nfile, myid, nodes, isnd,  &
& is1, is2, is3, ibs3, idst, cgjr, nplwex, nbnod, gdcr )
!-----------------------------------------------------------------------
!    all to all communication to unify wavefunctions
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension nfile(*)
#if CCREAL4
real*4    cgjr(nplwex,*)
real*4    gdcr(nplwex,*)
#else
dimension cgjr(nplwex,*)
dimension gdcr(nplwex,*)
#endif


!--- sending data : gdcr

nsd = is3*nbnod
nrc = nplwex*ibs3
#if CCREAL4
call mespasr(120+isnd,gdcr,cgjr,nsd,nrc,idst,0)  
#else
call mespasd(120+isnd,gdcr,cgjr,nsd,nrc,idst,0)  
#endif
!c---------myid.lt.idst: send & recv, if not empty
!          if (myid.lt.idst) then
!            if (nsd.ne.0) call cdsend(120,gdcr,nsd,idst,0)
!            if (nrc.ne.0) call cdrecv(130,cgjr,nrc,0)
!c---------myid.gt.idst: recv & send, if not empty
!          else
!            if (nrc.ne.0) call cdrecv(120,cgjr,nrc,0)
!            if (nsd.ne.0) call cdsend(130,gdcr,nsd,idst,0)
!          endif

return
end




subroutine unif3wv( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, idstnd, ncgcnt, ncgdsp, nspnod, npnnbnx, ltimecnt )
!-----------------------------------------------------------------------
!    unify wavefunctions decomposed by components  << version 3 >>
!
!    by alldgather
!
! (input)
!    gdcr : wavefunctions decomposed by components
!
! (output)
!    cgjr : unified wavefunctions
!-----------------------------------------------------------------------
use outfile
use all2all_out
implicit real*8 ( a-h, o-z )
dimension nfile(*)
#if CCREAL4
real*4    cgjr(nplwex,nband)
real*4    gdcr(nspnod)
real*4    bufcr(npnnbnx)
#else
dimension cgjr(nplwex,nband)
dimension gdcr(nspnod)
dimension bufcr(npnnbnx)
#endif
dimension iod(nbnod), iodg(nband), nbndsp(nodes), nbncnt(nodes)
dimension nplcnt(nodes), npldsp(nodes)
dimension ncgcnt(nodes), ncgdsp(nodes)
dimension idstnd(nodes)
logical ltimecnt


!    --- to convert G decomposition to band decomposition ---
#if CCREAL4
call rgdtobd( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, idstnd, nspnod, npnnbnx, ltimecnt )

call allrgatherv(gdcr,nplwex*nbnod,cgjr,ncgcnt,ncgdsp)

if( ltimecnt ) then
    ct = timecnt()
    tc(4) = tc(4) + ct - ct0
    if( la2aout ) then
        if(loutfile(1)) write(nfile(1),*) '  real*4 unify w.f. (by alldgatherv)',  &
&                         ' : cpu-time :', tc(4)
        if(loutfile(2)) write(nfile(2),*) '  real*4 unify w.f. (by alldgatherv)',  &
&                         ' : cpu-time :', tc(4)
        tc(4) = 0.d0
    end if
    ct0 = ct
endif
#else
call gdtobd2( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, idstnd, ltimecnt )

call alldgatherv(gdcr,nplwex*nbnod,cgjr,ncgcnt,ncgdsp)

if( ltimecnt ) then
    ct = timecnt()
    tc(4) = tc(4) + ct - ct0
    if( la2aout ) then
        if(loutfile(1)) write(nfile(1),*) '         unify w.f. (by alldgatherv)',  &
&                         ' : cpu-time :', tc(4)
        if(loutfile(2)) write(nfile(2),*) '         unify w.f. (by alldgatherv)',  &
&                         ' : cpu-time :', tc(4)
        tc(4) = 0.d0
    endif
    ct0 = ct
endif
#endif


return
end




subroutine pwpunifylc( nfile, myid, nodes,  &
& hdiag, vexc, vhar, vhshdp, vext, mshnod, glocal, vlocud,  &
& mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d )
!-----------------------------------------------------------------------
!   unify local potential/density
!-----------------------------------------------------------------------
use commwv_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: hdiag(*), vexc(*), vhar(*), vhshdp(*), vext(*)
real*8  :: glocal(*), vlocud(*)
integer :: mshnod, mftnod(*), mftdsp(*), mfd2ft(*), ntotfd, nd1vks(*)
integer :: kfft0d


!if( .not.lnoncollinear ) then
    call pwpunifylc2( nfile, myid, nodes,  &
& hdiag, vexc, vhar, vhshdp, vext, mshnod, glocal,  &
& mftnod, mftdsp, mfd2ft, ntotfd, nd1vks )
!else
!    !-----noncollinear magnetism
!    call ncpwpunifylc2( nfile, myid, nodes,  &
!& hdiag, vhar, vhshdp, vext, mshnod, glocal, vlocud,  &
!& mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d )
!end if


return
end





subroutine pwpunifylc2( nfile, myid, nodes,  &
& hdiag, vexc, vhar, vhshdp, vext, mshnod, glocal,  &
& mftnod, mftdsp, mfd2ft, ntotfd, nd1vks )
!-----------------------------------------------------------------------
!   unify local potential/density
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
real*8  :: hdiag(*), vexc(*), vhar(*), vhshdp(*), vext(*), glocal(*)
integer :: mshnod, mftnod(*), mftdsp(*), mfd2ft(*), ntotfd, nd1vks(*)


hdiag(1:mshnod) = vhar(1:mshnod) + vhshdp(1:mshnod) + vext(1:mshnod)  &
&               + vexc(1:mshnod)

call unifylc( nfile, myid, nodes,  &
& hdiag, mshnod, glocal, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks )


return
end




subroutine unifylc( nfile, myid, nodes,  &
& hdiag, mshnod, glocal, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks )
!-----------------------------------------------------------------------
!   unify local potential/density
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
real*8  :: hdiag(*), glocal(*)
integer :: mshnod, mftnod(*), mftdsp(*), mfd2ft(*), ntotfd, nd1vks(*)


call alldgatherv(hdiag,mshnod,glocal,mftnod,mftdsp)


return
end




subroutine distlc( nfile, myid, nodes,  &
& hdiag, mshnod, glocal, mftdsp )
!-----------------------------------------------------------------------
!   store local values in local variables
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
real*8  :: hdiag(*), glocal(*)
integer :: mshnod, mftdsp(*)

!------declare local variables
integer :: mdsp, m


mdsp = mftdsp(myid+1)
do m = 1, mshnod
   hdiag(m) = glocal(m+mdsp)
enddo


return
end




subroutine unifyeg( nfile, myid, nodes, ct0,  &
& eig, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, iod, iodg,  &
& prod, prodr, ltimecnt, leig_start, leig_end )
!-----------------------------------------------------------------------
!   unify eigenvalues
!-----------------------------------------------------------------------
use outfile
implicit real*8 ( a-h, o-z )
dimension nfile(*)
dimension eig(*)
dimension iod(*), iodg(*), nbndsp(*), nbncnt(*)
dimension prod(*), prodr(*)
logical ltimecnt, leig_start, leig_end
real*8  :: tc = 0.d0
save tc


if( leig_start ) then
if( ltimecnt ) then
    tc = 0.d0
end if
end if

do i = 1, nbnod
   ib = iod(i)
   prod(i) = eig(ib)
enddo
call alldgatherv(prod,nbnod,prodr,nbncnt,nbndsp)
do i = 1, nband
   ib = iodg(i)
   eig(ib) = prodr(i)
enddo

if( ltimecnt ) then
    ct = timecnt()
    tc = tc + ct - ct0
    ct0 = ct
end if

if( leig_end ) then
if( ltimecnt ) then
    if(loutfile(1)) write(nfile(1),*) '                          unify eig ',  &
&                         ' : cpu-time :', tc
    if(loutfile(2)) write(nfile(2),*) '                          unify eig ',  &
&                         ' : cpu-time :', tc
end if
end if


return
end




subroutine bdtogd( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iodg, idstnd, lflag, ltimecnt )
!-----------------------------------------------------------------------
!    all to all communication 
!    to convert band decomposition to G decomposition
!
! (input)
!    cgjr : wavefunctions decomposed by bands
!
! (output)
!    gdcr : wavefunctions decomposed by components
!-----------------------------------------------------------------------
use commwv_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: ct0
integer :: nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt(*), nbndsp(*)
real*8  :: cgjr(*)
integer :: npnod1, npnod2, npnod, nplcnt(*), npldsp(*)
real*8  :: gdcr(2*npnod,*)
real*8  :: bufcr(*)
integer :: iodg(*), idstnd(*)
logical :: lflag, ltimecnt


if( .not.lnoncollinear ) then
    call bdtogd2( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iodg, idstnd, lflag, ltimecnt )
else
    !-----noncollinear magnetism
    call bdtogd2( nfile, myid, nodes, ct0,  &
& cgjr, nplwex*2, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iodg, idstnd, lflag, ltimecnt )
end if


return
end




subroutine bdtogd2( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iodg, idstnd, lflag, ltimecnt )
!-----------------------------------------------------------------------
!    all to all communication 
!    to convert band decomposition to G decomposition
!
! (input)
!    cgjr : wavefunctions decomposed by bands
!
! (output)
!    gdcr : wavefunctions decomposed by components
!-----------------------------------------------------------------------
use outfile
use all2all_out
implicit real*8 ( a-h, o-z )
dimension nfile(*)
dimension cgjr(nplwex,*)
dimension gdcr(2*npnod,*)
dimension bufcr(*)
dimension iodg(*), nbndsp(*), nbncnt(*)
dimension nplcnt(*), npldsp(*)
dimension idstnd(*)
logical   lflag, ltimecnt


!--- within myself
ir1 = 2*( npldsp(myid+1) + 1 ) - 1
ir2 = 2*( npldsp(myid+1) + nplcnt(myid+1) )
ir3 = ir2 - ir1 + 1
ib0 = nbndsp(myid+1)
do i = 1, nbnod
   do ig = ir1, ir2
      gdcr(ig-ir1+1,i+ib0) = cgjr(ig,i)
   enddo
enddo

do isnd = 1, nodes - 1
!--- destination node
   idst = idstnd(isnd)

!--- set No. of sending data (components)
   is1 = 2*( npldsp(idst+1) + 1 ) - 1
   is2 = 2*( npldsp(idst+1) + nplcnt(idst+1) )
   is3 = is2 - is1 + 1

!--- set No. of recieving data (bands)
   ibs3 = nbncnt(idst+1)
   ib0  = nbndsp(idst+1)

   call bdtog2( nfile, myid, nodes,  &
& is1, is2, is3, ibs3, idst,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod,  &
& bufcr, gdcr(1,ib0+1), npnod1, npnod2, npnod, isnd )

!---------Internode synchronization
!         call gsync
enddo

if( lflag ) then
!--- store data in correct order
    call bdtog3( nfile, myid, nodes,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, cgjr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iodg, idstnd )
endif


if( ltimecnt ) then
    ct = timecnt()
    tc(1) = tc(1) + ct - ct0
    if( la2aout ) then
        if(loutfile(1)) write(nfile(1),*) '               all to all (bd-to-gd)',  &
&                         ' : cpu-time :', tc(1)
        if(loutfile(2)) write(nfile(2),*) '               all to all (bd-to-gd)',  &
&                         ' : cpu-time :', tc(1)
        tc(1) = 0.d0
    end if
    ct0 = ct
endif


return
end




subroutine bdtog2( nfile, myid, nodes,  &
& is1, is2, is3, ibs3, idst,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod,  &
& bufcr, gdcr, npnod1, npnod2, npnod, isnd )
!-----------------------------------------------------------------------
!    all to all communication 
!    to convert band decomposition to G decomposition
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension nfile(*)
dimension cgjr(nplwex,*)
dimension gdcr(2*npnod,*)
dimension bufcr(is3,*)


!--- set sending data
do i = 1, nbnod
   do ig = is1, is2
      bufcr(ig-is1+1,i) = cgjr(ig,i)
   enddo
enddo

nsd = is3*nbnod
nrc = 2*npnod*ibs3
call mespasd(120+isnd,bufcr,gdcr,nsd,nrc,idst,0)  
!c---------myid.lt.idst: send & recv, if not empty
!          if (myid.lt.idst) then
!            if (nsd.ne.0) call cdsend(120,bufcr,nsd,idst,0)
!            if (nrc.ne.0) call cdrecv(130,gdcr,nrc,0)
!c---------myid.gt.idst: recv & send, if not empty
!          else
!            if (nrc.ne.0) call cdrecv(120,gdcr,nrc,0)
!            if (nsd.ne.0) call cdsend(130,bufcr,nsd,idst,0)
!          endif

return
end




subroutine bdtog3( nfile, myid, nodes,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, cgjr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iodg, idstnd )
!-----------------------------------------------------------------------
!    all to all communication 
!    to convert band decomposition to G decomposition
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension nfile(*)
dimension cgjr(2*npnod,*)
dimension gdcr(2*npnod,*)
dimension iodg(*), nbndsp(*), nbncnt(*)
dimension nplcnt(*), npldsp(*)
dimension idstnd(*)


do i = 1, nband
   ib = iodg(i)
   do ig = 1, 2*npnod
      cgjr(ig,ib) = gdcr(ig,i)
   enddo
enddo
do ib = 1, nband
   do ig = 1, 2*npnod
      gdcr(ig,ib) = cgjr(ig,ib)
   enddo
enddo


return
end




subroutine gdtobd( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, idstnd, ltimecnt )
!-----------------------------------------------------------------------
!    all to all communication 
!    to convert G decomposition to band decomposition
!
! (input)
!    gdcr : wavefunctions decomposed by components
!
! (output)
!    gdcr : wavefunctions decomposed by bands
!-----------------------------------------------------------------------
use commwv_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: ct0
integer :: nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt(*), nbndsp(*)
real*8  :: cgjr(*)
integer :: npnod1, npnod2, npnod, nplcnt(*), npldsp(*)
real*8  :: gdcr(*)
real*8  :: bufcr(*)
integer :: iod(*), iodg(*), idstnd(*)
logical :: ltimecnt


if( .not.lnoncollinear ) then
    call gdtobd2( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, idstnd, ltimecnt )
else
    !-----noncollinear magnetism
    call gdtobd2( nfile, myid, nodes, ct0,  &
& cgjr, nplwex*2, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, idstnd, ltimecnt )
end if


return
end




subroutine gdtobd2( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, idstnd, ltimecnt )
!-----------------------------------------------------------------------
!    all to all communication 
!    to convert G decomposition to band decomposition
!
! (input)
!    gdcr : wavefunctions decomposed by components
!
! (output)
!    gdcr : wavefunctions decomposed by bands
!-----------------------------------------------------------------------
use outfile
use all2all_out
implicit real*8 ( a-h, o-z )
dimension nfile(*)
dimension cgjr(nbnod,*)
dimension gdcr(*)
dimension bufcr(*)
dimension iod(*), iodg(*), nbndsp(*), nbncnt(*)
dimension nplcnt(*), npldsp(*)
dimension idstnd(*)
logical   ltimecnt


!--- within myself
ir1 = 2*( npldsp(myid+1) + 1 ) - 1
ir2 = 2*( npldsp(myid+1) + nplcnt(myid+1) )
ir3 = ir2 - ir1 + 1
call gdtob1( nfile, myid, nodes,  &
& ir1, ir2, ir3, cgjr, nband, nbnod1, nbnod2, nbnod,  &
& gdcr, npnod1, npnod2, npnod, iod )

do isnd = 1, nodes - 1
!--- destination node
   idst = idstnd(isnd)

!--- set No. of sending data (bands)
   is1 = nbndsp(idst+1) + 1
   is2 = nbndsp(idst+1) + nbncnt(idst+1)
   is3 = is2 - is1 + 1

!--- set No. of recieving data (components)
   ir1 = 2*( npldsp(idst+1) + 1 ) - 1
   ir2 = 2*( npldsp(idst+1) + nplcnt(idst+1) )
   ir3 = ir2 - ir1 + 1

   if( ir3 == 0 ) ir1 = 1
   call gdtob2( nfile, myid, nodes,  &
& is1, is2, is3, ir3, idst,  &
& cgjr(1,ir1), nband, nbnod1, nbnod2, nbnod,  &
& bufcr, gdcr, npnod1, npnod2, npnod, iodg, isnd )

!---------Internode synchronization
!         call gsync
enddo

!--- store data in correct order
call gdtob3( nfile, myid, nodes, gdcr, cgjr, nplwex, nbnod )


if( ltimecnt ) then
    ct = timecnt()
    tc(2) = tc(2) + ct - ct0
    if( la2aout ) then
        if(loutfile(1)) write(nfile(1),*) '               all to all (gd-to-bd)',  &
&                         ' : cpu-time :', tc(2)
        if(loutfile(2)) write(nfile(2),*) '               all to all (gd-to-bd)',  &
&                         ' : cpu-time :', tc(2)
        tc(2) = 0.d0
    end if
    ct0 = ct
endif


return
end




subroutine gdtob1( nfile, myid, nodes,  &
& ir1, ir2, ir3, cgjr, nband, nbnod1, nbnod2, nbnod,  &
& gdcr, npnod1, npnod2, npnod, iod )
!-----------------------------------------------------------------------
!    all to all communication 
!    to convert G decomposition to band decomposition
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension nfile(*)
dimension cgjr(nbnod,*)
dimension gdcr(2*npnod,*)
dimension iod(*)


!--- within myself
do i = 1, nbnod
   ib = iod(i)
   do ig = ir1, ir2
      cgjr(i,ig) = gdcr(ig-ir1+1,ib)
   enddo
enddo


return
end




subroutine gdtob2( nfile, myid, nodes,  &
& is1, is2, is3, ibs3, idst,  &
& cgjr, nband, nbnod1, nbnod2, nbnod,  &
& bufcr, gdcr, npnod1, npnod2, npnod, iodg, isnd )
!-----------------------------------------------------------------------
!    all to all communication 
!    to convert G decomposition to band decomposition
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension nfile(*)
dimension cgjr(nbnod,*)
dimension gdcr(2*npnod,*)
dimension bufcr(is3,*)
dimension iodg(*)


!--- set sending data
do i = is1, is2
   ib = iodg(i)
   do ig = 1, 2*npnod
      bufcr(i-is1+1,ig) = gdcr(ig,ib)
   enddo
enddo

nsd = is3*2*npnod
nrc = nbnod*ibs3
call mespasd(120+isnd,bufcr,cgjr,nsd,nrc,idst,0)  
!c---------myid.lt.idst: send & recv, if not empty
!          if (myid.lt.idst) then
!            if (nsd.ne.0) call cdsend(120,bufcr,nsd,idst,0)
!            if (nrc.ne.0) call cdrecv(130,cgjr,nrc,0)
!c---------myid.gt.idst: recv & send, if not empty
!          else
!            if (nrc.ne.0) call cdrecv(120,cgjr,nrc,0)
!            if (nsd.ne.0) call cdsend(130,bufcr,nsd,idst,0)
!          endif

return
end




subroutine gdtob3( nfile, myid, nodes, gdcr, cgjr, nplwex, nbnod )
!-----------------------------------------------------------------------
!    all to all communication 
!    to convert G decomposition to band decomposition
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension nfile(*)
dimension cgjr(nbnod,*)
dimension gdcr(nplwex,*)


do i = 1, nbnod
   do ig = 1, nplwex
      gdcr(ig,i) = cgjr(i,ig)
   enddo
enddo


return
end




#if CCREAL4
subroutine rgdtobd( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, idstnd, nspnod, npnnbnx, ltimecnt )
!-----------------------------------------------------------------------
!    all to all communication 
!    to convert G decomposition to band decomposition
!
!   << for real*4 variables >>
!
! (input)
!    gdcr : wavefunctions decomposed by components
!
! (output)
!    gdcr : wavefunctions decomposed by bands
!-----------------------------------------------------------------------
use outfile
use all2all_out
implicit real*8 ( a-h, o-z )
dimension nfile(*)
real*4    cgjr(nbnod,nplwex)
real*4    gdcr(nspnod)
real*4    bufcr(npnnbnx)
dimension iod(nbnod), iodg(nband), nbndsp(nodes), nbncnt(nodes)
dimension nplcnt(nodes), npldsp(nodes)
dimension idstnd(nodes)
logical   ltimecnt


!--- within myself
ir1 = 2*( npldsp(myid+1) + 1 ) - 1
ir2 = 2*( npldsp(myid+1) + nplcnt(myid+1) )
ir3 = ir2 - ir1 + 1
call rgdtob1( nfile, myid, nodes,  &
& ir1, ir2, ir3, cgjr, nplwex, nband, nbnod1, nbnod2, nbnod,  &
& gdcr, npnod1, npnod2, npnod, iod )

do isnd = 1, nodes - 1
!--- destination node
   idst = idstnd(isnd)

!--- set No. of sending data (bands)
   is1 = nbndsp(idst+1) + 1
   is2 = nbndsp(idst+1) + nbncnt(idst+1)
   is3 = is2 - is1 + 1

!--- set No. of recieving data (components)
   ir1 = 2*( npldsp(idst+1) + 1 ) - 1
   ir2 = 2*( npldsp(idst+1) + nplcnt(idst+1) )
   ir3 = ir2 - ir1 + 1

   if( ir3 == 0 ) ir1 = 1
   call rgdtob2( nfile, myid, nodes,  &
& is1, is2, is3, ir3, idst,  &
& cgjr(1,ir1), nplwex-ir1+1, nband, nbnod1, nbnod2, nbnod,  &
& bufcr, gdcr, npnod1, npnod2, npnod, iodg, isnd )

!---------Internode synchronization
!         call gsync
enddo

!--- store data in correct order
call rgdtob3( nfile, myid, nodes, gdcr, cgjr, nplwex, nbnod )


if( ltimecnt ) then
    ct = timecnt()
    tc(3) = tc(3) + ct - ct0
    if( la2aout ) then
        if(loutfile(1)) write(nfile(1),*) '        real*4 all to all (gd-to-bd)',  &
&                         ' : cpu-time :', tc(3)
        if(loutfile(2)) write(nfile(2),*) '        real*4 all to all (gd-to-bd)',  &
&                         ' : cpu-time :', tc(3)
        tc(3) = 0.d0
    end if
    ct0 = ct
endif


return
end




subroutine rgdtob1( nfile, myid, nodes,  &
& ir1, ir2, ir3, cgjr, nplwex, nband, nbnod1, nbnod2, nbnod,  &
& gdcr, npnod1, npnod2, npnod, iod )
!-----------------------------------------------------------------------
!    all to all communication 
!    to convert G decomposition to band decomposition
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension nfile(*)
real*4    cgjr(nbnod,nplwex)
real*4    gdcr(2*npnod,nband)
dimension iod(nbnod)


!--- within myself
do i = 1, nbnod
   ib = iod(i)
   do ig = ir1, ir2
      cgjr(i,ig) = gdcr(ig-ir1+1,ib)
   enddo
enddo


return
end




subroutine rgdtob2( nfile, myid, nodes,  &
& is1, is2, is3, ibs3, idst,  &
& cgjr, nplwex, nband, nbnod1, nbnod2, nbnod,  &
& bufcr, gdcr, npnod1, npnod2, npnod, iodg, isnd )
!-----------------------------------------------------------------------
!    all to all communication 
!    to convert G decomposition to band decomposition
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension nfile(*)
real*4    cgjr(nbnod,nplwex)
real*4    gdcr(2*npnod,nband)
real*4    bufcr(is3,2*npnod)
dimension iodg(nband)


!--- set sending data
do i = is1, is2
   ib = iodg(i)
   do ig = 1, 2*npnod
      bufcr(i-is1+1,ig) = gdcr(ig,ib)
   enddo
enddo

nsd = is3*2*npnod
nrc = nbnod*ibs3
call mespasr(120+isnd,bufcr,cgjr,nsd,nrc,idst,0)  
!c---------myid.lt.idst: send & recv, if not empty
!          if (myid.lt.idst) then
!            if (nsd.ne.0) call crsend(120,bufcr,nsd,idst,0)
!            if (nrc.ne.0) call crrecv(130,cgjr,nrc,0)
!c---------myid.gt.idst: recv & send, if not empty
!          else
!            if (nrc.ne.0) call crrecv(120,cgjr,nrc,0)
!            if (nsd.ne.0) call crsend(130,bufcr,nsd,idst,0)
!          endif

return
end




subroutine rgdtob3( nfile, myid, nodes, gdcr, cgjr, nplwex,nbnod )
!-----------------------------------------------------------------------
!    all to all communication 
!    to convert G decomposition to band decomposition
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension nfile(*)
real*4    cgjr(nbnod,nplwex)
real*4    gdcr(nplwex,nbnod)


do i = 1, nbnod
   do ig = 1, nplwex
      gdcr(ig,i) = cgjr(i,ig)
   enddo
enddo


return
end
#endif




subroutine slm_bdtog3( nfile, myid, nodes,  &
& vbyplw, bufcr, multiple, nband, npnod, npnodx, iodg )
!-----------------------------------------------------------------------
!    all to all communication 
!    to convert band decomposition to G decomposition
!-----------------------------------------------------------------------
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: multiple, nband
integer :: npnod, npnodx
real*8  :: vbyplw(multiple*npnodx,nband)
real*8  :: bufcr(multiple*npnodx,nband)
integer :: iodg(nband)

!------declare local variables
integer :: i, ig, ib


do i = 1, nband
   ib = iodg(i)
   bufcr(1:multiple*npnod,ib) = vbyplw(1:multiple*npnod,i)
enddo

vbyplw(1:multiple*npnod,1:nband) = bufcr(1:multiple*npnod,1:nband)


return
end




subroutine slm_bdtovgd_org( nfile, myid, nodes, ct0,  &
& vbybnd, vbyplw, bufcr,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& multiple, nofplw, npnod1, npnod2, npnod, nplcnt, npldsp, npnodx,  &
& iodg, ioag, idstnd, lflag, ltimecnt )
!-----------------------------------------------------------------------
!    all to all communication 
!    to convert band decomposition to atom decomposition
!
! (input)
!    vbybnd : variables decomposed by bands
!
! (output)
!    vbyplw : variables decomposed by planewaves or atoms
!-----------------------------------------------------------------------
use outfile
use all2all_out
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
real*8  :: ct0
integer :: multiple, nband, nbnod1, nbnod2, nbnod, nbncnt(nodes), nbndsp(nodes)
integer :: nofplw, npnod1, npnod2, npnod, nplcnt(nodes), npldsp(nodes), npnodx
real*8  :: vbybnd(multiple,nofplw,nbnod)
real*8  :: vbyplw(multiple,npnodx,nband)
real*8  :: bufcr(*)
integer :: iodg(nband), ioag(nofplw), idstnd(nodes)
logical :: lflag, ltimecnt

!------declare local variables
integer :: ir1, ir2, ir3, ib0, i, ig, isnd, idst, is1, is2, is3, ibs3, igg
real*8  :: timecnt, ct


!--- within myself
ir1 = npldsp(myid+1) + 1
ir2 = npldsp(myid+1) + nplcnt(myid+1)
ir3 = ir2 - ir1 + 1
ib0 = nbndsp(myid+1)
do i = 1, nbnod
   do ig = ir1, ir2
      igg = ioag(ig)
      vbyplw(1:multiple,ig-ir1+1,i+ib0) = vbybnd(1:multiple,igg,i)
   enddo
enddo

do isnd = 1, nodes - 1
!--- destination node
   idst = idstnd(isnd)

!--- set No. of sending data (components)
   is1 = npldsp(idst+1) + 1
   is2 = npldsp(idst+1) + nplcnt(idst+1)
   is3 = is2 - is1 + 1

!--- set No. of recieving data (bands)
   ibs3 = nbncnt(idst+1)
   ib0  = nbndsp(idst+1)

   call slm_bdtovg2( nfile, myid, nodes,  &
& is1, is2, is3, ibs3, idst,  &
& vbybnd, multiple, nband, nbnod1, nbnod2, nbnod,  &
& vbyplw(1,1,ib0+1), nofplw, npnod1, npnod2, npnod, npnodx,  &
& isnd, ioag, bufcr )

!---------Internode synchronization
!         call gsync
enddo

if( lflag ) then
!--- store data in correct order
    call slm_bdtog3( nfile, myid, nodes,  &
& vbyplw, bufcr, multiple, nband, npnod, npnodx, iodg )
endif


if( ltimecnt ) then
    ct = timecnt()
    tc(5) = tc(5) + ct - ct0
    if( la2aout ) then
        if(loutfile(1)) write(nfile(1),*) '           all to all (slm-bd-to-gd)',  &
&                         ' : cpu-time :', tc(5)
        if(loutfile(2)) write(nfile(2),*) '           all to all (slm-bd-to-gd)',  &
&                         ' : cpu-time :', tc(5)
        tc(5) = 0.d0
    end if
    ct0 = ct
endif


return
end




subroutine slm_bdtovg2( nfile, myid, nodes,  &
& is1, is2, is3, ibs3, idst,  &
& vbybnd, multiple, nband, nbnod1, nbnod2, nbnod,  &
& vbyplw, nofplw, npnod1, npnod2, npnod, npnodx, isnd, ioag, bufcr )
!-----------------------------------------------------------------------
!    all to all communication 
!    to convert band decomposition to G decomposition
!-----------------------------------------------------------------------
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: is1, is2, is3, ibs3, idst
integer :: multiple, nband, nbnod1, nbnod2, nbnod
integer :: nofplw, npnod1, npnod2, npnod, npnodx
real*8  :: vbybnd(multiple,nofplw,nbnod)
real*8  :: vbyplw(multiple,npnodx,ibs3)
integer :: isnd, ioag(nofplw)
real*8  :: bufcr(multiple,is3,nbnod)

!------declare local variables
integer :: i, ig, nsd, nrc, igg


!--- set sending data
do i = 1, nbnod
   do ig = is1, is2
      igg = ioag(ig)
      bufcr(1:multiple,ig-is1+1,i) = vbybnd(1:multiple,igg,i)
   enddo
enddo

nsd = multiple*is3*nbnod
nrc = multiple*npnod*ibs3
call mespasd(120+isnd,bufcr,vbyplw,nsd,nrc,idst,0)  
!c---------myid.lt.idst: send & recv, if not empty
!          if (myid.lt.idst) then
!            if (nsd.ne.0) call cdsend(120,bufcr,nsd,idst,0)
!            if (nrc.ne.0) call cdrecv(130,gdcr,nrc,0)
!c---------myid.gt.idst: recv & send, if not empty
!          else
!            if (nrc.ne.0) call cdrecv(120,gdcr,nrc,0)
!            if (nsd.ne.0) call cdsend(130,bufcr,nsd,idst,0)
!          endif

return
end




subroutine slm_gdtob1( nfile, myid, nodes,  &
& ir1, ir2, ir3, worka, nband, nbnod1, nbnod2, nbnod,  &
& vbyplw, multiple, npnod1, npnod2, npnod, iod )
!-----------------------------------------------------------------------
!    all to all communication 
!    to convert G decomposition to band decomposition
!-----------------------------------------------------------------------
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: ir1, ir2, ir3
integer :: nband, nbnod1, nbnod2, nbnod
integer :: multiple, npnod1, npnod2, npnod
real*8  :: worka(nbnod,*)
real*8  :: vbyplw(multiple*npnod,*)
integer :: iod(nbnod)

!------declare local variables
integer :: i, ib, ig


!--- within myself
do i = 1, nbnod
   ib = iod(i)
   do ig = ir1, ir2
      worka(i,ig) = vbyplw(ig-ir1+1,ib)
   enddo
enddo


return
end




subroutine slm_gdtob2( nfile, myid, nodes,  &
& is1, is2, is3, ibs3, idst,  &
& worka, nband, nbnod1, nbnod2, nbnod,  &
& bufcr, vbyplw, multiple, npnod1, npnod2, npnod, iodg, isnd )
!-----------------------------------------------------------------------
!    all to all communication 
!    to convert G decomposition to band decomposition
!-----------------------------------------------------------------------
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: is1, is2, is3, ibs3, idst
integer :: nband, nbnod1, nbnod2, nbnod
integer :: multiple, npnod1, npnod2, npnod
real*8  :: worka(nbnod,*)
real*8  :: vbyplw(multiple*npnod,*)
real*8  :: bufcr(is3,*)
integer :: iodg(nband), isnd

!------declare local variables
integer :: i, ib, ig, nsd, nrc


!--- set sending data
do i = is1, is2
   ib = iodg(i)
   do ig = 1, multiple*npnod
      bufcr(i-is1+1,ig) = vbyplw(ig,ib)
   enddo
enddo

nsd = is3*multiple*npnod
nrc = nbnod*ibs3
call mespasd(120+isnd,bufcr,worka,nsd,nrc,idst,0)  
!c---------myid.lt.idst: send & recv, if not empty
!          if (myid.lt.idst) then
!            if (nsd.ne.0) call cdsend(120,bufcr,nsd,idst,0)
!            if (nrc.ne.0) call cdrecv(130,worka,nrc,0)
!c---------myid.gt.idst: recv & send, if not empty
!          else
!            if (nrc.ne.0) call cdrecv(120,worka,nrc,0)
!            if (nsd.ne.0) call cdsend(130,bufcr,nsd,idst,0)
!          endif

return
end




subroutine slm_vgdtobd_org( nfile, myid, nodes, ct0,  &
& vbyplw, vbybnd, worka, bufcr,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& multiple, nofplw, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, ioag, idstnd, ltimecnt )
!-----------------------------------------------------------------------
!    all to all communication 
!    to convert G decomposition to band decomposition
!
! (input)
!    vbyplw : variables decomposed by planewaves or atoms
!
! (output)
!    vbybnd : variables decomposed by bands
!
! Note:
!    worka : work array with the same saze as vbybnd
!-----------------------------------------------------------------------
use outfile
use all2all_out
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
real*8  :: ct0
integer :: nband, nbnod1, nbnod2, nbnod, nbncnt(nodes), nbndsp(nodes)
integer :: multiple, nofplw, npnod1, npnod2, npnod, nplcnt(nodes), npldsp(nodes)
real*8  :: vbyplw(multiple*npnod,nband)
real*8  :: vbybnd(multiple*nofplw,nbnod)
real*8  :: worka(nbnod,multiple*nofplw)
real*8  :: bufcr(*)
integer :: iod(nbnod), iodg(nband), ioag(nofplw), idstnd(nodes)
logical :: ltimecnt

!------declare local variables
integer :: ir1, ir2, ir3, isnd, idst, is1, is2, is3
real*8  :: timecnt, ct


!--- within myself
ir1 = multiple*( npldsp(myid+1) ) + 1
ir2 = multiple*( npldsp(myid+1) + nplcnt(myid+1) )
ir3 = ir2 - ir1 + 1
call slm_gdtob1( nfile, myid, nodes,  &
& ir1, ir2, ir3, worka, nband, nbnod1, nbnod2, nbnod,  &
& vbyplw, multiple, npnod1, npnod2, npnod, iod )

do isnd = 1, nodes - 1
!--- destination node
   idst = idstnd(isnd)

!--- set No. of sending data (bands)
   is1 = nbndsp(idst+1) + 1
   is2 = nbndsp(idst+1) + nbncnt(idst+1)
   is3 = is2 - is1 + 1

!--- set No. of recieving data (components)
   ir1 = multiple*( npldsp(idst+1) ) + 1
   ir2 = multiple*( npldsp(idst+1) + nplcnt(idst+1) )
   ir3 = ir2 - ir1 + 1

   if( ir3 == 0 ) ir1 = 1
   call slm_gdtob2( nfile, myid, nodes,  &
& is1, is2, is3, ir3, idst,  &
& worka(1,ir1), nband, nbnod1, nbnod2, nbnod,  &
& bufcr, vbyplw, multiple, npnod1, npnod2, npnod, iodg, isnd )

!---------Internode synchronization
!         call gsync
enddo

!--- store data in correct order
call slm_vgdtob3( nfile, myid, nodes,  &
& vbybnd, worka, multiple, nofplw, nbnod, ioag )


if( ltimecnt ) then
    ct = timecnt()
    tc(6) = tc(6) + ct - ct0
    if( la2aout ) then
        if(loutfile(1)) write(nfile(1),*) '           all to all (slm-gd-to-bd)',  &
&                         ' : cpu-time :', tc(6)
        if(loutfile(2)) write(nfile(2),*) '           all to all (slm-gd-to-bd)',  &
&                         ' : cpu-time :', tc(6)
        tc(6) = 0.d0
    end if
    ct0 = ct
endif


return
end




subroutine slm_vgdtob3( nfile, myid, nodes,  &
& vbybnd, worka, multiple, nofplw, nbnod, ioag )
!-----------------------------------------------------------------------
!    all to all communication 
!    to convert G decomposition to band decomposition
!-----------------------------------------------------------------------
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: multiple, nofplw, nbnod
integer :: ioag(nofplw)
real*8  :: worka(nbnod,multiple,nofplw)
real*8  :: vbybnd(multiple,nofplw,nbnod)

!------declare local variables
integer :: i, ig, ia


do ig = 1, nofplw
   ia = ioag(ig)
   do i = 1, nbnod
      vbybnd(1:multiple,ia,i) = worka(i,1:multiple,ig)
   enddo
enddo


return
end




subroutine dstata( nfile, myid, nodes, idstnd )
!-----------------------------------------------------------------------
!    set destination nodes for all-to-all communication
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension nfile(*)
dimension idstnd(*)


imyid = myid + 1

do i = 1, min(imyid-1,nodes+1-imyid)
   idstnd(imyid-2+i) = i
enddo
do i = imyid+1,nodes+1-imyid
   idstnd(imyid-2+i) = i
enddo
if( imyid.ne.nodes ) then
    do i = nodes+2-imyid, nodes-1
       if( i.ne.imyid ) then
           idstnd(imyid+i-nodes-1) = i
       endif
    enddo
    if( imyid.ge.2 .and. imyid.le.nodes/2 ) then
        idstnd(2*(imyid-1)) = nodes
    else if( imyid.ge.nodes/2+1 ) then
        idstnd(2*(imyid-nodes/2)-1) = nodes
    endif
  else
    do i = 2, nodes/2
       idstnd(2*(i-1)) = i
    enddo
    do i = nodes/2+1, nodes-1
       idstnd(2*(i-nodes/2)-1) = i
    enddo
endif
do i = 1, nodes-1
   idstnd(i) = idstnd(i) - 1
enddo


return
end




subroutine cprhcg( nfile, myid, nodes, cgjr, rhcr, nplwex, nbnod )
!-----------------------------------------------------------------------
!    copy rhcr to cgjr
!-----------------------------------------------------------------------
use commwv_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: nplwex, nbnod
real*8  :: cgjr(*), rhcr(*)

!------declare local variables
integer :: ndat

if( .not.lnoncollinear ) then
    ndat = nbnod*nplwex
else
    !-----noncollinear magnetism
    ndat = nbnod*nplwex*2
end if
cgjr(1:ndat) = rhcr(1:ndat)


return
end




subroutine cpgdrh( nfile, myid, nodes, rhcr, gdcr, npnod, nband )
!-----------------------------------------------------------------------
!    copy gdcr to rhcr
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension nfile(*)
dimension gdcr(*)
dimension rhcr(*)

!------declare local variables
integer :: ndat

ndat = 2*npnod*nband
rhcr(1:ndat) = gdcr(1:ndat)

return
end




subroutine cpgdg4( nfile, myid, nodes, rhcr, gdcr, npnod, nband )
!-----------------------------------------------------------------------
!    copy gdcr (real*8) to rhcr (real*4)
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension nfile(*)
dimension gdcr(*)
real*4    rhcr(*)

!------declare local variables
integer :: ndat

ndat = 2*npnod*nband
rhcr(1:ndat) = gdcr(1:ndat)

return
end




subroutine cpgdg47( nfile, myid, nodes,  &
& rhcr, gdcr, npnod, npnod7, nband )
!-----------------------------------------------------------------------
!    copy gdcr (real*8) to rhcr (real*4)
!-----------------------------------------------------------------------
use commwv_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: npnod, npnod7, nband
real*8  :: gdcr(2*npnod,*)
real*4  :: rhcr(2*npnod7,*)

!------declare local variables
integer :: ndat


if( .not.lnoncollinear .or. nodes > 1 ) then
    ndat = 2*npnod7
    rhcr(1:ndat,1:nband) = gdcr(1:ndat,1:nband)
else
    !-----noncollinear magnetism & nodes == 1
    ndat = 2*npnod7/2
    rhcr(1:ndat,1:nband) = gdcr(1:ndat,1:nband)
    rhcr(ndat+1:ndat*2,1:nband) = gdcr(npnod+1:npnod+ndat,1:nband)
end if

return
end




subroutine cpcgrd( nfile, myid, nodes, cgjr, gdcr, nplwex, nbnod, iod )
!-----------------------------------------------------------------------
!    copy rhcr to cgjr
!-----------------------------------------------------------------------
use commwv_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: nplwex, nbnod, iod(*)
real*8  :: cgjr(*), gdcr(*)

if( .not.lnoncollinear ) then
    call cpcgrd2( nfile, myid, nodes, cgjr, gdcr, nplwex, nbnod, iod )
else
    !-----noncollinear magnetism
    call cpcgrd2( nfile, myid, nodes, cgjr, gdcr, nplwex*2, nbnod, iod )
end if

return
end


subroutine cpcgrd2( nfile, myid, nodes, cgjr, gdcr, nplwex, nbnod, iod )
use commwv_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: nplwex, nbnod, iod(*)
real*8  :: cgjr(nplwex,*), gdcr(nplwex,*)

!------declare local variables
integer :: i, ib

do i = 1, nbnod
   ib = iod(i)
   gdcr(1:nplwex,i) = cgjr(1:nplwex,ib)
end do

return
end




subroutine lcalldgatherv( nfile, myid, nodes, idroot,  &
& buff, n, buffr, nbncnt, nbndsp, idstnd )
!-----------------------------------------------------------------------
!     alldgather within local nodes:
!     node id's are given by idstnd
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension nfile(*)
dimension buff(*), buffr(*)
dimension nbndsp(*), nbncnt(*)
dimension idstnd(*)

!--- within myself
ib0 = nbndsp(myid+1)
do i = 1, n
   buffr(i+ib0) = buff(i)
enddo

do isnd = 1, nodes - 1
!--- destination node
   idst = idstnd(isnd)

!--- set No. of recieving data
   ib0 = nbndsp(idst+1)

   nsd = n
   nrc = nbncnt(idst+1)
   if( nrc > 0 ) then
       ib01 = ib0 + 1
     else
       ib01 = 1
   end if
   call mespasd(400+isnd,buff,buffr(ib01),nsd,nrc,idst+idroot,0)

!---------Internode synchronization
!         call gsync
enddo

return
end




subroutine lcdsum( nfile, myid, nodes, idroot,  &
& buffs, n, temp, buffr, idstnd )
!-----------------------------------------------------------------------
!     take summation within local nodes:
!     node id's are given by idstnd
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension nfile(*)
dimension buffs(*), temp(*), buffr(*)
dimension idstnd(*)


!--- within myself
do i = 1, n
   temp(i) = buffs(i)
enddo

do isnd = 1, nodes - 1
!--- destination node
   idst = idstnd(isnd)

   nsd = n
   nrc = n
   call mespasd(300+isnd,buffs,buffr,nsd,nrc,idst+idroot,0)

   do i = 1, n
      temp(i) = temp(i) + buffr(i)
   enddo

!---------Internode synchronization
!         call gsync
enddo

do i = 1, n
   buffs(i) = temp(i)
enddo


return
end




subroutine rdtocd( nfile, myid, nodes, idroot,  &
& ar, mtrx, n, n1, n2, workarray, work3, nbncnt, nbndsp, idstnd )
!-----------------------------------------------------------------------
!    all to all communication 
!    to convert row decomposition to column decomposition
!
! (input)
!    ar : elements decomposed by row
!
! (output)
!    ar : elements decomposed by column
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension nfile(*)
dimension ar(mtrx,*)
dimension workarray(n1:n2,*), work3(*)
dimension nbndsp(*), nbncnt(*)
dimension idstnd(*)


if( n1 > n2 ) return

!--- within myself
do i = n1, n2
do j = n1, n2
   workarray(j,i) = ar(j,i)
enddo
enddo

do isnd = 1, nodes - 1
!--- destination node
   idst = idstnd(isnd)

!--- set No. of sending data (columns)
   is1 = nbndsp(idst+1) + 1
   is2 = nbndsp(idst+1) + nbncnt(idst+1)
   is3 = is2 - is1 + 1

!--- set No. of recieving data (row)
   ibs3 = nbncnt(idst+1)
   ib0  = nbndsp(idst+1)

   call rdtoc2( nfile, myid, nodes, idroot,  &
& is1, is2, is3, ibs3, idst, isnd,  &
& ar, mtrx, n, n1, n2, workarray(n1,ib0+1), work3 )

!---------Internode synchronization
!         call gsync
enddo


do i = 1, n
do j = n1, n2
   ar(j,i) = workarray(j,i)
enddo
enddo


return
end




subroutine rdtoc2( nfile, myid, nodes, idroot,  &
& is1, is2, is3, ibs3, idst, isnd,  &
& ar, mtrx, n, n1, n2, workarray, bufcr )
!-----------------------------------------------------------------------
!    all to all communication 
!    to convert row decomposition to column decomposition
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension nfile(*)
dimension ar(mtrx,*)
dimension workarray(n2-n1+1,*)
dimension bufcr(is3,*)


!--- set sending data
do i = n1, n2
   do j = is1, is2
      bufcr(j-is1+1,i-n1+1) = ar(j,i)
   enddo
enddo

nsd = is3*(n2-n1+1)
nrc = ibs3*(n2-n1+1)
call mespasd(500+isnd,bufcr,workarray,nsd,nrc,idst+idroot,0)  
!c---------myid.lt.idst: send & recv, if not empty
!          if (myid.lt.idst) then
!            if (nsd.ne.0) call cdsend(120,bufcr,nsd,idst+idroot,0)
!            if (nrc.ne.0) call cdrecv(130,gdcr,nrc,0)
!c---------myid.gt.idst: recv & send, if not empty
!          else
!            if (nrc.ne.0) call cdrecv(120,gdcr,nrc,0)
!            if (nsd.ne.0) call cdsend(130,bufcr,nsd,idst+idroot,0)
!          endif

return
end




