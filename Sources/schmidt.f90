



module ifcholesky
!-----------------------------------------------------------------------
! type declaration of shared variables in eigen.f
!-----------------------------------------------------------------------
implicit none

logical :: lstop = .true.
logical :: lcholesky = .true.

save

end module




module schmidt_timecnt
!-----------------------------------------------------------------------
! type declaration of shared variables in eigen.f
!-----------------------------------------------------------------------
implicit none

real*8,  dimension(11) :: tc = 0.d0
real*8  :: ctt0

save

end module




subroutine setlstop( nfile, myid, nodes, lstop_ )
!-----------------------------------------------------------------------
! set lcholesky_
!-----------------------------------------------------------------------
use ifcholesky
implicit none
integer :: nfile(*), myid, nodes
logical :: lstop_

lstop = lstop_

return
end




subroutine getlcholesky( nfile, myid, nodes, lcholesky_ )
!-----------------------------------------------------------------------
! set lcholesky_
!-----------------------------------------------------------------------
use ifcholesky
implicit none
integer :: nfile(*), myid, nodes
logical :: lcholesky_

lcholesky_ = lcholesky

return
end




subroutine schmidt( nfile, myid, nodes, ct0,  &
& cgjr_, iod_, iodg_, ltimecnt, lout, leig_start, leig_end )
!-----------------------------------------------------------------------
!   Gram-Schmidt orthonormalizaion by Cholesky decomposition
!
!   Wavefunctions are supposed to be decomposed by components.
!-----------------------------------------------------------------------
use param
use param_atom
use pwlda_pw
use pwlda_variables
use pwlda_proc
use pwlda_atom
implicit none
integer :: nfile(*), myid, nodes
real*8  :: ct0
real*8  :: cgjr_(2*npnod,nband)
integer :: iod_(nbnod), iodg_(nband)
logical :: ltimecnt, lout, leig_start, leig_end


!if( .not.lnoncollinear ) then
    call schmidt2( nfile, myid, nodes, ct0,  &
& lvand, cgjr_, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod,  &
& nbncnt, nbndsp, iod_, iodg_, prod, prodr, pdbuf,  &
& ltimecnt, lout, leig_start, leig_end, dmtrxr, nbxxxx,  &
& node_c, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, bijr, bm1r,  &
& rvol, ntype, nhk1_nat, nhk2_nat, natom, lvandi )
!else
!    !-----noncollinear magnetism
!    call schmidt_k2( nfile, myid, nodes, ct0,  &
!& lvand, cgjr_, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod,  &
!& nbncnt, nbndsp, iod_, iodg_, prod, prodr, prods, pdbuf,  &
!& ltimecnt, lout, leig_start, leig_end, dmtrxr, dmtrxi, nbxxxx,  &
!& node_c, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd,  &
!& bijr, biji, bm1r, bm1i,  &
!& rvol, ntype, nhk1_nat, nhk2_nat, natom, lvandi )
!end if


return
end




subroutine schmidt2( nfile, myid, nodes, ct0,  &
& lvand, cgjr, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod,  &
& nbncnt, nbndsp, iod, iodg, prod, prod2, dbuf,  &
& ltimecnt, lout, leig_start, leig_end,  &
& dmtrxr, nbxxxx,  &
& node_c, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, bijr, bm1r,  &
& rvol, ntype, nhk1, nhk2, natom, lvandi )
!-----------------------------------------------------------------------
!   Gram-Schmidt orthonormalizaion by Cholesky decomposition
!
!   Wavefunctions are supposed to be decomposed by components.
!-----------------------------------------------------------------------
use outfile
use ifcholesky
use schmidt_timecnt
implicit none
integer :: nfile(*), myid, nodes
real*8  :: ct0
logical :: lvand
integer :: npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod
real*8  :: cgjr(2*npnod,nband)
integer :: nbncnt(*), nbndsp(*)
integer :: iod(nbnod), iodg(nband)
real*8  :: prod(*), prod2(*), dbuf(*)
logical :: ltimecnt, lout, leig_start, leig_end
integer :: nbxxxx
real*8  :: dmtrxr(nbxxxx,*)
integer :: node_c, lbncnt(*), lbndsp(*), mbncnt(*), mbndsp(*), jdstnd(*)
real*8  :: bijr(*), bm1r(*)
real*8  :: rvol
integer :: ntype, natom
integer :: nhk1(ntype), nhk2(ntype)
logical :: lvandi(ntype)

!------declare local variables
integer :: i, j, ib, jb, npnodex, ig
real*8  :: tsum, cscr, prodjb, ar, prodr, dbuf1r
integer :: root, ierr
real*8  :: zero = 1.d-30
real*8  :: ct, timecnt
save zero


if( ltimecnt ) then
!          ct0 = timecnt()
if( leig_start ) then
    do i = 1, 11
       tc(i) = 0.d0
    end do
    ctt0 = ct0
end if
end if
!--- set S-matrix ------------------------------------------------------
do jb = 1, nband
do ib = 1, nbnod2 - nbnod1 + 1
   dmtrxr(ib,jb) = 0.d0
end do
end do

npnodex = 2*npnod
! --- This is serial code. ---
!      do i = 1, nband
!         do j = 1, i
!            produc = 0.d0
!            do m = 1, nmm
!               produc = produc + eigv(m,j)*eigv(m,i)
!            end do
!            dmtrxr(i,j) = produc*rdelv
!         end do
!      end do
! ----------------------------

root = 0
do i = 1, nband
#if PCHOLESKY
   if( i.gt.nbndsp(root+1)+nbncnt(root+1) ) root = root + 1
#endif
   do j = 1, i
      prod(j) = 0.d0
   end do
!CDIR OUTERUNROLL=5
!OCL UNROLL(9)
   do j = 1, i
      tsum = 0.d0
      do ig = 1, npnodex
!               prod(j) = prod(j) + cgjr(ig,j) * cgjr(ig,i)
         tsum = tsum + cgjr(ig,j) * cgjr(ig,i)
      end do
      prod(j) = tsum
#ifndef SX7
   end do
   do j = 1, i
#endif
      prod(j) = prod(j)*2.d0
   end do
   if( myid.eq.0 ) then
       do j = 1, i
          prod(j) = prod(j) - cgjr(1,j) * cgjr(1,i)
       end do
   end if
!   if( lvand ) then
!       do j = 1, i
!          call calcsc_atm( nfile, myid, nodes,  &
!& cscr, j, i, rvol, ntype, nhk1, nhk2, natom, lvandi )
!          prod(j) = prod(j) + cscr
!       end do
!   end if
#if PCHOLESKY
   call dsum(prod,i,prod2,root)
   if( myid.eq.root ) then
       do j = 1, i
          dmtrxr(i-nbnod1+1,j) = prod(j)
       end do
   end if
#else
   call gdsum(prod,i,prod2)
   do j = 1, i
      dmtrxr(i,j) = prod(j)
   end do
#endif
end do
!if( lvand ) call smat_vand( nfile, myid, nodes,  &
!& dmtrxr, nband, nbnod1, nbnod2, nbnod, bijr, bm1r, nbxxxx,  &
!& rvol, ntype, nhk1, nhk2, natom, lvandi )

if( ltimecnt ) then
    ct = timecnt()
    tc(1) = tc(1) + ct - ct0
    ct0 = ct
end if


!--- Inverse matrix of lower triangular matrix
!---     obtained by Cholesky decomposition of symmetric matrix
call icholesky( nfile, myid, nodes,  &
& dmtrxr, nband, nbnod1, nbnod2, nbxxxx, prod, prod2, nbncnt, nbndsp,  &
& node_c, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, bijr, bm1r )

!--- check whether the matrix is positive definite or not
do ib = 1, nband
   lcholesky = dmtrxr(ib,ib) > zero
   if( .not.lcholesky ) then
       if(loutfile(1)) write(nfile(1),*) '*** the matrix <phi_i|phi_j> is not positive definite.'
       if(loutfile(2)) write(nfile(2),*) '*** the matrix <phi_i|phi_j> is not positive definite.'

       !---emergency stop
       if( lstop ) call fstop( nfile, myid, nodes, 'in schmidt' )

       return
   end if
end do

if( ltimecnt ) then
    ct = timecnt()
    tc(2) = tc(2) + ct - ct0
    ct0 = ct
end if


!---TDDFT-FSSH: check wavefunction exchange
call tddft_fssh_check_trimat( nfile, myid, nodes,  &
& dmtrxr, nband, nbnod, nbxxxx, prod, prod2 )


#if PCHOLESKY
root = nodes - 1
#endif
do ib = nband, 1, -1
#if PCHOLESKY
   if( ib.le.nbndsp(root+1) ) root = root - 1
   if( myid.eq.root ) then
       do jb = 1, ib
          prod(jb) = dmtrxr(ib-nbnod1+1,jb)
       end do
   end if
   call dbcast(prod,ib,root)
#else
       do jb = 1, ib
          prod(jb) = dmtrxr(ib,jb)
       end do
#endif

   do ig = 1, npnodex
      dbuf(ig) = 0.d0
   end do
!OCL UNROLL(9)
   do jb = 1, ib
      prodjb = prod(jb)
      do ig = 1, npnodex
!               dbuf(ig) = dbuf(ig) + prod(jb)*cgjr(ig,jb)
         dbuf(ig) = dbuf(ig) + prodjb*cgjr(ig,jb)
      end do
   end do
   do ig = 1, npnodex
      cgjr(ig,ib) = dbuf(ig)
   end do
end do

!if( lvand ) call trans_slmi( nfile, myid, nodes,  &
!& dmtrxr, nband, nbnod1, nbnod2, nbnod, iod, iodg, nbxxxx,  &
!& ntype, nhk1, nhk2, natom, lvandi )


if( ltimecnt ) then
    ct = timecnt()
    tc(3) = tc(3) + ct - ct0
    ct0 = ct

    if( leig_end ) then
    do i = 1, 2
    if( loutfile(i) ) then
        write(nfile(i),*) ' Gram-Schmidt       '
        write(nfile(i),*) '                        set S-matrix',  &
&                     ' : cpu-time :', tc(1)
        write(nfile(i),*) '                    Cholesky decomp.',  &
&                     ' :          :', tc(2)
        write(nfile(i),*) '                            updating',  &
&                     ' :          :', tc(3)
        write(nfile(i),*) '                               total',  &
&                     ' : cpu-time :', ct-ctt0
    end if
    end do
    end if
end if

!check
if( lout ) then
if(loutfile(1)) write(nfile(1),*) ' check ortho-normalization'
if(loutfile(2)) write(nfile(2),*) ' check ortho-normalization'

do i = 1, nband
do j = i, nband
   ar = 0.d0
   do ig = 1, npnodex
      ar = ar + cgjr(ig,j) * cgjr(ig,i)
   end do
   if( myid.eq.0 ) then
       prodr = ar*2.d0 - cgjr(1,j) * cgjr(1,i)
     else
       prodr = ar*2.d0
   end if
!   if( lvand ) then
!       call calcsc_atm( nfile, myid, nodes,  &
!& cscr, i, j, rvol, ntype, nhk1, nhk2, natom, lvandi )
!       prodr = prodr + cscr
!   end if
   call gdsum(prodr,1,dbuf1r)
   if( i.eq.j ) prodr = prodr - 1.d0
   if( abs(prodr).gt.1.d-10 ) then
   if(loutfile(1)) write(nfile(1),*) 'orthonormalization error:',i, j, prodr
   if(loutfile(2)) write(nfile(2),*) 'orthonormalization error:',i, j, prodr
   end if
end do
end do

!-----Finalize the parallel environment
call end_parallel(ierr)
stop
end if

return
end




subroutine icholesky( nfile, myid, nodes,  &
& dmtrxr, nband, nbnod1, nbnod2, nbxxxx, prod, prod2,nbncnt,nbndsp,  &
& node_c, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, bijr, bm1r )
!-----------------------------------------------------------------------
!    Cholesky decomposition of symetric matrix
!                       and
!    Inverse matrix of lower triangular matrix
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension nfile(*)
dimension nbncnt(*), nbndsp(*)
dimension prod(*), prod2(*)
dimension dmtrxr(nbxxxx,*)
dimension lbncnt(*), lbndsp(*)
dimension mbncnt(*), mbndsp(*)
dimension jdstnd(*)
dimension bijr(*), bm1r(*)


!#if PCHOLESKY
!!      --- Cholesky decomposition ----------------------
!    call rdecmp( dmtrxr, nband, nbnod1, nbnod2, prod, nbxxxx,      &
!&                myid, nodes )
!!      --- inverse matrix  -----------------------------
!    call rinvcd( dmtrxr, nband, nbnod1, nbnod2, prod, prod2,  &
!&                nbncnt, nbndsp, nbxxxx, myid, nodes )
!#else
!if( node_c.eq.1 ) then
!      --- Cholesky decomposition ----------------------
    call rsdecmp( dmtrxr, nband, nbxxxx )
!      --- inverse matrix  -----------------------------
    call rsinvcd( dmtrxr, prod, nband, nbxxxx )
!  else
!    call vcholesky( nfile, myid, nodes,  &
!& dmtrxr, nband, nbxxxx, prod, prod2,  &
!& node_c, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, bijr, bm1r )
!end if
!#endif

return
end




subroutine rsdecmp( sr, n, nsize )
!-----------------------------------------------------------------------
!    Cholesky decomposition of symetric matrix ( sr )         '97/05/16
!                                                            F. Shimojo
!
!  ( input )
!     n        ......  order of matrix
!     sr(i,j)  ......  real      part of matrix
!       *** Diagonal and lower triangular elements have to be given. ***
!  ( output )
!     sr(i,j) ( i>=j )  ......  decomposed lower triangular matrix
!-----------------------------------------------------------------------
implicit real*8(a-h,o-z)
dimension  sr(nsize,*)
real*8  :: zero = 1.d-30
save zero

do i = 1, n
! ----- diagonal part -----
   do k = 1, i - 1
      sr(i,i) = sr(i,i) - sr(i,k)*sr(i,k)
   end do
   if( sr(i,i).gt.zero ) then
   sr(i,i) = sqrt(sr(i,i))
! ----- lower part -----
!             if( i-1.ge.n-i ) then
!                 do j = i + 1, n
!                 do k = 1, i - 1
!                    sr(j,i) = sr(j,i) - sr(j,k)*sr(i,k)
!                 end do
!                 end do
!               else
           do k = 1, i - 1
              srik = sr(i,k)
              do j = i + 1, n
!                       sr(j,i) = sr(j,i) - sr(j,k)*sr(i,k)
              sr(j,i) = sr(j,i) - sr(j,k)*srik
              end do
           end do
!             end if

       sriir = 1.d0/sr(i,i)
       do j = i + 1, n
          sr(j,i) = sr(j,i)*sriir
       end do
     else
       do j = i, n
          sr(j,i) = 0.d0
       end do
   end if
end do

! ----- upper part -----
do i = 1, n - 1
   do j = i + 1, n
      sr(i,j) = 0.d0
   end do
end do


return
end




subroutine rsinvcd( sr, workr, n, nsize )
!-----------------------------------------------------------------------
!    Inverse matrix of lower triangular matrix
!        obtained by Cholesky decomposition of symmetric matrix
!                                                            '97/05/16
!                                                            F. Shimojo
!
!  ( input )
!     n        ......  order of matrix
!     sr(i,j)  ......  real      part of matrix
!       *** Diagonal and lower triangular elements have to be given. ***
!  ( output )
!     sr(i,j) ( i>=j )  ......  inverse lower triangular matrix
!-----------------------------------------------------------------------
implicit real*8(a-h,o-z)
dimension  sr(nsize,*)
dimension  workr(*)
real*8  :: zero = 1.d-30
save zero


do i = 1, n
   if( abs(sr(i,i)).gt.zero ) then
       sriir = 1.d0/sr(i,i)
       do j = i + 1, n
          sr(j,i) = sr(j,i)*sriir
       end do
       sr(i,i) = sriir
   end if
end do

do j = n - 1, 1, -1
   do i = j + 1, n
      workr(i) = 0.d0
   end do
   do k = j + 1, n
      srkj = sr(k,j)
      do i = j + 1, n
!            do k = j + 1, i
!               workr(i) = workr(i) + sr(i,k)*sr(k,j)
         workr(i) = workr(i) + sr(i,k)*srkj
      end do
   end do
   do i = j + 1, n
      sr(i,j) = - workr(i)
   end do
end do


return
end




subroutine tmpdvn( nfile, myid, nodes,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp )
!-----------------------------------------------------------------------
!   index for decomposition
!-----------------------------------------------------------------------
dimension nfile(*)
dimension nbncnt(*), nbndsp(*)

if( nband.ge.nodes ) then

do inode = 0, nodes - 1
   call divno2( nband, nodes, inode, nbnod1, nbnod2, nbnod )
   nbncnt(inode+1) = nbnod
end do
nbndsp(1) = 0
do i = 2, nodes
   nbndsp(i) = nbndsp(i-1) + nbncnt(i-1)
end do

call divno2( nband, nodes, myid, nbnod1, nbnod2, nbnod )

else

do inode = 0, nband - 1
   nbncnt(inode+1) = 1
end do
do inode = nband, nodes - 1
   nbncnt(inode+1) = 0
end do
nbndsp(1) = 0
do i = 2, nodes
   nbndsp(i) = nbndsp(i-1) + nbncnt(i-1)
end do

if( myid .le. nband - 1 ) then
    nbnod1 = myid + 1
    nbnod2 = myid + 1
    nbnod  = 1
  else
    nbnod1 = nband + 1
    nbnod2 = nband
    nbnod  = 0
end if

end if

return
end




subroutine divno2( ndata, node, myid, ndnod1, ndnod2, ndnod )
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

ndnod  = ndata/node
ndnod1 = ndnod*myid
msrest = mod(ndata,node) - 1
if( myid.le.msrest ) then
    ndnod1 = ndnod1 + myid
    ndnod  = ndnod + 1
  else
    ndnod1 = ndnod1 + msrest + 1
end if
ndnod1 = ndnod1 + 1
ndnod2 = ndnod1 + ndnod - 1

return
end




