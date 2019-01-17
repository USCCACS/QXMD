



module engrad_variables
!-----------------------------------------------------------------------
! type declaration for engrad.f90
!-----------------------------------------------------------------------
implicit none

logical :: lnoncollinear

save

end module




subroutine set_lnoncollinear_in_engrad( nfile, myid, nodes, lnoncollinear_ )
!-----------------------------------------------------------------------
!     allocate memory for planewave-decomposition variables
!-----------------------------------------------------------------------
use engrad_variables
implicit none
integer :: nfile(*), myid, nodes
logical :: lnoncollinear_

lnoncollinear = lnoncollinear_

return
end subroutine




subroutine ekinetic( nfile, myid, nodes,  &
& ekint, ekib, occ, cgjr, rhcr, nplwex, nplw,  &
& nband, nbnod1, nbnod2, nbnod, iod, dkgnrm, dkrecnrm, iflag )
!-----------------------------------------------------------------------
!    kinetic energy & HC products
!
!    iflag = 0    : obtain HC products
!    iflag = else : obtain both HC products and kinetic energy
!
!-----------------------------------------------------------------------
use engrad_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: ekint, ekib(*), occ(*)
integer :: nplwex, nplw
real*8  :: cgjr(*), rhcr(*)
integer :: nband, nbnod1, nbnod2, nbnod, iod(*), iflag
real*8  :: dkgnrm(*), dkrecnrm(*)


!if( .not.lnoncollinear ) then
    call ekinetic2( nfile, myid, nodes,  &
& ekint, ekib, occ, cgjr, rhcr, nplwex, nplw,  &
& nband, nbnod1, nbnod2, nbnod, iod, dkgnrm, iflag )
!else
!    !-----noncollinear magnetism
!    call ncekinetic2( nfile, myid, nodes,  &
!& ekint, ekib, occ, cgjr, rhcr, nplwex, nplw,  &
!& nband, nbnod1, nbnod2, nbnod, iod, dkrecnrm, iflag )
!end if


return
end subroutine




subroutine ekinetic2( nfile, myid, nodes,  &
& ekint, ekib, occ, cgjr, rhcr, nplwex, nplw,  &
& nband, nbnod1, nbnod2, nbnod, iod, dkgnrm, iflag )
!-----------------------------------------------------------------------
!    kinetic energy & HC products
!
!    iflag = 0    : obtain HC products
!    iflag = else : obtain both HC products and kinetic energy
!
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
real*8  :: ekint, ekib(*), occ(*)
integer :: nplwex, nplw
real*8  :: cgjr(nplwex,*), rhcr(nplwex,*)
integer :: nband, nbnod1, nbnod2, nbnod, iod(*), iflag
real*8  :: dkgnrm(*)

!------declare local variables
integer :: i, ib, ig
real*8  :: ekin, ekinib, dbuf1r


!--- clear gradients
!      do i = 1, nbnod
!         do ig = 1, nplwex
!            rhcr(ig,i) = 0.d0
!         end do
!      end do

do i = 1, nbnod
   rhcr(1:nplwex,i) = dkgnrm(1:nplwex)*cgjr(1:nplwex,i)
end do

if( iflag.eq.0 ) return
ekin = 0.d0
do i = 1, nbnod
   ib = iod(i)
   ekinib = 0.d0
   do ig = 1, nplwex
      ekinib = ekinib + rhcr(ig,i)*cgjr(ig,i)
   end do
   ekinib = 2.d0*ekinib - rhcr(1,i)*cgjr(1,i)
   ekin = ekin + occ(ib)*ekinib
   ekib(i) = 2.d0*ekinib
end do
call gdsum(ekin,1,dbuf1r)

ekint = ekint + ekin


return
end




subroutine enonloc( nfile, myid, nodes,  &
& lspin, nspin, lpaw, enclt, enclvt, occ, glocal, cgjr, rhcr, hcsr,  &
& nplwex, nplw, nband, nbnod1, nbnod2, nbnod,  &
& iod, mfd2ft, ntotfd, nd1vks, kfft0, rvol,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh,  &
& eigr, eigi, apk, apki, scwr, scwi,  &
& ntype, nhk1, nhk2, iatoit,  &
& lkbpp_r, lvand_r, lkbpp_g, lvand_g, lkbppi, lvandi, lking,  &
& rdelv, iflag, lmax, mxl, lchk, lclno,  &
& nylmmx, natom, ycos, ysin, nplwcs, tc, ltimecnt )
!-----------------------------------------------------------------------
!    nonlocal pseudopotential energy and local/nonlocal HC products
!
!    iflag = 0    : obtain HC products
!    iflag = else : obtain both HC products and energy
!
!-----------------------------------------------------------------------
use engrad_variables
use ncmagne_variables
implicit none
integer :: nfile(*), myid, nodes
logical :: lspin
integer :: nspin
logical :: lpaw
real*8  :: enclt, enclvt
real*8  :: occ(*), glocal(*)
integer :: nplwex, nplw, nband, nbnod1, nbnod2, nbnod
real*8  :: cgjr(nplwex,*), rhcr(nplwex,*), hcsr(*)
integer :: iod(*)
integer :: mfd2ft(*), ntotfd, nd1vks(*), kfft0
real*8  :: rvol
real*8  :: fft3x(*), fft3y(*)
complex*16 :: fftwork(*)
integer :: ngenh
integer :: ijkgd(-ngenh:ngenh)
real*8  :: eigr(*), eigi(*), apk(*), apki(*), scwr(*), scwi(*)
integer :: ntype, nhk1(*), nhk2(*), iatoit(*)
logical :: lkbpp_r, lvand_r, lkbpp_g, lvand_g
logical :: lkbppi(ntype), lvandi(ntype), lking(ntype)
real*8  :: rdelv
integer :: iflag, lmax(*), mxl, lclno(*)
logical :: lchk(0:mxl,*)
integer :: nylmmx, natom
real*8  :: ycos(*), ysin(*)
integer :: nplwcs
real*8  :: tc(*)
logical :: ltimecnt


!if( .not.lnoncollinear ) then
    call enonloc2( nfile, myid, nodes,  &
& lspin, nspin, lpaw, enclt, enclvt, occ, glocal, cgjr, rhcr, hcsr,  &
& nplwex, nplw, nband, nbnod1, nbnod2, nbnod,  &
& iod, mfd2ft, ntotfd, nd1vks, kfft0, rvol,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, eigr, apk,  &
& ntype, nhk1, nhk2, iatoit,  &
& lkbpp_r, lvand_r, lkbpp_g, lvand_g, lkbppi, lvandi, lking,  &
& rdelv, iflag, lmax, mxl, lchk, lclno,  &
& nylmmx, natom, ycos, ysin, nplwcs, tc, ltimecnt )
!else
!    !-----noncollinear magnetism
!    call ncenonloc2( nfile, myid, nodes,  &
!& lspin, nspin, lpaw, enclt, enclvt, occ, glocal, cgjr, rhcr, hcsr,  &
!& nplwex, nplw, nband, nbnod1, nbnod2, nbnod,  &
!& iod, mfd2ft, ntotfd, nd1vks, kfft0, rvol,  &
!& fft3x, fft3y, fftwork, ncijkg,  &
!& eigr, eigi, eig2r, eig2i, apk, apki, scwr, scwi,  &
!& ntype, nhk1, nhk2, iatoit,  &
!& lkbpp_r, lvand_r, lkbpp_g, lvand_g, lkbppi, lvandi, lking,  &
!& rdelv, iflag, lmax, mxl, lchk, lclno,  &
!& nylmmx, natom, ycos, ysin, nplwcs, tc, ltimecnt )
!end if


return
end




subroutine enonloc2( nfile, myid, nodes,  &
& lspin, nspin, lpaw, enclt, enclvt, occ, glocal, cgjr, rhcr, hcsr,  &
& nplwex, nplw, nband, nbnod1, nbnod2, nbnod,  &
& iod, mfd2ft, ntotfd, nd1vks, kfft0, rvol,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, eigr, apk,  &
& ntype, nhk1, nhk2, iatoit,  &
& lkbpp_r, lvand_r, lkbpp_g, lvand_g, lkbppi, lvandi, lking,  &
& rdelv, iflag, lmax, mxl, lchk, lclno,  &
& nylmmx, natom, ycos, ysin, nplwcs, tc, ltimecnt )
!-----------------------------------------------------------------------
!    nonlocal pseudopotential energy and local/nonlocal HC products
!
!    iflag = 0    : obtain HC products
!    iflag = else : obtain both HC products and energy
!
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
logical :: lspin
integer :: nspin
logical :: lpaw
real*8  :: enclt, enclvt
real*8  :: occ(*), glocal(*)
integer :: nplwex, nplw, nband, nbnod1, nbnod2, nbnod
real*8  :: cgjr(nplwex,*), rhcr(nplwex,*), hcsr(*)
integer :: iod(*)
integer :: mfd2ft(*), ntotfd, nd1vks(*), kfft0
real*8  :: rvol
real*8  :: fft3x(*), fft3y(*)
complex*16 :: fftwork(*)
integer :: ngenh
integer :: ijkgd(-ngenh:ngenh)
real*8  :: apk(*), eigr(0:*)
integer :: ntype, nhk1(*), nhk2(*), iatoit(*)
logical :: lkbpp_r, lvand_r, lkbpp_g, lvand_g
logical :: lkbppi(ntype), lvandi(ntype), lking(ntype)
real*8  :: rdelv
integer :: iflag, lmax(*), mxl, lclno(*)
logical :: lchk(0:mxl,*)
integer :: nylmmx, natom
real*8  :: ycos(*), ysin(*)
integer :: nplwcs
real*8  :: tc(*)
logical :: ltimecnt

!------declare local variables
integer :: i, ib, nhcstt, j
real*8  :: encl, enclib, enclv, enclw, dbuf1r
real*8  :: ct0, ct, timecnt


!-----( H(local potenital) + H(nonlocal pp in real space) )*|phi>
encl  = 0.d0
do i = 1, nbnod
   ib = iod(i)

   !-----set pointer for hcsr
!   if( lvand_r ) then
!       nhcstt = nplwex*(i-1)+1
!       !-----set slmir
!       call setbackslm( nfile, myid, nodes,  &
!& i, ntype, nhk1, nhk2, natom, lvandi )
!     else
       nhcstt = 1
!   end if

   call hcnonloc( nfile, myid, nodes,  &
& enclib, glocal, cgjr(1,i), rhcr(1,i), hcsr(nhcstt), nplwex, nplw,  &
& nband, nbnod1, nbnod2, nbnod, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, eigr, apk,  &
& ntype, nhk1, nhk2, natom, iatoit,  &
& lkbpp_r, lvand_r, lkbppi, lvandi, lking, rdelv, iflag,  &
& lmax, mxl, lchk, lclno, nylmmx, tc, ltimecnt, .false. )

!check
!dbuf1r = 0.d0
!do j = 1, nplwex
!   dbuf1r = dbuf1r + cgjr(j,i)**2
!end do
!dbuf1r = 2.d0*dbuf1r - cgjr(1,i)**2
!if( myid == 0 ) write(*,*) i, dbuf1r
!dbuf1r = 0.d0
!do j = 1, nplwex
!   dbuf1r = dbuf1r + rhcr(j,i)**2
!end do
!dbuf1r = 2.d0*dbuf1r - rhcr(1,i)**2
!if( myid == 0 ) write(*,*) i, dbuf1r
!dbuf1r = 0.d0
!do j = 1, nplwex
!   dbuf1r = dbuf1r + hcsr(nhcstt+j-1)**2
!end do
!dbuf1r = 2.d0*dbuf1r - hcsr(nhcstt)**2
!if( myid == 0 ) write(*,*) i, dbuf1r
!check
   !-----norm-conserving pp. in real space
   if( lkbpp_r ) then
        if( iflag.ne.0 ) encl = encl + occ(ib)*enclib

        !------store sumr
        call savesumr( nfile, myid, nodes,  &
& i, ntype, nhk1, nhk2, natom, lkbppi, lking, nylmmx, nspin )
   end if

end do
if( lkbpp_r .and. iflag.ne.0 ) then
    call gdsum(encl,1,dbuf1r)
    encl = encl/dble(kfft0)
    enclt = enclt + encl
end if


if( ltimecnt ) then
    ct0  = timecnt()
end if

!-----( H(nonlocal pp in reciprocal space) )*|phi>
enclv = 0.d0

!-----norm-conseving pseudopotentials in KB form
if( lkbpp_g ) then

    call hckbpp_g( nfile, myid, nodes,  &
& cgjr, rhcr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod,  &
& ntype, nhk1, nhk2, natom, lkbppi, lking, ycos, ysin, nplwcs )

    !------store slmir
    call savesumr_g( nfile, myid, nodes,  &
& nbnod, ntype, nhk1, nhk2, natom, lkbppi, lking, lspin, nspin )

    if( iflag.ne.0 ) then
        encl = 0.d0
        call enkbpp_g( nfile, myid, nodes,  &
& encl, nband, nbnod1, nbnod2, nbnod, occ, iod,  &
& ntype, nhk1, nhk2, natom, lkbppi, lking, rvol )

        enclv = enclv + encl
    end if

end if


!-----ultrasoft pseudopotentials
!if( lvand_g .or. lvand_r ) then

!    if( lvand_g )  &
!&       call hcvand_g( nfile, myid, nodes,  &
!& lspin, nspin, cgjr, rhcr, hcsr, nplwex, nplw,  &
!& nband, nbnod1, nbnod2, nbnod,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking, rvol,  &
!& ycos, ysin, nplwcs )

!    if( iflag.ne.0 ) then
!        encl  = 0.d0
!        enclw = 0.d0
!        call envand_g( nfile, myid, nodes,  &
!& lpaw, nspin, encl, enclw, nband, nbnod1, nbnod2, nbnod, occ, iod,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking, rvol )
!
!        enclv = enclv + encl + enclw
!    end if

!end if

if( ltimecnt ) then
    ct = timecnt()
    tc(2) = tc(2) + ct - ct0
    ct0 = ct
end if

enclt = enclt + enclv
enclvt = enclvt + enclv


!check
!      call get_worldqm( myid_, nodes_ )
!      if( myid_ == 0 ) then
!      write(*,*) 'encl:', enclt, enclv
!      do ib = 1, nbnod
!         sum1 = 0.d0
!         sum2 = 0.d0
!         nhcstt = nplwex*(ib-1)+1
!         do ig = 1, nplwex
!            sum1 = sum1 + rhcr(ig,ib)*rhcr(ig,ib)
!            sum2 = sum2 + hcsr(nhcstt+ig-1)*hcsr(nhcstt+ig-1)
!         end do
!         sum1 = 2.d0*sum1 - rhcr(1,ib)*rhcr(1,ib)
!         sum2 = 2.d0*sum2 - hcsr(nhcstt)*hcsr(nhcstt)
!         write(*,*) ib, sum1, sum2
!      end do
!      end if
!      call fstop( nfile, myid, nodes, 'enonloc' )


return
end




subroutine hcnonloc( nfile, myid, nodes,  &
& enclib, glocal, cgjr, rhcr, hcsr, nplwex, nplw,  &
& nband, nbnod1, nbnod2, nbnod, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, eigr, apk,  &
& ntype, nhk1, nhk2, natom, iatoit,  &
& lkbpp_r, lvand_r, lkbppi, lvandi, lking, rdelv, iflag,  &
& lmax, mxl, lchk, lclno, nylmmx, tc, ltimecnt, lcslmir )
!-----------------------------------------------------------------------
!    nonlocal pseudopotential energy and local/nonlocal HC products
!
!    iflag = 0    : obtain HC products
!    iflag = else : obtain both HC products and energy
!
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
real*8  :: enclib, glocal(*)
integer :: nplwex, nplw
real*8  :: cgjr(nplwex), rhcr(nplwex), hcsr(*)
integer :: nband, nbnod1, nbnod2, nbnod
integer :: mfd2ft(*), ntotfd, nd1vks(*), kfft0
real*8  :: fft3x(*), fft3y(*)
complex*16 :: fftwork(*)
integer :: ngenh
integer :: ijkgd(-ngenh:ngenh)
real*8  :: apk(*), eigr(0:*)
integer :: ntype, natom
integer :: nhk1(*), nhk2(*), iatoit(*)
logical :: lkbpp_r, lvand_r
logical :: lkbppi(ntype), lvandi(ntype), lking(ntype)
real*8  :: rdelv
integer :: iflag, lmax(*), mxl, lclno(*)
logical :: lchk(0:mxl,*)
integer :: nylmmx
real*8  :: tc(*)
logical :: ltimecnt
logical :: lcslmir

!------declare local variables
integer :: iloopcnt, nml, inv, kfft1, kfft2, kfft3, ifd, ift, ig, ijk
integer :: ierrft
real*8  :: ct0, ct, timecnt


#ifndef VECTOR
iloopcnt = ntotfd
#else
iloopcnt = kfft0
#endif

if( ltimecnt ) then
    ct0  = timecnt()
end if

call wvg2r2( nfile, myid, nodes,  &
& cgjr, eigr, nplwex, nplw, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh )

if( ltimecnt ) then
    ct = timecnt()
    tc(1) = tc(1) + ct - ct0
    ct0 = ct
end if

!   --- constants for fft ---
nml = 1
inv = 2
kfft1 = nd1vks(1)
kfft2 = nd1vks(2)
kfft3 = nd1vks(3)


!-----clear gradients from nonlocal and local potentials
fft3x(1:kfft0) = 0.d0
fft3y(1:kfft0) = 0.d0

!--- gradients from local potentials
#ifndef VECTOR
do ifd = 1, ntotfd
   ift = mfd2ft(ifd)
   fft3x(ift) = glocal(ifd)*eigr(ifd)
end do
#else
fft3x(1:kfft0) = glocal(1:kfft0)*eigr(1:kfft0)
#endif


if( lkbpp_r ) then

    !-----norm-conserving pp. in real space
    apk(1:iloopcnt) = 0.d0
    call calnlc( nfile, myid, nodes,  &
& ntype, lmax, mxl, lchk, lclno, apk, eigr, nhk1, nhk2,  &
& lkbppi, lking, iatoit, rdelv, nylmmx, nd1vks  )

    if( iflag.ne.0 ) then
        !--- nonlocal pp energy
        enclib = 0.d0
        do ifd = 1, iloopcnt
           enclib = enclib + apk(ifd)*eigr(ifd)
        end do
    end if

    if( ltimecnt ) then
        ct = timecnt()
        tc(2) = tc(2) + ct - ct0
        ct0 = ct
    end if

    !--- gradients from nonlocal potentials
#ifndef VECTOR
    do ifd = 1, ntotfd
       ift = mfd2ft(ifd)
       fft3x(ift) = fft3x(ift) + apk(ifd)
    end do
#else
    fft3x(1:kfft0) = fft3x(1:kfft0) + apk(1:kfft0)
#endif

end if


!if( lvand_r ) then
!
!    !-----ultrasoft pp. in real space
!    apk(1:iloopcnt) = 0.d0
!    call hcvand_ib_r( nfile, myid, nodes,  &
!& apk, fft3y, eigr, ntype, nhk1, nhk2, natom, lvandi, lking,  &
!& iatoit, rdelv, lcslmir )
!
!    if( ltimecnt ) then
!        ct = timecnt()
!        tc(2) = tc(2) + ct - ct0
!        ct0 = ct
!    end if
!
!    !--- gradients from ultrasoft potentials
!#ifndef VECTOR
!    do ifd = 1, ntotfd
!       ift = mfd2ft(ifd)
!       fft3x(ift) = fft3x(ift) + apk(ifd)
!    end do
!#else
!    fft3x(1:kfft0) = fft3x(1:kfft0) + apk(1:kfft0)
!#endif
!
!    !-----store SC products
!    apk(1:iloopcnt) = fft3y(1:iloopcnt)
!    fft3y(1:kfft0) = 0.d0
!
!end if


call rfft3( nml, fft3x, fft3y, fftwork,  &
&                kfft1, kfft2, kfft3, kfft0, ierrft )

do ig = 0, nplw
   ijk  = ijkgd(ig)
   rhcr(2*ig+1) = rhcr(2*ig+1) + fft3x(ijk)
   rhcr(2*ig+2) = rhcr(2*ig+2) + fft3y(ijk)
end do


!if( lvand_r ) then
!
!    fft3x(1:kfft0) = 0.d0
!    fft3y(1:kfft0) = 0.d0
!#ifndef VECTOR
!    do ifd = 1, ntotfd
!       ift = mfd2ft(ifd)
!       fft3x(ift) = apk(ifd)
!    end do
!#else
!    fft3x(1:kfft0) = apk(1:kfft0)
!#endif
!    call rfft3( nml, fft3x, fft3y, fftwork,  &
!&                    kfft1, kfft2, kfft3, kfft0, ierrft )
!
!    do ig = 0, nplw
!       ijk  = ijkgd(ig)
!       hcsr(2*ig+1) = hcsr(2*ig+1) + fft3x(ijk)
!       hcsr(2*ig+2) = hcsr(2*ig+2) + fft3y(ijk)
!    end do
!
!end if


if( ltimecnt ) then
    ct = timecnt()
    tc(3) = tc(3) + ct - ct0
    ct0 = ct
end if


return
end




subroutine pwpgfd2ft( nfile, myid, nodes,  &
& glocal, apk, mfd2ft, ntotfd, kfft0d )
!-----------------------------------------------------------------------
!     convert FD meshes to FFT meshes
!-----------------------------------------------------------------------
use engrad_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: glocal(*), apk(*)
integer :: mfd2ft(*), ntotfd, kfft0d


!if( .not.lnoncollinear ) then
    call gfd2ft( nfile, myid, nodes,  &
& glocal, apk, mfd2ft, ntotfd, kfft0d )
!else
!    !-----noncollinear magnetism
!    call ncgfd2ft( nfile, myid, nodes,  &
!& glocal, apk, mfd2ft, ntotfd, kfft0d )
!end if


return
end




subroutine gfd2ft( nfile, myid, nodes,  &
& glocal, apk, mfd2ft, ntotfd, kfft0d )
!-----------------------------------------------------------------------
!     convert FD meshes to FFT meshes
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
real*8  :: glocal(*), apk(*)
integer :: mfd2ft(*), ntotfd, kfft0d

!------declare local variables
integer :: ifd, ift


apk(1:ntotfd) = glocal(1:ntotfd)
glocal(1:kfft0d) = 0.d0
do ifd = 1, ntotfd
   ift = mfd2ft(ifd)
   glocal(ift) = apk(ifd)
end do


return
end




subroutine pwplocal_on_coarse( nfile, myid, nodes,  &
& glocal, thrhgr, nplw5ex, nplw5, nplw2ex, nplw2,  &
& kfft1d, kfft2d, kfft3d, kfft0d, mfd2ft, ntotfd, ijkgd, ngenh,  &
& kfft1,  kfft2,  kfft3,  kfft0,  mfd2ft_c, ntotfd_c, ijkg,  &
& fft3x, fft3y, fftwork )
!-----------------------------------------------------------------------
!     glocal on dense grid -> coarse grid
!-----------------------------------------------------------------------
use engrad_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: glocal(*)
integer :: nplw5ex, nplw5, nplw2ex, nplw2
real*8  :: thrhgr(nplw5ex)
integer :: kfft1d, kfft2d, kfft3d, kfft0d
integer :: ntotfd
integer :: mfd2ft(ntotfd)
integer :: ngenh
integer :: ijkgd(-ngenh:ngenh)
integer :: kfft1, kfft2, kfft3, kfft0
integer :: ntotfd_c
integer :: mfd2ft_c(ntotfd_c)
integer :: ijkg(-nplw2:nplw2)
real*8  :: fft3x(*), fft3y(*)
complex*16 :: fftwork(*)


!if( .not.lnoncollinear ) then
    call local_on_coarse( nfile, myid, nodes,  &
& glocal, thrhgr, nplw5ex, nplw5, nplw2ex, nplw2,  &
& kfft1d, kfft2d, kfft3d, kfft0d, mfd2ft, ntotfd, ijkgd, ngenh,  &
& kfft1,  kfft2,  kfft3,  kfft0,  mfd2ft_c, ntotfd_c, ijkg,  &
& fft3x, fft3y, fftwork )
!else
!    !-----noncollinear magnetism
!    call nclocal_on_coarse( nfile, myid, nodes,  &
!& glocal, thrhgr, nplw5ex, nplw5, nplw2ex, nplw2,  &
!& kfft1d, kfft2d, kfft3d, kfft0d, mfd2ft, ntotfd, ijkgd, ngenh,  &
!& kfft1,  kfft2,  kfft3,  kfft0,  mfd2ft_c, ntotfd_c, ijkg,  &
!& fft3x, fft3y, fftwork )
!end if


return
end




subroutine local_on_coarse( nfile, myid, nodes,  &
& glocal, thrhgr, nplw5ex, nplw5, nplw2ex, nplw2,  &
& kfft1d, kfft2d, kfft3d, kfft0d, mfd2ft, ntotfd, ijkgd, ngenh,  &
& kfft1,  kfft2,  kfft3,  kfft0,  mfd2ft_c, ntotfd_c, ijkg,  &
& fft3x, fft3y, fftwork )
!-----------------------------------------------------------------------
!     glocal on dense grid -> coarse grid
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
real*8  :: glocal(*)
integer :: nplw5ex, nplw5, nplw2ex, nplw2
real*8  :: thrhgr(nplw5ex)
integer :: kfft1d, kfft2d, kfft3d, kfft0d
integer :: ntotfd
integer :: mfd2ft(ntotfd)
integer :: ngenh
integer :: ijkgd(-ngenh:ngenh)
integer :: kfft1, kfft2, kfft3, kfft0
integer :: ntotfd_c
integer :: mfd2ft_c(ntotfd_c)
integer :: ijkg(-nplw2:nplw2)
real*8  :: fft3x(*), fft3y(*)
complex*16 :: fftwork(*)

!------declare local variables
integer :: nml, inv, ifd, ift, ig, ijk, ijkm
integer :: ierrft


!   --- constants for fft ---
nml = 1
inv = 2

fft3x(1:kfft0d) = 0.d0
fft3y(1:kfft0d) = 0.d0
#ifndef VECTOR
do ifd = 1, ntotfd
   ift = mfd2ft(ifd)
   fft3x(ift) = glocal(ifd)
end do
#else
fft3x(1:kfft0d) = glocal(1:kfft0d)
#endif


if( kfft0d /= kfft0 ) then
    call rfft3( nml, fft3x, fft3y, fftwork,  &
&                    kfft1d, kfft2d, kfft3d, kfft0d, ierrft )

    do ig = 0, nplw2
       ijk  = ijkgd(ig)
       thrhgr(2*ig+1) = fft3x(ijk)
       thrhgr(2*ig+2) = fft3y(ijk)
    end do

    fft3x(1:kfft0) = 0.d0
    fft3y(1:kfft0) = 0.d0
    do ig = 0, nplw2
       ijk  = ijkg(ig)
       ijkm = ijkg(-ig)
       fft3x(ijk)  = thrhgr(2*ig+1)
       fft3y(ijk)  = thrhgr(2*ig+2)
       fft3x(ijkm) =   thrhgr(2*ig+1)
       fft3y(ijkm) = - thrhgr(2*ig+2)
    end do

    call rfft3( inv, fft3x, fft3y, fftwork,  &
&                    kfft1, kfft2, kfft3, kfft0, ierrft )
end if


#ifndef VECTOR
do ifd = 1, ntotfd_c
   ift = mfd2ft_c(ifd)
   glocal(ifd) = fft3x(ift)
end do
#else
glocal(1:kfft0) = fft3x(1:kfft0)
#endif


return
end




subroutine hmatrix( nfile, myid, nodes,  &
& eig, cgjr, rhcr, npnod1, npnod2, npnod,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& iod, bx0r, dmtrxr, dmtrxi, nbxxxx, prod, prodr )
!-----------------------------------------------------------------------
!    Hamiltonian matrix 
!-----------------------------------------------------------------------
use engrad_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: eig(*)
integer :: npnod1, npnod2, npnod
integer :: nband, nbnod1, nbnod2, nbnod, nbncnt(*), nbndsp(*)
real*8  :: cgjr(2*npnod,*), rhcr(2*npnod,*)
integer :: iod(*)
real*8  :: bx0r(nbnod,*)
integer :: nbxxxx
real*8  :: dmtrxr(nbxxxx,*), dmtrxi(*)
real*8  :: prod(*), prodr(*)


!if( .not.lnoncollinear ) then
    call hmatrix2( nfile, myid, nodes,  &
& eig, cgjr, rhcr, npnod1, npnod2, npnod,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& iod, bx0r, dmtrxr, nbxxxx, prod, prodr )
!else
!    !-----noncollinear magnetism
!    call hmatrix_k2( nfile, myid, nodes,  &
!& eig, cgjr, rhcr, npnod1, npnod2, npnod,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& iod, bx0r, dmtrxr, dmtrxi, nbxxxx, prod, prodr )
!end if


return
end




subroutine hmatrix2( nfile, myid, nodes,  &
& eig, cgjr, rhcr, npnod1, npnod2, npnod,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& iod, bx0r, dmtrxr, nbxxxx, prod, prodr )
!-----------------------------------------------------------------------
!    Hamiltonian matrix 
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
real*8  :: eig(*)
integer :: npnod1, npnod2, npnod
integer :: nband, nbnod1, nbnod2, nbnod, nbncnt(*), nbndsp(*)
real*8  :: cgjr(2*npnod,*), rhcr(2*npnod,*)
integer :: iod(*)
real*8  :: bx0r(nbnod,*)
integer :: nbxxxx
real*8  :: dmtrxr(nbxxxx,*)
real*8  :: prod(*), prodr(*)

!------declare local variables
integer :: ib, jb, ig, myid_, nodes_
real*8  :: prod1


do ib = 1, nband
   do jb = 1, ib
      prod1 = 0.0
      do ig = 1, 2*npnod
         prod1 = prod1 + rhcr(ig,ib)*cgjr(ig,jb)
      end do
      if( myid.eq.0 ) then
          prod1 = 2.d0*prod1 - rhcr(1,ib)*cgjr(1,jb)
        else
          prod1 = 2.d0*prod1
      end if
      prod(jb) = prod1
   end do
   call gdsum(prod,ib,prodr)
   do jb = 1, ib
#if PCHOLESKY
      if( jb.ge.nbnod1 .and. jb.le.nbnod2 )  &
&         dmtrxr(jb-nbnod1+1,ib) = prod(jb)
      if( ib.ge.nbnod1 .and. ib.le.nbnod2 )  &
&         dmtrxr(ib-nbnod1+1,jb) = prod(jb)
#else
      dmtrxr(jb,ib) = prod(jb)
      dmtrxr(ib,jb) = prod(jb)
#endif
   end do
   eig(ib) = prod(ib)
end do

!check
!      call get_worldqm( myid_, nodes_ )
!      if( myid_ == 0 ) then
!      do ib = 1, nband
!      do jb = 1, nband
!         write(*,*) ib, jb, dmtrxr(ib,jb) !, dmtrxi(ib,jb)
!      end do
!      end do
!      end if
!      call fstop( nfile, myid, nodes, 'hmatrix' )

!-8/22/00#if PCHOLESKY
!-8/22/00      do i = 1, nbnod
!-8/22/00         ib = iod(i)
!-8/22/00         prod(i) = dmtrxr(i,ib)
!-8/22/00      end do
!-8/22/00      call alldgatherv(prod,nbnod,prodr,nbncnt,nbndsp)
!-8/22/00      do ib = 1, nband
!-8/22/00c         ib = iodg(i)
!-8/22/00         eig(ib) = prodr(ib)
!-8/22/00      end do
!-8/22/00#else
!-8/22/00c--- unify H matrix
!-8/22/00      if( myid.eq.0 ) then
!-8/22/00          i0 = 0
!-8/22/00          do i = 1, nbnod
!-8/22/00             do jb = 1, nband
!-8/22/00                dmtrxr(i+i0,jb) = bx0r(i,jb)
!-8/22/00             end do
!-8/22/00          end do
!-8/22/00          do node_id = 1, nodes - 1
!-8/22/00             i0 = i0 + nbncnt(node_id)
!-8/22/00             call rcvin1( i0, nbncnt(node_id+1), nband, bx0r, dmtrxr )
!-8/22/00c             nrc = nbncnt(node_id+1)*nband
!-8/22/00c             call cdrecv(100,bx0r,nrc,0)
!-8/22/00c             do i = 1, nbncnt(node_id+1)
!-8/22/00c                do jb = 1, nband
!-8/22/00c                   dmtrxr(i+i0,jb) = bx0r(i,jb)
!-8/22/00c                end do
!-8/22/00c             end do
!-8/22/00             call gsync
!-8/22/00          end do
!-8/22/00        else
!-8/22/00          do node_id = 1, nodes - 1
!-8/22/00             if( myid.eq.node_id ) then
!-8/22/00                 nsd = nbnod*nband
!-8/22/00                 call cdsend(100,bx0r,nsd,0,0)
!-8/22/00             end if
!-8/22/00             call gsync
!-8/22/00          end do
!-8/22/00      end if
!-8/22/00      call dbcast(dmtrxr,nband*nband,0)
!-8/22/00      do ib = 1, nband
!-8/22/00         eig(ib) = dmtrxr(ib,ib)
!-8/22/00      end do
!-8/22/00#endif


return
end


!-8/22/00      subroutine rcvin1( i0, nbncnt, nband, dmtrxr, dmtrx2 )
!-8/22/00c-----------------------------------------------------------------------
!-8/22/00      implicit real*8 ( a-h, o-z )
!-8/22/00      dimension dmtrxr(nbncnt,*)
!-8/22/00      dimension dmtrx2(nband,*)
!-8/22/00
!-8/22/00             nrc = nbncnt*nband
!-8/22/00             call cdrecv(100,dmtrxr,nrc,0)
!-8/22/00             do i = 1, nbncnt
!-8/22/00                do jb = 1, nband
!-8/22/00                   dmtrx2(i+i0,jb) = dmtrxr(i,jb)
!-8/22/00                end do
!-8/22/00             end do
!-8/22/00
!-8/22/00      return
!-8/22/00      end




subroutine gradres( nfile, myid, nodes,  &
& lvand, eig, occ, cgjr, rhcr, hcsr, npnod1, npnod2, npnod,  &
& nband, zanst1, zanst2, prod, prodr )
!-----------------------------------------------------------------------
!    gradients and residuals
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
logical :: lvand
real*8  :: eig(*), occ(*)
integer :: npnod1, npnod2, npnod
integer :: nband
real*8  :: cgjr(2*npnod,*), rhcr(2*npnod,*), hcsr(*)
real*8  :: zanst1, zanst2
real*8  :: prod(*), prodr(*)


!if( lvand ) then
!    call gradresv( nfile, myid, nodes,  &
!& eig, occ, cgjr, rhcr, hcsr, npnod1, npnod2, npnod,  &
!& nband, zanst1, zanst2, prod, prodr )
!  else
    call gradresnv( nfile, myid, nodes,  &
& eig, occ, cgjr, rhcr, npnod1, npnod2, npnod,  &
& nband, zanst1, zanst2, prod, prodr )
!end if


return
end




subroutine gradresnv( nfile, myid, nodes,  &
& eig, occ, cgjr, rhcr, npnod1, npnod2, npnod,  &
& nband, zanst1, zanst2, prod, prodr )
!-----------------------------------------------------------------------
!    gradients and residuals
!-----------------------------------------------------------------------
use engrad_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: eig(*), occ(*)
integer :: npnod1, npnod2, npnod
integer :: nband
real*8  :: cgjr(2*npnod,*), rhcr(2*npnod,*)
real*8  :: zanst1, zanst2
real*8  :: prod(*), prodr(*)

!------declare local variables
logical :: lgammak
integer :: ib, ig
real*8  :: prod1, rhcrig, zansa1, zansa2


call get_lgammak( lgammak )

do ib = 1, nband
   prod1 = 0.0
   do ig = 1, 2*npnod
      rhcrig = eig(ib)*cgjr(ig,ib) - rhcr(ig,ib)
      prod1 = prod1 + rhcrig*rhcrig
   end do
   if( lgammak .and. .not.lnoncollinear ) then
       if( myid.eq.0 ) then
           rhcrig = eig(ib)*cgjr(1,ib) - rhcr(1,ib)
       else
           rhcrig = 0.d0
       end if
       prod1 = 2.d0*prod1 - rhcrig*rhcrig
   end if
   prod(ib) = prod1
end do
call gdsum(prod,nband,prodr)

zansa1 = 0.d0
zansa2 = 0.d0
do ib = 1, nband
   zansa1 = max( zansa1, prod(ib) )
   zansa2 = zansa2 + prod(ib)*occ(ib)
end do

zanst1 = max( zanst1, zansa1 )
zanst2 = zanst2 + zansa2


return
end




subroutine subgs3( nfile, myid, nodes,  &
& ib, ibt, ibb, cgjr, dgjr, nplwex, nplw, nplw7ex,  &
& nband, prod, prodr, ifnorm,  &
& lvand, lvand_r, lvand_g, rvol, rdelv, ntype, nhk1, nhk2, natom,  &
& lvandi, lking, iatoit, ycos, ysin, nplwcs, eigr,  &
& mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh )
!-----------------------------------------------------------------------
!    orthogonalize w.f. to previously determined w.f.'s
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: ib, ibt, ibb
integer :: nplwex, nplw, nplw7ex, nband
#if CCREAL4
real*4  :: cgjr(nplw7ex,nband)
#else
real*8  :: cgjr(nplw7ex,nband)
#endif
real*8  :: dgjr(nplwex)
real*8  :: prod(*), prodr(*)
integer :: ifnorm
logical :: lvand, lvand_r, lvand_g
real*8  :: rvol, rdelv
integer :: ntype, natom
integer :: nhk1(ntype), nhk2(ntype)
logical :: lvandi(ntype), lking(ntype)
integer :: iatoit(natom)
integer :: nplwcs
real*8  :: ycos(*), ysin(*)
real*8  :: eigr(0:*)
integer :: mfd2ft(*), ntotfd, nd1vks(*), kfft0
real*8  :: fft3x(*), fft3y(*)
complex*16 :: fftwork(*)
integer :: ngenh
integer :: ijkgd(-ngenh:ngenh)

!------declare local variables
integer :: jb, ig
real*8  :: prod1, cscr, prod1r


!if( lvand_g ) call calsl2( nfile, myid, nodes,  &
!& dgjr, nplw, nplwex,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking, ycos, ysin, nplwcs )

!if( lvand_r ) call calsl2_r( nfile, myid, nodes,  &
!& eigr, ntype, nhk1, nhk2, natom, lvandi, lking, iatoit, rdelv )


prod(1:ib) = 0.d0
do jb = 1, ib
   do ig = 1, nplw7ex
      prod(jb) = prod(jb) + cgjr(ig,jb)*dgjr(ig)
   end do
end do
do jb = 1, ib
   prod(jb) = 2.d0*prod(jb) - cgjr(1,jb)*dgjr(1)
end do

!if( lvand ) then
!    do jb = 1, ib
!       call calcs3( nfile, myid, nodes,  &
!& prodr(jb), jb, rvol, ntype, nhk1, nhk2, natom, lvandi )
!    end do
!    prod(1:ib) = prod(1:ib) + prodr(1:ib)
!end if

do jb = 1, ib
   dgjr(1:nplw7ex) = dgjr(1:nplw7ex) - prod(jb)*cgjr(1:nplw7ex,jb)
end do


!-----normalization, if needed
if( ifnorm.eq.0 ) return

!if( lvand ) then

!    if( lvand_g ) call calsl2( nfile, myid, nodes,  &
!& dgjr, nplw, nplwex,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking, ycos, ysin, nplwcs )

!    if( lvand_r ) then
!        !-----set wavefunction in real space
!        call wvg2r2( nfile, myid, nodes,  &
!& dgjr, eigr, nplwex, nplw, mfd2ft, ntotfd, nd1vks, kfft0,  &
!& fft3x, fft3y, fftwork, ijkgd, ngenh )
!        call calsl2_r( nfile, myid, nodes,  &
!& eigr, ntype, nhk1, nhk2, natom, lvandi, lking, iatoit, rdelv )
!    end if

!    call setslm( nfile, myid, nodes,  &
!& ibt, ibb, ntype, nhk1, nhk2, natom, lvandi )
!
!end if

prod1 = 0.d0
do ig = 1, nplwex
   prod1 = prod1 + dgjr(ig)*dgjr(ig)
end do
prod1 = 2.d0*prod1 - dgjr(1)*dgjr(1)

!if( lvand ) then
!    call calcsc( nfile, myid, nodes,  &
!& cscr, ibb, ibb, rvol, ntype, nhk1, nhk2, natom, lvandi )
!    prod1 = prod1 + cscr
!end if

prod1r = 1.d0/sqrt(prod1)
dgjr(1:nplwex) = dgjr(1:nplwex)*prod1r
cgjr(1:nplw7ex,ibt) = dgjr(1:nplw7ex)

!if( lvand ) then
!    call cbyslm( nfile, myid, nodes,  &
!& prod1r, ibt, ibb, ntype, nhk1, nhk2, natom, lvandi )
!end if


return
end




subroutine subgs4( nfile, myid, nodes,  &
& ib, gdcr, dgjr, nplwex, nplw, &
& lvand, lvand_r, lvand_g, rvol, rdelv, ntype, nhk1, nhk2, natom,  &
& lvandi, lking, iatoit, ycos, ysin, nplwcs, eigr )
!-----------------------------------------------------------------------
!    orthogonalize w.f. to the current w.f.
!-----------------------------------------------------------------------
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: ib, nplwex, nplw
real*8  :: dgjr(nplwex), gdcr(nplwex)
logical :: lvand, lvand_r, lvand_g
real*8  :: rvol, rdelv
integer :: ntype
integer :: natom
integer, dimension(ntype) :: nhk1, nhk2
logical, dimension(ntype) :: lvandi, lking
integer, dimension(natom) :: iatoit
real*8,  dimension(*) :: ycos, ysin
integer :: nplwcs
real*8, dimension(0:*) :: eigr

!-----declare local variables
real*8  :: hhx, cscr
integer :: ig


!if( lvand_g ) call calsl2( nfile, myid, nodes,  &
!& dgjr, nplw, nplwex,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking, ycos, ysin, nplwcs )

!if( lvand_r ) call calsl2_r( nfile, myid, nodes,  &
!& eigr, ntype, nhk1, nhk2, natom, lvandi, lking, iatoit, rdelv )


hhx = 0.d0
do ig = 1, nplwex
   hhx = hhx + gdcr(ig)*dgjr(ig)
end do
hhx = 2.d0*hhx - gdcr(1)*dgjr(1)
!if( lvand ) then
!    call calcs3( nfile, myid, nodes,  &
!& cscr, ib, rvol, ntype, nhk1, nhk2, natom, lvandi )
!    hhx = hhx + cscr
!end if
do ig = 1, nplwex
   dgjr(ig) = dgjr(ig) - hhx*gdcr(ig)
end do


return
end




subroutine chgr2g( nfile, myid, nodes,  &
& rho, rhgr, nplw5ex, nplw5, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh )
!-----------------------------------------------------------------------
!    charge density transformation from r-space to g-space
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
real*8  :: rho(*), rhgr(*)
integer :: nplw5ex, nplw5, mfd2ft(*), ntotfd, nd1vks(*), kfft0
real*8  :: fft3x(*), fft3y(*)
complex*16 :: fftwork(*)
integer :: ngenh
integer :: ijkgd(-ngenh:ngenh)

!------declare local variables
integer :: nml, inv, kfft1, kfft2, kfft3, ig, ijk, ifd, ift, ierrft


!   --- constants for fft ---
nml = 1
inv = 2
kfft1 = nd1vks(1)
kfft2 = nd1vks(2)
kfft3 = nd1vks(3)

fft3x(1:kfft0) = 0.d0
fft3y(1:kfft0) = 0.d0
do ifd = 1, ntotfd
   ift = mfd2ft(ifd)
   fft3x(ift) = rho(ifd)
end do

call rfft3( nml, fft3x, fft3y, fftwork,  &
&                kfft1, kfft2, kfft3, kfft0, ierrft )

do ig = 0, nplw5
   ijk  = ijkgd(ig)
   rhgr(2*ig+1) = fft3x(ijk)
   rhgr(2*ig+2) = fft3y(ijk)
end do


return
end




subroutine chgg2r( nfile, myid, nodes,  &
& rho, rhgr, nplw5ex, nplw5, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh )
!-----------------------------------------------------------------------
!    charge density transformation from g-space to r-space
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
real*8  :: rho(*), rhgr(*)
integer :: nplw5ex, nplw5, mfd2ft(*), ntotfd, nd1vks(*), kfft0
real*8  :: fft3x(*), fft3y(*)
complex*16 :: fftwork(*)
integer :: ngenh
integer :: ijkgd(-ngenh:ngenh)

!------declare local variables
integer :: nml, inv, kfft1, kfft2, kfft3, ig, ijk, ifd, ift, ierrft


!   --- constants for fft ---
nml = 1
inv = 2
kfft1 = nd1vks(1)
kfft2 = nd1vks(2)
kfft3 = nd1vks(3)

fft3x(1:kfft0) = 0.d0
fft3y(1:kfft0) = 0.d0

do ig = 0, nplw5
   ijk  = ijkgd(ig)
   fft3x(ijk)  = rhgr(2*ig+1)
   fft3y(ijk)  = rhgr(2*ig+2)
end do
do ig = 1, nplw5
   ijk = ijkgd(-ig)
   fft3x(ijk)  =  rhgr(2*ig+1)
   fft3y(ijk)  = -rhgr(2*ig+2)
end do

call rfft3( inv, fft3x, fft3y, fftwork,  &
&                kfft1, kfft2, kfft3, kfft0, ierrft )

!  --- convert FFT -> FD meshes ---
do ifd = 1, ntotfd
   ift = mfd2ft(ifd)
   rho(ifd) = fft3x(ift)
end do


return
end




subroutine chkchg( nfile, myid, nodes, rho, rdelv, nel, ntotfd )
!-----------------------------------------------------------------------
!     check No. of electrons and charge density
!-----------------------------------------------------------------------
use outfile
implicit none
integer :: nfile(*), myid, nodes
real*8  :: rho(*)
real*8  :: rdelv
integer :: nel, ntotfd

!------declare local variables
integer :: ifd
real*8  :: tel, telm


tel = 0.d0
telm = 0.d0
do ifd = 1, ntotfd
   tel = tel + rho(ifd)
   if( rho(ifd) < 0.d0 ) telm = telm + rho(ifd)
end do
tel  = tel  * rdelv
telm = telm * rdelv
if(loutfile(1)) write(nfile(1),*) ' No. of electrons :', tel, telm
if(loutfile(2)) write(nfile(2),*) ' No. of electrons :', tel, telm
tel = dble(nel)/tel
do ifd = 1, ntotfd
   rho(ifd) = rho(ifd)*tel
end do


return
end




subroutine chkchg2( nfile, myid, nodes,  &
& rho, rdelv, nel, ntotfd, tel )
!-----------------------------------------------------------------------
!     check No. of electrons and charge density
!-----------------------------------------------------------------------
use outfile
implicit real*8 ( a-h, o-z )
dimension nfile(*)
dimension rho(*)


tel = 0.d0
do ifd = 1, ntotfd
   tel = tel + rho(ifd)
end do
tel = tel * rdelv
if(loutfile(1)) write(nfile(1),*) ' No. of electrons :',tel
if(loutfile(2)) write(nfile(2),*) ' No. of electrons :',tel
tel = dble(nel)/tel
do ifd = 1, ntotfd
   rho(ifd) = rho(ifd)*tel
end do


return
end




subroutine chkcud( nfile, myid, nodes,  &
& rho, rdelv, diffud, lfixud, ntotfd )
!-----------------------------------------------------------------------
!     check spin charge density
!-----------------------------------------------------------------------
use outfile
implicit real*8 ( a-h, o-z )
dimension nfile(*)
dimension rho(*)
logical   lfixud


tel = 0.d0
do ifd = 1, ntotfd
   tel = tel + rho(ifd)
end do
tel = tel * rdelv
if(loutfile(1)) write(nfile(1),*) ' diff. up & down  :',tel
if(loutfile(2)) write(nfile(2),*) ' diff. up & down  :',tel
if( lfixud .and. abs(tel).gt.1.d-10 ) then
    tel = dble(diffud)/tel
    do ifd = 1, ntotfd
       rho(ifd) = rho(ifd)*tel
    end do
end if


return
end




subroutine chkcud2( nfile, myid, nodes,  &
& rho, rdelv, diffud, lfixud, ntotfd, tel )
!-----------------------------------------------------------------------
!     check spin charge density
!-----------------------------------------------------------------------
use outfile
implicit real*8 ( a-h, o-z )
dimension nfile(*)
dimension rho(*)
logical   lfixud


tel = 0.d0
do ifd = 1, ntotfd
   tel = tel + rho(ifd)
end do
tel = tel * rdelv
if(loutfile(1)) write(nfile(1),*) ' diff. up & down  :',tel
if(loutfile(2)) write(nfile(2),*) ' diff. up & down  :',tel
if( lfixud .and. abs(tel).gt.1.d-10 ) then
    telcor = dble(diffud)/tel
    do ifd = 1, ntotfd
       rho(ifd) = rho(ifd)*telcor
    end do
end if


return
end




subroutine wvg2r( nfile, myid, nodes,  &
& cgjr, nplwex, nplw, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, lflag )
!-----------------------------------------------------------------------
!    wavefunction transformation from g-space to r-space
!
! (input)
!     cgjr   ...... wavefunctions in g space
!
! (output)
!     fft3y  ...... real wavefunctions in r space (on FD  mesh)
!     fft3x  ...... real wavefunctions in r space (on FFT mesh)
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: nplwex, nplw
real*8  :: cgjr(*)
integer :: mfd2ft(*), ntotfd, nd1vks(*), kfft0
real*8  :: fft3x(*), fft3y(*)
complex*16 :: fftwork(*)
integer :: ngenh
integer :: ijkgd(-ngenh:ngenh)
logical :: lflag

!-----declare local variables
integer :: nml, inv, kfft1, kfft2, kfft3, ig, ijk, ifd, ift, ierrft


!   --- constants for fft ---
nml = 1
inv = 2
kfft1 = nd1vks(1)
kfft2 = nd1vks(2)
kfft3 = nd1vks(3)

fft3x(1:kfft0) = 0.d0
fft3y(1:kfft0) = 0.d0

ig = 0
ijk = ijkgd(ig)
fft3x(ijk) = cgjr(2*ig+1)
fft3y(ijk) = 0.d0
do ig = 1, nplw
   ijk  = ijkgd(ig)
   fft3x(ijk)  =  cgjr(2*ig+1)
   fft3y(ijk)  =  cgjr(2*ig+2)
end do
do ig = 1, nplw
   ijk = ijkgd(-ig)
   fft3x(ijk)  =  cgjr(2*ig+1)
   fft3y(ijk)  = -cgjr(2*ig+2)
end do

call rfft3( inv, fft3x, fft3y, fftwork,  &
&                kfft1, kfft2, kfft3, kfft0, ierrft )


if( lflag ) then
    ! --- convert FFT -> FD meshes ---
    do ifd = 1, ntotfd
       ift = mfd2ft(ifd)
       fft3y(ifd) = fft3x(ift)
    end do
end if


return
end




subroutine wvg2r2( nfile, myid, nodes,  &
& cgjr, eigr, nplwex, nplw, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh )
!-----------------------------------------------------------------------
!    wavefunction transformation from g-space to r-space
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: nplwex, nplw
real*8  :: cgjr(*), eigr(0:*)
integer :: mfd2ft(*), ntotfd, nd1vks(*), kfft0
real*8  :: fft3x(*), fft3y(*)
complex*16 :: fftwork(*)
integer :: ngenh
integer :: ijkgd(-ngenh:ngenh)

!-----declare local variables
logical :: lflag


#ifndef VECTOR
lflag = .true.
#else
lflag = .false.
#endif

call wvg2r( nfile, myid, nodes,  &
& cgjr, nplwex, nplw, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, lflag )

eigr(0) = 0.d0
#ifndef VECTOR
eigr(1:ntotfd) = fft3y(1:ntotfd)
#else
eigr(1:kfft0) = fft3x(1:kfft0)
#endif


return
end




subroutine wvr2g( nfile, myid, nodes,  &
& cgjr, nplwex, nplw, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, lflag )
!-----------------------------------------------------------------------
!    wavefunction transformation from g-space to r-space
!
! (input)
!     fft3x  ...... real wavefunctions in r space (on FFT mesh)
!
! (output)
!     cgjr   ...... wavefunctions in g space
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: nplwex, nplw
real*8  :: cgjr(*)
integer :: mfd2ft(*), ntotfd, nd1vks(*), kfft0
real*8  :: fft3x(*), fft3y(*)
complex*16 :: fftwork(*)
integer :: ngenh
integer :: ijkgd(-ngenh:ngenh)
logical :: lflag

!-----declare local variables
integer :: nml, inv, kfft1, kfft2, kfft3, ig, ijk, ifd, ift, ierrft


!   --- constants for fft ---
nml = 1
inv = 2
kfft1 = nd1vks(1)
kfft2 = nd1vks(2)
kfft3 = nd1vks(3)


if( lflag ) then
    ! --- convert FFT -> FD meshes ---
    do ifd = 1, ntotfd
       ift = mfd2ft(ifd)
       fft3x(ift) = fft3y(ifd)
    end do
end if

fft3y(1:kfft0) = 0.d0

call rfft3( nml, fft3x, fft3y, fftwork,  &
&                kfft1, kfft2, kfft3, kfft0, ierrft )

ig = 0
ijk = ijkgd(ig)
cgjr(2*ig+1) = fft3x(ijk)
cgjr(2*ig+2) = 0.d0
do ig = 1, nplw
   ijk  = ijkgd(ig)
   cgjr(2*ig+1) = fft3x(ijk)
   cgjr(2*ig+2) = fft3y(ijk)
end do


return
end




subroutine harting( nfile, myid, nodes,  &
& vhar, rhgr, nplw5ex, nplw5, recnrmex,  &
& mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh )
!-----------------------------------------------------------------------
!      Hartree potential calculation in g space
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension nfile(*)
dimension vhar(*)
dimension rhgr(*)
dimension recnrmex(0:*)
dimension mfd2ft(*), nd1vks(*)
dimension fft3x(*), fft3y(*)
complex*16 fftwork(*)
dimension ijkgd(-ngenh:ngenh)


!   --- constants for fft ---
nml = 1
inv = 2
kfft1 = nd1vks(1)
kfft2 = nd1vks(2)
kfft3 = nd1vks(3)

   do ik=1, kfft0
      fft3x(ik) = 0.d0
      fft3y(ik) = 0.d0
   end do

   ig = 0
   ijk  = ijkgd(ig)
   p4prec = 2.d0*recnrmex(ig)
   fft3x(ijk)  = p4prec*rhgr(2*ig+1)
#ifndef VECTOR
   do ig = 1, nplw5
      ijk  = ijkgd(ig)
      ijkm = ijkgd(-ig)
      p4prec = 2.d0*recnrmex(ig)
      fft3x(ijk)  = p4prec*rhgr(2*ig+1)
      fft3y(ijk)  = p4prec*rhgr(2*ig+2)
      fft3x(ijkm) =  fft3x(ijk)
      fft3y(ijkm) = -fft3y(ijk)
   end do
#else
   do ig = 1, nplw5
      ijk  = ijkgd(ig)
      p4prec = 2.d0*recnrmex(ig)
      fft3x(ijk)  = p4prec*rhgr(2*ig+1)
      fft3y(ijk)  = p4prec*rhgr(2*ig+2)
   end do
   do ig = 1, nplw5
      ijk = ijkgd(-ig)
      p4prec = 2.d0*recnrmex(ig)
      fft3x(ijk)  =  p4prec*rhgr(2*ig+1)
      fft3y(ijk)  = -p4prec*rhgr(2*ig+2)
   end do
#endif

   call rfft3( inv, fft3x, fft3y, fftwork,  &
&                   kfft1, kfft2, kfft3, kfft0, ierrft )

!      --- convert FFT -> FD meshes ---
   do ifd = 1, ntotfd
      ift = mfd2ft(ifd)
      vhar(ifd) = fft3x(ift)
   end do


return
end




