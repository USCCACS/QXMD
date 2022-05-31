



module eigen_variables
!-----------------------------------------------------------------------
! type declaration of shared variables in eigen.f
!-----------------------------------------------------------------------
implicit none

logical :: lnoncollinear
real(8) :: roughfeig, roughfzan

real*8,  dimension(11) :: tc = 0.d0
real*8  :: ct00
real*8,  dimension(11) :: tc2 = 0.d0

real*8  :: av_intcvg, av_ditavg, counteig

logical :: lout_nodecheck = .false.

logical :: lemin = .true.

save

end module




subroutine set_parameters_in_eigen( nfile, myid, nodes,  &
& lnoncollinear_, roughfeig_, roughfzan_ )
!-----------------------------------------------------------------------
! set lnoncollinear
!-----------------------------------------------------------------------
use eigen_variables
implicit none
integer :: nfile(*), myid, nodes
logical :: lnoncollinear_
real(8) :: roughfeig_, roughfzan_

lnoncollinear = lnoncollinear_
roughfeig = roughfeig_
roughfzan = roughfzan_

return
end subroutine




subroutine eigen( nfile, myid, nodes,  &
& eig, gdcr, cgjr, rhcr, hcsr, nplwex, nplw, nplw7ex,  &
& dkgnrm, prcd, ekib,  &
& gnk, hnk, sck, nband, nbnod1, nbnod2, nbnod, iod,  &
& glocal, mfd2ft, ntotfd, nd1vks, kfft0, hfv, hfv_c,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh,  &
& eigr, eigi, apk, apki, scwr, scwi, rdelv, rvol,  &
& ntype, nhk1, nhk2, natom, iatoit, lvand,  &
& lkbpp_r, lvand_r, lkbpp_g, lvand_g, lkbppi, lvandi, lking,  &
& prod, prodr, w1r, w1i, iter, istdot, tolera, threinn, wdinner, lortho,  &
& bzan, nriter, ihyoji, itermx, noband, iteremx,  &
& lhcunt, ltimecnt, leig_start, leig_end,  &
& lmax, mxl, lchk, lclno, nylmmx, ycos, ysin, nplwcs, method,  &
& hfrhcr, jhybrid, lefield_islts )
!-----------------------------------------------------------------------
!   solving eigenvalue problem                       since  '94/08/02
!   by conjugate gradient iteration method
!-----------------------------------------------------------------------
use outfile
use eigen_variables
use ncmagne_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: nplwex, nplw, nplw7ex
real*8  :: eig(*), dkgnrm(*)
real*8  :: gdcr(*), rhcr(*), hcsr(*)
#if CCREAL4
real*4  :: cgjr(*)
#else
real*8  :: cgjr(*)
#endif
real*8  :: prcd(*), ekib(*), gnk(*), hnk(*), sck(*)
integer :: nband, nbnod1, nbnod2, nbnod
integer :: iod(*)
real*8  :: glocal(*), hfv(*), hfv_c(*)
integer :: mfd2ft(*), nd1vks(*), ntotfd, kfft0
real*8  :: fft3x(*), fft3y(*)
complex*16 :: fftwork(*)
integer :: ngenh
integer :: ijkgd(-ngenh:ngenh)
real*8  :: eigr(0:*), eigi(*), apk(*), apki(*), scwr(*), scwi(*)
real*8  :: prod(*), prodr(*), w1r(*), w1i(*)
real*8  :: rdelv, rvol
integer :: nhk1(*), nhk2(*), iatoit(*), ntype, natom
logical :: lvand, lkbpp_r, lvand_r, lkbpp_g, lvand_g
logical :: lkbppi(ntype), lvandi(ntype), lking(ntype)
integer :: iter, istdot
real*8  :: tolera, threinn, wdinner
integer :: ihyoji, itermx, noband, iteremx
logical :: lortho
logical :: lhcunt, ltimecnt, leig_start, leig_end
real*8  :: bzan(*)
integer :: nriter(*)
integer :: lmax(*), lclno(*)
integer :: mxl
logical :: lchk(0:mxl,*)
real*8  :: ycos(*), ysin(*)
integer :: nylmmx, nplwcs
integer :: method  ! 1:line minimization, else:BKL
real*8  :: hfrhcr(*)
integer :: jhybrid
logical :: lefield_islts

!------declare local variables
!integer :: method = 2  ! 1:line minimization, else:BKL
!save method


!ncmagneif: if( .not.lnoncollinear ) then
!
!    hybridif: if( jhybrid == 0 .and. .not.lefield_islts ) then
!
!        if( method == 1 ) then
!            !--- local iterations by line minimization
!            call eigen_lm( nfile, myid, nodes,  &
!& eig, gdcr, cgjr, rhcr, hcsr, nplwex, nplw, nplw7ex,  &
!& dkgnrm, prcd, ekib,  &
!& gnk, hnk, sck, nband, nbnod1, nbnod2, nbnod, iod,  &
!& glocal, mfd2ft, ntotfd, nd1vks, kfft0,  &
!& fft3x, fft3y, fftwork, ijkgd, ngenh, eigr, apk, rdelv, rvol,  &
!& ntype, nhk1, nhk2, natom, iatoit, lvand,  &
!& lkbpp_r, lvand_r, lkbpp_g, lvand_g, lkbppi, lvandi, lking,  &
!& prod, prodr, iter, istdot, tolera, threinn, wdinner, lortho,  &
!& bzan, nriter, ihyoji, itermx, noband, iteremx,  &
!& lhcunt, ltimecnt, leig_start, leig_end,  &
!& lmax, mxl, lchk, lclno, nylmmx, ycos, ysin, nplwcs )
!        else
            !--- local iterations by BKL, PRB 42, 1394 (1990)
            call eigen_BKL( nfile, myid, nodes,  &
& eig, gdcr, cgjr, rhcr, hcsr, nplwex, nplw, nplw7ex,  &
& dkgnrm, prcd, ekib,  &
& gnk, hnk, sck, nband, nbnod1, nbnod2, nbnod, iod,  &
& glocal, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, eigr, apk, rdelv, rvol,  &
& ntype, nhk1, nhk2, natom, iatoit, lvand,  &
& lkbpp_r, lvand_r, lkbpp_g, lvand_g, lkbppi, lvandi, lking,  &
& prod, prodr, iter, istdot, tolera, threinn, wdinner, lortho,  &
& bzan, nriter, ihyoji, itermx, noband, iteremx,  &
& lhcunt, ltimecnt, leig_start, leig_end,  &
& lmax, mxl, lchk, lclno, nylmmx, ycos, ysin, nplwcs )
!        end if

!    else hybridif

        !--- local iterations by line minimization
!        call hfeigen_lm( nfile, myid, nodes,  &
!& eig, gdcr, cgjr, rhcr, hfrhcr, hcsr, nplwex, nplw, nplw7ex,  &
!& dkgnrm, prcd, ekib,  &
!& gnk, hnk, sck, nband, nbnod1, nbnod2, nbnod, iod,  &
!& glocal, mfd2ft, ntotfd, nd1vks, kfft0, hfv, hfv_c,  &
!& fft3x, fft3y, fftwork, ijkgd, ngenh, eigr, apk, rdelv, rvol,  &
!& ntype, nhk1, nhk2, natom, iatoit, lvand,  &
!& lkbpp_r, lvand_r, lkbpp_g, lvand_g, lkbppi, lvandi, lking,  &
!& prod, prodr, iter, istdot, tolera, threinn, wdinner, lortho,  &
!& bzan, nriter, ihyoji, itermx, noband, iteremx,  &
!& lhcunt, ltimecnt, leig_start, leig_end,  &
!& lmax, mxl, lchk, lclno, nylmmx, ycos, ysin, nplwcs )

!    end if hybridif
!
!else ncmagneif
    !-----noncollinear magnetism

    !--- local iterations by line minimization
!    call nceigen_lm( nfile, myid, nodes,  &
!& eig, gdcr, cgjr, rhcr, hcsr, nplwex, nplw, nplw7ex,  &
!& dkdkgnrm, prcd, ekib,  &
!& gnk, hnk, sck, nband, nbnod1, nbnod2, nbnod, iod,  &
!& glocal, mfd2ft, ntotfd, nd1vks, kfft0,  &
!& fft3x, fft3y, fftwork, ncijkg,  &
!& eigr, eigi, eig2r, eig2i, apk, apki, scwr, scwi, rdelv, rvol,  &
!& ntype, nhk1, nhk2, natom, iatoit, lvand,  &
!& lkbpp_r, lvand_r, lkbpp_g, lvand_g, lkbppi, lvandi, lking,  &
!& prod, prodr, w1r, w1i, iter, istdot, tolera, threinn, wdinner, lortho,  &
!& bzan, nriter, ihyoji, itermx, noband, iteremx,  &
!& lhcunt, ltimecnt, leig_start, leig_end,  &
!& lmax, mxl, lchk, lclno, nylmmx, ycos, ysin, nplwcs )

!end if ncmagneif


return
end




subroutine eigen_BKL( nfile, myid, nodes,  &
& eig, gdcr, cgjr, rhcr, hcsr, nplwex, nplw, nplw7ex,  &
& dkgnrm, prcd, ekib,  &
& gnk, hnk, sck, nband, nbnod1, nbnod2, nbnod, iod,  &
& glocal, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, eigr, apk, rdelv, rvol,  &
& ntype, nhk1, nhk2, natom, iatoit, lvand,  &
& lkbpp_r, lvand_r, lkbpp_g, lvand_g, lkbppi, lvandi, lking,  &
& prod, prodr, iter, istdot, tolera, threinn, wdinner, lortho,  &
& bzan, nriter, ihyoji, itermx, noband, iteremx,  &
& lhcunt, ltimecnt, leig_start, leig_end,  &
& lmax, mxl, lchk, lclno, nylmmx, ycos, ysin, nplwcs )
!-----------------------------------------------------------------------
!   solving eigenvalue problem                       since  '94/08/02
!   by conjugate gradient iteration method
!   with BKL, PRB 42, 1394 (1990)
!-----------------------------------------------------------------------
use outfile
use eigen_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: nplwex, nplw, nplw7ex
real*8  :: eig(*)
real*8  :: gdcr(nplwex,*), rhcr(nplwex,*), hcsr(*)
#if CCREAL4
real*4  :: cgjr(nplw7ex,*)
#else
real*8  :: cgjr(nplw7ex,*)
#endif
real*8  :: dkgnrm(*), prcd(*), ekib(*), gnk(*), hnk(*), sck(*)
integer :: nband, nbnod1, nbnod2, nbnod, iod(*)
real*8  :: glocal(*)
integer :: mfd2ft(*), ntotfd, nd1vks(*), kfft0
real*8  :: fft3x(*), fft3y(*)
complex*16 :: fftwork(*)
integer :: ngenh
integer :: ijkgd(-ngenh:ngenh)
real*8  :: eigr(0:*), apk(*)
real*8  :: rdelv, rvol
integer :: ntype, nhk1(*), nhk2(*), natom, iatoit(*)
logical :: lvand
logical :: lkbpp_r, lvand_r, lkbpp_g, lvand_g
logical :: lkbppi(ntype), lvandi(ntype), lking(ntype)
real*8  :: prod(*), prodr(*)
integer :: iter, istdot
real*8  :: tolera, threinn, wdinner
logical :: lortho
real*8  :: bzan(*)
integer :: nriter(*)
integer :: ihyoji, itermx, noband, iteremx
logical :: lhcunt, ltimecnt, leig_start, leig_end
integer :: lmax(*), mxl, lclno(*), nylmmx, nplwcs
logical :: lchk(0:mxl,*)
real*8  :: ycos(*), ysin(*)

!------declare local variables
integer :: i, ib, ibb, jb, ig, nhcstt, izcnt
real*8  :: enclib, eigib
real*8  :: zansa, zansab, z000re, x000re, eigb, dinnerx, zansap, beta, betk
real*8  :: cai, cbi, hhi, cscr, hhx
real*8  :: h11, h12, h21, h22, enew, dtnew
real*8  :: cona0, cona1, cona2, cona3, cona4
real*8  :: sinrt0, sinrt1, sinrt2, sinrt, optm2, dinner
real*8  :: dtex2, alpsq2, alpsq, prod1
logical :: lreset
integer :: ii, intcvg, iitavg, ibuf(2), ibufr(2)
real*8  :: ditavg
real*8  :: ct0, ct, timecnt
save z000re, x000re


if( ltimecnt ) then
    ct0  = timecnt()
    if( leig_start ) then
        tc(:) = 0.d0
        ct00 = ct0
    end if
end if
if( leig_start ) then
    av_intcvg = 0.d0
    av_ditavg = 0.d0
    counteig  = 0.d0
end if

do ibb = 1, nbnod
   ib = iod(ibb)

   !-----set preconditioner
   call prpcg( nfile, myid, nodes,  &
&              ibb, prcd, ekib, dkgnrm, nplwex )

   if( ltimecnt ) then
       ct = timecnt()
       tc(3) = tc(3) + ct - ct0
       ct0 = ct
   end if
!=======================================================================
!     global iteration starts.
!=======================================================================

!--- if lortho=.false., skip the orthogonalization 
if( lortho ) then

   if( ib.ge.2 ) then

!       if( lvand_r ) then
!           !-----set wavefunction in real space
!           call wvg2r2( nfile, myid, nodes,  &
!& gdcr(1,ibb), eigr, nplwex, nplw, mfd2ft, ntotfd, nd1vks, kfft0,  &
!& fft3x, fft3y, fftwork, ijkgd, ngenh )
!       end if

       !---  orthogonalization
       call subgs3( nfile, myid, nodes,  &
& ib-1, ib, ibb, cgjr, gdcr(1,ibb), nplwex, nplw, nplw7ex,  &
& nband, prod, prodr, 1,  &
& lvand, lvand_r, lvand_g, rvol, rdelv, ntype, nhk1, nhk2, natom,  &
& lvandi, lking, iatoit, ycos, ysin, nplwcs, eigr,  &
& mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh )

   end if

   if( ltimecnt ) then
       ct = timecnt()
       tc(1) = tc(1) + ct - ct0
       ct0 = ct
   end if

end if


!--- if lhcunt=.true., skip the HC product calculation
if( .not.lhcunt .or. lortho ) then
!---  calculate H * x --------------------------------------------------
   !--- kinetic energy & HC products
   rhcr(1:nplwex,ibb) = dkgnrm(1:nplwex)*gdcr(1:nplwex,ibb)
   !-----clear hcsr
!   if( lvand ) then
!       ig = nplwex*(ibb-1)
!       hcsr(ig+1:ig+nplwex) = 0.d0
!   end if

   !-----norm-conseving pseudopotentials in KB form
   if( lkbpp_g ) then
       call hckbpp_ib_g( nfile, myid, nodes,  &
& gdcr(1,ibb), rhcr(1,ibb), nplwex, nplw,  &
& ntype, nhk1, nhk2, natom, lkbppi, lking, ycos, ysin, nplwcs )
   end if

   !-----set pointer for hcsr
!   if( lvand ) then
!       nhcstt = nplwex*(ibb-1) + 1
!     else
       nhcstt = 1
!   end if

   !-----ultrasoft pseudopotentials
!   if( lvand_g ) then
!       call hcvand_ib_g( nfile, myid, nodes,  &
!& gdcr(1,ibb), rhcr(1,ibb), hcsr(nhcstt), nplwex, nplw,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking, rvol,  &
!& ycos, ysin, nplwcs )
!   end if

   !---  local term & nonlocal term in real space
   call hcnonloc( nfile, myid, nodes,  &
& enclib, glocal, gdcr(1,ibb), rhcr(1,ibb), hcsr(nhcstt), nplwex, nplw,  &
& nband, nbnod1, nbnod2, nbnod, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, eigr, apk,  &
& ntype, nhk1, nhk2, natom, iatoit,  &
& lkbpp_r, lvand_r, lkbppi, lvandi, lking, rdelv, 0,  &
& lmax, mxl, lchk, lclno, nylmmx, tc(8), .false., .true. )

   if( ltimecnt ) then
       ct = timecnt()
       tc(2) = tc(2) + ct - ct0
       ct0 = ct
   end if
!-----------------------------------------------------------------------
end if

!---  expectation value : eig = x*H*x ----------------------------------
   eigib = 0.d0
   do ig = 1, nplwex
      eigib = eigib + gdcr(ig,ibb)*rhcr(ig,ibb)
   end do
   eigib = 2.d0*eigib - gdcr(1,ibb)*rhcr(1,ibb)
   eig(ib) = eigib
!-----------------------------------------------------------------------

!---  gradient : gnk = -( H - eig )*x = eig*x - rk ---------------------
!   if( lvand ) then
!       ig = nplwex*(ibb-1)
!!                gnk(ig) = eig(ib)*( gdcr(ig,ibb) + hcsr(ig,ibb) )
!       gnk(1:nplwex) = eig(ib)*( gdcr(1:nplwex,ibb)  &
!&                    + hcsr(ig+1:ig+nplwex) ) - rhcr(1:nplwex,ibb)
!     else
       gnk(1:nplwex) = eig(ib)*gdcr(1:nplwex,ibb) - rhcr(1:nplwex,ibb)
!   end if
!---  zansa = < gnk*gnk >  ---------------------------------------------
   zansa = 0.d0
   do ig = 1, nplwex
      zansa = zansa + gnk(ig)*gnk(ig)
   end do
   zansa = 2.d0*zansa - gnk(1)*gnk(1)
!-----------------------------------------------------------------------
   if( ltimecnt ) then
       ct = timecnt()
       tc(3) = tc(3) + ct - ct0
       ct0 = ct
   end if


iter = 0
zansab = zansa
z000re = zansa
x000re = 1.d+10
eigb   = eig(ib)
izcnt  = 0
lreset = .true.
dinnerx= 0.d0
iterdo: do
!-----------------------------------------------------------------------
!     local iteration starts.
!-----------------------------------------------------------------------
iter = iter + 1

!--- preconditioning
   gnk(1:nplwex) = prcd(1:nplwex)*gnk(1:nplwex)
   zansap = 0.d0
   do ig = 1, nplwex
      zansap = zansap + gnk(ig)*gnk(ig)
   end do
   zansap = 2.d0*zansap - gnk(1)*gnk(1)
!-----------------------------------------------------------------------
   if( ltimecnt ) then
       ct = timecnt()
       tc(4) = tc(4) + ct - ct0
       ct0 = ct
   end if

!--- conjugate gradient direction --------------------------------------
   if( lreset ) then
       beta = 0.d0
       hnk(1:nplwex) = 0.d0
     else
       beta = zansap/betk
   end if
   lreset = .false.
   hnk(1:nplwex) = gnk(1:nplwex) + beta*hnk(1:nplwex)
   betk = zansap
   if( ltimecnt ) then
       ct = timecnt()
       tc(3) = tc(3) + ct - ct0
       ct0 = ct
   end if


!   if( lvand_r ) then
!       !-----set wavefunction in real space
!       call wvg2r2( nfile, myid, nodes,  &
!& hnk, eigr, nplwex, nplw, mfd2ft, ntotfd, nd1vks, kfft0,  &
!& fft3x, fft3y, fftwork, ijkgd, ngenh )
!   end if

   !-----orthogonalization
   call subgs3( nfile, myid, nodes,  &
& ib-1, ib, ibb, cgjr, hnk, nplwex, nplw, nplw7ex,  &
& nband, prod, prodr, 0,  &
& lvand, lvand_r, lvand_g, rvol, rdelv, ntype, nhk1, nhk2, natom,  &
& lvandi, lking, iatoit, ycos, ysin, nplwcs, eigr,  &
& mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh )

!   if( lvand_r ) then
!   if( ib > 1 ) then
!       !-----set wavefunction in real space
!       call wvg2r2( nfile, myid, nodes,  &
!& hnk, eigr, nplwex, nplw, mfd2ft, ntotfd, nd1vks, kfft0,  &
!& fft3x, fft3y, fftwork, ijkgd, ngenh )
!   end if
!   end if

   !-----orthogonalization
   call subgs4( nfile, myid, nodes,  &
& ib, gdcr(1,ibb), hnk, nplwex, nplw,  &
& lvand, lvand_r, lvand_g, rvol, rdelv, ntype, nhk1, nhk2, natom,  &
& lvandi, lking, iatoit, ycos, ysin, nplwcs, eigr )

   if( ltimecnt ) then
       ct = timecnt()
       tc(5) = tc(5) + ct - ct0
       ct0 = ct
   end if
!-----------------------------------------------------------------------

!--- new wave function -------------------------------------------------

!---  cai = x*H*h = h*H*x
!---  cbi = h*H*h
!---  hhi = h*S*h
!---  hhx = x*S*h
   cai = 0.d0
   do ig = 1, nplwex
      cai = cai + hnk(ig)*rhcr(ig,ibb)
   end do
   cai = 2.d0*cai - hnk(1)*rhcr(1,ibb)
!   cai = cai * 2.d0

   if( ltimecnt ) then
       ct = timecnt()
       tc(3) = tc(3) + ct - ct0
       ct0 = ct
   end if

!---  gnk = H * h
   gnk(1:nplwex) = dkgnrm(1:nplwex)*hnk(1:nplwex)
!   if( lvand ) then
!       sck(1:nplwex) = 0.d0
!   end if

   !-----norm-conseving pseudopotentials in KB form
   if( lkbpp_g ) then
       call hckbpp_ib_g( nfile, myid, nodes,  &
& hnk, gnk, nplwex, nplw,  &
& ntype, nhk1, nhk2, natom, lkbppi, lking, ycos, ysin, nplwcs )
   end if

   !-----ultrasoft pseudopotentials
!   if( lvand_g ) then
!       call hcvand_ib_g( nfile, myid, nodes,  &
!& hnk, gnk, sck, nplwex, nplw,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking, rvol,  &
!& ycos, ysin, nplwcs )
!   end if

   if( lkbpp_g .or. lvand_g ) then
   if( ltimecnt ) then
       ct = timecnt()
       tc(11) = tc(11) + ct - ct0
       ct0 = ct
   end if
   end if

   !---  local term & nonlocal term in real space
   call hcnonloc( nfile, myid, nodes,  &
& enclib, glocal, hnk, gnk, sck, nplwex, nplw,  &
& nband, nbnod1, nbnod2, nbnod, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, eigr, apk,  &
& ntype, nhk1, nhk2, natom, iatoit,  &
& lkbpp_r, lvand_r, lkbppi, lvandi, lking, rdelv, 0,  &
& lmax, mxl, lchk, lclno, nylmmx, tc(8), ltimecnt, .true. )

   if( ltimecnt ) then
       ct = timecnt()
       tc(6) = tc(6) + ct - ct0
       ct0 = ct
   end if

   cbi = 0.d0
   do ig = 1, nplwex
      cbi = cbi + hnk(ig)*gnk(ig)
   end do
   cbi = 2.d0*cbi - hnk(1)*gnk(1)

   hhi = 0.d0
   do ig = 1, nplwex
      hhi = hhi + hnk(ig)*hnk(ig)
   end do
   hhi = 2.d0*hhi - hnk(1)*hnk(1)
!   if( lvand ) then
!       call calcs4( nfile, myid, nodes,  &
!& cscr, rvol, ntype, nhk1, nhk2, natom, lvandi )
!       hhi = hhi + cscr
!   end if

   ifhhi: if( hhi > 0.d0 ) then
   !=== BKL method =====================================================

   hhx = 0.d0
!   do ig = 1, nplwex
!      hhx = hhx + hnk(ig)*gdcr(ig,ibb)
!   end do
!   hhx = 2.d0*hhx - hnk(1)*gdcr(1,ibb)
!   if( lvand ) then
!       call calcs3( nfile, myid, nodes,  &
!& cscr, ib, rvol, ntype, nhk1, nhk2, natom, lvandi )
!       hhx = hhx + cscr
!   end if
!   !check
!   if(loutfile(1)) write(nfile(1),*) 'hhx:', hhx


   h11 = eig(ib)
   h12 = cai / sqrt(hhi)
   h21 = cai / sqrt(hhi)
   h22 = cbi / hhi
   enew= 0.5d0*(h11+h22-sqrt((h11+h22)**2-4.d0*(h11*h22-h12*h21)))
!   alpsq = h12/sqrt(h12*h12 + (h11-enew)**2))
   dtnew = (enew-h11)/h12 / sqrt(hhi)
   cona3 = hhi
   cona4 = 2.d0*hhx
   alpsq2 = 1.d0/( 1.d0 + dtnew*(cona4 + dtnew*cona3) )
   alpsq  = sqrt( alpsq2 )

!------ check direction of new w.f.
!------ dinner = inner product between old and new w.f.'s
!------ if dinner < threinn, w.f. is not updated.
!if( dtnew > 0.d0 ) then
!    dinner = alpsq*(1.d0 + dtnew*hhx)
!    if( dinner < threinn .or. dinner < dinnerx*wdinner ) then
!    !------ if dinner < dinnerx*wdinner, w.f. is not updated.
!        alpsq = 1.d0
!        dtnew = 0.d0
!        enew  = eig(ib)
!    end if
!    dinnerx = max( dinnerx, dinner )
!end if
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!--- new eigenvalue ----------------------------------------------------
   eig(ib) = enew

   else ifhhi
   !=== line minimization ==============================================

   hhx = 0.d0
   do ig = 1, nplwex
      hhx = hhx + hnk(ig)*gdcr(ig,ibb)
   end do
   hhx = 2.d0*hhx - hnk(1)*gdcr(1,ibb)
!   if( lvand ) then
!       call calcs3( nfile, myid, nodes,  &
!& cscr, ib, rvol, ntype, nhk1, nhk2, natom, lvandi )
!       hhx = hhx + cscr
!   end if


!c---  optimized time step : dtnew
cona0 = eig(ib)
cona1 = cai
cona2 = cbi
cona3 = hhi
cona4 = 2.d0*hhx
if( lemin ) then
    !-----Energy minimization
    sinrt0 = cona0*cona4 - cona1
    sinrt1 = cona0*cona3 - cona2
    sinrt2 = cona1*cona3 - cona2*cona4
    sinrt  = sinrt1*sinrt1 - sinrt2*sinrt0
    if( sinrt > 0.d0 ) then
!        optm1 = ( -sinrt1 + sqrt( sinrt ) ) / sinrt2
        optm2 = ( -sinrt1 - sqrt( sinrt ) ) / sinrt2
        dtnew = optm2
    else
        dtnew = -1.d0
!        write(*,*) 'no answer !?'
    end if
  else
    !-----Residual minimization
!    if( lvand ) then
!        nhcstt = nplwex*(ibb-1) + 1
!      else
        nhcstt = 1
        hcsr(1:nplwex) = 0.d0
         sck(1:nplwex) = 0.d0
!        call fstop( nfile, myid, nodes, 'not supported yet: RMM for NCPP' )
!    end if
    call rmtopt( nfile, myid, nodes, &
& cona0, cona1, cona2, cona3, cona4, rhcr(1,ibb), gdcr(1,ibb), hcsr(nhcstt), &
& gnk, hnk, sck, nplwex, dtnew, .true. )
end if

!------ check direction of new w.f.
!------ dinner = inner product between old and new w.f.'s
!------ if dinner < threinn, w.f. is not updated.
if( dtnew > 0.d0 ) then
    dinner = sqrt( 1.d0/(1.d0+dtnew*dtnew*hhi) )
    if( dinner < threinn ) dtnew = 0.d0
    !------ if dinner < dinnerx*wdinner, w.f. is not updated.
    if( dinner < dinnerx*wdinner ) dtnew = 0.d0
    dinnerx = max( dinnerx, dinner )
end if
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!--- new eigenvalue ----------------------------------------------------
   dtex2 = dtnew*dtnew
   alpsq2 = 1.d0/( 1.d0 + dtnew*cona4 + dtex2*cona3 )
   eig(ib) = ( cona0 + dtnew*cona1 + dtex2*cona2)*alpsq2

   alpsq  = sqrt( alpsq2 )

   !====================================================================
   end if ifhhi

!--- new eigv = eigv + dt*hnk ------------------------------------------
   gdcr(1:nplwex,ibb) = alpsq*( gdcr(1:nplwex,ibb) + dtnew*hnk(1:nplwex) )
   cgjr(1:nplw7ex,ib)  = gdcr(1:nplw7ex,ibb)

!   if( lvand ) then
!
!       call updtslm( nfile, myid, nodes,  &
!& ib, ibb, ntype, nhk1, nhk2, natom, lvandi, alpsq, dtnew )
!
!   end if

!--- new H * x product -------------------------------------------------
   rhcr(1:nplwex,ibb) = alpsq*( rhcr(1:nplwex,ibb) + dtnew*gnk(1:nplwex) )
   !---check eigenenergy
!   eigib = 0.d0
!   do ig = 1, nplwex
!      eigib = eigib + gdcr(ig,ibb)*rhcr(ig,ibb)
!   end do
!   eigib = 2.d0*eigib - gdcr(1,ibb)*rhcr(1,ibb)
!   if(loutfile(1)) write(nfile(1),*) 'eig:', eig(ib), eigib
!   eig(ib) = eigib
!--- new S * x product -------------------------------------------------
!   if( lvand ) then
!       ig = nplwex*(ibb-1)
!!                hcsr(ig,ibb) = alpsq*( hcsr(ig,ibb) + dtnew*sck(ig) )
!       hcsr(ig+1:ig+nplwex) = alpsq*( hcsr(ig+1:ig+nplwex) + dtnew*sck(1:nplwex) )
!   end if
!--- new gradient ------------------------------------------------------
!   if( lvand ) then
!       ig = nplwex*(ibb-1)
!!                gnk(ig) = eig(ib)*( gdcr(ig,ibb) + hcsr(ig,ibb) )
!       gnk(1:nplwex) = eig(ib)*( gdcr(1:nplwex,ibb)  &
!&                    + hcsr(ig+1:ig+nplwex) ) - rhcr(1:nplwex,ibb)
!   else
       gnk(1:nplwex) = eig(ib)*gdcr(1:nplwex,ibb) - rhcr(1:nplwex,ibb)
!   end if
!--- new zansa ---------------------------------------------------------
   zansa = 0.d0
   do ig = 1, nplwex
      zansa = zansa + gnk(ig)*gnk(ig)
   end do
   zansa = 2.d0*zansa - gnk(1)*gnk(1)
!-----------------------------------------------------------------------

if( mod(iter,ihyoji).eq.0 ) then
    if(loutfile(1)) write(nfile(1),'(1x,i5,a2,i4,a1,f15.10,4es12.4)')  &
&           ib,' (',iter,')',eig(ib), eig(ib) - eigb, zansa, dtnew !, dinner
    if(loutfile(2)) write(nfile(2),'(1x,i5,a2,i4,a1,f15.10,4es12.4)')  &
&           ib,' (',iter,')',eig(ib), eig(ib) - eigb, zansa, dtnew !, dinner
end if

   if( ltimecnt ) then
       ct = timecnt()
       tc(3) = tc(3) + ct - ct0
       ct0 = ct
   end if

!--- convergence criterion ---------------------------------------------
   if(  &
&      ib.le.noband .and. iter  .ge. itermx      .or.  &
&      ib.gt.noband .and. iter  .ge. iteremx     .or.  &
&      zansa .lt. tolera                         .or.  &
&      abs(eigb - eig(ib)) .lt. tolera           .or.  &
!&      zansa .lt. z000re*3.d-01                  .or.  &
!&      iter.gt.1 .and. abs(eigb - eig(ib)).lt.x000re*3.d-01  &
&      zansa .lt. z000re*roughfzan               .or.  &
&      iter.gt.1 .and. abs(eigb - eig(ib)).lt.x000re*roughfeig  &
&                                                       ) exit iterdo
   if( iter.eq.1 ) x000re = abs(eigb - eig(ib))

   if( zansa.gt.zansab ) then
       izcnt = izcnt + 1
       if( izcnt.gt.20 ) then
           izcnt = 0
           lreset = .true.
       end if
   end if
   zansab = zansa
   eigb   = eig(ib)
end do iterdo
bzan(ibb)   = zansa
nriter(ibb) = iter
!--------------------------------------------------------
if( lortho ) then
!check ---------------------------------------
!cc      write(*,*) ' check ortho-normalization'
do jb = 1, ib
   prod1 = 0.0
   do ig = 1, nplw7ex
      prod1 = prod1 + cgjr(ig,ib)*cgjr(ig,jb)
   end do
   prod1 = 2.d0*prod1 - cgjr(1,ib)*cgjr(1,jb)
!   if( lvand ) then
!       call calcsc_dup( nfile, myid, nodes,  &
!& cscr, ib, jb, rvol, ntype, nhk1, nhk2, natom, lvandi )
!       prod1 = prod1 + cscr
!   end if
   if( ib.eq.jb ) prod1 = prod1 - 1.d0
   if( abs(prod1).gt.1.d-10 ) then
   if(loutfile(1)) write(nfile(1),*) 'orthonormalization error:',ib,jb,prod1
   if(loutfile(2)) write(nfile(2),*) 'orthonormalization error:',ib,jb,prod1
   end if
end do
!check ---------------------------------------
  else
!check ---------------------------------------
!cc      write(*,*) ' check normalization'
   prod1 = 0.0
   do ig = 1, nplwex
      prod1 = prod1 + gdcr(ig,ibb)*gdcr(ig,ibb)
   end do
   prod1 = 2.d0*prod1 - gdcr(1,ibb)*gdcr(1,ibb)
!   if( lvand ) then
!       call calcsc( nfile, myid, nodes,  &
!& cscr, ibb, ibb, rvol, ntype, nhk1, nhk2, natom, lvandi )
!       prod1 = prod1 + cscr
!   end if
   if( abs(prod1-1.d0).gt.1.d-10 ) then
!cc         if( myid.eq.0 ) then
!       if(loutfile(1))
                       write(nfile(1),*) 'normalization error: iam=',myid,ib,prod1
       if(loutfile(2)) write(nfile(2),*) 'normalization error: iam=',myid,ib,prod1
!cc         end if
   end if
!check ---------------------------------------
end if
if( ltimecnt ) then
    ct = timecnt()
    tc(7) = tc(7) + ct - ct0
    ct0 = ct
end if
!-----------------------------------------------------------------------
!     local iteration end
!-----------------------------------------------------------------------
if( istdot.eq.1 ) then
    if(loutfile(1)) write(nfile(1),'(1x,i5,a2,i4,a1,f20.15,2es13.4)')  &
&       ib,' (',iter,')',eig(ib), eig(ib) - eigb, zansa
    if(loutfile(2)) write(nfile(2),'(1x,i5,a2,i4,a1,f20.15,2es13.4)')  &
&       ib,' (',iter,')',eig(ib), eig(ib) - eigb, zansa
end if
!=======================================================================
!     end of loop for band
!=======================================================================
end do
if( istdot.eq.2 ) then
    do ibb = 1, (nbnod-1)/4 + 1
       do ii = 1, 2
       if( loutfile(ii) ) then
          write(nfile(ii),'(3x,4(i4,a1,i3,e11.4))')  &
&         ( iod(ib), ':',nriter(ib), bzan(ib),  &
&                ib = 1+(ibb-1)*4, min( nbnod, ibb*4 ) )
       end if
       end do
    end do
end if
if( istdot.eq.3 ) then
    intcvg = 0
    iitavg = 0
    do ib = 1, nbnod
       if( bzan(ib).gt.tolera ) intcvg = intcvg + 1
       iitavg = iitavg + nriter(ib)
    end do
    ibuf(1) = intcvg
    ibuf(2) = iitavg
    call gisum(ibuf,2,ibufr)
    intcvg = ibuf(1)
    ditavg = ibuf(2)
    ditavg = ditavg/dble(nband)

    av_intcvg = av_intcvg + intcvg
    av_ditavg = av_ditavg + ditavg
    counteig  = counteig  + 1.d0

    if( leig_end ) then
        av_intcvg = av_intcvg/counteig
        av_ditavg = av_ditavg/counteig
        do ii = 1, 2
        if( loutfile(ii) ) then
!                 write(nfile(ii),'(1x,a33,f7.1)') 
!     &                  ' *** average No. of local it.   :', ditavg
!                 write(nfile(ii),'(1x,a33,i6)')
!     &                  ' *** No. of not converged bands :', intcvg
           write(nfile(ii),'(1x,a30,f6.1,i5)')  &
&     '  # of it. / not conv. bands :', av_ditavg, nint(av_intcvg)
        end if
        end do
    end if
end if


if( leig_end ) then
if( ltimecnt ) then
    ct  = timecnt()
    do i = 1, 2
    if( loutfile(i) ) then
        write(nfile(i),*) '         orthogonalization for w.f. ',  &
&                     ' : cpu-time :', tc(1)
        write(nfile(i),*) '               HC products for w.f. ',  &
&                     ' :          :', tc(2)
        write(nfile(i),*) '                    preconditioning ',  &
&                     ' :          :', tc(4)
        write(nfile(i),*) '         orthogonalization for c.g. ',  &
&                     ' :          :', tc(5)
        write(nfile(i),*) '               HC products for c.g. ',  &
&                     ' :    wvg2r :', tc(8)
        write(nfile(i),*) '               HC products for c.g. ',  &
&                     ' :   calnlc :', tc(9)
        write(nfile(i),*) '               HC products for c.g. ',  &
&                     ' :    rfft3 :', tc(10)
        write(nfile(i),*) '               HC products for c.g. ',  &
&                     ' :    total :', tc(6)
        if( lkbpp_g .or. lvand_g ) then
        write(nfile(i),*) '               HC products for c.g. ',  &
&                     ' : nonlocal :', tc(11)
        end if
        write(nfile(i),*) '          ortho-normalization check ',  &
&                     ' :          :', tc(7)
        write(nfile(i),*) '                      new w.f. etc. ',  &
&                     ' :          :', tc(3)
        write(nfile(i),*) '                          sub total ',  &
&                     ' :          :', ct-ct00
    end if
    end do
    ct0 = ct
    call gsync

!#ifdef IFC
    !-----check CPU time of each node
!    if( lout_nodecheck ) call chkdnodes( nfile, myid, nodes, ct-ct00, 'in eigen' )
!#endif

    ct  = timecnt()
    if(loutfile(1)) write(nfile(1),*) '            waiting time for myid=0 ',  &
&                     ' :          :', ct-ct0
    if(loutfile(2)) write(nfile(2),*) '            waiting time for myid=0 ',  &
&                     ' :          :', ct-ct0
end if
end if


return
end




subroutine rmtopt( nfile, myid, nodes,  &
& cona0, cona1, cona2, cona3, cona4, hwf, wf, swf, hgr, gr, sgr, nplwex,  &
& rmopt, lgammak )
!-----------------------------------------------------------------------
!     optimized time step by residual minimisation
!-----------------------------------------------------------------------
use outfile
implicit none
integer :: nfile(*), myid, nodes
real*8  :: cona0, cona1, cona2, cona3, cona4
integer :: nplwex
real*8,  dimension(nplwex) :: hwf, wf, swf, hgr, gr, sgr
real*8  :: rmopt
logical :: lgammak

!-----declare local variables
real*8  :: cca0, cca1, cca2, ccb0, ccb1, ccb2, ccc0, ccc1, ccc2
integer :: ig
real*8  :: rmoptx, rmmtol, dlambd, alpsq2, rmhh, rmhs, rmss, eilam, &
& resid2, hoptn2, drmhh, drmhs, drmss, dalps, ddlam, deilam, drrdl, &
& rmopt1, drrdl1, rmopt2, drrdl2, dzero
logical :: lrmopt
integer :: irmopt


dzero = 1.d-05

cca0  = 0.d0
cca1  = 0.d0
cca2  = 0.d0
ccb0  = 0.d0
ccb1  = 0.d0
ccb2  = 0.d0
ccc0  = 0.d0
ccc1  = 0.d0
ccc2  = 0.d0
do ig = 1, nplwex
   cca0 = cca0 + hwf(ig)*hwf(ig)
   cca1 = cca1 + hwf(ig)*hgr(ig)
   cca2 = cca2 + hgr(ig)*hgr(ig)
   ccb0 = ccb0 + hwf(ig)*( wf(ig) + swf(ig) )
   ccb1 = ccb1 + hgr(ig)*( wf(ig) + swf(ig) ) + hwf(ig)*( gr(ig) + sgr(ig) )
   ccb2 = ccb2 + hgr(ig)*( gr(ig) + sgr(ig) )
   ccc0 = ccc0 + ( wf(ig) + swf(ig) )*( wf(ig) + swf(ig) )
   ccc1 = ccc1 + ( wf(ig) + swf(ig) )*( gr(ig) + sgr(ig) )
   ccc2 = ccc2 + ( gr(ig) + sgr(ig) )*( gr(ig) + sgr(ig) )
end do
if( lgammak ) then
    ig = 1
    cca0 = 2.d0*cca0 - hwf(ig)*hwf(ig)
    cca1 = 2.d0*cca1 - hwf(ig)*hgr(ig)
    cca2 = 2.d0*cca2 - hgr(ig)*hgr(ig)
    ccb0 = 2.d0*ccb0 - hwf(ig)*( wf(ig) + swf(ig) )
    ccb1 = 2.d0*ccb1 - hgr(ig)*( wf(ig) + swf(ig) ) - hwf(ig)*( gr(ig) + sgr(ig) )
    ccb2 = 2.d0*ccb2 - hgr(ig)*( gr(ig) + sgr(ig) )
    ccc0 = 2.d0*ccc0 - ( wf(ig) + swf(ig) )*( wf(ig) + swf(ig) )
    ccc1 = 2.d0*ccc1 - ( wf(ig) + swf(ig) )*( gr(ig) + sgr(ig) )
    ccc2 = 2.d0*ccc2 - ( gr(ig) + sgr(ig) )*( gr(ig) + sgr(ig) )
end if
cca1  = cca1 * 2.d0
ccb0  = ccb0 * 2.d0
ccb1  = ccb1 * 2.d0
ccb2  = ccb2 * 2.d0
ccc1  = ccc1 * 2.d0


   rmopt = 0.d0
   rmoptx = 1.d0
   rmmtol = 1.d-04
   lrmopt = .false.
   irmopt = 0
do
   irmopt = irmopt + 1
   dlambd = cona0 + ( cona1 + cona2*rmopt )*rmopt
   alpsq2 = 1.d0/( 1.d0 + ( cona4 + cona3*rmopt )*rmopt )
   rmhh   = cca0 + ( cca1 + cca2*rmopt )*rmopt
   rmhs   = ccb0 + ( ccb1 + ccb2*rmopt )*rmopt
   rmss   = ccc0 + ( ccc1 + ccc2*rmopt )*rmopt
   eilam  = dlambd*alpsq2
   resid2 = ( rmhh + ( - rmhs + eilam*rmss )*eilam )*alpsq2
!        write(*,*) dlambd,alpsq2
!        write(*,*) rmhh,rmhs,rmss

   hoptn2 = rmopt*2.d0
   drmhh = cca1 + cca2*hoptn2
   drmhs = ccb1 + ccb2*hoptn2
   drmss = ccc1 + ccc2*hoptn2
   dalps = cona4 + cona3*hoptn2
   ddlam = cona1 + cona2*hoptn2
   deilam = ( -dlambd*dalps*alpsq2 + ddlam )*alpsq2
   drrdl = alpsq2*( -dalps*resid2 &
&              + drmhh - deilam*rmhs - eilam*drmhs &
&              +  ( 2.d0*deilam*rmss + eilam*drmss )*eilam )
!         write(*,*) eilam, resid2
!         return

!         write(*,'(i4,f8.4,10e14.6)') irmopt, rmopt, resid2, drrdl
   if( irmopt > 1000 ) exit
   if( irmopt.eq.1 ) then
       if( drrdl.gt.0.d0 ) rmoptx = -rmoptx
       rmopt1 = rmopt
       drrdl1 = drrdl
       rmopt = rmoptx
       cycle
   else if( lrmopt ) then
       if( drrdl*drrdl1.lt.0.d0 ) then
           rmopt2 = rmopt
           drrdl2 = drrdl
       else if( drrdl*drrdl2.lt.0.d0 ) then
           rmopt1 = rmopt
           drrdl1 = drrdl
       end if
       if( abs(rmopt2-rmopt1).gt.rmmtol ) then
           rmopt = rmopt1 - drrdl1*(rmopt2-rmopt1)/(drrdl2-drrdl1)
           if( rmopt.ge.max(rmopt1, rmopt2) .or. &
&              rmopt.le.min(rmopt1, rmopt2) ) &
&              rmopt = 0.5d0*( rmopt1 + rmopt2 )
           cycle
       end if
       exit
   else
       lrmopt = drrdl*drrdl1.lt.0.d0
       if( lrmopt ) then
           rmopt2 = rmopt
           drrdl2 = drrdl
         else
           rmopt1 = rmopt
           drrdl1 = drrdl
           rmopt  = rmopt + rmoptx
       end if
       cycle
   end if
end do

   if( .not.lrmopt .or. abs(rmopt).lt.dzero ) then
!       if(loutfile(1))
                       write(nfile(1),*) ' *** RMM warning : possibility of ', &
&                 'no optimized time step.'
       if(loutfile(2)) write(nfile(2),*) ' *** RMM warning : possibility of ', &
&                 'no optimized time step.'
       rmopt = 1.d0
   end if


return
end




subroutine prpcg( nfile, myid, nodes,  &
& ib, prcd, ekib, dkgnrm, nplwex )
!-----------------------------------------------------------------------
!    modified version of M.P.Teter, M.C.Payne and D.C.Allen preconditioner
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension nfile(*)
dimension prcd(*), ekib(*), dkgnrm(*)

    do ig = 1, nplwex
       ekib15 = 1.5d0*ekib(ib)
       xcst = dkgnrm(ig)/ekib15
       tesu = 27.d0 + xcst*( 18.d0 + xcst*( 12.d0 + 8.d0*xcst ) )
       xcst2 = xcst*xcst
       xcst4 = xcst2*xcst2
       tesu = tesu/( tesu + 16.d0*xcst4 )
       prcd(ig)  = tesu * 2.d0/ekib15
    end do

return
end




subroutine untryt( nfile, myid, nodes, ct0,  &
& gdcr, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& prod, prodr, prods, betk, w1r, w1i, w2r, ew1, tmpr, iblock,  &
& ltimecnt, leig_start, leig_end,  &
& bijr, biji, dmtrxr, dmtrxi, nbxxxx,  &
& node_r, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, bm1r, bm1i, bx0r )
!-----------------------------------------------------------------------
!   subspace rotation (Unitary transformation)
!-----------------------------------------------------------------------
use outfile
use eigen_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: ct0
integer :: npnod1, npnod2, npnod
integer :: nband, nbnod1, nbnod2, nbnod, nbncnt(*), nbndsp(*)
real*8  :: gdcr(2*npnod,nband)
real*8  :: prod(*), prodr(*), prods(*), betk(*)
real*8  :: w1r(*), w1i(*), w2r(*), ew1(*)
integer :: iblock
real*8  :: tmpr(iblock,*)
logical :: ltimecnt, leig_start, leig_end
integer :: nbxxxx
real*8  :: dmtrxr(nbxxxx,*), dmtrxi(*), bijr(nbxxxx,*), biji(*)
integer :: node_r
integer :: lbncnt(*), lbndsp(*), mbncnt(*), mbndsp(*), jdstnd(*)
real*8  :: bm1r(*), bm1i(*), bx0r(*)


!if( .not.lnoncollinear ) then
    call untryt2( nfile, myid, nodes, ct0,  &
& gdcr, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& prod, prodr, w1r, w1i, ew1, tmpr, iblock,  &
& ltimecnt, leig_start, leig_end,  &
& bijr, dmtrxr, nbxxxx,  &
& node_r, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, bm1r, bx0r )
!else
!    !-----noncollinear magnetism
!    call untryt_k2( nfile, myid, nodes, ct0,  &
!& gdcr, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& prod, prodr, prods, betk, w1r, w1i, w2r, ew1, tmpr, iblock,  &
!& ltimecnt, leig_start, leig_end,  &
!& bijr, biji, dmtrxr, dmtrxi, nbxxxx,  &
!& node_r, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, bm1r, bm1i, bx0r )
!end if


return
end




subroutine untryt2( nfile, myid, nodes, ct0,  &
& gdcr, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& prod, prodr, w1r, w1i, ew1, tmpr, iblock,  &
& ltimecnt, leig_start, leig_end,  &
& bijr, dmtrxr, nbxxxx,  &
& node_r, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, bm1r, bx0r )
!-----------------------------------------------------------------------
!   subspace rotation (Unitary transformation)
!-----------------------------------------------------------------------
use outfile
use eigen_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: ct0
integer :: npnod1, npnod2, npnod
integer :: nband, nbnod1, nbnod2, nbnod, nbncnt(*), nbndsp(*)
real*8  :: gdcr(2*npnod,nband)
real*8  :: prod(*), prodr(*)
real*8  :: w1r(*), w1i(*), ew1(*)
integer :: iblock
real*8  :: tmpr(iblock,*)
logical :: ltimecnt, leig_start, leig_end
integer :: nbxxxx
real*8  :: dmtrxr(nbxxxx,*), bijr(nbxxxx,*), bm1r(*), bx0r(*)
integer :: node_r
integer :: lbncnt(*), lbndsp(*), mbncnt(*), mbndsp(*), jdstnd(*)

!------declare local variables
integer :: ib, jb, jb1, ifvec, ifsrt, i
logical :: lgamma
character*5 :: checkc
real*8  :: epsimr = 1.0d-15
real*8  :: ct, timecnt
save epsimr


if( leig_start ) then
if( ltimecnt ) then
    tc2(1) = 0.d0
    tc2(2) = 0.d0
    tc2(3) = 0.d0
    tc2(4) = 0.d0
end if
end if

do ib = 1, nband
#if PCHOLESKY
   do jb = max(nbnod1, ib), min(nbnod2, nband)
      jb1 = jb - nbnod1 + 1
      bijr(jb1,ib) = dmtrxr(jb1,ib)
   end do
#else
   bijr(1:nband,ib) = dmtrxr(1:nband,ib)
#endif
end do
!----------- eigenvalues and eigenvectors ------------------------------
ifvec = 1
ifsrt = 2
checkc = 'OK   '
!#if PCHOLESKY
!call REIGQR( nfile, myid, nodes,  &
!&            bijr, ew1, nband, nbnod1, nbnod2, nbxxxx,  &
!&            epsimr, ifvec, ifsrt, checkc, w1r, w1i,  &
!&            prod, prodr, nbncnt, nbndsp )
!#else
!if( node_r.eq.1 ) then
    call RsEIGQR( bijr, ew1, nband, nbxxxx,  &
&            epsimr, ifvec, ifsrt, checkc, w1r, w1i, prod )
!  else
!    call vREIGQR( nfile, myid, nodes,  &
!& bijr, ew1, nband, nbxxxx,  &
!& epsimr, ifvec, ifsrt, checkc, w1r, w1i, prod, prodr,  &
!& node_r, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, bx0r, bm1r )
!end if
!#endif
if( checkc.ne.'OK   ' ) then
    if(loutfile(1)) write(nfile(1),*)  &
&           ' warning : error-1 in M(R)EIGQR( in untryt2 )', checkc
    if(loutfile(2)) write(nfile(2),*)  &
&           ' warning : error-1 in M(R)EIGQR( in untryt2 )', checkc
    call fstop( nfile, myid, nodes, 'stop in untryt2' )
end if
!-----------------------------------------------------------------------
if( ltimecnt ) then
    ct = timecnt()
    tc2(1) = tc2(1) + ct - ct0
    ct0 = ct
end if


!---TDDFT-FSSH: check wavefunction exchange
call tddft_fssh_check_unitary( nfile, myid, nodes, &
& bijr, nband, nbnod, nbxxxx, prod, prodr )


call untrycore( nfile, myid, nodes,  &
& gdcr, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& prod, prodr, w1i, tmpr, iblock,  &
& bijr, nbxxxx )


if( ltimecnt ) then
    ct = timecnt()
    tc2(2) = tc2(2) + ct - ct0
    ct0 = ct
end if

if( leig_end ) then
if( ltimecnt ) then
    call get_lgamma( lgamma )
    do i = 1, 2
    if( loutfile(i) ) then
        write(nfile(i),*) '                              reigqr',  &
&                     ' : cpu-time :', tc2(1)
        write(nfile(i),*) '                      transformation',  &
&                     ' :          :', tc2(2)
      if( .not.lgamma ) then
        write(nfile(i),*) '                              meigqr',  &
&                     ' : cpu-time :', tc2(3)
        write(nfile(i),*) '              complex transformation',  &
&                     ' :          :', tc2(4)
      end if
    end if
    end do
end if
end if

!---  new eigenvalues
!      do ib = 1, nband
!         eigib = 0.d0
!#if PCHOLESKY
!         do jb = 1, nbnod2 - nbnod1 + 1
!            w1i(jb) = bijr(jb,ib)
!         end do
!         call alldgatherv(w1i,nbnod2-nbnod1+1,prod,nbncnt,nbndsp)
!#else
!            do jb = 1, nband
!               prod(jb) = bijr(jb,ib)
!            end do
!#endif
!#if PCHOLESKY
!         do jb = nbnod1, nbnod2
!            jb1 = jb - nbnod1 + 1
!#else
!         do jb = 1, nband
!            jb1 = jb
!#endif
!            dnegr = 0.d0
!            do kb = 1, nband
!               dnegr = dnegr + prod(kb)*dmtrxr(jb1,kb)
!            end do
!            eigib = eigib + dnegr*prod(jb)
!         end do
!#if PCHOLESKY
!      call gdsum(eigib,1,prod1)
!#endif
!      end do


return
end




subroutine untryhc( nfile, myid, nodes, ct0,  &
& gdcr, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& prod, prodr, w1r, w1i, tmpr, iblock,  &
& ltimecnt, leig_start, leig_end,  &
& bijr, biji, nbxxxx )
!-----------------------------------------------------------------------
!   subspace rotation (Unitary transformation)
!-----------------------------------------------------------------------
use outfile
use eigen_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: ct0
integer :: npnod1, npnod2, npnod
integer :: nband, nbnod1, nbnod2, nbnod, nbncnt(*), nbndsp(*)
real*8  :: gdcr(2*npnod,nband)
real*8  :: prod(*), prodr(*)
real*8  :: w1r(*), w1i(*)
integer :: iblock
real*8  :: tmpr(iblock,*)
logical :: ltimecnt, leig_start, leig_end
integer :: nbxxxx
real*8  :: bijr(nbxxxx,*), biji(*)


!if( .not.lnoncollinear ) then
    call untryhc2( nfile, myid, nodes, ct0,  &
& gdcr, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& prod, prodr, w1i, tmpr, iblock,  &
& ltimecnt, leig_start, leig_end,  &
& bijr, nbxxxx )
!else
!    !-----noncollinear magnetism
!    call untryhc_k2( nfile, myid, nodes, ct0,  &
!& gdcr, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& prod, prodr, w1r, w1i, tmpr, iblock,  &
!& ltimecnt, leig_start, leig_end,  &
!& bijr, biji, nbxxxx )
!end if


return
end




subroutine untryhc2( nfile, myid, nodes, ct0,  &
& gdcr, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& prod, prodr, w1i, tmpr, iblock,  &
& ltimecnt, leig_start, leig_end,  &
& bijr, nbxxxx )
!-----------------------------------------------------------------------
!   subspace rotation (Unitary transformation)
!-----------------------------------------------------------------------
use outfile
use eigen_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: ct0
integer :: npnod1, npnod2, npnod
integer :: nband, nbnod1, nbnod2, nbnod, nbncnt(*), nbndsp(*)
real*8  :: gdcr(2*npnod,nband)
real*8  :: prod(*), prodr(*)
real*8  :: w1i(*)
integer :: iblock
real*8  :: tmpr(iblock,*)
logical :: ltimecnt, leig_start, leig_end
integer :: nbxxxx
real*8  :: bijr(nbxxxx,*)

!------declare local variables
integer :: i
logical :: lgamma
real*8  :: ct, timecnt


if( leig_start ) then
if( ltimecnt ) then
    tc2(5) = 0.d0
    tc2(6) = 0.d0
end if
end if


call untrycore( nfile, myid, nodes,  &
& gdcr, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& prod, prodr, w1i, tmpr, iblock,  &
& bijr, nbxxxx )


if( ltimecnt ) then
    ct = timecnt()
    tc2(5) = tc2(5) + ct - ct0
    ct0 = ct
end if

if( leig_end ) then
if( ltimecnt ) then
    call get_lgamma( lgamma )
    do i = 1, 2
    if( loutfile(i) ) then
        write(nfile(i),*) '                transformation of HC',  &
&                     ' :          :', tc2(5)
      if( .not.lgamma ) then
        write(nfile(i),*) '        complex transformation of HC',  &
&                     ' :          :', tc2(6)
      end if
    end if
    end do
end if
end if


return
end




subroutine untrycore( nfile, myid, nodes,  &
& gdcr, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& prod, prodr, w1i, tmpr, iblock,  &
& bijr, nbxxxx )
!-----------------------------------------------------------------------
!   subspace rotation (Unitary transformation)
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: npnod1, npnod2, npnod
integer :: nband, nbnod1, nbnod2, nbnod, nbncnt(*), nbndsp(*)
real*8  :: gdcr(2*npnod,nband)
real*8  :: prod(*), prodr(*)
real*8  :: w1i(*)
integer :: iblock
real*8  :: tmpr(iblock,*)
integer :: nbxxxx
real*8  :: bijr(nbxxxx,*)

!------declare local variables
integer :: ib, jb
integer :: nmaxx, ncblck, icblck, icb1, icb2
real*8  :: prodjb


!--- transform HC products ---
nmaxx = 2*npnod
call gimax(nmaxx)
ncblck = nmaxx/iblock
do icblck = 1, ncblck + 1
   icb1 = (icblck-1)*iblock + 1
   icb2 =  icblck   *iblock
   icb2 = min( icb2, 2*npnod )

   do ib = 1, nband
#if PCHOLESKY
      do jb = 1, nbnod2 - nbnod1 + 1
         w1i(jb) = bijr(jb,ib)
      end do
      call alldgatherv(w1i,nbnod2-nbnod1+1,prod,nbncnt,nbndsp)
#else
      prod(1:nband) = bijr(1:nband,ib)
#endif

      tmpr(1:iblock,ib) = 0.d0
!OCL UNROLL(9)
      do jb = 1, nband
         prodjb = prod(jb)
         tmpr(1:icb2-icb1+1,ib) = tmpr(1:icb2-icb1+1,ib) + prodjb*gdcr(icb1:icb2,jb)
      end do
   end do
!OCL UNROLL(9)
   gdcr(icb1:icb2,1:nband) = tmpr(1:icb2-icb1+1,1:nband)

end do


return
end




