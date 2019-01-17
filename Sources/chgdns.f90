



subroutine pwpchg( nfile, myid, nodes, ltimecnt, lfssh_updtcg )
!-----------------------------------------------------------------------
!    charge density
!-----------------------------------------------------------------------
use param
implicit none
integer :: nfile(*), myid, nodes
logical :: ltimecnt, lfssh_updtcg


!if( .not.lnoncollinear ) then
    call pwpchg2( nfile, myid, nodes, ltimecnt, lfssh_updtcg )
!else
!    !-----noncollinear magnetism
!    call ncpwpchg2( nfile, myid, nodes, ltimecnt, lfssh_updtcg )
!end if


return
end




subroutine pwpchg2( nfile, myid, nodes, ltimecnt, lfssh_updtcg )
!-----------------------------------------------------------------------
!    charge density
!-----------------------------------------------------------------------
use param
use param_atom
use pwlda_pp
use pwlda_atom
use pwlda_variables
use pwlda_proc
use pwlda_pw
use pwlda_grid
implicit none
integer :: nfile(*), myid, nodes
logical :: ltimecnt, lfssh_updtcg

!------declare local variables
integer :: nspin, kvec, ifd
real*8  :: densup, densdn
real*8  :: csh0, timecnt
logical :: leig_end


if( ltimecnt ) then
    csh0 = timecnt()
end if
do nspin = 1, nspnmx
   if( lspin ) then
       call ldocck( nspin, occ, wegud, nband, nkpnt )
       if( lgamma ) then
           call stspud( nspin, gdcr, gdcrsv, nspnod, lcgjsv )
!           if( lvand ) call ldslmi( nfile, myid, nodes, nspin, lvandi, .false. )
       end if
   end if

   dupif3: if( nspin == 1 .or. nspin == 2 .and. .not.lduplicate .or. lfssh_updtcg ) then

   if( .not.lgamma ) then
       glocal(:) = 0.d0
   end if

call set_all2all_out( .false. )
kvecdo7: do kvec = 1, nknod
!   if( .not.lgamma ) then
!       call set_kvec( kvec )
!       call get_wfk( gdcr, nspin )
!       call get_plwk( nplw, nplwex, nspnod )
!       call get_indexk(  &
!& npnod1, npnod2, npnod, nplcnt, npldsp, nodes )
!!             if( lvand ) then
!!                 call ldslmi_k( nfile, myid, nodes, lvandi, nspin, .false. )
!!             end if
!   end if

   leig_end   = kvec == nknod

   !   --- copy gdcr to rhcr
   call cpgdrh( nfile, myid, nodes, rhcr, gdcr, npnod, nband )
   !   --- to convert G decomposition to band decomposition ---
   call set_all2all_out( leig_end )
   call gdtobd( nfile, myid, nodes, csh0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& rhcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iods, iodsg, idstnd, ltimecnt )
   !--- charge density on coarse grid (from soft part) --------------------
   if( lgamma ) then
       call softchg( nfile, myid, nodes,  &
& glocal, occ, rhcr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod,  &
& iods, rvol, mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, ijkg, nplw2,  &
& fft3x, fft3y, fftwork )
   else
!       call softchg_k( nfile, myid, nodes,  &
!& apk, occ, rhcr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod,  &
!& iods, rvol,  &
!& mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, ijkg, nplw2ex, nplw2,  &
!& fft3x, fft3y, fftwork, eigr, eigi )
!       !---symmetry operation
!       call softchg_k_sym( nfile, myid, nodes,  &
!& apk, nplw2ex, nplw2, mfd2ft_c, ntotfd_c, nd1vks_c, kfft0,  &
!& nga, ngb, ngc, nplw5, fft3x, fft3y, fftwork, ijkg )
!       !---add apk to glocal
!       call dadd_b_to_a( glocal, apk, ntotfd )
   end if

end do kvecdo7

   call outchg( nfile, myid, nodes, csh0,  &
& glocal, occ, nband, nbnod1, nbnod2, nbnod,  &
& rvol, mfd2ft, ntotfd, nd1vks, kfft0d, ijkgd, nplw5,  &
& mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, ijkg, nplw2ex, nplw2,  &
& fft3x, fft3y, fftwork, rdelv, nel, rhgr, nplw5ex, nplw5,  &
& nplw3, ngboxa, ngboxb, ngboxc,  &
& kfft1b, kfft2b, kfft3b, kfft0b, ijkgb,  &
& eigr, lvand, ntype, natom, iatoit, ioa, lvandi,  &
& lspin, nspin, ltimecnt )

!   if( .not.lgamma ) then
!
!       !---unify glocal across k-point decomposition domain
!       call unify_sumn( glocal, ntotfd, apk )
!
!       if( .not.lspin ) then
!       if( myid.eq.0 ) then
!           if( ltimecnt ) then
!               !--- check No. of electrons and charge density
!               call chkchg( nfile, myid, nodes, glocal, rdelv, nel, ntotfd )
!           end if
!           if( lvand ) then
!               !--- charge density transformation from r- to g-spaces
!               call chgr2g( nfile, myid, nodes,  &
!& glocal, rhgr, nplw5ex, nplw5, mfd2ft, ntotfd, nd1vks, kfft0d,  &
!& fft3x, fft3y, fftwork, ijkgd, nplw5 )
!           else
!               call chgr2g( nfile, myid, nodes,  &
!& glocal, rhgr, nplw2ex, nplw2, mfd2ft_c, ntotfd_c, nd1vks_c, kfft0,  &
!& fft3x, fft3y, fftwork, ijkg, nplw2 )
!           end if
!       end if
!       end if
!   end if

   if( lgamma ) then
   if( .not.lspin ) then
   if( ltimecnt ) then
   if( myid.eq.0 ) then
       !--- check No. of electrons and charge density
       call chkchg( nfile, myid, nodes, glocal, rdelv, nel, ntotfd )
   end if
   end if
   end if
   end if

   if( lspin ) then
       if( myid == 0 ) then
           if( .not.lduplicate .or. lfssh_updtcg ) then
               call exrhcr( glocal, hlocal, ntotfd )
           else
               call dcopy_a_to_b( glocal, hlocal, ntotfd )
           end if
       end if
   end if

   end if dupif3
end do


if( lspin ) then
if( myid == 0 ) then
!--- glocal :   up-spin electron density -> total density
!--- hlocal : down-spin electron density -> spin  density
    do ifd = 1, ntotfd
       densup = glocal(ifd)
       densdn = hlocal(ifd)
       glocal(ifd) = densup + densdn
       hlocal(ifd) = densup - densdn
    end do

    !--- check charge density and spin density
!    call checkchg( nfile, myid, nodes,  &
!& glocal, hlocal, mfd2ft, ntotfd, nd1vks, kfft0d, fft3x, fft3y )

!    if( lvand ) then
!        !--- charge density transformation from r- to g-spaces
!        call chgr2g( nfile, myid, nodes,  &
!& glocal, rhgr, nplw5ex, nplw5, mfd2ft, ntotfd, nd1vks, kfft0d,  &
!& fft3x, fft3y, fftwork, ijkgd, nplw5 )
!      else
    !--- convert local-order grid (mfd2ft_c)
    !---     to global-order grid (mfd2ft)
        !--- charge density transformation from r- to g-spaces
        call chgr2g( nfile, myid, nodes,  &
& hlocal, rhgr, nplw2ex, nplw2, mfd2ft_c, ntotfd_c, nd1vks_c, kfft0,  &
& fft3x, fft3y, fftwork, ijkg, nplw2 )
        !--- charge density transformation from g- to r-spaces
        call chgg2r( nfile, myid, nodes,  &
& hlocal, rhgr, nplw5ex, nplw5, mfd2ft, ntotfd, nd1vks, kfft0d,  &
& fft3x, fft3y, fftwork, ijkgd, nplw5 )
        !--- charge density transformation from r- to g-spaces
        call chgr2g( nfile, myid, nodes,  &
& glocal, rhgr, nplw2ex, nplw2, mfd2ft_c, ntotfd_c, nd1vks_c, kfft0,  &
& fft3x, fft3y, fftwork, ijkg, nplw2 )
!    end if

    if( ltimecnt ) then
        !--- check No. of electrons and charge density
        call chkchg( nfile, myid, nodes,  &
& glocal, rdelv, nel, ntotfd )
!        call chkcud( nfile, myid, nodes,  &
!& hlocal, rdelv, diffud, lfixud, ntotfd )
    end if
end if
end if
!--- store local charge density
!      call distlc( nfile, myid, nodes,
!     & tmpn, mshnod(1), glocal, mftdsp )
if( lspin ) then
!          call distlc( nfile, myid, nodes,
!     & tmpl, mshnod(1), hlocal, mftdsp )
    call dscatterv( hlocal, mftnod, mftdsp, tmpl, mshnod, 0 )
end if
!--- all to all communication ------------------------------------------
!--- to obtain charge density in spatial decomposition
!      call atachg( nfile, myid, nodes, csh0,
!     & tmpn, glocal, mftnod, mftdsp, mfd2ft, ntotfd, dbuf,
!     & idstnd, ltimecnt )
!          call chknel( nfile, myid, nodes,                                   &
!     & tmpn, mshnod(1), rdelv, nel, noddatx, ltimecnt )
if( lspin ) then
    call chknud( nfile, myid, nodes,                                   &
& tmpl, mshnod(1), rdelv, diffud, lfixud, noddatx, ltimecnt )
end if


return
end




subroutine outchg( nfile, myid, nodes, ct0,  &
& glocal, occ, nband, nbnod1, nbnod2, nbnod,  &
& rvol, mfd2ft, ntotfd, nd1vks, kfft0d, ijkgd, ngenh,  &
& mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, ijkg, nplw2ex, nplw2,  &
& fft3x, fft3y, fftwork, rdelv, nel, rhgr, nplw5ex, nplw5,  &
& nplw3, ngboxa, ngboxb, ngboxc,  &
& kfft1b, kfft2b, kfft3b, kfft0b, ijkgb,  &
& eigr, lvand, ntype, natom, iatoit, ioa, lvandi,  &
& lspin, nspin, ltimecnt )
!-----------------------------------------------------------------------
!    output charge density in each node
!-----------------------------------------------------------------------
use outfile
implicit none
integer :: nfile(*), myid, nodes
real*8  :: ct0
integer :: nband, nbnod1, nbnod2, nbnod
real*8  :: glocal(*)
real*8  :: occ(*)
real*8  :: rvol
integer :: mfd2ft(*), ntotfd, nd1vks(*), kfft0d
integer :: ngenh
integer :: ijkgd(-ngenh:ngenh)
integer :: mfd2ft_c(*), ntotfd_c, nd1vks_c(*), kfft0, nplw2ex, nplw2
integer :: ijkg(-nplw2:nplw2)
real*8  :: fft3x(*), fft3y(*)
complex*16 :: fftwork(*)
real*8  :: rdelv
integer :: nel
real*8  :: rhgr(*)
integer :: nplw5ex, nplw5
integer :: nplw3
integer :: ngboxa(0:nplw3), ngboxb(0:nplw3), ngboxc(0:nplw3)
integer :: kfft1b, kfft2b, kfft3b, kfft0b
integer :: ijkgb(-nplw3:nplw3)
real*8  :: eigr(0:*)
logical :: lvand
integer :: ntype, natom, iatoit(natom), ioa(*)
logical :: lvandi(ntype)
logical :: lspin
integer :: nspin
logical :: ltimecnt

!------declare local variables
integer :: ifd
real*8  :: tel, telm
real*8  :: ct, timecnt


rhgr(1:nplw5ex) = 0.d0


if( ltimecnt ) then
    tel  = 0.d0
    telm = 0.d0
    do ifd = 1, ntotfd_c
       tel = tel + glocal(ifd)
       if( glocal(ifd) < 0.d0 ) telm = telm + glocal(ifd)
    end do
    tel = tel * rdelv * dble(kfft0d)/dble(kfft0)
    call unify_sum1( tel )
    telm = telm * rdelv * dble(kfft0d)/dble(kfft0)
    call unify_sum1( telm )
    if(loutfile(1)) write(nfile(1),*) ' No. of el. from soft part :', tel, telm
    if(loutfile(2)) write(nfile(2),*) ' No. of el. from soft part :', tel, telm
    call gsync
    ct = timecnt()
    if(loutfile(1)) write(nfile(1),*) '               c.d. from soft part ',  &
&                         ' : cpu-time :', ct-ct0
    if(loutfile(2)) write(nfile(2),*) '               c.d. from soft part ',  &
&                         ' : cpu-time :', ct-ct0
    ct0 = ct
end if


if( .not.lspin .or. lspin.and.lvand ) then
    !--- charge density transformation from r- to g-spaces
    if( myid.eq.0 ) then
        call chgr2g( nfile, myid, nodes,  &
& glocal, rhgr, nplw2ex, nplw2, mfd2ft_c, ntotfd_c, nd1vks_c, kfft0,  &
& fft3x, fft3y, fftwork, ijkg, nplw2 )
    end if
end if


!if( lvand ) then
!
!    !-----c.d. on dense grid (from hard part, i.e. Q functions)
!    call chgfrq( nfile, myid, nodes,  &
!& glocal, nplw3, ngboxa, ngboxb, ngboxc,  &
!& kfft1b, kfft2b, kfft3b, kfft0b, fft3x, fft3y, fftwork, ijkgb,  &
!& nband, occ, eigr, ntotfd, kfft0d, rvol,  &
!& ntype, natom, iatoit, ioa, lvandi, nspin )
!
!    if( ltimecnt ) then
!        tel = 0.d0
!        telm = 0.d0
!        do ifd = 1, ntotfd
!           tel = tel + glocal(ifd)
!           if( glocal(ifd) < 0.d0 ) telm = telm + glocal(ifd)
!        end do
!        tel = tel * rdelv
!        call unify_sum1( tel )
!        telm = telm * rdelv
!        call unify_sum1( telm )
!        if(loutfile(1)) write(nfile(1),*) ' No. of el. from hard part :', tel, telm
!        if(loutfile(2)) write(nfile(2),*) ' No. of el. from hard part :', tel, telm
!    end if
!
!    if( myid.eq.0 ) then
!        !--- hard part of c.d. transformation from r- to g-spaces
!        call chgr2g_sum( nfile, myid, nodes,  &
!& glocal, rhgr, nplw5ex, nplw5, mfd2ft, ntotfd, nd1vks, kfft0d,  &
!& fft3x, fft3y, fftwork, ijkgd, ngenh )
!
!        !--- charge density transformation from g- to r-spaces
!        call chgg2r( nfile, myid, nodes,  &
!& glocal, rhgr, nplw5ex, nplw5, mfd2ft, ntotfd, nd1vks, kfft0d,  &
!& fft3x, fft3y, fftwork, ijkgd, ngenh )
!    end if
!
!    if( ltimecnt ) then
!        call gsync
!        ct = timecnt()
!        if(loutfile(1)) write(nfile(1),*) '               c.d. from hard part ',  &
!&                         ' : cpu-time :', ct-ct0
!        if(loutfile(2)) write(nfile(2),*) '               c.d. from hard part ',  &
!&                         ' : cpu-time :', ct-ct0
!        ct0 = ct
!    end if
!
!end if


return
end




subroutine softchg( nfile, myid, nodes,  &
& glocal, occ, cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod,  &
& iod, rvol, mfd2ft, ntotfd, nd1vks, kfft0, ijkgd, ngenh,  &
& fft3x, fft3y, fftwork )
!-----------------------------------------------------------------------
!    charge density from soft part of pp.
!
!    output: glocal ...... charge density only on node=0
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: nplwex, nplw, nband, nbnod1, nbnod2, nbnod, iod(*)
real*8  :: glocal(*)
real*8  :: occ(*)
real*8  :: cgjr(nplwex,*)
real*8  :: rvol
integer :: mfd2ft(*), ntotfd, nd1vks(*), kfft0
integer :: ngenh
integer :: ijkgd(-ngenh:ngenh)
real*8  :: fft3x(*), fft3y(*)
complex*16 :: fftwork(*)

!------declare local variables
integer :: ifd, ift, i, ib
real*8  :: wrvolr


#ifndef VECTOR
glocal(1:ntotfd) = 0.d0
#else
fft3y(1:kfft0) = 0.d0
#endif

do i = 1, nbnod
   ib = iod(i)
   wrvolr = occ(ib)/rvol

#ifndef VECTOR
   call wvg2r( nfile, myid, nodes,  &
& cgjr(1,i), nplwex, nplw, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, .true. )

   glocal(1:ntotfd) = glocal(1:ntotfd) + fft3y(1:ntotfd)*fft3y(1:ntotfd)*wrvolr
#else
   call wvg2r( nfile, myid, nodes,  &
& cgjr(1,i), nplwex, nplw, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, glocal, fftwork, ijkgd, ngenh, .false. )

   fft3y(1:kfft0) = fft3y(1:kfft0) + fft3x(1:kfft0)*fft3x(1:kfft0)*wrvolr
#endif

end do

#ifdef VECTOR
! --- convert FFT -> FD meshes ---
do ifd = 1, ntotfd
   ift = mfd2ft(ifd)
   glocal(ifd) = fft3y(ift)
end do
#endif

!--- global sum of charge density in r-space
call dsum(glocal,ntotfd,fft3x,0)


return
end




subroutine pwpinichg( nfile, myid, nodes )
!-----------------------------------------------------------------------
! return charge density to input charge density
!-----------------------------------------------------------------------
use param
use param_atom
use pwlda_pp
use pwlda_atom
use pwlda_variables
use pwlda_proc
use pwlda_pw
use pwlda_grid
implicit none
integer :: nfile(*), myid, nodes

!------declare local variables


!--- return rhgr to input charge density
rhgr(1:nplw5ex) = rinr(1:nplw5ex)

call chgg2r( nfile, myid, nodes,  &
& glocal, rhgr, nplw5ex, nplw5, mfd2ft, ntotfd, nd1vks, kfft0d,  &
& fft3x, fft3y, fftwork, ijkgd, nplw5 )

!--- store local charge density
call distlc( nfile, myid, nodes, rho, mshnod(1), glocal, mftdsp )


!if( lnoncollinear ) then
!    !-----noncollinear magnetism
!    call pwpinimagne( nfile, myid, nodes, mshnod(1) )
!end if


return
end




module mxchgg_variables
!-----------------------------------------------------------------------
! type declaration of variables for charge density mixing
!-----------------------------------------------------------------------
implicit none

integer :: nplwcex, imxchg, itratn
real*8,  allocatable, dimension(:,:) :: unr, dfmr

real*8,  allocatable, dimension(:,:) :: aij, betakl, wrk
real*8,  allocatable, dimension(:)   :: cmk
integer, allocatable, dimension(:)   :: ip

real*8,  allocatable, dimension(:,:) :: dmu, dfi
real*8,  allocatable, dimension(:)   :: wn

integer :: lplchex = 0    ! nplwcex
integer :: kntchg  = 0

real*8  :: factq1
real*8  :: g2min, g2max

save

end module




subroutine set_g2ming2max_in_chgdns( g2min_, g2max_ )
!-----------------------------------------------------------------------
use mxchgg_variables
implicit none
real*8  :: g2min_, g2max_

g2min = g2min_
g2max = g2max_

return
end




subroutine imxchgg_change( nfile, myid, nodes, imxchg_ )
!-----------------------------------------------------------------------
use mxchgg_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: imxchg_

imxchg  = imxchg_

return
end




subroutine mxchgg_alloc( nfile, myid, nodes, &
& alloc_mem, lspin, nplwcex_, imxchg_, itratn_, imxspin_, nmxspin_,  &
& cmetric, spinmetric, wkerker, pwscale, lnoncollinear )
!-----------------------------------------------------------------------
!     allocate memory for charge density mixing
!-----------------------------------------------------------------------
use mxchgg_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
real*8  :: alloc_mem
logical :: lspin
integer :: nplwcex_, imxchg_, itratn_, imxspin_, nmxspin_
real*8  :: cmetric, spinmetric, wkerker, pwscale
logical :: lnoncollinear

!-----declare local variables
integer :: status, ierror
real*8  :: small, bbb
real*8  :: the_mem
logical :: licall = .true.
save licall


if( licall ) then

    imxchg  = imxchg_
    itratn  = itratn_

    if( imxchg > 0 ) then
        !------allocate memory
        allocate( aij(itratn,itratn), betakl(itratn,itratn),  &
   & wrk(itratn,itratn), cmk(itratn), ip(itratn),  &
   & stat=status )

        the_mem = 8.d0 * ( size(aij) + size(betakl) + size(wrk) + size(cmk) ) &
   & + 4.d0 * size(ip)

        if( imxchg == 5 .or. imxchg == 6 ) then
            allocate( wn(itratn), stat=status )
            the_mem = the_mem + 8.d0 * ( size(wn) )
        end if

        !------error trap
        call check_alloc( nfile, myid, nodes,  &
   & status, alloc_mem, the_mem, 'mxchgg_alloc', .true. )
    end if

    if( imxchg == 5 .or. imxchg == 6 ) then
        wn(1:itratn) = 1.d0
    end if

end if


nplwcex = nplwcex_
if( imxchg > 0 .and. nplwcex > lplchex  ) then

    !-----if already allocated, deallocate arrays
    if( allocated(unr) ) then

        the_mem = 8.d0 * ( size(unr) + size(dfmr) )

        !-----deallocate arrays
        deallocate( unr, dfmr, stat=status )

        if( allocated(dmu) ) then
            the_mem = the_mem + 8.d0 * ( size(dmu) + size(dfi) )
            !-----deallocate arrays
            deallocate( dmu, dfi, stat=status )
        end if

        !------error trap
        call check_dealloc( nfile, myid, nodes,  &
     & status, alloc_mem, the_mem, 'mxchgg_alloc(2)', .true. )

    end if


    !-----allocate arrays
    lplchex    = nplwcex* pwscale

    allocate( unr(lplchex,itratn), dfmr(lplchex,itratn),  &
   & stat=status )

    the_mem = 8.d0 * ( size(unr) + size(dfmr) )

    if( imxchg == 4 .or. imxchg == 5 .or. imxchg == 6 ) then
        allocate( dmu(lplchex,2), dfi(lplchex,2),  &
   & stat=status )

        the_mem = the_mem + 8.d0 * ( size(dmu) + size(dfi) )
    end if

    !------error trap
    call check_alloc( nfile, myid, nodes,  &
   & status, alloc_mem, the_mem, 'mxchgg_alloc(2)', .true. )

end if


!--- factq1 : parameter for the metric in charge density mixing
small=1.0d-5
bbb = g2max - cmetric*g2min
if( bbb.gt.small ) then
    factq1 = (cmetric - 1.d0)*g2max*g2min/bbb
  else
    factq1 = 0.d0
end if


if( lspin ) call mxchggud_alloc( nfile, myid, nodes, &
& alloc_mem, lspin, nplwcex_, imxspin_, nmxspin_, spinmetric, wkerker,  &
& g2min, g2max, pwscale )

!if( lnoncollinear ) call ncmxchg_alloc( nfile, myid, nodes, &
!& alloc_mem, imxspin_, nmxspin_, spinmetric, wkerker, nplwcex_,  &
!& g2min, g2max, pwscale )

licall = .false.


return
end subroutine




subroutine pwpmxchgg( nfile, myid, nodes, ct0,  &
& lreset, aslh_, bslh_, amxspin_, bmxspin_, lmixhold, ltimecnt, lresetud )
!-----------------------------------------------------------------------
!     mixing charge density
!-----------------------------------------------------------------------
use param
implicit none
integer :: nfile(*), myid, nodes
real*8  :: ct0
logical :: lreset, lmixhold, ltimecnt
real*8  :: aslh_, bslh_, amxspin_, bmxspin_
logical :: lresetud


!if( .not.lnoncollinear ) then
    call pwpmxchgg2( nfile, myid, nodes, ct0,  &
& lreset, lresetud, aslh_, bslh_, amxspin_, bmxspin_, lmixhold, ltimecnt )
!else
!    !-----noncollinear magnetism
!    call ncpwpmxchgg2( nfile, myid, nodes, ct0,  &
!& lreset, lresetud, aslh_, bslh_, amxspin_, bmxspin_, lmixhold, ltimecnt )
!end if


return
end subroutine




subroutine pwpmxchgg2( nfile, myid, nodes, ct0,  &
& lreset, lresetud, aslh_, bslh_, amxspin_, bmxspin_, lmixhold, ltimecnt )
!-----------------------------------------------------------------------
!     mixing charge density
!-----------------------------------------------------------------------
use outfile
use param
use param_atom
use pwlda_pp
use pwlda_atom
use pwlda_variables
use pwlda_proc
use pwlda_pw
use pwlda_grid
implicit none
integer :: nfile(*), myid, nodes
real*8  :: ct0
logical :: lreset, lmixhold, ltimecnt
real*8  :: aslh_, bslh_, amxspin_, bmxspin_
logical :: lresetud

!------declare local variables
integer :: ii
real*8  :: ct, timecnt


call mxchgg( nfile, myid, nodes,  &
& lreset, rhgr, rinr, nplw5ex, nplw5, dkgnrm, rvol, aslh_, bslh_,  &
& glocal, apk, eigr, mfd2ft, ntotfd, nd1vks, kfft0d,  &
& fft3x, fft3y, fftwork, ijkgd, nplw5, rdelv, nel, lmixhold )
!--- store local charge density
!          call distlc( nfile, myid, nodes,
!     & rho, mshnod(1), glocal, mftdsp )
call dscatterv( glocal, mftnod, mftdsp, rho, mshnod, 0 )

if( lspin ) then
    rhoud(1:mshnod(1)) = tmpl(1:mshnod(1))
!              call unifylc( nfile, myid, nodes,
!     & rhoud, mshnod(1), hlocal,
!     & mftnod, mftdsp, mfd2ft, ntotfd, nd1vks )
    call dgatherv( rhoud, mshnod, hlocal, mftnod, mftdsp, 0)
!--- charge density transformation from r- to g-spaces
    if( myid == 0 ) then
        call chgr2g( nfile, myid, nodes,  &
& hlocal, rhgr(nplw5ex+1), nplw5ex, nplw5, mfd2ft, ntotfd, nd1vks, kfft0d,  &
& fft3x, fft3y, fftwork, ijkgd, nplw5 )
        !-----Check symmetry of charge density
!        call csymcd( nfile, myid, nodes, ltimecnt,  &
!& rhgr(nplw5ex+1), nd1vks, nplw5, nplw5ex, kfft0d, nga, ngb, ngc, ijkgd, 2 )
    end if
    call mxchggud( nfile, myid, nodes,  &
& lresetud, rhgr(nplw5ex+1), rinr(nplw5ex+1), nplw5ex, nplw5, dkgnrm, &
& rvol, amxspin_, bmxspin_,  &
& hlocal, mfd2ft, ntotfd, nd1vks, kfft0d,  &
& fft3x, fft3y, fftwork, ijkgd, nplw5, rdelv, nel, lmixhold )
    call dscatterv( hlocal, mftnod, mftdsp, rhoud, mshnod, 0 )
end if

ct = timecnt()
if( ltimecnt ) then
    do ii = 1, 2
       if( loutfile(ii) ) then
           write(nfile(ii),*) ' Charge mixing                 ',  &
&                             ' : cpu-time :', ct-ct0
!              write(nfile(ii),*) '                               ',
!     &                           ' : com-time :', t_comm
       end if
    end do
end if
ct0 = ct


return
end subroutine




subroutine mxchgg( nfile, myid, nodes,  &
& lreset, rhgr, rinr, nplw5ex, nplw5, dkgnrm, rvol, aslh, bslh,  &
& glocal, apk, eigr, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, rdelv, nel, lmixhold )
!-----------------------------------------------------------------------
!     mixing charge density in reciprocal space
!
!  (input)
!       rinr   ...... input charge density
!       rhgr   ...... output charge density
!
!  (output)
!       rhgr   ...... charge density for next SCF cycle
!
!-----------------------------------------------------------------------
use mxchgg_variables
implicit none
integer :: nfile(*), myid, nodes
logical :: lreset
integer :: nplw5ex, nplw5
real*8  :: rhgr(nplw5ex), rinr(nplw5ex)
real*8  :: dkgnrm(nplwcex)
real*8  :: rvol, aslh, bslh
integer :: ntotfd, kfft0
real*8  :: glocal(ntotfd), apk(*), eigr(*)
integer :: mfd2ft(ntotfd), nd1vks(3)
real*8  :: fft3x(kfft0), fft3y(kfft0)
complex*16 :: fftwork(kfft0)
integer :: ngenh, nel
integer :: ijkgd(-ngenh:ngenh)
real*8  :: rdelv
logical :: lmixhold   ! .true. = hold matrix in charge mixing

!-----declare local variables
logical :: lvarweit

if( imxchg == 0 ) then

    !--- No mixing
    call mxchgg_none( nfile, myid, nodes,  &
& lreset, rhgr, rinr, nplw5ex, nplw5,  &
& glocal, mfd2ft, ntotfd, nd1vks, kfft0, fft3x, fft3y, fftwork, ijkgd, ngenh )

else if( imxchg == 1 ) then

    !--- Pulay mixing
    call mxchgg_pulay( nfile, myid, nodes,  &
& lreset, rhgr, rinr, nplw5ex, nplw5, dkgnrm, rvol, aslh, bslh, factq1,  &
& glocal, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, rdelv, nel, lmixhold, &
& kntchg, itratn, nplwcex, unr, dfmr, lplchex, aij, betakl, wrk, cmk, ip,  &
& 1.d0 )

!else if( imxchg == 2 ) then

    !--- Anderson mixing
!    call mxchgg_anderson( nfile, myid, nodes,  &
!& lreset, rhgr, rinr, nplw5ex, nplw5, dkgnrm, rvol, aslh, bslh, factq1,  &
!& glocal, mfd2ft, ntotfd, nd1vks, kfft0,  &
!& fft3x, fft3y, fftwork, ijkgd, ngenh, rdelv, nel, lmixhold, &
!& kntchg, itratn, nplwcex, unr, dfmr, lplchex, aij, betakl, wrk, cmk, ip,  &
!& 1.d0 )

else if( imxchg == 3 ) then

    !--- Simple mixing
    call mxchgg_simple( nfile, myid, nodes,  &
& lreset, rhgr, rinr, nplw5ex, nplw5, aslh, &
& glocal, mfd2ft, ntotfd, nd1vks, kfft0, &
& fft3x, fft3y, fftwork, ijkgd, ngenh, rdelv, nel )

!else if( imxchg == 4 ) then

    !--- Broyden's second mixing by Srivastava
!    call mxchgg_srivastava( nfile, myid, nodes,  &
!& lreset, rhgr, rinr, nplw5ex, nplw5, dkgnrm, rvol, aslh, bslh, factq1,  &
!& glocal, mfd2ft, ntotfd, nd1vks, kfft0,  &
!& fft3x, fft3y, fftwork, ijkgd, ngenh, rdelv, nel, lmixhold, &
!& kntchg, itratn, nplwcex, unr, dfmr, lplchex, aij, betakl, wrk, cmk, ip,  &
!& dmu, dfi, 1.d0 )

!else if( imxchg == 5 .or. imxchg == 6 ) then

!    lvarweit = imxchg == 6
!    !--- Broyden's second mixing by Johnson
!    call mxchgg_johnson( nfile, myid, nodes,  &
!& lreset, rhgr, rinr, nplw5ex, nplw5, dkgnrm, rvol, aslh, bslh, factq1,  &
!& glocal, mfd2ft, ntotfd, nd1vks, kfft0,  &
!& fft3x, fft3y, fftwork, ijkgd, ngenh, rdelv, nel, lmixhold, &
!& kntchg, itratn, nplwcex, unr, dfmr, lplchex, aij, betakl, wrk, cmk, ip,  &
!& dmu, dfi, wn, 1.d0, lvarweit )

end if

!---check negative charge density
!call check_negative_rho( nfile, myid, nodes,  &
!& rhgr, rinr, nplw5ex, nplw5,  &
!& glocal, apk, eigr, mfd2ft, ntotfd, nd1vks, kfft0, &
!& fft3x, fft3y, fftwork, ijkgd, ngenh, rdelv, nel )


return
end subroutine




subroutine mxchgg_none( nfile, myid, nodes,  &
& lreset, rhgr, rinr, nplw5ex, nplw5,  &
& glocal, mfd2ft, ntotfd, nd1vks, kfft0, fft3x, fft3y, fftwork, ijkgd, ngenh )
!-----------------------------------------------------------------------
!     Do not mix charge density
!
!  (input)
!       rinr   ...... input charge density
!       rhgr   ...... output charge density
!
!  (output)
!       rhgr   ...... charge density for next SCF cycle
!
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
logical :: lreset
integer :: nplw5ex, nplw5
real*8  :: rhgr(nplw5ex), rinr(nplw5ex)
integer :: ntotfd, kfft0
real*8  :: glocal(ntotfd)
integer :: mfd2ft(ntotfd), nd1vks(3)
real*8  :: fft3x(kfft0), fft3y(kfft0)
complex*16 :: fftwork(kfft0)
integer :: ngenh
integer :: ijkgd(-ngenh:ngenh)


!-----broadcast to unify rhgr
call dbcast(rhgr,nplw5ex,0)

!--- store rinr as input charge density
rinr(1:nplw5ex) = rhgr(1:nplw5ex)

!---to keep consistency with symmetry operations
if( myid == 0 ) then
    !--- charge density transformation from g- to r-spaces
      call chgg2r( nfile, myid, nodes,  &
& glocal, rhgr, nplw5ex, nplw5, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh )
end if

lreset = .false.


return
end subroutine




subroutine mxchgg_pulay( nfile, myid, nodes,  &
& lreset, rhgr, rinr, nplw5ex, nplw5, dkgnrm, rvol, aslh, bslh, factq1,  &
& glocal, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, rdelv, nel, lmixhold, &
& kntchg, itratn, nplwcex, unr, dfmr, lplchex, aij, betakl, wrk, cmk, ip,  &
& wkerker )
!-----------------------------------------------------------------------
!     mixing charge density in reciprocal space
!
!                modified Broyden's method ( PRB38(1988)12807 )
!                Pulay's method ( Chem.Phys.Lett.73(1980)393 ) 
!
!  (input)
!       rinr   ...... input charge density
!       rhgr   ...... output charge density
!
!  (output)
!       rhgr   ...... charge density for next SCF cycle
!
!-----------------------------------------------------------------------
use outfile
implicit none
integer :: nfile(*), myid, nodes
logical :: lreset
integer :: nplw5ex, nplw5
real*8  :: rhgr(nplw5ex), rinr(nplw5ex)
real*8  :: dkgnrm(nplwcex)
real*8  :: rvol, aslh, bslh, factq1
integer :: ntotfd, kfft0
real*8  :: glocal(ntotfd)
integer :: mfd2ft(ntotfd), nd1vks(3)
real*8  :: fft3x(kfft0), fft3y(kfft0)
complex*16 :: fftwork(kfft0)
integer :: ngenh, nel
integer :: ijkgd(-ngenh:ngenh)
real*8  :: rdelv
logical :: lmixhold   ! .true. = hold matrix in charge mixing

integer :: kntchg, itratn, nplwcex, lplchex
real*8  :: unr(lplchex,itratn), dfmr(lplchex,itratn)
real*8  :: aij(itratn,itratn), betakl(itratn,itratn), wrk(itratn,itratn), cmk(itratn)
integer :: ip(itratn)
real*8  :: wkerker

!-----declare local variables
integer :: ig, i1, icntr, ii, jj, ierr, nhi, i, j
real*8  :: dfmsq, factq, fslhig, amax, bbnb
!integer :: icount = 0, ix, iy, iz, ijk, kfft1, kfft2, kfft3
!save icount
real*8  :: small = 1.d-15
save small


node0if: if( myid == 0 ) then

      kntchg = kntchg + 1
!      aslh = 0.8d0
!      bslh = 0.64d0

!--- Pulay mixing ------------------------------------------------------
      if( kntchg.gt.itratn .or. lreset ) kntchg = 1
      resetdo: do
       unr(1:nplwcex,kntchg) = rinr(1:nplwcex)
      dfmr(1:nplwcex,kntchg) = rhgr(1:nplwcex) - rinr(1:nplwcex)
      do i1 = 1, kntchg
         dfmsq = 0.d0
         do ig = 3, nplwcex
            factq = ( dkgnrm(ig) + factq1 )/dkgnrm(ig)
            dfmsq = dfmsq + dfmr(ig,i1)*dfmr(ig,kntchg) * factq
         end do
         ig = 3
         factq = ( dkgnrm(ig) + factq1 )/dkgnrm(ig)
         dfmsq = 2.d0*dfmsq + dfmr(1,i1)*dfmr(1,kntchg) * factq
         dfmsq = dfmsq*rvol
         aij(i1,kntchg) = dfmsq
         aij(kntchg,i1) = dfmsq
      end do

      ierr = 0
      kntchgif: if( kntchg.eq.1 ) then
          !--- trial by Kerker mixing ---
          do ig = 1, nplwcex
             fslhig = dkgnrm(ig)/( dkgnrm(ig) + bslh )
             fslhig = aslh*( 1.d0 - wkerker + wkerker*fslhig )
             rhgr(ig) = rinr(ig) + fslhig*dfmr(ig,kntchg)
          end do
      else kntchgif
          !---
          i1 = kntchg
          icntr = i1/2 + mod(i1,2)
!                  amax = max( abs(aij(icntr,icntr)), 1.d0 )
          amax = abs(aij(icntr,icntr))
          !---error trap
          if( amax < small ) then
              kntchg = 1
              cycle resetdo
          end if
          aij(1:i1,1:i1) = aij(1:i1,1:i1)/amax
          !---
          CALL INVERS( aij, betakl, i1, WRK, IP, itratn, ierr )
          ierrif: if( ierr == 0 ) then
              !---
                 aij(1:i1,1:i1) =    aij(1:i1,1:i1)*amax
              betakl(1:i1,1:i1) = betakl(1:i1,1:i1)/amax
              !---
              bbnb = 0.d0
              do ii = 1, i1
                 cmk(ii) = 0.d0
                 do jj = 1, i1
                    bbnb    = bbnb    + betakl(ii,jj)
                    cmk(ii) = cmk(ii) + betakl(ii,jj)
                 end do
              end do
              cmk(1:i1) = cmk(1:i1) / bbnb
              rhgr(1:nplwcex) = 0.d0
              do ii = 1, i1
                 do ig = 1, nplwcex
                    fslhig = dkgnrm(ig)/( dkgnrm(ig) + bslh )
                    fslhig = aslh*( 1.d0 - wkerker + wkerker*fslhig )
                    rhgr(ig) = rhgr(ig)  &
&                            + cmk(ii)*( unr(ig,ii) + fslhig*dfmr(ig,ii) )
                 end do
              end do
          else ierrif
              !--- ierr /= 0
!              kntchg = 1
              nhi = 1    !--- max( 1, kntchg - 1 )
              kntchg = max( 1, kntchg - nhi )
              if(loutfile(1)) write(nfile(1),*) '*** information(Pulay) : reset kntchg =', kntchg
              if(loutfile(2)) write(nfile(2),*) '*** information(Pulay) : reset kntchg =', kntchg
              do ii = 1, kntchg - 1
                 unr(1:nplwcex,ii) =  unr(1:nplwcex,ii+nhi)
                dfmr(1:nplwcex,ii) = dfmr(1:nplwcex,ii+nhi)
              end do
              do j = 1, kntchg - 1
              do i = 1, kntchg - 1
                 aij(i,j) = aij(i+nhi,j+nhi)*amax
              end do
              end do
          end if ierrif
      end if kntchgif
      if( ierr == 0 ) exit
      end do resetdo

      if( kntchg.eq.itratn ) then
          do i = 1, kntchg - 1
              unr(1:nplwcex,i) =  unr(1:nplwcex,i+1)
             dfmr(1:nplwcex,i) = dfmr(1:nplwcex,i+1)
          end do
          do j = 1, kntchg - 1
          do i = 1, kntchg - 1
             aij(i,j) = aij(i+1,j+1)
          end do
          end do

          if( lmixhold ) kntchg = kntchg - 1
      end if
!-----------------------------------------------------------------------


!--- charge density transformation from g- to r-spaces
      call chgg2r( nfile, myid, nodes,  &
& glocal, rhgr, nplw5ex, nplw5, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh )

!!--- check No. of electrons and charge density
!      call chkchg( nfile, myid, nodes, glocal, rdelv, nel, ntotfd )
!    icount = icount + 1
!    if( mod(icount-1,10) == 0 ) then
!
!        kfft1 = nd1vks(1)
!        kfft2 = nd1vks(2)
!        kfft3 = nd1vks(3)
!
!        iz = kfft3/2
!        iy = kfft2/2
!        do ix = 0, kfft1-1
!           ijk = 1 + ix + kfft1*iy + kfft1*kfft2*iz
!           write(100+icount-1,'(i4, es14.6)') ix, fft3x(ijk)
!        end do
!
!    end if

end if node0if


!-----broadcast to unify rhgr
call dbcast(rhgr,nplw5ex,0)

!--- store rinr as input charge density
rinr(1:nplw5ex) = rhgr(1:nplw5ex)

lreset = .false.


return
end




subroutine mxchgg_simple( nfile, myid, nodes,  &
& lreset, rhgr, rinr, nplw5ex, nplw5, aslh, &
& glocal, mfd2ft, ntotfd, nd1vks, kfft0, &
& fft3x, fft3y, fftwork, ijkgd, ngenh, rdelv, nel )
!-----------------------------------------------------------------------
!     simple mixing charge density in reciprocal space
!
!  (input)
!       rinr   ...... input charge density
!       rhgr   ...... output charge density
!
!  (output)
!       rhgr   ...... charge density for next SCF cycle
!
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
logical :: lreset
integer :: nplw5ex, nplw5
real*8  :: rhgr(nplw5ex), rinr(nplw5ex)
real*8  :: aslh
integer :: ntotfd, kfft0
real*8  :: glocal(ntotfd)
integer :: mfd2ft(ntotfd), nd1vks(3)
real*8  :: fft3x(kfft0), fft3y(kfft0)
complex*16 :: fftwork(kfft0)
integer :: ngenh, nel
integer :: ijkgd(-ngenh:ngenh)
real*8  :: rdelv
!integer :: kntchg

!-----declare local variables
integer :: ig
real*8  :: ph


node0if: if( myid == 0 ) then

!      kntchg = kntchg + 1
!      ph  =  ( ph + phi ) / 2.d0
      ph = aslh
      rhgr(1:nplw5ex) =  ( 1.d0 - ph ) * rinr(1:nplw5ex) + ph * rhgr(1:nplw5ex)

!--- charge density transformation from g- to r-spaces
      call chgg2r( nfile, myid, nodes,  &
& glocal, rhgr, nplw5ex, nplw5, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh )

!!--- check No. of electrons and charge density
!      call chkchg( nfile, myid, nodes, glocal, rdelv, nel, ntotfd )

end if node0if


!-----broadcast to unify rhgr
call dbcast(rhgr,nplw5ex,0)

!--- store rinr as input charge density
rinr(1:nplw5ex) = rhgr(1:nplw5ex)

lreset = .false.


return
end




module mxchggud_variables
!-----------------------------------------------------------------------
! type declaration of variables for charge density mixing
!-----------------------------------------------------------------------
implicit none

integer :: nplwcex, imxchg, itratn
real*8,  allocatable, dimension(:,:) :: unr, dfmr

real*8,  allocatable, dimension(:,:) :: aij, betakl, wrk
real*8,  allocatable, dimension(:)   :: cmk
integer, allocatable, dimension(:)   :: ip

real*8,  allocatable, dimension(:,:) :: dmu, dfi
real*8,  allocatable, dimension(:)   :: wn

integer :: lplchex = 0    ! nplwcex
integer :: kntchg  = 0

real*8  :: factq1
real*8  :: wkerker

save

end module




subroutine imxchggud_change( nfile, myid, nodes, imxchg_ )
!-----------------------------------------------------------------------
use mxchggud_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: imxchg_

imxchg  = imxchg_

return
end




subroutine mxchggud_alloc( nfile, myid, nodes, &
& alloc_mem, lspin, nplwcex_, imxspin_, nmxspin_, spinmetric, wkerker_,  &
& g2min, g2max, pwscale )
!-----------------------------------------------------------------------
!     allocate memory for charge density mixing
!-----------------------------------------------------------------------
use mxchggud_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
real*8  :: alloc_mem
logical :: lspin
integer :: nplwcex_, imxspin_, nmxspin_
real*8  :: spinmetric, wkerker_, pwscale
real*8  :: g2min, g2max

!-----declare local variables
integer :: status, ierror
real*8  :: small, bbb
real*8  :: the_mem
logical :: licall = .true.
save licall


if( licall ) then

    imxchg  = imxspin_
    itratn  = nmxspin_
    wkerker = wkerker_

    if( imxchg > 0 ) then
        !------allocate memory
        allocate( aij(itratn,itratn), betakl(itratn,itratn),  &
   & wrk(itratn,itratn), cmk(itratn), ip(itratn),  &
   & stat=status )

        the_mem = 8.d0 * ( size(aij) + size(betakl) + size(wrk) + size(cmk) ) &
   & + 4.d0 * size(ip)

        if( imxchg == 5 .or. imxchg == 6 ) then
            allocate( wn(itratn), stat=status )
            the_mem = the_mem + 8.d0 * ( size(wn) )
        end if

        !------error trap
        call check_alloc( nfile, myid, nodes,  &
   & status, alloc_mem, the_mem, 'mxchggud_alloc', .true. )
    end if

    if( imxchg == 5 .or. imxchg == 6 ) then
        wn(1:itratn) = 1.d0
    end if

end if


nplwcex = nplwcex_
if( imxchg > 0 .and. nplwcex > lplchex  ) then

    !-----if already allocated, deallocate arrays
    if( allocated(unr) ) then

        the_mem = 8.d0 * ( size(unr) + size(dfmr) )

        !-----deallocate arrays
        deallocate( unr, dfmr, stat=status )

        if( allocated(dmu) ) then
            the_mem = the_mem + 8.d0 * ( size(dmu) + size(dfi) )
            !-----deallocate arrays
            deallocate( dmu, dfi, stat=status )
        end if

        !------error trap
        call check_dealloc( nfile, myid, nodes,  &
     & status, alloc_mem, the_mem, 'mxchggud_alloc(2)', .true. )

    end if


    !-----allocate arrays
    lplchex    = nplwcex* pwscale

    allocate( unr(lplchex,itratn), dfmr(lplchex,itratn),  &
   & stat=status )

    the_mem = 8.d0 * ( size(unr) + size(dfmr) )

    if( imxchg == 4 .or. imxchg == 5 .or. imxchg == 6  ) then
        allocate( dmu(lplchex,2), dfi(lplchex,2),  &
   & stat=status )

        the_mem = the_mem + 8.d0 * ( size(dmu) + size(dfi) )
    end if

    !------error trap
    call check_alloc( nfile, myid, nodes,  &
   & status, alloc_mem, the_mem, 'mxchggud_alloc(2)', .true. )

end if


!--- factq1 : parameter for the metric in spin-density mixing
small=1.0d-5
bbb = g2max - spinmetric*g2min
if( bbb.gt.small ) then
    factq1 = (spinmetric - 1.d0)*g2max*g2min/bbb
  else
    factq1 = 0.d0
end if


licall = .false.


return
end subroutine




subroutine mxchggud( nfile, myid, nodes,  &
& lreset, rhgr, rinr, nplw5ex, nplw5, dkgnrm, rvol, aslh, bslh,  &
& glocal, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, rdelv, nel, lmixhold )
!-----------------------------------------------------------------------
!     mixing spin-charge density in reciprocal space
!
!  (input)
!       rinr   ...... input charge density
!       rhgr   ...... output charge density
!
!  (output)
!       rhgr   ...... charge density for next SCF cycle
!
!-----------------------------------------------------------------------
use mxchggud_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
logical :: lreset
integer :: nplw5ex, nplw5
real*8,  dimension(nplw5ex) :: rhgr, rinr
real*8,  dimension(nplwcex) :: dkgnrm
real*8  :: rvol, aslh, bslh
integer :: ntotfd, kfft0
real*8,  dimension(ntotfd) :: glocal
integer :: mfd2ft(ntotfd), nd1vks(3)
real*8,  dimension(kfft0) :: fft3x, fft3y
complex*16, dimension(kfft0) :: fftwork
integer :: ngenh, nel
integer :: ijkgd(-ngenh:ngenh)
real*8  :: rdelv
logical :: lmixhold   ! .true. = hold matrix in charge mixing

!-----declare local variables
logical :: lvarweit


if( imxchg == 0 ) then

    !--- No mixing
    call mxchgg_none( nfile, myid, nodes,  &
& lreset, rhgr, rinr, nplw5ex, nplw5,  &
& glocal, mfd2ft, ntotfd, nd1vks, kfft0, fft3x, fft3y, fftwork, ijkgd, ngenh )

else if( imxchg == 1 ) then

    !--- Pulay mixing
    call mxchgg_pulay( nfile, myid, nodes,  &
& lreset, rhgr, rinr, nplw5ex, nplw5, dkgnrm, rvol, aslh, bslh, factq1,  &
& glocal, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, rdelv, nel, lmixhold, &
& kntchg, itratn, nplwcex, unr, dfmr, lplchex, aij, betakl, wrk, cmk, ip,  &
& wkerker )

!else if( imxchg == 2 ) then

    !--- Anderson mixing
!    call mxchgg_anderson( nfile, myid, nodes,  &
!& lreset, rhgr, rinr, nplw5ex, nplw5, dkgnrm, rvol, aslh, bslh, factq1,  &
!& glocal, mfd2ft, ntotfd, nd1vks, kfft0,  &
!& fft3x, fft3y, fftwork, ijkgd, ngenh, rdelv, nel, lmixhold, &
!& kntchg, itratn, nplwcex, unr, dfmr, lplchex, aij, betakl, wrk, cmk, ip,  &
!& wkerker )

else if( imxchg == 3 ) then

    !--- Simple mixing
    call mxchgg_simple( nfile, myid, nodes,  &
& lreset, rhgr, rinr, nplw5ex, nplw5, aslh, &
& glocal, mfd2ft, ntotfd, nd1vks, kfft0, &
& fft3x, fft3y, fftwork, ijkgd, ngenh, rdelv, nel )

!else if( imxchg == 4 ) then

    !--- Broyden's second mixing by Srivastava
!    call mxchgg_srivastava( nfile, myid, nodes,  &
!& lreset, rhgr, rinr, nplw5ex, nplw5, dkgnrm, rvol, aslh, bslh, factq1,  &
!& glocal, mfd2ft, ntotfd, nd1vks, kfft0,  &
!& fft3x, fft3y, fftwork, ijkgd, ngenh, rdelv, nel, lmixhold, &
!& kntchg, itratn, nplwcex, unr, dfmr, lplchex, aij, betakl, wrk, cmk, ip,  &
!& dmu, dfi, wkerker )

!else if( imxchg == 5 .or. imxchg == 6 ) then

!    lvarweit = imxchg == 6
!    !--- Broyden's second mixing by Johnson
!    call mxchgg_johnson( nfile, myid, nodes,  &
!& lreset, rhgr, rinr, nplw5ex, nplw5, dkgnrm, rvol, aslh, bslh, factq1,  &
!& glocal, mfd2ft, ntotfd, nd1vks, kfft0,  &
!& fft3x, fft3y, fftwork, ijkgd, ngenh, rdelv, nel, lmixhold, &
!& kntchg, itratn, nplwcex, unr, dfmr, lplchex, aij, betakl, wrk, cmk, ip,  &
!& dmu, dfi, wn, wkerker, lvarweit )

end if


return
end subroutine




subroutine chknel( nfile, myid, nodes,                             &
& rho, mshnod, rdelv, nel, nmt0, ltimecnt )
!-----------------------------------------------------------------------
!     check No. of electrons
!-----------------------------------------------------------------------
use outfile
implicit none
integer :: nfile(*), myid, nodes
integer :: mshnod, nel, nmt0
real*8  :: rdelv
real*8  :: rho(nmt0)
logical :: ltimecnt

!-----declare local variables
integer :: m1
real*8  :: tel, telm, dbuf1r
!      real*8  :: tolera = 1.d-10
real*8  :: tolera = 1.d-01
save tolera


tel  = 0.d0
telm = 0.d0
do m1 = 1, mshnod
   tel = tel + rho( m1 )
   if( rho( m1 ) < 0.d0 ) telm = telm + rho( m1 )
end do
call gdsum(tel,1,dbuf1r)
tel = tel * rdelv
call gdsum(telm,1,dbuf1r)
telm = telm * rdelv
if( ltimecnt ) then
    if(loutfile(1)) write(nfile(1),*) ' No. of electrons (2) =', tel, telm
    if(loutfile(2)) write(nfile(2),*) ' No. of electrons (2) =', tel, telm
end if
if( abs(dble(nel)-tel) > tolera ) then
    if(loutfile(1)) write(nfile(1),*) ' error in chknel: No. of electrons =', tel
    if(loutfile(2)) write(nfile(2),*) ' error in chknel: No. of electrons =', tel
end if
tel = dble(nel)/tel
do m1 = 1, mshnod
   rho( m1 ) = rho( m1 )*tel
end do


return
end




subroutine chknud( nfile, myid, nodes,                             &
& rhoud, mshnod, rdelv, diffud, lfixud, nmt0, ltimecnt )
!-----------------------------------------------------------------------
!     check No. of electrons
!-----------------------------------------------------------------------
use outfile
implicit none

integer :: myid, nodes
integer, dimension(*) :: nfile

integer :: mshnod, nmt0
real*8  :: rdelv, diffud
real*8, dimension(nmt0) :: rhoud
logical :: lfixud
logical :: ltimecnt

!-----declare local variables
integer :: m1
real*8  :: tel, dbuf1r
!      real*8  :: tolera = 1.d-10
real*8  :: tolera = 1.d-01
save tolera


tel = 0.d0
do m1 = 1, mshnod
   tel = tel + rhoud( m1 )
end do
call gdsum(tel,1,dbuf1r)
tel = tel * rdelv
if( ltimecnt ) then
    if(loutfile(1)) write(nfile(1),*) ' diff. up & down  :',tel
    if(loutfile(2)) write(nfile(2),*) ' diff. up & down  :',tel
end if
if( lfixud .and. abs(diffud-tel) > tolera ) then
    if(loutfile(1)) write(nfile(1),*) ' error in chknud: diff. up & down  =',tel
    if(loutfile(2)) write(nfile(2),*) ' error in chknud: diff. up & down  =',tel
end if
if( lfixud .and. abs(tel) > 1.d-10 ) then
    tel = diffud/tel
    do m1 = 1, mshnod
       rhoud( m1 ) = rhoud( m1 )*tel
    end do
end if


return
end




