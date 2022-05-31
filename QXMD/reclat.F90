



module reclat_variables
!-----------------------------------------------------------------------
! type declaration of variables in reclat.f
!-----------------------------------------------------------------------
implicit none

real*8,  dimension(3,3) :: cross
real*8,  dimension(3)   :: rk1d
real*8,  dimension(3,3) :: crsbox
real*8  :: ecutbox      ! cutoff energies for BOX
logical :: lnoncollinear

save

end module




subroutine reclat_ini( nfile, myid, nodes,  &
& nplw5, nplw, nplw2, nplw3, nplwc, nplw7, nplwcs,  &
& nplw5ex, nplwex, nplw2ex, nplw3ex, nplwcex, nplw7ex, rvol,  &
& nd1v, nd1vks, hcell, rba, ecut, ecutdens, ecutsoft, ecutc, ecutorth,  &
& lvand, rctmax, kfft1b, kfft2b, kfft3b, kfft0b,  &
& pwscale, lvshape )
!-----------------------------------------------------------------------
!     set reciprocal vectors & count the number of plane waves
!-----------------------------------------------------------------------
use reclat_variables
implicit none

integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: nplw5, nplw, nplw2, nplw3, nplwc, nplw7, nplwcs  ! the number of plane waves
integer :: nplw5ex, nplwex, nplw2ex, nplw3ex, nplwcex, nplw7ex
real*8  :: rvol          ! volume of supercell
integer, dimension(3) :: nd1v     ! the number of FD  meshes
integer, dimension(3) :: nd1vks   ! the number of FFT meshes
real*8,  dimension(3,3) :: hcell  ! supercell vectors
real*8,  dimension(3,3) :: rba
real*8  :: ecut          ! cutoff energy for wavefunctions
real*8  :: ecutdens      ! cutoff energy for charge density
real*8  :: ecutsoft      ! cutoff energy for soft part of density
real*8  :: ecutc         ! cutoff energy for charge density mixing
real*8  :: ecutorth      ! cutoff energy for orthogonalization in CG iteration
logical :: lvand
real*8  :: rctmax
integer :: kfft1b, kfft2b, kfft3b, kfft0b
real*8  :: pwscale
logical :: lvshape

!-----declare local variables
integer :: i, j


call get_lnoncollinear( lnoncollinear )

!--- set supercell
do i = 1, 3
do j = 1, 3
   rba(j,i) = hcell(j,i)*dble(nd1vks(i))/dble(nd1v(i))
end do
end do


nplw5 = 0
call strvec( nfile, myid, nodes,  &
& nplw5, nplw, nplw2, nplwc, nplw7, nplwcs,  &
& ecutdens, ecut, ecutsoft, ecutc, ecutorth,  &
& rba, rvol, cross, rk1d, .true., pwscale, lvshape, lnoncollinear )

!-----get box-grid
!if( lvand ) then
!    call get_vandbox( nfile, myid, nodes,  &
!& nplw3, kfft1b, kfft2b, kfft3b, kfft0b, rctmax, nd1vks, ecutdens,  &
!& nplw5, cross, rk1d, crsbox, ecutbox )
!  else
    nplw3 = 0
!end if


!-----allocate memory for variables related to nplw3
call nplw3_alloc( nfile, myid, nodes, nplw3 )


nplwex  = 2*( nplw  + 1 )
nplw5ex = 2*( nplw5 + 1 )
nplw2ex = 2*( nplw2 + 1 )
nplw3ex = 2*( nplw3 + 1 )
nplwcex = 2*( nplwc + 1 )
nplw7ex = 2*( nplw7 + 1 )


!-----copy rba to datadump.f90
call rba_in_datadump( rba )

!-----copy rba to wannier.f90
!call rba_in_wannier( rba )

!-----copy rba to efield.f90
!call rba_in_efield( rba )


return
end subroutine




subroutine get_ecutbox( nfile, myid, nodes,  &
& ecutbox, ecutdens, crsbox, kfft1b, kfft2b, kfft3b )
!-----------------------------------------------------------------------
!     determine cutoff energy for BOX : ecutbox
!-----------------------------------------------------------------------
implicit none

integer :: myid, nodes
integer, dimension(*) :: nfile
real*8  :: ecutbox, ecutdens
real*8,  dimension(3,3) :: crsbox
integer :: kfft1b, kfft2b, kfft3b

!-----declare local variables
integer, dimension(3) :: IIP1 = (/ 2, 3, 1 /)
integer, dimension(3) :: IIP2 = (/ 3, 1, 2 /)
integer, dimension(3) :: IIP3 = (/ 1, 2, 3 /)
integer, dimension(3) :: nd1
integer :: i, ip1, ip2, ip3
real*8  :: a1x, a1y, a1z, aaa, pk3, gmax
save IIP1, IIP2, IIP3


nd1(1) = ( kfft1b - 1 )/2
nd1(2) = ( kfft2b - 1 )/2
nd1(3) = ( kfft3b - 1 )/2

ecutbox = ecutdens
do i = 1, 3
   IP1 = IIP1(i)
   IP2 = IIP2(i)
   IP3 = IIP3(i)
   A1X = crsbox(2,IP1)*crsbox(3,IP2) - crsbox(3,IP1)*crsbox(2,IP2)
   A1Y = crsbox(3,IP1)*crsbox(1,IP2) - crsbox(1,IP1)*crsbox(3,IP2)
   A1Z = crsbox(1,IP1)*crsbox(2,IP2) - crsbox(2,IP1)*crsbox(1,IP2)
   AAA = A1X*A1X + A1Y*A1Y + A1Z*A1Z
   AAA = 1.0D0/SQRT(AAA)
   A1X = A1X*AAA
   A1Y = A1Y*AAA
   A1Z = A1Z*AAA
   PK3 = crsbox(1,IP3)*A1X + crsbox(2,IP3)*A1Y + crsbox(3,IP3)*A1Z

   gmax    = pk3*( dble(nd1(i)) + 0.99d0 )
   ecutbox = min( ecutbox, gmax*gmax )
end do


return
end subroutine




subroutine reclat( nfile, myid, nodes,  &
& nd1vks, rba, ecut, ecutdens, ecutsoft, ecutorth,  &
& nplw5, nplw, nplw2, nplw3, nplwc, nplw7, nplwcs,  &
& nplw5ex, nplwex, nplw2ex, nplw3ex, nplwcex, nplw7ex, kfft0d,  &
& nga, ngb, ngc, recnrm, ijkgd, ngenh, lplh,  &
& kfft1,  kfft2,  kfft3,  kfft0, ijkg,  &
& kfft1b, kfft2b, kfft3b, kfft0b,  &
& gboxx, gboxy, gboxz, ngboxa, ngboxb, ngboxc, abgvc3, ijkgb,  &
& kmax1, kmax2, kmax3, kmax1d, kmax2d, kmax3d,  &
& kmax1b, kmax2b, kmax3b, kmax1cs, kmax2cs, kmax3cs )
!-----------------------------------------------------------------------
!     set reciprocal vectors for plane waves
!-----------------------------------------------------------------------
use outfile
use reclat_variables
implicit real*8 ( a-h, o-z )

integer :: myid, nodes
integer, dimension(*) :: nfile
integer, dimension(3) :: nd1vks   ! the number of FFT meshes
real*8,  dimension(3,3) :: rba
real*8  :: ecut          ! cutoff energy for wavefunctions
real*8  :: ecutdens      ! cutoff energy for charge density
real*8  :: ecutsoft      ! cutoff energy for soft part of density
real*8  :: ecutorth      ! cutoff energy for orthogonalization in CG iteration
real*8,  dimension(0:ngenh) :: recnrm
integer, dimension(0:ngenh) :: nga, ngb, ngc
integer, dimension(-ngenh:ngenh) :: ijkgd
integer :: kfft1,  kfft2,  kfft3,  kfft0
integer, dimension(-nplw2:nplw2) :: ijkg
integer :: kfft1b, kfft2b, kfft3b, kfft0b
real*8,  dimension(0:nplw3) :: gboxx, gboxy, gboxz
integer, dimension(0:nplw3) :: ngboxa, ngboxb, ngboxc
real*8,  dimension(0:nplw3) :: abgvc3
integer, dimension(-nplw3:nplw3) :: ijkgb


!-----------------------------------------------------------------------
call set_g2ming2max_in_chgdns( recnrm(1), recnrm(nplwc) )

!-----------------------------------------------------------------------
if( .not.lnoncollinear ) then
    nplw0 = nplw
else
    nplw0 = nplw/2
end if
kmax1 = 0
kmax2 = 0
kmax3 = 0
do ig = 1, nplw0
   kmax1 = max( kmax1, iabs(nga(ig)) )
   kmax2 = max( kmax2, iabs(ngb(ig)) )
   kmax3 = max( kmax3, iabs(ngc(ig)) )
end do
kmax1cs = kmax1
kmax2cs = kmax2
kmax3cs = kmax3
do ig = nplw0 + 1, nplwcs
   kmax1cs = max( kmax1cs, iabs(nga(ig)) )
   kmax2cs = max( kmax2cs, iabs(ngb(ig)) )
   kmax3cs = max( kmax3cs, iabs(ngc(ig)) )
end do
kmax1s = kmax1cs
kmax2s = kmax2cs
kmax3s = kmax3cs
do ig = nplwcs + 1, nplw2
   kmax1s = max( kmax1s, iabs(nga(ig)) )
   kmax2s = max( kmax2s, iabs(ngb(ig)) )
   kmax3s = max( kmax3s, iabs(ngc(ig)) )
end do
kmax1d = kmax1s
kmax2d = kmax2s
kmax3d = kmax3s
do ig = nplw2 + 1, nplw5
   kmax1d = max( kmax1d, iabs(nga(ig)) )
   kmax2d = max( kmax2d, iabs(ngb(ig)) )
   kmax3d = max( kmax3d, iabs(ngc(ig)) )
end do


!-----allocate memory for reciprocal vectors
!call set_kmax_k( nfile, myid, nodes, kmax1, kmax2, kmax3 )


!--- store FFT-mesh No. 

!-----for electronic charge density
kfft1d = nd1vks(1)
kfft2d = nd1vks(2)
kfft3d = nd1vks(3)
kfft0d = kfft1d*kfft2d*kfft3d

ijkgd(0) = 1
do m = 1, nplw5
   nd1 = nga(m)
   nd2 = ngb(m)
   nd3 = ngc(m)
   if(nd1.lt.0) nd1=nd1+kfft1d
   if(nd2.lt.0) nd2=nd2+kfft2d
   if(nd3.lt.0) nd3=nd3+kfft3d
   ijkgd(m) = 1+nd1+(nd2+nd3*kfft2d)*kfft1d

   nd1 = -nga(m)
   nd2 = -ngb(m)
   nd3 = -ngc(m)
   if(nd1.lt.0) nd1=nd1+kfft1d
   if(nd2.lt.0) nd2=nd2+kfft2d
   if(nd3.lt.0) nd3=nd3+kfft3d
   ijkgd(-m) = 1+nd1+(nd2+nd3*kfft2d)*kfft1d
end do


!-----for soft part of electronic charge density
if( nplw2 == nplw5 ) then
    kfft1 = kfft1d
    kfft2 = kfft2d
    kfft3 = kfft3d
  else
    kfft1 = 2*kmax1s + 1
    kfft2 = 2*kmax2s + 1
    kfft3 = 2*kmax3s + 1
end if

!--- check FFT mesh
call chkftm( nfile, myid, kfft1, kfft2, kfft3, .false. )
kfft0 = kfft1*kfft2*kfft3

ijkg(0) = 1
do m = 1, nplw2
   nd1 = nga(m)
   nd2 = ngb(m)
   nd3 = ngc(m)
   if(nd1.lt.0) nd1=nd1+kfft1
   if(nd2.lt.0) nd2=nd2+kfft2
   if(nd3.lt.0) nd3=nd3+kfft3
   ijkg(m) = 1+nd1+(nd2+nd3*kfft2)*kfft1

   nd1 = -nga(m)
   nd2 = -ngb(m)
   nd3 = -ngc(m)
   if(nd1.lt.0) nd1=nd1+kfft1
   if(nd2.lt.0) nd2=nd2+kfft2
   if(nd3.lt.0) nd3=nd3+kfft3
   ijkg(-m) = 1+nd1+(nd2+nd3*kfft2)*kfft1
end do

!-----noncollinear magnetism
!if( lnoncollinear ) then
!    call set_ncijkg( nfile, myid, nodes,  &
!& nplw, kfft1, kfft2, kfft3 )
!end if


!--- for k-point sampling
!call set_ijkgk( nfile, myid, nodes,  &
!& kfft1, kfft2, kfft3 )


!-----for box-grid mesh for ultrasoft pseudopotentials
if( nplw3 > 0 ) then
    nplw3 = 0
    gboxx(nplw3)  = 0.d0
    gboxy(nplw3)  = 0.d0
    gboxz(nplw3)  = 0.d0
    ngboxa(nplw3) = 0
    ngboxb(nplw3) = 0
    ngboxc(nplw3) = 0
    do ig = 1, nplw5
       G1 = crsbox(1,1)*dble(nga(ig)) + crsbox(1,2)*dble(ngb(ig))  &
&         + crsbox(1,3)*dble(ngc(ig))
       G2 = crsbox(2,1)*dble(nga(ig)) + crsbox(2,2)*dble(ngb(ig))  &
&         + crsbox(2,3)*dble(ngc(ig))
       G3 = crsbox(3,1)*dble(nga(ig)) + crsbox(3,2)*dble(ngb(ig))  &
&         + crsbox(3,3)*dble(ngc(ig))
       ggs = G1*G1 + G2*G2 + G3*G3
       if( ggs <= ecutbox  ) then
           nplw3 = nplw3 + 1
           gboxx(nplw3)  = G1
           gboxy(nplw3)  = G2
           gboxz(nplw3)  = G3
           ngboxa(nplw3) = nga(ig)
           ngboxb(nplw3) = ngb(ig)
           ngboxc(nplw3) = ngc(ig)
       end if
    end do

    km1 = 0
    km2 = 0
    km3 = 0
    do ig = 1, nplw3
       km1 = max( km1, iabs(ngboxa(ig)) )
       km2 = max( km2, iabs(ngboxb(ig)) )
       km3 = max( km3, iabs(ngboxc(ig)) )
    end do
    kmax1b = km1
    kmax2b = km2
    kmax3b = km3

    abgvc3(0) = 0.d0
    do ig = 1, nplw3
       ggs = gboxx(ig)*gboxx(ig) + gboxy(ig)*gboxy(ig)  &
&          + gboxz(ig)*gboxz(ig)
       abgvc3(ig)  = sqrt(ggs)
    end do

    ijkgb(0) = 1
    do m = 1, nplw3
       nd1 = ngboxa(m)
       nd2 = ngboxb(m)
       nd3 = ngboxc(m)
       if(nd1.lt.0) nd1=nd1+kfft1b
       if(nd2.lt.0) nd2=nd2+kfft2b
       if(nd3.lt.0) nd3=nd3+kfft3b
       ijkgb(m) = 1+nd1+(nd2+nd3*kfft2b)*kfft1b

       nd1 = -ngboxa(m)
       nd2 = -ngboxb(m)
       nd3 = -ngboxc(m)
       if(nd1.lt.0) nd1=nd1+kfft1b
       if(nd2.lt.0) nd2=nd2+kfft2b
       if(nd3.lt.0) nd3=nd3+kfft3b
       ijkgb(-m) = 1+nd1+(nd2+nd3*kfft2b)*kfft1b
    end do

  else
    kfft1b = 0
    kfft2b = 0
    kfft3b = 0
    kfft0b = 0
    kmax1b = 0
    kmax2b = 0
    kmax3b = 0
end if


!--- check and print out data ------------------------------------------
do i = 1, 2
if( loutfile(i) ) then
write(nfile(i), 605) rba(1,1), rba(2,1), rba(3,1),  &
&                    rba(1,2), rba(2,2), rba(3,2),  &
&                    rba(1,3), rba(2,3), rba(3,3)
write(nfile(i), 600) cross(1,1), cross(2,1), cross(3,1),  &
&                    cross(1,2), cross(2,2), cross(3,2),  &
&                    cross(1,3), cross(2,3), cross(3,3)
write(nfile(i), 625) ecut, ecutdens, ecutsoft, ecutorth,  &
&                  kmax1,  kmax2,  kmax3,  nplw,  &
&                  kmax1s, kmax2s, kmax3s, nplw2,  &
&                  kmax1d, kmax2d, kmax3d, nplw5, nplwc, nplw7,  &
&                  kfft1,  kfft2,  kfft3,  kfft0,  &
&                  kfft1d, kfft2d, kfft3d, kfft0d

if( nplw3.gt.0 ) then
!    dabai1 = crsbox(1,1) / cross(1,1)
!    dabai2 = crsbox(2,2) / cross(2,2)
!    dabai3 = crsbox(3,3) / cross(3,3)
    dabai1 = vecratio(crsbox(1,1),cross(1,1))
    dabai2 = vecratio(crsbox(1,2),cross(1,2))
    dabai3 = vecratio(crsbox(1,3),cross(1,3))
    write(nfile(i), 630) 'Information for box-grid FFT.          '
    write(nfile(i), 605)  &
&         rba(1,1)/dabai1, rba(2,1)/dabai1, rba(3,1)/dabai1,  &
&         rba(1,2)/dabai2, rba(2,2)/dabai2, rba(3,2)/dabai2,  &
&         rba(1,3)/dabai3, rba(2,3)/dabai3, rba(3,3)/dabai3
    write(nfile(i), 600)  &
&        cross(1,1)*dabai1, cross(2,1)*dabai1, cross(3,1)*dabai1,  &
&        cross(1,2)*dabai2, cross(2,2)*dabai2, cross(3,2)*dabai2,  &
&        cross(1,3)*dabai3, cross(2,3)*dabai3, cross(3,3)*dabai3
    write(nfile(i), 635) ecutbox,  &
&                  kmax1b, kmax2b, kmax3b, nplw3,  &
&                  kfft1b, kfft2b, kfft3b, kfft0b
end if
end if
end do
605 format(/1x,'supercell for the plane wave method'/  &
&  5x,'a  = (',es13.5,',',es13.5,',',es13.5, ')' /  &
&  5x,'b  = (',es13.5,',',es13.5,',',es13.5, ')' /  &
&  5x,'c  = (',es13.5,',',es13.5,',',es13.5, ')' /)
600 format(1x,'unit reciprocal cell'/  &
&  5x,'a* = (',es13.5,',',es13.5,',',es13.5, ')' /  &
&  5x,'b* = (',es13.5,',',es13.5,',',es13.5, ')' /  &
&  5x,'c* = (',es13.5,',',es13.5,',',es13.5, ')' /)
625 format(1x, 'energy cut off'/  &
& 5x,' ecut(w.f.)=', es16.8,'    ecut(c.d.) =',es16.8 /  &
& 5x,' ecut(soft)=', es16.8,'    ecut(orth) =',es16.8 /  &
& 5x,' kmax1 =', i4,', kmax2 =', i4,', kmax3 =', i4,  &
& ',   nplw =', i7 /  &
& 5x,' kmax1s=', i4,', kmax2s=', i4,', kmax3s=', i4,  &
& ',   nplw2=', i7 /  &
& 5x,' kmax1d=', i4,', kmax2d=', i4,', kmax3d=', i4,  &
& ',   nplw5=', i7 /  &
& 47x, 'nplwc=', i7 /  &
& 47x, 'nplw7=', i7 /  &
& 1x,'fft descreated mesh'/  &
& 5x,' kfft1 =', i4,', kfft2 =', i4,', kfft3 =', i4,  &
& ',   kfft0 =',i9 /  &
& 5x,' kfft1d=', i4,', kfft2d=', i4,', kfft3d=', i4,  &
& ',   kfft0d=',i9 /)
630 format(/a41)
635 format(1x, 'energy cut off'/  &
& 5x,'    ecut(box ) =',es16.8 /  &
& 5x,' kmax1b=', i4,', kmax2b=', i4,', kmax3b=', i4,  &
& ',   nplw3=', i7 /  &
& 1x,'fft descreated mesh'/  &
& 5x,' kfft1b=', i4,', kfft2b=', i4,', kfft3b=', i4,  &
& ',   kfft0b=',i9 /)


return
end




module strvec_variables
!-----------------------------------------------------------------------
! type declaration of allocatable variables in strvec
!-----------------------------------------------------------------------
implicit none

integer :: ifftk0 = 0
integer, allocatable, dimension(:) :: ktmp1, ktmp2, ktmp3
real*8,  allocatable, dimension(:) :: gs1, gs2, gs3, gss

!--- for sorting ------------------------------
integer, allocatable, dimension(:) :: key, ix, dstkey, x1, x2
real*8 :: dkeymx = 1000000000.d0
save ifftk0, dkeymx

end module




subroutine strvec( nfile, myid, nodes,  &
& nplw5, nplw, nplw2, nplwc, nplw7, nplwcs,  &
& ecutdens, ecut, ecutsoft, ecutc, ecutorth,  &
& rba, rvol, cross, rk1d, lflag, pwscale, lvshape, lnoncollinear )
!-----------------------------------------------------------------------
!     set reciprocal vectors for plane waves
!
!     lflag = .true.  : small cell
!     lflag = .false. : large cell in the double-grid method
!-----------------------------------------------------------------------
use outfile
use strvec_variables
implicit real*8 ( a-h, o-z )

integer :: myid
integer, dimension(*) :: nfile
integer :: nplw5, nplw, nplw2, nplwc, nplw7, nplwcs
real*8  :: ecutdens, ecut, ecutsoft, ecutc, ecutorth
real*8, dimension(3,3) :: rba
real*8  :: rvol
real*8, dimension(3,3) :: cross
real*8, dimension(3)   :: rk1d
logical :: lflag
real*8  :: pwscale
logical :: lvshape, lnoncollinear

!------declare local variables
integer :: KMAXM0, KMAXM
real*8, dimension(3,3) :: b
integer :: NUM
!      integer, allocatable, dimension(:) :: ktmp1, ktmp2, ktmp3
!      real*8,  allocatable, dimension(:) :: gs1, gs2, gs3, gss
integer :: status
integer, dimension(3) :: IIP1 = (/ 2, 3, 1 /)
integer, dimension(3) :: IIP2 = (/ 3, 1, 2 /)
integer, dimension(3) :: IIP3 = (/ 1, 2, 3 /)
integer, dimension(3) :: kmaxm1
integer, dimension(24) :: tmp1, tmp2, tmp3
integer :: NUMg
integer :: nplw0
save IIP1, IIP2, IIP3


KMAXM0 = 1
!-----------------------------------------------------------------------
!   ... reciprocal lattice ...
CALL RCIPRL( rba, B, rvol )

dpi = 2.d0*acos(-1.d0)
do i=1, 3
do j=1, 3
    cross(j,i) = dpi*B(j,i)
end do
end do


GXFFT  = sqrt( ecutdens )/dpi
GXFFT2 = GXFFT*GXFFT
!--- determination of maximum integer : KMAXM --------------------------
do i = 1, 3
   RK1D(i) = B(1,i)*B(1,i) + B(2,i)*B(2,i) + B(3,i)*B(3,i)
end do

KMAXM = 0
do i = 1, 3
   IP1 = IIP1(i)
   IP2 = IIP2(i)
   IP3 = IIP3(i)
   A1X = B(2,IP1)*B(3,IP2) - B(3,IP1)*B(2,IP2)
   A1Y = B(3,IP1)*B(1,IP2) - B(1,IP1)*B(3,IP2)
   A1Z = B(1,IP1)*B(2,IP2) - B(2,IP1)*B(1,IP2)
   AAA = A1X*A1X + A1Y*A1Y + A1Z*A1Z
   AAA = 1.0D0/SQRT(AAA)
   A1X = A1X*AAA
   A1Y = A1Y*AAA
   A1Z = A1Z*AAA
   PK3 = B(1,IP3)*A1X + B(2,IP3)*A1Y + B(3,IP3)*A1Z

   DMMZ  = GXFFT/PK3
   kmaxm1(i) = int(ABS(DMMZ))
   KMAXM = max( KMAXM, int(ABS(DMMZ)) )
end do
!-----------------------------------------------------------------------


!      ifftk = 2*(KMAXM + 1)
!      ifftk = ifftk*ifftk*ifftk/2

ifftk1 = 2*(kmaxm1(1) + 2)
ifftk2 = 2*(kmaxm1(2) + 2)
ifftk3 = 2*(kmaxm1(3) + 2)
ifftk = ifftk1*ifftk2*ifftk3/2

!-----check array sizes
if( ifftk > ifftk0 ) then

    !-----if already allocated, deallocate arrays
    if( allocated(gs1) ) then
        deallocate( gs1, gs2, gs3, gss, ktmp1, ktmp2, ktmp3, key,  &
& ix, dstkey, x1, x2, stat=status )
        dealocmem = 8.d0 * ifftk0 * 4.d0 + 4.d0 * ifftk0 * 8.d0
        if(loutfile(1)) write(nfile(1),*) nint(dealocmem/1d6),' MB, deallocated (strvec)'
        if(loutfile(2)) write(nfile(2),*) nint(dealocmem/1d6),' MB, deallocated (strvec)'
    end if

    !------allocate memory
    ifftk0 = ifftk * pwscale

    allocate( gs1(ifftk0), gs2(ifftk0), gs3(ifftk0), gss(ifftk0),  &
& ktmp1(ifftk0), ktmp2(ifftk0), ktmp3(ifftk0), key(ifftk0),  &
& ix(ifftk0), dstkey(ifftk0), x1(ifftk0), x2(ifftk0),  &
& stat=status )

    !------error trap
    status = abs(status)
    call gimax(status)
    if( status /= 0 ) call fstop( nfile, myid, nodes,  &
& 'memory allocation error in strvec' )

    alocmem = 8.d0 * ifftk0 * 4.d0 + 4.d0 * ifftk0 * 8.d0
    if(loutfile(1)) write(nfile(1),*) nint(alocmem/1d6),' MB, allocated (strvec)'
    if(loutfile(2)) write(nfile(2),*) nint(alocmem/1d6),' MB, allocated (strvec)'

end if


numg = 0
do i = KMAXM0 + 1, KMAXM + 1
   ih  = i - 1
   do j = 1, i
      jh  = j - 1
      do k = 1, j
         kh  = k - 1

         num = 0
         kji = 1
         if( k == 1 ) kji = 2
         if( k == j ) kji = kji + 2
         if( i == j ) kji = kji + 4

         num = num + 1
         tmp1(num) =  ih
         tmp2(num) =  jh
         tmp3(num) =  kh
         if( kji /= 7 ) then
             num = num + 1
             tmp1(num) =  kh
             tmp2(num) =  ih
             tmp3(num) =  jh
             num = num + 1
             tmp1(num) =  jh
             tmp2(num) =  kh
             tmp3(num) =  ih
         end if
         if( kji /= 4 ) then
             num = num + 1
             tmp1(num) = -ih 
             tmp2(num) =  jh 
             tmp3(num) =  kh
             num = num + 1
             tmp1(num) =  kh 
             tmp2(num) = -ih 
             tmp3(num) =  jh
             num = num + 1
             tmp1(num) =  jh 
             tmp2(num) =  kh 
             tmp3(num) = -ih
             if( kji /= 7 .and. kji /= 6 ) then
             if( kji /= 2 ) then
                 num = num + 1
                 tmp1(num) =  ih 
                 tmp2(num) = -jh 
                 tmp3(num) =  kh
                 num = num + 1
                 tmp1(num) =  kh 
                 tmp2(num) =  ih 
                 tmp3(num) = -jh
                 num = num + 1
                 tmp1(num) = -jh 
                 tmp2(num) =  kh 
                 tmp3(num) =  ih

                 num = num + 1
                 tmp1(num) =  ih 
                 tmp2(num) =  jh 
                 tmp3(num) = -kh
                 num = num + 1
                 tmp1(num) = -kh 
                 tmp2(num) =  ih 
                 tmp3(num) =  jh
                 num = num + 1
                 tmp1(num) =  jh 
                 tmp2(num) = -kh 
                 tmp3(num) =  ih
             end if
             if( kji /= 3 .and. kji /= 5 ) then
                 num = num + 1
                 tmp1(num) =  ih 
                 tmp2(num) =  kh 
                 tmp3(num) =  jh
                 num = num + 1
                 tmp1(num) =  jh 
                 tmp2(num) =  ih 
                 tmp3(num) =  kh
                 num = num + 1
                 tmp1(num) =  kh 
                 tmp2(num) =  jh 
                 tmp3(num) =  ih

                 num = num + 1
                 tmp1(num) = -ih 
                 tmp2(num) =  kh 
                 tmp3(num) =  jh
                 num = num + 1
                 tmp1(num) =  jh 
                 tmp2(num) = -ih 
                 tmp3(num) =  kh
                 num = num + 1
                 tmp1(num) =  kh 
                 tmp2(num) =  jh 
                 tmp3(num) = -ih
                 if( kji /= 2 ) then
                     num = num + 1
                     tmp1(num) =  ih 
                     tmp2(num) =  kh 
                     tmp3(num) = -jh
                     num = num + 1
                     tmp1(num) = -jh 
                     tmp2(num) =  ih 
                     tmp3(num) =  kh
                     num = num + 1
                     tmp1(num) =  kh 
                     tmp2(num) = -jh 
                     tmp3(num) =  ih

                     num = num + 1
                     tmp1(num) =  ih 
                     tmp2(num) = -kh 
                     tmp3(num) =  jh
                     num = num + 1
                     tmp1(num) =  jh 
                     tmp2(num) =  ih 
                     tmp3(num) = -kh
                     num = num + 1
                     tmp1(num) = -kh 
                     tmp2(num) =  jh 
                     tmp3(num) =  ih
                 end if
             end if
             end if
         end if
         do mkk = 1, num
            if( abs(tmp1(mkk)) .le. kmaxm1(1)+1 .and.  &
&               abs(tmp2(mkk)) .le. kmaxm1(2)+1 .and.  &
&               abs(tmp3(mkk)) .le. kmaxm1(3)+1      ) then
                numg = numg + 1
                ktmp1(numg) = tmp1(mkk) 
                ktmp2(numg) = tmp2(mkk) 
                ktmp3(numg) = tmp3(mkk) 
            end if
         end do
      end do
   end do
end do
num = numg
!-----------------------------------
do mkk = 1, num
   d1h = ktmp1(mkk)
   d2h = ktmp2(mkk)
   d3h = ktmp3(mkk)
   gs1(mkk) = cross(1,1)*d1h + cross(1,2)*d2h + cross(1,3)*d3h
   gs2(mkk) = cross(2,1)*d1h + cross(2,2)*d2h + cross(2,3)*d3h
   gs3(mkk) = cross(3,1)*d1h + cross(3,2)*d2h + cross(3,3)*d3h
!*                ( gs1, gs2, gs3 ) is reciprocal lattice vector

   gss(mkk) = gs1(mkk)*gs1(mkk) + gs2(mkk)*gs2(mkk)  &
&           + gs3(mkk)*gs3(mkk)
end do


!--- sorting by length of vector : gss ------------------------------
!     ix     : data pointer of record

ddmx = 0.d0
do mkk = 1, num
   ddmx = max( ddmx, abs(gss(mkk)) )
end do

fctmax = dkeymx/ddmx
do mkk = 1, num
   key(mkk) = fctmax*gss(mkk)
   ix(mkk)  = mkk
end do

call vsgar( num, key, ix, ifftk, dstkey, x1, x2, x2 )
call vsoe(  num, key, ix, ifftk )
!-----------------------------------------------------------------------


!-----count number of plane waves for charge density
nplw5 = 0
do i = 1, num
   if( gss(i).le.ecutdens ) nplw5 = nplw5 + 1
end do

!if( lflag ) then

    !-----count number of plane waves for wave functions
    nplw  = 0
    do i = 1, num
       if( gss(i).le.ecut ) nplw = nplw + 1
    end do

    !-----count number of plane waves for soft part of charge density
    nplw2  = 0
    do i = 1, num
       if( gss(i).le.ecutsoft ) nplw2 = nplw2 + 1
    end do

    !-----count number of plane waves for charge density mixing
    nplwc  = 0
    do i = 1, num
       if( gss(i).le.ecutc ) nplwc = nplwc + 1
    end do

    !-----count number of plane waves for orthogonalization in CG iteration
    nplw7  = 0
    do i = 1, num
       if( gss(i).le.ecutorth ) nplw7 = nplw7 + 1
    end do

    nplw0 = nplw
    if( lnoncollinear ) then
        nplw  = nplw * 2
        nplw7 = nplw7 * 2
    end if

    !-----allocate memory for reciprocal vectors
    call set_rciprl( nfile, myid, nodes,  &
& nplw5, nplw, nplw2, nplwc,  &
& gs1, gs2, gs3, gss, ktmp1, ktmp2, ktmp3, ix, ifftk )

    !-----set number of plave waves for ycos, ysin
    nplwcs = nplw0

    !-----set reciprocal vectors for k point sampling
!    call strveck( nfile, myid, nodes,  &
!& ecut, ecutorth, nplw, nplw2, nplw7, nplwcs,  &
!& b, gs1, gs2, gs3, gss, ktmp1, ktmp2, ktmp3,  &
!& key, ix, ifftk, dstkey, x1, x2, dkeymx )

!  else
!
!    !-----allocate memory for reciprocal vectors
!    call double_grid_alloc( nfile, myid, nodes,  &
!& nplw5, gs1, gs2, gs3, gss, ktmp1, ktmp2, ktmp3, ix, ifftk )
!
!end if


if( .not.lvshape ) then

    !------deallocate arrays
    deallocate( gs1, gs2, gs3, gss, ktmp1, ktmp2, ktmp3, key,  &
& ix, dstkey, x1, x2, stat=status )
    dealocmem = 8.d0 * ifftk0 * 4.d0 + 4.d0 * ifftk0 * 8.d0
    if(loutfile(1)) write(nfile(1),*) nint(dealocmem/1d6),' MB, deallocated (strvec)'
    if(loutfile(2)) write(nfile(2),*) nint(dealocmem/1d6),' MB, deallocated (strvec)'

    ifftk0 = 0

end if


return
end




module chkftm_param
!-----------------------------------------------------------------------
! type declaration and initialization of variables for chkftm
!
! For low-resolution calculations in van der Waals DFT,
! FFT mesh must be even.
!-----------------------------------------------------------------------
implicit none

#ifdef LIBFFTW3
integer, parameter :: npwrmx = 264
integer, dimension(npwrmx) :: npwrnm = &
 (/   1,    2,    3,    4,    5,    6,    7,    8,    9,   10,  &
&    12,   14,   15,   16,   18,   20,   21,   24,   25,   27,  &
&    28,   30,   32,   35,   36,   40,   42,   45,   48,   49,  &
&    50,   54,   56,   60,   63,   64,   70,   72,   75,   80,  &
&    81,   84,   90,   96,   98,  100,  105,  108,  112,  120,  &
&   125,  126,  128,  135,  140,  144,  147,  150,  160,  162,  &
&   168,  175,  180,  189,  192,  196,  200,  210,  216,  224,  &
&   225,  240,  243,  245,  250,  252,  256,  270,  280,  288,  &
&   294,  300,  315,  320,  324,  336,  343,  350,  360,  375,  &
&   378,  384,  392,  400,  405,  420,  432,  441,  448,  450,  &
&   480,  486,  490,  500,  504,  512,  525,  540,  560,  567,  &
&   576,  588,  600,  625,  630,  640,  648,  672,  675,  686,  &
&   700,  720,  729,  735,  750,  756,  768,  784,  800,  810,  &
&   840,  864,  875,  882,  896,  900,  945,  960,  972,  980,  &
&  1000, 1008, 1024, 1029, 1050, 1080, 1120, 1125, 1134, 1152,  &
&  1176, 1200, 1215, 1225, 1250, 1260, 1280, 1296, 1323, 1344,  &
&  1350, 1372, 1400, 1440, 1458, 1470, 1500, 1512, 1536, 1568,  &
&  1575, 1600, 1620, 1680, 1701, 1715, 1728, 1750, 1764, 1792,  &
&  1800, 1875, 1890, 1920, 1944, 1960, 2000, 2016, 2025, 2048,  &
&  2058, 2100, 2160, 2187, 2205, 2240, 2250, 2268, 2304, 2352,  &
&  2400, 2401, 2430, 2450, 2500, 2520, 2560, 2592, 2625, 2646,  &
&  2688, 2700, 2744, 2800, 2835, 2880, 2916, 2940, 3000, 3024,  &
&  3072, 3087, 3125, 3136, 3150, 3200, 3240, 3360, 3375, 3402,  &
&  3430, 3456, 3500, 3528, 3584, 3600, 3645, 3675, 3750, 3780,  &
&  3840, 3888, 3920, 3969, 4000, 4032, 4050, 4096, 4116, 4200,  &
&  4320, 4374, 4375, 4410, 4480, 4500, 4536, 4608, 4704, 4725,  &
&  4800, 4802, 4860, 4900 /)
#else
#ifdef LIBFFTW
integer, parameter :: npwrmx = 225
integer, dimension(npwrmx) :: npwrnm = &
 (/   1,    2,    3,    4,    5,    6,    7,    8,    9,   10,  &
&    11,   12,   13,   14,   15,   16,   18,   20,   21,   22,  &
&    24,   25,   26,   27,   28,   30,   32,   33,   35,   36,  &
&    39,   40,   42,   44,   45,   48,   49,   50,   52,   54,  &
&    55,   56,   60,   63,   64,   65,   66,   70,   72,   75,  &
&    77,   78,   80,   81,   84,   88,   90,   91,   96,   98,  &
&    99,  100,  104,  105,  108,  110,  112,  117,  120,  125,  &
&   126,  128,  130,  132,  135,  140,  144,  147,  150,  154,  &
&   156,  160,  162,  165,  168,  175,  176,  180,  182,  189,  &
&   192,  195,  196,  198,  200,  208,  210,  216,  220,  224,  &
&   225,  231,  234,  240,  243,  245,  250,  252,  256,  260,  &
&   264,  270,  273,  275,  280,  288,  294,  297,  300,  308,  &
&   312,  315,  320,  324,  325,  330,  336,  343,  350,  351,  &
&   352,  360,  364,  375,  378,  384,  385,  390,  392,  396,  &
&   400,  405,  416,  420,  432,  440,  441,  448,  450,  455,  &
&   462,  468,  480,  486,  490,  495,  500,  504,  512,  520,  &
&   525,  528,  539,  540,  546,  550,  560,  567,  576,  585,  &
&   588,  594,  600,  616,  624,  625,  630,  637,  640,  648,  &
&   650,  660,  672,  675,  686,  693,  700,  702,  704,  720,  &
&   728,  729,  735,  750,  756,  768,  770,  780,  784,  792,  &
&   800,  810,  819,  825,  832,  840,  864,  875,  880,  882,  &
&   891,  896,  900,  910,  924,  936,  945,  960,  972,  975,  &
&   980,  990, 1000, 1008, 1024 /)
#else
#ifdef SSL2VP
integer, parameter :: npwrmx = 60
integer, dimension(npwrmx) :: npwrnm = &
        (/    1,    2,    3,    4,    5,    6,    7,    8,    9,  &
&            10,   12,   14,   15,   16,   18,   20,   21,   24,  &
&            28,   30,   35,   36,   40,   42,   45,   48,   56,  &
&            60,   63,   70,   72,   80,   84,   90,  105,  112,  &
&           120,  126,  140,  144,  168,  180,  210,  240,  252,  &
&           280,  315,  336,  360,  420,  504,  560,  630,  720,  &
&           840, 1008, 1260, 1680, 2520, 5040  / )
#else
#ifdef ASLSX
integer, parameter :: npwrmx = 109
integer, dimension(npwrmx) :: npwrnm = &
  (/  2,    3,    4,    5,    6,    8,    9,   10,   12,   15,  &
&    16,   18,   20,   24,   25,   27,   30,   32,   36,   40,  &
&    45,   48,   50,   54,   60,   64,   72,   75,   80,   81,  &
&    90,   96,  100,  108,  120,  125,  128,  135,  144,  150,  &
&   160,  162,  180,  192,  200,  216,  225,  240,  243,  250,  &
&   256,  270,  288,  300,  320,  324,  360,  375,  384,  400,  &
&   405,  432,  450,  480,  486,  500,  512,  540,  576,  600,  &
&   625,  640,  648,  675,  720,  729,  750,  768,  800,  810,  &
&   864,  900,  960,  972, 1000, 1024, 1080, 1125, 1152, 1200,  &
&  1215, 1250, 1280, 1296, 1350, 1440, 1458, 1500, 1536, 1600,  &
&  1620, 1728, 1800, 1875, 1920, 1944, 2000, 2025, 2048 /)
#else
#endif
#endif
#endif
#endif

integer :: imev0 = 1

save


end module




subroutine chkftm( nfile, myid, kfft1, kfft2, kfft3, leven )
!-----------------------------------------------------------------------
!     check No. of FFT meshes
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid
integer :: kfft1, kfft2, kfft3
logical :: leven


call chkftm_org( nfile, myid, kfft1, kfft2, kfft3, leven )


return
end




subroutine chkftm_org( nfile, myid, kfft1, kfft2, kfft3, leven )
!-----------------------------------------------------------------------
!     check No. of FFT meshes
!-----------------------------------------------------------------------
use outfile
use chkftm_param
implicit real*8 ( a-h, o-z )
integer, dimension(*) :: nfile
logical :: leven

integer, dimension(3) :: kfftdm
real*8,  dimension(3) :: ratio
integer :: imev


if( leven ) then
    imev = imev0
  else
    imev = 1
end if

!--- original ratio
ratio(1) = 1.d0
ratio(2) = dble(kfft2)/dble(kfft1)
ratio(3) = dble(kfft3)/dble(kfft1)

kfftdm(1) = kfft1
kfftdm(2) = kfft2
kfftdm(3) = kfft3
ik = 1
do i = 1, npwrmx
#ifdef SSL2VP
   if( npwrnm(i) >= kfftdm(ik) .and. mod(npwrnm(i),2) == 0 .and. &
&      mod(npwrnm(i),imev) == 0 ) then
#else
   if( npwrnm(i) >= kfftdm(ik) .and. mod(npwrnm(i),imev) == 0 ) then
#endif
       kfftdm(ik) = npwrnm(i)
       exit
   end if
end do
do ik = 2, 3
   do i = 2, npwrmx
!            if( dble(npwrnm(i))/dble(kfftdm(1)) .ge. ratio(ik) ) then
!                dif1 = dble(npwrnm(i-1))/dble(kfftdm(1)) - ratio(ik)
!                dif2 = dble(npwrnm(i  ))/dble(kfftdm(1)) - ratio(ik)
!                if( npwrnm(i-1).lt.kfftdm(ik) .or.
!     &              abs(dif1).gt.abs(dif2) ) then
!                    kfftdm(ik) = npwrnm(i)
!                  else
!                    kfftdm(ik) = npwrnm(i-1)
!                end if
      if( npwrnm(i) >= kfftdm(ik) .and. mod(npwrnm(i),imev) == 0 ) then
          kfftdm(ik) = npwrnm(i)
          exit
      end if
   end do
end do
if( kfft2.eq.kfft1 ) kfftdm(2) = kfftdm(1)
if( kfft3.eq.kfft1 ) kfftdm(3) = kfftdm(1)

if( kfftdm(1).ne.kfft1 .or. kfftdm(2).ne.kfft2 .or.  &
&   kfftdm(3).ne.kfft3 ) then
    do i = 1, 2
    if( loutfile(i) ) then
           write(nfile(i),'(1x,a34,3i5,a5,3i5)')  &
&       '*** FFT meshes have been changed: ', kfft1, kfft2, kfft3,  &
&           ' --->', kfftdm(1), kfftdm(2), kfftdm(3)
    end if
    end do
end if
kfft1 = kfftdm(1)
kfft2 = kfftdm(2)
kfft3 = kfftdm(3)


return
end




subroutine rciprl( a, b, v )
!-----------------------------------------------------------------------
!     volume and reciprocal lattice vectors
!-----------------------------------------------------------------------
!  ( input )
!     a     ...... supercell vectors
!  ( output )
!     b     ...... primitive reciprocal lattice vectors / 2*pai
!     v     ...... volume of supercell
!-----------------------------------------------------------------------
!   Caution !!   a_i * b_j .eq. delta_ij  .not. 2*pai*delta_ij
!-----------------------------------------------------------------------
implicit none
real*8,  dimension(3,3) :: a, b
real*8  :: v

!------declare local variables
real*8  :: vr


b(1,1)  = a(2,2)*a(3,3) - a(3,2)*a(2,3)
b(2,1)  = a(3,2)*a(1,3) - a(1,2)*a(3,3)
b(3,1)  = a(1,2)*a(2,3) - a(2,2)*a(1,3)
v  = a(1,1)*b(1,1) + a(2,1)*b(2,1) + a(3,1)*b(3,1)
vr = 1.0d0/v


b(1,1) = b(1,1)*vr
b(2,1) = b(2,1)*vr
b(3,1) = b(3,1)*vr

b(1,2) = (a(2,3)*a(3,1) - a(3,3)*a(2,1))*vr
b(2,2) = (a(3,3)*a(1,1) - a(1,3)*a(3,1))*vr
b(3,2) = (a(1,3)*a(2,1) - a(2,3)*a(1,1))*vr

b(1,3) = (a(2,1)*a(3,2) - a(3,1)*a(2,2))*vr
b(2,3) = (a(3,1)*a(1,2) - a(1,1)*a(3,2))*vr
b(3,3) = (a(1,1)*a(2,2) - a(2,1)*a(1,2))*vr

v = abs( v )


return
end subroutine
