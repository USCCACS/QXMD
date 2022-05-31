



module pcc_variables
!-----------------------------------------------------------------------
! type declaration of shared variables in pcc.f
!-----------------------------------------------------------------------
implicit none

!-----variables for partial core correction (PCC)
real*8, allocatable, dimension(:,:) :: rhocr, rhocra   ! in r space
real*8, allocatable, dimension(:,:) :: drhoc, drhoca
real*8, allocatable, dimension(:)   :: rmxcor, dlcor
real*8, allocatable, dimension(:,:) :: rhocg, rhocgp   ! in g space
real*8, allocatable, dimension(:)   :: rhmur
real*8, allocatable, dimension(:)   :: tmprhocg, tmprhocga

save

end module




subroutine pcc_variables_alloc( nfile, myid, nodes,  &
& alloc_mem, ntype, mx1 )
!-----------------------------------------------------------------------
!     allocate memory for variables in vxc.f
!-----------------------------------------------------------------------
use pcc_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
real*8  :: alloc_mem
integer :: ntype
integer :: mx1

!-----declare local variables
integer :: status
real*8  :: the_mem


!------allocate memory
 allocate( rhocr(0:mx1,ntype), rhocra(0:mx1,ntype),  &
& drhoc(0:mx1,ntype), drhoca(0:mx1,ntype),  &
& rmxcor(ntype), dlcor(ntype),  &
& stat=status )

the_mem =  &
& + 8.d0 * ( (mx1+1)*ntype*4 + ntype*2 )

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'pcc_variables_alloc', .true. )


return
end subroutine




subroutine pcc_variables_alloc2( nfile, myid, nodes,  &
& alloc_mem, ntype, nplw5, lstress, pwscale )
!-----------------------------------------------------------------------
!     allocate memory for partial core correction in g space
!-----------------------------------------------------------------------
use pcc_variables
use planewave_decomp_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
real*8  :: alloc_mem
integer :: ntype
integer :: nplw5
logical :: lstress
real*8  :: pwscale

!-----declare local variables
integer :: status
real*8  :: the_mem
integer :: ngenh = 0
integer :: ngenhnod = 0
save ngenh, ngenhnod


if( nplw5+1  > ngenh     .or. &
&   nplw5nod > ngenhnod ) then

    !-----if already allocated, deallocate arrays
    if( allocated(rhocg) ) then

        the_mem = 8.d0 * ( size(rhocg) + size(rhmur) )

        !------deallocate memory
        deallocate( rhocg, rhmur, stat=status )

        if( allocated(rhocgp) ) then
            the_mem = the_mem + 8.d0 * size(rhocgp)
            deallocate( rhocgp, stat=status )
        end if

        !------error trap
        call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'pcc_variables_alloc', .true. )

    end if

    ngenh    = (nplw5+1)* pwscale
    ngenhnod = nplw5nod * pwscale
    !------allocate memory
    allocate( rhocg(ngenhnod,ntype), rhmur(2*ngenh), stat=status )

    if( status == 0 .and. lstress )  &
&       allocate( rhocgp(ngenhnod,ntype), stat=status )

    the_mem = 8.d0 * ( size(rhocg) + size(rhmur) )
    if( lstress ) the_mem = the_mem + 8.d0 * size(rhocgp)

    !------error trap
    call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'pcc_variables_alloc', .true. )

end if


return
end subroutine




subroutine pccchg( nfile, myid, nodes,  &
& rhocore, lclust, rdel, rdelg, hcell, hci, lorthrhmbc, lvacuum,  &
& ntype, nhk1, nhk2, ratm, ltimecnt, lpcci, rpcc, lpking,  &
& mshxyz, mulx1, mulx2, muly1, muly2, mulz1, mulz2,  &
& mshnx, mshny, mshnz, mshnod,  &
& mshx1, mshy1, mshz1, mshx, mshy, mshz )
!-----------------------------------------------------------------------
!    partial core charge for exchange-correlation energy
!-----------------------------------------------------------------------
use outfile
use pcc_variables
use nlpp_variables
implicit real*8 ( a-h, o-z )
dimension mshxyz(mulx1:mulx2,muly1:muly2,mulz1:mulz2)
dimension mshnx(mshnod), mshny(mshnod), mshnz(mshnod)
dimension rhocore(mshnod)
dimension rpcc(*), nfile(*), rdel(3), nhk1(*), nhk2(*),  &
&         ratm(3,*), hcell(3,3)
logical :: lorthrhmbc
real*8,  dimension(3,3) :: hci
real*8,  dimension(3,3) :: rdelg
logical   lvacuum(3)
logical   lpcci(*), lpking(*)
logical   lclust, ltimecnt
dimension ndv3nl(3)


!if( .not.lpcc ) return
if( ltimecnt ) then
    ct0 = timecnt()
end if

if( lclust ) then
!--- for atomic cluster calculation ------------------------------------
    ! --- coordinate of the corner
    ixx = mshx1
    xmn = dble(ixx-1)*rdel(1)
    iyy = mshy1
    ymn = dble(iyy-1)*rdel(2)
    izz = mshz1
    zmn = dble(izz-1)*rdel(3)
    do it = 1, ntype
    if( lpcci(it) .and. lpking(it) ) then
        rmxnl2 = rmxcor(it)*rmxcor(it)
        dlxrr  = dlcor(it)
        do icf = 1, 3
           ndv3nl(icf) = rmxcor(it)/rdel(icf) + 1.d0
        end do
        do i = nhk1(it), nhk2(it)
           xx = ratm(1,i) - xmn
           yy = ratm(2,i) - ymn
           zz = ratm(3,i) - zmn
           nxx = (xx+0.5d0*rdel(1))/rdel(1) + 1.d0
           nyy = (yy+0.5d0*rdel(2))/rdel(2) + 1.d0
           nzz = (zz+0.5d0*rdel(3))/rdel(3) + 1.d0
           do ix = max(1,nxx-ndv3nl(1)), min(mshx,nxx+ndv3nl(1))
           do iy = max(1,nyy-ndv3nl(2)), min(mshy,nyy+ndv3nl(2))
           do iz = max(1,nzz-ndv3nl(3)), min(mshz,nzz+ndv3nl(3))
              m1 = mshxyz(ix,iy,iz)
              if( m1.gt.0 ) then
                  ixx = mshx1 + ix - 1
                  xm = dble(ixx-1)*rdel(1)
                  iyy = mshy1 + iy - 1
                  ym = dble(iyy-1)*rdel(2)
                  izz = mshz1 + iz - 1
                  zm = dble(izz-1)*rdel(3)

                  xm = xm - ratm(1,i)
                  ym = ym - ratm(2,i)
                  zm = zm - ratm(3,i)
                  r2 = xm*xm + ym*ym + zm*zm
                  if( r2.le.rmxnl2 ) then
                      m = r2/(2.d0*dlxrr)
                      m = 2*m
                      d = 0.5d0*( r2/dlxrr - dble(m) )
                      rhps = d*( (d-1.d0)*rhocra(m,it)  &
&                    + rhocr(m+2,it) - rhocr(m,it) ) + rhocr(m,it)
                      rhocore( m1 ) = rhocore( m1 ) + rhps
                  end if
              end if
           end do
           end do
           end do
        end do
    end if
    end do
else
!--- for bulk calculation ----------------------------------------------
    if( lvacuum(1) ) then
        ndsx1 = 0
      else
        ndsx1 = 1
    end if
    if( lvacuum(2) ) then
        ndsx2 = 0
      else
        ndsx2 = 1
    end if
    if( lvacuum(3) ) then
        ndsx3 = 0
      else
        ndsx3 = 1
    end if

orthoif: if( lorthrhmbc ) then

!------orthorhombic super cell -----------------------------------
    ! --- coordinate of the corner
    ixx = mshx1
    xmn = dble(ixx-1)*rdel(1)
    iyy = mshy1
    ymn = dble(iyy-1)*rdel(2)
    izz = mshz1
    zmn = dble(izz-1)*rdel(3)
    do it = 1, ntype
    if( lpcci(it) .and. lpking(it) ) then
        rmxnl2 = rmxcor(it)*rmxcor(it)
        dlxrr  = dlcor(it)
        do icf = 1, 3
           ndv3nl(icf) = rmxcor(it)/rdel(icf) + 1.d0
        end do
        do i = nhk1(it), nhk2(it)
           do ixxx = 0, ndsx1
              sxm = ratm(1,i) - sign(dble(ixxx),ratm(1,i)-0.5d0)
              sxm = sxm*hcell(1,1)
              xx = sxm - xmn
              nxx = (xx+0.5d0*rdel(1))/rdel(1) + 1.d0
           do iyyy = 0, ndsx2
              sym = ratm(2,i) - sign(dble(iyyy),ratm(2,i)-0.5d0)
              sym = sym*hcell(2,2)
              yy = sym - ymn
              nyy = (yy+0.5d0*rdel(2))/rdel(2) + 1.d0
           do izzz = 0, ndsx3
              szm = ratm(3,i) - sign(dble(izzz),ratm(3,i)-0.5d0)
              szm = szm*hcell(3,3)
              zz = szm - zmn
              nzz = (zz+0.5d0*rdel(3))/rdel(3) + 1.d0
           do ix = max(1,nxx-ndv3nl(1)), min(mshx,nxx+ndv3nl(1))
           do iy = max(1,nyy-ndv3nl(2)), min(mshy,nyy+ndv3nl(2))
           do iz = max(1,nzz-ndv3nl(3)), min(mshz,nzz+ndv3nl(3))
              m1 = mshxyz(ix,iy,iz)
              if( m1.gt.0 ) then
                  ixx = mshx1 + ix - 1
                  xm = dble(ixx-1)*rdel(1)
                  iyy = mshy1 + iy - 1
                  ym = dble(iyy-1)*rdel(2)
                  izz = mshz1 + iz - 1
                  zm = dble(izz-1)*rdel(3)

                  xm = xm - sxm
                  ym = ym - sym
                  zm = zm - szm
                  r2 = xm*xm + ym*ym + zm*zm
                  if( r2.le.rmxnl2 ) then
                      m = r2/(2.d0*dlxrr)
                      m = 2*m
                      d = 0.5d0*( r2/dlxrr - dble(m) )
                      rhps = d*( (d-1.d0)*rhocra(m,it)  &
&                    + rhocr(m+2,it) - rhocr(m,it) ) + rhocr(m,it)
                      rhocore( m1 ) = rhocore( m1 ) + rhps
                  end if
              end if
           end do
           end do
           end do
           end do
           end do
           end do
        end do
    end if
    end do

else orthoif

!------non-orthorhombic super cell -------------------------------
!    dn1v1 = hcell(1,1)/rdel(1)
!    dn1v2 = hcell(2,2)/rdel(2)
!    dn1v3 = hcell(3,3)/rdel(3)
    dn1v1 = vecratio(hcell(1,1),rdelg(1,1))
    dn1v2 = vecratio(hcell(1,2),rdelg(1,2))
    dn1v3 = vecratio(hcell(1,3),rdelg(1,3))
    ixx = mshx1
    xmn = dble(ixx-1)/dn1v1
    iyy = mshy1
    ymn = dble(iyy-1)/dn1v2
    izz = mshz1
    zmn = dble(izz-1)/dn1v3
    do it = 1, ntype
    if( lpcci(it) .and. lpking(it) ) then
        rmxnl2 = rmxcor(it)*rmxcor(it)
        dlxrr  = dlcor(it)
        do icf = 1, 3
          xc = sqrt( hci(1,icf)*hci(1,icf) + hci(2,icf)*hci(2,icf)  &
&                  + hci(3,icf)*hci(3,icf) ) * rmxcor(it)
!          ndv3nl(icf) = xc*hcell(icf,icf)/rdel(icf) + 1.d0
          ndv3nl(icf) = xc*vecratio(hcell(1,icf),rdelg(1,icf)) + 1.d0
        end do
        do i = nhk1(it), nhk2(it)
           do ixxx = 0, ndsx1
              q1 = ratm(1,i) - sign(dble(ixxx),ratm(1,i)-0.5d0)
              xx = q1 - xmn
              nxx = ( xx*dn1v1 + 0.5d0 ) + 1.d0
           do iyyy = 0, ndsx2
              q2 = ratm(2,i) - sign(dble(iyyy),ratm(2,i)-0.5d0)
              yy = q2 - ymn
              nyy = ( yy*dn1v2 + 0.5d0 ) + 1.d0
           do izzz = 0, ndsx3
              q3 = ratm(3,i) - sign(dble(izzz),ratm(3,i)-0.5d0)
              zz = q3 - zmn
              nzz = ( zz*dn1v3 + 0.5d0 ) + 1.d0
              sxm = hcell(1,1)*q1 + hcell(1,2)*q2 + hcell(1,3)*q3
              sym = hcell(2,1)*q1 + hcell(2,2)*q2 + hcell(2,3)*q3
              szm = hcell(3,1)*q1 + hcell(3,2)*q2 + hcell(3,3)*q3
           do ix = max(1,nxx-ndv3nl(1)), min(mshx,nxx+ndv3nl(1))
           do iy = max(1,nyy-ndv3nl(2)), min(mshy,nyy+ndv3nl(2))
           do iz = max(1,nzz-ndv3nl(3)), min(mshz,nzz+ndv3nl(3))
              m1 = mshxyz(ix,iy,iz)
              if( m1.gt.0 ) then
                  ixx = mshx1 + ix - 1
                  d1  = dble(ixx-1)
                  iyy = mshy1 + iy - 1
                  d2  = dble(iyy-1)
                  izz = mshz1 + iz - 1
                  d3  = dble(izz-1)
                xm = rdelg(1,1)*d1 + rdelg(1,2)*d2 + rdelg(1,3)*d3
                ym = rdelg(2,1)*d1 + rdelg(2,2)*d2 + rdelg(2,3)*d3
                zm = rdelg(3,1)*d1 + rdelg(3,2)*d2 + rdelg(3,3)*d3

                  xm = xm - sxm
                  ym = ym - sym
                  zm = zm - szm
                  r2 = xm*xm + ym*ym + zm*zm
                  if( r2.le.rmxnl2 ) then
                      m = r2/(2.d0*dlxrr)
                      m = 2*m
                      d = 0.5d0*( r2/dlxrr - dble(m) )
                      rhps = d*( (d-1.d0)*rhocra(m,it)  &
&                    + rhocr(m+2,it) - rhocr(m,it) ) + rhocr(m,it)
                      rhocore( m1 ) = rhocore( m1 ) + rhps
                  end if
              end if
           end do
           end do
           end do
           end do
           end do
           end do
        end do
    end if
    end do

end if orthoif
!--- end ---------------------------------------------------------------
end if

if( ltimecnt ) then
    call gsync
    ct = timecnt()
    if(loutfile(1)) write(nfile(1),*) ' partial core charge in real space: cpu-time :', ct-ct0
    if(loutfile(2)) write(nfile(2),*) ' partial core charge in real space: cpu-time :', ct-ct0
    ct0 = ct
end if


return
end




subroutine pccchg_g( nfile, myid, nodes, ltimecnt,  &
& rhocore, ntype, nhk, nhk1, nhk2, iatoit, iatmpt, natom, nion,  &
& lpcci, lpking, nplw5, nplw5ex, nplw, kfft1d, kfft2d, kfft3d, kfft0d,  &
& nga, ngb, ngc, ijkgd, thrhgr, rvol,  &
& mshnod, mftnod, mftdsp, mfd2ft, ntotfd,  &
& fft3x, fft3y, fftwork, lstress,  & !& ycos, ysin, nplwcs,
& DCQ1, DCQ2, DCQ3, DSQ1, DSQ2, DSQ3, kmax1, kmax2, kmax3 )
!-----------------------------------------------------------------------
!    partial core charge calculated in reciprocal space
!-----------------------------------------------------------------------
use outfile
use planewave_decomp_variables
use pcc_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
logical :: ltimecnt
integer :: ntype
integer, dimension(ntype) :: nhk
integer, dimension(ntype) :: nhk1, nhk2
integer :: nion, natom
integer, dimension(natom) :: iatoit
integer, dimension(nion)  :: iatmpt
logical, dimension(ntype) :: lpcci, lpking
integer :: nplw5, nplw5ex, nplw
integer :: kfft1d, kfft2d, kfft3d, kfft0d
integer, dimension(0:nplw5) :: nga, ngb, ngc
integer, dimension(-nplw5:nplw5) :: ijkgd
real*8,  dimension(nplw5ex) :: thrhgr
real*8  :: rvol
integer :: mshnod, ntotfd
integer, dimension(*) :: mftnod, mftdsp, mfd2ft
real*8,  dimension(mshnod) :: rhocore
real*8,  dimension(*) :: fft3x, fft3y
complex*16,  dimension(*) :: fftwork
logical :: lstress
!real*8,  dimension(0:nplwcs,natom) :: ycos, ysin
integer :: kmax1, kmax2, kmax3
real*8,  dimension(-kmax1:kmax1,natom) :: DCQ1, DSQ1
real*8,  dimension(-kmax2:kmax2,natom) :: DCQ2, DSQ2
real*8,  dimension(-kmax3:kmax3,natom) :: DCQ3, DSQ3

!-----declare local variables
integer :: it, ig1, ig2, k, ia, ig, igg, k1m, k2m, k3m, ijk, ijkm
real*8  :: rvolr, DCSN, DCC12, DSC12, ycosiK, ysineK, DGGG, phar, phai
integer :: nml, inv, ifd, ift, ierrft
integer :: myid_, nodes_, nkd_
real*8  :: ct0, ct, timecnt


if( ltimecnt ) then
    ct0 = timecnt()
end if

fft3x(1:kfft0d) = 0.d0
fft3y(1:kfft0d) = 0.d0

rvolr = 1.d0/rvol

if( myid_pw == 0 ) then
    ig  = 0
    ijk = ijkgd(ig)
    do it = 1, ntype
    if( lpcci(it) .and. .not.lpking(it) ) then
        DCSN = dble( nhk(it) )*rhocg(1,it)
        fft3x(ijk) = fft3x(ijk) + DCSN*rvolr
    end if
    end do
    ig1 = 2
  else
    ig1 = 1
end if

thrhgr(1:nplw5ex) = 0.d0
do it = 1, ntype
if( lpcci(it) .and. .not.lpking(it) ) then
    do ig = 1, 2*nplw5nod
       rhmur(ig) = 0.d0
    end do
    do k = nhk1(it), nhk2(it)
       ia = iatmpt(k)

!      do ig = 1, nplwcs
!         rhmur(2*ig-1) = rhmur(2*ig-1) + ycos(ig,ia)
!         rhmur(2*ig  ) = rhmur(2*ig  ) + ysin(ig,ia)
!      end do
!      do ig = nplwcs + 1, nplw5
!      do ig = 1, nplw5
       do igg = ig1, nplw5nod
          ig = igg + nplw5nod1 - 2
          K1M = nga(ig)
          K2M = ngb(ig)
          K3M = ngc(ig)
          DCC12 = DCQ1(K1M,ia)*DCQ2(K2M,ia) -  &
&                 DSQ1(K1M,ia)*DSQ2(K2M,ia)
          DSC12 = DSQ1(K1M,ia)*DCQ2(K2M,ia) +  &
&                 DCQ1(K1M,ia)*DSQ2(K2M,ia)

          ycosiK = DCC12*DCQ3(K3M,ia) - DSC12*DSQ3(K3M,ia)
          ysineK = DSC12*DCQ3(K3M,ia) + DCC12*DSQ3(K3M,ia)

          rhmur(2*igg-1) = rhmur(2*igg-1) + ycosiK
          rhmur(2*igg  ) = rhmur(2*igg  ) + ysineK
       end do
    end do

!    do ig = 1, nplw5
    do igg = ig1, nplw5nod
       DGGG =  rhocg(igg,it)*rvolr
       phar =  rhmur(2*igg-1)*DGGG
       phai = -rhmur(2*igg  )*DGGG
       thrhgr(2*igg-1) = thrhgr(2*igg-1) + phar
       thrhgr(2*igg  ) = thrhgr(2*igg  ) + phai
    end do

end if
end do
if( nodes_pw > 1 ) then
    !-----set communicator
    call get_worldpw( myid_, nodes_ )
!    !-----global sum
!    call gdsum(thrhgr,2*nplw5,rhmur)
    do igg = 2*ig1-1, 2*nplw5nod
       rhmur(igg) = thrhgr(igg)
    end do
    call alldgatherv(rhmur,2*nplw5nod,thrhgr,nplw5excnt,nplw5exdsp)

    ig = 0
    ijk = ijkgd(ig)
    call ibcast(ijk,1,0)
    call dbcast(fft3x(ijk),1,0)

    !-----set communicator
    call get_worldkd( myid_, nodes_, nkd_ )
end if
!-----global sum
call gdsum(thrhgr(3),2*nplw5,rhmur)
do ig = 1, nplw5
   ijk  = ijkgd(ig)
   fft3x(ijk)  = fft3x(ijk)  + thrhgr(2*ig+1)
   fft3y(ijk)  = fft3y(ijk)  + thrhgr(2*ig+2)
end do
do ig = 1, nplw5
   ijkm = ijkgd(-ig)
   fft3x(ijkm) = fft3x(ijkm) + thrhgr(2*ig+1)
   fft3y(ijkm) = fft3y(ijkm) - thrhgr(2*ig+2)
end do

!        ........................
!        ... pcc(g) to pcc(r) ...
!        ........................
!          .................
!          ... constants ...
!          .................
nml = 1
inv = 2

call rfft3( inv, fft3x, fft3y, fftwork,  &
&                kfft1d, kfft2d, kfft3d, kfft0d, ierrft )

!      --- convert FFT -> FD meshes ---
do ifd = 1, ntotfd
   ift = mfd2ft(ifd)
   fft3y(ifd) = fft3x(ift)
end do

call distlc( nfile, myid, nodes,  &
& fft3x, mshnod, fft3y, mftdsp )

rhocore(1:mshnod) = rhocore(1:mshnod) + fft3x(1:mshnod)

!-----for stress calculation
!if( lstress ) then
!    call save_vext( nfile, myid, nodes, fft3x, mshnod )
!end if

if( ltimecnt ) then
    call gsync
    ct = timecnt()
    if(loutfile(1)) write(nfile(1),*) ' partial core charge in reciprocal space: cpu-time :', ct-ct0
    if(loutfile(2)) write(nfile(2),*) ' partial core charge in reciprocal space: cpu-time :', ct-ct0
    ct0 = ct
end if


return
end




subroutine pccfrc( nfile, myid, nodes,  &
& vexc, fpcc, lclust, rdel, rdelg, rdelv, hcell, hci, lorthrhmbc,  &
& lvacuum,  &
& ntype, natom, nhk1, nhk2, ratm, ltimecnt, dbuf, dbufr,  &
& lpcci, rpcc, lpking,  &
& mshxyz, mulx1, mulx2, muly1, muly2, mulz1, mulz2,  &
& mshnx, mshny, mshnz, mshnod,  &
& mshx1, mshy1, mshz1, mshx, mshy, mshz )
!-----------------------------------------------------------------------
!    partial core charge for exchange-correlation energy
!-----------------------------------------------------------------------
use outfile
use pcc_variables
implicit real*8 ( a-h, o-z )
dimension mshxyz(mulx1:mulx2,muly1:muly2,mulz1:mulz2)
dimension mshnx(mshnod), mshny(mshnod), mshnz(mshnod)
dimension fpcc(3,natom)
dimension vexc(*), rpcc(*), nfile(*), rdel(3), nhk1(*), nhk2(*),  &
&         ratm(3,*), hcell(3,3)
logical :: lorthrhmbc
real*8,  dimension(3,3) :: hci
real*8,  dimension(3,3) :: rdelg
logical   lvacuum(3)
logical   lpcci(*), lpking(*)
logical   lclust, ltimecnt
dimension ndv3nl(3)
real*8,  dimension(3*natom) :: dbuf, dbufr


if( ltimecnt ) then
    ct0 = timecnt()
end if

if( lclust ) then
!--- for atomic cluster calculation ------------------------------------
    ! --- coordinate of the corner
    ixx = mshx1
    xmn = dble(ixx-1)*rdel(1)
    iyy = mshy1
    ymn = dble(iyy-1)*rdel(2)
    izz = mshz1
    zmn = dble(izz-1)*rdel(3)
    do it = 1, ntype
    if( lpcci(it) .and. lpking(it) ) then
        rmxnl2 = rmxcor(it)*rmxcor(it)
        dlxrr  = dlcor(it)
        do icf = 1, 3
           ndv3nl(icf) = rmxcor(it)/rdel(icf) + 1.d0
        end do
        do i = nhk1(it), nhk2(it)
           xx = ratm(1,i) - xmn
           yy = ratm(2,i) - ymn
           zz = ratm(3,i) - zmn
           nxx = (xx+0.5d0*rdel(1))/rdel(1) + 1.d0
           nyy = (yy+0.5d0*rdel(2))/rdel(2) + 1.d0
           nzz = (zz+0.5d0*rdel(3))/rdel(3) + 1.d0
           do ix = max(1,nxx-ndv3nl(1)), min(mshx,nxx+ndv3nl(1))
           do iy = max(1,nyy-ndv3nl(2)), min(mshy,nyy+ndv3nl(2))
           do iz = max(1,nzz-ndv3nl(3)), min(mshz,nzz+ndv3nl(3))
              m1 = mshxyz(ix,iy,iz)
              if( m1.gt.0 ) then
                  ixx = mshx1 + ix - 1
                  xm = dble(ixx-1)*rdel(1)
                  iyy = mshy1 + iy - 1
                  ym = dble(iyy-1)*rdel(2)
                  izz = mshz1 + iz - 1
                  zm = dble(izz-1)*rdel(3)

                  xm = xm - ratm(1,i)
                  ym = ym - ratm(2,i)
                  zm = zm - ratm(3,i)
                  r2 = xm*xm + ym*ym + zm*zm
                  if( r2.le.rmxnl2 ) then
                      m = r2/(2.d0*dlxrr)
                      m = 2*m
                      d = 0.5d0*( r2/dlxrr - dble(m) )
                      rhps = d*( (d-1.d0)*drhoca(m,it)  &
&                    + drhoc(m+2,it) - drhoc(m,it) ) + drhoc(m,it)
                      vexcm1r = vexc(m1)*rhps * rdelv
                      fpcc(1,i) = fpcc(1,i) + vexcm1r * xm
                      fpcc(2,i) = fpcc(2,i) + vexcm1r * ym
                      fpcc(3,i) = fpcc(3,i) + vexcm1r * zm
                  end if
              end if
           end do
           end do
           end do
        end do
    end if
    end do
else
!--- for bulk calculation ----------------------------------------------
    if( lvacuum(1) ) then
        ndsx1 = 0
      else
        ndsx1 = 1
    end if
    if( lvacuum(2) ) then
        ndsx2 = 0
      else
        ndsx2 = 1
    end if
    if( lvacuum(3) ) then
        ndsx3 = 0
      else
        ndsx3 = 1
    end if

orthoif: if( lorthrhmbc ) then

!------orthorhombic super cell -----------------------------------
    ! --- coordinate of the corner
    ixx = mshx1
    xmn = dble(ixx-1)*rdel(1)
    iyy = mshy1
    ymn = dble(iyy-1)*rdel(2)
    izz = mshz1
    zmn = dble(izz-1)*rdel(3)
    do it = 1, ntype
    if( lpcci(it) .and. lpking(it) ) then
        rmxnl2 = rmxcor(it)*rmxcor(it)
        dlxrr  = dlcor(it)
        do icf = 1, 3
           ndv3nl(icf) = rmxcor(it)/rdel(icf) + 1.d0
        end do
        do i = nhk1(it), nhk2(it)
           do ixxx = 0, ndsx1
              sxm = ratm(1,i) - sign(dble(ixxx),ratm(1,i)-0.5d0)
              sxm = sxm*hcell(1,1)
              xx = sxm - xmn
              nxx = (xx+0.5d0*rdel(1))/rdel(1) + 1.d0
           do iyyy = 0, ndsx2
              sym = ratm(2,i) - sign(dble(iyyy),ratm(2,i)-0.5d0)
              sym = sym*hcell(2,2)
              yy = sym - ymn
              nyy = (yy+0.5d0*rdel(2))/rdel(2) + 1.d0
           do izzz = 0, ndsx3
              szm = ratm(3,i) - sign(dble(izzz),ratm(3,i)-0.5d0)
              szm = szm*hcell(3,3)
              zz = szm - zmn
              nzz = (zz+0.5d0*rdel(3))/rdel(3) + 1.d0
           do ix = max(1,nxx-ndv3nl(1)), min(mshx,nxx+ndv3nl(1))
           do iy = max(1,nyy-ndv3nl(2)), min(mshy,nyy+ndv3nl(2))
           do iz = max(1,nzz-ndv3nl(3)), min(mshz,nzz+ndv3nl(3))
              m1 = mshxyz(ix,iy,iz)
              if( m1.gt.0 ) then
                  ixx = mshx1 + ix - 1
                  xm = dble(ixx-1)*rdel(1)
                  iyy = mshy1 + iy - 1
                  ym = dble(iyy-1)*rdel(2)
                  izz = mshz1 + iz - 1
                  zm = dble(izz-1)*rdel(3)

                  xm = xm - sxm
                  ym = ym - sym
                  zm = zm - szm
                  r2 = xm*xm + ym*ym + zm*zm
                  if( r2.le.rmxnl2 ) then
                      m = r2/(2.d0*dlxrr)
                      m = 2*m
                      d = 0.5d0*( r2/dlxrr - dble(m) )
                      rhps = d*( (d-1.d0)*drhoca(m,it)  &
&                    + drhoc(m+2,it) - drhoc(m,it) ) + drhoc(m,it)
                      vexcm1r = vexc(m1)*rhps * rdelv
                      fpcc(1,i) = fpcc(1,i) + vexcm1r * xm
                      fpcc(2,i) = fpcc(2,i) + vexcm1r * ym
                      fpcc(3,i) = fpcc(3,i) + vexcm1r * zm
                  end if
              end if
           end do
           end do
           end do
           end do
           end do
           end do
        end do
    end if
    end do

else orthoif

!------non-orthorhombic super cell -------------------------------
!    dn1v1 = hcell(1,1)/rdel(1)
!    dn1v2 = hcell(2,2)/rdel(2)
!    dn1v3 = hcell(3,3)/rdel(3)
    dn1v1 = vecratio(hcell(1,1),rdelg(1,1))
    dn1v2 = vecratio(hcell(1,2),rdelg(1,2))
    dn1v3 = vecratio(hcell(1,3),rdelg(1,3))
    ixx = mshx1
    xmn = dble(ixx-1)/dn1v1
    iyy = mshy1
    ymn = dble(iyy-1)/dn1v2
    izz = mshz1
    zmn = dble(izz-1)/dn1v3
    do it = 1, ntype
    if( lpcci(it) .and. lpking(it) ) then
        rmxnl2 = rmxcor(it)*rmxcor(it)
        dlxrr  = dlcor(it)
        do icf = 1, 3
          xc = sqrt( hci(1,icf)*hci(1,icf) + hci(2,icf)*hci(2,icf)  &
&                  + hci(3,icf)*hci(3,icf) ) * rmxcor(it)
!          ndv3nl(icf) = xc*hcell(icf,icf)/rdel(icf) + 1.d0
          ndv3nl(icf) = xc*vecratio(hcell(1,icf),rdelg(1,icf)) + 1.d0
        end do
        do i = nhk1(it), nhk2(it)
           do ixxx = 0, ndsx1
              q1 = ratm(1,i) - sign(dble(ixxx),ratm(1,i)-0.5d0)
              xx = q1 - xmn
              nxx = ( xx*dn1v1 + 0.5d0 ) + 1.d0
           do iyyy = 0, ndsx2
              q2 = ratm(2,i) - sign(dble(iyyy),ratm(2,i)-0.5d0)
              yy = q2 - ymn
              nyy = ( yy*dn1v2 + 0.5d0 ) + 1.d0
           do izzz = 0, ndsx3
              q3 = ratm(3,i) - sign(dble(izzz),ratm(3,i)-0.5d0)
              zz = q3 - zmn
              nzz = ( zz*dn1v3 + 0.5d0 ) + 1.d0
              sxm = hcell(1,1)*q1 + hcell(1,2)*q2 + hcell(1,3)*q3
              sym = hcell(2,1)*q1 + hcell(2,2)*q2 + hcell(2,3)*q3
              szm = hcell(3,1)*q1 + hcell(3,2)*q2 + hcell(3,3)*q3
           do ix = max(1,nxx-ndv3nl(1)), min(mshx,nxx+ndv3nl(1))
           do iy = max(1,nyy-ndv3nl(2)), min(mshy,nyy+ndv3nl(2))
           do iz = max(1,nzz-ndv3nl(3)), min(mshz,nzz+ndv3nl(3))
              m1 = mshxyz(ix,iy,iz)
              if( m1.gt.0 ) then
                  ixx = mshx1 + ix - 1
                  d1  = dble(ixx-1)
                  iyy = mshy1 + iy - 1
                  d2  = dble(iyy-1)
                  izz = mshz1 + iz - 1
                  d3  = dble(izz-1)
                xm = rdelg(1,1)*d1 + rdelg(1,2)*d2 + rdelg(1,3)*d3
                ym = rdelg(2,1)*d1 + rdelg(2,2)*d2 + rdelg(2,3)*d3
                zm = rdelg(3,1)*d1 + rdelg(3,2)*d2 + rdelg(3,3)*d3

                  xm = xm - sxm
                  ym = ym - sym
                  zm = zm - szm
                  r2 = xm*xm + ym*ym + zm*zm
                  if( r2.le.rmxnl2 ) then
                      m = r2/(2.d0*dlxrr)
                      m = 2*m
                      d = 0.5d0*( r2/dlxrr - dble(m) )
                      rhps = d*( (d-1.d0)*drhoca(m,it)  &
&                    + drhoc(m+2,it) - drhoc(m,it) ) + drhoc(m,it)
                      vexcm1r = vexc(m1)*rhps * rdelv
                      fpcc(1,i) = fpcc(1,i) + vexcm1r * xm
                      fpcc(2,i) = fpcc(2,i) + vexcm1r * ym
                      fpcc(3,i) = fpcc(3,i) + vexcm1r * zm
                  end if
              end if
           end do
           end do
           end do
           end do
           end do
           end do
        end do
    end if
    end do

end if orthoif
!--- end ---------------------------------------------------------------
end if

!do i = 1, nhk2(ntype)
!   dbuf(3*i-2) = fpcc(1,i)
!   dbuf(3*i-1) = fpcc(2,i)
!   dbuf(3*i-0) = fpcc(3,i)
!end do
!
!call gdsum(dbuf,3*nhk2(ntype),dbufr)
!
!do i = 1, nhk2(ntype)
!   fpcc(1,i) = dbuf(3*i-2)*rdelv
!   fpcc(2,i) = dbuf(3*i-1)*rdelv
!   fpcc(3,i) = dbuf(3*i-0)*rdelv
!end do

if( ltimecnt ) then
    call gsync
    ct = timecnt()
    if(loutfile(1)) write(nfile(1),*) '      by partial core correction ',  &
&                             ':          :', ct - ct0
    if(loutfile(2)) write(nfile(2),*) '      by partial core correction ',  &
&                             ':          :', ct - ct0
    ct0 = ct
end if


return
end




subroutine pccfrc_g( nfile, myid, nodes,  &
& glocal, fpcc, strpcc, lcstress, ntype, nhk1, nhk2, iatoit,  &
& iatmpt, natom, nion, lpcci, lpking, nplw5, nplw5ex, nplw,  &
& nga, ngb, ngc, gx, gy, gz, thrhgr,  & !ycos, ysin, nplwcs,  &
& kmax1, kmax2, kmax3, DCQ1, DCQ2, DCQ3, DSQ1, DSQ2, DSQ3,  &
& kfft1d, kfft2d, kfft3d, kfft0d, mfd2ft, ntotfd, ijkgd, ngenh,  &
& fft3x, fft3y, fftwork, ltimecnt )
!-----------------------------------------------------------------------
!    fpcc : force by partial core correction in reciprocal space
!-----------------------------------------------------------------------
use outfile
use planewave_decomp_variables
use pcc_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: ntotfd
real*8,  dimension(ntotfd) :: glocal
integer :: ntype
logical :: lcstress
integer, dimension(ntype) :: nhk1, nhk2
integer :: nion, natom
real*8,  dimension(3,natom) :: fpcc
real*8,  dimension(3,3)     :: strpcc
integer, dimension(natom) :: iatoit
integer, dimension(nion)  :: iatmpt
logical, dimension(ntype) :: lpcci, lpking
integer :: nplw5, nplw5ex, nplw
integer, dimension(0:nplw5) :: nga, ngb, ngc
real*8,  dimension(0:nplw5) :: gx, gy, gz
real*8,  dimension(nplw5ex) :: thrhgr
!real*8,  dimension(0:nplwcs,natom) :: ycos, ysin
integer :: kmax1, kmax2, kmax3
real*8,  dimension(-kmax1:kmax1,natom) :: DCQ1, DSQ1
real*8,  dimension(-kmax2:kmax2,natom) :: DCQ2, DSQ2
real*8,  dimension(-kmax3:kmax3,natom) :: DCQ3, DSQ3
integer :: kfft1d, kfft2d, kfft3d, kfft0d
integer :: ngenh
integer, dimension(ntotfd) :: mfd2ft
integer, dimension(-ngenh:ngenh) :: ijkgd
real*8,  dimension(*) :: fft3x, fft3y
complex*16,  dimension(*) :: fftwork
logical :: ltimecnt

!-----declare local variables
real*8,  dimension(6) :: bufst, bufstr
integer :: myid_, nodes_, nkd_
integer :: nml, inv, ifd, ift, ierrft
integer :: ig1, it, k, ia, ig, igg, k1m, k2m, k3m, ijk
real*8  :: DCC12, DSC12, zcosk, zsink, ftmp, stmp, stmp2
real*8  :: ct0, ct, timecnt


if( ltimecnt ) then
    ct0 = timecnt()
end if


!   --- constants for fft ---
nml = 1
inv = 2
fft3x(1:kfft0d) = 0.d0
fft3y(1:kfft0d) = 0.d0
do ifd = 1, ntotfd
   ift = mfd2ft(ifd)
   fft3x(ift) = glocal(ifd)
end do

call rfft3( nml, fft3x, fft3y, fftwork,  &
&                kfft1d, kfft2d, kfft3d, kfft0d, ierrft )

do ig = 0, nplw5
   ijk  = ijkgd(ig)
   thrhgr(2*ig+1) = fft3x(ijk)
   thrhgr(2*ig+2) = fft3y(ijk)
end do


if( myid_pw == 0 ) then
    ig1 = 2
  else
    ig1 = 1
end if

!-----set communicator
call get_worldpw( myid_, nodes_ )

do it = 1, ntype
if( lpcci(it) .and. .not.lpking(it) ) then

   do k = nhk1(it), nhk2(it)
      ia = iatmpt(k)
!      do ig = 1, nplwcs
!         ftmp = rhocg(ig,it)*( thrhgr(2*ig+1)*ysin(ig,ia)  &
!&                            + thrhgr(2*ig+2)*ycos(ig,ia) )
!         fpcc(1,ia) = fpcc(1,ia) + ftmp*gx(ig)
!         fpcc(2,ia) = fpcc(2,ia) + ftmp*gy(ig)
!         fpcc(3,ia) = fpcc(3,ia) + ftmp*gz(ig)
!      end do
!      do ig = nplwcs + 1, nplw5
!      do ig = 1, nplw5
      bufst(1:3) = 0.d0
      do igg = ig1, nplw5nod
         ig = igg + nplw5nod1 - 2
         K1M = nga(ig)
         K2M = ngb(ig)
         K3M = ngc(ig)
         DCC12 = DCQ1(K1M,ia)*DCQ2(K2M,ia) -  &
&                DSQ1(K1M,ia)*DSQ2(K2M,ia)
         DSC12 = DSQ1(K1M,ia)*DCQ2(K2M,ia) +  &
&                DCQ1(K1M,ia)*DSQ2(K2M,ia)
         zcosk = DCC12*DCQ3(K3M,ia) - DSC12*DSQ3(K3M,ia)
         zsink = DSC12*DCQ3(K3M,ia) + DCC12*DSQ3(K3M,ia)
         ftmp = rhocg(igg,it)*( thrhgr(2*ig+1)*zsink  &
&                             + thrhgr(2*ig+2)*zcosk )
!         fpcc(1,ia) = fpcc(1,ia) + ftmp*gx(ig)
!         fpcc(2,ia) = fpcc(2,ia) + ftmp*gy(ig)
!         fpcc(3,ia) = fpcc(3,ia) + ftmp*gz(ig)
         bufst(1) = bufst(1) + ftmp*gx(ig)
         bufst(2) = bufst(2) + ftmp*gy(ig)
         bufst(3) = bufst(3) + ftmp*gz(ig)
      end do
      if( nodes_pw > 1 ) then
          !-----global sum
          call gdsum(bufst,3,bufstr)
      end if
      fpcc(1,ia) = fpcc(1,ia) + bufst(1)*2.d0
      fpcc(2,ia) = fpcc(2,ia) + bufst(2)*2.d0
      fpcc(3,ia) = fpcc(3,ia) + bufst(3)*2.d0
   end do

end if
end do

!-----set communicator
call get_worldkd( myid_, nodes_, nkd_ )


!--- stress calculation
if( lcstress ) then

    !-----set communicator
    call get_worldpw( myid_, nodes_ )

    bufst(1:6) = 0.d0

   do it = 1, ntype
   if( lpcci(it) .and. .not.lpking(it) ) then

   do k = nhk1(it), nhk2(it)
      ia = iatmpt(k)

!      do ig = 1, nplwcs
!         stmp = rhocgp(ig,it)*( thrhgr(2*ig+1)*ycos(ig,ia)  &
!&                             - thrhgr(2*ig+2)*ysin(ig,ia) )
!         bufst(1) = bufst(1) + stmp*gx(ig)*gx(ig)
!         bufst(2) = bufst(2) + stmp*gy(ig)*gy(ig)
!         bufst(3) = bufst(3) + stmp*gz(ig)*gz(ig)
!         bufst(4) = bufst(4) + stmp*gy(ig)*gz(ig)
!         bufst(5) = bufst(5) + stmp*gz(ig)*gx(ig)
!         bufst(6) = bufst(6) + stmp*gx(ig)*gy(ig)
!      end do
!      do ig = nplwcs + 1, nplw5
!      do ig = 1, nplw5
      do igg = ig1, nplw5nod
         ig = igg + nplw5nod1 - 2
         K1M = nga(ig)
         K2M = ngb(ig)
         K3M = ngc(ig)
         DCC12 = DCQ1(K1M,ia)*DCQ2(K2M,ia) -  &
&                DSQ1(K1M,ia)*DSQ2(K2M,ia)
         DSC12 = DSQ1(K1M,ia)*DCQ2(K2M,ia) +  &
&                DCQ1(K1M,ia)*DSQ2(K2M,ia)
         zcosk = DCC12*DCQ3(K3M,ia) - DSC12*DSQ3(K3M,ia)
         zsink = DSC12*DCQ3(K3M,ia) + DCC12*DSQ3(K3M,ia)
         stmp2 = thrhgr(2*ig+1)*zcosk - thrhgr(2*ig+2)*zsink
         ftmp = rhocg(igg,it) *stmp2
         stmp = rhocgp(igg,it)*stmp2
         bufst(1) = bufst(1) + stmp*gx(ig)*gx(ig) + ftmp
         bufst(2) = bufst(2) + stmp*gy(ig)*gy(ig) + ftmp
         bufst(3) = bufst(3) + stmp*gz(ig)*gz(ig) + ftmp
         bufst(4) = bufst(4) + stmp*gy(ig)*gz(ig)
         bufst(5) = bufst(5) + stmp*gz(ig)*gx(ig)
         bufst(6) = bufst(6) + stmp*gx(ig)*gy(ig)
      end do

   end do
   !---contribution from g = 0 term
   if( myid_pw == 0 ) then
       ftmp = ( nhk2(it) - nhk1(it) + 1 )*thrhgr(1)*rhocg(1,it) * 0.5d0
       bufst(1) = bufst(1) + ftmp
       bufst(2) = bufst(2) + ftmp
       bufst(3) = bufst(3) + ftmp
   end if

   end if
   end do
   if( nodes_pw > 1 ) then
       !-----global sum
       call gdsum(bufst,6,bufstr)
   end if

   !-----set communicator
   call get_worldkd( myid_, nodes_, nkd_ )

   call gdsum(bufst,6,bufstr)
   strpcc(1,1) = strpcc(1,1) - bufst(1)*2.d0
   strpcc(2,2) = strpcc(2,2) - bufst(2)*2.d0
   strpcc(3,3) = strpcc(3,3) - bufst(3)*2.d0
   strpcc(2,3) = strpcc(2,3) - bufst(4)*2.d0
   strpcc(3,1) = strpcc(3,1) - bufst(5)*2.d0
   strpcc(1,2) = strpcc(1,2) - bufst(6)*2.d0
   strpcc(3,2) = strpcc(2,3)
   strpcc(1,3) = strpcc(3,1)
   strpcc(2,1) = strpcc(1,2)

end if

if( ltimecnt ) then
    ct = timecnt()
    if(loutfile(1)) write(nfile(1),*) '    PCC in reciprocal space ',  &
&                         ':          :', ct - ct0
    if(loutfile(2)) write(nfile(2),*) '    PCC in reciprocal space ',  &
&                         ':          :', ct - ct0
    ct0 = ct
end if


return
end




subroutine pccstr( nfile, myid, nodes,  &
& vexc, strpcc, rdel, rdelg, rdelv, hcell, hci, lorthrhmbc,  &
& lvacuum, ntype, natom, nhk1, nhk2, ratm, ltimecnt,  &
& lpcci, rpcc, lpking,  &
& mshxyz, mulx1, mulx2, muly1, muly2, mulz1, mulz2,  &
& mshnx, mshny, mshnz, mshnod,  &
& mshx1, mshy1, mshz1, mshx, mshy, mshz )
!-----------------------------------------------------------------------
!    partial core charge for exchange-correlation energy
!-----------------------------------------------------------------------
use outfile
use pcc_variables
implicit real*8 ( a-h, o-z )
dimension mshxyz(mulx1:mulx2,muly1:muly2,mulz1:mulz2)
dimension mshnx(mshnod), mshny(mshnod), mshnz(mshnod)
dimension strpcc(3,3)
dimension vexc(*), rpcc(*), nfile(*), rdel(3), nhk1(*), nhk2(*),  &
&         ratm(3,*), hcell(3,3)
logical :: lorthrhmbc
real*8,  dimension(3,3) :: hci
real*8,  dimension(3,3) :: rdelg
logical   lvacuum(3)
logical   lpcci(*), lpking(*)
logical   ltimecnt
dimension ndv3nl(3)
dimension dbuf(6), dbufr(6)


if( ltimecnt ) then
    ct0 = timecnt()
end if

!--- for bulk calculation ----------------------------------------------
    if( lvacuum(1) ) then
        ndsx1 = 0
      else
        ndsx1 = 1
    end if
    if( lvacuum(2) ) then
        ndsx2 = 0
      else
        ndsx2 = 1
    end if
    if( lvacuum(3) ) then
        ndsx3 = 0
      else
        ndsx3 = 1
    end if

orthoif: if( lorthrhmbc ) then

!------orthorhombic super cell -----------------------------------
    ! --- coordinate of the corner
    ixx = mshx1
    xmn = dble(ixx-1)*rdel(1)
    iyy = mshy1
    ymn = dble(iyy-1)*rdel(2)
    izz = mshz1
    zmn = dble(izz-1)*rdel(3)
    do it = 1, ntype
    if( lpcci(it) .and. lpking(it) ) then
        rmxnl2 = rmxcor(it)*rmxcor(it)
        dlxrr  = dlcor(it)
        do icf = 1, 3
           ndv3nl(icf) = rmxcor(it)/rdel(icf) + 1.d0
        end do
        do i = nhk1(it), nhk2(it)
           do ixxx = 0, ndsx1
              sxm = ratm(1,i) - sign(dble(ixxx),ratm(1,i)-0.5d0)
              sxm = sxm*hcell(1,1)
              xx = sxm - xmn
              nxx = (xx+0.5d0*rdel(1))/rdel(1) + 1.d0
           do iyyy = 0, ndsx2
              sym = ratm(2,i) - sign(dble(iyyy),ratm(2,i)-0.5d0)
              sym = sym*hcell(2,2)
              yy = sym - ymn
              nyy = (yy+0.5d0*rdel(2))/rdel(2) + 1.d0
           do izzz = 0, ndsx3
              szm = ratm(3,i) - sign(dble(izzz),ratm(3,i)-0.5d0)
              szm = szm*hcell(3,3)
              zz = szm - zmn
              nzz = (zz+0.5d0*rdel(3))/rdel(3) + 1.d0
           do ix = max(1,nxx-ndv3nl(1)), min(mshx,nxx+ndv3nl(1))
           do iy = max(1,nyy-ndv3nl(2)), min(mshy,nyy+ndv3nl(2))
           do iz = max(1,nzz-ndv3nl(3)), min(mshz,nzz+ndv3nl(3))
              m1 = mshxyz(ix,iy,iz)
              if( m1.gt.0 ) then
                  ixx = mshx1 + ix - 1
                  xm = dble(ixx-1)*rdel(1)
                  iyy = mshy1 + iy - 1
                  ym = dble(iyy-1)*rdel(2)
                  izz = mshz1 + iz - 1
                  zm = dble(izz-1)*rdel(3)

                  xm = xm - sxm
                  ym = ym - sym
                  zm = zm - szm
                  r2 = xm*xm + ym*ym + zm*zm
                  if( r2.le.rmxnl2 ) then
                      m = r2/(2.d0*dlxrr)
                      m = 2*m
                      d = 0.5d0*( r2/dlxrr - dble(m) )
                      rhps = d*( (d-1.d0)*drhoca(m,it)  &
&                    + drhoc(m+2,it) - drhoc(m,it) ) + drhoc(m,it)
                      vxcrhc = vexc(m1)*rhps
                      strpcc(1,1) = strpcc(1,1) + vxcrhc*xm*xm
                      strpcc(2,2) = strpcc(2,2) + vxcrhc*ym*ym
                      strpcc(3,3) = strpcc(3,3) + vxcrhc*zm*zm
                      strpcc(2,3) = strpcc(2,3) + vxcrhc*ym*zm
                      strpcc(3,1) = strpcc(3,1) + vxcrhc*zm*xm
                      strpcc(1,2) = strpcc(1,2) + vxcrhc*xm*ym
                  end if
              end if
           end do
           end do
           end do
           end do
           end do
           end do
        end do
    end if
    end do

else orthoif

!------non-orthorhombic super cell -------------------------------
!    dn1v1 = hcell(1,1)/rdel(1)
!    dn1v2 = hcell(2,2)/rdel(2)
!    dn1v3 = hcell(3,3)/rdel(3)
    dn1v1 = vecratio(hcell(1,1),rdelg(1,1))
    dn1v2 = vecratio(hcell(1,2),rdelg(1,2))
    dn1v3 = vecratio(hcell(1,3),rdelg(1,3))
    ixx = mshx1
    xmn = dble(ixx-1)/dn1v1
    iyy = mshy1
    ymn = dble(iyy-1)/dn1v2
    izz = mshz1
    zmn = dble(izz-1)/dn1v3
    do it = 1, ntype
    if( lpcci(it) .and. lpking(it) ) then
        rmxnl2 = rmxcor(it)*rmxcor(it)
        dlxrr  = dlcor(it)
        do icf = 1, 3
          xc = sqrt( hci(1,icf)*hci(1,icf) + hci(2,icf)*hci(2,icf)  &
&                  + hci(3,icf)*hci(3,icf) ) * rmxcor(it)
!          ndv3nl(icf) = xc*hcell(icf,icf)/rdel(icf) + 1.d0
          ndv3nl(icf) = xc*vecratio(hcell(1,icf),rdelg(1,icf)) + 1.d0
        end do
        do i = nhk1(it), nhk2(it)
           do ixxx = 0, ndsx1
              q1 = ratm(1,i) - sign(dble(ixxx),ratm(1,i)-0.5d0)
              xx = q1 - xmn
              nxx = ( xx*dn1v1 + 0.5d0 ) + 1.d0
           do iyyy = 0, ndsx2
              q2 = ratm(2,i) - sign(dble(iyyy),ratm(2,i)-0.5d0)
              yy = q2 - ymn
              nyy = ( yy*dn1v2 + 0.5d0 ) + 1.d0
           do izzz = 0, ndsx3
              q3 = ratm(3,i) - sign(dble(izzz),ratm(3,i)-0.5d0)
              zz = q3 - zmn
              nzz = ( zz*dn1v3 + 0.5d0 ) + 1.d0
              sxm = hcell(1,1)*q1 + hcell(1,2)*q2 + hcell(1,3)*q3
              sym = hcell(2,1)*q1 + hcell(2,2)*q2 + hcell(2,3)*q3
              szm = hcell(3,1)*q1 + hcell(3,2)*q2 + hcell(3,3)*q3
           do ix = max(1,nxx-ndv3nl(1)), min(mshx,nxx+ndv3nl(1))
           do iy = max(1,nyy-ndv3nl(2)), min(mshy,nyy+ndv3nl(2))
           do iz = max(1,nzz-ndv3nl(3)), min(mshz,nzz+ndv3nl(3))
              m1 = mshxyz(ix,iy,iz)
              if( m1.gt.0 ) then
                  ixx = mshx1 + ix - 1
                  d1  = dble(ixx-1)
                  iyy = mshy1 + iy - 1
                  d2  = dble(iyy-1)
                  izz = mshz1 + iz - 1
                  d3  = dble(izz-1)
                xm = rdelg(1,1)*d1 + rdelg(1,2)*d2 + rdelg(1,3)*d3
                ym = rdelg(2,1)*d1 + rdelg(2,2)*d2 + rdelg(2,3)*d3
                zm = rdelg(3,1)*d1 + rdelg(3,2)*d2 + rdelg(3,3)*d3

                  xm = xm - sxm
                  ym = ym - sym
                  zm = zm - szm
                  r2 = xm*xm + ym*ym + zm*zm
                  if( r2.le.rmxnl2 ) then
                      m = r2/(2.d0*dlxrr)
                      m = 2*m
                      d = 0.5d0*( r2/dlxrr - dble(m) )
                      rhps = d*( (d-1.d0)*drhoca(m,it)  &
&                    + drhoc(m+2,it) - drhoc(m,it) ) + drhoc(m,it)
                      vxcrhc = vexc(m1)*rhps
                      strpcc(1,1) = strpcc(1,1) + vxcrhc*xm*xm
                      strpcc(2,2) = strpcc(2,2) + vxcrhc*ym*ym
                      strpcc(3,3) = strpcc(3,3) + vxcrhc*zm*zm
                      strpcc(2,3) = strpcc(2,3) + vxcrhc*ym*zm
                      strpcc(3,1) = strpcc(3,1) + vxcrhc*zm*xm
                      strpcc(1,2) = strpcc(1,2) + vxcrhc*xm*ym
                  end if
              end if
           end do
           end do
           end do
           end do
           end do
           end do
        end do
    end if
    end do

end if orthoif
!--- end ---------------------------------------------------------------

dbuf(1) = strpcc(1,1)
dbuf(2) = strpcc(2,2)
dbuf(3) = strpcc(3,3)
dbuf(4) = strpcc(2,3)
dbuf(5) = strpcc(3,1)
dbuf(6) = strpcc(1,2)

call gdsum(dbuf,6,dbufr)

strpcc(1,1) = dbuf(1)*rdelv
strpcc(2,2) = dbuf(2)*rdelv
strpcc(3,3) = dbuf(3)*rdelv
strpcc(2,3) = dbuf(4)*rdelv
strpcc(3,1) = dbuf(5)*rdelv
strpcc(1,2) = dbuf(6)*rdelv
strpcc(3,2) = strpcc(2,3)
strpcc(1,3) = strpcc(3,1)
strpcc(2,1) = strpcc(1,2)

if( ltimecnt ) then
    call gsync
    ct = timecnt()
    if(loutfile(1)) write(nfile(1),*) '      by partial core correction ',  &
&                             ':          :', ct - ct0
    if(loutfile(2)) write(nfile(2),*) '      by partial core correction ',  &
&                             ':          :', ct - ct0
    ct0 = ct
end if


return
end




subroutine subpcc( nfile, myid, nodes,                             &
&                  mxrhcr, ntype, lpcci, rpcc, rhoc, aaa, bbb,  &
&                  ecutdens, lpking, rpking, gpkgmax, gpkgexct,  &
&                  r, pi, it, mesh,  &
&                  vlocli, xitgrd, rr, nvlcl,  &
& betaq, betar, amatrx, vmatrx, ivmatrx, noofdt )
!-----------------------------------------------------------------------
!     ... partial core correction ...
!-----------------------------------------------------------------------
use outfile
use pcc_variables
IMPLICIT REAL*8 (A-H,O-Z)

integer :: myid, nodes
integer, dimension(*) :: nfile

logical   lpcci(*)
dimension rpcc(*), rhoc(*), r(*)
logical   lpking(*)
dimension rpking(*), gpkgmax(*), gpkgexct(*)

integer :: nvlcl
real*8,  dimension(0:nvlcl)  :: vlocli, xitgrd, rr
real*8,  dimension(0:noofdt) :: betaq
real*8 :: gdel
real*8,  dimension(0:noofdt) :: betar
real*8 :: rspdel
real*8,  dimension(0:noofdt,0:noofdt) :: amatrx
real*8,  dimension(noofdt,noofdt+1)   :: vmatrx
integer, dimension(noofdt)            :: ivmatrx

!-----variables for spline interpolation : splins, splinv
integer, parameter :: mm = 50, mm1 = mm + 1, ndv = 2
integer, parameter :: lm = 7,  km = lm + 1, kqm = mm + km
real*8,  dimension(mm)       :: xs, ys
real*8,  dimension(kqm)      :: qs
real*8,  dimension(mm,mm1)   :: bs
real*8,  dimension(mm,0:ndv) :: as
real*8,  dimension(mm)       :: ws
integer, dimension(mm)       :: ip
!-----------------------------------------------------------------------
character(1), dimension(0:9) :: num =  &
&             (/ '0','1','2','3','4','5','6','7','8','9' /)
integer :: iunit
save ls, num
DATA ls / 5 /


!-----if ifout == 1, output data to ckeck
ifout = 0
!-----------------------------------------------------------------------
!--- store partial core charge by table
   do i = 1, mesh
      rhoc(i) = rhoc(i)/( 4.d0*pi*r(i)*r(i) )
   end do
   irec = 0
   do i = mesh-1, 1, -1
      if( rhoc(i) > 1.d-05 ) then
          rhcrmx = r(i+1)
          irec = i
          exit
      end if
   end do
   if( irec == 0 ) then
       lpcci(it) = .false.
       return
   end if
rhcrmx = min( rhcrmx, rpcc(it)*2.d0 )
!-----------------------------------------------------------------------
!     interpolation  by  spline
   n = 1
   xs(n) = 0.d0
   ys(n) = rhoc(1)
   rdels = rhcrmx/dble(mm-2)
   rdels = rdels*( 1.d0 - 1.d-10 )
   rbs  = rdels
   do i = 1, mesh
      if( r(i).gt.rbs ) then
          n = n + 1
          rbs = rbs + rdels
          xs(n) = r(i)
          ys(n) = rhoc(i)
          if( n >= mm ) exit
!          if( r(i) > rhcrmx ) exit
      end if
   end do
!   n = n + 1
!   xs(n) = xs(n-1) + rdels
!   ys(n) = ys(n-1)
   ll   = ls
   call splins( ll, n, xs, ys, qs, bs, as, ws, ip,  &
&               kqm, mm, mm1, ndv )

   dlxrr = rhcrmx*rhcrmx/dble(mxrhcr-2)
   xrr = 0.d0
   do ir = 1, mxrhcr - 2
      xrr = xrr + dlxrr
      xr  = sqrt(xrr)
      call splinv( ll, n, xr, sgc, qs, as, ws, kqm, mm, ndv, 0 )
      rhocr(ir,it) = sgc
      call splinv( ll, n, xr, dsg, qs, as, ws, kqm, mm, ndv, 1 )
      drhoc(ir,it) = dsg/xr
   end do
!---correction at the cutoff distance
   call splinv( ll, n, rhcrmx, sgc0, qs, as, ws, kqm, mm, ndv, 0 )
   call splinv( ll, n, rhcrmx, dsg0, qs, as, ws, kqm, mm, ndv, 1 )
   call splinv( ll, n, rhcrmx, d2s0, qs, as, ws, kqm, mm, ndv, 2 )
   nalp = 2
   alpp = 8.d0*log(10.d0)/rhcrmx**nalp
   xrr = 0.d0
   do ir = 1, mxrhcr - 2
      xrr = xrr + dlxrr
      xr  = sqrt(xrr)
      xrc = xr - rhcrmx
      fac  = exp(-alpp*(-xrc)**nalp)
      dfac = alpp*nalp*(-xrc)**(nalp-1) * fac
      sgc = rhocr(ir,it) - ( sgc0 + xrc*dsg0 + 0.5d0*xrc*xrc*d2s0 )*fac
      dsg = drhoc(ir,it)*xr - ( dsg0 + xrc*d2s0 )*fac  &
&         - ( sgc0 + xrc*dsg0 + 0.5d0*xrc*xrc*d2s0 )*dfac
      rhocr(ir,it) = sgc
      drhoc(ir,it) = dsg/xr
   end do


   rhocr(0,       it) = aaa
   rhocr(mxrhcr-1,it) = 0.d0
   rhocr(mxrhcr  ,it) = 0.d0
   drhoc(0,       it) = bbb
   drhoc(mxrhcr-1,it) = 0.d0
   drhoc(mxrhcr  ,it) = 0.d0
   rmxcor(it)  = rhcrmx
   dlcor(it)   = dlxrr
   do ir = 0, mxrhcr - 2
      rhocra(ir,it) = 2.D0*( rhocr(ir,it) + rhocr(ir+2,it)  &
&                             - 2.d0*rhocr(ir+1,it) )
      drhoca(ir,it) = 2.D0*( drhoc(ir,it) + drhoc(ir+2,it)  &
&                             - 2.d0*drhoc(ir+1,it) )
   end do
!-----------------------------------------------------------------------
!check
   elcore = 0.d0
   xrr = 0.d0
   do ir = 1, mxrhcr
      xrr = xrr + dlxrr
      xr  = sqrt(xrr)
      elcore = elcore + 2.d0*pi*xr*rhocr(ir,it)
   end do
!   if(loutfile(1)) write(nfile(1),*) ' check : radius, No. of core el.:',  &
!&                        rmxcor(it), elcore*dlxrr
!   if(loutfile(2)) write(nfile(2),*) ' check : radius, No. of core el.:',  &
!&                        rmxcor(it), elcore*dlxrr


   if( ifout == 1 .and. loutfile(1) ) then
       call allocate_unit_number( iunit )
!----- output original partial core charge
       open( iunit, file=( 'data/qm_pcc_org.it='//num(it) ),  &
&            status='unknown' )
       do i = 1, mesh
          write(iunit,'(2e14.6)') r(i), rhoc(i)
       end do
       close(iunit)

!----- output interpolated partial core charge
       open( iunit, file=( 'data/qm_pcc_int.it='//num(it) ),  &
&            status='unknown' )
       xrr = 0.d0
       do ir = 0, mxrhcr
          xr  = sqrt(xrr)
          write(iunit,'(3e14.6)') xr, rhocr(ir,it), drhoc(ir,it)
          xrr = xrr + dlxrr
       end do
       close(iunit)
       call deallocate_unit_number( iunit )
   end if

   !----- return core charge * 4*pi*r^2
   do i = 1, mesh
      rhoc(i) = rhoc(i)*( 4.d0*pi*r(i)*r(i) )
   end do

   if( .not.lpking(it) ) return
!--- Gaussian filtering ------------------------------------------------
   lltk = 0
   iitk = it

   deltr = rmxcor(it)/dble(nvlcl)
   rmxnl2 = rmxcor(it)*rmxcor(it)
   xrr = 0.d0
   do ir = 0, nvlcl
      rr(ir) = xrr
      r1 = xrr
      r2 = r1*r1
      if( r2.le.rmxnl2 ) then
          m = r2/(2.d0*dlxrr)
          m = 2*m
          d = 0.5d0*( r2/dlxrr - dble(m) )
          rhps = d*( (d-1.d0)*rhocra(m,it)  &
&                    + rhocr(m+2,it) - rhocr(m,it) ) + rhocr(m,it)
!                rhps = d*( (d-1.d0)*drhoca(m,it)
!     &                     + drhoc(m+2,it) - drhoc(m,it) ) + drhoc(m,it)
        else
          rhps = 0.d0
      end if
      vlocli(ir) = rhps*xrr*xrr
      xrr = xrr + deltr
   end do
       dkmax4 = sqrt(ecutdens) * gpkgmax(it)
       dkmax  = sqrt(ecutdens) * gpkgexct(it)
       eg4   = dkmax4*dkmax4
       egmax = dkmax*dkmax
       rmxnlc     = rmxcor(it) * rpking(it)
       rmxcor(it) = rmxnlc
!            --- Gaussian filtering ------------------------------
       call kingsm( nfile, myid, nodes,                            &
& lltk, iitk, deltr, egmax, eg4, .false., rmxnlc,  &
& vlocli ,xitgrd, rr, nvlcl, betaq, gdel, betar, rspdel, noofdt,  &
& amatrx, vmatrx, ivmatrx )
!            -----------------------------------------------------
       dlxrr  = rmxcor(it)*rmxcor(it)/dble(mxrhcr-2)
!        --- interpolation by  spline ---
       rskip = rmxcor(it)/dble(mm-2)
       n  = 1
       nr = 0
       xs(n) = rspdel*dble(nr)
       ys(n) = betar(nr)
       xrr = rskip
       do nr = 1, noofdt
          xnrr = rspdel*dble(nr)
          if( xnrr.ge.xrr ) then
              n = n + 1
              xs(n) = xnrr
              ys(n) = betar(nr)
              xrr = xrr + rskip
              if( n >= mm-1 ) exit
              if( xnrr > rmxcor(it) ) exit
          end if
       end do
       n = n + 1
       xs(n) = xs(n-1) + rskip
       ys(n) = ys(n-1)
       ll   = ls
       call splins( ll, n, xs, ys, qs, bs, as, ws, ip,  &
&               kqm, mm, mm1, ndv )

       xrr = 0.d0
       do ir = 1, mxrhcr - 2
          xrr = xrr + dlxrr
          xr  = sqrt(xrr)
        call splinv( ll, n, xr, sgc, qs, as, ws, kqm, mm, ndv, 0 )
          rhocr(ir,it) = sgc
        call splinv( ll, n, xr, dsg, qs, as, ws, kqm, mm, ndv, 1 )
          drhoc(ir,it) = dsg/xr
       end do
       rhocr(0,       it) = betar(0)
       rhocr(mxrhcr-2,it) = betar(noofdt)
       rhocr(mxrhcr-1,it) = 0.d0
       rhocr(mxrhcr  ,it) = 0.d0
       drhoc(0,       it) = 2.d0*drhoc(1,it) - drhoc(2,it)
       drhoc(mxrhcr-2,it) = drhoc(mxrhcr-3,it)
       drhoc(mxrhcr-1,it) = 0.d0
       drhoc(mxrhcr  ,it) = 0.d0
       dlcor(it)   = dlxrr
       do ir = 0, mxrhcr - 2
          rhocra(ir,it) = 2.D0*( rhocr(ir,it) + rhocr(ir+2,it)  &
&                             - 2.d0*rhocr(ir+1,it) )
          drhoca(ir,it) = 2.D0*( drhoc(ir,it) + drhoc(ir+2,it)  &
&                             - 2.d0*drhoc(ir+1,it) )
       end do
!check
   elcore = 0.d0
   xrr = 0.d0
   do ir = 1, mxrhcr
      xrr = xrr + dlxrr
      xr  = sqrt(xrr)
      elcore = elcore + 2.d0*pi*xr*rhocr(ir,it)
   end do
!   if(loutfile(1)) write(nfile(1),*) ' check : after filtering        :',  &
!&                        rmxcor(it), elcore*dlxrr
!   if(loutfile(2)) write(nfile(2),*) ' check : after filtering        :',  &
!&                        rmxcor(it), elcore*dlxrr


   if( ifout == 1 .and. loutfile(1) ) then
       call allocate_unit_number( iunit )
!----- output filterred partial core charge
       open( iunit, file=( 'data/qm_pcc_fil.it='//num(it) ),  &
&            status='unknown' )
       do nr = 0, noofdt
          write(iunit,'(2e14.6)') rspdel*dble(nr), betar(nr)
       end do
       close(iunit)

!----- output interpolated partial core charge
       open( iunit, file=( 'data/qm_pcc_fil_int.it='//num(it) ),  &
&            status='unknown' )
       xrr = 0.d0
       do ir = 0, mxrhcr
          xr  = sqrt(xrr)
          write(iunit,'(3e14.6)') xr, rhocr(ir,it), drhoc(ir,it)
          xrr = xrr + dlxrr
       end do
       close(iunit)
       call deallocate_unit_number( iunit )
   end if



!         rmxnl2 = rmxcor(it)*rmxcor(it)
!         do i = 1, mesh
!            r2 = r(i)*r(i)
!            if( r2.le.rmxnl2 ) then
!                m = r2/(2.d0*dlxrr)
!                m = 2*m
!                d = 0.5d0*( r2/dlxrr - dble(m) )
!                rhps = d*( (d-1.d0)*rhocra(m,it)
!     &                     + rhocr(m+2,it) - rhocr(m,it) ) + rhocr(m,it)
!!                rhps = d*( (d-1.d0)*drhoca(m,it)
!!     &                     + drhoc(m+2,it) - drhoc(m,it) ) + drhoc(m,it)
!              else
!                rhps = 0.d0
!            end if
!            rhoc(i) = rhps*( 4.d0*pi*r(i)*r(i) )
!         end do


return
end




subroutine subpcc_g( nfile, myid, nodes,                           &
& ntype, lpcci, lpking,  &
& recnrm, recnrmex, nplw5, vlocli, xitgrd, rr, nvlcl, lstress,  &
& lvshape )
!-----------------------------------------------------------------------
!     Tables for partial core charge density in reciprocal space
!-----------------------------------------------------------------------
use outfile
use planewave_decomp_variables
use pcc_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: ntype
logical, dimension(ntype) :: lpcci, lpking
integer :: nplw5
real*8,  dimension(0:nplw5) :: recnrm
real*8,  dimension(0:nplw5) :: recnrmex
integer :: nvlcl
real*8,  dimension(0:nvlcl)  :: vlocli, xitgrd, rr
logical :: lstress
logical :: lvshape

!-----declare local variables
integer, parameter :: tmpsize = 4096
integer :: it, i, ir, m
real*8  :: rhpsmx, xrdel, xrr, r2, d, rhps
integer :: ntnod1, ntnod2, ntnod, ik
real*8  :: small, g1maxq, dlqab, pi4, ak1, cl
real*8  :: akdr, sindr, cosdr, akrr, singr, cosgr, sint, cost, qbesl0
integer :: ig, igg, ig1
real*8  :: xr, dlxrr, akrrr, qbesl1
real*8  :: alocmem
integer :: status


!-----if not allocated, allocate arrays
if( .not.allocated(tmprhocg) ) then
    !------allocate memory
    allocate( tmprhocg(0:tmpsize), tmprhocga(0:tmpsize), stat=status )

    !------error trap
    status = abs(status)
    call gimax(status)
    if( status /= 0 ) call fstop( nfile, myid, nodes,  &
& 'memory allocation error in subpcc_g' )

    alocmem = 8.d0 * ( size(tmprhocg) + size(tmprhocga) )
    if(loutfile(1)) write(nfile(1),*) nint(alocmem/1d0),' B, allocated (subpcc_g)'
    if(loutfile(2)) write(nfile(2),*) nint(alocmem/1d0),' B, allocated (subpcc_g)'
end if


itdo: do it = 1, ntype
if( lpcci(it) .and. .not.lpking(it) ) then
!=======================================================================

    !-----------------------------------------------------------------------
    !     set cutoff length
    rhpsmx = rmxcor(it)
    xrdel  = rhpsmx/dble(nvlcl)
    dlxrr  = dlcor(it)

    do ir = 0, nvlcl
       xrr = xrdel*dble(ir)
       r2  = xrr*xrr
       m = r2/(2.d0*dlxrr)
       m = 2*m
       d = 0.5d0*( r2/dlxrr - dble(m) )
       rhps = d*( (d-1.d0)*rhocra(m,it)  &
&                        + rhocr(m+2,it) - rhocr(m,it) ) + rhocr(m,it)

       rr(ir) = xrr
       vlocli(ir) = rhps * r2
    end do

    small=1.0d-5
    !-----------------------------------------------------------------------
    !  create table
    g1maxq = sqrt(recnrm(nplw5))
    dlqab = g1maxq*(1.d0+10.d0/dble(tmpsize))/dble(tmpsize)
    pi4 = 4.d0*acos(-1.d0)

    tmprhocg(0:tmpsize) = 0.d0

    !-----distribute to nodes
    call divnod( tmpsize+1, nodes, myid, ntnod1, ntnod2, ntnod )

!         do ik = 0, tmpsize
    do ik = ntnod1-1, ntnod2-1
       ak1 = dlqab * dble(ik)
       if( ak1.le.small ) then
!          -----------------------------------------
           xitgrd(0:nvlcl) = vlocli(0:nvlcl)
    !---                           BSL0( ak1*rr(ir) ) = 1
           call INTGB3( nvlcl, xrdel, xitgrd, CL )
!          -----------------------------------------
       else

            akdr = ak1*xrdel
           sindr = sin(akdr)
           cosdr = cos(akdr)
            akrr = 0.d0
           singr = 0.d0
           cosgr = 1.d0
           ir = 0
           xitgrd(ir) = vlocli(ir)
           do ir = 1, nvlcl
              akrr = akrr + akdr
              sint = singr*cosdr + cosgr*sindr
              cost = cosgr*cosdr - singr*sindr
              singr = sint
              cosgr = cost
              qbesl0 = singr/akrr
              xitgrd(ir) = vlocli(ir) * qbesl0
           end do
           call INTGB3( nvlcl, xrdel, xitgrd, CL )

       end if

       tmprhocg(ik) = CL

    end do
    !-----unify tmprhocg
    call gdsum(tmprhocg,tmpsize+1,tmprhocga)

    do ik = 0, tmpsize - 2
       tmprhocga(ik) = 2.D0*( tmprhocg(ik) + tmprhocg(ik+2)  &
&                      - 2.d0*tmprhocg(ik+1) )
    end do
    do ik = tmpsize - 1, tmpsize
       tmprhocga(ik) = 0.d0
    end do

    !--- set rhocg
!    do ig = 0, nplw5
    do igg = 1, nplw5nod
       ig = igg + nplw5nod1 - 2
       ak1 = sqrt(recnrm(ig))
       m = ak1/(2.d0*dlqab)
       m = 2*m
       d = 0.5d0*( ak1/dlqab - dble(m) )
       rhps = d*( (d-1.d0)*tmprhocga(m) + tmprhocg(m+2)  &
&                        - tmprhocg(m) ) + tmprhocg(m)

       rhocg(igg,it) = pi4*rhps

    end do
    !-----------------------------------------------------------------------

    stressif: if( lstress ) then
    !-----------------------------------------------------------------------
    ! for stress calculation
    do ir = 1, nvlcl
       xr = rr(ir)
       vlocli(ir) = xr*vlocli(ir)
    end do

    tmprhocg(0:tmpsize) = 0.d0

    !-----distribute to nodes
    call divnod( tmpsize+1, nodes, myid, ntnod1, ntnod2, ntnod )

!         do ik = 0, tmpsize
    do ik = ntnod1-1, ntnod2-1
       ak1 = dlqab * dble(ik)
       if( ak1.le.small ) then
!          -----------------------------------------
           CL = 0.d0
!          -----------------------------------------
       else

            akdr = ak1*xrdel
           sindr = sin(akdr)
           cosdr = cos(akdr)
            akrr = 0.d0
           singr = 0.d0
           cosgr = 1.d0
           ir = 0
           xitgrd(ir) = 0.d0
           do ir = 1, nvlcl
              akrr = akrr + akdr
              akrrr= 1.d0/akrr
              sint = singr*cosdr + cosgr*sindr
              cost = cosgr*cosdr - singr*sindr
              singr = sint
              cosgr = cost
              qbesl1 = ( singr*akrrr - cosgr )*akrrr
              xitgrd(ir) = vlocli(ir) * qbesl1
           end do
           call INTGB3( nvlcl, xrdel, xitgrd, CL )

       end if

       tmprhocg(ik) = -CL

    end do
    !-----unify tmprhocg
    call gdsum(tmprhocg,tmpsize+1,tmprhocga)

    do ik = 0, tmpsize - 2
       tmprhocga(ik) = 2.D0*( tmprhocg(ik) + tmprhocg(ik+2)  &
&                      - 2.d0*tmprhocg(ik+1) )
    end do
    do ik = tmpsize - 1, tmpsize
       tmprhocga(ik) = 0.d0
    end do

    !--- set rhocgp
    if( myid_pw == 0 ) then
        rhocgp(1,it) = 0.d0
        ig1 = 2
      else
        ig1 = 1
    end if
!    do ig = 1, nplw5
    do igg = ig1, nplw5nod
       ig = igg + nplw5nod1 - 2
       ak1 = sqrt(recnrm(ig))
       m = ak1/(2.d0*dlqab)
       m = 2*m
       d = 0.5d0*( ak1/dlqab - dble(m) )
       rhps = d*( (d-1.d0)*tmprhocga(m) + tmprhocg(m+2)  &
&                        - tmprhocg(m) ) + tmprhocg(m)

       rhocgp(igg,it) = pi4*rhps/ak1
    end do
    !-----------------------------------------------------------------------
    end if stressif

!=======================================================================
end if
end do itdo


if( .not.lvshape ) then
    alocmem = 8.d0 * ( size(tmprhocg) + size(tmprhocga) )

    !------deallocate memory
    deallocate( tmprhocg, tmprhocga, stat=status )

    !------error trap
    status = abs(status)
    call gimax(status)
    if( status /= 0 ) call fstop( nfile, myid, nodes,  &
& 'memory deallocation error in subpcc_g' )

    if(loutfile(1)) write(nfile(1),*) nint(alocmem/1d0),' B, deallocated (subpcc_g)'
    if(loutfile(2)) write(nfile(2),*) nint(alocmem/1d0),' B, deallocated (subpcc_g)'
end if


return
end
