



subroutine dtdump( myid, nodes, v, savedata, dbuf, dbufr,  &
&               xx, nd1v, mshmx, mshmy, mshmz, nmt0, nmt3, fname1,  &
&               mshxyz, mulx1, mulx2, muly1, muly2, mulz1, mulz2,  &
&               mshnx, mshny, mshnz, mshnod,  &
&               mshx1, mshy1, mshz1, mshx, mshy, mshz,  &
&   meshx1, meshx, meshy1, meshy, meshz1, meshz, npx, npy, npz )
!-----------------------------------------------------------------------
!     dump data to file
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
real      savedata(mshmx,mshmy,mshmz)
dimension nd1v(3), v(mshnod)
dimension dbuf(nmt0), dbufr(nmt0)
dimension xx(0:nmt3)
dimension mshxyz(mulx1:mulx2,muly1:muly2,mulz1:mulz2)
dimension mshnx(mshnod), mshny(mshnod), mshnz(mshnod)
dimension meshx1(0:npx-1), meshx(0:npx-1),  &
&         meshy1(0:npy-1), meshy(0:npy-1),  &
&         meshz1(0:npz-1), meshz(0:npz-1)
character fname1*30


do iz = 1, nd1v(3)
do iy = 1, nd1v(2)
do ix = 1, nd1v(1)
   savedata( ix, iy, iz ) = 0.0
end do
end do
end do

if( myid.eq.0 ) then
    do m1 = 1, mshnod
!---  local mesh
       ixg = mshnx(m1)
       iyg = mshny(m1)
       izg = mshnz(m1)
!---  global mesh
       ix = mshx1 + ixg - 1
       iy = mshy1 + iyg - 1
       iz = mshz1 + izg - 1
       savedata( ix, iy, iz ) = v(m1)
    end do
    do node_id = 1, nodes - 1
       node_x=node_id/(npy*npz)
       node_y=mod(node_id/npz,npy)
       node_z=mod(node_id,npz)
       nrc = meshx(node_x)*meshy(node_y)*meshz(node_z)
       call cdrecv(100,dbuf,nrc,0)
       i = 0
       do izg = 1, meshz(node_z)
       do iyg = 1, meshy(node_y)
       do ixg = 1, meshx(node_x)
          i = i + 1
          ix = meshx1(node_x) + ixg - 1
          iy = meshy1(node_y) + iyg - 1
          iz = meshz1(node_z) + izg - 1
          savedata( ix, iy, iz ) = dbuf(i)
       end do
       end do
       end do
       call gsync
    end do
!--- open file and dump data
    call allocate_unit_number( nfile )
    open(nfile,file=fname1,status='unknown')
    write(nfile,'(3i5)') nd1v(1), nd1v(2), nd1v(3)
    do iz = 1, nd1v(3)
    do iy = 1, nd1v(2)
    do ix = 1, nd1v(1)
       write(nfile,'(es12.5)') savedata( ix, iy, iz )
    end do
    end do
    end do
    close(nfile)
    call deallocate_unit_number( nfile )
  else
    xx(0) = 0.d0
    do i = 1, mshnod
       xx(i) = v(i)
    end do
    i = 0
    do iz = 1, mshz
    do iy = 1, mshy
    do ix = 1, mshx
       i = i + 1
       m1 = mshxyz( ix, iy, iz )
       dbuf(i) = xx(m1)
    end do
    end do
    end do
    nrc = mshx*mshy*mshz
    do node_id = 1, nodes - 1
       if( myid.eq.node_id ) then
           call cdsend(100,dbuf,nrc,0,0)
       end if
       call gsync
    end do
end if


return
end




module datadump_variables
!-----------------------------------------------------------------------
! type declaration of variables in reclat.f
!-----------------------------------------------------------------------
implicit none

real*8,  dimension(3,3) :: rba, b
save

end module




subroutine rba_in_datadump( rba_ )
!-----------------------------------------------------------------------
!     copy rba from reclat.f
!-----------------------------------------------------------------------
use datadump_variables
implicit none
real*8,  dimension(3,3) :: rba_
real*8  :: rvol
integer :: i, j


do i = 1, 3
do j = 1, 3
   rba(j,i) = rba_(j,i)
end do
end do

!----- transpose matrix of b = inverse of rba
call rciprl( rba, b, rvol )


return
end





subroutine wrgeom( nfile, myid, nodes,  &
& nstep, lclust, ratm, frc, vatm, nhk, ntype, natom,  &
& nd1v, nd1vks, hcell )
!-----------------------------------------------------------------------
!     output supercell, atomic coordinates and atomic forces
!-----------------------------------------------------------------------
use outfile
use datadump_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: nstep
logical :: lclust
integer :: ntype
integer :: natom
real*8,  dimension(3,natom) :: ratm
real*8,  dimension(3,natom) :: frc
real*8,  dimension(3,natom) :: vatm
integer, dimension(ntype) :: nhk
integer, dimension(3) :: nd1v
integer, dimension(3) :: nd1vks
real*8,  dimension(3,3) :: hcell

!-----declare local variables
integer :: i, j, ia, it
logical :: loutedg
real*8  :: dltca, dltcb, dltcc, angalf, angbet, anggam
real*8  :: atconfx, x, y, z
real*8, dimension(3,3) :: rbarec = 0.d0
save rbarec


if( loutfile(1) ) then

    !-----supercell vectors
    call caledge( loutedg, rba, rbarec,  &
&             dltca, dltcb, dltcc, angalf, angbet, anggam )
    if( loutedg ) then
        write(nfile(6),'(i7,3es14.7,3f11.6)') nstep,  &
&             dltca, dltcb, dltcc, angalf, angbet, anggam
        write(nfile(27),'(i7,9es15.7)') nstep, rba(1:3,1:3)
    end if

end if


!--- atomic configuration ---

if( lclust ) then

    !-----if cluster, transform (real) ratm to scaled form
    do i = 1, 3
    do ia = 1, natom
       vatm(i,ia) = b(1,i)*ratm(1,ia) + b(2,i)*ratm(2,ia)  &
&                 + b(3,i)*ratm(3,ia)
    end do
    end do

  else

    !-----if bulk,
    !----- transform (scaled) ratm to real form,
    do i = 1, 3
    do ia = 1, natom
       vatm(i,ia) = hcell(i,1)*ratm(1,ia) + hcell(i,2)*ratm(2,ia)  &
&                 + hcell(i,3)*ratm(3,ia)
    end do
    end do
    !----- and transform them to scaled form for rba
    do ia = 1, natom
       x = vatm(1,ia)
       y = vatm(2,ia)
       z = vatm(3,ia)
       vatm(1,ia) = b(1,1)*x + b(2,1)*y + b(3,1)*z
       vatm(2,ia) = b(1,2)*x + b(2,2)*y + b(3,2)*z
       vatm(3,ia) = b(1,3)*x + b(2,3)*y + b(3,3)*z
    end do

end if


if( loutfile(1) ) then

    !-----atomic scaled coordinates
    atconfx = 0.1d0
    write(nfile(4),'(100i7)') nstep,  &
&              ntype, ( nhk(it), it = 1, ntype )
    write(nfile(4),'(es14.7)') atconfx
    write(nfile(4),'(9f8.5)')  &
&              ( ( min(vatm(i,ia)/atconfx,9.99999d0), i = 1, 3 ), ia = 1, natom )

end if


!--- atomic real forces ---

atconfx = 0.d0
do i  = 1, 3
do ia = 1, natom
   atconfx = max( atconfx, abs(frc(i,ia)) )
end do
end do
atconfx = atconfx / 9.999d0


if( loutfile(1) ) then

    write(nfile(8),'(100i7)') nstep,  &
&              ntype, ( nhk(it), it = 1, ntype )
    write(nfile(8),'(es14.7)') atconfx
    !--- avoid floating divide by zero
    if( abs(atconfx) < 1.d-30 ) atconfx = 1.d0
    write(nfile(8),'(9f8.5)')  &
&        ( ( frc(i,ia)/atconfx, i = 1, 3 ), ia = 1, natom )

end if


return
end




subroutine wvdump( nfile, myid, nodes,  &
& dname_eig, nstep, ibstt1, ibstt2,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, iod,  &
& mfd2ft, ntotfd, nd1vks, kfft0, rvol,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, apk, eigr, nspin, lspin,  &
& kfft1b, kfft2b, kfft3b, kfft0b, ijkgb, nplw3, ngboxa, ngboxb, ngboxc,  &
& ntype, natom, iatoit, ioa, lvand, lvandi,  &
& xmin, xmax, ymin, ymax, zmin, zmax,  &
& lcompress, ndigit )
!-----------------------------------------------------------------------
!    dump wavefunction in r space
!
!      output to file 'eig.dat.ib=mmmm.nnnnnn',
!      where nnnnnn / mmmm are nstep / band index.
!
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
character(*) :: dname_eig
integer :: nstep
integer :: ibstt1, ibstt2
integer :: nplwex, nplw, nband, nbnod1, nbnod2, nbnod, iod(*)
real*8  :: cgjr(nplwex,*)
integer :: mfd2ft(*), ntotfd, nd1vks(*), kfft0
real*8  :: rvol
real*8  :: fft3x(*), fft3y(*)
complex*16 :: fftwork(*)
integer :: ngenh, ijkgd(-ngenh:ngenh)
real*8  :: apk(*), eigr(*)
integer :: nspin
logical :: lspin
integer :: nplw3
integer :: ngboxa(0:nplw3), ngboxb(0:nplw3), ngboxc(0:nplw3)
integer :: kfft1b, kfft2b, kfft3b, kfft0b
integer :: ijkgb(-nplw3:nplw3)
integer :: ntype, natom
integer :: iatoit(natom)
integer :: ioa(*)
logical :: lvand, lvandi(ntype)
real*8  :: xmin, xmax, ymin, ymax, zmin, zmax
logical :: lcompress
integer :: ndigit

!-----declare local variables
integer :: kfft1, kfft2, kfft3, ib, ir, ifd, ift, ibb
integer :: numch2, icount, i1, i2, i3, i123, numch3
real*8  :: sqrvol, ch2mx, plus, ch2mxx
logical :: larea
integer :: meshx1, meshx2
integer :: meshy1, meshy2
integer :: meshz1, meshz2
integer :: iunit
character(50) :: fname
integer :: digit
!real*8,  parameter :: facch2 = 999.d0
!logical, parameter :: lcompress  = .true.
real*8  :: facch2
integer :: ierror


facch2 = 10d0**ndigit - 1d0

kfft1 = nd1vks(1)
kfft2 = nd1vks(2)
kfft3 = nd1vks(3)

!-----choose area
if( xmin > xmax ) then
    meshx1 = 0
    meshx2 = kfft1-1
  else
    meshx1 = nint(dble(kfft1)*xmin)
    meshx2 = nint(dble(kfft1)*xmax)
    if( meshx1 < 0       ) meshx1 = 0
    if( meshx2 > kfft1-1 ) meshx2 = kfft1-1
end if
if( ymin > ymax ) then
    meshy1 = 0
    meshy2 = kfft2-1
  else
    meshy1 = nint(dble(kfft2)*ymin)
    meshy2 = nint(dble(kfft2)*ymax)
    if( meshy1 < 0       ) meshy1 = 0
    if( meshy2 > kfft2-1 ) meshy2 = kfft2-1
end if
if( zmin > zmax ) then
    meshz1 = 0
    meshz2 = kfft3-1
  else
    meshz1 = nint(dble(kfft3)*zmin)
    meshz2 = nint(dble(kfft3)*zmax)
    if( meshz1 < 0       ) meshz1 = 0
    if( meshz2 > kfft3-1 ) meshz2 = kfft3-1
end if
larea = meshx1 /= 0 .or. meshx2 /= kfft1-1 .or.  &
&       meshy1 /= 0 .or. meshy2 /= kfft2-1 .or.  &
&       meshz1 /= 0 .or. meshz2 /= kfft3-1

do ib = 1, nband
!      do ibb = 1, nbnod
!         ib = iod(ibb)
   bandif: if( ib.ge.ibstt1 .and. ib.le.ibstt2 ) then

!    if( lvand ) then
!        !--- contributions from hard part
!        call chgfrq_ibjb( nfile, myid, nodes,  &
!& ib, ib, apk, nplw3, ngboxa, ngboxb, ngboxc,  &
!& kfft1b, kfft2b, kfft3b, kfft0b, fft3x, fft3y, fftwork, ijkgb,  &
!& eigr, ntotfd, rvol,  &
!& ntype, natom, iatoit, ioa, lvandi )
!
!        !--- convert FD -> FFT meshes ---
!        do ir = 1, kfft0
!           eigr(ir) = 0.d0
!        end do
!        do ifd = 1, ntotfd
!           ift = mfd2ft(ifd)
!           eigr(ift) = apk(ifd)
!        end do
!    end if


   ibb = ib - nbnod1 + 1
   nodeif: if( ib.ge.nbnod1 .and. ibb.le.nbnod ) then

       call wvg2r( nfile, myid, nodes,  &
& cgjr(1,ibb), nplwex, nplw, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, .false. )

!--- sqrvol is a factor to get the normalized w.f. in the r space
       sqrvol = 1.d0/sqrt(rvol)
       do ir = 1, kfft0
          fft3x(ir) = fft3x(ir) * sqrvol
       end do

!       if( lvand ) then
           !--- charge density for USPP
!           do ir = 1, kfft0
!              fft3x(ir) = fft3x(ir) * fft3x(ir) + eigr(ir)
!           end do
!       end if

       digit = log(dble(nband)+0.5d0)/log(10.d0) + 1.d0
       fname = dname_eig
      !!! fname = ( fname(1:len_trim(fname)) // '.ib=' )
       fname = ( fname(1:len_trim(fname)) // '.' )
       call get_dumpfname( fname, ib, digit )
       if( lspin ) then
           if( nspin.eq.1 ) then
               fname = ( fname(1:len_trim(fname)) // 'u' )
             else
               fname = ( fname(1:len_trim(fname)) // 'd' )
           end if
       end if

       fname = ( fname(1:len_trim(fname)) // '.' )
       call get_dumpfname( fname, nstep, 6 )


       call allocate_unit_number( iunit )
       open( iunit, file=fname, status='unknown',  &
&                   action='write', iostat=ierror )
       errif: if( ierror == 0 ) then

       !-----dump data
       ch2mx = 0.d0
       plus = 0
       do ir = 1, kfft0
          ch2mx = max( ch2mx, abs(fft3x(ir)) )
          plus = plus + sign(1.d0, fft3x(ir))
       end do
       if( plus < 0.d0 ) ch2mx = -ch2mx


       modeif: if( lcompress ) then

       !--- compress format
       if( .not.larea ) then
           write(iunit,'(3i5)') 0, 0, 0
           write(iunit,'(3i5)') kfft1, kfft2, kfft3
         else
           write(iunit,'(3i5)') -1, -1, -1
           write(iunit,'(3i5)') 0, 0, 0
           write(iunit,'(3i5)') kfft1, kfft2, kfft3
           write(iunit,'(6i5)') meshx1, meshy1, meshz1,  &
&                 meshx2-meshx1+1, meshy2-meshy1+1, meshz2-meshz1+1
       end if
           write(iunit,'(2es15.6)') ch2mx / facch2

           ch2mxx = facch2/ch2mx
!           numch2 = nint( fft3x(1)*ch2mxx )
           i3 = meshz1
           i2 = meshy1
           i1 = meshx1
           i123=1+i1+(i2+i3*kfft2)*kfft1
           numch2 = nint( fft3x(i123)*ch2mxx )
           icount = 0
!           do i3=0, kfft3-1
!           do i2=0, kfft2-1
!           do i1=0, kfft1-1
           do i3 = meshz1, meshz2
           do i2 = meshy1, meshy2
           do i1 = meshx1, meshx2

              i123=1+i1+(i2+i3*kfft2)*kfft1
              numch3 = nint( fft3x(i123)*ch2mxx )

              if( numch3 /= numch2 ) then

                 !--- output icount
                 write(iunit,'(i0)') icount

                 !--- output numch2
                 write(iunit,'(i0)') numch2

                 icount = 1
                 numch2 = numch3

              else

                  icount = icount + 1

              end if

           end do
           end do
           end do

                 !--- output icount
                 write(iunit,'(i0)') icount

                 !--- output numch2
                 write(iunit,'(i0)') numch2

       else modeif

       !--- normal format
!!!             write(iunit,'(i7)') nstep
       if( .not.larea ) then
           write(iunit,'(3i5)') kfft1, kfft2, kfft3
         else
           write(iunit,'(3i5)') -1, -1, -1
           write(iunit,'(3i5)') kfft1, kfft2, kfft3
           write(iunit,'(6i5)') meshx1, meshy1, meshz1,  &
&                 meshx2-meshx1+1, meshy2-meshy1+1, meshz2-meshz1+1
       end if
       write(iunit,'(2es15.6)') ch2mx / facch2
!             write(iunit,'(2es15.6)') ch2mx, facch2

       ch2mxx = facch2/ch2mx
!       do i3=0, kfft3-1
!       do i2=0, kfft2-1
!       do i1=0, kfft1-1
       do i3 = meshz1, meshz2
       do i2 = meshy1, meshy2
       do i1 = meshx1, meshx2
          i123=1+i1+(i2+i3*kfft2)*kfft1
          numch2 = nint( fft3x(i123)*ch2mxx )
          write(iunit,'(i0)') numch2
       end do
       end do
       end do

       end if modeif

       close(iunit)
       call deallocate_unit_number( iunit )

       end if errif

   end if nodeif
   end if bandif
end do


return
end




subroutine get_dumpfname( fname, nstep, digit )
!-----------------------------------------------------------------------
implicit none
character(*) :: fname
integer :: nstep
integer :: digit

integer :: iu
integer :: ibuf, i, itenth
character(1), dimension(0:9) :: num =  &
& (/ '0','1','2','3','4','5','6','7','8','9' /)


ibuf = nstep
do i = digit, 1, -1
   itenth = 10**(i-1)
   iu     = ibuf/itenth
   ibuf   = ibuf - iu*itenth

   fname = ( fname(1:len_trim(fname)) // num(iu) )
end do


return
end




subroutine get_frmt( frmt, digit )
!-----------------------------------------------------------------------
! Chenge frmt = '(i7******' to '(i[digit]******'
!-----------------------------------------------------------------------
implicit none
character(*) :: frmt
integer :: digit

character(9) :: c

!write(c,*) digit
write(c(1:9),'(i9)') digit
frmt = '(i' // trim(adjustl(c)) // frmt(4:len_trim(frmt))

return
end




subroutine remd_wrgeom( nfile, myid, nodes, iogpsz,  &
& nstep, ifmd, x, v, is, n, dt, ntype, a, h, acon, lfixion,  &
& sxog, syog, szog, pit1, fack, volume,  &
& lmdout, locoor, lovelo, loforc, ioskip, ioskipcoor, ioskipvelo, ioskipforc,  &
& ncbonds, lstress, nskip_stress, nskip_stress_out, latomic, nskip_atomic,  &
& wepot, wstrs, wstrsk, nmdmax__, ntot,  &
& lrestrict_area, xio_min, xio_max, yio_min, yio_max, zio_min, zio_max,  &
& ierror )
!----------------------------------------------------------------------c
!     output atomic configurations etc.
!----------------------------------------------------------------------c
implicit none
integer :: nfile(*), myid, nodes
integer :: iogpsz
integer :: nstep
integer :: ifmd
integer :: n
real*8,  dimension(3,n) :: x
real*8,  dimension(3,n) :: v
integer, dimension(n) :: is
real*8  :: dt
integer :: ntype
real*8,  dimension(3*n) :: a
real*8,  dimension(3,3,0:1) :: h
real*8,  dimension(ntype) :: fack
real*8  :: volume
real*8,  dimension(ntype) :: acon
logical, dimension(ntype) :: lfixion
real*8  :: sxog, syog, szog
real*8,  dimension(3,3) :: pit1
logical :: lmdout, locoor, lovelo, loforc
integer :: ioskip, ioskipcoor, ioskipvelo, ioskipforc
integer :: ierror
integer :: ncbonds
!real*8,  dimension(ncbonds) :: cblambda, cblambdav    ! possibly not allocated
logical :: lstress, latomic
integer :: nskip_stress, nskip_stress_out, nskip_atomic
integer :: nmdmax__
real*8  :: wepot(nmdmax__), wstrs(6,nmdmax__), wstrsk(6,nmdmax__)
integer, dimension(0:ntype) :: ntot
logical :: lrestrict_area
real*8  :: xio_min, xio_max, yio_min, yio_max, zio_min, zio_max

!-----declare local variables
logical :: loutedg
real*8  :: dltca, dltcb, dltcc, angalf, angbet, anggam
real*8  :: atconfx
integer :: i, i3, ix, iy, iz, digit
real*8, parameter :: prau = 2.94216d4
real*8, dimension(3,3) :: rbarec = 0.d0
logical :: lcstress, lcatomic
integer, parameter :: ireci = 36, irecd = 3
integer :: isr(ireci), ii, nrec, i2, it
real*8  :: xr(3,irecd)
real*8  :: vgt(3,3), veig(3), veigv(3,3)
logical :: larea
character(50) :: frmt
!---for group I/O
integer :: nfiles, n1, n1alloc, isnd
integer, allocatable :: ibuf(:)
real(8), allocatable :: abuf(:)
save rbarec


digit = log(dble(nstep)+0.5d0)/log(10.d0) + 1.d0

!-----check if stresses are calculated
lcstress = lstress .and. mod(nstep,nskip_stress) == 0
!-----check if atomic stresses & energies are calculated
lcatomic = latomic .and. mod(nstep,nskip_atomic) == 0

!-----internal stress tensor
if( lcstress ) then
    if( ifmd >= 2 ) then
        !-----Kinetic energy contribution in MD
        call get_ekinp( nfile, myid, nodes,  &
& pit1, v, is, n, ntype, h, fack, volume,  &
& lcatomic, wstrsk, nmdmax__, ntot )
    end if
end if

if( myid == 0 ) then

    !-----MD cell vectors
    call caledge( loutedg, h(1,1,0), rbarec,  &
&             dltca, dltcb, dltcc, angalf, angbet, anggam )
    if( loutedg ) then
        frmt = '(i7,3es14.7,3f11.6)'
        if( digit > 7 ) call get_frmt( frmt, digit )
        write(nfile(4),frmt) nstep,  &
&             dltca, dltcb, dltcc, angalf, angbet, anggam
        frmt = '(i7,9es15.7)'
        if( digit > 7 ) call get_frmt( frmt, digit )
        write(nfile(9),frmt) nstep, h(1:3,1:3,0)
    end if

!    if( lcstress ) then
    if( nskip_stress_out > 0 ) then
    if( lcstress .and. mod(nstep,nskip_stress_out) == 0 ) then
        !-----internal stress tensor
        frmt = '(i7,10es16.8)'
        if( digit > 7 ) call get_frmt( frmt, digit )
        write(nfile(16),frmt) nstep,  &
&         pit1(1,1)*prau, pit1(2,2)*prau, pit1(3,3)*prau,  &
&         pit1(2,3)*prau, pit1(3,1)*prau, pit1(1,2)*prau

        !----- diagonalize stress tensor
        vgt(1:3,1:3) = pit1(1:3,1:3)
        call eigen3(vgt,veig,veigv)
        frmt = '(i7,3es16.8,9es12.4)'
        if( digit > 7 ) call get_frmt( frmt, digit )
        write(nfile(22),frmt) nstep,  &
&         veig(1)*prau, veig(2)*prau, veig(3)*prau, veigv(1:3,1:3)
    end if
    end if

!    if( ifmd >= 1 .and. ncbonds > 0 ) then
!        call out_cblambda( nfile, myid, nodes, nstep, digit )
!    end if

    if( ifmd == 10 ) then
        call out_hugoniot( nfile, myid, nodes, nstep, digit )
    end if

end if


if( lmdout ) then
!if( mod(nstep,ioskip) == 0 .or. ifmd == 0 ) then
if( mod(nstep,1) == 0 .or. ifmd == 0 ) then

    if(    locoor .and. mod(nstep,ioskipcoor) == 0  &
&     .or. lovelo .and. mod(nstep,ioskipvelo) == 0  &
&     .or. loforc .and. mod(nstep,ioskipforc) == 0 ) then

    ioif: if( mod(myid,iogpsz) == 0 ) then
        nfiles=min(iogpsz,nodes-myid)

        !--- write my node data
        !-----Atomic species
!        if( .not.lrestrict_area ) then
            frmt = '(i7,i7)'
            if( digit > 7 ) call get_frmt( frmt, digit )
            write(nfile(5),frmt) nstep, n
            write(nfile(5),'(36i2)') ( is(i), i = 1, n )
!        else
!            !---precalculation of the number of atoms
!            nrec = 0
!            do i = 1, n
!               if( larea( x(1,i), sxog, syog, szog,  &
!& xio_min, xio_max, yio_min, yio_max, zio_min, zio_max ) ) then
!                   nrec = nrec + 1
!               end if
!            end do
!
!            frmt = '(i7,i7)'
!            if( digit > 7 ) call get_frmt( frmt, digit )
!            write(nfile(5),frmt) nstep, nrec
!            ii = 0; i2 = 0
!            do i = 1, n
!               if( larea( x(1,i), sxog, syog, szog,  &
!& xio_min, xio_max, yio_min, yio_max, zio_min, zio_max ) ) then
!                   ii = ii + 1; i2 = i2 + 1
!                   isr(ii) = is(i)
!                   if( ii == ireci .or. i2 == nrec ) then
!                       write(nfile(5),'(36i2)') isr(1:ii)
!                       ii = 0
!                   end if
!               end if
!            end do
!        end if


        n1alloc = 0
        do ii=1, nfiles-1
           call cirecvs(myid+ii,n1,1,myid+ii,0)

           if( n1 > n1alloc ) then
               if( allocated(ibuf) ) deallocate(ibuf)
               n1alloc = n1
               allocate(ibuf(n1alloc))
           end if

           call cirecvs(myid+ii,ibuf,n1,myid+ii,0)

           write(nfile(5),'(2i7)') n1
           write(nfile(5),'(36i2)') ( ibuf(i), i = 1, n1 )

        end do
        if( allocated(ibuf) ) deallocate(ibuf)

    else ioif

        isnd=(myid/iogpsz)*iogpsz
!        if( .not.lrestrict_area ) then
            call cisend(myid,n,1,isnd,0)
            call cisend(myid,is,n,isnd,0)
!        else
!            !---precalculation of the number of atoms
!            nrec = 0
!            do i = 1, n
!               if( larea( x(1,i), sxog, syog, szog,  &
!& xio_min, xio_max, yio_min, yio_max, zio_min, zio_max ) ) then
!                   nrec = nrec + 1
!               end if
!            end do
!            call cisend(myid,nrec,1,isnd,0)
!
!            allocate(ibuf(nrec))
!
!            ii = 0
!            do i = 1, n
!               if( larea( x(1,i), sxog, syog, szog,  &
!& xio_min, xio_max, yio_min, yio_max, zio_min, zio_max ) ) then
!                   ii = ii + 1
!                   ibuf(ii) = is(i)
!               end if
!            end do
!            call cisend(myid,ibuf,nrec,isnd,0)
!
!            deallocate(ibuf)
!        end if

    end if ioif

    end if


    !-----atomic scaled coordinates
    if( locoor .and. mod(nstep,ioskipcoor) == 0 ) then

        atconfx = 0.1d0

        ioif2: if( mod(myid,iogpsz) == 0 ) then
            nfiles=min(iogpsz,nodes-myid)

            !--- write my node data
!            if( .not.lrestrict_area ) then

                frmt = '(i7,9i7)'
                if( digit > 7 ) call get_frmt( frmt, digit )
                write(nfile(6),frmt) nstep, n
                write(nfile(6),'(es14.7)') atconfx
                write(nfile(6),'(9f8.5)')  &
&           ( min( (x(1,i)+sxog - floor(x(1,i)+sxog))/atconfx, 9.99999d0 ),  &
&             min( (x(2,i)+syog - floor(x(2,i)+syog))/atconfx, 9.99999d0 ),  &
&             min( (x(3,i)+szog - floor(x(3,i)+szog))/atconfx, 9.99999d0 ), i = 1, n )

!            else
!
!                frmt = '(i7,9i7)'
!                if( digit > 7 ) call get_frmt( frmt, digit )
!                write(nfile(6),frmt) nstep, nrec
!                write(nfile(6),'(es14.7)') atconfx
!                ii = 0; i2 = 0
!                do i = 1, n
!                   if( larea( x(1,i), sxog, syog, szog,  &
!& xio_min, xio_max, yio_min, yio_max, zio_min, zio_max ) ) then
!                       ii = ii + 1; i2 = i2 + 1
!                       xr(1,ii) = min( (x(1,i)+sxog - floor(x(1,i)+sxog))/atconfx, 9.99999d0 )
!                       xr(2,ii) = min( (x(2,i)+syog - floor(x(2,i)+syog))/atconfx, 9.99999d0 )
!                       xr(3,ii) = min( (x(3,i)+szog - floor(x(3,i)+szog))/atconfx, 9.99999d0 )
!                       if( ii == irecd .or. i2 == nrec ) then
!                           write(nfile(6),'(9f8.5)') ( xr(1:3,i3), i3 = 1, ii )
!                           ii = 0
!                       end if
!                   end if
!                end do
!
!            end if


            n1alloc = 0
            do ii=1, nfiles-1
               call cirecvs(myid+ii,n1,1,myid+ii,0)

               if( n1*3 > n1alloc ) then
                   if( allocated(abuf) ) deallocate(abuf)
                   n1alloc = n1*3
                   allocate(abuf(n1alloc))
               end if

               call cdrecvs(myid+ii,abuf,n1*3,myid+ii,0)

               write(nfile(6),'(10i7)') n1
               write(nfile(6),'(9f8.5)') abuf(1:3*n1)

            end do
            if( allocated(abuf) ) deallocate(abuf)

        else ioif2

            isnd=(myid/iogpsz)*iogpsz

!            if( .not.lrestrict_area ) then
                call cisend(myid,n,1,isnd,0)

                allocate(abuf(3*n))

                do i = 1, n
                   abuf(3*i-2) = min( (x(1,i)+sxog - floor(x(1,i)+sxog))/atconfx, 9.99999d0 )
                   abuf(3*i-1) = min( (x(2,i)+syog - floor(x(2,i)+syog))/atconfx, 9.99999d0 )
                   abuf(3*i  ) = min( (x(3,i)+szog - floor(x(3,i)+szog))/atconfx, 9.99999d0 )
                end do

                call cdsend(myid,abuf,n*3,isnd,0)

                deallocate(abuf)

!            else
!
!                call cisend(myid,nrec,1,isnd,0)
!
!                allocate(abuf(3*nrec))
!
!                ii = 0
!                do i = 1, n
!                   if( larea( x(1,i), sxog, syog, szog,  &
!& xio_min, xio_max, yio_min, yio_max, zio_min, zio_max ) ) then
!                       ii = ii + 1
!                       abuf(3*ii-2) = min( (x(1,i)+sxog - floor(x(1,i)+sxog))/atconfx, 9.99999d0 )
!                       abuf(3*ii-1) = min( (x(2,i)+syog - floor(x(2,i)+syog))/atconfx, 9.99999d0 )
!                       abuf(3*ii  ) = min( (x(3,i)+szog - floor(x(3,i)+szog))/atconfx, 9.99999d0 )
!                   end if
!                end do
!
!                call cdsend(myid,abuf,nrec*3,isnd,0)
!
!                deallocate(abuf)
!
!            end if

        end if ioif2

    end if


    !--- atomic scaled velocities ---
    if( ifmd >= 2 ) then
        if( lovelo .and. mod(nstep,ioskipvelo) == 0 ) then

            atconfx = 0.d0
            do i = 1, n
               atconfx = max( atconfx, abs(v(1,i)), abs(v(2,i)), abs(v(3,i)) )
            end do
            atconfx = atconfx / 9.999d0

            ioif3: if( mod(myid,iogpsz) == 0 ) then
                nfiles=min(iogpsz,nodes-myid)

                !--- write my node data
!                if( .not.lrestrict_area ) then

                    frmt = '(i7,9i7)'
                    if( digit > 7 ) call get_frmt( frmt, digit )
                    write(nfile(7),frmt) nstep, n
                    write(nfile(7),'(es14.7)') atconfx/dt
                    !--- avoid floating divide by zero
                    if( abs(atconfx) < 1.d-30 ) atconfx = 1.d0
                    write(nfile(7),'(9f8.5)') (( v(ix,i)/atconfx, ix = 1, 3 ), i = 1, n )

!                else
!
!                    frmt = '(i7,9i7)'
!                    if( digit > 7 ) call get_frmt( frmt, digit )
!                    write(nfile(7),frmt) nstep, nrec
!                    write(nfile(7),'(es14.7)') atconfx/dt
!                    !--- avoid floating divide by zero
!                    if( abs(atconfx) < 1.d-30 ) atconfx = 1.d0
!                    ii = 0; i2 = 0
!                    do i = 1, n
!                       if( larea( x(1,i), sxog, syog, szog,  &
!& xio_min, xio_max, yio_min, yio_max, zio_min, zio_max ) ) then
!                           ii = ii + 1; i2 = i2 + 1
!                           xr(1:3,ii) = v(1:3,i)/atconfx
!                           if( ii == irecd .or. i2 == nrec ) then
!                               write(nfile(7),'(9f8.5)') ( xr(1:3,i3), i3 = 1, ii )
!                               ii = 0
!                           end if
!                       end if
!                    end do
!
!                end if


                n1alloc = 0
                do ii=1, nfiles-1
                   call cirecvs(myid+ii,n1,1,myid+ii,0)

                   if( n1*3+1 > n1alloc ) then
                       if( allocated(abuf) ) deallocate(abuf)
                       n1alloc = n1*3+1
                       allocate(abuf(n1alloc))
                   end if

                   call cdrecvs(myid+ii,abuf,3*n1+1,myid+ii,0)

                   write(nfile(7),'(10i7)') n1
                   write(nfile(7),'(es14.7)') abuf(3*n1+1)/dt
                   write(nfile(7),'(9f8.5)') abuf(1:3*n1)

                end do
                if( allocated(abuf) ) deallocate(abuf)

            else ioif3

                isnd=(myid/iogpsz)*iogpsz

!                if( .not.lrestrict_area ) then

                    call cisend(myid,n,1,isnd,0)

                    allocate(abuf(3*n+1))

                    abuf(3*n+1) = atconfx
                    !--- avoid floating divide by zero
                    if( abs(atconfx) < 1.d-30 ) atconfx = 1.d0
                    do i = 1, n
                       abuf(3*i-2) = v(1,i)/atconfx
                       abuf(3*i-1) = v(2,i)/atconfx
                       abuf(3*i  ) = v(3,i)/atconfx
                    end do

                    call cdsend(myid,abuf,3*n+1,isnd,0)

                    deallocate(abuf)

!                else
!
!                    call cisend(myid,nrec,1,isnd,0)
!
!                    allocate(abuf(3*nrec+1))
!
!                    abuf(3*nrec+1) = atconfx
!                    !--- avoid floating divide by zero
!                    if( abs(atconfx) < 1.d-30 ) atconfx = 1.d0
!                    ii = 0
!                    do i = 1, n
!                       if( larea( x(1,i), sxog, syog, szog,  &
!& xio_min, xio_max, yio_min, yio_max, zio_min, zio_max ) ) then
!                           ii = ii + 1; i2 = i2 + 1
!                           abuf(3*ii-2) = v(1,i)/atconfx
!                           abuf(3*ii-1) = v(2,i)/atconfx
!                           abuf(3*ii  ) = v(3,i)/atconfx
!                       end if
!                    end do
!
!                    call cdsend(myid,abuf,3*nrec+1,isnd,0)
!
!                    deallocate(abuf)
!
!                end if

            end if ioif3

        end if
    end if


    !--- atomic scaled forces ---
    if( loforc .and. mod(nstep,ioskipforc) == 0 ) then

        atconfx = 0.d0
!        do i3  = 1, 3*n
!           i=(i3+2)/3
!           atconfx = max( atconfx, abs( a(i3)/acon(is(i)) ) )
!        end do
        do i  = 1, n
           it = is(i)
           if( .not.lfixion(it) ) then
               atconfx = max( atconfx, abs( a(i*3-2)/acon(it) ),  &
&         abs( a(i*3-1)/acon(it) ), abs( a(i*3  )/acon(it) ) )
           end if
        end do
        atconfx = atconfx / 9.999d0

        ioif4: if( mod(myid,iogpsz) == 0 ) then
            nfiles=min(iogpsz,nodes-myid)

            !--- write my node data
!            if( .not.lrestrict_area ) then

                frmt = '(i7,9i7)'
                if( digit > 7 ) call get_frmt( frmt, digit )
                write(nfile(8),frmt) nstep, n
                write(nfile(8),'(es14.7)') atconfx
                !--- avoid floating divide by zero
                if( abs(atconfx) < 1.d-30 ) atconfx = 1.d0
!                write(nfile(8),'(9f8.5)')  &
!&                ( a(i3)/acon(is((i3+2)/3))/atconfx, i3 = 1, 3*n )
                ii = 0; i2 = 0
                do i = 1, n
                   ii = ii + 1; i2 = i2 + 1
                   it = is(i)
                   if( .not.lfixion(it) ) then
                       xr(1,ii) = a(3*i-2)/acon(it)/atconfx
                       xr(2,ii) = a(3*i-1)/acon(it)/atconfx
                       xr(3,ii) = a(3*i  )/acon(it)/atconfx
                   else
                       xr(1:3,ii) = 0.d0
                   end if
                   if( ii == irecd .or. i2 == n ) then
                       write(nfile(8),'(9f8.5)') ( xr(1:3,i3), i3 = 1, ii )
                       ii = 0
                   end if
                end do

!            else
!
!                frmt = '(i7,9i7)'
!                if( digit > 7 ) call get_frmt( frmt, digit )
!                write(nfile(8),frmt) nstep, nrec
!                write(nfile(8),'(es14.7)') atconfx
!                !--- avoid floating divide by zero
!                if( abs(atconfx) < 1.d-30 ) atconfx = 1.d0
!                ii = 0; i2 = 0
!                do i = 1, n
!                   if( larea( x(1,i), sxog, syog, szog,  &
!& xio_min, xio_max, yio_min, yio_max, zio_min, zio_max ) ) then
!                       ii = ii + 1; i2 = i2 + 1
!                       it = is(i)
!                       if( .not.lfixion(it) ) then
!                           xr(1,ii) = a(3*i-2)/acon(it)/atconfx
!                           xr(2,ii) = a(3*i-1)/acon(it)/atconfx
!                           xr(3,ii) = a(3*i  )/acon(it)/atconfx
!                       else
!                           xr(1:3,ii) = 0.d0
!                       end if
!                       if( ii == irecd .or. i2 == nrec ) then
!                           write(nfile(8),'(9f8.5)') ( xr(1:3,i3), i3 = 1, ii )
!                           ii = 0
!                       end if
!                   end if
!                end do
!
!            end if


            n1alloc = 0
            do ii=1, nfiles-1
               call cirecvs(myid+ii,n1,1,myid+ii,0)

               if( n1*3+1 > n1alloc ) then
                   if( allocated(abuf) ) deallocate(abuf)
                   n1alloc = n1*3+1
                   allocate(abuf(n1alloc))
               end if

               call cdrecvs(myid+ii,abuf,3*n1+1,myid+ii,0)

               write(nfile(8),'(10i7)') n1
               write(nfile(8),'(es14.7)') abuf(3*n1+1)
               write(nfile(8),'(9f8.5)') abuf(1:3*n1)

            end do
            if( allocated(abuf) ) deallocate(abuf)

        else ioif4

            isnd=(myid/iogpsz)*iogpsz

!            if( .not.lrestrict_area ) then

                call cisend(myid,n,1,isnd,0)

                allocate(abuf(3*n+1))

                abuf(3*n+1) = atconfx
                !--- avoid floating divide by zero
                if( abs(atconfx) < 1.d-30 ) atconfx = 1.d0
                do i = 1, n
                   it = is(i)
                   if( .not.lfixion(it) ) then
                       abuf(3*i-2:3*i) = a(3*i-2:3*i)/acon(it)/atconfx
                   else
                       abuf(3*i-2:3*i) = 0.d0
                   end if
                end do

                call cdsend(myid,abuf,3*n+1,isnd,0)

                deallocate(abuf)

!            else
!
!                call cisend(myid,nrec,1,isnd,0)
!
!                allocate(abuf(3*nrec+1))
!
!                abuf(3*nrec+1) = atconfx
!                !--- avoid floating divide by zero
!                if( abs(atconfx) < 1.d-30 ) atconfx = 1.d0
!                ii = 0
!                do i = 1, n
!                   if( larea( x(1,i), sxog, syog, szog,  &
!& xio_min, xio_max, yio_min, yio_max, zio_min, zio_max ) ) then
!                       ii = ii + 1
!                       it = is(i)
!                       if( .not.lfixion(it) ) then
!                           abuf(3*ii-2:3*ii) = a(3*i-2:3*i)/acon(it)/atconfx
!                       else
!                           abuf(3*ii-2:3*ii) = 0.d0
!                       end if
!                   end if
!                end do
!
!                call cdsend(myid,abuf,3*nrec+1,isnd,0)
!
!                deallocate(abuf)
!
!            end if

        end if ioif4

    end if

end if
end if


!if( lcatomic ) then
!    if( .not.lrestrict_area ) then
!
!        !--- atomic energies
!        call out_f85( nfile, myid, nodes, iogpsz,  &
!& nfile(19), nstep, digit, n, wepot )
!
!        if( lcstress ) then
!            !--- atomic stresses from potential energy
!            call out_f85( nfile, myid, nodes, iogpsz,  &
!& nfile(20), nstep, digit, 6*n, wstrs )
!            if( ifmd >= 2 ) then
!                !--- atomic stresses from kinetic energy
!                call out_f85( nfile, myid, nodes, iogpsz,  &
!& nfile(21), nstep, digit, 6*n, wstrsk )
!            end if
!        end if
!
!    else
!
!        !--- atomic energies
!        call out_f85_area( nfile, myid, nodes, iogpsz,  &
!& nfile(19), nstep, digit, 1, n, wepot,  &
!& nrec, x, sxog, syog, szog,  &
!& xio_min, xio_max, yio_min, yio_max, zio_min, zio_max )
!
!        if( lcstress ) then
!            !--- atomic stresses from potential energy
!            call out_f85_area( nfile, myid, nodes, iogpsz,  &
!& nfile(20), nstep, digit, 6, n, wstrs,  &
!& nrec, x, sxog, syog, szog,  &
!& xio_min, xio_max, yio_min, yio_max, zio_min, zio_max )
!            if( ifmd >= 2 ) then
!                !--- atomic stresses from kinetic energy
!                call out_f85_area( nfile, myid, nodes, iogpsz,  &
!& nfile(21), nstep, digit, 6, n, wstrsk,  &
!& nrec, x, sxog, syog, szog,  &
!& xio_min, xio_max, yio_min, yio_max, zio_min, zio_max )
!            end if
!        end if
!
!    end if
!end if


return
end




subroutine caledge( loutedg, rba, rbarec,  &
&                   dltca, dltcb, dltcc, angalf, angbet, anggam )
!-----------------------------------------------------------------------
!     calculate supercell edges
!
! (output)
!    dltca,  dltcb,  dltcc  : lengths of edges
!    angalf, angbet, anggam : angles between edges
!-----------------------------------------------------------------------
implicit real*8(a-h,o-z)
real*8, parameter :: angpi = 57.29577951308d0
logical :: loutedg
real*8, dimension(3,3) :: rba
real*8, dimension(3,3) :: rbarec


loutedg = .false.
do i = 1, 3
do j = 1, 3
   loutedg = loutedg .or. abs(rbarec(j,i)-rba(j,i)).gt.1.d-12
end do
end do

dltca =sqrt(rba(1,1)*rba(1,1)+rba(2,1)*rba(2,1)+rba(3,1)*rba(3,1))
dltcb =sqrt(rba(1,2)*rba(1,2)+rba(2,2)*rba(2,2)+rba(3,2)*rba(3,2))
dltcc =sqrt(rba(1,3)*rba(1,3)+rba(2,3)*rba(2,3)+rba(3,3)*rba(3,3))

cosbc = rba(1,2)*rba(1,3) + rba(2,2)*rba(2,3) + rba(3,2)*rba(3,3)
cosbc = cosbc / ( dltcb*dltcc )
angalf = acos(cosbc) * angpi
cosbc = rba(1,3)*rba(1,1) + rba(2,3)*rba(2,1) + rba(3,3)*rba(3,1)
cosbc = cosbc / ( dltcc*dltca )
angbet = acos(cosbc) * angpi
cosbc = rba(1,1)*rba(1,2) + rba(2,1)*rba(2,2) + rba(3,1)*rba(3,2)
cosbc = cosbc / ( dltca*dltcb )
anggam = acos(cosbc) * angpi

do i = 1, 3
do j = 1, 3
  rbarec(j,i) = rba(j,i)
end do
end do


return
end
