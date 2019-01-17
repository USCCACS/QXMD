



subroutine stress( nfile, myid, nodes,  &
& str, strloc, strnlc, strclm, strgga, strpcc, ltimecnt, nstep,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, occ, rvol, iod,  &
& rhgr, nplw5ex, nplw5, recnrm, gx, gy, gz,  &
& ehar, vexc, strexc,  &
& rdel, rdelg, rdelv, hcell, hci, lorthrhmbc, lvacuum,  &
& ntype, natom, nhk1, nhk2, ratm, mx1, encl, t_comm,  &
& lpcc, lpcc_r, lpcc_g, lpcci, rpcc, lpking,  &
& mshxyz, mulx1, mulx2, muly1, muly2, mulz1, mulz2,  &
& mshnx, mshny, mshnz, mshnod, mshnew,  &
& mshx1, mshy1, mshz1, mshx, mshy, mshz,  &
& lspin, nspnmx, gdcrsv, nspnod, lcgjsv, wegud, ldftd,  &
& strhfs, exhfsr, jhybrid, ncscale, strefd, lefield )
!-----------------------------------------------------------------------
!    internal stress tensor from kinetic E. & Hartree E. & Exc
!-----------------------------------------------------------------------
use outfile
implicit none
integer :: nfile(*), myid, nodes
real*8,  dimension(3,3) :: str
real*8,  dimension(3,3) :: strloc, strnlc, strclm, strgga, strpcc
logical :: ltimecnt
integer :: nstep
integer :: nplwex, nplw, nband, nbnod1, nbnod2, nbnod
real*8  :: cgjr(nplwex,*)
real*8  :: occ(*)
real*8  :: rvol
integer :: iod(*)
integer :: nplw5ex, nplw5
real*8  :: rhgr(nplw5ex)
real*8  :: recnrm(0:nplw5)
real*8  :: gx(0:nplw5), gy(0:nplw5), gz(0:nplw5)
real*8  :: ehar, vexc(*), strexc
real*8  :: rdel(3), rdelg(3,3), rdelv
real*8  :: hcell(3,3), hci(3,3)
logical :: lorthrhmbc
logical :: lvacuum(3)
integer :: ntype, natom, nhk1(*), nhk2(*)
real*8  :: ratm(3,*)
integer :: mx1
real*8  :: encl
real*8  :: t_comm
logical :: lpcc, lpcc_r, lpcc_g, lpcci(*), lpking(*)
real*8  :: rpcc(*)
integer :: mulx1, mulx2, muly1, muly2, mulz1, mulz2
integer :: mshxyz(mulx1:mulx2,muly1:muly2,mulz1:mulz2)
integer :: mshnod, mshnew
integer :: mshnx(mshnod), mshny(mshnod), mshnz(mshnod)
integer :: mshx1, mshy1, mshz1, mshx, mshy, mshz
logical :: lspin, lcgjsv
integer :: nspnmx, nspnod
real*8  :: gdcrsv(*), wegud(*)
logical :: ldftd
real*8  :: strhfs(3,3), exhfsr
integer :: jhybrid, ncscale
real*8  :: strefd(3,3)
logical :: lefield

!-----declare local variables
integer :: i, j, nf
real*8, parameter :: prau = 1.47108d4
real*8  :: strparts(3,3)


t_comm = 0.d0
str(1:3,1:3) = 0.d0


!--- Stress from kinetic E.
call strkine( nfile, myid, nodes,  &
& strparts, cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod,  &
& rvol, iod, occ, gx, gy, gz,  &
& lspin, nspnmx, gdcrsv, nspnod, lcgjsv, wegud, ncscale )

if( ltimecnt ) then
    do nf = 1, 2
    if( loutfile(nf) ) then
       write(nfile(nf),*) ' '
       write(nfile(nf),*) 'Stress by kinetic E. in GPa'
       do j = 1, 3
          write(nfile(nf),721) ( strparts(j,i)*prau, i = 1, 3 )
       end do
    end if
    end do
end if
721   format(5x,3e16.8)

str(1:3,1:3) = str(1:3,1:3) + strparts(1:3,1:3)


!--- Stress from Hartree E.
call strhart( nfile, myid, nodes,  &
& strparts, rhgr, nplw5ex, nplw5, recnrm, gx, gy, gz )

strparts(1,1) = strparts(1,1) + ehar/rvol
strparts(2,2) = strparts(2,2) + ehar/rvol
strparts(3,3) = strparts(3,3) + ehar/rvol

if( ltimecnt ) then
    do nf = 1, 2
    if( loutfile(nf) ) then
       write(nfile(nf),*) ' '
       write(nfile(nf),*) 'Stress by Hartree E. in GPa'
       do j = 1, 3
          write(nfile(nf),721) ( strparts(j,i)*prau, i = 1, 3 )
       end do
    end if
    end do
end if

str(1:3,1:3) = str(1:3,1:3) + strparts(1:3,1:3)


!--- Stress from exchange-correlation E.
strparts(1:3,1:3) = 0.d0
do i = 1, 3
   strparts(i,i) = strparts(i,i) + strexc
end do
strparts(1:3,1:3) = strparts(1:3,1:3) + strgga(1:3,1:3)
strparts(1:3,1:3) = strparts(1:3,1:3) * (-1.d0) / rvol

    !---symmetry operations on stress
!    call symk_stress( nfile, myid, nodes, strparts )

if( ltimecnt ) then
    do nf = 1, 2
    if( loutfile(nf) ) then
       write(nfile(nf),*) ' '
       write(nfile(nf),*) 'Stress by Exc in GPa'
       do j = 1, 3
          write(nfile(nf),721) ( strparts(j,i)*prau, i = 1, 3 )
       end do
    end if
    end do
end if

str(1:3,1:3) = str(1:3,1:3) + strparts(1:3,1:3)


!--- Stress from exchange-correlation E.
if( lpcc ) then
    strparts(1:3,1:3) = 0.d0
    if( lpcc_r ) then
        call pccstr( nfile, myid, nodes,  &
& vexc, strparts, rdel, rdelg, rdelv, hcell, hci, lorthrhmbc,  &
& lvacuum, ntype, natom, nhk1, nhk2, ratm, .false.,  &
& lpcci, rpcc, lpking,  &
& mshxyz, mulx1, mulx2, muly1, muly2, mulz1, mulz2,  &
& mshnx, mshny, mshnz, mshnod,  &
& mshx1, mshy1, mshz1, mshx, mshy, mshz )
    end if
    strparts(1:3,1:3) = strparts(1:3,1:3) + strpcc(1:3,1:3)
    strparts(1:3,1:3) = strparts(1:3,1:3) * (-1.d0) / rvol

    !---symmetry operations on stress
!    call symk_stress( nfile, myid, nodes, strparts )

    if( ltimecnt ) then
        do nf = 1, 2
        if( loutfile(nf) ) then
           write(nfile(nf),*) ' '
           write(nfile(nf),*) 'Stress by PCC in GPa'
           do j = 1, 3
              write(nfile(nf),721) ( strparts(j,i)*prau, i = 1, 3 )
           end do
        end if
        end do
    end if

    str(1:3,1:3) = str(1:3,1:3) + strparts(1:3,1:3)
end if

!--- Stress from local pseudopotentials
strloc(1:3,1:3) = strloc(1:3,1:3) * (-1.d0) / rvol

if( ltimecnt ) then
    do nf = 1, 2
    if( loutfile(nf) ) then
       write(nfile(nf),*) ' '
       write(nfile(nf),*)  &
&                        'Stress by local pseudopotentials in GPa'
       do j = 1, 3
          write(nfile(nf),721) ( strloc(j,i)*prau, i = 1, 3 )
       end do
    end if
    end do
end if

str(1:3,1:3) = str(1:3,1:3) + strloc(1:3,1:3)


!--- Stress from nonlocal pseudopotentials
do i = 1, 3
   strnlc(i,i) = strnlc(i,i) + encl
do j = 1, 3
   strnlc(j,i) = strnlc(j,i) * (-1.d0) / rvol
end do
end do

if( ltimecnt ) then
    do nf = 1, 2
    if( loutfile(nf) ) then
       write(nfile(nf),*) ' '
       write(nfile(nf),*)  &
&                     'Stress by nonlocal pseudopotentials in GPa'
       do j = 1, 3
          write(nfile(nf),721) ( strnlc(j,i)*prau, i = 1, 3 )
       end do
    end if
    end do
end if

str(1:3,1:3) = str(1:3,1:3) + strnlc(1:3,1:3)


!--- Stress from direct coulomb
strclm(1:3,1:3) = strclm(1:3,1:3) * (-1.d0) / rvol

if( ltimecnt ) then
    do nf = 1, 2
    if( loutfile(nf) ) then
       write(nfile(nf),*) ' '
       write(nfile(nf),*) 'Stress by direct Coulomb in GPa'
       do j = 1, 3
          write(nfile(nf),721) ( strclm(j,i)*prau, i = 1, 3 )
       end do
    end if
    end do
end if

str(1:3,1:3) = str(1:3,1:3) + strclm(1:3,1:3)


!--- Stress from short-ranged HF energy
if( jhybrid /= 0 ) then
    do i = 1, 3
       strhfs(i,i) = strhfs(i,i) - exhfsr
    end do
    strhfs(1:3,1:3) = strhfs(1:3,1:3) * (-1.d0) / rvol

    if( ltimecnt ) then
        do nf = 1, 2
        if( loutfile(nf) ) then
           write(nfile(nf),*) ' '
           write(nfile(nf),*) 'Stress by short-ranged HF energy in GPa'
           do j = 1, 3
              write(nfile(nf),721) ( strhfs(j,i)*prau, i = 1, 3 )
           end do
        end if
        end do
    end if

    str(1:3,1:3) = str(1:3,1:3) + strhfs(1:3,1:3)

end if


!--- Stress from uniform electric field
if( lefield ) then
    do i = 1, 3
!       strefd(i,i) = strefd(i,i) - exhfsr
    end do
    strefd(1:3,1:3) = strefd(1:3,1:3) * (-1.d0) / rvol

    if( ltimecnt ) then
        do nf = 1, 2
        if( loutfile(nf) ) then
           write(nfile(nf),*) ' '
           write(nfile(nf),*) 'Stress by uniform electric field in GPa'
           do j = 1, 3
              write(nfile(nf),721) ( strefd(j,i)*prau, i = 1, 3 )
           end do
        end if
        end do
    end if

    str(1:3,1:3) = str(1:3,1:3) + strefd(1:3,1:3)

end if


!-----Stress from an empirical correction for vdW : DFT-D
if( ldftd ) call stress_dftd( nfile, myid, nodes, &
str, rvol, prau, ltimecnt )


!--- in [GPa]
!      do i = 1, 3
!      do j = 1, 3
!         str(j,i) = str(j,i)*prau
!      end do
!      end do
!--- Total stress
if( ltimecnt ) then
    do nf = 1, 2
    if( loutfile(nf) ) then
        write(nfile(nf),*) ' '
        write(nfile(nf),*) 'Total internal stress tensor in GPa'
        do j = 1, 3
           write(nfile(nf),721) ( str(j,i)*prau, i = 1, 3 )
        end do
    end if
    end do
end if

if( loutfile(1) ) then
    write(nfile(14),1155) nstep,  &
& str(1,1)*prau, str(2,2)*prau, str(3,3)*prau,  &
& str(2,3)*prau, str(3,1)*prau, str(1,2)*prau
1155     format(i6,10es16.8)
end if


return
end




subroutine strkine( nfile, myid, nodes,  &
& strparts, cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod,  &
& rvol, iod, occ, gx, gy, gz,  &
& lspin, nspnmx, gdcrsv, nspnod, lcgjsv, wegud, ncscale )
!-----------------------------------------------------------------------
!    Stress from kinetic E.
!-----------------------------------------------------------------------
use ncmagne_variables
implicit none
integer :: nfile(*), myid, nodes
real*8  :: strparts(3,3)
integer :: nplwex, nplw, nband, nbnod1, nbnod2, nbnod
integer :: ncscale
real*8  :: cgjr(nplwex,ncscale,*)
real*8  :: rvol
integer :: iod(*)
real*8  :: occ(*)
real*8  :: gx(0:*), gy(0:*), gz(0:*)
logical :: lspin, lcgjsv
integer :: nspnmx, nspnod
real*8  :: gdcrsv(*), wegud(*)

!-----declare local variables
integer :: nspin, i, j, ib, ig, mx
real*8  :: st11, st22, st33, st23, st31, st12, cgtmp2, occ2
real*8  :: bufst(6), bufstr(6)
logical :: lgamma


call get_lgamma( lgamma )

bufst(1:6) = 0.d0
!if( lgamma ) then

    if( ncscale == 1 ) then
        call strkine2( nfile, myid, nodes,  &
& bufst, cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod,  &
& iod, occ, gx, gy, gz,  &
& lspin, nspnmx, gdcrsv, nspnod, lcgjsv, wegud, ncscale )
    else
        !-----noncollinear magnetism
        call strkine2( nfile, myid, nodes,  &
& bufst, cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod,  &
& iod, occ, dkgx, dkgy, dkgz,  &
& lspin, nspnmx, gdcrsv, nspnod, lcgjsv, wegud, ncscale )
    end if

!else
!
!    if( lspin ) then
!        call strkine_k( nfile, myid, nodes,  &
!& bufst, nspnmx, wegud, nband, nbnod, iod, ncscale )
!    else
!        call strkine_k( nfile, myid, nodes,  &
!& bufst, nspnmx, occ, nband, nbnod, iod, ncscale )
!    end if
!
!end if


call gdsum(bufst,6,bufstr)
call unify_sumn(bufst,6,bufstr)
strparts(1,1) = bufst(1)
strparts(2,2) = bufst(2)
strparts(3,3) = bufst(3)
strparts(2,3) = bufst(4)
strparts(3,1) = bufst(5)
strparts(1,2) = bufst(6)
strparts(3,2) = strparts(2,3)
strparts(1,3) = strparts(3,1)
strparts(2,1) = strparts(1,2)

strparts(1:3,1:3) = strparts(1:3,1:3) * 2.d0 / rvol


return
end




subroutine strkine2( nfile, myid, nodes,  &
& bufst, cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod,  &
& iod, occ, dkgx, dkgy, dkgz,  &
& lspin, nspnmx, gdcrsv, nspnod, lcgjsv, wegud, ncscale )
!-----------------------------------------------------------------------
!    Stress from kinetic E.
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
real*8  :: bufst(6)
integer :: nplwex, nplw, nband, nbnod1, nbnod2, nbnod
integer :: ncscale
real*8  :: cgjr(nplwex,ncscale,*)
integer :: iod(*)
real*8  :: occ(*)
real*8  :: dkgx(0:*), dkgy(0:*), dkgz(0:*)
logical :: lspin, lcgjsv
integer :: nspnmx, nspnod
real*8  :: gdcrsv(*), wegud(*)

!-----declare local variables
integer :: nspin, i, j, ib, ig, mx
real*8  :: st11, st22, st33, st23, st31, st12, cgtmp2, occ2, spdeg


if( ncscale == 1 ) then
    spdeg = 2.d0
  else
    spdeg = 1.d0
end if
do nspin = 1, nspnmx

   if( lspin ) then
       call stspud( nspin, cgjr, gdcrsv, nspnod, lcgjsv )
       call ldocc( nspin, occ, wegud, nband )
   end if

do i = 1, nbnod
   ib = iod(i)

   st11 = 0.d0
   st22 = 0.d0
   st33 = 0.d0
   st23 = 0.d0
   st31 = 0.d0
   st12 = 0.d0
   do mx = 1, ncscale
   do ig = 1, nplw
      cgtmp2 = cgjr(2*ig+1,mx,i)*cgjr(2*ig+1,mx,i)  &
&            + cgjr(2*ig+2,mx,i)*cgjr(2*ig+2,mx,i)
      st11 = st11 + dkgx(ig)*dkgx(ig)*cgtmp2
      st22 = st22 + dkgy(ig)*dkgy(ig)*cgtmp2
      st33 = st33 + dkgz(ig)*dkgz(ig)*cgtmp2
      st23 = st23 + dkgy(ig)*dkgz(ig)*cgtmp2
      st31 = st31 + dkgz(ig)*dkgx(ig)*cgtmp2
      st12 = st12 + dkgx(ig)*dkgy(ig)*cgtmp2
   end do
   end do
   occ2 = occ(ib)*spdeg
   bufst(1) = bufst(1) + st11*occ2
   bufst(2) = bufst(2) + st22*occ2
   bufst(3) = bufst(3) + st33*occ2
   bufst(4) = bufst(4) + st23*occ2
   bufst(5) = bufst(5) + st31*occ2
   bufst(6) = bufst(6) + st12*occ2

end do

end do


return
end




subroutine strhart( nfile, myid, nodes,  &
& strparts, rhgr, nplw5ex, nplw5, recnrm, gx, gy, gz )
!-----------------------------------------------------------------------
!    Stress from Hartree E.
!-----------------------------------------------------------------------
implicit none

integer :: myid, nodes
integer, dimension(*) :: nfile
real*8, dimension(3,3) :: strparts
integer :: nplw5ex
integer :: nplw5
real*8, dimension(nplw5ex) :: rhgr
real*8, dimension(0:nplw5) :: recnrm
real*8, dimension(0:nplw5) :: gx, gy, gz

real*8, parameter :: pi8 = 3.14159265358979323846d0 * 8.d0
real*8  :: p4prec, vr, vi, buff
real*8  :: sxx, syy, szz, syz, szx, sxy
integer :: ig


sxx = 0.0d0
syy = 0.0d0
szz = 0.0d0
syz = 0.0d0
szx = 0.0d0
sxy = 0.0d0
do ig = 1, nplw5
   p4prec = pi8/recnrm(ig)
   vr  = p4prec*rhgr(2*ig+1)
   vi  = p4prec*rhgr(2*ig+2)
   buff = vr*vr + vi*vi
   sxx = sxx + gx(ig)*gx(ig)*buff
   syy = syy + gy(ig)*gy(ig)*buff
   szz = szz + gz(ig)*gz(ig)*buff
   syz = syz + gy(ig)*gz(ig)*buff
   szx = szx + gz(ig)*gx(ig)*buff
   sxy = sxy + gx(ig)*gy(ig)*buff
end do
strparts(1,1) = -sxx*2.d0/pi8
strparts(2,2) = -syy*2.d0/pi8
strparts(3,3) = -szz*2.d0/pi8
strparts(2,3) = -syz*2.d0/pi8
strparts(3,1) = -szx*2.d0/pi8
strparts(1,2) = -sxy*2.d0/pi8
strparts(3,2) = strparts(2,3)
strparts(1,3) = strparts(3,1)
strparts(2,1) = strparts(1,2)


return
end
