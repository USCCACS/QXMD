



subroutine pwpxc( nfile, myid, nodes,  &
& eexc, evexc, ecnl, ecor, evcor, lcstress, lvdw_pre, ltimecnt, t_comm )
!-----------------------------------------------------------------------
!    exchange & correlation   potential : vexc
!-----------------------------------------------------------------------
use param
implicit none
integer :: nfile(*), myid, nodes
real*8  :: eexc, evexc, ecnl, ecor, evcor
logical :: lcstress, lvdw_pre, ltimecnt
real*8  :: t_comm


!if( .not.lnoncollinear ) then
    call pwpxc2( nfile, myid, nodes,  &
& eexc, evexc, ecnl, ecor, evcor, lcstress, lvdw_pre, ltimecnt, t_comm )
!else
!    !-----noncollinear magnetism
!    call ncpwpxc2( nfile, myid, nodes,  &
!& eexc, evexc, ecnl, ecor, evcor, lcstress, lvdw_pre, ltimecnt, t_comm )
!end if


return
end




subroutine pwpxc2( nfile, myid, nodes,  &
& eexc, evexc, ecnl, ecor, evcor, lcstress, lvdw_pre, ltimecnt, t_comm )
!-----------------------------------------------------------------------
!    exchange & correlation   potential : vexc
!-----------------------------------------------------------------------
!    exchange energy       : eexc
!    exchange potential    : evexc
!    correlation energy    : ecor
!    correlation potential : evcor
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
real*8  :: eexc, evexc, ecnl, ecor, evcor
logical :: lcstress, lvdw_pre, ltimecnt
real*8  :: t_comm


ecnl = 0.d0
if( jgga == 1 ) then

    !--- LDA ( perdew-zunger + ceperley-alder )
    if( lspin ) then

        call xclsda( rho, rhoud, rhocore, mshnod(1), rdelv,  &
& vlocud, eeexc, noddatx, eexc, evexc, ecor, evcor )

      else

        call xclda( rho, rhocore, mshnod(1), rdelv,  &
&               vexc, eeexc, eexc, evexc, ecor, evcor )

    end if

 else

    !--- GGA
    spinif: if( lspin ) then

!        if( jgga == 2 ) then

!            hybridif: if( jhybrid == 0 ) then

                !--- GGA ( PBE )
                call xcpbe_spin( nfile, myid, nodes,  &
& vlocud, eeexc, eexc, evexc, ecor, evcor,  &
& rho, rhoud, rhocore, mshnod(1), rdelv,  &
& x, rk, xx, vk, tmpk, tmpl, tmpm,  &
& tmpn, hdiag, vexc, thrhgr, nplw5ex, nplw5,  &
& glocal, apk, eigr, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d,  &
& fft3x, fft3y, fftwork, ijkgd, nplw5, gx, gy, gz, t_comm,  &
& lcstress, strgga )

!            else hybridif
!
!                !--- HSE hybrid functional with PBE
!                call xcpbe_HSE_spin( nfile, myid, nodes,  &
!& vlocud, eeexc, eexc, evexc, ecor, evcor,  &
!& rho, rhoud, rhocore, mshnod(1), rdelv,  &
!& x, rk, xx, vk, tmpk, tmpl, tmpm,  &
!& tmpn, hdiag, vexc, thrhgr, nplw5ex, nplw5,  &
!& glocal, apk, eigr, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d,  &
!& fft3x, fft3y, fftwork, ijkgd, nplw5, gx, gy, gz, t_comm,  &
!& lcstress, strgga, hmixing, hrange )
!
!            end if hybridif

!          else
            !--- GGA ( RPBE )
!            call xcrpbe_spin( nfile, myid, nodes,  &
!& vlocud, eeexc, eexc, evexc, ecor, evcor,  &
!& rho, rhoud, rhocore, mshnod(1), rdelv,  &
!& x, rk, xx, vk, tmpk, tmpl, tmpm,  &
!& tmpn, hdiag, vexc, thrhgr, nplw5ex, nplw5,  &
!& glocal, apk, eigr, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d,  &
!& fft3x, fft3y, fftwork, ijkgd, nplw5, gx, gy, gz, t_comm,  &
!& lcstress, strgga )
!        end if

    else spinif

!        if( jgga == 2 ) then

!            hybridif2: if( jhybrid == 0 ) then

                if( .not.lvdw .or. lvdw_pre ) then
                    !--- GGA ( PBE )
                    call xcpbe( nfile, myid, nodes,  &
& vexc, eeexc, eexc, evexc, ecor, evcor, rho, rhocore,  &
& mshnod(1), rdelv, x, rk, xx, vk, thrhgr, nplw5ex, nplw5,  &
& glocal, apk, eigr, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d,  &
& fft3x, fft3y, fftwork, ijkgd, nplw5, gx, gy, gz, t_comm,  &
& lcstress, strgga )
!                  else
!                    !--- van der Waals DFT
!                    call xcvdw( nfile, myid, nodes,  &
!& vexc, eeexc, eexc, evexc, ecor, ecnl, evcor, rho, rhocore,  &
!& rdelv, x, rk, xx, vk, thrhgr, nplw5ex, nplw5,  &
!& glocal, apk, eigr, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d,  &
!& fft3x, fft3y, fftwork, ijkgd, nplw5, gx, gy, gz, t_comm,  &
!& lcstress, strgga, lorthrhmbc, &
!& hcell, nd1v, mshglb, nd1vks(1), nd1vks(2), nd1vks(3), lvacuum,  &
!& mshxyz(mulpit(1)),  &
!& mulx1(1), mulx2(1), muly1(1), muly2(1), mulz1(1), mulz2(1),  &
!& mshnx(mulpms(1)), mshny(mulpms(1)), mshnz(mulpms(1)), mshnod(1),  &
!& mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1), ltimecnt )
                end if

!            else hybridif2
!
!                !--- HSE hybrid functional with PBE
!                call xcpbe_HSE( nfile, myid, nodes,  &
!& vexc, eeexc, eexc, evexc, ecor, evcor, rho, rhocore,  &
!& mshnod(1), rdelv, x, rk, xx, vk, thrhgr, nplw5ex, nplw5,  &
!& glocal, apk, eigr, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d,  &
!& fft3x, fft3y, fftwork, ijkgd, nplw5, gx, gy, gz, t_comm,  &
!& lcstress, strgga, hmixing, hrange )
!
!            end if hybridif2

!          else

            !--- GGA ( RPBE )
!            call xcrpbe( nfile, myid, nodes,  &
!& vexc, eeexc, eexc, evexc, ecor, evcor, rho, rhocore,  &
!& mshnod(1), rdelv, x, rk, xx, vk, thrhgr, nplw5ex, nplw5,  &
!& glocal, apk, eigr, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d,  &
!& fft3x, fft3y, fftwork, ijkgd, nplw5, gx, gy, gz, t_comm,  &
!& lcstress, strgga )
!        end if

    end if spinif

end if


return
end




subroutine xclsda( rho, rhoud, rhocore, mshnod, rdelv,  &
&               vexc, eeexc, noddatx, eexc, evexc, ecor, evcor )
!-----------------------------------------------------------------------
! ( LSDA )
!
!  ( input )
!       rho   : total charge density
!     rhoud   : difference between up and down densities
!
!  ( output )
!    exchange & correlation   potential : vexc
!    exchange energy                    : eexc
!    exchange potential                 : evexc
!    correlation energy                 : ecor
!    correlation potential              : evcor
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
parameter( pi4 = 12.56637061436d0 )
dimension rho(mshnod), rhoud(mshnod)
dimension rhocore(mshnod)
dimension vexc(mshnod,2)
dimension eeexc(mshnod)
dimension vx(2), vc(2)
dimension dbuf(4), dbufr(4)


! --- check to avoid numerical instability ---
!do ir = 1, mshnod
!   if( rho(ir).lt.0.d0 ) then
!       rhoud(ir) = 0.d0
!   else if( 0.5d0*( rho(ir) + rhoud(ir) ).lt.0.d0 ) then
!       rhoud(ir) = - rho(ir)
!   else if( 0.5d0*( rho(ir) - rhoud(ir) ).lt.0.d0 ) then
!       rhoud(ir) = rho(ir)
!   end if
!end do
do m1 = 1, mshnod
   if( rhocore(m1).ge.0.d0 ) then
       rx = rho(m1) + rhocore(m1)
     else
       rx = rho(m1)
   end if
   if( rx.ge.1.d-30 ) then
       zeta = rhoud(m1)/rx
       if( abs(zeta) > 1.d0 ) then
           zeta = sign(1.d0,rhoud(m1))
       end if
       vexc(m1,1) = zeta
   end if
end do

eexc  = 0.d0
evexc = 0.d0
ecor  = 0.d0
evcor = 0.d0
do m1 = 1, mshnod
   if( rhocore(m1).ge.0.d0 ) then
       rhot = rho(m1) + rhocore(m1)
     else
       rhot = rho(m1)
   end if
   if( rhot.lt.1.d-30 ) then
       vexc(m1,1) = 0.d0
       vexc(m1,2) = 0.d0
       eeexc(m1)  = 0.d0
     else
       zeta = vexc(m1,1)
       rs = ( 3.D0/( pi4*rhot ) )**( 1.D0/3.D0 )
       call vxclsda( rs, zeta, ex, vx, ec, vc, 0 )
       vexc(m1,1) = 2.d0*( vx(1) + vc(1) )
       vexc(m1,2) = 2.d0*( vx(2) + vc(2) )
       eeexc(m1)  = 2.d0*( ex + ec )
       rhoup = 0.5d0*( rho(m1) + rhoud(m1) )
       rhodn = 0.5d0*( rho(m1) - rhoud(m1) )
       eexc  = eexc  + ex * rhot
       evexc = evexc + vx(1) * rhoup + vx(2) * rhodn
       ecor  = ecor  + ec * rhot
       evcor = evcor + vc(1) * rhoup + vc(2) * rhodn
   end if
end do
dbuf(1) = eexc
dbuf(2) = evexc
dbuf(3) = ecor
dbuf(4) = evcor
call gdsum(dbuf,4,dbufr)
eexc  = dbuf(1) * 2.d0 * rdelv
evexc = dbuf(2) * 2.d0 * rdelv
ecor  = dbuf(3) * 2.d0 * rdelv
evcor = dbuf(4) * 2.d0 * rdelv


return
end




module pbe_variables
!-----------------------------------------------------------------------
! type declaration of shared variables in vxc.f
!-----------------------------------------------------------------------
implicit none

real*8,  parameter :: c1b3 = 1.d0/3.d0
real*8,  parameter :: c4b3 = 4.d0/3.d0
real*8,  parameter :: c5b3 = 5.d0/3.d0
real*8,  parameter :: c7b3 = 7.d0/3.d0
real*8,  parameter :: c1b6m = -1.d0/6.d0
real*8,  parameter :: c8b3   =  8.d0/3.d0
real*8,  parameter :: c16b3  = 16.d0/3.d0
real*8,  parameter :: c2p4b3 = 2.5198420997897463295344212145565d0 ! = 2^(4/3)
real*8,  parameter :: c6sqpi = 0.1184723887405400601343407529603d0 ! = 2^(1/3)/(6*sqrt(pi))
real*8,  parameter :: sqpi   = 1.7724538509055160272981674833411d0 ! = sqrt(pi)

!-----variables for PBE exchange functional
real*8  :: uk = 0.8040d0

real*8  :: rs0 = 0.62035049090d0
!--- parameter for exchange energy
real*8  :: ax = -0.738558766382022405884230032680836d0
real*8  :: bx =  0.16162045967d0
real*8  :: cx =  0.16162045967d0*2.d0
real*8  :: um =  0.2195149727645171d0
!      data uk / 0.8040d0 /

!--- parameter for correlation energy
real*8  :: a0    = 0.0621814d0
real*8  :: alp01 = 0.21370d0
real*8  :: bet01 = 7.5957d0
real*8  :: bet02 = 3.5876d0
real*8  :: bet03 = 1.6382d0
real*8  :: bet04 = 0.49294d0
real*8  :: p0    = 1.d0
real*8  :: betapw = 0.06672455060314922d0
real*8  :: t0 = 0.2519289703424d0
real*8  :: h1 = 0.03109069086965489503494086371273d0

!--- parameter for correlation energy with spin polarization
real*8  :: a1    = 0.0310907d0
real*8  :: alp11 = 0.20548d0
real*8  :: bet11 = 14.1189d0
real*8  :: bet12 = 6.1977d0
real*8  :: bet13 = 3.3662d0
real*8  :: bet14 = 0.62517d0
real*8  :: p1    = 1.d0
real*8  :: aa    = 0.0337738d0
real*8  :: alpa1 = 0.11125d0
real*8  :: beta1 = 10.3570d0
real*8  :: beta2 = 3.6231d0
real*8  :: beta3 = 0.88026d0
real*8  :: beta4 = 0.49671d0
real*8  :: pa    = 1.d0
real*8  :: f1   = 1.92366105094d0
real*8  :: fd20 = 0.5848223622608d0
real*8  :: eta = 1.d-12


!--- parameter for relativistic correction
real*8  :: BB = 0.0140D0

save

end module




subroutine xcpbe_spin( nfile, myid, nodes,  &
& vlocud, eeexc, eexc, evexc, ecor, evcor,  &
& rho, rhoud, prho, mshnod, rdelv,  &
& x, dxrho, dyrho, dzrho, dxrhod, dyrhod, dzrhod,  &
& drhot, drhot1, drhot2, thrhgr, nplw5ex, nplw5,  &
& glx, gly, glz, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, gx, gy, gz, t_comm,  &
& lcstress, strgga )
!-----------------------------------------------------------------------
! Generalized gradient corrected exchange correlation potential by PBE
!     for spin polarized case
!              ^^^^^^^^^
!  input
!       rho       : total charge density
!     rhoud       : difference between up and down densities
!  output
!    exchange & correlation   potential : vlocud
!    exchange energy                    : eexc
!    exchange potential                 : evexc
!    correlation energy                 : ecor
!    correlation potential              : evcor
!-----------------------------------------------------------------------
use pbe_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: mshnod
real*8  :: vlocud(mshnod,2), eeexc(mshnod), eexc, evexc, ecor, evcor
real*8  :: rho(mshnod), rhoud(mshnod), prho(mshnod)
real*8  :: rdelv
real*8  :: x(*), dxrho(*), dyrho(*), dzrho(*), dxrhod(*), dyrhod(*), dzrhod(*)
real*8  :: drhot(*), drhot1(*), drhot2(*)
integer :: nplw5ex, nplw5
real*8  :: thrhgr(*)
real*8  :: glx(*), gly(*), glz(*)
integer :: mftnod(*), mftdsp(*), mfd2ft(*), ntotfd, nd1vks(*), kfft0
real*8  :: fft3x(*), fft3y(*)
complex*16 :: fftwork(*)
integer :: ngenh
integer, dimension(-ngenh:ngenh) :: ijkgd
real*8, dimension(0:ngenh) :: gx, gy, gz
real*8  :: t_comm
logical :: lcstress
real*8, dimension(3,3) :: strgga

!-----declare local variables
!integer, parameter :: c1b3 = 1.d0/3.d0
!integer, parameter :: c4b3 = 4.d0/3.d0
!integer, parameter :: c7b3 = 7.d0/3.d0
!integer, parameter :: c1b6m = -1.d0/6.d0
real*8, dimension(6) :: dbuf, dbufr
integer :: ir, m1, ijk
real*8  :: rx, zeta, rhohf, abdrho, dxud, dyud, dzud, rhot1, rhot2
real*8  :: rx1, ex, vxup, vxdn, rho3, rs, adrho, adrbro, ex0, s, s2, ul, fsbb, fst, fss
real*8  :: dfsds, dfxdn, dfxddn
real*8  :: zeta2, zeta4, zetp1, zetm1, zetp2, zetm2, fzeta, rs1h
real*8  :: g0bb, g0b2, g0ln, g0rs1, g0rs, g1bb, g1b2, g1ln, g1rs1, g1rs, gabb, gab2, galn, gars1, gars
real*8  :: ecrsz, ec, dg0dr1, dg0drs, dg1dr1, dg1drs, dgadr1, dgadrs, drsdn, decdn, dfdz
real*8  :: zeta3, zet34f, zet4df, decdz, gzeta, gzetar, gzeta2, gzeta3, t, t2, t4
real*8  :: dfcdn, dfcdz, vcup, vcdn, h2, g3beta, eecrsz, arsz, arsz2
real*8  :: h0bs, h0bb, h0b2, h0ln, h0trsz, ech, h0bsa, h0bba, h0bst, h0bbt, dh0c1, dh0c, dh0da, dh0dt
real*8  :: dadcc, dadn, dtdn, dh0dn, dtddn, dh0ddn, dfcddn, dadzc, dadz, dh0dg, dgdz, dh0dz
real*8  :: dxrhi, dyrhi, dzrhi, dxrhdi, dyrhdi, dzrhdi, drhot1i, drhot2i, rhoup, rhodn
!save rs0, ax, bx, um, &
!&    a0, alp01, bet01, bet02, bet03, bet04, p0,  &
!&    a1, alp11, bet11, bet12, bet13, bet14, p1,  &
!&    aa, alpa1, beta1, beta2, beta3, beta4, pa,  &
!&    betapw, t0, h1,  &
!&    f1, fd20, eta
!
!data rs0 / 0.62035049090d0 /
!!--- parameter for exchange energy
!data ax / -0.738558766382022405884230032680836d0 /
!data bx /  0.16162045967d0 /
!data um / 0.2195149727645171d0 /
!!      data uk / 0.8040d0 /
!
!!--- parameter for correlation energy
!data a0, alp01, bet01, bet02, bet03, bet04, p0  &
!&  / 0.0621814d0, 0.21370d0,  &
!&    7.5957d0, 3.5876d0, 1.6382d0, 0.49294d0, 1.d0 /
!data a1, alp11, bet11, bet12, bet13, bet14, p1  &
!&  / 0.0310907d0, 0.20548d0,  &
!&   14.1189d0, 6.1977d0, 3.3662d0, 0.62517d0, 1.d0 /
!data aa, alpa1, beta1, beta2, beta3, beta4, pa  &
!&  / 0.0337738d0, 0.11125d0,  &
!&   10.3570d0, 3.6231d0, 0.88026d0, 0.49671d0, 1.d0 /
!data betapw / 0.06672455060314922d0 /
!data t0 / 0.2519289703424d0 /
!data h1 / 0.03109069086965489503494086371273d0 /
!data f1, fd20 / 1.92366105094d0,  0.5848223622608d0 /
!data eta / 1.d-12 /


! --- check to avoid numerical instability ---
!do ir = 1, mshnod
!   if( rho(ir).lt.0.d0 ) then
!       rhoud(ir) = 0.d0
!   else if( 0.5d0*( rho(ir) + rhoud(ir) ).lt.0.d0 ) then
!       rhoud(ir) = - rho(ir)
!   else if( 0.5d0*( rho(ir) - rhoud(ir) ).lt.0.d0 ) then
!       rhoud(ir) = rho(ir)
!   end if
!end do
do m1 = 1, mshnod
   if( prho(m1).ge.0.d0 ) then
       rx = rho(m1) + prho(m1)
     else
       rx = rho(m1)
   end if
   if( rx.ge.1.d-30 ) then
       zeta = rhoud(m1)/rx
       if( abs(zeta) > 1.d0 ) then
           zeta = sign(1.d0,rhoud(m1))
       end if
     else
       zeta = 0.d0
   end if
!cc         rhohf = 0.5d0*( rho(m1) + prho(m1) )
   rhohf = 0.5d0*rx
   vlocud(m1,1) = rhohf*( 1.d0 + zeta )
   vlocud(m1,2) = rhohf*( 1.d0 - zeta )
end do


!  --- gradient of up-spin charge density ---
do ir = 1, mshnod
   x(ir) = vlocud(ir,1)
end do

!--- unify charge density
!      call unifylc( nfile, myid, nodes,
!     & x, mshnod, glx, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks )

call graden( nfile, myid, nodes,  &
& dxrho, dyrho, dzrho, mshnod, x, glx, gly, glz,  &
& thrhgr, nplw5ex, nplw5,  &
& mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, gx, gy, gz )


!  --- gradient of down-spin charge density ---
do ir = 1, mshnod
   x(ir) = vlocud(ir,2)
end do

!--- unify charge density
!      call unifylc( nfile, myid, nodes,
!     & x, mshnod, glx, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks )

call graden( nfile, myid, nodes,  &
& dxrhod, dyrhod, dzrhod, mshnod, x, glx, gly, glz,  &
& thrhgr, nplw5ex, nplw5,  &
& mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, gx, gy, gz )


eexc  = 0.d0
evexc = 0.d0
ecor  = 0.d0
evcor = 0.d0
!=======================================================================
do ijk = 1, mshnod
   abdrho = dxrho(ijk)*dxrho(ijk) + dyrho(ijk)*dyrho(ijk)  &
&         + dzrho(ijk)*dzrho(ijk)
   drhot1(ijk) = sqrt(abdrho)

   abdrho = dxrhod(ijk)*dxrhod(ijk) + dyrhod(ijk)*dyrhod(ijk)  &
&         + dzrhod(ijk)*dzrhod(ijk)
   drhot2(ijk) = sqrt(abdrho)

   dxud = dxrho(ijk) + dxrhod(ijk)
   dyud = dyrho(ijk) + dyrhod(ijk)
   dzud = dzrho(ijk) + dzrhod(ijk)
   abdrho = dxud*dxud + dyud*dyud + dzud*dzud
   drhot(ijk) = sqrt(abdrho)
end do

do ijk = 1, mshnod

   rx = vlocud(ijk,1) + vlocud(ijk,2)
   if( rx.lt.1.d-30 ) then
       vlocud(ijk,1) = 0.d0
       vlocud(ijk,2) = 0.d0
       drhot1(ijk) = 0.d0
       drhot2(ijk) = 0.d0
       drhot(ijk)  = 0.d0
       eeexc(ijk)  = 0.d0
   else

   rhot1 = vlocud(ijk,1)
   rhot2 = vlocud(ijk,2)
!--- exchange terms ----------------------------------------------------

!  --- for up-spin -------------
rx1   = 2.d0*rhot1
!  -----------------------------
if( rx1 .lt. 1.d-30 ) then
    ex   = 0.d0
    vxup = 0.d0
    drhot1(ijk) = 0.d0
else


rho3 = rx1**c1b3
rs   = rs0/rho3

adrho  = 2.d0*drhot1(ijk)
adrbro = adrho / rx1

ex0  = ax * rho3

s  = bx * adrbro / rho3
if( s.lt.1.d-15 ) then
    ex   = ex0 * rx1
    vxup = c4b3*ex0
    drhot1(ijk) = 0.d0
  else
    s2 = s*s

    ul = um/uk
    fsbb = 1.d0 + s2*ul
    fsbb = 1.d0/fsbb
    fst  = - uk*fsbb
    fss  = 1.d0 + uk + fst

!     --- exchange energy ---
    ex = ex0*fss * rx1


    dfsds  = - 2.d0*ul*s * fst*fsbb

    dfxdn  = c4b3*ex0 * ( fss - s*dfsds )
    dfxddn = ax * bx * dfsds

    vxup = dfxdn
    drhot1(ijk) = dfxddn/drhot1(ijk)

end if
end if




!  --- for down-spin -----------
rx1   = 2.d0*rhot2
!  -----------------------------
if( rx1 .lt. 1.d-30 ) then
    vxdn = 0.d0
    drhot2(ijk) = 0.d0
else


rho3 = rx1**c1b3
rs   = rs0/rho3

adrho  = 2.d0*drhot2(ijk)
adrbro = adrho / rx1

ex0  = ax * rho3

s  = bx * adrbro / rho3
if( s.lt.1.d-15 ) then
    ex   = ex + ex0 * rx1
    vxdn = c4b3*ex0
    drhot2(ijk) = 0.d0
  else
    s2 = s*s

    ul = um/uk
    fsbb = 1.d0 + s2*ul
    fsbb = 1.d0/fsbb
    fst  = - uk*fsbb
    fss  = 1.d0 + uk + fst

!     --- exchange energy ---
    ex = ex + ex0*fss * rx1


    dfsds  = - 2.d0*ul*s * fst*fsbb

    dfxdn  = c4b3*ex0 * ( fss - s*dfsds )
    dfxddn = ax * bx * dfsds

    vxdn = dfxdn
    drhot2(ijk) = dfxddn/drhot2(ijk)

end if
end if

!     --- exchange energy ---
ex = 0.5d0 * ex




!--- correlation terms -------------------------------------------------
!      rx    = rhot1   + rhot2

rho3 = rx**c1b3
rs   = rs0/rho3

zeta = ( rhot1 - rhot2 )/rx
zeta2 = zeta*zeta
zeta4 = zeta2*zeta2
if( zeta.gt.0.d0 ) then
    zetp1 = 1.d0 + zeta
    zetm1 = 1.d0 - zeta
    if( zetm1.lt.1.d-08 ) zetm1 = 2.d0 * rhot2/rhot1
  else
    zetm1 = 1.d0 - zeta
    zetp1 = 1.d0 + zeta
    if( zetp1.lt.1.d-08 ) zetp1 = 2.d0 * rhot1/rhot2
end if
zetp1 = zetp1**c1b3
zetm1 = zetm1**c1b3
zetp2 = zetp1*zetp1
zetm2 = zetm1*zetm1
fzeta = f1*( zetp2*zetp2 + zetm2*zetm2 - 2.d0 )

rs1h = dsqrt(rs)


g0bb = rs1h*( bet01 + rs1h*( bet02 + rs1h*( bet03 + rs1h*bet04 )))
g0b2 = 1.d0/(a0*g0bb)
if( abs(g0b2).gt.1.d-10 ) then
    g0ln = log( 1.d0 + g0b2 )
  else
    g0ln = g0b2
end if
g0rs1 = 1.d0 + alp01*rs
g0rs  = -a0*g0rs1*g0ln

g1bb = rs1h*( bet11 + rs1h*( bet12 + rs1h*( bet13 + rs1h*bet14 )))
g1b2 = 1.d0/(a1*g1bb)
if( abs(g1b2).gt.1.d-10 ) then
    g1ln = log( 1.d0 + g1b2 )
  else
    g1ln = g1b2
end if
g1rs1 = 1.d0 + alp11*rs
g1rs  = -a1*g1rs1*g1ln

gabb = rs1h*( beta1 + rs1h*( beta2 + rs1h*( beta3 + rs1h*beta4 )))
gab2 = 1.d0/(aa*gabb)
if( abs(gab2).gt.1.d-10 ) then
    galn = log( 1.d0 + gab2 )
  else
    galn = gab2
end if
gars1 = 1.d0 + alpa1*rs
gars  = -aa*gars1*galn

ecrsz = g0rs + fzeta*( -gars*fd20*( 1.d0 - zeta4 )  &
&                      + ( g1rs - g0rs )*zeta4     )
!  --- correlation energy ---
ec = ecrsz




   dg0dr1 = 0.5d0*bet01/rs1h + bet02  &
&         + rs1h*(1.5d0*bet03 + 2.d0*bet04*rs1h)
   dg0drs = -a0*alp01*g0ln  &
&         + g0rs1*dg0dr1 / ( ( 1.d0 + g0b2 )*g0bb*g0bb )

   dg1dr1 = 0.5d0*bet11/rs1h + bet12  &
&         + rs1h*(1.5d0*bet13 + 2.d0*bet14*rs1h)
   dg1drs = -a1*alp11*g1ln  &
&         + g1rs1*dg1dr1 / ( ( 1.d0 + g1b2 )*g1bb*g1bb )

   dgadr1 = 0.5d0*beta1/rs1h + beta2  &
&         + rs1h*(1.5d0*beta3 + 2.d0*beta4*rs1h)
   dgadrs = -aa*alpa1*galn  &
&         + gars1*dgadr1 / ( ( 1.d0 + gab2 )*gabb*gabb )

   drsdn = - rs * c1b3

   decdn = drsdn*( dg0drs + fzeta*( -dgadrs*fd20*( 1.d0 - zeta4 )  &
&                         + ( dg1drs - dg0drs )*zeta4     ) )



   dfdz = c4b3 * f1*( zetp1 - zetm1 )

   zeta3 = zeta2*zeta
   zet34f = 4.d0*zeta3*fzeta
   zet4df = zeta4*dfdz
   decdz  = gars*fd20*( zet34f - dfdz + zet4df )  &
&                      + ( g1rs - g0rs )*( zet34f + zet4df )


!==== gradient terms ================
gzeta = 0.5d0*( zetp2 + zetm2 )
gzetar= 1.d0/gzeta
gzeta2= gzeta*gzeta
gzeta3= gzeta2*gzeta

adrho = drhot(ijk)
adrbro = adrho / rx

t  = t0*adrbro
t2 = t*t / (gzeta2*rho3)
if( t2.lt.1.d-30 ) then
    dfcdn = ec + decdn
    dfcdz = decdz
    vcup = dfcdn + dfcdz*(  1.d0 - zeta )
    vcdn = dfcdn + dfcdz*( -1.d0 - zeta )
    drhot(ijk) = 0.d0
else

h2 = betapw/h1

g3beta = gzeta3*betapw
eecrsz = exp(-ecrsz/(gzeta3*h1))
arsz   = h2 / ( eecrsz - 1.d0 )
arsz2  = arsz*arsz

t4 = t2*t2

h0bs =             t2 + arsz *t4 
h0bb = 1.d0 + arsz*t2 + arsz2*t4 
h0b2 = h0bs/h0bb
if( abs(h0b2).gt.1.d-10 ) then
    h0ln = log( 1.d0 + h2*h0b2 )
  else
    h0ln = h2*h0b2
end if
h0trsz = gzeta3*h1*h0ln

ech  = h0trsz

!  --- correlation energy by gradient terms ---
ec = ec + ech




   h0bsa =                t4
   h0bba = t2 + 2.d0*arsz*t4
   h0bst = 1.d0 + 2.d0*arsz *t2
   h0bbt = arsz + 2.d0*arsz2*t2

   dh0c1 = h2*h0bs + h0bb
   dh0c  = g3beta/dh0c1
   dh0da = dh0c*( h0bsa - h0b2*h0bba )
   dh0dt = dh0c*( h0bst - h0b2*h0bbt )

   dadcc= eecrsz*arsz2/(g3beta)
   dadn = dadcc * decdn
   dtdn = - c7b3 * t2

   dh0dn = dh0da*dadn + dh0dt*dtdn


   dfcdn = ec + decdn + dh0dn



   dtddn = 2.d0*t2/adrbro

   dh0ddn = dh0dt*dtddn

   dfcddn = dh0ddn
   drhot(ijk) = dfcddn/drhot(ijk)


if( abs(zeta).lt.1.d-30 ) then
    vcup = dfcdn
    vcdn = dfcdn
  else

   dadzc= dh0da*dadcc
   dadz = dadzc * decdz
   dh0dg = ( 3.d0*( h0trsz - dadzc * ecrsz )  &
&           - 2.d0*t2*dh0dt                  )* gzetar
   dgdz = c1b3 * ( ( ( 1.d0 + zeta )**2 + eta )**c1b6m  &
&                - ( ( 1.d0 - zeta )**2 + eta )**c1b6m )
   dh0dz = dh0dg * dgdz

   dfcdz = decdz + dh0dz + dadz

    vcup = dfcdn + dfcdz*(  1.d0 - zeta )
    vcdn = dfcdn + dfcdz*( -1.d0 - zeta )
end if

end if



   vlocud(ijk,1) = vxup + vcup
   vlocud(ijk,2) = vxdn + vcdn

   eexc  = eexc + ex
   ecor  = ecor + rx*ec
   eeexc(ijk) = ex/rx + ec

   end if
end do
!=======================================================================


!   --- gradient of d_fxc(up)/d_dn ---
if( lcstress ) then
    strgga(1:3,1:3) = 0.d0
    do ijk = 1, mshnod
       dxrhi  = dxrho(ijk)
       dyrhi  = dyrho(ijk)
       dzrhi  = dzrho(ijk)
       dxrhdi = dxrhod(ijk)
       dyrhdi = dyrhod(ijk)
       dzrhdi = dzrhod(ijk)
       drhot1i = drhot1(ijk) + drhot(ijk)
       dxrho(ijk) = drhot1i*dxrhi + drhot(ijk)*dxrhdi
       dyrho(ijk) = drhot1i*dyrhi + drhot(ijk)*dyrhdi
       dzrho(ijk) = drhot1i*dzrhi + drhot(ijk)*dzrhdi
       drhot2i = drhot2(ijk) + drhot(ijk)
       dxrhod(ijk) = drhot(ijk)*dxrhi + drhot2i*dxrhdi
       dyrhod(ijk) = drhot(ijk)*dyrhi + drhot2i*dyrhdi
       dzrhod(ijk) = drhot(ijk)*dzrhi + drhot2i*dzrhdi
       strgga(1,1)=strgga(1,1)+dxrho(ijk)*dxrhi+dxrhod(ijk)*dxrhdi
       strgga(2,2)=strgga(2,2)+dyrho(ijk)*dyrhi+dyrhod(ijk)*dyrhdi
       strgga(3,3)=strgga(3,3)+dzrho(ijk)*dzrhi+dzrhod(ijk)*dzrhdi
       strgga(2,3)=strgga(2,3)+dyrho(ijk)*dzrhi+dyrhod(ijk)*dzrhdi
       strgga(3,1)=strgga(3,1)+dzrho(ijk)*dxrhi+dzrhod(ijk)*dxrhdi
       strgga(1,2)=strgga(1,2)+dxrho(ijk)*dyrhi+dxrhod(ijk)*dyrhdi
    end do
  else
    do ijk = 1, mshnod
       dxrhi  = dxrho(ijk)
       dyrhi  = dyrho(ijk)
       dzrhi  = dzrho(ijk)
       dxrhdi = dxrhod(ijk)
       dyrhdi = dyrhod(ijk)
       dzrhdi = dzrhod(ijk)
       drhot1i = drhot1(ijk) + drhot(ijk)
       dxrho(ijk) = drhot1i*dxrhi + drhot(ijk)*dxrhdi
       dyrho(ijk) = drhot1i*dyrhi + drhot(ijk)*dyrhdi
       dzrho(ijk) = drhot1i*dzrhi + drhot(ijk)*dzrhdi
       drhot2i = drhot2(ijk) + drhot(ijk)
       dxrhod(ijk) = drhot(ijk)*dxrhi + drhot2i*dxrhdi
       dyrhod(ijk) = drhot(ijk)*dyrhi + drhot2i*dyrhdi
       dzrhod(ijk) = drhot(ijk)*dzrhi + drhot2i*dzrhdi
    end do
end if

call ddnrtg( nfile, myid, nodes,  &
& dxrho, dyrho, dzrho, mshnod, glx, gly, glz,  &
& thrhgr, nplw5ex, nplw5,  &
& mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, gx, gy, gz )


!   --- gradient of d_fxc(down)/d_dn ---
call ddnrtg( nfile, myid, nodes,  &
& dxrhod, dyrhod, dzrhod, mshnod, glx, gly, glz,  &
& thrhgr, nplw5ex, nplw5,  &
& mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, gx, gy, gz )


do ijk = 1, mshnod
   vlocud(ijk,1) = vlocud(ijk,1) - dxrho(ijk)
   vlocud(ijk,2) = vlocud(ijk,2) - dxrhod(ijk)
   rhoup = 0.5d0*( rho(ijk) + rhoud(ijk) )
   rhodn = 0.5d0*( rho(ijk) - rhoud(ijk) )
   evexc = evexc + rhoup*vlocud(ijk,1) + rhodn*vlocud(ijk,2)
end do


dbuf(1) = eexc
dbuf(2) = evexc
dbuf(3) = ecor
dbuf(4) = evcor
call gdsum(dbuf,4,dbufr)
eexc  = dbuf(1) * 2.d0 * rdelv
evexc = dbuf(2) * 2.d0 * rdelv
ecor  = dbuf(3) * 2.d0 * rdelv
evcor = dbuf(4) * 2.d0 * rdelv
do ijk = 1, mshnod
   vlocud(ijk,1) = 2.d0 * vlocud(ijk,1)
   vlocud(ijk,2) = 2.d0 * vlocud(ijk,2)
   eeexc(ijk)    = 2.d0 * eeexc(ijk)
end do


if( lcstress ) then
    dbuf(1) = strgga(1,1)
    dbuf(2) = strgga(2,2)
    dbuf(3) = strgga(3,3)
    dbuf(4) = strgga(2,3)
    dbuf(5) = strgga(3,1)
    dbuf(6) = strgga(1,2)
    call gdsum(dbuf,6,dbufr)
    strgga(1,1) = -dbuf(1) * 2.d0 * rdelv
    strgga(2,2) = -dbuf(2) * 2.d0 * rdelv
    strgga(3,3) = -dbuf(3) * 2.d0 * rdelv
    strgga(2,3) = -dbuf(4) * 2.d0 * rdelv
    strgga(3,1) = -dbuf(5) * 2.d0 * rdelv
    strgga(1,2) = -dbuf(6) * 2.d0 * rdelv
    strgga(3,2) = strgga(2,3)
    strgga(1,3) = strgga(3,1)
    strgga(2,1) = strgga(1,2)
end if


return
end




subroutine xclda( rho, rhocore, mshnod, rdelv,  &
&                 vexc, eeexc, eexc, evexc, ecor, evcor )
!-----------------------------------------------------------------------
! ( LDA )
!
!  ( output )
!    exchange & correlation   potential : vexc
!    exchange energy                    : eexc
!    exchange potential                 : evexc
!    correlation energy                 : ecor
!    correlation potential              : evcor
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
parameter( pi4 = 12.56637061436d0 )
dimension rho(mshnod)
dimension rhocore(mshnod)
dimension vexc(mshnod)
dimension eeexc(mshnod)
dimension dbuf(4), dbufr(4)


eexc  = 0.d0
evexc = 0.d0
ecor  = 0.d0
evcor = 0.d0
do m1 = 1, mshnod
   if( rhocore(m1).ge.0.d0 ) then
       rhot = rho(m1) + rhocore(m1)
     else
       rhot = rho(m1)
   end if
   if( rhot.lt.1.d-30 ) then
        vexc(m1) = 0.d0
       eeexc(m1) = 0.d0
     else
       rs = ( 3.D0/( pi4*rhot ) )**( 1.D0/3.D0 )
       call vxc( rs, ex, vx, ec, vc )
       vexc(m1)  = 2.d0*( vx + vc )
       eeexc(m1) = 2.d0*( ex + ec )
       eexc  = eexc  + ex * rhot
       evexc = evexc + vx * rho(m1)
       ecor  = ecor  + ec * rhot
       evcor = evcor + vc * rho(m1)
   end if
end do
dbuf(1) = eexc
dbuf(2) = evexc
dbuf(3) = ecor
dbuf(4) = evcor
call gdsum(dbuf,4,dbufr)
eexc  = dbuf(1) * 2.d0 * rdelv
evexc = dbuf(2) * 2.d0 * rdelv
ecor  = dbuf(3) * 2.d0 * rdelv
evcor = dbuf(4) * 2.d0 * rdelv


return
end




subroutine xcpbe( nfile, myid, nodes,  &
& vexc, eeexc, eexc, evexc, ecor, evcor, rho, prho, mshnod, rdelv,  &
& x, dxrho, dyrho, dzrho, thrhgr, nplw5ex, nplw5,  &
& glx, gly, glz, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, gx, gy, gz, t_comm,  &
& lcstress, strgga )
!-----------------------------------------------------------------------
! Generalized gradient corrected exchange correlation potential by PBE
!     for spin unpolarized case
!              ^^^^^^^^^^^
!-----------------------------------------------------------------------
use pbe_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: mshnod
real*8  :: vexc(mshnod), eeexc(mshnod), eexc, evexc, ecor, evcor
real*8  :: rho(mshnod), prho(mshnod), rdelv
real*8  :: x(*), dxrho(*), dyrho(*), dzrho(*)
real*8  :: thrhgr(*)
integer :: nplw5ex, nplw5
real*8  :: glx(*), gly(*), glz(*)
integer :: mftnod(*), mftdsp(*), mfd2ft(*), ntotfd, nd1vks(*), kfft0
real*8,  dimension(*) :: fft3x, fft3y
complex*16, dimension(*) :: fftwork
integer :: ngenh
integer :: ijkgd(-ngenh:ngenh)
real*8, dimension(0:ngenh) :: gx, gy, gz
real*8  :: t_comm
logical :: lcstress
real*8  :: strgga(3,3)

!-----declare local variables
!integer, parameter :: c1b3 = 1.d0/3.d0
!integer, parameter :: c4b3 = 4.d0/3.d0
!integer, parameter :: c7b3 = 7.d0/3.d0
integer :: ir, ijk
real*8  :: dbuf(6), dbufr(6)
real*8  :: rx, rho3, rs, drho, adrho, adrbro, ex0, ex, s, s2, dfxdn, dfxddn
real*8  :: ul, fsbb, fst, fss, dfsds
real*8  :: rs1h, g0bb, g0b2, g0ln, g0rs1, g0rs, ecrsz, ec, dg0dr1, dg0drs
real*8  :: drsdn, decdn, t, t2, t4, dfcdn, dfcddn, h2, eecrsz, arsz, arsz2
real*8  :: h0bs, h0bb, h0b2, h0ln, h0trsz, ech, h0bsa, h0bba, h0bst, h0bbt
real*8  :: dh0c1, dh0c, dh0da, dh0dt, dadn, dtdn, dh0dn, dtddn, dh0ddn
real*8  :: fxcdn, fxcddn, fxcddx, fxcddy, fxcddz
integer :: myid_, nodes_, nkd_
!save rs0, ax, bx, um, &
!&    a0, alp01, bet01, bet02, bet03, bet04, p0,  &
!&    betapw, t0, h1
!
!data rs0 / 0.62035049090d0 /
!!--- parameter for exchange energy
!data ax / -0.738558766382022405884230032680836d0 /
!data bx /  0.16162045967d0 /
!data um / 0.2195149727645171d0 /
!!      data uk / 0.8040d0 /
!
!!--- parameter for correlation energy
!data a0, alp01, bet01, bet02, bet03, bet04, p0  &
!&  / 0.0621814d0, 0.21370d0,  &
!&    7.5957d0, 3.5876d0, 1.6382d0, 0.49294d0, 1.d0 /
!data betapw / 0.06672455060314922d0 /
!data t0 / 0.2519289703424d0 /
!data h1 / 0.03109069086965489503494086371273d0 /


if( lcstress ) then
    strgga(1:3,1:3) = 0.d0
end if

!  --- gradient of charge density ---
do ir = 1, mshnod
   if( prho(ir).ge.0.d0 ) then
       x(ir) = rho(ir) + prho(ir)
     else
       x(ir) = rho(ir)
   end if
end do


!--- unify charge density
!      call unifylc( nfile, myid, nodes,
!     & x, mshnod, glx, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks )

call graden( nfile, myid, nodes,  &
& dxrho, dyrho, dzrho, mshnod, x, glx, gly, glz,  &
& thrhgr, nplw5ex, nplw5,  &
& mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, gx, gy, gz )

eexc  = 0.d0
evexc = 0.d0
ecor  = 0.d0
evcor = 0.d0
!=======================================================================
do ijk = 1, mshnod

if( prho(ijk).ge.0.d0 ) then
    rx = rho(ijk) + prho(ijk)
  else
    rx = rho(ijk)
end if
if( rx.lt.1.d-30 ) then
    vexc(ijk)  = 0.d0
    dxrho(ijk) = 0.d0
    dyrho(ijk) = 0.d0
    dzrho(ijk) = 0.d0
    eeexc(ijk) = 0.d0
else

rho3 = rx**c1b3
rs   = rs0/rho3

drho = dxrho(ijk)*dxrho(ijk) + dyrho(ijk)*dyrho(ijk)  &
&    + dzrho(ijk)*dzrho(ijk)
adrho = sqrt(drho)
adrbro = adrho / rx

!--- exchange terms ----------------------------------------------------
ex0  = ax * rho3

s  = bx * adrbro / rho3
if( s.lt.1.d-15 ) then
    ex = ex0
    dfxdn  = c4b3*ex0
    dfxddn = 0.d0
  else
    s2 = s*s

    ul = um/uk
    fsbb = 1.d0 + s2*ul
    fsbb = 1.d0/fsbb
    fst  = - uk*fsbb
    fss  = 1.d0 + uk + fst

!     --- exchange energy ---
    ex = ex0*fss


    dfsds  = - 2.d0*ul*s * fst*fsbb

    dfxdn  = c4b3*ex0 * ( fss - s*dfsds )
    dfxddn = ax * bx * dfsds / adrho
end if

!--- correlation terms -------------------------------------------------

rs1h = dsqrt(rs)

g0bb = rs1h*( bet01 + rs1h*( bet02 + rs1h*( bet03 + rs1h*bet04 )))
g0b2 = 1.d0/(a0*g0bb)
if( abs(g0b2).gt.1.d-10 ) then
    g0ln = log( 1.d0 + g0b2 )
  else
    g0ln = g0b2
end if
g0rs1 = 1.d0 + alp01*rs
g0rs  = -a0*g0rs1*g0ln

ecrsz = g0rs
!  --- correlation energy ---
ec    = ecrsz


   dg0dr1 = 0.5d0*bet01/rs1h + bet02  &
&         + rs1h*(1.5d0*bet03 + 2.d0*bet04*rs1h)
   dg0drs = -a0*alp01*g0ln  &
&         + g0rs1*dg0dr1 / ( ( 1.d0 + g0b2 )*g0bb*g0bb )

   drsdn = - rs*c1b3

   decdn = drsdn * dg0drs

!==== gradient terms ================
t  = t0*adrbro
t2 = t*t / rho3
if( t2.lt.1.d-30 ) then
    dfcdn  = ec + decdn
    dfcddn = 0.d0
else

h2 = betapw/h1

eecrsz = exp(-ecrsz/h1)
arsz   = h2 / ( eecrsz - 1.d0 )
arsz2  = arsz*arsz

t4 = t2*t2

h0bs =             t2 + arsz *t4 
h0bb = 1.d0 + arsz*t2 + arsz2*t4 
h0b2 = h0bs/h0bb
if( abs(h0b2).gt.1.d-10 ) then
    h0ln = log( 1.d0 + h2*h0b2 )
  else
    h0ln = h2*h0b2
end if
h0trsz = h1*h0ln

ech  = h0trsz

!  --- correlation energy by gradient terms ---
ec = ec + ech


   h0bsa =                t4
   h0bba = t2 + 2.d0*arsz*t4
   h0bst = 1.d0 + 2.d0*arsz *t2
   h0bbt = arsz + 2.d0*arsz2*t2

   dh0c1 = h2*h0bs + h0bb
   dh0c  = betapw/dh0c1
   dh0da = dh0c*( h0bsa - h0b2*h0bba )
   dh0dt = dh0c*( h0bst - h0b2*h0bbt )

   dadn = eecrsz*arsz2/betapw * decdn
   dtdn = - c7b3 * t2

   dh0dn = dh0da*dadn + dh0dt*dtdn


   dfcdn = ec + decdn + dh0dn


   dtddn = 2.d0*t2/adrbro

   dh0ddn = dh0dt*dtddn

   dfcddn = dh0ddn / adrho
end if
!-----------------------------------------------------------------------

!        exc(ijk) = ex  + ec
  fxcdn  = dfxdn  + dfcdn
  fxcddn = dfxddn + dfcddn

  eexc = eexc + rx*ex
  ecor = ecor + rx*ec
!      --- for check of consistency with total energy ---
!        evexc = evexc + rho(ijk)*fxcdn

  vexc(ijk) = fxcdn
  eeexc(ijk) = ex + ec


  if( lcstress ) then
      fxcddx = fxcddn*dxrho(ijk)
      fxcddy = fxcddn*dyrho(ijk)
      fxcddz = fxcddn*dzrho(ijk)
      strgga(1,1) = strgga(1,1) + fxcddx*dxrho(ijk)
      strgga(2,2) = strgga(2,2) + fxcddy*dyrho(ijk)
      strgga(3,3) = strgga(3,3) + fxcddz*dzrho(ijk)
      strgga(1,2) = strgga(1,2) + fxcddx*dyrho(ijk)
      strgga(2,3) = strgga(2,3) + fxcddy*dzrho(ijk)
      strgga(3,1) = strgga(3,1) + fxcddz*dxrho(ijk)
      dxrho(ijk) = fxcddx
      dyrho(ijk) = fxcddy
      dzrho(ijk) = fxcddz
    else
      dxrho(ijk) = fxcddn*dxrho(ijk)
      dyrho(ijk) = fxcddn*dyrho(ijk)
      dzrho(ijk) = fxcddn*dzrho(ijk)
  end if

end if
end do
!=======================================================================

call ddnrtg( nfile, myid, nodes,  &
& dxrho, dyrho, dzrho, mshnod, glx, gly, glz,  &
& thrhgr, nplw5ex, nplw5,  &
& mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, gx, gy, gz )

    do ijk = 1, mshnod
       vexc(ijk) = vexc(ijk) - dxrho(ijk)
       evexc = evexc + rho(ijk)*vexc(ijk)
    end do

dbuf(1) = eexc
dbuf(2) = evexc
dbuf(3) = ecor
dbuf(4) = evcor
call gdsum(dbuf,4,dbufr)


eexc  = dbuf(1) * 2.d0 * rdelv
evexc = dbuf(2) * 2.d0 * rdelv
ecor  = dbuf(3) * 2.d0 * rdelv
evcor = dbuf(4) * 2.d0 * rdelv
do ijk = 1, mshnod
    vexc(ijk) = 2.d0 * vexc(ijk)
   eeexc(ijk) = 2.d0 * eeexc(ijk)
end do


if( lcstress ) then
    dbuf(1) = strgga(1,1)
    dbuf(2) = strgga(2,2)
    dbuf(3) = strgga(3,3)
    dbuf(4) = strgga(2,3)
    dbuf(5) = strgga(3,1)
    dbuf(6) = strgga(1,2)
    call gdsum(dbuf,6,dbufr)
    strgga(1,1) = -dbuf(1) * 2.d0 * rdelv
    strgga(2,2) = -dbuf(2) * 2.d0 * rdelv
    strgga(3,3) = -dbuf(3) * 2.d0 * rdelv
    strgga(2,3) = -dbuf(4) * 2.d0 * rdelv
    strgga(3,1) = -dbuf(5) * 2.d0 * rdelv
    strgga(1,2) = -dbuf(6) * 2.d0 * rdelv
    strgga(3,2) = strgga(2,3)
    strgga(1,3) = strgga(3,1)
    strgga(2,1) = strgga(1,2)
end if


return
end




subroutine graden( nfile, myid, nodes,  &
& dxrho, dyrho, dzrho, mshnod, x, glx, gly, glz,  &
& thrhgr, nplw5ex, nplw5,  &
& mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, gx, gy, gz )
!-----------------------------------------------------------------------
!     gradient of charge density
!
!  input
!      glx(ir) : charge density
!
!  output
!     ( dxrho(ir), dyrho(ir), dzrho(ir) ) : gradient of charge density
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes
integer :: mshnod
real*8,  dimension(*) :: dxrho, dyrho, dzrho
real*8,  dimension(*) :: x
real*8,  dimension(*) :: glx, gly, glz
real*8,  dimension(*) :: thrhgr
integer :: nplw5ex, nplw5
integer, dimension(*) :: mftnod, mftdsp, mfd2ft
integer :: ntotfd, kfft0
integer, dimension(3) :: nd1vks
real*8,  dimension(*) :: fft3x, fft3y
complex*16, dimension(*) :: fftwork
integer :: ngenh
integer, dimension(-ngenh:ngenh) :: ijkgd
real*8, dimension(0:ngenh) :: gx, gy, gz


!--- gather charge density to node 0
call dgatherv(x,mshnod,glx,mftnod,mftdsp,0)

node1 = min( 1, nodes-1 )
node2 = min( 2, nodes-1 )


!   --- constants for fft ---
nml = 1
inv = 2
kfft1 = nd1vks(1)
kfft2 = nd1vks(2)
kfft3 = nd1vks(3)


if( myid==0 ) then
!   --- rho(r) -->  rho(g) ---
   do ik=1, kfft0
      fft3x(ik) = 0.d0
      fft3y(ik) = 0.d0
   end do
   do ifd = 1, ntotfd
      ift = mfd2ft(ifd)
      fft3x(ift) = glx(ifd)
   end do

   call rfft3( nml, fft3x, fft3y, fftwork,  &
&                   kfft1, kfft2, kfft3, kfft0, ierrft )

   do ig = 0, nplw5
      ijk  = ijkgd(ig)
      thrhgr(2*ig+1) = fft3x(ijk)
      thrhgr(2*ig+2) = fft3y(ijk)
   end do
end if

if( myid==0 ) then
    do idist = 1, node2
       call cdsend(100+idist,thrhgr,2*nplw5+2,idist,0)
    end do
  else
    do idist = 1, node2
       if( myid==idist ) then
           call cdrecv(100+idist,thrhgr,2*nplw5+2,0)
       end if
    end do
end if


!   --- gradient of density ---
! - x -
node0if: if( myid.eq.0 ) then
do ik=1, kfft0
   fft3x(ik) = 0.d0
   fft3y(ik) = 0.d0
end do
   ig = 0
   ijk = ijkgd(ig)
   fft3x(ijk) = 0.d0
   fft3y(ijk) = 0.d0
#ifndef VECTOR
do ig = 1, nplw5
   ijk  = ijkgd(ig)
   ijkm = ijkgd(-ig)
   fft3x(ijk)  =  -gx(ig)*thrhgr(2*ig+2)
   fft3y(ijk)  =   gx(ig)*thrhgr(2*ig+1)
   fft3x(ijkm) =  fft3x(ijk)
   fft3y(ijkm) = -fft3y(ijk)
end do
#else
do ig = 1, nplw5
   ijk  = ijkgd(ig)
   fft3x(ijk)  =  -gx(ig)*thrhgr(2*ig+2)
   fft3y(ijk)  =   gx(ig)*thrhgr(2*ig+1)
end do
do ig = 1, nplw5
   ijk = ijkgd(-ig)
   fft3x(ijk)  =  -gx(ig)*thrhgr(2*ig+2)
   fft3y(ijk)  =  -gx(ig)*thrhgr(2*ig+1)
end do
#endif
!     --- i*gx*rho(g) -->  d rho(r)/dx ---
call rfft3( inv, fft3x, fft3y, fftwork,  &
&                kfft1, kfft2, kfft3, kfft0, ierrft )

!     --- convert FFT -> FD meshes ---
do ifd = 1, ntotfd
   ift = mfd2ft(ifd)
   glx(ifd) = fft3x(ift)
end do
else node0if
do ifd = 1, ntotfd
   glx(ifd) = 0.d0
end do
end if node0if

!     --- store local variables
!      call distlc( nfile, myid, nodes,
!     & dxrho, mshnod, glx, mftdsp )


! - y -
node1if: if( myid.eq.node1 ) then
do ik=1, kfft0
   fft3x(ik) = 0.d0
   fft3y(ik) = 0.d0
end do
   ig = 0
   ijk = ijkgd(ig)
   fft3x(ijk) = 0.d0
   fft3y(ijk) = 0.d0
#ifndef VECTOR
do ig = 1, nplw5
   ijk  = ijkgd(ig)
   ijkm = ijkgd(-ig)
   fft3x(ijk)  =  -gy(ig)*thrhgr(2*ig+2)
   fft3y(ijk)  =   gy(ig)*thrhgr(2*ig+1)
   fft3x(ijkm) =  fft3x(ijk)
   fft3y(ijkm) = -fft3y(ijk)
end do
#else
do ig = 1, nplw5
   ijk  = ijkgd(ig)
   fft3x(ijk)  =  -gy(ig)*thrhgr(2*ig+2)
   fft3y(ijk)  =   gy(ig)*thrhgr(2*ig+1)
end do
do ig = 1, nplw5
   ijk = ijkgd(-ig)
   fft3x(ijk)  =  -gy(ig)*thrhgr(2*ig+2)
   fft3y(ijk)  =  -gy(ig)*thrhgr(2*ig+1)
end do
#endif
!     --- i*gy*rho(g) -->  d rho(r)/dy ---
call rfft3( inv, fft3x, fft3y, fftwork,  &
&                kfft1, kfft2, kfft3, kfft0, ierrft )

!     --- convert FFT -> FD meshes ---
do ifd = 1, ntotfd
   ift = mfd2ft(ifd)
   gly(ifd) = fft3x(ift)
end do
else node1if
do ifd = 1, ntotfd
   gly(ifd) = 0.d0
end do
end if node1if

!     --- store local variables
!      call distlc( nfile, myid, nodes,
!     & dyrho, mshnod, gly, mftdsp )


! - z -
node2if: if( myid.eq.node2 ) then
do ik=1, kfft0
   fft3x(ik) = 0.d0
   fft3y(ik) = 0.d0
end do
   ig = 0
   ijk = ijkgd(ig)
   fft3x(ijk) = 0.d0
   fft3y(ijk) = 0.d0
#ifndef VECTOR
do ig = 1, nplw5
   ijk  = ijkgd(ig)
   ijkm = ijkgd(-ig)
   fft3x(ijk)  =  -gz(ig)*thrhgr(2*ig+2)
   fft3y(ijk)  =   gz(ig)*thrhgr(2*ig+1)
   fft3x(ijkm) =  fft3x(ijk)
   fft3y(ijkm) = -fft3y(ijk)
end do
#else
do ig = 1, nplw5
   ijk  = ijkgd(ig)
   fft3x(ijk)  =  -gz(ig)*thrhgr(2*ig+2)
   fft3y(ijk)  =   gz(ig)*thrhgr(2*ig+1)
end do
do ig = 1, nplw5
   ijk = ijkgd(-ig)
   fft3x(ijk)  =  -gz(ig)*thrhgr(2*ig+2)
   fft3y(ijk)  =  -gz(ig)*thrhgr(2*ig+1)
end do
#endif
!     --- i*gz*rho(g) -->  d rho(r)/dz ---
call rfft3( inv, fft3x, fft3y, fftwork,  &
&                kfft1, kfft2, kfft3, kfft0, ierrft )

!     --- convert FFT -> FD meshes ---
do ifd = 1, ntotfd
   ift = mfd2ft(ifd)
   glz(ifd) = fft3x(ift)
end do
else node2if
do ifd = 1, ntotfd
   glz(ifd) = 0.d0
end do
end if node2if

!     --- store local variables
!      call distlc( nfile, myid, nodes,
!     & dzrho, mshnod, glz, mftdsp )

call dscatterv( glx,mftnod,mftdsp,dxrho,mshnod,0 )
call dscatterv( gly,mftnod,mftdsp,dyrho,mshnod,node1 )
call dscatterv( glz,mftnod,mftdsp,dzrho,mshnod,node2 )


return
end




subroutine ddnrtg( nfile, myid, nodes,  &
& dxrho, dyrho, dzrho, mshnod, glx, gly, glz,  &
& thrhgr, nplw5ex, nplw5,  &
& mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, gx, gy, gz )
!-----------------------------------------------------------------------
!     gradient of d_fxc/d_dn
!
!  input
!      dxrho, dyrho, dzrho
!
!  output
!      dxrho :   gradient of d_fxc/d_dn
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes
integer :: mshnod
real*8,  dimension(*) :: dxrho, dyrho, dzrho
real*8,  dimension(*) :: glx, gly, glz
real*8,  dimension(*) :: thrhgr
integer :: nplw5ex, nplw5
integer, dimension(*) :: mftnod, mftdsp, mfd2ft
integer :: ntotfd, kfft0
integer, dimension(3) :: nd1vks
real*8,  dimension(*) :: fft3x, fft3y
complex*16, dimension(*) :: fftwork
integer :: ngenh
integer, dimension(-ngenh:ngenh) :: ijkgd
real*8, dimension(0:ngenh) :: gx, gy, gz


node1 = min( 1, nodes-1 )
node2 = min( 2, nodes-1 )
!-----gather dxrho to key nodes
call dgatherv(dxrho,mshnod,glx,mftnod,mftdsp,0)
call dgatherv(dyrho,mshnod,gly,mftnod,mftdsp,node1)
call dgatherv(dzrho,mshnod,glz,mftnod,mftdsp,node2)


!   --- constants for fft ---
nml = 1
inv = 2
kfft1 = nd1vks(1)
kfft2 = nd1vks(2)
kfft3 = nd1vks(3)


do ig = 0, nplw5
   thrhgr(2*ig+1) = 0.d0
   thrhgr(2*ig+2) = 0.d0
end do

!   --- dxrho(r) -->  dxrho(g) ---

!--- unify
!      call unifylc( nfile, myid, nodes,
!     & dxrho, mshnod, glx, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks )

node0if: if( myid == 0 ) then
   do ik=1, kfft0
      fft3x(ik) = 0.d0
      fft3y(ik) = 0.d0
   end do
   do ifd = 1, ntotfd
      ift = mfd2ft(ifd)
      fft3x(ift) = glx(ifd)
   end do

   call rfft3( nml, fft3x, fft3y, fftwork,  &
&                   kfft1, kfft2, kfft3, kfft0, ierrft )

   do ig = 0, nplw5
      ijk  = ijkgd(ig)
      thrhgr(2*ig+1) = -gx(ig)*fft3y(ijk)
      thrhgr(2*ig+2) =  gx(ig)*fft3x(ijk)
   end do
end if node0if


!   --- dyrho(r) -->  dyrho(g) ---

!--- unify
!      call unifylc( nfile, myid, nodes,
!     & dyrho, mshnod, gly, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks )

node1if: if( myid == node1 ) then
   do ik=1, kfft0
      fft3x(ik) = 0.d0
      fft3y(ik) = 0.d0
   end do
   do ifd = 1, ntotfd
      ift = mfd2ft(ifd)
      fft3x(ift) = gly(ifd)
   end do

   call rfft3( nml, fft3x, fft3y, fftwork,  &
&                   kfft1, kfft2, kfft3, kfft0, ierrft )

   do ig = 0, nplw5
      ijk  = ijkgd(ig)
      thrhgr(2*ig+1) = thrhgr(2*ig+1) - gy(ig)*fft3y(ijk)
      thrhgr(2*ig+2) = thrhgr(2*ig+2) + gy(ig)*fft3x(ijk)
   end do
end if node1if


!   --- dzrho(r) -->  dzrho(g) ---

!--- unify
!      call unifylc( nfile, myid, nodes,
!     & dzrho, mshnod, glz, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks )

node2if: if( myid == node2 ) then
   do ik=1, kfft0
      fft3x(ik) = 0.d0
      fft3y(ik) = 0.d0
   end do
   do ifd = 1, ntotfd
      ift = mfd2ft(ifd)
      fft3x(ift) = glz(ifd)
   end do

   call rfft3( nml, fft3x, fft3y, fftwork,  &
&                   kfft1, kfft2, kfft3, kfft0, ierrft )

   do ig = 0, nplw5
      ijk  = ijkgd(ig)
      thrhgr(2*ig+1) = thrhgr(2*ig+1) - gz(ig)*fft3y(ijk)
      thrhgr(2*ig+2) = thrhgr(2*ig+2) + gz(ig)*fft3x(ijk)
   end do
end if node2if


!   --- gradient of d_fxc/d_dn in r-space ---
if( myid <= node2 ) then
do ik=1, kfft0
   fft3x(ik) = 0.d0
   fft3y(ik) = 0.d0
end do
   ig = 0
   ijk = ijkgd(ig)
   fft3x(ijk) = thrhgr(2*ig+1)
   fft3y(ijk) = thrhgr(2*ig+2)
do ig = 1, nplw5
   ijk  = ijkgd(ig)
   fft3x(ijk) =  thrhgr(2*ig+1)
   fft3y(ijk) =  thrhgr(2*ig+2)
end do
do ig = 1, nplw5
   ijkm = ijkgd(-ig)
   fft3x(ijkm) =  thrhgr(2*ig+1)
   fft3y(ijkm) = -thrhgr(2*ig+2)
end do

call rfft3( inv, fft3x, fft3y, fftwork,  &
&                kfft1, kfft2, kfft3, kfft0, ierrft )

!     --- convert FFT -> FD meshes ---
if( myid == 0 ) then
do ifd = 1, ntotfd
   ift = mfd2ft(ifd)
   glx(ifd) = fft3x(ift)
   gly(ifd) = 0.d0
   glz(ifd) = 0.d0
end do
else if( myid == node1 ) then
do ifd = 1, ntotfd
   ift = mfd2ft(ifd)
   gly(ifd) = fft3x(ift)
   glx(ifd) = 0.d0
   glz(ifd) = 0.d0
end do
else if( myid == node2 ) then
do ifd = 1, ntotfd
   ift = mfd2ft(ifd)
   glz(ifd) = fft3x(ift)
   glx(ifd) = 0.d0
   gly(ifd) = 0.d0
end do
end if
end if

call dscatterv( glx,mftnod,mftdsp,dxrho,mshnod,0 )
call dscatterv( gly,mftnod,mftdsp,dyrho,mshnod,node1 )
call dscatterv( glz,mftnod,mftdsp,dzrho,mshnod,node2 )
do i = 1, mshnod
   dxrho(i) = dxrho(i) + dyrho(i) + dzrho(i)
end do

!     --- store local variables
!      call distlc( nfile, myid, nodes,
!     & dxrho, mshnod, glx, mftdsp )


return
end




subroutine xcpbe_nscforce( nfile, myid, nodes,  &
& xexc, rho, prho, mshnod, rdelv,  &
& x, dxrho, dyrho, dzrho, delrho, ddelrho, thrhgr, nplw5ex, nplw5,  &
& glx, gly, glz, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, gx, gy, gz )
!-----------------------------------------------------------------------
! Generalized gradient corrected exchange correlation potential by PBE
!     for spin unpolarized case
!              ^^^^^^^^^^^
!-----------------------------------------------------------------------
use pbe_variables
implicit none
!parameter( c1b3 = 1.d0/3.d0 )
!parameter( c4b3 = 4.d0/3.d0 )
!parameter( c7b3 = 7.d0/3.d0 )
integer :: nfile(*), myid, nodes
integer :: mshnod
real*8,  dimension(mshnod) :: xexc, rho, prho, x, dxrho, dyrho, dzrho
real*8  :: delrho(mshnod), ddelrho(mshnod,3)
real*8  :: rdelv
integer :: nplw5ex, nplw5
real*8,  dimension(*) :: thrhgr, glx, gly, glz
integer, dimension(*) :: mftnod, mftdsp, mfd2ft, nd1vks
integer :: ntotfd, kfft0
real*8,  dimension(*) :: fft3x, fft3y
complex*16, dimension(*) :: fftwork
integer :: ngenh
integer, dimension(-ngenh:ngenh) :: ijkgd
real*8, dimension(0:ngenh) :: gx, gy, gz

!-----declare local variables
integer :: ir, ijk
real*8  :: rx, rho3, rs, drho, adrho, adrhor, adrbro, ex0, s, s2, axbx, rxr
real*8  :: ul, fsbb, fst, fss, dfsds0, dfsds, d2fsds2, dsdn, d2sdn2, dsddn, d2sdnddn
real*8  :: d2fxdn2, dfxddn, d2fxddn2, d2fxdnddn
real*8  :: dnx, dny, dnz, dotn, xcn1, xcn2

!save rs0, ax, bx, um, &
!&    a0, alp01, bet01, bet02, bet03, bet04, p0,  &
!&    betapw, t0, h1
!
!data rs0 / 0.62035049090d0 /
!!--- parameter for exchange energy
!data ax / -0.738558766382022405884230032680836d0 /
!data bx /  0.16162045967d0 /
!data um / 0.2195149727645171d0 /
!!      data uk / 0.8040d0 /
!
!!--- parameter for correlation energy
!data a0, alp01, bet01, bet02, bet03, bet04, p0  &
!&  / 0.0621814d0, 0.21370d0,  &
!&    7.5957d0, 3.5876d0, 1.6382d0, 0.49294d0, 1.d0 /
!data betapw / 0.06672455060314922d0 /
!data t0 / 0.2519289703424d0 /
!data h1 / 0.03109069086965489503494086371273d0 /


!  --- gradient of charge density ---
do ir = 1, mshnod
   if( prho(ir).ge.0.d0 ) then
       x(ir) = rho(ir) + prho(ir)
     else
       x(ir) = rho(ir)
   end if
end do


!--- unify charge density
!      call unifylc( nfile, myid, nodes,
!     & x, mshnod, glx, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks )

call graden( nfile, myid, nodes,  &
& dxrho, dyrho, dzrho, mshnod, x, glx, gly, glz,  &
& thrhgr, nplw5ex, nplw5,  &
& mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, gx, gy, gz )

call graden( nfile, myid, nodes,  &
& ddelrho(1,1), ddelrho(1,2), ddelrho(1,3), mshnod, delrho, glx, gly, glz,  &
& thrhgr, nplw5ex, nplw5,  &
& mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, gx, gy, gz )

axbx = ax * bx
!=======================================================================
do ijk = 1, mshnod

if( prho(ijk).ge.0.d0 ) then
    rx = rho(ijk) + prho(ijk)
  else
    rx = rho(ijk)
end if
if( rx.lt.1.d-30 ) then
    xexc(ijk)  = 0.d0
    dxrho(ijk) = 0.d0
    dyrho(ijk) = 0.d0
    dzrho(ijk) = 0.d0
else

rho3 = rx**c1b3
rs   = rs0/rho3

drho = dxrho(ijk)*dxrho(ijk) + dyrho(ijk)*dyrho(ijk)  &
&    + dzrho(ijk)*dzrho(ijk)
adrho = sqrt(drho)
rxr   = 1.d0 / rx
adrbro = adrho * rxr

!--- exchange terms ----------------------------------------------------
ex0  = ax * rho3

s  = bx * adrbro / rho3
if( s.lt.1.d-15 ) then
    d2fxdn2  = c4b3*ex0 * rxr * c1b3
    xexc(ijk)  = delrho(ijk) * d2fxdn2
    dxrho(ijk) = 0.d0
    dyrho(ijk) = 0.d0
    dzrho(ijk) = 0.d0
  else
    adrhor = 1.d0/adrho

    s2 = s*s

    ul = um/uk
    fsbb = 1.d0 + s2*ul
    fsbb = 1.d0/fsbb
    fst  = - uk*fsbb
    fss  = 1.d0 + uk + fst

    dfsds0  = - 2.d0*ul * fst*fsbb
    dfsds   = dfsds0 * s
    d2fsds2 = dfsds0 * ( 1.d0 - 4.d0*fsbb*ul*s2 )

    dsdn   = -c4b3*s * rxr
    d2sdn2 = -c7b3*dsdn * rxr
    dsddn  = bx / (rx*rho3)
    d2sdnddn = -c4b3*dsddn * rxr

    d2fxdn2  = ex0 * ( c4b3*(fss/(3.d0*rx) + 2.d0*dfsds*dsdn)  &
&                      + rx*(d2fsds2*dsdn*dsdn + dfsds*d2sdn2) )
    dfxddn   = axbx * dfsds * adrhor
    d2fxddn2 = axbx * d2fsds2*dsddn
    d2fxdnddn= ex0 * ( c4b3*dfsds*dsddn +  &
&                      rx*( d2fsds2*dsdn*dsddn + dfsds*d2sdnddn ) )

    dnx = dxrho(ijk) * adrhor
    dny = dyrho(ijk) * adrhor
    dnz = dzrho(ijk) * adrhor

    dotn = dnx*ddelrho(ijk,1) + dny*ddelrho(ijk,2) + dnz*ddelrho(ijk,3)

    xexc(ijk)  = delrho(ijk) * d2fxdn2 + d2fxdnddn * dotn

    xcn1 = delrho(ijk) * d2fxdnddn
    xcn2 = (d2fxddn2 - dfxddn)*dotn
    dxrho(ijk) = (xcn1 + xcn2)*dnx + dfxddn*ddelrho(ijk,1)
    dyrho(ijk) = (xcn1 + xcn2)*dny + dfxddn*ddelrho(ijk,2)
    dzrho(ijk) = (xcn1 + xcn2)*dnz + dfxddn*ddelrho(ijk,3)
end if

!--- correlation terms -------------------------------------------------

!rs1h = dsqrt(rs)
!
!g0bb = rs1h*( bet01 + rs1h*( bet02 + rs1h*( bet03 + rs1h*bet04 )))
!g0b2 = 1.d0/(a0*g0bb)
!if( abs(g0b2).gt.1.d-10 ) then
!    g0ln = log( 1.d0 + g0b2 )
!  else
!    g0ln = g0b2
!end if
!g0rs1 = 1.d0 + alp01*rs
!g0rs  = -a0*g0rs1*g0ln
!
!ecrsz = g0rs
!!  --- correlation energy ---
!ec    = ecrsz
!
!
!   dg0dr1 = 0.5d0*bet01/rs1h + bet02  &
!&         + rs1h*(1.5d0*bet03 + 2.d0*bet04*rs1h)
!   dg0drs = -a0*alp01*g0ln  &
!&         + g0rs1*dg0dr1 / ( ( 1.d0 + g0b2 )*g0bb*g0bb )
!
!   drsdn = - rs*c1b3
!
!   decdn = drsdn * dg0drs
!
!!==== gradient terms ================
!t  = t0*adrbro
!t2 = t*t / rho3
!if( t2.lt.1.d-30 ) then
!    dfcdn  = ec + decdn
!    dfcddn = 0.d0
!else
!
!h2 = betapw/h1
!
!eecrsz = exp(-ecrsz/h1)
!arsz   = h2 / ( eecrsz - 1.d0 )
!arsz2  = arsz*arsz
!
!t4 = t2*t2
!
!h0bs =             t2 + arsz *t4 
!h0bb = 1.d0 + arsz*t2 + arsz2*t4 
!h0b2 = h0bs/h0bb
!if( abs(h0b2).gt.1.d-10 ) then
!    h0ln = log( 1.d0 + h2*h0b2 )
!  else
!    h0ln = h2*h0b2
!end if
!h0trsz = h1*h0ln
!
!ech  = h0trsz
!
!!  --- correlation energy by gradient terms ---
!ec = ec + ech
!
!
!
!   h0bsa =                t4
!   h0bba = t2 + 2.d0*arsz*t4
!   h0bst = 1.d0 + 2.d0*arsz *t2
!   h0bbt = arsz + 2.d0*arsz2*t2
!
!   dh0c1 = h2*h0bs + h0bb
!   dh0c  = betapw/dh0c1
!   dh0da = dh0c*( h0bsa - h0b2*h0bba )
!   dh0dt = dh0c*( h0bst - h0b2*h0bbt )
!
!   dadn = eecrsz*arsz2/betapw * decdn
!   dtdn = - c7b3 * t2
!
!   dh0dn = dh0da*dadn + dh0dt*dtdn
!
!
!   dfcdn = ec + decdn + dh0dn
!
!
!
!   dtddn = 2.d0*t2/adrbro
!
!   dh0ddn = dh0dt*dtddn
!
!   dfcddn = dh0ddn / adrho
!end if
!-----------------------------------------------------------------------

end if
end do
!=======================================================================

call ddnrtg( nfile, myid, nodes,  &
& dxrho, dyrho, dzrho, mshnod, glx, gly, glz,  &
& thrhgr, nplw5ex, nplw5,  &
& mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, gx, gy, gz )

do ijk = 1, mshnod
   xexc(ijk) = xexc(ijk) - dxrho(ijk)
end do

do ijk = 1, mshnod
    xexc(ijk) = 2.d0 * xexc(ijk)
end do


return
end




SUBROUTINE VXC( RS, EX, VX, EC, VC )
!-----------------------------------------------------------------------
!     one - electron potential 
!     energy units are a.u. ( = Ryd * 0.5 )
!-----------------------------------------------------------------------
IMPLICIT REAL*8 (A-H,O-Z)                                         
!S-- Slater exchange ---------------------------------------------------
!S    VXC   = -3.D0 * 0.6108870575D0 / RS * 0.5
!-----------------------------------------------------------------------
!--- Gunarsson-Lundqvist exchange correlation potential. ---------------
!      BETA  =  1.D0 + 0.0545D0 * RS * LOG( 1.D0 + 11.4D0/RS )
!      VXC   = - 0.6108870575D0*BETA / RS
!-----------------------------------------------------------------------

DATA  A1 / -0.4582D0 /
DATA  B1,B2,B3 / -0.1423D0, 1.0529, 0.3334 /
DATA  C1,C2,C3,C4 / -0.0480D0, 0.0311, -0.0116, 0.0020 /
!      DATA  BB / 0.0140D0 /
!--- Ceperly and Alder as parametrized by Perdew and Zunger ------------
!       EX    ...... exchange energy
!       VX    ...... exchange potential ( non-relativistic )
!       RCORR ...... relativistic correction fot exchange potential
!       EC    ...... correlation energy
!       VC    ...... correlation potential
!-----------------------------------------------------------------------
EX    =  A1/RS
VX    =  4.0D0*EX/3.0D0
!C      BETA  =  BB/RS
!C      YETA  =  SQRT( 1.0D0 + BETA*BETA )
!C      ZETA  =  LOG( BETA + YETA )
!C      RCORR =  -0.5D0 + 1.5D0*ZETA/BETA/YETA
!C      VXR   =  VX * RCORR
IF( RS.GE.1.0D0 ) THEN
    BUNB = 1.0D0 + B2*SQRT(RS) + B3*RS
    EC   = B1/BUNB
    VC   = EC*(1.0D0 + (B2*SQRT(RS)/6.0D0 + B3*RS/3.0D0)/BUNB )
ELSE
    EC = C1 + C2*LOG(RS) + C3*RS + C4*RS*LOG(RS)
    VC = EC - (C2 + (C3 + C4)*RS + C4*RS*LOG(RS))/3.0D0
end if
!C      VXC   = VXR + VC
!      VXC   = VX + VC
!-----------------------------------------------------------------------

RETURN
END




SUBROUTINE vxclsda( RS, zeta, EX, VX, EC, VC, iflag )
!-----------------------------------------------------------------------
!     one - electron potential ( LSDA )
!     energy units are a.u. ( = Ryd * 0.5 )
!
!     iflag  : 1 -> relativistic correction / else -> not
!       RCORRE...... relativistic correction for exchange energy
!       RCORR ...... relativistic correction for exchange potential
!-----------------------------------------------------------------------
IMPLICIT REAL*8 (A-H,O-Z)                                         
dimension vx(2), vc(2)
save A1,B1,B2,B3,C1,C2,C3,C4,Ap1,Bp1,Bp2,Bp3,Cp1,Cp2,Cp3,Cp4,BB,  &
&    fc1,fc2
!S-- Slater exchange ---------------------------------------------------
!S    VXC   = -3.D0 * 0.6108870575D0 / RS * 0.5
!-----------------------------------------------------------------------
!--- Gunarsson-Lundqvist exchange correlation potential. ---------------
!      BETA  =  1.D0 + 0.0545D0 * RS * LOG( 1.D0 + 11.4D0/RS )
!      VXC   = - 0.6108870575D0*BETA / RS
!-----------------------------------------------------------------------
DATA  A1 / -0.4582D0 /
DATA  B1,B2,B3 / -0.1423D0, 1.0529d0, 0.3334d0 /
DATA  C1,C2,C3,C4 / -0.0480D0, 0.0311d0, -0.0116d0, 0.0020d0 /
DATA  Ap1 / -0.5773D0 /
DATA  Bp1,Bp2,Bp3 / -0.0843D0, 1.3981d0, 0.2611d0 /
DATA  Cp1,Cp2,Cp3,Cp4 / -0.0269D0, 0.01555d0, -0.0048d0, 0.0007d0/
DATA  BB / 0.0140D0 /

data fc1, fc2 / 1.92366d0, 2.56488d0 /
!--- Ceperly and Alder as parametrized by Perdew and Zunger ------------
!       EX    ...... exchange energy
!       VX    ...... exchange potential ( non-relativistic )
!       RCORR ...... relativistic correction fot exchange potential
!       EC    ...... correlation energy
!       VC    ...... correlation potential
!-----------------------------------------------------------------------
EXu    =  A1/RS
VXu    =  4.0D0*EXu/3.0D0
IF( RS.GE.1.0D0 ) THEN
    sqrs = sqrt(rs)
    BUNB = 1.0D0 + B2*sqrs + B3*RS
    ECu   = B1/BUNB
    VCu   = ECu*(1.0D0 + (B2*sqrs/6.0D0 + B3*RS/3.0D0)/BUNB )
ELSE
    dlnrs = log(rs)
    ECu = C1 + C2*dlnrs + C3*RS + C4*RS*dlnrs
    VCu = ECu - (C2 + (C3 + C4)*RS + C4*RS*dlnrs)/3.0D0
end if
if( abs(zeta).lt.1.d-30 ) then
    EX    = exu
    VX(1) = vxu
    VX(2) = vxu
    EC    = ecu
    VC(1) = vcu
    VC(2) = vcu
 else
!-----------------------------------------------------------------------
    EXpp   =  Ap1/RS
    VXp    =  4.0D0*EXpp/3.0D0
    IF( RS.GE.1.0D0 ) THEN
!c             sqrs = sqrt(rs)
       BUNB = 1.0D0 + Bp2*sqrs + Bp3*RS
       ECp   = Bp1/BUNB
       VCp   = ECp*(1.0D0 + (Bp2*sqrs/6.0D0 + Bp3*RS/3.0D0)/BUNB )
    ELSE
!             dlnrs = log(rs)
       ECp = Cp1 + Cp2*dlnrs + Cp3*RS + Cp4*RS*dlnrs
       VCp = ECp - (Cp2 + (Cp3 + Cp4)*RS + Cp4*RS*dlnrs)/3.0D0
    end if

    zetp13 = ( 1.d0 + zeta )**(1.d0/3.d0)
    zetm13 = ( 1.d0 - zeta )**(1.d0/3.d0)
    zetp43 = ( 1.d0 + zeta )*zetp13
    zetm43 = ( 1.d0 - zeta )*zetm13
    fzeta  = fc1*( zetp43 + zetm43 - 2.d0 )
    dfzeta = fc2*( zetp13 - zetm13 )

    EX    = exu + fzeta*( expp - exu )
    vx1   =  fzeta*( vxp - vxu )
    vx2   = dfzeta*( expp - exu )
   VX(1) = vxu + vx1 + vx2*(  1.d0 - zeta )
   VX(2) = vxu + vx1 + vx2*( -1.d0 - zeta )

   EC    = ecu + fzeta*( ecp - ecu )
   vc1   =  fzeta*( vcp - vcu )
   vc2   = dfzeta*( ecp - ecu )
   VC(1) = vcu + vc1 + vc2*(  1.d0 - zeta )
   VC(2) = vcu + vc1 + vc2*( -1.d0 - zeta )

end if


!      if( iflag.eq.1 ) then
!c   --- relativistic corrections ---------------------------------------
!          betar =  BB/rs
!          YETA  =  SQRT( 1.0D0 + betar*betar )
!          ZETA  =  LOG( betar + YETA )
!          R2    =  ( YETA - ZETA/betar )/betar
!          RCORRE=  1.D0 - 1.5D0*R2*R2
!          RCORR =  -0.5D0 + 1.5D0*ZETA/betar/YETA
!
!          ex    = ex * RCORRE
!          vx(1) = vx(1) * RCORR
!          vx(2) = vx(2) * RCORR
!c   --------------------------------------------------------------------
!      end if

RETURN
END




subroutine normal_pbe( trho, drho, d2rho, ri, ex, vx, ec, vc, iflag )
!-----------------------------------------------------------------------
! Generalized gradient corrected exchange correlation potential by PBE
!     for spin unpolarized case
!              ^^^^^^^^^^^
! input
!     trho   : density
!     drho   : derivative of density          d trho/ dr
!     d2rho  : second derivative of density  d^2trho/ dr^2
!     iflag  : 1 -> relativistic correction / else -> not
!
!       RCORRE...... relativistic correction for exchange energy
!       RCORR ...... relativistic correction for exchange potential
!
! output
!     ex     : exchange energy
!     vx     : exchange potential
!     ec     : correlation energy
!     vc     : correlation potential
!
!   ***  energy units are hartree ( = Ryd * 0.5 )
!-----------------------------------------------------------------------
use pbe_variables
implicit real*8 (a-h,o-z)   
!parameter( c1b3 = 1.d0/3.d0 )
!parameter( c4b3 = 4.d0/3.d0 )
!parameter( c7b3 = 7.d0/3.d0 )
!save rs0, ax, bx, um, &
!&    a0, alp01, bet01, bet02, bet03, bet04, p0,  &
!&    betapw, BB, t0, h1
!
!!--- parameter for relativistic correction
!DATA  BB / 0.0140D0 /
!
!data rs0 / 0.62035049090d0 /
!!--- parameter for exchange energy
!data ax / -0.738558766382022405884230032680836d0 /
!data bx /  0.16162045967d0 /
!data um / 0.2195149727645171d0 /
!!      data uk / 0.8040d0 /
!
!!--- parameter for correlation energy
!data a0, alp01, bet01, bet02, bet03, bet04, p0  &
!&  / 0.0621814d0, 0.21370d0,  &
!&    7.5957d0, 3.5876d0, 1.6382d0, 0.49294d0, 1.d0 /
!data betapw / 0.06672455060314922d0 /
!data t0 / 0.2519289703424d0 /
!data h1 / 0.03109069086965489503494086371273d0 /


if( trho .lt. 1.d-30 ) then
    ex = 0.d0
    vx = 0.d0
    ec = 0.d0
    vc = 0.d0
    return
end if

rho3 = trho**c1b3
rs   = rs0/rho3


adrho = abs(drho)
adrbro = adrho / trho

!--- exchange terms ----------------------------------------------------
ex0  = ax * rho3

s  = bx * adrbro / rho3
if( s.lt.1.d-15 ) then
    ex = ex0
    vx = c4b3*ex0
  else
    s2 = s*s

    ul = um/uk
    fsbb = 1.d0 + s2*ul
    fsbb = 1.d0/fsbb
    fst  = - uk*fsbb
    fs   = 1.d0 + uk + fst

!     --- exchange energy ---
    ex = ex0*fs


    dfsds  = - 2.d0*ul*s * fst*fsbb

    dfxdn  = c4b3*ex0 * ( fs - s*dfsds )
    dfxddn = ax * bx * dfsds

!=== from here, for olny sherical symmetric case ===
    d2fsds = - 2.d0*ul*( fst + 2.d0*dfsds*s )*fsbb

    drbro  = drho /trho
    d2rbro = d2rho/drho
    dsdr = s*( d2rbro - c4b3 * drbro )

    dfxdr = ax*bx*d2fsds*dsdr

!     --- exchange potential ---
    vx = dfxdn - ( dfxdr + 2.d0*dfxddn/ri )*sign(1.d0,drho)

!          write(*,'(10e15.7)') ri, ex0
!          write(*,'(10e15.7)') ri, dfxdn, -dfxdr, dfxddn
!          write(*,'(10e15.7)') ri, fs, dfsds, d2fsds
!          write(*,'(10e15.7)') ri, trho, drho, d2rho
!          write(*,'(10e15.7)') ri, s, d2rho/drho, 4.d0*drho/(3.d0*trho)

!===       end  for olny sherical symmetric case ===

end if

if( iflag.eq.1 ) then
!   --- relativistic corrections ---------------------------------------
    betar =  BB/rs
    YETA  =  SQRT( 1.0D0 + betar*betar )
    ZETA  =  LOG( betar + YETA )
    R2    =  ( YETA - ZETA/betar )/betar
    RCORRE=  1.D0 - 1.5D0*R2*R2
    RCORR =  -0.5D0 + 1.5D0*ZETA/betar/YETA

    ex = ex * RCORRE
    vx = vx * RCORR
!   --------------------------------------------------------------------
end if




!--- correlation terms -------------------------------------------------

rs1h = dsqrt(rs)

g0bb = rs1h*( bet01 + rs1h*( bet02 + rs1h*( bet03 + rs1h*bet04 )))
g0b2 = 1.d0/(a0*g0bb)
if( abs(g0b2).gt.1.d-10 ) then
    g0ln = log( 1.d0 + g0b2 )
  else
    g0ln = g0b2
end if
g0rs1 = 1.d0 + alp01*rs
g0rs  = -a0*g0rs1*g0ln

ecrsz = g0rs
!  --- correlation energy ---
ec    = ecrsz


   dg0dr1 = 0.5d0*bet01/rs1h + bet02  &
&         + rs1h*(1.5d0*bet03 + 2.d0*bet04*rs1h)
   dg0drs = -a0*alp01*g0ln  &
&         + g0rs1*dg0dr1 / ( ( 1.d0 + g0b2 )*g0bb*g0bb )

   drsdn = - rs * c1b3

   decdn = drsdn * dg0drs


!==== gradient terms ================
t  = t0*adrbro
t2 = t*t / rho3
if( t2.lt.1.d-30 ) then
    vc = ec + decdn
else

h2 = betapw/h1

eecrsz = exp(-ecrsz/h1)
arsz   = h2 / ( eecrsz - 1.d0 )
arsz2  = arsz*arsz

t4 = t2*t2

h0bs =             t2 + arsz *t4 
h0bb = 1.d0 + arsz*t2 + arsz2*t4 
h0b2 = h0bs/h0bb
if( abs(h0b2).gt.1.d-10 ) then
    h0ln = log( 1.d0 + h2*h0b2 )
  else
    h0ln = h2*h0b2
end if
h0trsz = h1*h0ln

ech  = h0trsz

!  --- correlation energy by gradient terms ---
ec = ec + ech



   h0bsa =                t4
   h0bba = t2 + 2.d0*arsz*t4
   h0bst = 1.d0 + 2.d0*arsz *t2
   h0bbt = arsz + 2.d0*arsz2*t2

   dh0c1 = h2*h0bs + h0bb
   dh0c  = betapw/dh0c1
   dh0da = dh0c*( h0bsa - h0b2*h0bba )
   dh0dt = dh0c*( h0bst - h0b2*h0bbt )

   dadn = eecrsz*arsz2/betapw * decdn
   dtdn = - c7b3 * t2

   dh0dn = dh0da*dadn + dh0dt*dtdn


   dfcdn = ec + decdn + dh0dn



   dtddn = 2.d0*t2/adrbro

   dh0ddn = dh0dt*dtddn

   dfcddn = dh0ddn


!=== from here, for olny sherical symmetric case ===

   h0bstt = 2.d0*arsz
   h0bsta = 2.d0*t2
   h0bbtt = 2.d0*arsz2
   h0bbta = 1.d0 + 4.d0*arsz*t2

   dh0c2  = -dh0dt/dh0c1
   dh0dtt = dh0c2*( h2*h0bst + h0bbt )  &
&  + dh0c*( h0bstt  &
&        - ( h0bst*h0bbt + h0bs*h0bbtt - h0b2*h0bbt*h0bbt )/h0bb )
   dh0dat = dh0c2*( h2*h0bsa + h0bba )  &
&  + dh0c*( h0bsta  &
&        - ( h0bsa*h0bbt + h0bs*h0bbta - h0b2*h0bbt*h0bba )/h0bb )

   dtdr   = t2*( 2.d0*d2rbro - c7b3 * drbro )
   dadr   = dadn * drbro

   dh0drt = dh0dtt*dtdr + dh0dat*dadr


   dtdrdn = dtddn*( d2rbro - c4b3 * drbro )

   dfcdr = dh0drt*dtddn + dh0dt*dtdrdn

!     --- correlation potential ---
    vc = dfcdn - ( dfcdr + 2.d0*dfcddn/ri )*sign(1.d0,drho)

!      if( iflag.eq.1 ) then
!      write(181,'(10e20.12)') t2, dfcdn, dh0dn, dh1dn,
!     &                        dfcdr, 2.d0*dfcddn/ri
!      write(181,'(20e14.6)') t2, 
!     &    dfcdr, dh0drt, dh1drt, dtddn,  dh0dt, dh1dt, dtdrdn,
!     &         dh1drs,  dsddn,           dh1ds,  dsdrdn
!      end if

!===       end  for olny sherical symmetric case ===

end if

return
end




