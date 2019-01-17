



module aspc_param
!-----------------------------------------------------------------------
! type declaration and initialization of variables for ASPC method
!-----------------------------------------------------------------------
implicit none

logical :: laspc_chg, laspc_wv                    ! for backward compatibility
logical :: lxlbomd
integer :: kxlbomd

save

end module




module aspc
!-----------------------------------------------------------------------
! type declaration and initialization of variables for ASPC method
!-----------------------------------------------------------------------
implicit none

integer :: kaspc    ! the total number of previous integration steps
real*8,  allocatable, dimension(:,:) :: coefaspc  ! for predictor
real*8,  allocatable, dimension(:)   :: omegaaspc ! for corrector
real*8,  allocatable, dimension(:,:,:) :: pspmat    ! <phi(t-j)|phi(t)>
real*8  :: aspc_corr  ! parameter for corrector


integer, parameter :: kprvmax = 9
real*8  :: dkappa, c0(kprvmax+1)

real*8  :: dkappa0(kprvmax) =  &
& (/2.d0, 1.53d0, 1.69d0, 1.75d0, 1.82d0, 1.84d0, 1.86d0, 1.88d0, 1.89d0 /)
real*8  :: alpha0(kprvmax) =  &
& (/0.d0, 190d0, 150d0, 57d0, 18d0, 5.5d0, 1.6d0, 0.44d0, 0.12d0 /)
real*8  :: c00(kprvmax+1,kprvmax) = reshape(  &
& (/    0d0,    0d0,    0d0,    0d0,    0d0,    0d0,    0d0,    0d0,    0d0,    0d0,  &
&      -2d0,    2d0,   -1d0,    0d0,    0d0,    0d0,    0d0,    0d0,    0d0,    0d0,  &
&      -2d0,    3d0,    0d0,   -1d0,    0d0,    0d0,    0d0,    0d0,    0d0,    0d0,  &
&      -3d0,    6d0,   -2d0,   -2d0,    1d0,    0d0,    0d0,    0d0,    0d0,    0d0,  &
&      -6d0,   14d0,   -8d0,   -3d0,    4d0,   -1d0,    0d0,    0d0,    0d0,    0d0,  &
&     -14d0,   36d0,  -27d0,   -2d0,   12d0,   -6d0,    1d0,    0d0,    0d0,    0d0,  &
&     -36d0,   99d0,  -88d0,   11d0,   32d0,  -25d0,    8d0,   -1d0,    0d0,    0d0,  &
&     -99d0,  286d0, -286d0,   78d0,   78d0,  -90d0,   42d0,  -10d0,    1d0,    0d0,  &
&    -286d0,  858d0, -936d0,  364d0,  168d0, -300d0,  184d0,  -63d0,   12d0,   -1d0   &
&  /), (/kprvmax+1,kprvmax/) )

save

end module




module xlbomd
!-----------------------------------------------------------------------
! type declaration and initialization of variables for XL-BOMD method
!-----------------------------------------------------------------------
implicit none

logical :: lwavxlbomd = .false.
logical :: lchgxlbomd = .false.

save

end module




subroutine set_coef_ASPC( nfile, myid, nodes, &
& ichest, ihest, laspc_chg_, laspc_wv_, nband, aspc_corr_, lxlbomd_, kxlbomd_ )
!-----------------------------------------------------------------------
!     Set coefficients for the ASPC method
!-----------------------------------------------------------------------
use outfile
use aspc_param
use aspc
implicit none

integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: ichest, ihest, nband
logical :: laspc_chg_, laspc_wv_
real*8  :: aspc_corr_
logical :: lxlbomd_
integer :: kxlbomd_

!-----declare local variables
integer :: k, j, i
real*8  :: a, b
integer :: status


laspc_chg = laspc_chg_
laspc_wv  = laspc_wv_

kaspc = max( ichest, ihest, 0 ) + 1

lxlbomd = lxlbomd_
kxlbomd = kxlbomd_
if( lxlbomd ) then
    dkappa = dkappa0(kxlbomd)
    c0(1:kxlbomd+1) = c00(1:kxlbomd+1,kxlbomd)*alpha0(kxlbomd)*1.d-03
    do i = 1, 2
       if(loutfile(i)) write(nfile(i),*) '*** Parameters for the XL-BOMD method'
       if(loutfile(i)) write(nfile(i),*) ' K =', kxlbomd
       if(loutfile(i)) write(nfile(i),*) ' kappa =', dkappa
       if(loutfile(i)) write(nfile(i),*) ' alpha(/1d-3) =', alpha0(kxlbomd)
       if(loutfile(i)) write(nfile(i),'(a,10f7.1)')  &
& '  c0 =', c0(1:kxlbomd+1)/(alpha0(kxlbomd)*1.d-03)
    end do
end if

!------allocate memory
allocate( coefaspc(kaspc,kaspc), stat=status )

coefaspc = 0.d0

if( laspc_chg .or. laspc_wv ) then

    !------allocate memory
    allocate( pspmat(nband,nband,0:kaspc-1), omegaaspc(kaspc), &
& stat=status )

!    if(loutfile(1)) write(nfile(1),*) '*** Coefficients for the ASPC method'
!    if(loutfile(2)) write(nfile(2),*) '*** Coefficients for the ASPC method'
    do k = 2, kaspc
    do j = 1, k
       a = 1.d0
       do i = k - j + 1, k - 1
          a = a * dble(i)
       end do
       b = 1.d0
       do i = k, k + j
          b = b * dble(i)
       end do
       coefaspc(j,k) = (-1)**(j+1) * dble(j*(2*k)*(2*k-1))*a/b
    
!       if(loutfile(1)) write(nfile(1),*) k,j,coefaspc(j,k)
!       if(loutfile(2)) write(nfile(2),*) k,j,coefaspc(j,k)
    end do
       omegaaspc(k) = dble(k)/dble(2*k-1)
!       if(loutfile(1)) write(nfile(1),*) k, omegaaspc(k)
!       if(loutfile(2)) write(nfile(2),*) k, omegaaspc(k)
    end do
end if

if( aspc_corr_ > 1.d-05 ) then
    aspc_corr = aspc_corr_
  else
    if( allocated(omegaaspc) ) then
        aspc_corr = omegaaspc( max( ihest,  1 ) + 1)
      else
        aspc_corr = 0.57d0
    end if
end if


return
end subroutine




subroutine set_mixprm_aspc( nfile, myid, nodes, &
& ncprvmx, aslh_aspc, bslh_aspc, aslh, bslh, laspc_exec )
!-----------------------------------------------------------------------
!     Set parameters for charge mixing in the ASPC method
!-----------------------------------------------------------------------
use aspc_param
use aspc
implicit none

integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: ncprvmx
real*8  :: aslh_aspc, bslh_aspc, aslh, bslh
logical :: laspc_exec

!-----declare local variables
!real*8  :: svaslh, svbslh
!logical :: licall = .true.
!save licall, svaslh, svbslh


!if( licall ) then
!    licall = .false.
!    svaslh = aslh
!    svbslh = bslh
!end if

if( laspc_exec ) then
    !---aslh = omegaaspc(ncprvmx+1)
    aslh = aslh_aspc
    bslh = bslh_aspc
!  else
!    aslh = svaslh
!    bslh = svbslh
end if

return
end subroutine




subroutine sbalin( nfile, myid, nodes,  &
& facc1, facc2, facc3, laspc_exec, t_comm, ltimecnt )
!-----------------------------------------------------------------------
!     estimate input wavefunctions by subspace alignment
!-----------------------------------------------------------------------
use param
use param_atom
use pwlda_atom
use pwlda_variables
use pwlda_proc
use pwlda_pw
use aspc
implicit none
integer :: nfile(*), myid, nodes
real*8  :: facc1, facc2, facc3
logical :: laspc_exec
real*8  :: t_comm
logical :: ltimecnt


!xlbomdif: if( lxlbomd ) then
!
!    if( laspc_exec ) then
!
!        !---Extended Lagrangian scheme XL-BOMD
!        call sbalinxlbomd( nfile, myid, nodes,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& gdcr, nplw, nplwex, npnod1, npnod2, npnod, nplcnt, npldsp,  &
!& iods, iodsg, ioa, ioag, idstnd, bufcr, prveig, nprvmx, nspnmx, prcd,  &
!& prod, prodr, prods, pdbuf,  &
!& dmtrxr, dmtrxi, bijr, biji, bm1r, bm1i, bx0r, bx0i, nbxxxx,  &
!& w1r, w1i, ew1, bzan, betk, tmpr, iblock,  &
!& node_c, node_r, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd,  &
!& t_comm, lspin, gdcrsv, nspnod, lcgjsv, rhcr, cgjr,  &
!& lvand, ntype, nhk1, nhk2, natom, nhk1_nat, nhk2_nat,  &
!& natnod1, natnod2, natnod, natcnt, natdsp,  &
!& lvandi, lking, lvand_r, lvand_g,  &
!& iatoit, mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, rvol, rdelv_c,  &
!& fft3x, fft3y, fftwork, ijkg, nplw2, eigr, eigi, ycos, ysin, nplwcs,  &
!& kxlbomd, dkappa, c0, ltimecnt )
!
!    else
!
!        if( laspc_wv ) then
!!----- To keep backward compatibility
!!      The following settings give the same results as sbalinold
!!        coefaspc(1,2) = 1.d0 + facc3
!!        coefaspc(2,2) = - facc3
!!        coefaspc(1,3) = 1.d0 + facc1
!!        coefaspc(2,3) = facc2  - facc1
!!        coefaspc(3,3) = - facc2
!
!            call sbalinaspc1( nfile, myid, nodes,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& gdcr, nplw, nplwex, npnod1, npnod2, npnod, nplcnt, npldsp,  &
!& iods, iodsg, ioa, ioag, idstnd, bufcr, keypwv, prveig, nprvmx, nspnmx, prcd,  &
!& ihest, nstep, nstep_ini, prod, prodr, prods, pdbuf, coefaspc, pspmat, kaspc,  &
!& dmtrxr, dmtrxi, bijr, biji, bm1r, bm1i, bx0r, bx0i, nbxxxx,  &
!& w1r, w1i, ew1, bzan, betk, tmpr, iblock,  &
!& node_c, node_r, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd,  &
!& t_comm, lspin, gdcrsv, nspnod, lcgjsv, rhcr, cgjr,  &
!& lvand, ntype, nhk1, nhk2, natom, nhk1_nat, nhk2_nat,  &
!& natnod1, natnod2, natnod, natcnt, natdsp,  &
!& lvandi, lking, lvand_r, lvand_g,  &
!& iatoit, mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, rvol, rdelv_c,  &
!& fft3x, fft3y, fftwork, ijkg, nplw2, eigr, eigi, ycos, ysin, nplwcs,  &
!& ltimecnt )
!
!        else
!
!            call sbalinold1( nfile, myid, nodes,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& gdcr, nplw, nplwex, npnod1, npnod2, npnod, nplcnt, npldsp,  &
!& iods, iodsg, ioa, ioag, idstnd, bufcr, keypwv, prveig, nprvmx, nspnmx, prcd,  &
!& ihest, nstep, nstep_ini, prod, prodr, prods, pdbuf, facc1, facc2, facc3,  &
!& dmtrxr, dmtrxi, bijr, biji, bm1r, bm1i, bx0r, bx0i, nbxxxx,  &
!& w1r, w1i, ew1, bzan, betk, tmpr, iblock,  &
!& node_c, node_r, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd,  &
!& t_comm, lspin, gdcrsv, nspnod, lcgjsv, rhcr, cgjr,  &
!& lvand, ntype, nhk1, nhk2, natom, nhk1_nat, nhk2_nat,  &
!& natnod1, natnod2, natnod, natcnt, natdsp,  &
!& lvandi, lking, lvand_r, lvand_g,  &
!& iatoit, mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, rvol, rdelv_c,  &
!& fft3x, fft3y, fftwork, ijkg, nplw2, eigr, eigi, ycos, ysin, nplwcs,  &
!& ltimecnt )
!
!        end if
!
!    end if
!
!else xlbomdif

!    if( laspc_wv ) then
!----- To keep backward compatibility
!      The following settings give the same results as sbalinold
!        coefaspc(1,2) = 1.d0 + facc3
!        coefaspc(2,2) = - facc3
!        coefaspc(1,3) = 1.d0 + facc1
!        coefaspc(2,3) = facc2  - facc1
!        coefaspc(3,3) = - facc2

        call sbalinaspc( nfile, myid, nodes,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, nplw, nplwex, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iods, iodsg, ioa, ioag, idstnd, bufcr, keypwv, prveig, nprvmx, nspnmx, prcd,  &
& ihest, nstep, nstep_ini, prod, prodr, prods, pdbuf, coefaspc, pspmat, kaspc,  &
& dmtrxr, dmtrxi, bijr, biji, bm1r, bm1i, bx0r, bx0i, nbxxxx,  &
& w1r, w1i, ew1, bzan, betk, tmpr, iblock,  &
& node_c, node_r, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd,  &
& t_comm, lspin, gdcrsv, nspnod, lcgjsv, rhcr, cgjr,  &
& lvand, ntype, nhk1, nhk2, natom, nhk1_nat, nhk2_nat,  &
& natnod1, natnod2, natnod, natcnt, natdsp,  &
& lvandi, lking, lvand_r, lvand_g,  &
& iatoit, mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, rvol, rdelv_c,  &
& fft3x, fft3y, fftwork, ijkg, nplw2, eigr, eigi, ycos, ysin, nplwcs,  &
& ltimecnt )

!    else
!
!        call sbalinold( nfile, myid, nodes,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& gdcr, nplw, nplwex, npnod1, npnod2, npnod, nplcnt, npldsp,  &
!& iods, iodsg, ioa, ioag, idstnd, bufcr, keypwv, prveig, nprvmx, nspnmx, prcd,  &
!& ihest, nstep, nstep_ini, prod, prodr, prods, pdbuf, facc1, facc2, facc3,  &
!& dmtrxr, dmtrxi, bijr, biji, bm1r, bm1i, bx0r, bx0i, nbxxxx,  &
!& w1r, w1i, ew1, bzan, betk, tmpr, iblock,  &
!& node_c, node_r, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd,  &
!& t_comm, lspin, gdcrsv, nspnod, lcgjsv, rhcr, cgjr,  &
!& lvand, ntype, nhk1, nhk2, natom, nhk1_nat, nhk2_nat,  &
!& natnod1, natnod2, natnod, natcnt, natdsp,  &
!& lvandi, lking, lvand_r, lvand_g,  &
!& iatoit, mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, rvol, rdelv_c,  &
!& fft3x, fft3y, fftwork, ijkg, nplw2, eigr, eigi, ycos, ysin, nplwcs,  &
!& ltimecnt )
!
!    end if

!end if xlbomdif


!-----Save predictor in the ASPC method
!if( laspc_exec ) then
!if( laspc_exec .and. .false. ) then !---not used
!    call save_pwv_in_aspc( nfile, myid, nodes,  &
!& gdcr, nplw, nplwex, npnod1, npnod2, npnod, nband, &
!& lspin, gdcrsv, nspnod, lcgjsv, nspnmx, lvand, lvandi, ntype, pgdcr, pgdcrsv )
!end if


return
end subroutine




subroutine sbalinaspc( nfile, myid, nodes,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, nplw, nplwex, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iods, iodsg, ioa, ioag, idstnd, bufcr, keypwv, prveig, nprvmx, nspnmx, prcd,  &
& ihest, nstep, nstep_ini, prod, prodr, prods, pdbuf, coefaspc, pspmat, kaspc,  &
& dmtrxr, dmtrxi, bijr, biji, bm1r, bm1i, bx0r, bx0i, nbxxxx,  &
& w1r, w1i, ew1, bzan, betk, tmpr, iblock,  &
& node_c, node_r, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd,  &
& t_comm, lspin, gdcrsv, nspnod, lcgjsv, rhcr, cgjr,  &
& lvand, ntype, nhk1, nhk2, natom, nhk1_nat, nhk2_nat,  &
& natnod1, natnod2, natnod, natcnt, natdsp,  &
& lvandi, lking, lvand_r, lvand_g,  &
& iatoit, mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, rvol, rdelv_c,  &
& fft3x, fft3y, fftwork, ijkg, nplw2, eigr, eigi, ycos, ysin, nplwcs,  &
& ltimecnt )
!-----------------------------------------------------------------------
!     estimate input wavefunctions by subspace alignment
!-----------------------------------------------------------------------
use outfile
use xlbomd
implicit none
integer :: nfile(*), myid, nodes
integer :: nband, nbnod1, nbnod2, nbnod, nbncnt(*), nbndsp(*)
integer :: nplw, nplwex, npnod1, npnod2, npnod, nplcnt(*), npldsp(*)
real*8  :: gdcr(2*npnod,*)
integer :: iods(*), iodsg(*), ioa(*), ioag(*), idstnd(*)
real*8  :: bufcr(*)
integer :: keypwv, nprvmx, nspnmx
real*8  :: prveig(2*npnod,nband,nprvmx,1:nspnmx)
real*8  :: prcd(2*npnod)
integer :: ihest, nstep, nstep_ini
real*8  :: prod(*), prodr(*), prods(*), bzan(*), betk(*), pdbuf(*)
integer :: kaspc
real*8  :: coefaspc(kaspc,kaspc), pspmat(nband,nband,0:kaspc-1)
integer :: nbxxxx
real*8  :: dmtrxr(nbxxxx,*), bijr(nbxxxx,*), bm1r(nbnod,*), bx0r(nbnod,*)
real*8  :: dmtrxi(*), biji(*), bm1i(*), bx0i(*)
real*8  :: w1r(*), w1i(*), ew1(*)
integer :: iblock
real*8  :: tmpr(iblock,*)
integer :: node_c, node_r
integer :: lbncnt(*), lbndsp(*), mbncnt(*), mbndsp(*), jdstnd(*)
real*8  :: t_comm
logical :: lspin, lcgjsv
real*8  :: gdcrsv(*)
integer :: nspnod
real*8  :: rhcr(*), cgjr(*)
logical :: lvand
integer :: ntype, nhk1(ntype), nhk2(ntype), natom
integer :: nhk1_nat(ntype), nhk2_nat(ntype)
integer :: natnod1, natnod2, natnod, natcnt(nodes), natdsp(nodes)
logical :: lvandi(ntype), lking(ntype), lvand_r, lvand_g
integer :: iatoit(natom)
integer :: mfd2ft_c(*), ntotfd_c, nd1vks_c(3), kfft0
real*8  :: rvol, rdelv_c
real*8  :: fft3x(*), fft3y(*)
complex*16 :: fftwork(*)
integer :: nplw2
integer :: ijkg(-nplw2:nplw2)
real*8  :: eigr(0:*), eigi(0:*)
real*8  :: ycos(*), ysin(*)
integer :: nplwcs
logical :: ltimecnt

!---declare local variables
!logical :: lsbalin
integer :: nnspin, nspin, i, ib, jb, iaspc
real*8  :: ctt0, ct0


!lsbalin = .true.

if( lspin ) then
    nnspin = 2
  else
    nnspin = 1
end if

if( lwavxlbomd ) then
!    do nspin = 1, nnspin
!       do i = 2, nprvmx
!          do ib = 1, nband
!             prveig(1:2*npnod,ib,i-1,nspin) = prveig(1:2*npnod,ib,i,nspin)
!          end do
!       end do
!
!       !----- save slmir
!       if( lvand ) then
!           call revsftprvslmir( nfile, myid, nodes, nspin )
!!           call setprvslmir( nfile, myid, nodes, 1, nspin )
!       end if
!    end do
!
!    keypwv = nprvmx - 1
!---The subspace rotation for the w.f. at the current time step to minimize the distance
!---to the w.f. at the previous steps seems to be the cause of a drift of the conserved quantity.
    keypwv = 0
    lwavxlbomd = .false.
end if


if( nstep.eq.nstep_ini+1 .or. keypwv.eq.0 ) then
    keypwv = 1
    do nspin = 1, nnspin
       if( lspin ) then
           call stspud( nspin, gdcr, gdcrsv, nspnod, lcgjsv )
!           if( lvand ) call ldslmi( nfile, myid, nodes, nspin, lvandi, .false. )
       end if

       !---set nspin in tddft-fssh.f90
       call set_nspin_in_tddft_fssh( nspin )

       !----- save wavefunctions
       do ib = 1, nband
          prveig(1:2*npnod,ib,1,nspin) = gdcr(1:2*npnod,ib)
       end do

       !----- save slmir
!       if( lvand ) then
!           call setprvslmir( nfile, myid, nodes, 1, nspin )
!       end if


!       if( lvand ) then
!
!           if( lvand_r ) then
!               !--- copy gdcr to rhcr
!               call cpgdrh( nfile, myid, nodes,  &
!& rhcr, gdcr, npnod, nband )
!               !--- to convert G decomposition to band decomposition
!               call gdtobd( nfile, myid, nodes, ctt0,  &
!& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& rhcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
!& iods, iodsg, idstnd, .false. )
!               call calslm_r( nfile, myid, nodes,  &
!& rhcr, nband, nbnod1, nbnod2, nbnod, nplw, nplwex,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking, iatoit,  &
!& mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, rvol, rdelv_c,  &
!& fft3x, fft3y, fftwork, ijkg, nplw2, eigr, eigi )
!
!               if( lvand_g ) then
!                   call calslm( nfile, myid, nodes,  &
!& rhcr, nband, nbnod1, nbnod2, nbnod, nplw, nplwex,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking, ycos, ysin, nplwcs )
!               end if
!
!               call slm_bdtogd( nfile, myid, nodes, ctt0,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& natom, natnod1, natnod2, natnod, natcnt, natdsp,  &
!& iodsg, ioag, idstnd, .false., .false. )
!
!           end if
!
!           if( lvand_g .and. .not.lvand_r ) then
!               call calslm_g( nfile, myid, nodes,  &
!& gdcr, nband, nplw, nplwex, npnod1, npnod2, npnod,  &
!& ntype, nhk1, nhk2, natom, natnod1, natnod2, natnod,  &
!& lvandi, lking, ycos, ysin, nplwcs, ioa )
!           end if
!
!           if(loutfile(1)) write(nfile(1),*) ' ( Gram-Schmidt )      '
!           if(loutfile(2)) write(nfile(2),*) ' ( Gram-Schmidt )      '
!!           call schmidt( nfile, myid, nodes, ct0,  &
!!& lvand, gdcr, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod,  &
!!& nbncnt, nbndsp, iods, iodsg, prod, prodr, pdbuf,  &
!!& .false., .false., .false., .false.,  &
!!& dmtrxr, nbxxxx,  &
!!& node_c, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, bijr, bm1r,  &
!!& rvol, ntype, nhk1_nat, nhk2_nat, natom, lvandi )
!           call schmidt( nfile, myid, nodes, ct0,  &
!& gdcr, iods, iodsg, .false., .false., .false., .false. )
!
!           if( lspin ) then
!               call svslmi( nfile, myid, nodes, nspin, lvandi, .false. )
!           end if
!       end if

    end do

    !----- save rhoij
!    if( lvand ) then
!        call setprvrhoij( nfile, myid, nodes, .true. )
!    end if

    return
end if


do i = 0, nprvmx
   if( ihest == i .or. nstep == nstep_ini+1+i .or. keypwv == i ) then
       iaspc = i + 1
       exit
   end if
end do

keypwv = min( keypwv + 1, nprvmx )

do nspin = 1, nnspin
   if( lspin ) then
       call stspud( nspin, gdcr, gdcrsv, nspnod, lcgjsv )
!       if( lvand ) call ldslmi( nfile, myid, nodes, nspin, lvandi, .false. )
   end if

!   if( lsbalin ) then
!       call roteig1( nfile, myid, nodes,  &
!& gdcr, prveig(1,1,1,nspin), nband, nbnod1, nbnod2, nbnod,  &
!& nbncnt, nbndsp, npnod1, npnod2, npnod, prod, prodr, prods,  &
!& lvand, rvol, nspin,  &
!& dmtrxr, dmtrxi, bijr, biji, bm1r, bm1i, bx0r, bx0i, nbxxxx,  &
!& w1r, w1i, ew1, bzan, betk, tmpr, iblock,  &
!& node_r, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, t_comm,  &
!& ltimecnt )
       call roteig( nfile, myid, nodes,  &
& gdcr, prveig(1,1,1,nspin), nband, nbnod1, nbnod2, nbnod,  &
& nbncnt, nbndsp, npnod1, npnod2, npnod, prod, prodr, prods,  &
& lvand, rvol, nspin, nprvmx, max(iaspc-1,1),  &
& dmtrxr, dmtrxi, bijr, biji, bm1r, bm1i, bx0r, bx0i, nbxxxx,  &
& w1r, w1i, ew1, bzan, betk, tmpr, iblock,  &
& node_r, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, t_comm, ltimecnt )
!   else
!       call aspcmat( nfile, myid, nodes,  &
!& gdcr, prveig(1,1,1,nspin), nband, nbnod1, nbnod2, nbnod,  &
!& nbncnt, nbndsp, npnod1, npnod2, npnod, prod, prodr,  &
!& lvand, rvol, nspin, nprvmx, pspmat, iaspc-1,  &
!& dmtrxr, nbxxxx, ltimecnt )
!   end if

   do ib = 1, nband
      if( iaspc > 1 ) then
          prcd(1:2*npnod) = coefaspc(1,iaspc)*gdcr(1:2*npnod,ib)

!        if( lsbalin ) then
          do i = 2, iaspc
             prcd(1:2*npnod) = prcd(1:2*npnod) + coefaspc(i,iaspc)*prveig(1:2*npnod,ib,i-1,nspin)
          end do
!        else
!          do i = 2, iaspc
!             gnk(1:2*npnod) = 0.d0
!             do jb = 1, nband
!                gnk(1:2*npnod) = gnk(1:2*npnod) + prveig(1:2*npnod,jb,i-1,nspin)*pspmat(ib,jb,i-1)
!             end do
!             prcd(1:2*npnod) = prcd(1:2*npnod) + coefaspc(i,iaspc)*gnk(1:2*npnod)
!          end do
!        end if
      end if

      do i = nprvmx - 1, 1, -1
         prveig(1:2*npnod,ib,i+1,nspin) = prveig(1:2*npnod,ib,i,nspin)
      end do
      prveig(1:2*npnod,ib,1,nspin) = gdcr(1:2*npnod,ib)

      if( iaspc > 1 ) then
          gdcr(1:2*npnod,ib) = prcd(1:2*npnod)
      end if
   end do

   !----- save slmir
!   if( lvand ) then
!       call sftprvslmir( nfile, myid, nodes, nspin )
!       call setprvslmir( nfile, myid, nodes, 1, nspin )
!   end if

end do


!----- save rhoij
!if( lvand ) then
!    call setprvrhoij( nfile, myid, nodes, .false. )
!end if


if( .not.lvand .and. ihest.eq.0 ) return

if(loutfile(1)) write(nfile(1),*) ' ( Gram-Schmidt )      '
if(loutfile(2)) write(nfile(2),*) ' ( Gram-Schmidt )      '
do nspin = 1, nnspin
   if( lspin ) then
       call stspud( nspin, gdcr, gdcrsv, nspnod, lcgjsv )
   end if

   !---set nspin in tddft-fssh.f90
   call set_nspin_in_tddft_fssh( nspin )

!   if( lvand ) then
!
!       if( lvand_r ) then
!           !--- copy gdcr to rhcr
!           call cpgdrh( nfile, myid, nodes,  &
!& rhcr, gdcr, npnod, nband )
!           !--- to convert G decomposition to band decomposition
!           call gdtobd( nfile, myid, nodes, ctt0,  &
!& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& rhcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
!& iods, iodsg, idstnd, .false. )
!           call calslm_r( nfile, myid, nodes,  &
!& rhcr, nband, nbnod1, nbnod2, nbnod, nplw, nplwex,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking, iatoit,  &
!& mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, rvol, rdelv_c,  &
!& fft3x, fft3y, fftwork, ijkg, nplw2, eigr, eigi )
!
!           if( lvand_g ) then
!               call calslm( nfile, myid, nodes,  &
!& rhcr, nband, nbnod1, nbnod2, nbnod, nplw, nplwex,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking, ycos, ysin, nplwcs )
!           end if
!
!           call slm_bdtogd( nfile, myid, nodes, ctt0,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& natom, natnod1, natnod2, natnod, natcnt, natdsp,  &
!& iodsg, ioag, idstnd, .false., .false. )
!
!       end if
!
!       if( lvand_g .and. .not.lvand_r ) then
!           call calslm_g( nfile, myid, nodes,  &
!& gdcr, nband, nplw, nplwex, npnod1, npnod2, npnod,  &
!& ntype, nhk1, nhk2, natom, natnod1, natnod2, natnod,  &
!& lvandi, lking, ycos, ysin, nplwcs, ioa )
!       end if
!
!   end if

!   call schmidt( nfile, myid, nodes, ct0,  &
!& lvand, gdcr, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod,  &
!& nbncnt, nbndsp, iods, iodsg, prod, prodr, pdbuf,  &
!& .false., .false., .false., .false.,  &
!& dmtrxr, nbxxxx,  &
!& node_c, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, bijr, bm1r,  &
!& rvol, ntype, nhk1_nat, nhk2_nat, natom, lvandi )
   call schmidt( nfile, myid, nodes, ct0,  &
& gdcr, iods, iodsg, .false., .false., .false., .false. )

!   if( lspin .and. lvand ) then
!       call svslmi( nfile, myid, nodes, nspin, lvandi, .false. )
!   end if

end do


return
end




subroutine roteig( nfile, myid, nodes,  &
& cgjr, rhcr, nband_, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& npnod1, npnod2, npnod, prod, prodr, prods,  &
& lvand, rvol, nspin, nprvmx, kpwvx,  &
& aijr, aiji, bijr, biji, bm1r, bm1i, bx0r, bx0i, nbxxxx,  &
& w1r, w1i, ew1, bzan, betk, tmpr, iblock_,  &
& node_r_, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, t_comm, ltimecnt )
!-----------------------------------------------------------------------
!     subspace alignment by  cgjr and rhcr
!-----------------------------------------------------------------------
use param
implicit none
integer :: nfile(*), myid, nodes
integer :: nband_, nbnod1, nbnod2, nbnod, nbncnt(*), nbndsp(*)
integer :: npnod1, npnod2, npnod
integer :: nprvmx, kpwvx
real*8  :: cgjr(2*npnod,*), rhcr(2*npnod,nband,nprvmx)
real*8  :: prod(*), prodr(*), prods(*)
logical :: lvand
real*8  :: rvol
integer :: nspin
integer :: nbxxxx
real*8  :: aijr(nbxxxx,*), bijr(nbxxxx,*), bm1r(nbnod,*), bx0r(nbnod,*)
real*8  :: aiji(nbxxxx,*), biji(nbxxxx,*), bm1i(nbnod,*), bx0i(nbnod,*)
real*8  :: w1r(*), w1i(*), ew1(*), bzan(*), betk(*)
integer :: iblock_
real*8  :: tmpr(iblock,*)
integer :: node_r_
integer :: lbncnt(*), lbndsp(*), mbncnt(*), mbndsp(*), jdstnd(*)
real*8  :: t_comm
logical :: ltimecnt


!if( .not.lnoncollinear ) then
    call roteig2( nfile, myid, nodes,  &
& cgjr, rhcr, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& npnod1, npnod2, npnod, prod, prodr,  &
& lvand, rvol, nspin, nprvmx, kpwvx,  &
& aijr, bijr, bm1r, bx0r, nbxxxx, w1r, w1i, ew1, tmpr, iblock,  &
& node_r, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, t_comm, ltimecnt )
!else
!    !-----noncollinear magnetism
!    call roteig_k2( nfile, myid, nodes,  &
!& cgjr, rhcr, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& npnod1, npnod2, npnod, prod, prodr, prods,  &
!& lvand, rvol, nspin, nprvmx, kpwvx,  &
!& aijr, aiji, bijr, biji, bm1r, bm1i, bx0r, bx0i, nbxxxx,  &
!& w1r, w1i, ew1, bzan, betk, tmpr, iblock,  &
!& node_r, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, t_comm, ltimecnt )
!end if


return
end




subroutine roteig2( nfile, myid, nodes,  &
& cgjr, rhcr, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& npnod1, npnod2, npnod, prod, prodr,  &
& lvand, rvol, nspin, nprvmx, kpwvx,  &
& aijr, bijr, bm1r, bx0r, nbxxxx, w1r, w1i, ew1, tmpr, iblock,  &
& node_r, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, t_comm, ltimecnt )
!-----------------------------------------------------------------------
!     subspace alignment by  cgjr and rhcr
!-----------------------------------------------------------------------
use outfile
implicit none
integer :: nfile(*), myid, nodes
integer :: nband, nbnod1, nbnod2, nbnod, nbncnt(*), nbndsp(*)
integer :: npnod1, npnod2, npnod
integer :: nprvmx, kpwvx
real*8  :: cgjr(2*npnod,*), rhcr(2*npnod,nband,nprvmx)
real*8  :: prod(*), prodr(*)
logical :: lvand
real*8  :: rvol
integer :: nspin
integer :: nbxxxx
real*8  :: aijr(nbxxxx,*), bijr(nbxxxx,*), bm1r(nbnod,*), bx0r(nbnod,*)
real*8  :: w1r(*), w1i(*), ew1(*)
integer :: iblock
real*8  :: tmpr(iblock,*)
integer :: node_r
integer :: lbncnt(*), lbndsp(*), mbncnt(*), mbndsp(*), jdstnd(*)
real*8  :: t_comm
logical :: ltimecnt

!------declare local variables
character(5) :: checkc
integer :: i, j, k, m, mshnod, root, ifvec, ifsrt, ib, jb
integer :: nmaxx, ncblck, icblck, icb1, icb2, keypwv
real*8  :: cscr, cosi, t1
logical :: lerror
real*8  :: tc(10), ct, ct0, timecnt
real*8, parameter :: epsimr = 1.0d-15


do i = 1, 10
   tc(i) = 0.d0
end do

if( ltimecnt ) then
    ct  = timecnt()
    ct0 = ct
end if

mshnod = 2*npnod

!!check
!      do i = 1, nband
!         do j = 1, nband
!            prod(j) = 0.d0
!            do m = 1, mshnod
!               prod(j) = prod(j) + rhcr(m,j,1)*rhcr(m,i,1)
!            end do
!            if( myid.eq.0 ) then
!                prod(j) = prod(j)*2.d0 - rhcr(1,j,1)*rhcr(1,i,1)
!              else
!                prod(j) = prod(j)*2.d0
!            end if
!         end do
!         if( lvand ) then
!             do j = 1, nband
!                call calpcspc( nfile, myid, nodes,
!     & cscr, i, j, rvol, 1, 1, nspin )
!                prod(j) = prod(j) + cscr
!             end do
!         end if
!         call gdsum(prod,nband,prodr)
!         lerror = .false.
!         errmx  = 0.d0
!         do j = 1, nband
!         if( i.eq.j ) prod(j) = prod(j) - 1.d0
!         lerror = lerror .or. abs(prod(j)).gt.1.d-10
!         errmx  = max( errmx, abs(prod(j)) )
!!         if( abs(prod(j)).gt.1.d-10 ) then
!!         if( myid.eq.0 ) then
!!             write(nfile(1),'(a39,2i5,e15.6)')
!!     &       'orthonormalization error in roteig:',i, j, prod(j)
!!         end if
!!         end if
!         end do
!         if( lerror ) then
!         if( myid.eq.0 ) then
!             write(nfile(1),'(a39,i5,e15.6)')
!     &       'orthonormalization error in roteig:',i, errmx
!         end if
!         end if
!      end do
!!check
!      do i = 1, nband
!         do j = 1, nband
!            prod(j) = 0.d0
!            do m = 1, mshnod
!               prod(j) = prod(j) + cgjr(m,j)*cgjr(m,i)
!            end do
!            if( myid.eq.0 ) then
!                prod(j) = prod(j)*2.d0 - cgjr(1,j)*cgjr(1,i)
!              else
!                prod(j) = prod(j)*2.d0
!            end if
!         end do
!         if( lvand ) then
!             do j = 1, nband
!                call calccscc( nfile, myid, nodes,
!     & cscr, i, j, rvol, nspin )
!                prod(j) = prod(j) + cscr
!             end do
!         end if
!         call gdsum(prod,nband,prodr)
!         lerror = .false.
!         errmx  = 0.d0
!         do j = 1, nband
!         if( i.eq.j ) prod(j) = prod(j) - 1.d0
!         lerror = lerror .or. abs(prod(j)).gt.1.d-10
!         errmx  = max( errmx, abs(prod(j)) )
!!         if( abs(prod(j)).gt.1.d-10 ) then
!!         if( myid.eq.0 ) then
!!             write(nfile(1),'(a39,2i5,e15.6)')
!!     &       'orthonormalization error in roteig:',i, j, prod(j)
!!         end if
!!         end if
!         end do
!         if( lerror ) then
!         if( myid.eq.0 ) then
!             write(nfile(1),'(a39,i5,e15.6)')
!     &       'orthonormalization error (2) in roteig:',i, errmx
!         end if
!         end if
!      end do
!!check
!--- calculate U-matrix : ( aijr, aiji ) -------------------------------
! --- This is serial code. ---
!      do ib=1, nband
!      do jb=1, nband
!         aijr(ib,jb) = 0.d0
!         do ig = 1, mshnod
!            aijr(ib,jb) = aijr(ib,jb) + cgjr(ig,jb)*rhcr(ig,ib)
!         end do
!      end do
!      end do
! ----------------------------
root = 0
do i = 1, nband
   if( i.gt.nbndsp(root+1)+nbncnt(root+1) ) root = root + 1
   do j = 1, nband
      prod(j) = 0.d0
      do m = 1, mshnod
         prod(j) = prod(j) + cgjr(m,j)*rhcr(m,i,1)
      end do
      if( myid.eq.0 ) then
          prod(j) = prod(j)*2.d0 - cgjr(1,j)*rhcr(1,i,1)
        else
          prod(j) = prod(j)*2.d0
      end if
   end do
   !---in order to stabilize the DC calculations 2012/06/08
!   if( lvand ) then
!       do j = 1, nband
!          call calpcscc( nfile, myid, nodes,  &
!& cscr, j, i, rvol, 1, nspin )
!          prod(j) = prod(j) + cscr
!       end do
!   end if
#if PCHOLESKY
   call dsum(prod,nband,prodr,root)
   if( myid.eq.root ) then
       do j = 1, nband
          aijr(i-nbnod1+1,j) = prod(j)
       end do
   end if
#else
   call gdsum(prod,nband,prodr)
       do j = 1, nband
          aijr(i,j) = prod(j)
       end do
#endif
end do
if( ltimecnt ) then
    ct = timecnt()
    tc(1) = tc(1) + ct - ct0
    ct0 = ct
end if
!-----------------------------------------------------------------------

!!---TDDFT-FSSH: set matrix of wavefunctions
!call tddft_fssh_set_matrix( nfile, myid, nodes,  &
!& aijr, nband, nbnod, nbxxxx, nspin, prod, prodr )


!--- calculate U(Hermite conjugate)*U ----------------------------------
! --- This is serial code. ---
!      do ib=1, nband
!      do jb=1, nband
!         bijr(ib,jb) = 0.d0
!         do kb=1, nband
!            bijr(ib,jb) = bijr(ib,jb) + aijr(kb,ib)*aijr(kb,jb)
!         end do
!      end do
!      end do
! ----------------------------
#if PCHOLESKY
root = 0
do i = 1, nband
   if( i.gt.nbndsp(root+1)+nbncnt(root+1) ) root = root + 1
   do j = 1, i
      prod(j) = 0.d0
      do k = 1, nbnod2 - nbnod1 + 1
         prod(j) = prod(j) + aijr(k,i)*aijr(k,j)
      end do
   end do
   call dsum(prod,i,prodr,root)
   if( myid.eq.root ) then
       do j = 1, i
          bijr(i-nbnod1+1,j) = prod(j)
       end do
   end if
end do
#else
do i = 1, nband
   do j = 1, i
      prod(j) = 0.d0
      do k = 1, nband
         prod(j) = prod(j) + aijr(k,i)*aijr(k,j)
      end do
   end do
   do j = 1, i
      bijr(i,j) = prod(j)
   end do
end do
#endif
if( ltimecnt ) then
    ct = timecnt()
    tc(2) = tc(2) + ct - ct0
    ct0 = ct
end if
!-----------------------------------------------------------------------


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
&           ' warning : error-1 in M(R)EIGQR( in sbalin )', checkc
    if(loutfile(2)) write(nfile(2),*)  &
&           ' warning : error-1 in M(R)EIGQR( in sbalin )', checkc
end if
!-----------------------------------------------------------------------
if( ltimecnt ) then
    ct = timecnt()
    tc(3) = tc(3) + ct - ct0
    ct0 = ct
end if


!--- Unitary matrix A' : ( bijr, biji ) --------------------------------
! --- This is serial code. ---
!      do ib=1, nband
!      do jb=1, nband
!         bm1r(ib,jb) =  aijr(jb,ib)
!         bx0r(ib,jb) =  bijr(ib,jb)
!      end do
!      end do
!      do ib=1, nband
!         cosi = 1.d0/sqrt(ew1(ib))
!      do jb=1, nband
!         bijr(ib,jb) = 0.d0
!         do kb=1, nband
!            bijr(ib,jb) = bijr(ib,jb) + bx0r(kb,ib)*bm1r(kb,jb)
!         end do
!         bijr(ib,jb) = bijr(ib,jb)*cosi
!      end do
!      end do
! ----------------------------
#if PCHOLESKY
do jb = 1, nband
do i = 1, nbnod2 - nbnod1 + 1
   bx0r(i,jb) =  bijr(i,jb)
end do
end do
!--- All to All communication : transpose matirix aijr
root = 0
do i = 1, nband
   if( i.gt.nbndsp(root+1)+nbncnt(root+1) ) root = root + 1
   do j = 1, nbnod2 - nbnod1 + 1
      w1i(j) = aijr(j,i)
   end do
   call dgatherv(w1i,nbnod2-nbnod1+1,prod,nbncnt,nbndsp,root)
   if( myid.eq.root ) then
       do j = 1, nband
          bm1r(i-nbnod1+1,j) = prod(j)
       end do
   end if
end do
if( ltimecnt ) then
    ct = timecnt()
    tc(4) = tc(4) + ct - ct0
    ct0 = ct
end if

root = 0
do i = 1, nband
   if( i.gt.nbndsp(root+1)+nbncnt(root+1) ) root = root + 1
   do j = 1, nband
      prod(j) = 0.d0
      do k = 1, nbnod2 - nbnod1 + 1
         prod(j) = prod(j) + bx0r(k,i)*bm1r(k,j)
      end do
   end do
   call dsum(prod,nband,prodr,root)
   if( myid.eq.root ) then
       cosi = 1.d0/sqrt(ew1(i))
       do j = 1, nband
          bijr(i-nbnod1+1,j) = prod(j)*cosi
       end do
   end if
end do
#else
do j = 2, nband
do k = 1, j - 1
   t1 = aijr(j,k)
   aijr(j,k) = aijr(k,j)
   aijr(k,j) = t1
end do
end do
do i = 1, nbnod
   ib = i + nbnod1 - 1
   do j = 1, nband
      prod(j) = 0.d0
      do k = 1, nband
         prod(j) = prod(j) + bijr(k,ib)*aijr(k,j)
      end do
   end do
   cosi = 1.d0/sqrt(ew1(ib))
   do j = 1, nband
      bx0r(i,j) = prod(j)*cosi
   end do
end do
#endif
if( ltimecnt ) then
    ct = timecnt()
    tc(2) = tc(2) + ct - ct0
    ct0 = ct
end if
!-----------------------------------------------------------------------


!---  A^{+} * A' : ( aijr, aiji ) --------------------------------------
! --- This is serial code. ---
!      do ib=1, nband
!      do jb=1, nband
!         bm1r(ib,jb) =  bx0r(jb,ib)
!      end do
!      end do
!      do ib=1, nband
!      do jb=1, nband
!         aijr(ib,jb) = 0.d0
!         do kb=1, nband
!            aijr(ib,jb) = aijr(ib,jb) + bm1r(kb,ib)*bijr(kb,jb)
!         end do
!      end do
!      end do
! ----------------------------
#if PCHOLESKY
!--- All to All communication : transpose matirix aijr
root = 0
do i = 1, nband
   if( i.gt.nbndsp(root+1)+nbncnt(root+1) ) root = root + 1
   do j = 1, nbnod2 - nbnod1 + 1
      w1i(j) = bx0r(j,i)
   end do
   call dgatherv(w1i,nbnod2-nbnod1+1,prod,nbncnt,nbndsp,root)
   if( myid.eq.root ) then
       do j = 1, nband
          bm1r(i-nbnod1+1,j) = prod(j)
       end do
   end if
end do
if( ltimecnt ) then
    ct = timecnt()
    tc(4) = tc(4) + ct - ct0
    ct0 = ct
end if

root = 0
do i = 1, nband
   if( i.gt.nbndsp(root+1)+nbncnt(root+1) ) root = root + 1
   do j = 1, nband
      prod(j) = 0.d0
      do k = 1, nbnod2 - nbnod1 + 1
         prod(j) = prod(j) + bm1r(k,i)*bijr(k,j)
      end do
   end do
   call dsum(prod,nband,prodr,root)
   if( myid.eq.root ) then
       do j = 1, nband
          aijr(i-nbnod1+1,j) = prod(j)
       end do
   end if
end do
#else
do jb = 1, nband
do i = 1, nbnod2 - nbnod1 + 1
   bm1r(i,jb) =  bijr(jb,i+nbnod1-1)
end do
end do

root = 0
do i = 1, nband
   if( i.gt.nbndsp(root+1)+nbncnt(root+1) ) root = root + 1
   do j = 1, nband
      prod(j) = 0.d0
      do k = 1, nbnod2 - nbnod1 + 1
         prod(j) = prod(j) + bm1r(k,i)*bx0r(k,j)
      end do
   end do
   call gdsum(prod,nband,prodr)
       do j = 1, nband
          aijr(i,j) = prod(j)
       end do
end do
#endif
if( ltimecnt ) then
    ct = timecnt()
    tc(2) = tc(2) + ct - ct0
    ct0 = ct
end if
!-----------------------------------------------------------------------


!c--- Unitary transformation --------------------------------------------
! --- This is serial code. ---
!      do ib = 1, nband
!         do ig = 1, mshnod
!            c2gr(ig,ib) = 0.0
!         end do
!         do jb = 1, nband
!            do ig = 1, mshnod
!               c2gr(ig,ib) = c2gr(ig,ib) + aijr(ib,jb)*rhcr(ig,jb)
!            end do
!         end do
!      end do
!      do ib = 1, nband
!         do ig = 1, mshnod
!            rhcr(ig,ib) = c2gr(ig,ib)
!         end do
!      end do
! ----------------------------
nmaxx = mshnod
call gimax(nmaxx)
ncblck = nmaxx/iblock
do keypwv = 1, min( nprvmx, kpwvx )

do icblck = 1, ncblck + 1
   icb1 = (icblck-1)*iblock + 1
   icb2 =  icblck   *iblock
   icb2 = min( icb2, mshnod )

   root = 0
   do ib = 1, nband
#if PCHOLESKY
   if( ib.gt.nbndsp(root+1)+nbncnt(root+1) ) root = root + 1
!             do jb = 1, nbnod2 - nbnod1 + 1
!                w1i(jb) = aijr(jb,ib)
!             end do
!             call alldgatherv(w1i,nbnod2-nbnod1+1,prod,nbncnt,nbndsp)
       if( myid.eq.root ) then
           do jb = 1, nband
              prod(jb) = aijr(ib-nbnod1+1,jb)
           end do
       end if
       call dbcast(prod,nband,root)
#else
       do jb = 1, nband
          prod(jb) = aijr(ib,jb)
       end do
#endif

       do i = 1, iblock
          tmpr(i,ib) = 0.d0
       end do
       do jb = 1, nband
       do i = icb1, icb2
          tmpr(i-icb1+1,ib) = tmpr(i-icb1+1,ib)  &
&                           + prod(jb)*rhcr(i,jb,keypwv)
       end do
       end do
   end do
   do ib = 1, nband
   do i = icb1, icb2
      rhcr(i,ib,keypwv) = tmpr(i-icb1+1,ib)
   end do
   end do

end do

!----- rotation of slmir
!---in order to stabilize the DC calculations 2012/06/08
!if( lvand ) then
!    call trans_prvslmi2( nfile, myid, nodes,  &
!& aijr, nband, nbnod1, nbnod2, nbnod, nbxxxx, keypwv, nspin, prod )
!end if


!!check
!      do i = 1, nband
!         do j = 1, nband
!            prod(j) = 0.d0
!            do m = 1, mshnod
!               prod(j) = prod(j) + rhcr(m,j,keypwv)*rhcr(m,i,keypwv)
!            end do
!            if( myid.eq.0 ) then
!                prod(j) = prod(j)*2.d0
!     &                  - rhcr(1,j,keypwv)*rhcr(1,i,keypwv)
!              else
!                prod(j) = prod(j)*2.d0
!            end if
!         end do
!         if( lvand ) then
!             do j = 1, nband
!                call calpcspc( nfile, myid, nodes,
!     & cscr, i, j, rvol, keypwv, keypwv, nspin )
!                prod(j) = prod(j) + cscr
!             end do
!         end if
!         call gdsum(prod,nband,prodr)
!         lerror = .false.
!         errmx  = 0.d0
!         do j = 1, nband
!         if( i.eq.j ) prod(j) = prod(j) - 1.d0
!         lerror = lerror .or. abs(prod(j)).gt.1.d-10
!         errmx  = max( errmx, abs(prod(j)) )
!!         if( abs(prod(j)).gt.1.d-10 ) then
!!         if( myid.eq.0 ) then
!!             write(nfile(1),'(a39,2i5,e15.6)')
!!     &       'orthonormalization error in roteig:',i, j, prod(j)
!!         end if
!!         end if
!         end do
!         if( lerror ) then
!         if( myid.eq.0 ) then
!             write(nfile(1),'(a39,i5,e15.6)')
!     &       'orthonormalization error (3) in roteig:',i, errmx
!         end if
!         end if
!      end do
!!check

end do


if( ltimecnt ) then
    ct = timecnt()
    tc(5) = tc(5) + ct - ct0
    ct0 = ct
end if
!-----------------------------------------------------------------------


if( ltimecnt ) then
do i = 1, 2
   if( loutfile(i) ) then
       write(nfile(i),*) '                        set U-matrix',  &
&                     ' : cpu-time :', tc(1)
       write(nfile(i),*) '                              reigqr',  &
&                     ' :          :', tc(3)
       write(nfile(i),*) '         matrix-by-matrix operations',  &
&                     ' :          :', tc(2)
       write(nfile(i),*) '               all-to-all operations',  &
&                     ' :          :', tc(4)
       write(nfile(i),*) '              Unitary transformation',  &
&                     ' :          :', tc(5)
       write(nfile(i),*) '                          sub total ',  &
&   ' :          :', tc(1)+tc(2)+tc(3)+tc(4)+tc(5)+tc(6)+tc(7)
   end if
end do
end if



return
end




subroutine chgest( nfile, myid, nodes,  &
& rhgr, rinr, nplw5ex, nplw5, prvrho, ncprvmx, nspnmx2, thrhgr,  &
& ichest, keypcd, keypmg, nstep, nstep_ini, facc1, facc2, facc3,  &
& lspin, rhoud, diffud, lfixud,  &
& glocal, rho, mshnod, rdelv, nel,  &
& mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh, laspc_exec )
!-----------------------------------------------------------------------
!     estimation of input charge
!-----------------------------------------------------------------------
use aspc_param
use aspc
implicit none
integer :: nfile(*), myid, nodes
integer :: nplw5ex, nplw5
real*8  :: rhgr(nplw5ex,*), rinr(nplw5ex,*)
integer :: ncprvmx, nspnmx2
real*8  :: prvrho(nplw5ex,ncprvmx,1:nspnmx2)
real*8  :: thrhgr(nplw5ex)
integer :: ichest, keypcd, keypmg, nstep, nstep_ini
real*8  :: facc1, facc2, facc3
logical :: lspin, lfixud
real*8  :: rhoud(*), diffud
real*8  :: glocal(*)
real*8  :: rho(*)
integer :: mshnod(*)
real*8  :: rdelv
integer :: nel
integer :: mftnod(*), mftdsp(*), mfd2ft(*), ntotfd, nd1vks(*), kfft0
real*8  :: fft3x(*), fft3y(*)
complex*16 :: fftwork(*)
integer :: ngenh, ijkgd(-ngenh:ngenh)
logical :: laspc_exec


!if( lxlbomd .and. laspc_exec ) then
!    !---Extended Lagrangian scheme XL-BOMD
!    call chgxlbomd( nfile, myid, nodes,  &
!& rhgr, nplw5ex, nplw5, prvrho, ncprvmx, nspnmx2, thrhgr,  &
!& kxlbomd, dkappa, c0 )
!
!else

!if( laspc_chg ) then

!----- To keep backward compatibility
!      The following settings give the same results as chgestold
!    coefaspc(1,2) = 1.d0 + facc3
!    coefaspc(2,2) = - facc3
!    coefaspc(1,3) = 1.d0 + facc1
!    coefaspc(2,3) = facc2  - facc1
!    coefaspc(3,3) = - facc2

    call chgestaspc( nfile, myid, nodes,  &
& rhgr, nplw5ex, nplw5, prvrho, ncprvmx, nspnmx2, thrhgr,  &
& ichest, keypcd, keypmg, nstep, nstep_ini, coefaspc, kaspc )

!else
!
!    call chgestold( nfile, myid, nodes,  &
!& rhgr, nplw5ex, nplw5, prvrho, ncprvmx, nspnmx2, thrhgr,  &
!& ichest, keypcd, keypmg, nstep, nstep_ini, facc1, facc2, facc3 )
!
!end if
!end if

call chgest2( nfile, myid, nodes,  &
& rhgr, rinr, nplw5ex, nplw5, nspnmx2, &
& lspin, rhoud, diffud, lfixud,  &
& glocal, rho, mshnod, rdelv, nel,  &
& mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh )


return
end




subroutine chgestaspc( nfile, myid, nodes,  &
& rhgr, nplw5ex, nplw5, prvrho, ncprvmx, nspnmx, thrhgr,  &
& ichest, keypcd, keypmg, nstep, nstep_ini, coefaspc, kaspc )
!-----------------------------------------------------------------------
!     estimation of input charge
!-----------------------------------------------------------------------
use xlbomd
implicit none
integer :: nfile(*), myid, nodes
integer :: nplw5ex, nplw5
real*8  :: rhgr(nplw5ex,*)
integer :: ncprvmx, nspnmx
real*8  :: prvrho(nplw5ex,ncprvmx,1:nspnmx)
real*8  :: thrhgr(nplw5ex)
integer :: ichest, keypcd, keypmg, nstep, nstep_ini
integer :: kaspc
real*8  :: coefaspc(kaspc,kaspc)

!------declare local variables
integer :: nspin, i, iaspc(4)


if( lchgxlbomd ) then
    do nspin = 1, nspnmx
    do i = 2, ncprvmx
       prvrho(1:nplw5ex,i-1,nspin) = prvrho(1:nplw5ex,i,nspin)
    end do
    end do

    keypcd = ncprvmx - 1
    if( nspnmx >= 2 ) keypmg = ncprvmx - 1
    lchgxlbomd = .false.
end if


do i = 0, ncprvmx
   if( ichest == i .or. nstep == nstep_ini+1+i .or. keypcd == i ) then
       iaspc(1) = i + 1
       exit
   end if
end do

keypcd = min( keypcd + 1, ncprvmx )

if( nspnmx >= 2 ) then
    do i = 0, ncprvmx
       if( ichest == i .or. nstep == nstep_ini+1+i .or. keypmg == i ) then
           iaspc(2:nspnmx) = i + 1
           exit
       end if
    end do

    keypmg = min( keypmg + 1, ncprvmx )
end if

!-----Extrapolation from previous charge denisties
do nspin = 1, nspnmx
   if( iaspc(nspin) > 1 ) then
       thrhgr(1:nplw5ex) = coefaspc(1,iaspc(nspin))*rhgr(1:nplw5ex,nspin)
       do i = 2, iaspc(nspin)
          thrhgr(1:nplw5ex) = thrhgr(1:nplw5ex)  &
&                           + coefaspc(i,iaspc(nspin))*prvrho(1:nplw5ex,i-1,nspin)
       end do
   end if

   do i = ncprvmx - 1, 1, -1
      prvrho(1:nplw5ex,i+1,nspin) = prvrho(1:nplw5ex,i,nspin)
   end do
   prvrho(1:nplw5ex,1,nspin) = rhgr(1:nplw5ex,nspin)

   if( iaspc(nspin) > 1 ) then
       rhgr(1:nplw5ex,nspin) = thrhgr(1:nplw5ex)
   end if
end do


return
end




subroutine chgest2( nfile, myid, nodes,  &
& rhgr, rinr, nplw5ex, nplw5, nspnmx, &
& lspin, rhoud, diffud, lfixud,  &
& glocal, rho, mshnod, rdelv, nel,  &
& mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh )
!-----------------------------------------------------------------------
!     estimation of input charge
!-----------------------------------------------------------------------
use outfile
use ncmagne_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: nplw5ex, nplw5, nspnmx
real*8  :: rhgr(nplw5ex,*), rinr(nplw5ex,*)
logical :: lspin, lfixud
real*8  :: rhoud(*), diffud
real*8  :: glocal(*)
real*8  :: rho(*)
integer :: mshnod(*)
real*8  :: rdelv
integer :: nel
integer :: mftnod(*), mftdsp(*), mfd2ft(*), ntotfd, nd1vks(*), kfft0
real*8  :: fft3x(*), fft3y(*)
complex*16 :: fftwork(*)
integer :: ngenh, ijkgd(-ngenh:ngenh)

!------declare local variables
integer :: m1
real*8  :: telcor, dbuf1r


!--- charge density transformation from g- to r-spaces
call chgg2r( nfile, myid, nodes,  &
& glocal, rhgr, nplw5ex, nplw5, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh )

!--- check No. of electrons and charge density
!      call chkchg( nfile, myid, nodes, glocal, rdelv, nel, ntotfd )
call chkchg2( nfile, myid, nodes,  &
& glocal, rdelv, nel, ntotfd, telcor )
rhgr(1:nplw5ex,1)  = rhgr(1:nplw5ex,1) * telcor

!--- store local charge density
call distlc( nfile, myid, nodes, rho, mshnod(1), glocal, mftdsp )

if( lspin ) then
    !--- spin charge density transformation from g- to r-spaces
    call chgg2r( nfile, myid, nodes,  &
& glocal, rhgr(1,2), nplw5ex, nplw5, mfd2ft, ntotfd, nd1vks, kfft0,  &
& fft3x, fft3y, fftwork, ijkgd, ngenh )

    !--- check No. of electrons and charge density
!          call chkcud( nfile, myid, nodes,
!     & glocal, rdelv, diffud, lfixud, ntotfd )
   call chkcud2( nfile, myid, nodes,  &
& glocal, rdelv, diffud, lfixud, ntotfd, telcor )
    if( lfixud .and. abs(telcor).gt.1.d-10 ) then
        telcor = dble(diffud)/telcor
        rhgr(1:nplw5ex,2)  = rhgr(1:nplw5ex,2) * telcor
    end if

    !--- store local charge density
    call distlc( nfile, myid, nodes, rhoud, mshnod(1), glocal, mftdsp )
end if


!--- noncollinear magnetism
!if( nspnmx == 4 ) then
!    !--- rhomx
!    call chgg2r( nfile, myid, nodes,  &
!& glocal, rhgr(1,2), nplw5ex, nplw5, mfd2ft, ntotfd, nd1vks, kfft0,  &
!& fft3x, fft3y, fftwork, ijkgd, ngenh )
!
!    !--- store local charge density
!    call distlc( nfile, myid, nodes, rhomx, mshnod(1), glocal, mftdsp )
!
!    !--- rhomy
!    call chgg2r( nfile, myid, nodes,  &
!& glocal, rhgr(1,3), nplw5ex, nplw5, mfd2ft, ntotfd, nd1vks, kfft0,  &
!& fft3x, fft3y, fftwork, ijkgd, ngenh )
!
!    !--- store local charge density
!    call distlc( nfile, myid, nodes, rhomy, mshnod(1), glocal, mftdsp )
!
!    !--- rhomz
!    call chgg2r( nfile, myid, nodes,  &
!& glocal, rhgr(1,4), nplw5ex, nplw5, mfd2ft, ntotfd, nd1vks, kfft0,  &
!& fft3x, fft3y, fftwork, ijkgd, ngenh )
!
!    !--- store local charge density
!    call distlc( nfile, myid, nodes, rhomz, mshnod(1), glocal, mftdsp )
!
!    !---check
!    telcor = 0.d0
!    do m1 = 1, mshnod(1)
!       telcor = telcor + rhomx( m1 )
!    end do
!    call gdsum(telcor,1,dbuf1r)
!    if(loutfile(1)) write(nfile(1),*) 'sum of rhomx', telcor * rdelv
!    if(loutfile(2)) write(nfile(2),*) 'sum of rhomx', telcor * rdelv
!
!    telcor = 0.d0
!    do m1 = 1, mshnod(1)
!       telcor = telcor + rhomy( m1 )
!    end do
!    call gdsum(telcor,1,dbuf1r)
!    if(loutfile(1)) write(nfile(1),*) 'sum of rhomy', telcor * rdelv
!    if(loutfile(2)) write(nfile(2),*) 'sum of rhomy', telcor * rdelv
!
!    telcor = 0.d0
!    do m1 = 1, mshnod(1)
!       telcor = telcor + rhomz( m1 )
!    end do
!    call gdsum(telcor,1,dbuf1r)
!    if(loutfile(1)) write(nfile(1),*) 'sum of rhomz', telcor * rdelv
!    if(loutfile(2)) write(nfile(2),*) 'sum of rhomz', telcor * rdelv
!end if


!--- store rinr as input charge density
rinr(1:nplw5ex,1:nspnmx) = rhgr(1:nplw5ex,1:nspnmx)


return
end




subroutine redpcd( nfile,  &
& fname1, lstart, ifmd, nstepCG, nstepMD,  &
& keypcd, prvrho, nplw5ex, nplw5, ncprvmx, lspin, nspnmx2, ngenh,  &
& lnoncollinear, keypmg,  &
& glocal, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d, fft3x )
!-----------------------------------------------------------------------
!     read charge densities at previous time steps
!-----------------------------------------------------------------------
use xlbomd
implicit none
integer :: nfile(*)
character(*) :: fname1
logical :: lstart, lspin
integer :: ifmd, nstepCG, nstepMD, keypcd
integer :: nplw5ex, nplw5, ncprvmx, nspnmx2, ngenh
real*8  :: prvrho(nplw5ex,ncprvmx,1:nspnmx2)
logical :: lnoncollinear
integer :: keypmg
real*8  :: glocal(*), fft3x(*)
integer :: mftnod(*), mftdsp(*), mfd2ft(*), ntotfd, nd1vks(*), kfft0d

!-----declare local variables
integer :: myid, nodes, nkd
integer :: i, j, m1, nplw5o, nmmo, keypcdo, nspnmxo, keypmgo, mg
real*8  :: dummy
integer :: iunit, nversion, istat, ierror
logical :: lnoncollinearo


if( .not.lstart .or. ifmd == 0 ) return
if( ifmd == 1 .and. nstepCG == 0 .or.  &
&   ifmd >= 2 .and. nstepMD == 0 ) return


!-----set communicator
call get_worldqm( myid, nodes )

ierror = 0
nspnmxo = 0
lnoncollinearo = .false.

if( myid.eq.0 ) then
    prvrho(1:nplw5ex,1:ncprvmx,1:nspnmx2) = 0.d0

    call allocate_unit_number( iunit )

    if( myid == 0 ) then
        write(nfile(1),*) 'open file(in redpcd): ',  &
&                         fname1(1:len_trim(fname1))
        write(nfile(2),*) 'open file(in redpcd): ',  &
&                         fname1(1:len_trim(fname1))
    end if

    open(iunit,iostat=istat,file=fname1,status='old',form='unformatted')
    if( istat == 0 ) read(iunit,iostat=istat) nplw5o
    nversion = 0
    if( istat == 0 .and. nplw5o < 0 ) then
        nversion = abs(nplw5o)
        if( istat == 0 ) read(iunit,iostat=istat) nplw5o
    end if
    if( istat == 0 ) then
        if( nplw5o.gt.ngenh ) then
            write(nfile(1),*) ' nplw5o.gt.ngenh', nplw5o, ngenh
            write(nfile(2),*) ' nplw5o.gt.ngenh', nplw5o, ngenh
            istat = 1
        end if
        if( nplw5o.ne.nplw5 ) then
            write(nfile(1),*) ' warning: nplw5o.ne.nplw5',nplw5o,nplw5
            write(nfile(2),*) ' warning: nplw5o.ne.nplw5',nplw5o,nplw5
        end if
        nmmo = 2*( nplw5o + 1 )
    end if
    if( istat == 0 ) read(iunit,iostat=istat) keypcdo
    if( istat == 0 ) then
        if( keypcdo > ncprvmx ) then
            write(nfile(1),*) ' warning: keypcdo > ncprvmx', keypcdo, ncprvmx
            write(nfile(2),*) ' warning: keypcdo > ncprvmx', keypcdo, ncprvmx
        end if
        keypcd = min( keypcdo, ncprvmx )
    end if
    if( nversion >= 3 ) then
        if( istat == 0 ) read(iunit,iostat=istat) nspnmxo
        if( istat == 0 ) read(iunit,iostat=istat) lnoncollinearo
        if( istat == 0 ) read(iunit,iostat=istat) keypmgo
        if( istat == 0 ) then
            if( keypmgo > ncprvmx ) then
                write(nfile(1),*) ' warning: keypmgo > ncprvmx', keypmgo, ncprvmx
                write(nfile(2),*) ' warning: keypmgo > ncprvmx', keypmgo, ncprvmx
            end if
            keypmg = min( keypmgo, ncprvmx )
        end if
    end if
    do i = 1, keypcd
       if( istat == 0 ) read(iunit,iostat=istat) ( prvrho(m1,i,1), m1 = 1, nmmo )
    end do
    do i = keypcd + 1, keypcdo
       if( istat == 0 ) read(iunit,iostat=istat) dummy
    end do

    if( istat == 0 ) then
        if( lspin ) then
            if( .not.lnoncollinearo ) then
                do i = 1, keypcd
                   if( istat == 0 ) read(iunit,iostat=istat) prvrho(1:nmmo,i,2)
                end do
                if( istat /= 0 ) then
                    prvrho(1:nmmo,1:ncprvmx,2) = 0.d0 !prvrho(1:nmmo,1:keypcd,1)
                    keypmg = 0
                end if
            else
                !-----from noncollinear magnetism
                prvrho(1:nmmo,1:ncprvmx,2) = 0.d0
                keypmg = 0
            end if
        else
            if( lnoncollinear ) then
                !-----noncollinear magnetism
                if( lnoncollinearo ) then
                    do mg = 2, 4
                       do i = 1, keypcd
                          if( istat == 0 ) read(iunit,iostat=istat) prvrho(1:nmmo,i,mg)
                       end do
                       do i = keypcd + 1, keypcdo
                          if( istat == 0 ) read(iunit,iostat=istat) dummy
                       end do
                    end do
                else 
                    prvrho(1:nmmo,1:keypcd,1:4) = 0.d0
                    keypmg = 0
                end if
            end if
        end if
    else
        ierror = 1
    end if
end if
!--- error trap
call gimax(ierror)


!---noncollinear magnetism
!    call ibcast(nspnmxo,1,0)
!call lbcast(lnoncollinearo,1,0)

!if( ierror == 0 ) then
!    if( lnoncollinear ) then
!        if( lnoncollinearo ) then
!            call redpmagne( nfile, myid,  &
!& fname1, iunit, prvrhom, keypmg, ncprvmx, mshnod,  &
!& glocal, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d, fft3x )
!        else
!            keypmg = 0
!        end if
!    end if
!end if


if( myid == 0 ) then
    close(iunit)
    call deallocate_unit_number( iunit )
end if

if( ierror == 0 ) then

    !--- broadcast charge density to all nodes
    call ibcast(keypcd,1,0)
    call ibcast(keypmg,1,0)
    call dbcast(prvrho,nplw5ex*keypcd*nspnmx2,0)

    if( myid == 0 ) then
        write(nfile(1),*) ' read charge densities ',  &
&                         'at previous time steps'
        write(nfile(2),*) ' read charge densities ',  &
&                         'at previous time steps'
    end if

else

    if( myid == 0 ) then
        write(nfile(1),*) 'error : in file pcd000 ...'
        write(nfile(2),*) 'error : in file pcd000 ...'
    end if
!       -----Finalize the parallel environment
!              call end_parallel(ierr)
!              stop
    keypcd = 0
    keypmg = 0

end if

!-----set communicator
call get_worldkd( myid, nodes, nkd )


return
end




!!================== cut from here ======================================
!
!subroutine redpmagne( nfile, myid,  &
!& fname1, iunit, prvrhom, keypmg, ncprvmx, mshnod,  &
!& glocal, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d, fft3x )
!!-----------------------------------------------------------------------
!!     save charge densities at previous time steps
!!-----------------------------------------------------------------------
!use xlbomd
!implicit none
!integer :: nfile(*), myid
!character(*) :: fname1
!integer :: iunit, keypmg, ncprvmx, mshnod
!real*8  :: prvrhom(mshnod,3,ncprvmx)
!real*8  :: glocal(*), fft3x(*)
!integer :: mftnod(*), mftdsp(*), mfd2ft(*), ntotfd, nd1vks(*), kfft0d
!
!!-----declare local variables
!integer :: nodes, myid_kd, nodes_kd, nkd, myid_lr, nodes_lr, myid_pw, nodes_pw
!integer :: i, j, kfft0do, keypmgo
!integer :: istat, ierro3
!
!
!!-----set communicator
!call get_worldlr( myid_lr, nodes_lr )
!
!!-----set communicator
!call get_worldpw( myid_pw, nodes_pw )
!
!!-----set communicator
!call get_worldkd( myid_kd, nodes_kd, nkd )
!
!nlrif: if( myid_lr == 0 ) then
!npwif: if( myid_pw == 0 ) then
!nkdif: if( nkd == 0 ) then
!
!    if( myid == 0 ) then
!        read(iunit,iostat=istat) kfft0do
!        if( istat == 0 ) then
!            if( kfft0do /= kfft0d .and. myid.eq.0 ) then
!                ierro3 = 1
!                write(nfile(1),*) 'error - different kfft0d in ', fname1
!                write(nfile(2),*) 'error - different kfft0d in ', fname1
!            end if
!        else
!            ierro3 = 1
!        end if
!        read(iunit,iostat=istat) keypmgo
!        if( istat == 0 ) then
!            if( keypmgo /= keypmg .and. myid.eq.0 ) then
!                write(nfile(1),*) 'warning - different keypmg in ', fname1
!                write(nfile(2),*) 'warning - different keypmg in ', fname1
!            end if
!            keypmg = min( keypmgo, ncprvmx )
!        else
!            ierro3 = 1
!        end if
!    end if
!    !--- error trap
!    call gimax(ierro3)
!    ier0: if( ierro3 == 0 ) then
!
!        call ibcast(keypmg,1,0)
!        ido: do i = 1, keypmg
!        do j = 1, 3
!           !---rhomx, y, z for j = 1, 2, 3
!           call rd_local_data( nfile, myid, nodes, iunit, ierro3,  &
!& prvrhom(1,j,i), mshnod, glocal, mftnod, mftdsp, mfd2ft, ntotfd, kfft0d, fft3x )
!           if( ierro3 /= 0 ) exit ido
!        end do
!        end do ido
!
!    end if ier0
!
!end if nkdif
!
!!-----set communicator
!call get_worldun( myid, nodes )
!
!!--- error trap
!call gimax(ierro3)
!if( ierro3 == 0 ) then
!    call ibcast(keypmg,1,0)
!    call dbcast(prvrhom,mshnod*3*keypmg,0)
!end if
!
!end if npwif
!
!
!if( nodes_pw > 1 ) then
!    !-----set communicator
!    call get_worldpw( myid, nodes )
!
!    !--- error trap
!    call gimax(ierro3)
!    if( ierro3 == 0 ) then
!        call ibcast(keypmg,1,0)
!        call dbcast(prvrhom,mshnod*3*keypmg,0)
!    end if
!
!end if
!
!end if nlrif
!
!if( nodes_lr > 1 ) then
!    !-----set communicator
!    call get_worldlr( myid, nodes )
!
!    !--- error trap
!    call gimax(ierro3)
!    if( ierro3 == 0 ) then
!        call ibcast(keypmg,1,0)
!        call dbcast(prvrhom,mshnod*3*keypmg,0)
!    end if
!
!end if
!
!!-----set communicator
!call get_worldqm( myid, nodes )
!
!
!return
!end
!
!!================== cut up to here ======================================




subroutine savpcd( nfile,  &
& fname1, ifmd, keypcd, prvrho, nplw5ex, nplw5, ncprvmx, nspnmx2,  &
& lnoncollinear, keypmg,  &
& glocal, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d, fft3x )
!-----------------------------------------------------------------------
!     save charge densities at previous time steps
!-----------------------------------------------------------------------
use xlbomd
implicit none
integer :: nfile(*)
character(*) :: fname1
integer :: ifmd, keypcd, nplw5ex, nplw5, ncprvmx, nspnmx2
real*8  :: prvrho(nplw5ex,ncprvmx,1:nspnmx2)
logical :: lnoncollinear
integer :: keypmg
real*8  :: glocal(*), fft3x(*)
integer :: mftnod(*), mftdsp(*), mfd2ft(*), ntotfd, nd1vks(*), kfft0d

!-----declare local variables
integer :: myid, nodes, myid_kd, nodes_kd, nkd
integer :: nspin, i, m1
integer :: iunit, nversion
real*8  :: ct, ct0, timecnt


if( ifmd.lt.1 ) return


!-----set communicator
call get_worldqm( myid, nodes )


ct0 = timecnt()

if( myid == 0 ) then
    call allocate_unit_number( iunit )

    if( myid == 0 ) then
        write(nfile(1),*) 'open file(in savpcd): ',  &
&                         fname1(1:len_trim(fname1))
        write(nfile(2),*) 'open file(in savpcd): ',  &
&                         fname1(1:len_trim(fname1))
    end if

    open(iunit, file=fname1, status='unknown', form='unformatted')
    nversion = -3
    write(iunit) nversion
    write(iunit) nplw5
    write(iunit) keypcd
    write(iunit) nspnmx2         ! <- nversion = -3
    write(iunit) lnoncollinear   ! <- nversion = -3
    write(iunit) keypmg          ! <- nversion = -3
    do nspin = 1, nspnmx2
    do i = 1, keypcd
       write(iunit) ( prvrho(m1,i,nspin), m1 = 1, nplw5ex )
    end do
    end do
end if


!-----set communicator
call get_worldkd( myid_kd, nodes_kd, nkd )


!!---non-collinear magnetism
!if( lnoncollinear ) then
!    call savpmagne( nfile, myid, myid_kd, nodes_kd,  &
!& iunit, prvrhom, keypmg, ncprvmx, mshnod,  &
!& glocal, mftnod, mftdsp, mfd2ft, ntotfd, kfft0d, fft3x )
!end if


if( myid == 0 ) then
    close(iunit)
    call deallocate_unit_number( iunit )
end if

ct = timecnt()
if( myid.eq.0 ) then
    write(nfile(1),*) '              savpcd :          :', ct - ct0
    write(nfile(2),*) '              savpcd :          :', ct - ct0
end if
ct0 = ct


!!-----set communicator
!call get_worldkd( myid, nodes, nkd )


return
end




!!================== cut from here ======================================
!
!subroutine savpmagne( nfile, myid, myid_kd, nodes_kd,  &
!& iunit, prvrhom, keypmg, ncprvmx, mshnod,  &
!& glocal, mftnod, mftdsp, mfd2ft, ntotfd, kfft0d, fft3x )
!!-----------------------------------------------------------------------
!!     save charge densities at previous time steps
!!-----------------------------------------------------------------------
!use xlbomd
!implicit none
!integer :: nfile(*), myid, myid_kd, nodes_kd
!integer :: iunit, keypmg, ncprvmx, mshnod
!real*8  :: prvrhom(mshnod,3,ncprvmx)
!real*8  :: glocal(*), fft3x(*)
!integer :: mftnod(*), mftdsp(*), mfd2ft(*), ntotfd, kfft0d
!
!!-----declare local variables
!integer :: i, j
!
!
!if( myid == 0 ) then
!    write(iunit) kfft0d
!    write(iunit) keypmg
!end if
!
!do i = 1, keypmg
!   do j = 1, 3
!      !---rhomx, y, z for j = 1, 2, 3
!      call sv_local_data( nfile, myid, iunit,  &
!& prvrhom(1,j,i), mshnod, glocal, mftnod, mftdsp, mfd2ft, ntotfd, kfft0d, fft3x )
!   end do
!end do
!
!
!return
!end
!
!!================== cut up to here ======================================




subroutine redpwv( nfile, iogpsz, ct0,  &
& fname1, lstart, ifmd, nstepCG, nstepMD, keypwv,  &
& cgjr, rhcr, bufcr, iodg, ioag, idstnd,  &
& nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& prveig, npnod1, npnod2, npnod, nplcnt, npldsp, nprvmx, nspnmx,  &
& prod, prodr, pdbuf, dmtrxr, nbxxxx,  &
& node_c, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, bijr, bm1r, ncscale )
!-----------------------------------------------------------------------
!     read wavefunctions at previous time steps
!-----------------------------------------------------------------------
use outfile
use xlbomd
implicit real*8 ( a-h, o-z )
dimension nfile(*)
integer :: iogpsz
character(*) :: fname1
dimension cgjr(nplwex,*)
dimension rhcr(nplwex,*)
dimension prveig(2*npnod,nband,nprvmx,1:nspnmx)
dimension nbncnt(*), nbndsp(*)
dimension nplcnt(*), npldsp(*)
dimension bufcr(*)
dimension iodg(*), ioag(*)
dimension idstnd(*)
dimension pdbuf(*)
dimension prod(*), prodr(*)
dimension dmtrxr(nbxxxx,*)
dimension lbncnt(*), lbndsp(*)
dimension mbncnt(*), mbndsp(*)
dimension jdstnd(*)
dimension bijr(*), bm1r(*)
logical   lstart, lsifreal, lrdslmir, lrdrhoij
integer :: ncscale

!------declare local variables
character(80) :: fname
integer :: iunit, istat = 0, digit
real*8  :: dummy
integer*2 :: saveigv
integer :: nplwoex2
logical :: ltimecnt
integer :: nfiles, ii, isnd, nsd, ibuf(10)


if( .not.lstart ) return
!if( .not.lstart .or. ifmd == 0 ) return
!if( ifmd == 1 .and. nstepCG == 0 .or.  &
!&   ifmd >= 2 .and. nstepMD == 0 ) return

ctt0 = ct0

!-----set communicator
call get_worldlr( myid_lr, nodes_lr )

!-----set communicator
call get_worldpw( myid_pw, nodes_pw )

!-----set communicator
call get_worldkd( myid, nodes, nkd )

!-----set communicator
call get_worldqmun( myid_qm_un, nodes_qm_un )


nlrif: if( myid_lr == 0 ) then
npwif: if( myid_pw == 0 ) then

if( myid.eq.0 ) then
    if(loutfile(1)) write(nfile(1),*) ' read wavefunctions ',  &
&                         'at previous time steps'
    if(loutfile(2)) write(nfile(2),*) ' read wavefunctions ',  &
&                         'at previous time steps'
end if

ierror = 0
lrdslmir = .false.
lrdrhoij = .false.
nversion = 0
if( myid.eq.0 ) then

  ioif: if( mod(myid_qm_un,iogpsz) == 0 ) then ! always .ture. for no divide-and-conquer MD
     nfiles=min(iogpsz,nodes_qm_un-myid_qm_un) ! always   1    for no divide-and-conquer MD

    fname = trim(fname1)

    if( myid == 0 ) then
        if(loutfile(1)) write(nfile(1),*) 'open file(in redpwv): ', trim(fname)
        if(loutfile(2)) write(nfile(2),*) 'open file(in redpwv): ', trim(fname)
    end if

    call allocate_unit_number( iunit )

    open(iunit,iostat=istat, file=fname, status='old', form='unformatted')
    if( istat == 0 ) read(iunit,iostat=istat) nbando
    if( istat == 0 ) then
        lrdslmir = nbando.lt.0
        if( lrdslmir ) then
            nversion = abs(nbando)
            read(iunit,iostat=istat) nbando
        end if
    end if
    if( istat == 0 ) read(iunit,iostat=istat) nplwo
    if( istat == 0 ) read(iunit,iostat=istat) keypwvo
    if( istat == 0 ) read(iunit,iostat=istat) lsifreal
    if( istat == 0 ) read(iunit,iostat=istat) nnspin
    if( nversion >= 2 ) then
        if( istat == 0 ) read(iunit,iostat=istat) lwavxlbomd, lchgxlbomd
    end if
    ncscaleo = 1
    if( nversion >= 3 ) then
        if( istat == 0 ) read(iunit,iostat=istat) ncscaleo
    end if
    if( nbando /= nband .or. nplwo /= nplw ) istat = 1
    if( istat == 0 ) then
        nplwoex = 2*( nplwo + 1 )
    else
       ierror = 1
    end if

    !-----for divide-and-conquer MD
    do ii=1, nfiles-1
       ibuf(:) = 0
       ibuf(1) = nbando
       ibuf(2) = nplwo
       ibuf(3) = keypwvo
       ibuf(4) = nnspin
       ibuf(5) = ncscaleo
       if( lrdslmir   ) ibuf(6) = 1
       if( lsifreal   ) ibuf(7) = 1
       if( lwavxlbomd ) ibuf(8) = 1
       if( lchgxlbomd ) ibuf(9) = 1
       ibuf(10) = ierror
       call cisend(myid_qm_un+ii,ibuf,10,myid_qm_un+ii,0)
    end do

  else ioif

    !-----for divide-and-conquer MD
    isnd=(myid_qm_un/iogpsz)*iogpsz
    call cirecvs(myid_qm_un,ibuf,10,isnd,0)
    nbando   = ibuf(1)
    nplwo    = ibuf(2)
    keypwvo  = ibuf(3)
    nnspin   = ibuf(4)
    ncscaleo = ibuf(5)
    lrdslmir   = ibuf(6) == 1
    lsifreal   = ibuf(7) == 1
    lwavxlbomd = ibuf(8) == 1
    lchgxlbomd = ibuf(9) == 1
    ierror   = ibuf(10)
    nplwoex = 2*( nplwo + 1 )

  end if ioif

end if

!-----set communicator
call get_worldkd( myid, nodes, nkd )

call gimax(ierror)
if( ierror.gt.0 ) go to 98
call ibcast(keypwvo,1,0)
call ibcast(nnspin,1,0)
call lbcast(lwavxlbomd,1,0)
call lbcast(lchgxlbomd,1,0)
call ibcast(ncscaleo,1,0)

if( ncscale /= ncscaleo ) then
    ierror = 1
    go to 98
end if

keypwv = min( keypwvo, nprvmx )

csh0 = timecnt()
ltimecnt = .true.
do nspin = 1, min( nnspin, nspnmx )
do ikey = 1, keypwv
   if( ncscale == 1 ) then
       call rdwfn22( nfile, myid, nodes, iogpsz, iunit, &
& rhcr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& cgjr, lsifreal, nplwoex, nband, nband, .true., dummy, ierror )
   else
       !---noncollinear magnetism
       if( myid.eq.0 ) nplwoex2 = nplwoex*ncscale
       call rdwfn22( nfile, myid, nodes, iogpsz, iunit, &
& rhcr, nplwex*ncscale, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& cgjr, lsifreal, nplwoex2, nband, nband, .false., dummy, ierror )
   end if
   if( ierror.gt.0 ) go to 98

   csh = timecnt()
   if( myid.eq.0 ) then
       if(loutfile(1)) write(nfile(1),*) '                read & distrib. w.f.',  &
&                         ' : cpu-time :', csh-csh0
       if(loutfile(2)) write(nfile(2),*) '                read & distrib. w.f.',  &
&                         ' : cpu-time :', csh-csh0
   end if
   csh0 = csh

!    --- convert band decomposition to G decomposition ---
   call bdtogd( nfile, myid, nodes, csh0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& prveig(1,1,ikey,nspin), bufcr,  &
& npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iodg, idstnd, .false., ltimecnt )

   ltimecnt = .false.
end do

!if( nspin == 1 .and. min(nnspin, nspnmx) == 2 .and. keypwv < keypwvo ) then
if( keypwv < keypwvo ) then
if( myid.eq.0 ) then
  ioif2: if( mod(myid_qm_un,iogpsz) == 0 ) then ! always .ture. for no divide-and-conquer MD
     nfiles=min(iogpsz,nodes_qm_un-myid_qm_un) ! always   1    for no divide-and-conquer MD

    do ii=0, nfiles-1
    do ikey = keypwv + 1, keypwvo
       if( lsifreal ) then
           do n = 1, nband
             read(iunit,iostat=ierror) dummy
             if( ierror /= 0 ) exit
           end do
       else
           do n = 1, nband
              read(iunit,iostat=ierror) efactr
              if( ierror /= 0 ) exit
              read(iunit,iostat=ierror) saveigv
              if( ierror /= 0 ) exit
           end do
       end if
       if( ierror /= 0 ) exit
    end do
    end do

  end if ioif2
end if
ierror = abs(ierror)
call gimax(ierror)
if( ierror > 0 ) then
    ierror = 0
    nnspin = 1
    exit
end if
end if

end do


!----- read previous slmir
call lbcast(lrdslmir,1,0)
if( lrdslmir ) then
!    call rdslmir( nfile, myid, nodes, iogpsz, iunit, &
!& keypwv, keypwvo, nnspin, lrdslmir,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& iodg, ioag, idstnd )
    lrdslmir = .false.
    if( lrdslmir ) then
    if( myid.eq.0 ) then
        if(loutfile(1)) write(nfile(1),*) 'read prvslmir : successful'
        if(loutfile(2)) write(nfile(2),*) 'read prvslmir : successful'
    end if
    end if
end if


!----- read previous rhoij
!call rdrhoij( nfile, myid, nodes, lrdrhoij, iunit, &
!& ioag, idstnd )
!if( lrdrhoij ) then
!if( myid.eq.0 ) then
!    if(loutfile(1)) write(nfile(1),*) 'read prvrhoij : successful'
!    if(loutfile(2)) write(nfile(2),*) 'read prvrhoij : successful'
!end if
!end if


98 continue
if( myid.eq.0 ) then
  ioif3: if( mod(myid_qm_un,iogpsz) == 0 ) then ! always .ture. for no divide-and-conquer MD
    close(iunit)
    call deallocate_unit_number( iunit )
  end if ioif3
end if
if( ierror == 0 ) then

!      ct = timecnt()
!      ct0 = ct
!      do i = 1, keypwv
!         if( myid.eq.0 ) then
!           write(nfile(1),*) ' Gram-Schmidt for w.f. at previous step',i
!           write(nfile(2),*) ' Gram-Schmidt for w.f. at previous step',i
!         end if
!         do nspin = 1, nnspin
!!--- Gram-Schmidt orthonormalizaion ------------------------------------
!             call schmidt( nfile, myid, nodes, ct0,
!     & prveig(1,1,i,nspin), npnod1, npnod2, npnod, nbnod1, nbnod2,
!     & nband,
!     & nbncnt, nbndsp, prod, prodr, pdbuf, .false., .false.,
!     & dmtrxr, nbxxxx,
!     & node_c, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, bijr, bm1r )
!!-----------------------------------------------------------------------
!         end do
!         ct = timecnt()
!         if( myid.eq.0 ) then
!          write(nfile(1),*) '              total : cpu-time :', ct - ct0
!          write(nfile(2),*) '              total : cpu-time :', ct - ct0
!         end if
!         ct0 = ct
!      end do

    if( nspnmx.eq.2 .and. nnspin.eq.1 ) then
        do ikey = 1, keypwv
           call cpgdrh( nfile, myid, nodes,  &
&          prveig(1,1,ikey,2), prveig(1,1,ikey,1), npnod, nband )
        end do
!        if( lrdslmir ) then
!            call cpprvslmir( nfile, myid, nodes, keypwv )
!        end if
!        call cpprvrhoij( nfile, myid, nodes )
    end if
!    call cpprvrhomm( nfile, myid, nodes )

else

    keypwv = 0
    if( myid.eq.0 ) then
        write(nfile(1),*) 'error : in files peg000 ...'
        write(nfile(2),*) 'error : in files peg000 ...'
    end if
!   -----Finalize the parallel environment
!          call end_parallel(ierr)
!          stop
end if
end if npwif


if( nodes_pw > 1 ) then
    !-----set communicator
    call get_worldpw( myid, nodes )
    !-----internode synchronization
    call gsync

    call ibcast(ierror,1,0)
    call ibcast(keypwv,1,0)
    call lbcast(lwavxlbomd,1,0)
    call lbcast(lchgxlbomd,1,0)
    if( ierror == 0 ) then
        do nspin = 1, nspnmx
           call dbcast(prveig(1,1,1,nspin),2*npnod*nband*keypwv,0)
        end do
        call lbcast(lrdslmir,1,0)
!        if( lrdslmir ) then
!            call bcastslmir( nfile, myid, nodes, keypwv )
!        end if
!        call bcastrhoij( nfile, myid, nodes )
    end if

end if

end if nlrif

if( nodes_lr > 1 ) then
    !-----set communicator
    call get_worldlr( myid, nodes )
    !-----internode synchronization
    call gsync

    call ibcast(ierror,1,0)
    call ibcast(keypwv,1,0)
    call lbcast(lwavxlbomd,1,0)
    call lbcast(lchgxlbomd,1,0)
    if( ierror == 0 ) then
        do nspin = 1, nspnmx
           call dbcast(prveig(1,1,1,nspin),2*npnod*nband*keypwv,0)
        end do
        call lbcast(lrdslmir,1,0)
!        if( lrdslmir ) then
!            call bcastslmir( nfile, myid, nodes, keypwv )
!        end if
!        call bcastrhoij( nfile, myid, nodes )
    end if

end if


!-----set communicator
call get_worldkd( myid, nodes, nkd )

ct = timecnt()
if(loutfile(1)) write(nfile(1),*) ' set w.f. at prev. t: cpu-time :', ct-ctt0
if(loutfile(2)) write(nfile(2),*) ' set w.f. at prev. t: cpu-time :', ct-ctt0
ct0 = ct


return
end




subroutine savpwv( nfile, iogpsz,  &
& fname1, ifmd, lsreal8, keypwv,  &
& cgjr, nplwex, nplw, nbnod1, nbnod2,  &
& nbnod, nbncnt, nbndsp, nband, prveig,  &
& rhcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, ioag, idstnd, bufcr, nprvmx, nspnmx, ncscale )
!-----------------------------------------------------------------------
!     save wavefunctions at previous time steps
!-----------------------------------------------------------------------
use outfile
use xlbomd
implicit real*8 ( a-h, o-z )
dimension nfile(*)
integer :: iogpsz
character(*) :: fname1
dimension cgjr(*)
dimension prveig(2*npnod,nband,nprvmx,1:nspnmx)
dimension rhcr(*)
dimension iod(*), iodg(*), ioag(*), nbndsp(*), nbncnt(*)
dimension nplcnt(*), npldsp(*)
dimension idstnd(*)
dimension bufcr(*)
logical   lsreal8
integer :: ncscale

!------declare local variables
character(80) :: fname
integer :: iunit, istat = 0, digit
logical :: ltimecnt
integer :: nfiles, ii, isnd, nsd, ibuf(10)


!if( ifmd.lt.1 ) return

ct00 = timecnt()

!-----set communicator
call get_worldlr( myid_lr, nodes_lr )

!-----set communicator
call get_worldpw( myid_pw, nodes_pw )

!-----set communicator
call get_worldkd( myid, nodes, nkd )

!-----set communicator
call get_worldqmun( myid_qm_un, nodes_qm_un )


nlrif: if( myid_lr == 0 ) then
npwif: if( myid_pw == 0 ) then
if( myid.eq.0 ) then

  ioif: if( mod(myid_qm_un,iogpsz) == 0 ) then ! always .ture. for no divide-and-conquer MD
!     nfiles=min(iogpsz,nodes_qm_un-myid_qm_un) ! always   1    for no divide-and-conquer MD

    fname = trim(fname1)

    call allocate_unit_number( iunit )

    if( myid == 0 ) then
        if(loutfile(1)) write(nfile(1),*) 'open file(in savpwv): ', trim(fname)
        if(loutfile(2)) write(nfile(2),*) 'open file(in savpwv): ', trim(fname)
    end if

    open(iunit, file=fname, status='unknown', form='unformatted')
    nversion = -3
    write(iunit) nversion    ! for backward compatibility
    write(iunit) nband
    write(iunit) nplw
    write(iunit) keypwv
    write(iunit) lsreal8
    write(iunit) nspnmx
    write(iunit) lwavxlbomd, lchgxlbomd   ! <- nversion = -2
    write(iunit) ncscale                  ! <- nversion = -3

  end if ioif

end if

!-----set communicator
call get_worldkd( myid, nodes, nkd )

ltimecnt = .true.
do nspin = 1, nspnmx
do i = 1, keypwv

   ct0 = timecnt()
!   --- copy prveig to rhcr
   call cpgdrh( nfile, myid, nodes, rhcr, prveig(1,1,i,nspin),  &
& npnod, nband )
!   --- to convert G decomposition to band decomposition ---
   call gdtobd( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& rhcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, idstnd, ltimecnt )

   !----- save previous wavefunctions
   call svwfn22( nfile, myid, nodes, iogpsz, iunit,  &
& rhcr, nplwex*ncscale, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, cgjr )

   ltimecnt = .false.
end do
end do

!----- save previous slmir
!call svslmir( nfile, myid, nodes, iogpsz, iunit, keypwv,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& iod, iodg, ioag, idstnd )

!----- save previous rhoij
!call svrhoij( nfile, myid, nodes, iunit,  &
!& ioag, idstnd )


if( myid.eq.0 ) then
  ioif2: if( mod(myid_qm_un,iogpsz) == 0 ) then ! always .ture. for no divide-and-conquer MD
    close(iunit)
    call deallocate_unit_number( iunit )
  end if ioif2
end if

end if npwif

if( nodes_pw > 1 ) then
    !-----set communicator
    call get_worldpw( myid, nodes )
    !-----internode synchronization
    call gsync
end if

end if nlrif

if( nodes_lr > 1 ) then
    !-----set communicator
    call get_worldlr( myid, nodes )
    !-----internode synchronization
    call gsync
end if


!-----set communicator
call get_worldkd( myid, nodes, nkd )

ct = timecnt()
if(loutfile(1)) write(nfile(1),*) '              savpwv :          :', ct- ct00
if(loutfile(2)) write(nfile(2),*) '              savpwv :          :', ct- ct00
ct0 = ct


return
end
