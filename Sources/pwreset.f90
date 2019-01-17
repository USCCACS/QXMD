



module pwreset_variables
!-----------------------------------------------------------------------
! type declaration of shared variables in pwreset.f
!-----------------------------------------------------------------------
implicit none

integer, dimension(3) :: nd1vks_old
integer :: nplw5o, nplwo, nplw2o, nplw3o, nplwco, nplw7o
integer :: npnod1o, npnod2o, npnodo, nspnodo
integer, allocatable, dimension(:) :: nplcnto, npldspo
integer :: nplw5exo, nplwexo, nplw2exo, nplw3exo, nplwcexo,  &
&          nplw7exo
integer :: kmax1o, kmax2o, kmax3o, kmax1do, kmax2do, kmax3do

integer, allocatable, dimension(:) :: ngao, ngbo, ngco
real*8  :: rvolo

!-----for noncollinear magnetism
integer :: mshnodo, ntotfdo
integer, allocatable, dimension(:) :: mftnodo, mftdspo, mfd2fto, ijkgdo

save


end module




subroutine pwlda_reset( nfile, myid, nodes, ltimecnt, ierror )
!-----------------------------------------------------------------------
!     reset LDA or GGA calculation with plane-wave method
!-----------------------------------------------------------------------
use outfile
use constants
use param
use param_atom
use pwlda_pp
use pwlda_atom
use pwlda_variables
use pwlda_proc
use pwlda_pw
use pwlda_grid
use ncmagne_variables
implicit none

integer :: myid, nodes
integer, dimension(*) :: nfile
logical :: ltimecnt
integer :: ierror

!-----declare local variables
logical :: lrecal
real*8  :: ct, ct0
!integer :: kmax1b, kmax2b, kmax3b
integer :: i
real*8  :: timecnt


!-----check if recalculation is necessary or not
call chkreset( nfile, myid, nodes, lrecal, h_bkup, hcell )
if( .not.lrecal ) return


ct0 = timecnt()

!-----not supported for cluster calculations and
!----- bulk calculations with vacuum
if( lclust .or. ldouble_grid_recip ) then
    do i = 1, 2
    if( loutfile(i) ) then
        write(nfile(i),*) ' *** not supported for cluster calculations'
        write(nfile(i),*) ' ***  and bulk calculations with vacuum'
        write(nfile(i),*) ' ***  in the variable-shape methods'
    end if
    end do
    ierror = 1001
end if

if(loutfile(1)) write(nfile(1),*) ' *** reset the plane-wave method'
if(loutfile(2)) write(nfile(2),*) ' *** reset the plane-wave method'


!----- re-read and re-check input data, if necessary
call reinput( nfile, myid, nodes, lclust )


!----- store previous variables
call reclat_save( nfile, myid, nodes,  &
& alloc_mem, nplw5, nplw, nplw2, nplw3, nplwc, nplw7,  &
& npnod1, npnod2, npnod, nplcnt, npldsp, nspnod,  &
& nplw5ex, nplwex, nplw2ex, nplw3ex, nplwcex, nplw7ex,  &
& nga, ngb, ngc, kmax1, kmax2, kmax3, kmax1d, kmax2d, kmax3d, rvol,  &
& pwscale, mshnod(1), mftnod, mftdsp, mfd2ft, ntotfd, ijkgd, nd1vks,  &
& lnoncollinear )


!--- set reciprocal vectors for plane waves
call reclat_ini( nfile, myid, nodes,  &
& nplw5, nplw, nplw2, nplw3, nplwc, nplw7, nplwcs,  &
& nplw5ex, nplwex, nplw2ex, nplw3ex, nplwcex, nplw7ex, rvol,  &
& nd1v, nd1vks, hcell, rba, ecut, ecutdens, ecutsoft, ecutc, ecutorth,  &
& lvand, rctmax, kfft1b, kfft2b, kfft3b, kfft0b,  &
& pwscale, lvshape )

volume = rvol

call reclat( nfile, myid, nodes,  &
& nd1vks, rba, ecut, ecutdens, ecutsoft, ecutorth,  &
& nplw5, nplw, nplw2, nplw3, nplwc, nplw7, nplwcs,  &
& nplw5ex, nplwex, nplw2ex, nplw3ex, nplwcex, nplw7ex, kfft0d,  &
& nga, ngb, ngc, recnrm, ijkgd, nplw5, nplw,  &
& kfft1,  kfft2,  kfft3,  kfft0, ijkg,  &
& kfft1b, kfft2b, kfft3b, kfft0b,  &
& gboxx, gboxy, gboxz, ngboxa, ngboxb, ngboxc, abgvc3, ijkgb,  &
& kmax1, kmax2, kmax3, kmax1d, kmax2d, kmax3d,  &
& kmax1b, kmax2b, kmax3b, kmax1cs, kmax2cs, kmax3cs )


!-----set indexes for planewave decomposition
call planewave_decomp( nfile, &
& nplw5ex, nplw5 )


!--- set local and non-local pseudopotential elements 
call setpp( nfile, myid, nodes,  &
& jgga, lrela, ntype, zatom, zv, lmax, lclno, lchk, ecutdens, ecutsoft,  &
& lpaw, lpaw_sym, lkbpp, lvand, lvfull, llocl, lkbppi, lvandi,  &
& lvflag, llocli,  &
& lnlpp_r, lnlpp_g, lkbpp_r, lkbpp_g, lvand_r, lvand_g,  &
& llclpp_r, llclpp_g, llking,  &
& tablc, tablca, dltlc, rmxlc, tbflc, tbflca, dltflc, rmxflc,  &
& lking, rking, gkgmax, gkgexct, rctflc,  &
& lpcc_r, lpcc_g, lpcci, rpcc, lpking, rpking, gpkgmax, gpkgexct, frac_rcomp_it,  &
& ms, aname, mxl, mx1, mx1loc, mxref, vlocli ,xitgrd, rr, nvlcl,  &
& nd1vks, rvol, nplw5, nplw, gx, gy, gz, recnrm,  &
& nplw3, gboxx, gboxy, gboxz, kfft1b, kfft2b, kfft3b, kfft0b,  &
& alloc_mem, lstress, ldouble_grid_recip, dgalpha, recnrmex,  &
& pwscale, lvshape, dkgx, dkgy, dkgz, dkrecnrm )

!-----set indexes for band decomposition
call band_decomp( nfile, myid, nodes,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, nbxxxx,  &
& iod, iodg, iods, iodsg )


!-----set indexes for G decomposition
call index_decomp( nfile, myid, nodes,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& nplwex, nplw, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& ncgcnt, ncgdsp, nspnod, lnoncollinear )

!-----set indexes for G decomposition related to nplw7
call index_decomp7( nfile, myid, nodes,  &
& nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& nplw7ex, nplw7, npnod71, npnod72, npnod7, nplcnt7, npldsp7,  &
& ncgcnt7, ncgdsp7, nspnod7, lnoncollinear )

!-----set indexes for G decomposition for k-point sampling
!call index_decomp_k( nfile, myid, nodes,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp )


!-----allocate memory for the plane wave method
call pwlda_alloc( nfile, myid, nodes )


!-----rearrange wavefunctions etc.
call rarrwv( nfile, myid, nodes, ltimecnt )


!-----allocate memory for the symmetry operations
!call symm_alloc( nfile, myid, nodes,  &
!& alloc_mem, nplw5, nga, ngb, ngc, nfftk, pwscale )

!-----allocate memory for the symmetry operations
!call kpointsym_alloc( nfile, myid, nodes,  &
!& alloc_mem, natom, nplwex, nplw, pwscale, nbnod, lnlpp_g, lnlpp_r, ncscale )


!-----allocate memory for uniform electric field for insulator
!if( lefield .and. .not.lsawtooth ) then
!    call efield_alloc( nfile, myid, nodes,  &
!& alloc_mem, nproc, nfftk, nplwex, npnod, nband, nbnod,  &
!& lvand, ntype, natom, natnod, pwscale, ncscale, nspnmx )
!    !-----prepare variables
!    call preset_efield( nfile, myid, nodes,  &
!& efield, lvand, ntype, lmax, lchk, lvandi, lvflag, mxl,  &
!& vlocli, xitgrd, rr, nvlcl )
!end if


!-----allocate memory for maximully localized wannier functions
!if( lwannier ) then
!    call wannier_alloc( nfile, myid, nodes,  &
!& alloc_mem, lwannier, nproc, nplwex, npnod, nband, nbnod,  &
!& lvand, ntype, natom, natnod, pwscale, ncscale )
!    !-----prepare variables
!    call preset_wannier( nfile, myid, nodes,  &
!& lwannier, lvand, ntype, lmax, lchk, lvandi, lvflag, mxl,  &
!& vlocli, xitgrd, rr, nvlcl )
!end if

!-----allocate memory for Mulliken population analysis
!if( lmulken ) then
!    call mulliken_alloc( nfile, myid, nodes,  &
!& alloc_mem, lmulken, ldecmpovp, lspdatom, lmdecmpdos, lpmchgat, leda,  &
!& nband, nbnod, nplw, npnod, rdecmp, lmulao, nmulao, mulaox, mxl,  &
!& ntype, nhk1, nhk2, lvand, pwscale )
!end if

!-----allocate memory for conductivity calculation
!if( lconduct ) then
!    call conductivity_alloc( nfile, myid, nodes,  &
!& alloc_mem, lconduct, ldcconduct, nband, kfft0d, nspnmx, pwscale )
!end if


!--- set atomic basis for Mulliken population analysis
!if( lmulken ) then
!
!    if( .not.lnoncollinear ) then
!        call setmulvec( nfile, myid, nodes,  &
!& nga, ngb, ngc, nplw, nplw5 )
!
!        ierror = 0
!        call setabsis( nfile, myid, nodes,  &
!& jgga, ntype, zatom, zv, lmax, lchk, lkbppi, lvandi, lpaw,  &
!& rdecmp, rxmulk, lmulao, lmulbsis, nmulao, mulaox, recnrm(nplw), gx, gy, gz,  &
!& nplw, npnod1, npnod2, npnod, ms, aname, mxl,  &
!& vlocli ,xitgrd, rr, mx1, ierror )
!    else
!        call setmulvec( nfile, myid, nodes,  &
!& ndkga, ndkgb, ndkgc, nplw, nplw )
!
!        ierror = 0
!        call setabsis( nfile, myid, nodes,  &
!& jgga, ntype, zatom, zv, lmax, lchk, lkbppi, lvandi, lpaw,  &
!& rdecmp, rxmulk, lmulao, lmulbsis, nmulao, mulaox, dkrecnrm(nplw),  &
!& dkgx, dkgy, dkgz, nplw, npnod1, npnod2, npnod, ms, aname, mxl,  &
!& vlocli ,xitgrd, rr, mx1, ierror )
!    end if
!
!    if( ierror /= 0 ) then
!        if(loutfile(1)) write(nfile(1),*) 'stop Mulliken analysis, ierror=',ierror
!        if(loutfile(2)) write(nfile(2),*) 'stop Mulliken analysis, ierror=',ierror
!    end if
!    lmulken = ierror == 0
!    leda = leda .and. lmulken
!    ierror = 0
!
!end if


!if( lvdw ) then
!    !----- set van der Waals kernel table : tbkernel
!    call setvdWkernel( nfile )

!    !----- set van der Waals kernel table : tbkernel
!    call setvdWfactorization( nfile,  &
!& recnrm, nplw5, gx, gy, gz, nd1vks, lstress )
!end if


!-----------------------------------------------------------------------
!-----get multigrid level

    call getmultg(  lclust, multg, nd1vks, nd2v, mulnd2 )
    call getsmultg( lclust, multg, nd1vks, muls, nd1vs,  &
&                   mulnd2, kulnd2, msrhmx, msrhmy, msrhmz )

    !-----allocate memory for multigrid
    call pwlda_mgrid_alloc( nfile, myid, nodes, npx, npy, npz,  &
& alloc_mem, nd1vks, nd2v, multg, multg, lspin, pwscale, lrtddft,  &
& ltddft_fssh, lnoncollinear, ncprvmx, lwell, lefield, lsawtooth )


!-----coefficients for numerical differentiation : cdv2, cdv1
!call coefnd( multg, nd2v, rdel, amcdv2, cdv2, cdv1,  &
!&            multg, ndv2x3 )

!do i = 1, 2
!if( loutfile(i) ) then
!    if( lclust ) then
!        write(nfile(i),*) ' '
!        write(nfile(i),*) ' << multigrid to solve Poisson eq. >>'
!    else
!        write(nfile(i),*) ' '
!        write(nfile(i),*) ' << grid points in r-space >>'
!    end if
!end if
!end do

!-----set mesh points
call setmsh( nfile, myid, nodes, myid, nodes,  &
& myx, myy, myz, nn, myparity, npx, npy, npz, nproc, lsphere,  &
& lclust, multg, 1, nd1vks, rdel, rdelv, rmax, cofmas, lvacuum,  &
& mshnod, mshnew, mshxyz, mulpit, mulx1, mulx2, muly1, muly2,  &
& mulz1, mulz2, mshnx, mshny, mshnz, mshx1, mshy1, mshz1,  &
& mshx, mshy, mshz, mulpms, mulpem, ismx, ismy, ismz,  &
& ismshx, ismshy, ismshz, mulpsx, mulpsy, mulpsz,  &
& mulyzs, mulzxs, mulxys, irmx, irmy, irmz, irmshx, irmshy, irmshz,  &
& mulprx, mulpry, mulprz, mulyzr, mulzxr, mulxyr,  &
& ismfpx, ismfpy, ismfpz, irmfpx, irmfpy, irmfpz, ncomct, mulnd2,  &
& multg, mulexf, muldatx, mulnpx, mulnpy, mulnpz,  &
& mulyzx, mulzxx, mulxyx, mulext,  &
& nmm, nmmcnt, nmmdsp, nmmrc1, nmmrc2,  &
& meshx1, meshx, meshy1, meshy, meshz1, meshz,  &
& jsmshx, jsmshy, jsmshz, nodyzh, nodzxh, nodxyh,  &
& idbuf, idbufr, mshndv, nhit1, mulgmx, mulgmy, mulgmz )


!-----pre-set FFT mesh on coarse grid
call cfftmsh( nfile, myid, nodes,  &
& rdel_c, rdelg_c, rdelv_c, ntotfd_c, nd1vks_c, mshglb,  &
& lclust, lsphere, rmax, cofmas, rba, rvol,  &
& kfft1, kfft2, kfft3, kfft0 )

!-----allocate memory for FFT mesh on coarse grid
call cfftmsh_alloc( nfile, myid, nodes,  &
& alloc_mem, lnlpp_r, kfft1, kfft2, kfft3, kfft0, pwscale,  &
& jhybrid, mxkpdupnm, nbnod, lefield_islts )


!-----set FFT mesh on dense grid
call fftmsh( nfile, myid, nodes,  &
& mftnod, mftdsp, mfd2ft, ntotfd, mftwrk, noddatx,  &
& mshglb, nd1vks(1), nd1vks(2), nd1vks(3),  &
& mshnx(mulpms(1)), mshny(mulpms(1)), mshnz(mulpms(1)), mshnod(1),  &
& mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1) )


!-----noncollinear magnetism
!call ncrarrmagne( nfile, myid, nodes, ltimecnt )


    !-----initialize Ewald method & shape denpendent terms
    call ewaldini( nfile, myid, nodes,  &
& ntype, zv, ecutdens, alloc_mem, mx1loc, llclpp_r, llclpp_g,  &
& tberf, tberfa, drcut, tbeff, tbeffa, drcutf,  &
& tablc, tablca, dltlc, rmxlc, tbflc, tbflca, dltflc, rmxflc,  &
& rctflc, llking, rlking, glkgmax, glkgexct, hcell, h_MD, gamma,  &
& rccc2, lhfull, vlocli ,xitgrd, rr, nvlcl, nion_nod,  &
& mshx(1), mshy(1), mshz(1),  &
& ldouble_grid_recip, nd1v, nd1vks, dgalpha, pwscale, lvshape )




!-----memory allocation for variables for the calculations
!-----in reciprocal space
call struc_atoms_alloc( nfile, myid, nodes,  &
& alloc_mem, natom, nplwcs, kmax1cs, kmax2cs, kmax3cs,  &
& kmax1cs, kmax2cs, kmax3cs, kmax1d, kmax2d, kmax3d,  &
& lnlpp_g, llclpp_g, lpcc_g, lsphexp, pwscale )


!-----for stress calcu by local pp.
call svvext_alloc( nfile, myid, nodes,  &
& alloc_mem, llclpp_r, llclpp_g, lstress, mshnod(1), pwscale )


if( lnlpp_g .or. lvand ) then

    !-----memory allocation for ultrasoft pp
!    if( lvand )  &
!&       call uspp_variables_alloc( nfile, myid, nodes,  &
!& alloc_mem, ntype, natom, natnod1, natnod2, natnod, natcnt, natdsp, ioag,  &
!& nband, nbnod, nbncnt, nbndsp,  &
!& kfft0b, nplw3, kmax1b, kmax2b, kmax3b, gboxx, gboxy, gboxz,  &
!& lstress, lspin, pwscale, jhybrid, mxkpdupnm, lefield, lsawtooth,  &
!& natom_alloc, natnod_alloc )

end if


if( lnlpp_g ) then

    !-----memory allocation for nonlocal pp in reciprocal space
    call nlg_variables_alloc( nfile, myid, nodes,  &
& alloc_mem, nplw, nplwex, nband, nbnod, nbncnt, nbndsp,  &
& ntype, lvandi, lstress, pwscale )

    !-----check if big memory allocation is allowed
    call nl_bigmem_alloc( nfile, myid, nodes,  &
& alloc_mem, natom, nplw, pwscale, natom_alloc )

end if


!-----check if big memory allocation is allowed for USPP
!if( lvand ) call uspp_bigmem_alloc( nfile, myid, nodes,  &
!& alloc_mem, natom, nplw3, pwscale, jhybrid, natnod_alloc )


ct = timecnt()
do i = 1, 2
if( loutfile(i) ) then
    write(nfile(i),*) ' pwlda_reset        : cpu-time :', ct - ct0
    write(nfile(i),*) ' '
end if
end do
ct0 = ct


return
end subroutine




subroutine chkreset( nfile, myid, nodes, lrecal, h_bkup, hcell )
!-----------------------------------------------------------------------
!     check if recalculation is necessary or not
!-----------------------------------------------------------------------
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
logical :: lrecal
real*8, dimension(3,3) :: hcell
real*8, dimension(3,3) :: h_bkup

!-----declare local variables
integer :: i, ix
real*8, parameter :: small = 1.d-13
real*8  :: dhmax


dhmax = 0.d0
do i = 1, 3
do ix = 1, 3
   dhmax = max( dhmax, abs(hcell(ix,i)-h_bkup(ix,i)) )
end do
end do
call gdmax(dhmax)
lrecal = dhmax > small


return
end subroutine




subroutine reinput( nfile, myid, nodes, lclust )
!-----------------------------------------------------------------------
!    re-read and re-check input data, if necessary
!-----------------------------------------------------------------------
implicit none

integer :: myid, nodes
integer, dimension(*) :: nfile
logical :: lclust


!----- read input variables from an input file
!      call reread_data( nfile, myid, nodes, lclust )


!----- read input variables for pseudopotential
!      call reread_pp( nfile, myid, nodes )


!----- check input variables
call recheck_data( nfile, myid, nodes, lclust )


!----- write input variables
!      call rewrite_data( nfile, myid, nodes, lclust )


return
end subroutine




subroutine recheck_data( nfile, myid, nodes, lclust )
!-----------------------------------------------------------------------
!    re-check input data
!-----------------------------------------------------------------------
use param
use param_atom
!      use constants
use pwreset_variables
implicit none

integer :: myid, nodes
integer, dimension(*) :: nfile
logical :: lclust

!-----declare local variables
integer :: i, j, ix


do i = 1, 3
   if( .not.lvacuum(i) ) vacuum(i) = hcell(i,i)
end do

lorthrhmbc = .false.

!----- store previous FFT meshes
do i = 1, 3
   nd1vks_old(i) = nd1vks(i)
end do

!----- get FD & FFT meshes
call get_mesh( nfile, myid, nodes,  &
& ecutdens, hcell, nd1v, nd1vks )

do ix = 1, 3
   cofmas(ix) = 0.5d0*( hcell(ix,1)+hcell(ix,2)+hcell(ix,3) )
end do


!----- the grid spacing : rdel
do i = 1, 3
   rdel(i)  = hcell(i,i)/dble(nd1v(i))
end do

do i = 1, 3
do j = 1, 3
   rdelg(j,i) = hcell(j,i)/dble(nd1v(i))
end do
end do


!----- the volume element for dense grid : rdelv
CALL RCIPRL( hcell, hci, volume )
rdelv = volume/dble( nd1v(1)*nd1v(2)*nd1v(3) )


return
end subroutine




subroutine reclat_save( nfile, myid, nodes,  &
& alloc_mem, nplw5, nplw, nplw2, nplw3, nplwc, nplw7,  &
& npnod1, npnod2, npnod, nplcnt, npldsp, nspnod,  &
& nplw5ex, nplwex, nplw2ex, nplw3ex, nplwcex, nplw7ex,  &
& nga, ngb, ngc, kmax1, kmax2, kmax3, kmax1d, kmax2d, kmax3d, rvol,  &
& pwscale, mshnod, mftnod, mftdsp, mfd2ft, ntotfd, ijkgd, nd1vks,  &
& lnoncollinear )
!-----------------------------------------------------------------------
!     store previous variables
!-----------------------------------------------------------------------
use pwreset_variables
implicit none

integer :: myid, nodes
integer, dimension(*) :: nfile
real*8  :: alloc_mem
integer :: nplw5, nplw, nplw2, nplw3, nplwc, nplw7
integer :: npnod1, npnod2, npnod, nspnod
integer, dimension(nodes) :: nplcnt, npldsp
integer :: nplw5ex, nplwex, nplw2ex, nplw3ex, nplwcex, nplw7ex
integer, dimension(0:nplw5) :: nga, ngb, ngc
integer :: kmax1, kmax2, kmax3, kmax1d, kmax2d, kmax3d
real*8  :: rvol
real*8  :: pwscale
integer :: mshnod, mftnod(nodes), mftdsp(nodes), mfd2ft(*), ntotfd
integer :: ijkgd(-nplw5:nplw5), nd1vks(3)
logical :: lnoncollinear

!-----declare local variables
integer :: nfftk
integer :: status
real*8  :: the_mem
integer :: ngenh = 0
integer :: lplnfftk = 0
logical :: licall = .true.
save ngenh, licall, lplnfftk


if( licall ) then

    !------allocate memory
    allocate( nplcnto(nodes), npldspo(nodes),  &
& stat=status )

    !-----noncollinear magnetism
    if( lnoncollinear ) then
        !------allocate memory
        allocate( mftnodo(nodes), mftdspo(nodes),  &
& stat=status )
    end if

end if


if( nplw5 > ngenh ) then

    !-----if already allocated, deallocate arrays
    if( allocated(ngao) ) then

        the_mem = 4.d0 * ( size(ngao) + size(ngbo) + size(ngco) )

        !------deallocate memory
        deallocate( ngao, ngbo, ngco, stat=status )

        if( allocated(ijkgdo) ) then
            the_mem = the_mem + 4.d0 * ( size(ijkgdo) )
            !------deallocate memory
            deallocate( ijkgdo, stat=status )
        end if

        !------error trap
        call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'reclat_save', .true. )

    end if

    ngenh = nplw5 * pwscale
    !------allocate memory
    allocate( ngao(0:ngenh), ngbo(0:ngenh), ngco(0:ngenh),  &
& stat=status )

    the_mem = 4.d0 * ( size(ngao) + size(ngbo) + size(ngco) )

    !-----noncollinear magnetism
    if( lnoncollinear ) then
        !------allocate memory
        allocate( ijkgdo(2*ngenh+1),  &
& stat=status )
        the_mem = the_mem + 4.d0 * ( size(ijkgdo) )
    end if

    !------error trap
    call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'reclat_save', .true. )

end if

!-----noncollinear magnetism
if( lnoncollinear ) then
    nfftk = nd1vks(1)*nd1vks(2)*nd1vks(3)
    if( nfftk   > lplnfftk ) then
        !-----if already allocated, deallocate arrays
        if( allocated(mfd2fto) ) then

            the_mem = 4.d0 * ( size(mfd2fto) )
            !------deallocate memory
            deallocate( mfd2fto, stat=status )

            !------error trap
            call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'reclat_save(2)', .true. )
        end if

        lplnfftk   = nfftk  * pwscale
        !------allocate memory
        allocate( mfd2fto(lplnfftk), stat=status )
        the_mem = 4.d0 * ( size(mfd2fto) )

        !------error trap
        call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'reclat_save(2)', .true. )

    end if
end if


nplw5o = nplw5
nplwo  = nplw
nplw2o = nplw2
nplw3o = nplw3
nplwco = nplwc
nplw7o = nplw7
npnod1o = npnod1
npnod2o = npnod2
npnodo = npnod
nspnodo= nspnod
nplw5exo = nplw5ex
nplwexo  = nplwex
nplw2exo = nplw2ex
nplw3exo = nplw3ex
nplwcexo = nplwcex
nplw7exo = nplw7ex
kmax1o = kmax1
kmax2o = kmax2
kmax3o = kmax3
kmax1do = kmax1d
kmax2do = kmax2d
kmax3do = kmax3d
rvolo = rvol

nplcnto(1:nodes) = nplcnt(1:nodes)
npldspo(1:nodes) = npldsp(1:nodes)

!-----store PW index
ngao(0:nplw5o) = nga(0:nplw5o)
ngbo(0:nplw5o) = ngb(0:nplw5o)
ngco(0:nplw5o) = ngc(0:nplw5o)


!-----noncollinear magnetism
!if( lnoncollinear ) then
!    mshnodo = mshnod
!    mftnodo(1:nodes) = mftnod(1:nodes)
!    mftdspo(1:nodes) = mftdsp(1:nodes)
!    call cpijkgdo( ijkgdo, ijkgd, nplw5o )
!    ntotfdo = ntotfd
!    mfd2fto(1:ntotfdo) = mfd2ft(1:ntotfdo)
!end if


!-----for k-point sampling
!call reclat_k_save( nfile, myid, nodes,  &
!& alloc_mem, nplw, pwscale )

licall = .false.


return
end subroutine




module wv_in_scratch
!-----------------------------------------------------------------------
! type declaration of shared variables in savewv_in_scratch and readwv_in_scratch
!-----------------------------------------------------------------------
implicit none

logical :: lalloc = .false.
real*8,  allocatable, dimension(:) ::  &
& sv_gdcr, sv_prveig, sv_gdcrsv, sv_rhgr, sv_prvrho

integer :: iunit
logical :: lopen = .false.
save

end module




subroutine savewv_in_scratch( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!     save wavefunctions : gdcr, gdcrsv, prveig
!       charge densities : rhgr, prvrho
!     in temporal variables or a scratch file
!-----------------------------------------------------------------------
use param
use param_atom
use pwlda_pp
use pwlda_atom
use pwlda_variables
use pwlda_proc
use pwlda_pw
use pwlda_grid
use pwreset_variables
use wv_in_scratch
implicit none

integer :: myid, nodes
integer, dimension(*) :: nfile

!-----declare local variables
integer :: status
integer :: ierror
integer :: i


ierror = 0
status = 0
!------allocate memory
if( lgamma ) then
    allocate( sv_gdcr(2*npnodo*nband), sv_prveig(2*npnodo*nband*nprvmx*nspnmx),  &
& stat=status )

    if( status == 0 ) then
        sv_gdcr(1:2*npnodo*nband) = gdcr(1:2*npnodo*nband)
        sv_prveig(1:2*npnodo*nband*nprvmx*nspnmx) = prveig(1:2*npnodo*nband*nprvmx*nspnmx)
    end if

    if( lspin .and. status == 0 ) then
        allocate( sv_gdcrsv(2*npnodo*nband),  &
& stat=status )

        if( status == 0 ) then
            sv_gdcrsv(1:2*npnodo*nband) = gdcrsv(1:2*npnodo*nband)
       end if
   end if
end if
if( status == 0 ) then
    allocate( sv_rhgr(nplw5exo*nspnmx2), sv_prvrho(nplw5exo*ncprvmx*nspnmx2),  &
& stat=status )

    if( status == 0 ) then
        sv_rhgr(1:nplw5exo*nspnmx2) = rhgr(1:nplw5exo*nspnmx2)
        sv_prvrho(1:nplw5exo*ncprvmx*nspnmx2) = prvrho(1:nplw5exo*ncprvmx*nspnmx2)
    end if
end if
lalloc = status == 0
if( lalloc ) return


call allocate_unit_number( iunit )

!-----Allocate I/0 file
open( iunit, status='scratch', form='unformatted',  &
&     iostat=ierror )

lopen = ierror == 0
if( .not.lopen ) return

if( lgamma ) then
   write(iunit) (   gdcr(i), i = 1, 2*npnodo*nband )
   write(iunit) ( prveig(i), i = 1, 2*npnodo*nband*nprvmx*nspnmx )
   if( lspin ) then
       write(iunit) ( gdcrsv(i), i = 1, 2*npnodo*nband )
   end if
end if
write(iunit) (   rhgr(i), i = 1, nplw5exo*nspnmx2 )
write(iunit) ( prvrho(i), i = 1, nplw5exo*ncprvmx*nspnmx2 )


return
end subroutine




subroutine readwv_in_scratch( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!     read wavefunctions : gdcr, gdcrsv, prveig
!       charge densities : rhgr, prvrho
!     in temporal variables or a scratch file
!-----------------------------------------------------------------------
use param
use param_atom
use pwlda_pp
use pwlda_atom
use pwlda_variables
use pwlda_proc
use pwlda_pw
use pwlda_grid
use pwreset_variables
use wv_in_scratch
implicit none

integer :: myid, nodes
integer, dimension(*) :: nfile

!-----declare local variables
integer :: ierror, istat, status
integer :: i


ierror = 0
if( lalloc ) then
    if( lgamma ) then
        gdcr(1:2*npnodo*nband) = sv_gdcr(1:2*npnodo*nband)
        prveig(1:2*npnodo*nband*nprvmx*nspnmx) = sv_prveig(1:2*npnodo*nband*nprvmx*nspnmx)
        if( lspin ) then
            gdcrsv(1:2*npnodo*nband) = sv_gdcrsv(1:2*npnodo*nband)
        end if
    end if

    rhgr(1:nplw5exo*nspnmx2) = sv_rhgr(1:nplw5exo*nspnmx2)
    prvrho(1:nplw5exo*ncprvmx*nspnmx2) = sv_prvrho(1:nplw5exo*ncprvmx*nspnmx2)

    !------deallocate memory
    if( lgamma ) then
        deallocate( sv_gdcr, sv_prveig, stat=status )
        if( lspin ) then
            deallocate( sv_gdcrsv, stat=status )
       end if
    end if
    deallocate( sv_rhgr, sv_prvrho, stat=status )
    lalloc = .false.
    return
end if


if( .not.lopen ) return


rewind(iunit)

istat = 0
if( lgamma ) then
    read(iunit,iostat=istat)  (   gdcr(i), i = 1, 2*npnodo*nband )

    if( istat == 0 )  &
&       read(iunit,iostat=istat)  &
&           ( prveig(i), i = 1, 2*npnodo*nband*nprvmx*nspnmx )

    if( lspin ) then
        if( istat == 0 )  &
&    read(iunit,iostat=istat) ( gdcrsv(i), i = 1, 2*npnodo*nband )
    end if
end if

if( istat == 0 )  &
& read(iunit,iostat=istat) (   rhgr(i), i = 1, nplw5exo*nspnmx2 )

if( istat == 0 )  &
& read(iunit,iostat=istat)  &
&                 ( prvrho(i), i = 1, nplw5exo*ncprvmx*nspnmx2 )


ierror = istat
close(iunit)
lopen = .false.

call deallocate_unit_number( iunit )


return
end subroutine




!!================== cut from here ======================================
!
!module magne_in_scratch
!!-----------------------------------------------------------------------
!! type declaration of shared variables in savewv_in_scratch and readwv_in_scratch
!!-----------------------------------------------------------------------
!implicit none
!
!logical :: lalloc = .false.
!real*8,  allocatable, dimension(:,:) :: sv_rhom
!real*8,  allocatable, dimension(:) :: sv_prvrhom
!
!integer :: iunit
!logical :: lopen = .false.
!save
!
!end module
!
!
!
!
!subroutine savemagne_in_scratch( nfile, myid, nodes, ierror )
!!-----------------------------------------------------------------------
!!     save magnetic moment : rhomx, rhomy, rhomz, prvrhom
!!     in temporal variables or a scratch file
!!-----------------------------------------------------------------------
!use param
!use param_atom
!use pwlda_pp
!use pwlda_atom
!use pwlda_variables
!use pwlda_proc
!use pwlda_pw
!use pwlda_grid
!use ncmagne_variables
!use pwreset_variables
!use magne_in_scratch
!implicit none
!
!integer :: myid, nodes
!integer, dimension(*) :: nfile
!
!!-----declare local variables
!integer :: status
!integer :: ierror
!integer :: i
!
!
!ierror = 0
!status = 0
!!------allocate memory
!allocate( sv_rhom(mshnodo,3), sv_prvrhom(mshnodo*ncprvmx*3),  &
!& stat=status )
!
!if( status == 0 ) then
!       sv_rhom(1:mshnodo,1)         =   rhomx(1:mshnodo)
!       sv_rhom(1:mshnodo,2)         =   rhomy(1:mshnodo)
!       sv_rhom(1:mshnodo,3)         =   rhomz(1:mshnodo)
!    sv_prvrhom(1:mshnodo*ncprvmx*3) = prvrhom(1:mshnodo*ncprvmx*3)
!end if
!lalloc = status == 0
!if( lalloc ) return
!
!
!call allocate_unit_number( iunit )
!
!!-----Allocate I/0 file
!open( iunit, status='scratch', form='unformatted',  &
!&     iostat=ierror )
!
!lopen = ierror == 0
!if( .not.lopen ) return
!
!write(iunit) rhomx(1:mshnodo)
!write(iunit) rhomy(1:mshnodo)
!write(iunit) rhomz(1:mshnodo)
!write(iunit) prvrhom(1:mshnodo*ncprvmx*3)
!
!
!return
!end subroutine
!
!
!
!
!subroutine readmagne_in_scratch( nfile, myid, nodes, ierror )
!!-----------------------------------------------------------------------
!!     read magnetic moment : rhomx, rhomy, rhomz, prvrhom
!!     in temporal variables or a scratch file
!!-----------------------------------------------------------------------
!use param
!use param_atom
!use pwlda_pp
!use pwlda_atom
!use pwlda_variables
!use pwlda_proc
!use pwlda_pw
!use pwlda_grid
!use ncmagne_variables
!use pwreset_variables
!use magne_in_scratch
!implicit none
!
!integer :: myid, nodes
!integer, dimension(*) :: nfile
!
!!-----declare local variables
!integer :: ierror, istat, status
!integer :: i
!
!
!ierror = 0
!if( lalloc ) then
!      rhomx(1:mshnodo)           =    sv_rhom(1:mshnodo,1)
!      rhomy(1:mshnodo)           =    sv_rhom(1:mshnodo,2)
!      rhomz(1:mshnodo)           =    sv_rhom(1:mshnodo,3)
!    prvrhom(1:mshnodo*ncprvmx*3) = sv_prvrhom(1:mshnodo*ncprvmx*3)
!
!    !------deallocate memory
!    deallocate( sv_rhom, sv_prvrhom, stat=status )
!    lalloc = .false.
!    return
!end if
!
!
!if( .not.lopen ) return
!
!
!rewind(iunit)
!
!istat = 0
!read(iunit,iostat=istat) rhomx(1:mshnodo)
!read(iunit,iostat=istat) rhomy(1:mshnodo)
!read(iunit,iostat=istat) rhomz(1:mshnodo)
!read(iunit,iostat=istat) prvrhom(1:mshnodo*ncprvmx*3)
!
!
!ierror = istat
!close(iunit)
!lopen = .false.
!
!call deallocate_unit_number( iunit )
!
!
!return
!end subroutine
!
!!================== cut up to here ======================================



subroutine rarrwv( nfile, myid, nodes, ltimecnt )
!-----------------------------------------------------------------------
!     rearrange wavefunctions : gdcr, gdcrsv, prveig
!            charge densities : rhgr, prvrho
!-----------------------------------------------------------------------
use param
use param_atom
use pwlda_pp
use pwlda_atom
use pwlda_variables
use pwlda_proc
use pwlda_pw
use pwlda_grid
use pwreset_variables
implicit none

integer :: myid, nodes
integer, dimension(*) :: nfile
logical :: ltimecnt

!-----declare local variables
integer :: nspin, i, j
real*8  :: ct0, c
integer :: nspnodx
integer :: kmx1, kmx2, kmx3
integer :: munit, munito
real*8  :: timecnt


!-----zero clear fft3x
fft3x = 0.d0
!-----zero clear fft3y
fft3y = 0.d0

kmx1 = max( kmax1o, kmax1 )
kmx2 = max( kmax2o, kmax2 )
kmx3 = max( kmax3o, kmax3 )


lgammaif: if( lgamma ) then

nspnodx = max( nspnodo, nspnod )
if( ltimecnt ) ct0 = timecnt()
!-----convert wavefuction : gdcr ---------------------------------------
if( lspin ) then
    nspin = 1
    call stspud( nspin, gdcr, gdcrsv, nspnodo, lcgjsv )
end if
do nspin = 1, nspnmx
!         if( lspin ) then
!   --- copy gdcr to rhcr
   call cpgdrh( nfile, myid, nodes, rhcr, gdcr, npnodo, nband )
!   --- to convert G decomposition to band decomposition ---
   call gdtobd( nfile, myid, nodes, ct0,  &
& cgjr, nplwexo, nplwo,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& rhcr, bufcr, npnod1o, npnod2o, npnodo, nplcnto, npldspo,  &
& iods, iodsg, idstnd, ltimecnt )
!         end if

   !-----convert rhcr to cgjr
   call storering( nfile, myid, nodes,  &
& rhcr, nplwo, nplwexo, ngao, ngbo, ngco,  &
& cgjr, nplw,  nplwex,  nga,  ngb,  ngc,  &
& fft3x, fft3y, kmx1, kmx2, kmx3, nbnod )

!    --- convert band decomposition to G decomposition ---
   call bdtogd( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iodg, idstnd, .false., ltimecnt )

   if( lspin ) then
       call stpsud( nspin, gdcr, gdcrsv, nspnodx, lcgjsv )
   end if
end do


!-----convert previous wavefuctions : prveig ---------------------------
munit  = 2*npnod*nband
munito = 2*npnodo*nband
if( munit > munito ) then
    !-----shift array
    do i = nprvmx*nspnmx, 2, -1
!       call shifta( prveig(munito*(i-1)+1), prveig(munit*(i-1)+1), munito )
       rhcr(1:munito) = prveig(munito*(i-1)+1:munito*(i-1)+munito)
       prveig(munit*(i-1)+1:munit*(i-1)+munito) = rhcr(1:munito)
    end do
    munito = munit
end if


if( ltimecnt ) ct0 = timecnt()
do nspin = 1, nspnmx
do j = 1, keypwv

   i = j + nprvmx*( nspin - 1 )
!   --- copy prveig to rhcr
   call cpgdrh( nfile, myid, nodes, rhcr, prveig(munito*(i-1)+1),  &
& npnodo, nband )
!   --- to convert G decomposition to band decomposition ---
   call gdtobd( nfile, myid, nodes, ct0,  &
& cgjr, nplwexo, nplwo,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& rhcr, bufcr, npnod1o, npnod2o, npnodo, nplcnto, npldspo,  &
& iods, iodsg, idstnd, ltimecnt )

   !-----convert rhcr to cgjr
   call storering( nfile, myid, nodes,  &
& rhcr, nplwo, nplwexo, ngao, ngbo, ngco,  &
& cgjr, nplw,  nplwex,  nga,  ngb,  ngc,  &
& fft3x, fft3y, kmx1, kmx2, kmx3, nbnod )

!    --- convert band decomposition to G decomposition ---
   call bdtogd( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& prveig(munit*(i-1)+1), bufcr,  &
& npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iodg, idstnd, .false., ltimecnt )

end do
end do

else lgammaif

    !-----for k-point sampling
!    call rarrwv_k( nfile, myid, nodes, ltimecnt,  &
!& nspnmx, nprvmx, keypwv, rhcr, cgjr, bufcr,  &
!& fft3x, fft3y, kmx1, kmx2, kmx3,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& iods, iodg, iodsg, idstnd )


end if lgammaif


!-----zero clear fft3x
fft3x = 0.d0
!-----zero clear fft3y
fft3y = 0.d0

kmx1 = max( kmax1do, kmax1d )
kmx2 = max( kmax2do, kmax2d )
kmx3 = max( kmax3do, kmax3d )

!-----convert charge density : rhgr ------------------------------------
munit  = nplw5ex
munito = nplw5exo
if( munit > munito ) then
    !-----shift array
    do i = nspnmx2, 2, -1
!       call shifta( rhgr(munito*(i-1)+1), rhgr(munit*(i-1)+1), munito )
       thrhgr(1:munito) = rhgr(munito*(i-1)+1:munito*(i-1)+munito)
       rhgr(munit*(i-1)+1:munit*(i-1)+munito) = thrhgr(1:munito)
    end do
    munito = munit
end if

do i = 1, nspnmx2
   !-----convert rhcr to cgjr
!   call storering( nfile, myid, nodes,  &
!& rhgr(munito*(i-1)+1), nplw5o, nplw5exo, ngao, ngbo, ngco,  &
!& rhgr(munit*(i-1)+1), nplw5,  nplw5ex,  nga,  ngb,  ngc,  &
!& fft3x, fft3y, kmx1, kmx2, kmx3, 1 )
   call storerhgr( nfile, myid, nodes,  &
& rhgr(munito*(i-1)+1), nplw5o, nplw5exo, ngao, ngbo, ngco,  &
& fft3x, fft3y, kmx1, kmx2, kmx3 )
   call restorerhgr( nfile, myid, nodes,  &
& rhgr(munit*(i-1)+1), nplw5,  nplw5ex,  nga,  ngb,  ngc,  &
& fft3x, fft3y, kmx1, kmx2, kmx3 )

   !-----correction by rvolo/rvol
   c = rvolo/rvol
   call abyc( rhgr(munit*(i-1)+1), nplw5ex, c )
end do

!-----convert previous charge density : prvrho -------------------------
munit  = nplw5ex
munito = nplw5exo
if( munit > munito ) then
    !-----shift array
    do i = ncprvmx*nspnmx2, 2, -1
!       call shifta( prvrho(munito*(i-1)+1), prvrho(munit*(i-1)+1), munito )
       thrhgr(1:munito) = prvrho(munito*(i-1)+1:munito*(i-1)+munito)
       prvrho(munit*(i-1)+1:munit*(i-1)+munito) = thrhgr(1:munito)
    end do
    munito = munit
end if

do nspin = 1, nspnmx2
do j = 1, keypcd

   i = j + ncprvmx*( nspin - 1 )
   !-----convert rhcr to cgjr
!   call storering( nfile, myid, nodes,  &
!& prvrho(munito*(i-1)+1), nplw5o, nplw5exo, ngao, ngbo, ngco,  &
!& prvrho(munit*(i-1)+1), nplw5,  nplw5ex,  nga,  ngb,  ngc,  &
!& fft3x, fft3y, kmx1, kmx2, kmx3, 1 )
   call storerhgr( nfile, myid, nodes,  &
& prvrho(munito*(i-1)+1), nplw5o, nplw5exo, ngao, ngbo, ngco,  &
& fft3x, fft3y, kmx1, kmx2, kmx3 )
   call restorerhgr( nfile, myid, nodes,  &
& prvrho(munit*(i-1)+1), nplw5,  nplw5ex,  nga,  ngb,  ngc,  &
& fft3x, fft3y, kmx1, kmx2, kmx3 )

   !-----correction by rvolo/rvol
   c = rvolo/rvol
   call abyc( prvrho(munit*(i-1)+1), nplw5ex, c )
end do
end do


return
end subroutine




!!================== cut from here ======================================
!
!subroutine ncrarrmagne( nfile, myid, nodes, ltimecnt )
!!-----------------------------------------------------------------------
!!     rearrange magnetic moment : rhomx, rhomy, rhomz, prvrhom
!!-----------------------------------------------------------------------
!use param
!use param_atom
!use pwlda_pp
!use pwlda_atom
!use pwlda_variables
!use pwlda_proc
!use pwlda_pw
!use pwlda_grid
!use ncmagne_variables
!use pwreset_variables
!implicit none
!integer :: myid, nodes
!integer, dimension(*) :: nfile
!logical :: ltimecnt
!
!!-----declare local variables
!integer :: kmx1, kmx2, kmx3, i, j, k, munit, munito
!
!
!if( .not.lnoncollinear ) return
!
!kmx1 = max( kmax1do, kmax1d )
!kmx2 = max( kmax2do, kmax2d )
!kmx3 = max( kmax3do, kmax3d )
!
!!---rhomx
!call ncrarrmagne2( nfile, myid, nodes,  &
!& rhomx, glocal, rinr, thrhgr, nplw5exo, nplw5o, ngao, ngbo, ngco,  &
!& mshnodo, mftnodo, mftdspo, mfd2fto, ntotfdo, nd1vks_old, ijkgdo,  &
!& nplw5ex, nplw5, nga, ngb, ngc,  &
!& mshnod(1), mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d, ijkgd,  &
!& kmx1, kmx2, kmx3, fft3x, fft3y, fftwork )
!
!!---rhomy
!call ncrarrmagne2( nfile, myid, nodes,  &
!& rhomy, glocal, rinr, thrhgr, nplw5exo, nplw5o, ngao, ngbo, ngco,  &
!& mshnodo, mftnodo, mftdspo, mfd2fto, ntotfdo, nd1vks_old, ijkgdo,  &
!& nplw5ex, nplw5, nga, ngb, ngc,  &
!& mshnod(1), mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d, ijkgd,  &
!& kmx1, kmx2, kmx3, fft3x, fft3y, fftwork )
!
!!---rhomz
!call ncrarrmagne2( nfile, myid, nodes,  &
!& rhomz, glocal, rinr, thrhgr, nplw5exo, nplw5o, ngao, ngbo, ngco,  &
!& mshnodo, mftnodo, mftdspo, mfd2fto, ntotfdo, nd1vks_old, ijkgdo,  &
!& nplw5ex, nplw5, nga, ngb, ngc,  &
!& mshnod(1), mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d, ijkgd,  &
!& kmx1, kmx2, kmx3, fft3x, fft3y, fftwork )
!
!!---prvrhom
!munit  = mshnod(1)
!munito = mshnodo
!if( munit > munito ) then
!    !-----shift array
!    do i = keypmg*3, 2, -1
!!       call shifta( prvrho(munito*(i-1)+1), prvrho(munit*(i-1)+1), munito )
!       dbuf(1:munito) = prvrhom(munito*(i-1)+1:munito*(i-1)+munito)
!       prvrhom(munit*(i-1)+1:munit*(i-1)+munito) = dbuf(1:munito)
!    end do
!    munito = munit
!end if
!do k = 1, keypmg
!do j = 1, 3
!   i = j + 3*(k-1)
!   dbuf(1:munito) = prvrhom(munito*(i-1)+1:munito*(i-1)+munito)
!   !---prvrhom
!   call ncrarrmagne2( nfile, myid, nodes,  &
!& dbuf, glocal, rinr, thrhgr, nplw5exo, nplw5o, ngao, ngbo, ngco,  &
!& mshnodo, mftnodo, mftdspo, mfd2fto, ntotfdo, nd1vks_old, ijkgdo,  &
!& nplw5ex, nplw5, nga, ngb, ngc,  &
!& mshnod(1), mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d, ijkgd,  &
!& kmx1, kmx2, kmx3, fft3x, fft3y, fftwork )
!   prvrhom(munit*(i-1)+1:munit*(i-1)+munit) = dbuf(1:munit)
!end do
!end do
!
!
!return
!end subroutine
!
!
!
!
!subroutine ncrarrmagne2( nfile, myid, nodes,  &
!& rhom, glocal, rinr, thrhgr, nplw5exo, nplw5o, ngao, ngbo, ngco,  &
!& mshnodo, mftnodo, mftdspo, mfd2fto, ntotfdo, nd1vks_old, ijkgdo,  &
!& nplw5ex, nplw5, nga, ngb, ngc,  &
!& mshnod, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d, ijkgd,  &
!& kmx1, kmx2, kmx3, fft3x, fft3y, fftwork )
!!-----------------------------------------------------------------------
!!     rearrange magnetic moment : rhomx, rhomy, rhomz, prvrhom
!!-----------------------------------------------------------------------
!implicit none
!integer :: myid, nodes
!integer, dimension(*) :: nfile
!real*8  :: rhom(*), glocal(*), rinr(*), thrhgr(*)
!integer :: nplw5exo, nplw5o, ngao(0:nplw5o), ngbo(0:nplw5o), ngco(0:nplw5o)
!integer :: mshnodo, mftnodo(*), mftdspo(*), mfd2fto(*), ntotfdo, nd1vks_old(*)
!integer :: ijkgdo(-nplw5o:nplw5o)
!integer :: nplw5ex, nplw5, nga(0:nplw5), ngb(0:nplw5), ngc(0:nplw5)
!integer :: mshnod, mftnod(*), mftdsp(*), mfd2ft(*), ntotfd, nd1vks(*), kfft0d
!integer :: ijkgd(-nplw5:nplw5)
!integer :: kmx1, kmx2, kmx3
!real*8  :: fft3x(*), fft3y(*)
!complex*16 :: fftwork(*)
!
!!-----declare local variables
!integer :: kfft0do
!
!
!!---unify rhom
!call alldgatherv(rhom,mshnodo,glocal,mftnodo,mftdspo)
!
!kfft0do = nd1vks_old(1)*nd1vks_old(2)*nd1vks_old(3)
!!--- charge density transformation from r- to g-spaces
!call chgr2g( nfile, myid, nodes,  &
!& glocal, rinr, nplw5exo, nplw5o, mfd2fto, ntotfdo, nd1vks_old, kfft0do,  &
!& fft3x, fft3y, fftwork, ijkgdo, nplw5o )
!
!call storering2( nfile, myid, nodes,  &
!& rinr, nplw5o, nplw5exo, ngao, ngbo, ngco,  &
!& rinr, nplw5,  nplw5ex,  nga,  ngb,  ngc,  &
!& fft3x, fft3y, kmx1, kmx2, kmx3, 1 )
!
!!--- charge density transformation from g- to r-spaces
!call chgg2r( nfile, myid, nodes,  &
!& glocal, rinr, nplw5ex, nplw5, mfd2ft, ntotfd, nd1vks, kfft0d,  &
!& fft3x, fft3y, fftwork, ijkgd, nplw5 )
!
!!--- store local charge density
!call distlc( nfile, myid, nodes,  &
!& rhom, mshnod, glocal, mftdsp )
!
!
!return
!end subroutine
!
!!================== cut up to here ======================================




subroutine storering( nfile, myid, nodes,  &
& rhcr, nplwo, nplwexo, ngao, ngbo, ngco,  &
& gdcr, nplw,  nplwex,  nga,  ngb,  ngc,  &
& fft3x, fft3y, kmx1, kmx2, kmx3, nbnod )
!-----------------------------------------------------------------------
!     convert rhcr to gdcr
!-----------------------------------------------------------------------
use param
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: nbnod
integer :: nplwo, nplwexo
real*8,  dimension(nplwexo,*) :: rhcr
integer, dimension(0:nplwo) :: ngao, ngbo, ngco
integer :: nplw, nplwex
real*8,  dimension(nplwex,*) :: gdcr
integer, dimension(0:nplw) :: nga, ngb, ngc
integer :: kmx1, kmx2, kmx3
real*8, dimension(-kmx1:kmx1,-kmx2:kmx2,-kmx3:kmx3) :: fft3x,fft3y


!if( .not.lnoncollinear ) then
    call storering2( nfile, myid, nodes,  &
& rhcr, nplwo, nplwexo, ngao, ngbo, ngco,  &
& gdcr, nplw,  nplwex,  nga,  ngb,  ngc,  &
& fft3x, fft3y, kmx1, kmx2, kmx3, nbnod )
!else
!    !-----noncollinear magnetism
!    if( lgamma ) then
!        call ncstorering2( nfile, myid, nodes,  &
!& rhcr, nplwo, nplwexo*2, ngao, ngbo, ngco,  &
!& gdcr, nplw,  nplwex*2,  nga,  ngb,  ngc,  &
!& fft3x, fft3y, kmx1, kmx2, kmx3, nbnod )
!    else
!        call ncstorering2_k( nfile, myid, nodes,  &
!& rhcr, nplwo, nplwexo*2, ngao, ngbo, ngco,  &
!& gdcr, nplw,  nplwex*2,  nga,  ngb,  ngc,  &
!& fft3x, fft3y, kmx1, kmx2, kmx3, nbnod )
!    end if
!end if


return
end subroutine




subroutine storering2( nfile, myid, nodes,  &
& rhcr, nplwo, nplwexo, ngao, ngbo, ngco,  &
& gdcr, nplw,  nplwex,  nga,  ngb,  ngc,  &
& fft3x, fft3y, kmx1, kmx2, kmx3, nbnod )
!-----------------------------------------------------------------------
!     convert rhcr to gdcr
!-----------------------------------------------------------------------
implicit none

integer :: myid, nodes
integer, dimension(*) :: nfile

integer :: nbnod
integer :: nplwo, nplwexo
real*8,  dimension(nplwexo,*) :: rhcr
integer, dimension(0:nplwo) :: ngao, ngbo, ngco
integer :: nplw, nplwex
real*8,  dimension(nplwex,*) :: gdcr
integer, dimension(0:nplw) :: nga, ngb, ngc
integer :: kmx1, kmx2, kmx3
real*8, dimension(-kmx1:kmx1,-kmx2:kmx2,-kmx3:kmx3) :: fft3x,fft3y

!-----declare local variables
integer :: ig, ib, ix, iy, iz


do ib = 1, nbnod

   fft3x(:,:,:) = 0.d0
   fft3y(:,:,:) = 0.d0
   do ig = 0, nplwo
      ix = ngao(ig)
      iy = ngbo(ig)
      iz = ngco(ig)
      fft3x(ix,iy,iz) = rhcr(2*ig+1,ib)
      fft3y(ix,iy,iz) = rhcr(2*ig+2,ib)
   end do

   do ig = 0, nplw
      ix = nga(ig)
      iy = ngb(ig)
      iz = ngc(ig)
      gdcr(2*ig+1,ib) = fft3x(ix,iy,iz)
      gdcr(2*ig+2,ib) = fft3y(ix,iy,iz)
   end do

end do


return
end subroutine




subroutine storerhgr( nfile, myid, nodes,  &
& rhcr, nplwo, nplwexo, ngao, ngbo, ngco,  &
& fft3x, fft3y, kmx1, kmx2, kmx3 )
!-----------------------------------------------------------------------
!     convert rhcr to gdcr
!-----------------------------------------------------------------------
implicit none

integer :: myid, nodes
integer, dimension(*) :: nfile

integer :: nplwo, nplwexo
real*8,  dimension(nplwexo) :: rhcr
integer, dimension(0:nplwo) :: ngao, ngbo, ngco
integer :: kmx1, kmx2, kmx3
real*8, dimension(-kmx1:kmx1,-kmx2:kmx2,-kmx3:kmx3) :: fft3x,fft3y

!-----declare local variables
integer :: ig, ix, iy, iz


fft3x(:,:,:) = 0.d0
fft3y(:,:,:) = 0.d0
do ig = 0, nplwo
   ix = ngao(ig)
   iy = ngbo(ig)
   iz = ngco(ig)
   fft3x(ix,iy,iz) = rhcr(2*ig+1)
   fft3y(ix,iy,iz) = rhcr(2*ig+2)
end do


return
end subroutine




subroutine restorerhgr( nfile, myid, nodes,  &
& gdcr, nplw,  nplwex,  nga,  ngb,  ngc,  &
& fft3x, fft3y, kmx1, kmx2, kmx3 )
!-----------------------------------------------------------------------
!     convert rhcr to gdcr
!-----------------------------------------------------------------------
implicit none

integer :: myid, nodes
integer, dimension(*) :: nfile

integer :: nplw, nplwex
real*8,  dimension(nplwex) :: gdcr
integer, dimension(0:nplw) :: nga, ngb, ngc
integer :: kmx1, kmx2, kmx3
real*8, dimension(-kmx1:kmx1,-kmx2:kmx2,-kmx3:kmx3) :: fft3x,fft3y

!-----declare local variables
integer :: ig, ix, iy, iz


do ig = 0, nplw
   ix = nga(ig)
   iy = ngb(ig)
   iz = ngc(ig)
   gdcr(2*ig+1) = fft3x(ix,iy,iz)
   gdcr(2*ig+2) = fft3y(ix,iy,iz)
end do


return
end subroutine


!subroutine shifta( a, b, n )
!!-----copy a to b
!implicit none
!integer :: n
!real*8,  dimension(n) :: a, b
!integer :: i
!
!do i = n, 1, -1
!   b(i) = a(i)
!end do
!
!return
!end subroutine


subroutine abyc( a, n, c )
!-----copy a to b
implicit none
integer :: n
real*8,  dimension(n) :: a
real*8  :: c
integer :: i

do i = 1, n
   a(i) = a(i)*c
end do

return
end subroutine
