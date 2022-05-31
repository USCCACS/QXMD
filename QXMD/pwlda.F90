



subroutine pwlda_set( nfile, myid, nodes,  &
& npx_, npy_, npz_, myid_md, npx_md, npy_md, npz_md, iogpsz_, ierror )
!-----------------------------------------------------------------------
!     initial setting for LDA or GGA calculation with plane-wave method
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
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes
integer :: npx_, npy_, npz_, myid_md, npx_md, npy_md, npz_md
integer :: iogpsz_
integer :: ierror


ct0  = timecnt()

iogpsz  = iogpsz_

!-----Start parallel environment and keep my node ID --
npx = npx_
npy = npy_
npz = npz_
myx=myid/(npy*npz)
myy=mod(myid/npz,npy)
myz=mod(myid,npz)
!---  Prepares a neighbor-node-ID table --------------------------------
call ntset( myx, myy, myz, nn, myparity, npx, npy, npz )


!-----------------------------------------------------------------------
!----- read and check input data
call input( nfile, myid, nodes, lclust )

if( lspin ) then
    nspnmx = 2
  else
    nspnmx = 1
end if


ct = timecnt()
if(loutfile(1)) write(nfile(1),*) ' read input data    : cpu-time :', ct - ct0
if(loutfile(2)) write(nfile(2),*) ' read input data    : cpu-time :', ct - ct0
ct0 = ct
!-----------------------------------------------------------------------


!--- set indices of k points for each node
!call setkp( nfile,  &
!& nkpnt, nknod1, nknod2, nknod, bzk, lgammak, lgamma, kpsymop )


!-----allocate memory for variables for pseudopotentials
call pwlda_pp_alloc( nfile, myid, nodes,  &
& alloc_mem, ntype, mxl, mx1, mx1loc, lpcc )


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

ct = timecnt()
if(loutfile(1)) write(nfile(1),*) ' set reciprocal vectors : cpu-time :',  &
&                      ct - ct0
if(loutfile(2)) write(nfile(2),*) ' set reciprocal vectors : cpu-time :',  &
&                      ct - ct0
ct0 = ct


!if( ldouble_grid_recip ) then

    !--- set reciprocal vectors for the double-grid method
!    call reclat_double_grid( nfile, myid, nodes,  &
!& alloc_mem, ecutlong )

    !--- initialize the double-grid method
!    call ini_double_grid( nfile, myid, nodes, alloc_mem,  &
!& dgalpha, ecutlong, idouble_grid_method, lvacuum, nd1vks,  &
!& alpha_ldouble_grid_recip, alpha_llcexchange )

!if( ldouble_grid_recip ) then
!
!    ct = timecnt()
!    if(loutfile(1)) write(nfile(1),*) ' set double-grid method : cpu-time :',  &
!&                      ct - ct0
!    if(loutfile(2)) write(nfile(2),*) ' set double-grid method : cpu-time :',  &
!&                      ct - ct0
!    ct0 = ct
!
!end if


!-----allocate memory for planewave-decomposition variables
call planewave_decomp_variables_alloc( nfile, alloc_mem )

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


!-----allocate memory for processor-dependent variables
call pwlda_proc_alloc( nfile, myid, nodes, alloc_mem )


!-----allocate memory for iod, iodg, iods, iodsg
call pwlda_iod_alloc( nfile, myid, nodes,  &
& alloc_mem, nband )

!-----set indexes for band decomposition: nbnod, nbncnt, nbndsp
call decompose( nfile, myid, nodes,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp )

!-----set indexes for band decomposition: iod, iodg, iods, iodsg
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


!-----set destination nodes for all-to-all communication
call dstata( nfile, myid, nodes, idstnd )


!-----allocate memory for variables for atoms
call pwlda_atom_alloc( nfile, myid, nodes )
!&, ntype, nhk1, nhk2, lkbpp_g, lkbpp_r, lvand_r, mxl, nylmmx,  &
!& nbnod, nspnmx, jhybrid, lefield ) !, lefield_islts )


!--- set No. of atoms in each node
call statom( nfile, myid, nodes,  &
& ntype, nhk1, nhk2, nhk1_nod, nhk2_nod, iatoit, nion_nod, natom )

!-----allocate memory for iatmpt_nod
call iatmpt_nod_alloc( nfile, myid, nodes,  &
& alloc_mem, ntype, nhk1, nhk2 )

!-----set indexes for atom decomposition: natnod, natcnt, natdsp
call decompose( nfile, myid, nodes,  &
& natom, natnod1, natnod2, natnod, natcnt, natdsp )

!-----set indexes for atom decomposition: nhk*_nat, ioa*
call atom_decomp( nfile, myid, nodes,  &
& natom, natnod1, natnod2, natnod, natcnt, natdsp, iatoit,  &
& ntype, nhk1, nhk2, nhk1_nat, nhk2_nat, ioa, ioag, ioas, ioasg )


!-----allocate memory for the plane wave method
call pwlda_alloc( nfile, myid, nodes )


!-----allocate memory for nlmisc
!call nlmisc_alloc( nfile, myid, nodes,  &
!& alloc_mem, nband, nbnod, ntype, nhk1, nhk2, natom,  &
!& natnod1, natnod2, natnod, natcnt, natdsp, nhk1_nat, nhk2_nat,  &
!& ioa, ioag, ioas, ioasg )

!-----allocate memory for prvslmir
call slmir_alloc( nfile, myid, nodes,  &
& alloc_mem, nband, nbnod, ntype, natom, natnod,  &
& lkbpp, lvandi, lvand, lspin, nprvmx, nspnmx, laspc )

!if( lpaw ) then
!    !-----allocate memory for prvrhoij
!    call prvrhoij_variables_alloc( nfile, myid, nodes,  &
!& alloc_mem, ntype, natom, natnod1, natnod2, natnod, natcnt, natdsp,  &
!& iatoit, ioa, lvandi, nspnmx,  &
!& lplusU, lplusU_at, plusU_l, plusU_U, plusU_J )
!end if


!-----allocate memory for TDDFT-FSSH
call tddft_fssh_alloc( nfile, myid, nodes,  &
& alloc_mem, nband, nbnod, noband, nspnmx, natom )


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


!-----allocate memory for Energy Density Analysis (EDA)
!if( leda ) then
!    call eda_alloc( nfile, myid, nodes,  &
!& alloc_mem, leda, ntype, nhk1, nhk2 )
!end if

!-----allocate memory for noncollinear magnetism
!if( lnoncollinear ) then
!    call ncmagne_alloc( nfile, myid, nodes,  &
!& alloc_mem, ntype, nhk1, nhk2, zatom )
!end if


!----- Allocate I/0 files
call open_files( nfile, myid, nodes, ierror )

call title_in_files( nfile, myid, nodes, ierror )


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
!        !-----noncollinear magnetism
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


!-----initial occupancies
if( lspin ) then
    call iniocc_spin( nfile, myid, nodes,  &
& nel, nocc, nband, occ, norder, wegud, nsordr, nkpnt, wbzk )
    egvlud(:) = 0.d0
  else
    call iniocc( nfile, myid, nodes,  &
& nel, nocc, nband, occ, norder, nkpnt, wbzk )
    eig(:) = 0.d0
end if



!if( lvdw ) then
!    !----- set van der Waals kernel table : tbkernel
!    call setvdWkernel( nfile )
!
!    !----- set van der Waals kernel table : tbkernel
!    call setvdWfactorization( nfile,  &
!& recnrm, nplw5, gx, gy, gz, nd1vks, lstress )
!end if


!-----setup DFT-D
if( ldftd ) call setdftd( nfile, myid, nodes, &
jgga, ntype, nhk1, nhk2, zatom, evdj, hrdev, audang, avogad )

!-----setup LR-TDDFT
!if( lrtddft ) call lr_tddft_setup( nfile, myid, nodes )


ct = timecnt()
if(loutfile(1)) write(nfile(1),*) ' set parameters     : cpu-time :', ct - ct0
if(loutfile(2)) write(nfile(2),*) ' set parameters     : cpu-time :', ct - ct0
ct0 = ct


!-----------------------------------------------------------------------
!-----get multigrid level
!gridif: if( ldouble_grid ) then

!    nultg = multg
!    nnd2v = nd2v
!    multg = 1
!
!    !-----allocate memory for multigrid
!    call pwlda_grid_alloc( nfile, myid, nodes,  &
!& alloc_mem, multg, nultg )
!
!    call getmultg(  lclust, multg, nd1vks, nd2v, mulnd2 )
!    call getsmultg( lclust, multg, nd1vks, muls, nd1vs,  &
!&                   mulnd2, kulnd2, msrhmx, msrhmy, msrhmz )
!
!    !-----allocate memory for multigrid
!    call pwlda_mgrid_alloc( nfile, myid, nodes, npx, npy, npz,  &
!& alloc_mem, nd1vks, nd2v, multg, nultg, lspin, pwscale, lrtddft,  &
!& ltddft_fssh, lnoncollinear, ncprvmx, lwell, lefield, lsawtooth )
!
!    !-----get multigrid level for second grid
!    call getmultg(  lclust, nultg, nd1v, nnd2v, nulnd2 )
!    call getsmultg( lclust, nultg, nd1v, lmuls, lnd1vs,  &
!&                   nulnd2, lulnd2, msrhmx, msrhmy, msrhmz )
!
!    !-----allocate memory for second multigrid
!    call poisson_sec_alloc( nfile, myid, nodes, npx, npy, npz,  &
!& alloc_mem, nd1v, nnd2v, nultg, nulnd2 )
!
!    !-----allocate memory for serial calculation
!    call poisson_serial_alloc( nfile, myid, nodes,  &
!& alloc_mem, lnd1vs, nultg, lmuls, lulnd2 )

!else gridif

    nultg = 1

    !-----allocate memory for multigrid
    call pwlda_grid_alloc( nfile, myid, nodes,  &
& alloc_mem, multg, nultg )

    call getmultg(  lclust, multg, nd1vks, nd2v, mulnd2 )
    call getsmultg( lclust, multg, nd1vks, muls, nd1vs,  &
&                   mulnd2, kulnd2, msrhmx, msrhmy, msrhmz )

    !-----allocate memory for multigrid
    call pwlda_mgrid_alloc( nfile, myid, nodes, npx, npy, npz,  &
& alloc_mem, nd1vks, nd2v, multg, multg, lspin, pwscale, lrtddft,  &
& ltddft_fssh, lnoncollinear, ncprvmx, lwell, lefield, lsawtooth )

    !-----allocate memory for serial calculation
!    call poisson_serial_alloc( nfile, myid, nodes,  &
!& alloc_mem, nd1vs, multg, muls, kulnd2 )

!end if gridif


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


!if( .not.ldouble_grid ) then

    !-----set mesh points for serial calculation
!    call setsmsh( nfile, myid, nodes,  &
!& lclust, lsphere, multg, muls, nd1vs, kulnd2, rdel, rdelv, rmax,  &
!& nhit1, cofmas, lvacuum, amcdv2, ndv2x3, multg )

!else

    !-----set mesh points for double grid method in real space
!    do i = 1, 2
!    if( loutfile(i) ) then
!        write(nfile(i),*) ' '
!        write(nfile(i),*) ' << multigrid in double-grid method',  &
!& ' in r-space >>'
!    end if
!    end do
!    call setmpo( nfile, myid, nodes, npx, npy, npz,  &
!& myx, myy, myz, nn, myparity,  &
!& lclust, lsphere, nultg, nd1v, nnd2v, lmuls, lnd1vs, nulnd2,  &
!& lulnd2, rdel, rdelv, nhit1, rmax, cofmas, disp_org, ndisp_org,  &
!& meshx1, meshx, meshy1, meshy, meshz1, meshz )

!end if


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


!-----for the cardinal B-spline in the double-grid method
!if( ldouble_grid_recip ) then
!    call setmsh_bsp( nfile, myid, nodes,  &
!& alloc_mem, myx, myy, myz, npx, npy, npz,  &
!& mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1) )
!end if
!-----------------------------------------------------------------------
ct = timecnt()
if(loutfile(1)) write(nfile(1),*) ' set mesh points    : cpu-time :', ct - ct0
if(loutfile(2)) write(nfile(2),*) ' set mesh points    : cpu-time :', ct - ct0
ct0 = ct


return
end




subroutine pwlda_set2( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
! initial setting for LDA or GGA calculation with plane-wave method
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
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes
integer :: ierror

!------declare local variables
logical :: lstart_
logical :: leig_start, leig_end


ct0  = timecnt()


lstart_ = lstart
!=======================================================================
!    restore atomic configuration
call rdacon( nfile, iogpsz, ct0,  &
& fname_ion, lstart, prevr, nstepCG, nstepMD,  &
& ratm, frc, ntype, nhk1, nhk2, nhk1r, nhk2r, lclust, ifmd,  &
& natom, hcell )
!=======================================================================


!=======================================================================
!--- set initial wave functions : cgjr
!if( lgamma ) then
    call rdwfns( nfile, iogpsz, ct0,  &
& fname_eig, lstart, cgjr, eig, egvlud, occ, wegud, nplwex, nplw,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, lspin, gdcrsv, nspnod, lcgjsv, lwfrand, ncscale, wvratio )
!else
!
!    if( lstart ) then
!        call rd_k_wfns( nfile,  &
!& fname_eigk, lstart, bzk, ierror,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, eig, egvlud,  &
!& rhcr, cgjr, nspnmx, lspin, ncscale, wvratio )
!        if( ierror /= 0 ) then
!            !---Try to read data of Gamma-point calculation
!            call rdwfns( nfile, iogpsz, ct0,  &
!& fname_eig, lstart, cgjr, eig, egvlud, occ, wegud, nplwex, nplw,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& gdcr, lspin, gdcrsv, nspnod, lcgjsv, lwfrand, ncscale, wvratio )
!            call copywfns( nfile, myid, nodes,  &
!& cgjr, nplwex, nplw, nbnod, gdcrsv, lspin, ncscale )
!            call copyeig( nfile, myid, nodes,  &
!& eig, egvlud, nband, lspin )
!        end if
!        ierror = 0
!      else
!        call st_k_wfns( nfile, lspin, lwfrand, lcgjsv,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, gdcr, cgjr, ncscale )
!    end if
!
!end if


!--- set initial wave packets in TDDFT
!call rdwfns_tddft( nfile, ct0,  &
!fname_tddft, lstart, cgjr, nplwex, nplw, nbnod, gdcrsv, lspin )

!--- set data for TDDFT-FSSH
if( lspin ) then
    call read_tddft_fssh( nfile, iogpsz, wegud, nstepMD,  &
& gdcr, rhcr, bufcr, iodg, ioag, idstnd,  &
& nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& npnod1, npnod2, npnod, nplcnt, npldsp, nspnmx )
  else
    call read_tddft_fssh( nfile, iogpsz, occ, nstepMD,  &
& gdcr, rhcr, bufcr, iodg, ioag, idstnd,  &
& nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& npnod1, npnod2, npnod, nplcnt, npldsp, nspnmx )
end if

!--- set data for uniform electric field calculation
!call read_efield_data( nfile )

!--- read wavefunctions at previous time steps : prveig
keypwv = 0
!if( lgamma ) then
    call redpwv( nfile, iogpsz, ct0,  &
& fname_peig, lstart, ifmd, nstepCG, nstepMD, keypwv,  &
& gdcr, rhcr, bufcr, iodg, ioag, idstnd,  &
& nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& prveig, npnod1, npnod2, npnod, nplcnt, npldsp, nprvmx, nspnmx,  &
& prod, prodr, pdbuf, dmtrxr, nbxxxx,  &
& node_c, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, bijr, bm1r, ncscale )
!else
!    call redpwv_k( nfile,  &
!& fname_peigk, lstart, ifmd, nstepCG, nstepMD, keypwv, bzk,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& gdcr, rhcr, iodg, ioag, idstnd, bufcr, nprvmx, nspnmx, ncscale )
!end if
!=======================================================================

!--- set occupancies
if( lstart ) then
if( .not.ltddft_fssh .or. lfssh_gsscf ) then

    if( lspin ) then
        call fermi_k_spin( nfile, myid, nodes,  &
& nband, nkpnt, noband, nel, egvlud, nsordr, wegud, wbzk, nocc,  &
& eig, norder, occ,  &
& entrpy, feneud, fermie, lfermi, tfermi, diffud, lfixud )
    else
        call fermi_k( nfile, myid, nodes,  &
& nband, nkpnt, noband, nel, eig, norder, occ, wbzk, nocc,  &
& entrpy, fermie, lfermi, tfermi )
    end if

    !---tetrahedron method for BZ-zone sampling
!    if( ltetra ) then
!        if( lspin ) then
!            call BZtetrahedron_spin( nfile, myid, nodes,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, nkpnt, noband,  &
!& nel, egvlud, nsordr, wegud, wbzk, eig, norder, occ,  &
!& entrpy, feneud, fermie, diffud, lfixud, ekib, prod, prodr )
!        else
!            call BZtetrahedron( nfile, myid, nodes,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, nkpnt, noband,  &
!& nel, eig, norder, occ, wbzk, entrpy, fermie, ekib, prod, prodr, lnoncollinear )
!        end if
!    end if

end if
end if


if( ltddft_fssh ) then
    !---TDDFT-FSSH  initial occupations
    if( lspin ) then
        call tddft_fssh_ini( nfile, myid, nodes,  &
nband, nbnod, nspnmx, egvlud, nsordr, wegud, entrpy )
      else
        call tddft_fssh_ini( nfile, myid, nodes,  &
nband, nbnod, nspnmx, eig, norder, occ, entrpy )
    end if
end if


!=======================================================================
keypcd = 0
keypmg = 0
!--- initial electron density : rho
call rdcdty( nfile, iogpsz, ct0,  &
& fname_cds, lstart, hcell, h_MD, lorthrhmbc, rho, rdelv, rdel, nel,  &
& multg, lclust, nd1v, ntype, nhk1, nhk2, zv, lclno, ratm,  &
& natom, mx1,  &
& mshnod, mshnx, mshny, mshnz, mshx1, mshy1, mshz1,  &
& mshx, mshy, mshz, mulpms, noddatx,  &
& lspin, nspnmx, rhoud, diffud, lfixud,  &
& glocal, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d,  &
& rhgr, rinr, nplw5ex, nplw5, fft3x, fft3y, fftwork, ijkgd, nplw5,  &
& nga, ngb, ngc, lnoncollinear, keypcd, ncprvmx )


!--- read charge densities at previous time steps
call redpcd( nfile,  &
& fname_pcds, lstart, ifmd, nstepCG, nstepMD,  &
& keypcd, prvrho, nplw5ex, nplw5, ncprvmx, lspin, nspnmx2, nplw5,  &
& lnoncollinear, keypmg,  &
& glocal, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d, fft3x )
!=======================================================================


!-----deallocate memory for pseudopotentials
if( .not.lpaw .and. .not.lvshape .and. .not.lsphexp ) then
    call psvariables_dealloc( nfile, myid, nodes,  &
& alloc_mem, mxl, ntype )
end if


!-----load hartree potential
call rdvhar( nfile,  &
& fname_hrt, lstart, vhar, mshnod(1), ldouble_grid,  &
& glocal, mftnod, mftdsp, mfd2ft, ntotfd, kfft0d, fft3x )


!-----------------------------------------------------------------------
!if( lclust ) then
!    !-----initialize local pseudopotential filtering
!    call localpini( nfile, myid, nodes,  &
!& ntype, zv, ecutdens, alloc_mem, mx1loc, llclpp_r, llclpp_g,  &
!& tberf, tberfa, drcut, tbeff, tbeffa, drcutf,  &
!& tablc, tablca, dltlc, rmxlc, tbflc, tbflca, dltflc, rmxflc,  &
!& rctflc, llking, rlking, glkgmax, glkgexct, hcell, gamma,  &
!& rccc2, lhfull, vlocli ,xitgrd, rr, nvlcl, nion_nod,  &
!& mshx(1), mshy(1), mshz(1) )
!  else
    !-----initialize Ewald method & shape denpendent terms
    call ewaldini( nfile, myid, nodes,  &
& ntype, zv, ecutdens, alloc_mem, mx1loc, llclpp_r, llclpp_g,  &
& tberf, tberfa, drcut, tbeff, tbeffa, drcutf,  &
& tablc, tablca, dltlc, rmxlc, tbflc, tbflca, dltflc, rmxflc,  &
& rctflc, llking, rlking, glkgmax, glkgexct, hcell, h_MD, gamma,  &
& rccc2, lhfull, vlocli ,xitgrd, rr, nvlcl, nion_nod,  &
& mshx(1), mshy(1), mshz(1),  &
& ldouble_grid_recip, nd1v, nd1vks, dgalpha, pwscale, lvshape )
!end if
!      if( lshdep ) call shpset( lvacuum, volume, cshdep )
!-----------------------------------------------------------------------
ct = timecnt()
if(loutfile(1)) write(nfile(1),*) ' initialize Ewald   : cpu-time :', ct - ct0
if(loutfile(2)) write(nfile(2),*) ' initialize Ewald   : cpu-time :', ct - ct0
ct0 = ct




!-----memory allocation for variables for the calculations

!-----in reciprocal space
call struc_atoms_alloc( nfile, myid, nodes,  &
& alloc_mem, natom, nplwcs, kmax1cs, kmax2cs, kmax3cs,  &
& kmax1cs, kmax2cs, kmax3cs, kmax1d, kmax2d, kmax3d,  &
& lnlpp_g, llclpp_g, lpcc_g, lsphexp, pwscale )

!-----for the double-grid method
!if( ldouble_grid_recip ) then
!    call struc_atoms_dg_alloc( nfile, myid, nodes,  &
!& alloc_mem, natom, nion_nod, idouble_grid_method, lvacuum )
!end if


!-----for stress calcu by local pp.
call svvext_alloc( nfile, myid, nodes,  &
& alloc_mem, llclpp_r, llclpp_g, lstress, mshnod(1), pwscale )


if( lnlpp_g .or. lvand ) then

    !-----memory allocation for nonlocal pp
    call nl_variables_alloc( nfile, myid, nodes,  &
& alloc_mem, nplw, nplwex, nband, nbnod, nbncnt, nbndsp,  &
& ntype, nhk1, nhk2, natom, natnod1, natnod2, natnod, natcnt, natdsp,  &
& nhk1_nat, nhk2_nat, lvandi, lstress, lspin, natom_alloc, natnod_alloc )

    !-----memory allocation for ultrasoft pp
!    if( lvand )  &
!&       call uspp_variables_alloc( nfile, myid, nodes,  &
!& alloc_mem, ntype, natom, natnod1, natnod2, natnod, natcnt, natdsp, ioag,  &
!& nband, nbnod, nbncnt, nbndsp,  &
!& kfft0b, nplw3, kmax1b, kmax2b, kmax3b, gboxx, gboxy, gboxz,  &
!& lstress, lspin, pwscale, jhybrid, mxkpdupnm, lefield, lsawtooth,  &
!& natom_alloc, natnod_alloc )

    !-----memory allocation for the PAW method
!    if( lpaw )  &
!&       call nlpaw_variables_alloc( nfile, myid, nodes,  &
!& alloc_mem, ntype, natom, nspnmx, lrtddft )

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


!-----initialize the calculation in reciprocal space
if( lnlpp_g .or. llclpp_g .or. lpcc_g .or. lsphexp ) then
    !-----set sine, cosine function for atoms
    call stscfn_pw( nfile, myid, nodes,  &
& ycos, ysin, DCQ1, DCQ2, DCQ3, DSQ1, DSQ2, DSQ3,  &
& kmax1cs, kmax2cs, kmax3cs, natom, nplwcs, nplw5, nga, ngb, ngc,  &
& lclust, ratm, hcell, nd1v, nd1vks, lnlpp_g, llclpp_g, lpcc_g, lsphexp )
end if

!-----for the double-grid method
!if( ldouble_grid_recip ) then
!    call stscfn_pw_dg( nfile, myid, nodes,  &
!& ratm, iatoit, iatmpt_nod, natom, nion_nod )
!end if

!-----for uniform electric field for insulator
!if( lefield .and. .not.lsawtooth ) then
!    call stscfn_pw_efield( nfile, myid, nodes,  &
!& lclust, lvand, hcell, ntype, natom, natnod, iatoit, ioa, lvandi, lking,  &
!& ratm, wrka )
!end if


!-----initialize the calculation for ultrasoft pseudopotential
!if( lvand ) then
!    call inivan( nfile, myid, nodes,  &
!& nplw, nplw3, kfft0d, kfft1b, kfft2b, kfft3b, kfft0b,  &
!& gboxx, gboxy, gboxz, ngboxa, ngboxb, ngboxc,  &
!& mshglb, nd1vks(1), nd1vks(2), nd1vks(3),  &
!& lclust, hcell, nd1v, rvol, lstress,  &
!& ntype, nhk1, nhk2, natom, iatoit, ioa, ioag, ratm, lvandi, lking,  &
!& lvand_r, lvand_g, ycos, ysin, nplwcs,  &
!& lvacuum, lmax, rdel_c, rdelg_c, hci, lorthrhmbc,  &
!& mshglb_c, kfft1, kfft2, kfft3, ylmr, ylmi, mxl, alloc_mem, jhybrid )
!
!   kvecdo1: do kvec = 1, nknod
!      if( .not.lgamma ) then
!          call set_kvec( kvec )
!          call get_wfk( cgjr, 1 )
!          call get_plwk( nplw, nplwex, nspnod )
!          if( lspin ) then
!              call get_wfk( gdcrsv, 2 )
!              lcgjsv = .true.
!          end if
!      end if
!
!   do nspin = 1, nspnmx
!      if( lspin ) then
!          call stspud( nspin, cgjr, gdcrsv, nspnod, lcgjsv )
!      end if
!
!      dupif: if( nspin == 1 .or. nspin == 2 .and. .not.lduplicate ) then
!      if( lgamma ) then
!
!          if( lvand_g ) then
!              call calslm( nfile, myid, nodes,  &
!& cgjr, nband, nbnod1, nbnod2, nbnod, nplw, nplwex,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking, ycos, ysin, nplwcs )
!          end if
!          if( lvand_r ) then
!              call calslm_r( nfile, myid, nodes,  &
!& cgjr, nband, nbnod1, nbnod2, nbnod, nplw, nplwex,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking, iatoit,  &
!& mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, rvol, rdelv_c,  &
!& fft3x, fft3y, fftwork, ijkg, nplw2, eigr, eigi )
!          end if
!
!      else
!
!          if( lvand_g ) then
!              call set_ycosysin( nfile, myid, nodes,  &
!& natom, ycos, ysin, nplwcs )
!              call calslm_k( nfile, myid, nodes,  &
!& cgjr, nband, nbnod1, nbnod2, nbnod, nplw, nplwex,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking )
!          end if
!          if( lvand_r ) then
!              call set_eikr( nfile, myid, nodes,  &
!& mshglb_c, kfft1, kfft2, kfft3 )
!              call set_eikl( nfile, myid, nodes,  &
!& ratm, natom, .false., lvand_r )
!              call calslm_r_k( nfile, myid, nodes,  &
!& cgjr, nband, nbnod1, nbnod2, nbnod, nplw, nplwex,  &
!& ntype, nhk1, nhk2, natom, lvandi, lking, iatoit,  &
!& mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, rvol, rdelv_c,  &
!& fft3x, fft3y, fftwork, ijkg, nplw2, eigr, eigi )
!          end if
!
!      end if
!      end if dupif
!
!      if( lspin ) then
!          call svslmi( nfile, myid, nodes, nspin, lvandi, .true. )
!      end if
!   end do
!
!      if( .not.lgamma ) then
!          call svslmi_k_both( nfile, myid, nodes, lvandi, lspin, .true. )
!      end if
!   end do kvecdo1
!
!-----< beta | wave packets > in TDDFT
!call calslm_in_tddft( nfile, myid, nodes,  &
!nspnmx, lvand_g, lvand_r, nband, nbnod1, nbnod2, nbnod, nplw, nplwex,  &
!ntype, nhk1, nhk2, natom, lvandi, lking, ycos, ysin, nplwcs, iatoit, &
!mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, rvol, rdelv_c, &
!fft3x, fft3y, fftwork, ijkg, nplw2, eigr, apk )

!    !--- obtain S^(1/2) in TDDFT
!    call tddft_s_half( nfile, myid, nodes, &
! nplwex, nplw, ntype, nhk1, nhk2, natom, lvand_r, lvand_g, lvandi, lking, &
! rvol, ycos, ysin, nplwcs, iatoit, rdelv_c )
!
!end if

if( lnlpp_g .or. llclpp_g .or. lpcc_g .or. lvand ) then
    ct = timecnt()
    if(loutfile(1)) write(nfile(1),*)  &
& ' initialize g-space and/or USPP calc. : cpu-time :', ct - ct0
    if(loutfile(2)) write(nfile(2),*)  &
& ' initialize g-space and/or USPP calc. : cpu-time :', ct - ct0
    ct0 = ct
end if


!--- Gram-Schmidt orthonormalizaion ------------------------------------
call set_all2all_out( .false. )
kvecdo2: do kvec = 1, nknod
!   if( .not.lgamma ) then
!       call set_kvec( kvec )
!       call get_wfk( cgjr, 1 )
!       call get_plwk( nplw, nplwex, nspnod )
!       call get_indexk(  &
!& npnod1, npnod2, npnod, nplcnt, npldsp, nodes )
!       if( lspin ) then
!           call get_wfk( gdcrsv, 2 )
!           lcgjsv = .true.
!       end if
!       if( lvand ) call ldslmi_k_both( nfile, myid, nodes, lvandi, lspin, .true. )
!   end if

do nspin = 1, nspnmx
   dupif2: if( nspin == 1 .or. nspin == 2 .and. .not.lduplicate ) then
   if( lspin ) then
       if( nspin.eq.1 )  &
&          call stspud( nspin, cgjr, gdcrsv, nspnod, lcgjsv )
!       if( lvand ) call ldslmi( nfile, myid, nodes, nspin, lvandi, .true. )
   end if

   !---set nspin in tddft-fssh.f90
   call set_nspin_in_tddft_fssh( nspin )

   leig_start = kvec == 1     .and. nspin == 1
   leig_end   = kvec == nknod .and. nspin == nspnmx
   leig_end   = leig_end .or. kvec == nknod .and. nspin == 1 .and. lduplicate
   call set_all2all_out( leig_end )
   ctt0 = timecnt()

   call bdtogd( nfile, myid, nodes, ctt0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iodg, idstnd, .false., .true. )

!   if( lvand ) call slm_bdtogd( nfile, myid, nodes, ctt0,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& natom, natnod1, natnod2, natnod, natcnt, natdsp,  &
!& iodg, ioag, idstnd, .false., .true. )

!   if( lgamma .or. lgammak(kvec+nknod1-1) ) then
       call schmidt( nfile, myid, nodes, ct0,  &
& gdcr, iod, iodg, .true., .false., leig_start, leig_end )
!   else
!       call schmidt_k( nfile, myid, nodes, ct0,  &
!& gdcr, iod, iodg, .true., .false., leig_start, leig_end )
!   end if
   end if dupif2

   if( lspin ) then
       call stpsud( nspin, gdcr, gdcrsv, nspnod, lcgjsv )
!       if( lvand ) call svslmi( nfile, myid, nodes, nspin, lvandi, .false. )
       if( nspin.eq.1 ) then
           if( .not.lduplicate ) then
               !--- copy gdcr to cgjr ---
               call cprhcg( nfile, myid, nodes, cgjr, gdcr, nplwex, nbnod )
           else
               !--- copy gdcrsv to gdcr ---
               call dcopy_a_to_b( gdcrsv, gdcr, nspnod )
           end if
       end if
   end if
end do

!   if( .not.lgamma ) then
!       call set_wfk( gdcr, 1 )
!       if( lspin ) then
!           call set_wfk( gdcrsv, 2 )
!           lcgjsv = .true.
!       end if
!       if( lvand ) call svslmi_k_both( nfile, myid, nodes, lvandi, lspin, .false. )
!   end if
end do kvecdo2
!-----------------------------------------------------------------------

!      !==== stop here =========================
!      call fstop( nfile, myid, nodes, 'pwlda_set' )
!      !========================================


!-----set pratm
call set_pratm( nfile, myid, nodes, ratm, pratm, natom )


!-----set time step for MD
if( ifmd.eq.0 ) then
    nstep = max( nstepCG, nstepMD )
else if( ifmd.eq.1 ) then
    if( .not.lstart_ .or. lstart_ .and. nstepCG == 0 ) then
        nstepCG = nstep_ini
      else
        nstep_ini = 0
    end if
    nstep = nstepCG
else
    if( .not.lstart_ .or. lstart_ .and. nstepMD == 0 ) then
        nstepMD = nstep_ini
      else
        nstep_ini = 0
    end if
    nstep = nstepMD
end if


return
end




subroutine check_param_in_qm( nfile, myid, nodes,  &
& nstep_, nstop_, ifmd_, ierror  )
!-----------------------------------------------------------------------
!     check basic parameters for MD calculation
!-----------------------------------------------------------------------
use constants
use param
use param_atom
use pwlda_pp
use pwlda_atom
use pwlda_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: nstep_, nstop_, ifmd_
integer :: ierror

!-----error trap
if( nstep_ /= nstep ) ierror = 1
!      if( nstop_ /= nstop ) ierror = ierror + 10
nstop = nstop_
if( ifmd_  /= ifmd  ) ierror = ierror + 100
if( ierror /= 0 ) then
    write(*,*) 'myid=',myid,' error in check_param_in_qm',  &
&   nstep_, nstop_, ifmd_, nstep, nstop, ifmd
end if

nstop = nstep + nstop

return
end subroutine




subroutine set_nstep_in_qm( nfile, myid, nodes, nstep_ )
!-----------------------------------------------------------------------
!     set time step
!-----------------------------------------------------------------------
use pwlda_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: nstep_

nstep = nstep_

return
end subroutine




subroutine pwlda( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
use param
implicit none
integer :: nfile(*), myid, nodes
integer :: ierror


if( .not.ltddft ) then

    !--- LDA or GGA calculation with plane-wave method
    call pwlda_md( nfile, myid, nodes, ierror )

  else

    !--- TDDFT
!    call pwlda_tddft( nfile, myid, nodes, ierror )

end if


return
end subroutine




subroutine pwlda_md( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!     LDA or GGA calculation with plane-wave method
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
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes
integer :: ierror

!-----declare local variables
logical :: lreset, lresetud, lretrn, ltimecnt, lortho, lcstress, lhcexe
logical :: leig_start, leig_end
character(30) :: fname1
real*8, dimension(10) :: tc
real*8, dimension(10) :: unwork, unbuff
logical :: lonsite
logical :: lmixhold   ! .true. = hold matrix in charge mixing
logical :: lchkzansa
logical :: lvdw_pre = .true.
logical :: laspc_exec, laspc_conv, laspc_rot
logical :: lconverge, lfssh_updtcg
logical :: ltransition
real*8,  allocatable, dimension(:) :: prevzansa2
real*8  :: paw_mix0, plusU_mix0
logical :: lcholesky
logical :: lexcitedstates
logical :: lhfstart, lhfend
logical :: lncinit
integer :: status
character(100) :: fmtHF, fmtKS
character(200) :: fmtprts
save lvdw_pre


ct0 = timecnt()
ct00 = ct0

!---allocate array for convergence stabilizer
if( .not.allocated(prevzansa2) ) then
    allocate( prevzansa2(0:iscfmx+1), stat=status )
end if

ltimecnt = ifmd == 0 .or. loutlog
!      ltimecnt = laspc
if( ltimecnt ) then
    istdot = 1
  else
!    istdot = 2  !      simplified output form in eigen
!    istdot = 3  ! more simplified output
    istdot = 4   !              no output
end if

ltransition = .false.
fsshcycle: do

!lonsite = .true.
!lonsite = .not.lstart .and. nstep == nstep_ini .or. lstart .and. ifmd == 0
lonsite = nstep == nstep_ini .or. lstart .and. ifmd == 0
  paw_mix0 =   paw_mix
plusU_mix0 = plusU_mix


hfmix = 1.d0
!if( jhybrid /= 0 .and. .not.lgamma ) then
!    hfrhcr(:,:) = 0.d0
!    do nspin = 1, nspnmx
!    do kvec = 1, nknod
!       call set_kvec( kvec )
!       call set_hfrhcr( hfrhcr(1,nspin), nspin )
!    end do
!    end do
!end if

if( lnoncollinear ) then
    lncinit = wvratio(1,1,1) + wvratio(1,2,1) < 0.1d0
else
    lncinit = .false.
end if

ctmd0  = ct0
!-----------------------------------------------------------------------
if( nstep > nstep_ini .or. lstart ) then
!---    itrial = 0
    lvdw_pre = .false.
end if
lreset   = .true.
lresetud = .true.
zansa2min = 1.d99
do ii = 1, 2
if( loutfile(ii) ) then
    write(nfile(ii),*) ' '
    write(nfile(ii),*) ' '
    write(nfile(ii),'(30a,a,i12,a,30a)') ('=',i=1,30),' step:',nstep,' ',  &
&                      ('=',i=1,30)
    write(nfile(ii),*) ' '
end if
end do


!-----ASPC method
laspc_exec = laspc .and. keypwv >= nprvmx .and. keypcd >= ncprvmx
laspc_exec = laspc_exec .and. nspnmx2 == 1  &
&       .or. laspc_exec .and. nspnmx2 >= 2 .and. keypmg >= ncprvmx
laspc_exec = laspc_exec .and. nstep > (nstep/max(nskp_aspc,1))*nskp_aspc + nprvmx - 1
laspc_conv = .false.
laspc_rot  = .false. .and. laspc_exec

!-----reset plane-wave method
call pwlda_reset( nfile, myid, nodes, ltimecnt, ierror )
if( ierror /= 0 ) return


lcstress = lstress .and.  &
&        ( ifmd.eq.0 .or. mod(nstep,nskip_stress).eq.0 )

lmixhold = lstart .and. ifmd >= 1 .and. nstep >= nstep_ini+1
if( lmixhold ) then
    aslh = aslh0
!    if( myid == 0 ) then
!        write(nfile(1),'(a,f7.4)') ' *** lmixhold = .true. / aslh =', aslh
!        write(nfile(2),'(a,f7.4)') ' *** lmixhold = .true. / aslh =', aslh
!    end if
  else
    aslh = aslh0    !* 0.5d0
    if(loutfile(1)) write(nfile(1),'(a,f7.4)') ' *** lmixhold = .false. / aslh =', aslh
    if(loutfile(2)) write(nfile(2),'(a,f7.4)') ' *** lmixhold = .false. / aslh =', aslh
end if
lstart   = .true.

lfssh_updtcg = .false.

if( ifmd >= 1 .and. nstep >= nstep_ini+1 ) then
!--- update atomic configuration ---------------------------------------
    call updt_atom_disp( nfile, myid, nodes,  &
& lclust, ratm, pratm, prevr, hcell,  &
& ntype, natom, nhk1, nhk2 )

    !-----set pratm
    call set_pratm( nfile, myid, nodes, ratm, pratm, natom )


    !-----initialize the calculation in reciprocal space
    if( lnlpp_g .or. llclpp_g .or. lpcc_g .or. lsphexp ) then
        !-----set sine, cosine function for atoms
        call stscfn_pw( nfile, myid, nodes,  &
& ycos, ysin, DCQ1, DCQ2, DCQ3, DSQ1, DSQ2, DSQ3,  &
& kmax1cs, kmax2cs, kmax3cs, natom, nplwcs, nplw5, nga, ngb, ngc,  &
& lclust, ratm, hcell, nd1v, nd1vks, lnlpp_g, llclpp_g, lpcc_g, lsphexp )
    end if

!    !-----for the double-grid method
!    if( ldouble_grid_recip ) then
!        call stscfn_pw_dg( nfile, myid, nodes,  &
!& ratm, iatoit, iatmpt_nod, natom, nion_nod )
!    end if

    !-----for uniform electric field for insulator
!    if( lefield .and. .not.lsawtooth ) then
!        call stscfn_pw_efield( nfile, myid, nodes,  &
!& lclust, lvand, hcell, ntype, natom, natnod, iatoit, ioa, lvandi, lking,  &
!& ratm, wrka )
!    end if


    !-----initialize the calculation for ultrasoft pseudopotential
!    if( lvand ) then
!        call inivan( nfile, myid, nodes,  &
!& nplw, nplw3, kfft0d, kfft1b, kfft2b, kfft3b, kfft0b,  &
!& gboxx, gboxy, gboxz, ngboxa, ngboxb, ngboxc,  &
!& mshglb, nd1vks(1), nd1vks(2), nd1vks(3),  &
!& lclust, hcell, nd1v, rvol, lstress,  &
!& ntype, nhk1, nhk2, natom, iatoit, ioa, ioag, ratm, lvandi, lking,  &
!& lvand_r, lvand_g, ycos, ysin, nplwcs,  &
!& lvacuum, lmax, rdel_c, rdelg_c, hci, lorthrhmbc,  &
!& mshglb_c, kfft1, kfft2, kfft3, ylmr, ylmi, mxl, alloc_mem, jhybrid )
!    end if

    !--- set <phi(t-dt)|phi(t)> matrix for TDDFT-FSSH
    call tddft_fssh_set_matrix( nfile, myid, nodes,  &
& gdcr, nplw, nplwex, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& lspin, gdcrsv, nspnod, lcgjsv, lvand, lvandi, ntype, rvol,  &
& dmtrxr, nbxxxx, prod, prodr, bijr )

    !--- get coefficients for extrapolation
    call setacc( prevr, facc1, facc2, facc3,  &
&                ntype, natom, nhk1, nhk2 )

        ct = timecnt()
        ct0 = ct
        if(loutfile(1)) write(nfile(1),*) ' Input charge estimation'
        if(loutfile(2)) write(nfile(2),*) ' Input charge estimation'

    !--- estimation of input charge
    call chgest( nfile, myid, nodes,  &
& rhgr, rinr, nplw5ex, nplw5, prvrho, ncprvmx, nspnmx2, thrhgr,  &
& ichest, keypcd, keypmg, nstep, nstep_ini, facc1, facc2, facc3,  &
& lspin, rhoud, diffud, lfixud,  &
& glocal, rho, mshnod, rdelv, nel,  &
& mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d,  &
& fft3x, fft3y, fftwork, ijkgd, nplw5, laspc_exec )

        ct = timecnt()
        if(loutfile(1)) write(nfile(1),'(1x,a20,f10.4)')  &
&                     ' cpu-time ( sec )  :', ct - ct0
        if(loutfile(2)) write(nfile(2),'(1x,a20,f10.4)')  &
&                     ' cpu-time ( sec )  :', ct - ct0
        ct0 = ct

        if(loutfile(1)) write(nfile(1),*) ' Subspace alignment'
        if(loutfile(2)) write(nfile(2),*) ' Subspace alignment'
        t_comm = 0.d0

    !--- subspace alignment
!    if( lgamma ) then
        call sbalin( nfile, myid, nodes,  &
& facc1, facc2, facc3, laspc_exec, t_comm, ltimecnt )
!      else
!        call sbalin_k( nfile, myid, nodes,  &
!& facc1, facc2, facc3, t_comm, ltimecnt )
!    end if

    ct = timecnt()
    if(loutfile(1)) write(nfile(1),'(1x,a20,f10.4)')  &
&                     ' cpu-time ( sec )  :', ct - ct0
    if(loutfile(2)) write(nfile(2),'(1x,a20,f10.4)')  &
&                     ' cpu-time ( sec )  :', ct - ct0
    ct0 = ct

    !---TDDFT-FSSH: main routine to determine new occupancies
!    if( ltddft_fssh ) then
        if( lspin ) then
            call tddft_fssh( nfile, myid, nodes,  &
& nband, nbnod, nspnmx, egvlud, nsordr, wegud, dtmd, &
& bijr, bm1r, bx0r, prod, prodr, prods, betk, bzan, w1r, w1i, ew1, nstep, &
& keypwv, keypcd, laspc_exec, ltransition, lexcitedstates )
          else
            call tddft_fssh( nfile, myid, nodes,  &
& nband, nbnod, nspnmx, eig, norder, occ, dtmd, &
& bijr, bm1r, bx0r, prod, prodr, prods, betk, bzan, w1r, w1i, ew1, nstep, &
& keypwv, keypcd, laspc_exec, ltransition, lexcitedstates )
        end if
        if( .not.lexcitedstates ) then
            !---Ground state!!! stop TDDFT-FSSH
            if( ltddft_fssh ) then
                if(loutfile(1)) write(nfile(1),*) '***** Ground state!!! stop TDDFT-FSSH.'
                if(loutfile(2)) write(nfile(2),*) '***** Ground state!!! stop TDDFT-FSSH.'
                ltddft_fssh = .false.
                lfssh_gsscf = .false.
                call set_ltddft_fssh( ltddft_fssh, lfssh_gsscf )
            end if
        end if
!    end if
!-----------------------------------------------------------------------
end if
!=======================================================================
!        minimization problem starts with given atomic positions.
!=======================================================================

if( lmixhold ) then
    lmixhold = keypwv >= nprvmx .and. keypcd >= ncprvmx
    if( .not.lmixhold ) then
        aslh = aslh0    !* 0.5d0
        if(loutfile(1)) write(nfile(1),'(a,f7.4)') ' *** lmixhold = .false. / aslh =', aslh
        if(loutfile(2)) write(nfile(2),'(a,f7.4)') ' *** lmixhold = .false. / aslh =', aslh
    end if
end if

!-----------------------------------------------------------------------
!    partial core charge for exchange-correlation energy
rhocore(1:mshnod(1)) = 0.d0
lpccif: if( lpcc ) then

    !-----in real space
    if( lpcc_r ) then
        call pccchg( nfile, myid, nodes,  &
& rhocore, lclust, rdel, rdelg, hcell, hci, lorthrhmbc, lvacuum,  &
& ntype, nhk1, nhk2, ratm, ltimecnt, lpcci, rpcc, lpking,  &
& mshxyz(mulpit(1)),  &
& mulx1(1), mulx2(1), muly1(1), muly2(1), mulz1(1), mulz2(1),  &
& mshnx(mulpms(1)), mshny(mulpms(1)), mshnz(mulpms(1)), mshnod(1),  &
& mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1) )
    end if

    !-----in reciprocal space
    if( lpcc_g ) then
        call pccchg_g( nfile, myid, nodes, ltimecnt,  &
& rhocore, ntype, nhk, nhk1_nod, nhk2_nod, iatoit, iatmpt_nod, natom,  &
& nion_nod, lpcci, lpking, nplw5, nplw5ex, nplw,  &
& nd1vks(1), nd1vks(2), nd1vks(3), kfft0d,  &
& nga, ngb, ngc, ijkgd, thrhgr, rvol,  &
& mshnod(1), mftnod, mftdsp, mfd2ft, ntotfd,  &
& fft3x, fft3y, fftwork, lstress,  & !& ycos, ysin, nplwcs,
& DCQ1, DCQ2, DCQ3, DSQ1, DSQ2, DSQ3, kmax1cs, kmax2cs, kmax3cs )
    end if

    !check
    tsum = 0.d0
    do i = 1, mshnod(1)
       tsum = tsum + rhocore(i)
    end do
    call gdsum(tsum,1,dbuf1r)
    if( ltimecnt ) then
        ct = timecnt()
        if(loutfile(1)) write(nfile(1),*) ' total core charge:', tsum * rdelv
        if(loutfile(2)) write(nfile(2),*) ' total core charge:', tsum * rdelv
        ct0 = ct
    end if
end if lpccif


if( ltimecnt ) then
    if(loutfile(1)) write(nfile(1),*) ' Local pseudo potential'
    if(loutfile(2)) write(nfile(2),*) ' Local pseudo potential'
end if
!-----------------------------------------------------------------------
!    local pseudopotential : vext
do i = 1, mshnod(1)
   vext(i) = 0.d0
end do
!--- set sine, cosine function for atoms
if( .not.lclust ) then
    call stscfn(  &
& ntype, ratm, iatoit, iatmpt_nod, natom, nion_nod,  &
& ldouble_grid_recip, nd1v, nd1vks )
end if

!-----in real space
if( llclpp_r ) then
    call vlocal( nfile, myid, nodes,  &
& lclust, lorthrhmbc, nd1v, ntype, nhk1, nhk2, nhk1_nod, nhk2_nod,  &
& zv, nel, ratm, iatoit, iatmpt_nod, nion_nod, natom,  &
& llking, mx1loc,  &
& hcell, rdel, lhfull, vext, tablc, tablca, dltlc, rmxlc,  &
& mshxyz(mulpit(1)),  &
& mulx1(1), mulx2(1), muly1(1), muly2(1), mulz1(1), mulz2(1),  &
& mshnx(mulpms(1)), mshny(mulpms(1)), mshnz(mulpms(1)), mshnod(1),  &
& mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1),  &
& tmpk, ltimecnt )
end if

!-----in reciprocal space
if( llclpp_g ) then
    call vlocal_g( nfile, myid, nodes,  &
& vext, ntype, nhk, nhk1_nod, nhk2_nod, iatoit, iatmpt_nod, natom,  &
& nion_nod, llking, nplw5, nplw5ex, nplw,  &
& nd1vks(1), nd1vks(2), nd1vks(3), kfft0d,  &
& nga, ngb, ngc, ijkgd, thrhgr, rvol,  &
& mshnod(1), mftnod, mftdsp, mfd2ft, ntotfd,  &
& fft3x, fft3y, fftwork, lstress,  & !& ycos, ysin, nplwcs,
& DCQ1, DCQ2, DCQ3, DSQ1, DSQ2, DSQ3, kmax1cs, kmax2cs, kmax3cs )
end if

!-----for the double-grid method
!if( ldouble_grid_recip ) then
!    call vlocal_dg( nfile, myid, nodes,  &
!& vext, zv, ntype, nhk, nhk1_nod, nhk2_nod, iatoit, iatmpt_nod,  &
!& natom, nion_nod,  &
!& mshnod(1), mftnod, mftdsp, mfd2ft, ntotfd, fft3x, fft3y, fftwork,  &
!& hcell, glocal, kfft0d, ltimecnt )
!end if

ct = timecnt()
if( ltimecnt ) then
    if(loutfile(1)) write(nfile(1),*) ' set local pp.   : cpu-time :', ct-ct0
    if(loutfile(2)) write(nfile(2),*) ' set local pp.   : cpu-time :', ct-ct0
end if
ct0 = ct

!      if( lshdep ) then
!!---  shape denpendent terms
!          do i = 1, mshnod(1)
!             vlshdp(i) = 0.d0
!          end do
!          call shvloc( nfile, myid, nodes,  &
!     & lvacuum, cshdep, nel, ntype, natom, nhk1, nhk2, zv, ratm,
!     & hcell, rdel, vlshdp, dpion, dpion2,
!     & mshxyz(mulpit(1)),
!     & mulx1(1), mulx2(1), muly1(1), muly2(1), mulz1(1), mulz2(1),
!     & mshnx(mulpms(1)), mshny(mulpms(1)), mshnz(mulpms(1)), mshnod(1),
!     & mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1),
!     & ltimecnt )
!          do i = 1, mshnod(1)
!             vext(i) = vext(i) + vlshdp(i)
!          end do
!      end if
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!    set parameters for calculation of norm-conserving pseudopotentials

!-----initialize real-space calculation
if( lkbpp_r ) then
    call setnlc( nfile, myid, nodes,  &
& lvacuum, lclust, ntype, nhk1, nhk2, lmax,lclno,lchk,lkbppi,lking,  &
& iatoit, rdel_c, rdelg_c, hcell, hci, lorthrhmbc, ratm,  &
& mshglb_c, kfft1, kfft2, kfft3,  &
& ylmr, ylmi, natom, mxl, nylmmx, mx1, alloc_mem )
end if

!-----initialize reciprocal-space calculation
if( lkbpp_g ) then
    call setnlc_g( nfile, myid, nodes,  &
& nplw, ntype, nhk1, nhk2, natom, lkbppi, lking,  &
& ycos, ysin, nplwcs, rvol )
end if

if( lkbpp ) then
   ct = timecnt()
   if( ltimecnt ) then
       if(loutfile(1)) write(nfile(1),*) ' set nonlocal pp.   : cpu-time :',ct-ct0
       if(loutfile(2)) write(nfile(2),*) ' set nonlocal pp.   : cpu-time :',ct-ct0
   end if
   ct0 = ct
end if
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!    set meshes around each atom to calculate noncollinear magnetism
!if( lnoncollinear ) then
!    call ncmagne_calweight( nfile, myid, nodes,  &
!& lclust, rdel, rdelg, hcell, hci, lorthrhmbc, lvacuum, nd1v, nd1vks, rdelv,  &
!& ntype, nhk1, nhk2, iatoit, ratm, dbuf,  &
!& mshxyz(mulpit(1)),  &
!& mulx1(1), mulx2(1), muly1(1), muly2(1), mulz1(1), mulz2(1),  &
!& mshnx(mulpms(1)), mshny(mulpms(1)), mshnz(mulpms(1)), mshnod(1),  &
!& mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1) )
!end if
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! well potential
!if( lwell ) then
!    call set_well_potential( nfile, myid, nodes,  &
!& vwell, lclust, rdel, rdelg, hcell, lorthrhmbc, nd1v, nd1vks, rdelv,  &
!& ntype, nhk1, nhk2, iatoit, ratm,  &
!& mshxyz(mulpit(1)),  &
!& mulx1(1), mulx2(1), muly1(1), muly2(1), mulz1(1), mulz2(1),  &
!& mshnx(mulpms(1)), mshny(mulpms(1)), mshnz(mulpms(1)), mshnod(1),  &
!& mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1) )
!
!    vext(1:mshnod(1)) = vext(1:mshnod(1)) + vwell(1:mshnod(1))
!end if


!-----------------------------------------------------------------------
! sawtooth electric field
!if( lefield ) then
!    if( lsawtooth ) then
!        call set_efield_potential( nfile, myid, nodes,  &
!& vefield, efield, lsawtooth_shape, lsawtooth_xyz,  &
!& lclust, rdel, rdelg, hcell, lorthrhmbc, nd1v, nd1vks, rdelv,  &
!& mshxyz(mulpit(1)),  &
!& mulx1(1), mulx2(1), muly1(1), muly2(1), mulz1(1), mulz2(1),  &
!& mshnx(mulpms(1)), mshny(mulpms(1)), mshnz(mulpms(1)), mshnod(1),  &
!& mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1) )
!
!        vext(1:mshnod(1)) = vext(1:mshnod(1)) + vefield(1:mshnod(1))
!    end if
!
!    !---ionic polarization
!    call cal_polion( nfile, myid, nodes, &
!& lclust, ntype, nhk1, nhk2, ratm, iatoit, natom, zv, hcell,  &
!& ltimecnt, lefield_islts, lsawtooth_shape )
!end if
!-----------------------------------------------------------------------

!=======================================================================
!                    outer iteration starts.
!=======================================================================
jcycl = 0
jstab = 0
eneold = 1.d+10
outer: do
ct = timecnt()
ct0 = ct

if( ltimecnt ) then
    if(loutfile(1)) write(nfile(1),*) ' Hartree potential'
    if(loutfile(2)) write(nfile(2),*) ' Hartree potential'
end if
!-----------------------------------------------------------------------
!      Hartree potential : vhar
!-----------------------------------------------------------------------

t_comm = 0.d0
!hartif: if( .not.lclust .and. .not.ldouble_grid ) then !==========

!-----calculation in g-space
    call harting( nfile, myid, nodes,  &
& glocal, rhgr, nplw5ex, nplw5, recnrmex,  &
& mfd2ft, ntotfd, nd1vks, kfft0d,  &
& fft3x, fft3y, fftwork, ijkgd, nplw5 )
    !-----store local hartree potential
    call distlc( nfile, myid, nodes,  &
& vhar, mshnod(1), glocal, mftdsp )

    !-----for the double-grid method
!    if( ldouble_grid_recip ) then
!
!        if( ltimecnt ) then
!            ehars = 0.d0
!            do i = 1, mshnod(1)
!               ehars = ehars + vhar(i)*rho(i)
!            end do
!            call gdsum(ehars,1,dbuf1r)
!            ehars = 0.5d0 * ehars * rdelv
!          if(loutfile(1)) write(nfile(1),*) '       short-range (r-space) :',ehars
!          if(loutfile(2)) write(nfile(2),*) '       short-range (r-space) :',ehars
!        end if
!
!        call harting_dg( nfile, myid, nodes,  &
!& vhar, rhgr, nplw5, ijkgd, nel,  &
!& mshnod(1), mftnod, mftdsp, mfd2ft, ntotfd, fft3x, fft3y, fftwork,  &
!& hcell, glocal, kfft0d, ltimecnt )
!
!        if( ltimecnt ) then
!            eharl = 0.d0
!            do i = 1, mshnod(1)
!               eharl = eharl + vhar(i)*rho(i)
!            end do
!            call gdsum(eharl,1,dbuf1r)
!            eharl = 0.5d0 * eharl * rdelv
!            if(loutfile(1)) write(nfile(1),*) '       long -range (r-space) :',eharl-ehars
!            if(loutfile(2)) write(nfile(2),*) '       long -range (r-space) :',eharl-ehars
!        end if
!
!    end if


    ct = timecnt()
    if( ltimecnt ) then
        if(loutfile(1)) write(nfile(1),*) ' Hartree potential : cpu-time :',ct-ct0
        if(loutfile(2)) write(nfile(2),*) ' Hartree potential : cpu-time :',ct-ct0
    end if
    ct0 = ct

!else hartif !=====================================================

!!-----calculation in r-space
!iter = 1000
!
!!--- tolerance for Poisson eq. [/el] -> [/Nel]
!ttolcg = tolcg*dble(nel)/rdelv
!
!
!!-----Poisson equation
!!dgridif: if( .not.ldouble_grid ) then
!
!   !-----set boundary condition
!   clustif: if( lclust ) then
!
!      !-----for cluster calculation
!      !-----boundary condition by multipole expansions
!      npnd0 = 3*(nd2v-1)*(nd2v+1) + 1
!      call cormul( nfile, myid, nodes,  &
!& lsphere, vk, rho, cdv2(npnd0,1), nel, nd2v, rdel, rdelv, cofmas,  &
!& rmax, muldatx,  noddatx, nd1vks,  &
!& mshxyz(mulpit(1)),  &
!& mulx1(1), mulx2(1), muly1(1), muly2(1), mulz1(1), mulz2(1),  &
!& mshnx(mulpms(1)), mshny(mulpms(1)), mshnz(mulpms(1)), mshnod(1),  &
!& mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1),  &
!& tmpk, tmpl, tmpm, tmpn, xx, xxx )
!
!      do i = 1, mshnod(1)
!         vk(i) = 8.d0*pi*rho(i) - vk(i)
!      end do
!
!      if( ltimecnt ) then
!         ct = timecnt()
!         if(loutfile(1)) write(nfile(1),*) '                 set boundary cond. ',  &
!&                         ' : cpu-time :', ct - ct0
!         if(loutfile(2)) write(nfile(2),*) '                 set boundary cond. ',  &
!&                         ' : cpu-time :', ct - ct0
!      end if
!
!   else clustif
!
!      !-----for bulk calculation
!      avchg8 = 8.d0*pi*dble(nel)/volume
!      do i = 1, mshnod(1)
!         vk(i) = 8.d0*pi*rho(i) - avchg8
!      end do
!
!   end if clustif
!
!
!   !---- solve Poisson eq. by multigrid method
!   call mcgrad( nfile, myid, nodes, nn, myparity, npx, npy, npz,  &
!& vk, vhar, iter, rk, x, xx, dbuf, dbufr, tmpk,  &
!& hdiag, tmpn, residal, ttolcg, weigrd, nhit1, multg, muls,  &
!& noddatx,  &
!& muldatx, mulext, nodexe, ltimecnt, .false., .true., lretrn,  &
!& t_comm,  &
!& nmm, mshnod, mshnew, mshxyz, mulpit, mulx1, mulx2, muly1, muly2,  &
!& mulz1, mulz2, mshnx, mshny, mshnz, mshx1, mshy1, mshz1,  &
!& mshx, mshy, mshz, mulpms, ismx, ismy, ismz,  &
!& ismshx, ismshy, ismshz, mulpsx, mulpsy, mulpsz,  &
!& mulyzs, mulzxs, mulxys, irmx, irmy, irmz, irmshx, irmshy, irmshz,  &
!& mulprx, mulpry, mulprz, mulyzr, mulzxr, mulxyr,  &
!& ismfpx, ismfpy, ismfpz, irmfpx, irmfpy, irmfpz, ncomct, mulnd2,  &
!& amcdv2, multg, mulexf, muldatx, mulnpx, mulnpy, mulnpz,  &
!& mulyzx, mulzxx, mulxyx, ndv2x3,  &
!& fmshx, msndx4, msndy4, msndz4, .true., xxx, nrstrct, prevres )
!
! norder_pois = nmm


!else dgridif


   !---- solve Poisson eq. by multigrid & double-grid method
!   call mcgrpo( nfile, myid, nodes, nn, myparity, npx, npy, npz,  &
!& vhar, rho, lclust, lsphere, nel, volume, nd1v, nd2v, cdv2,  &
!& rdel, rdelv, cofmas, disp_org, rmax,  &
!& lmuls, norder_pois, iter, residal, ttolcg, weigrd, nhit1,  &
!& nultg, ltimecnt, .false., .true., lretrn, t_comm, ct0,  &
!& mshnod(1), mshxyz(mulpit(1)),  &
!& mulx1(1), mulx2(1), muly1(1), muly2(1), mulz1(1), mulz2(1),  &
!& mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1) )

!end if dgridif


!ct = timecnt()
!if( ltimecnt ) then
!    do i = 1, 2
!    if( loutfile(i) ) then
!        write(nfile(i),*) '                              total ',  &
!&                         ' : cpu-time :', ct - ct0
!        write(nfile(i),*) '                                    ',  &
!&                         ' : com-time :', t_comm
!
!!              write(nfile(i),*) ' Poisson eq. has been solved',
!!     &                          ' by conjugate gradient method.'
!        write(nfile(i),*) ' Poisson eq. has been solved',  &
!&                         ' by multigrid method.'
!        write(nfile(i),*) '  order of matrix  :', norder_pois
!        write(nfile(i),*) '  No. of iteration :', iter
!        write(nfile(i),*) '  residual         :', residal
!!           else
!!              write(nfile(i),*) ' Poisson eq. (multigrid method)',
!!     &                      ' : cpu-time :', ct - ct0
!!              write(nfile(i),*) '                               ',
!!     &                      ' : com-time :', t_comm
!!              write(nfile(i),*) '  residual         :', residal
!   end if
!end do
!end if
!ct0 = ct

!end if hartif !===================================================

if( lfssh_updtcg .and. ltddft_nscforce ) then
    call dexchange( vhar, vhar_out, mshnod(1) )
    delrho(1:mshnod(1)) = rho(1:mshnod(1)) - delrho(1:mshnod(1))
    !---xc contribution to NSC force
    call xc_nscforce( nfile, myid, nodes )
    if( ltimecnt ) then
!        call qm_gsync
        ct = timecnt()
        if(loutfile(1)) write(nfile(1),*) ' xc_nscforce : cpu-time :', ct-ct0
        if(loutfile(2)) write(nfile(2),*) ' xc_nscforce : cpu-time :', ct-ct0
        ct0 = ct
    end if
    exit outer
end if


do i = 1, mshnod(1)
   vhshdp(i) = 0.d0
end do
!      if( lshdep ) then
!!--- shape dependent terms
!          call shdhar( nfile, myid, nodes,  &
!     & lvacuum, cshdep, nel, rdel, rdelv, rho, vhshdp, dprho, dprho2,
!     & mshxyz(mulpit(1)),
!     & mulx1(1), mulx2(1), muly1(1), muly2(1), mulz1(1), mulz2(1),
!     & mshnx(mulpms(1)), mshny(mulpms(1)), mshnz(mulpms(1)), mshnod(1),
!     & mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1),
!     & ltimecnt )
!!          do i = 1, mshnod(1)
!!             vhar(i) = vhar(i) + vhshdp(i)
!!          end do
!!---
!          ehshdp = 0.d0
!          elshdp = 0.d0
!          do i = 1, mshnod(1)
!             ehshdp = ehshdp + vhshdp(i)*rho(i)
!             elshdp = elshdp + vlshdp(i)*rho(i)
!          end do
!          dbuf(1) = ehshdp
!          dbuf(2) = elshdp
!          call gdsum(dbuf,2,dbufr)
!          ehshdp = 0.5d0 * dbuf(1) * rdelv
!          elshdp =         dbuf(2) * rdelv
!          ct = timecnt()
!          if( ltimecnt .and. myid.eq.0 ) then
!              do ii = 1, 2
!                 write(nfile(ii),*) ' * Dipole moment'
!                 write(nfile(ii),'((6x,a10,3e19.10))')
!     &         'atomic    ',(dpion(i), i = 1, 3),
!     &         'electronic',(-dprho(i), i = 1, 3),
!     &         'total     ',(dpion(i)-dprho(i), i = 1, 3)
!              write(nfile(ii),*) '              shape dependent terms ',
!     &                      ' : cpu-time :', ct - ct0
!              end do
!          end if
!          ct0 = ct
!      end if


!--  Hartree energy : ehar
ehar = 0.d0
do i = 1, mshnod(1)
   ehar = ehar + ( vhar(i) + vhshdp(i) )*rho(i)
end do
call gdsum(ehar,1,dbuf1r)
ehar = 0.5d0 * ehar * rdelv
!-----------------------------------------------------------------------

!check ---------------------------------------
!      do i = 1, mshnod(1)
!         x(i) = vhar(i)
!      end do
!      mul = 1
!      mul2x = 2*(npx+1)*(mul-1) + 1
!      mul2y = 2*(npy+1)*(mul-1) + 1
!      mul2z = 2*(npz+1)*(mul-1) + 1
!      call bdcopy( x, dbuf, dbufr, mulext, noddatx, t_comm,
!     &    mshnod(mul), ismx(mul2x), ismy(mul2y), ismz(mul2z),
!     &    ismshx(mulpsx(mul)), ismshy(mulpsy(mul)), ismshz(mulpsz(mul)),
!     &    mulyzs(mul), mulzxs(mul), mulxys(mul),
!     &    irmx(mul2x), irmy(mul2y), irmz(mul2z),
!     &    irmshx(mulprx(mul)), irmshy(mulpry(mul)), irmshz(mulprz(mul)),
!     &    mulyzr(mul), mulzxr(mul), mulxyr(mul), ncomct(1,mul),
!     &    npx, npy, npz, nn, myparity )
!c  --- matrix by vector operation ---
!      npnd0 = 3*(nd2v-1)*(nd2v+1) + 1
!      call mbyv_p( x, rk, mulext, muldatx, amcdv2(npnd0,1), nd2v,
!     &       mshxyz(mulpit(1)),
!     &       mulx1(1), mulx2(1), muly1(1), muly2(1), mulz1(1), mulz2(1),
!     &       mshnx(mulpms(1)), mshny(mulpms(1)), mshnz(mulpms(1)),
!     &       mshnod(1), mshnew(1), xx, nodexe )
!
!c      do i = 1, mshnod(1)
!c         if( abs(rk(i)-vk(i)).gt.1.d-06 ) then
!c         write(nfile(1),'(i5,10e16.8)') i,x(i),divl,vk(i)
!c         write(nfile(2),'(i5,10e16.8)') i,x(i),divl,vk(i)
!c         end if
!c      end do
!
!      cd12 = 0.d0
!      do i = 1, mshnod(1)
!         rki = vk(i) - rk(i)
!         cd12 = cd12 + rki*rki
!      end do
!      call gdsum(cd12,1,dbuf1r)
!      if( myid.eq.0 ) then
!          write(nfile(1),*) '  zansa            :', cd12
!          write(nfile(2),*) '  zansa            :', cd12
!      end if
!c--- end of check -------------------------------
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!    exchange & correlation   potential : vexc
!-----------------------------------------------------------------------
!    exchange energy       : eexc
!    exchange potential    : evexc
!    correlation energy    : ecor
!    correlation potential : evcor
!-----------------------------------------------------------------------
call pwpxc( nfile, myid, nodes,  &
& eexc, evexc, ecnl, ecor, evcor, lcstress, lvdw_pre, ltimecnt, t_comm )

if( ltimecnt ) then
    call qm_gsync
    ct = timecnt()
    if(loutfile(1)) write(nfile(1),*) ' Exchange correlation potential : cpu-time :', ct-ct0
    if(loutfile(2)) write(nfile(2),*) ' Exchange correlation potential : cpu-time :', ct-ct0
    ct0 = ct
end if


!-----------------------------------------------------------------------
!   local charge densities and energies in the PAW method
!if( lpaw ) then
!    call pwponsite( nfile, myid, nodes, lfssh_updtcg, lonsite, ltimecnt )
!
!    ct = timecnt()
!    if( ltimecnt ) then
! if(loutfile(1)) write(nfile(1),*) ' Onsite el. density & energies : cpu-time :',  &
!&                 ct-ct0
! if(loutfile(2)) write(nfile(2),*) ' Onsite el. density & energies : cpu-time :',  &
!&                 ct-ct0
!    end if
!    ct0 = ct
!
!end if


!-----------------------------------------------------------------------
!    kinetic energy : ekin
!    external potenital energy : eext
!-----------------------------------------------------------------------
elcl = 0.d0
do i = 1, mshnod(1)
   elcl = elcl + vext(i)*rho(i)
end do
call gdsum(elcl,1,dbuf1r)
elcl = elcl * rdelv


ctt0 = timecnt()
t_cpu2 = 0.d0

do nspin = 1, nspnmx
   if( lspin ) then
       call ldvloc( nspin, vlocud, vexc, mshnod(1) )
       if( lgamma ) then
           call ldocc( nspin, occ, wegud, nband )
!           if( lvand ) call ldslmi( nfile, myid, nodes, nspin, lvandi, .false. )
       end if
   end if

!--- unify local potential
call pwpunifylc( nfile, myid, nodes,  &
& hdiag, vexc, vhar, vhshdp, vext, mshnod(1), glocal, vlocud,  &
& mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d )

if( ltimecnt ) then
    ctt = timecnt()
    t_cpu2 = t_cpu2 + ctt - ctt0
    ctt0 = ctt
end if


!if( lvand ) then
    !----- screened Dij in ultrasoft pp
!    call caldij( nfile, myid, nodes,  &
!& lspin, nspin, lpaw, nplw3,  &
!& gboxx, gboxy, gboxz, ngboxa, ngboxb, ngboxc,  &
!& kfft1b, kfft2b, kfft3b, kfft0b, fft3x, fft3y, fftwork, ijkgb,  &
!& nband, occ, glocal, eigr, ntotfd, kfft0d, rvol,  &
!& ntype, natom, iatoit, ioa, ioag, lvandi, lking, bufatm, .false., .false. )
!    call calscrdij( nfile, myid, nodes,  &
!& lspin, nspin, lpaw, nplw3, ngboxa, ngboxb, ngboxc,  &
!& kfft1b, kfft2b, kfft3b, kfft0b, fft3x, fft3y, fftwork, ijkgb,  &
!& glocal, eigr, kfft0d, rvol,  &
!!& glocal, eigr, ntotfd, rvol,  &
!& ntype, natom, iatoit, ioa, ioag, lvandi, lking, bufatm )
!end if


#ifdef VECTOR
!-----convert glocal on FD mesh to on FFT mesh
call pwpgfd2ft( nfile, myid, nodes,  &
& glocal, apk, mfd2ft, ntotfd, kfft0d )
#endif

!----- glocal on dense grid -> coarse grid
call pwplocal_on_coarse( nfile, myid, nodes,  &
& glocal, thrhgr, nplw5ex, nplw5, nplw2ex, nplw2,  &
& nd1vks(1), nd1vks(2), nd1vks(3), kfft0d,  &
& mfd2ft, ntotfd, ijkgd, nplw5,  &
& kfft1,  kfft2,  kfft3,  kfft0,  mfd2ft_c, ntotfd_c, ijkg,  &
& fft3x, fft3y, fftwork )


   if( lspin ) then
#ifndef VECTOR
       call exrhcr( glocal, hlocal, ntotfd_c )
#else
       call exrhcr( glocal, hlocal, kfft0 )
#endif
   end if
end do


1000 continue
t_cpu1 = 0.d0
t_cpu3 = 0.d0
t_cpu4 = 0.d0
t_cpu5 = 0.d0
t_cpu6 = 0.d0
do i = 1, 10
   tc(i) = 0.d0
end do
ekin = 0.d0
encl = 0.d0
enclv = 0.d0
!-----zero clear onsite energies in the PAW method
EKIN1 = 0.d0
EHAR1 = 0.d0
EEXT1 = 0.d0
EEXCE1= 0.d0
ECORE1= 0.d0
edc1  = 0.d0
EplusU   = 0.d0
EplusUdc = 0.d0
edc_scissors = 0.d0
!---HF energy
exhfsr = 0.d0
!---electric field energy E*P
eefield = 0d0
!---fully relativistic SO energy
!call fullrela_clear


!---uniform electric field for insulator
!call periodic_efield_polarization( nfile, myid, nodes, ltimecnt, eefield )

zansa1 = 0.d0
zansa2 = 0.d0
call set_all2all_out( .false. )
spindo: do nspin = 1, nspnmx
   if( lspin ) then
       call ldocck( nspin, occ, wegud, nband, nkpnt )
       if( lgamma ) then
           call stspud( nspin, gdcr, gdcrsv, nspnod, lcgjsv )
!           if( lvand ) call ldslmi( nfile, myid, nodes, nspin, lvandi, .false. )
       end if
!       if( lvand ) then
!           call lddij( nfile, myid, nodes,  &
!& nspin, ntype, nhk1, nhk2, natom, lvandi, rvol )
!       end if
   end if


kvecdo3: do kvec = 1, nknod
!   if( .not.lgamma ) then
!       call set_kvec( kvec )
!       call get_wfk( gdcr, nspin )
!!       if( jhybrid /= 0 ) call get_hfrhcr( hfrhcr(1,nspin), nspin )
!       call get_plwk( nplw, nplwex, nspnod )
!       call get_indexk(  &
!& npnod1, npnod2, npnod, nplcnt, npldsp, nodes )
!       if( lvand ) then
!           call ldslmi_k( nfile, myid, nodes, lvandi, nspin, .false. )
!       end if
!   end if

   leig_end   = kvec == nknod .and. nspin == nspnmx
   if( ltimecnt ) ctt0 = timecnt()

!   --- copy gdcr to rhcr
   call cpgdrh( nfile, myid, nodes, rhcr, gdcr, npnod, nband )
!   --- to convert G decomposition to band decomposition ---
   call set_all2all_out( leig_end )
   call gdtobd( nfile, myid, nodes, ctt0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& rhcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, idstnd, ltimecnt )
!   --- copy rhcr to cgjr
   call cprhcg( nfile, myid, nodes, cgjr, rhcr, nplwex, nbnod )

!   if( lvand ) call slm_gdtobd( nfile, myid, nodes, ctt0,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& natom, natnod1, natnod2, natnod, natcnt, natdsp,  &
!& iod, iodg, ioag, idstnd, ltimecnt )
   call set_all2all_out( .false. )


!--- kinetic energy & HC products
!if( lgamma ) then
    call ekinetic( nfile, myid, nodes,  &
& ekin, ekib, occ, cgjr, rhcr, nplwex, nplw,  &
& nband, nbnod1, nbnod2, nbnod, iod, dkgnrm, dkrecnrm, 1 )
!else
!    call ekinetic_k( nfile, myid, nodes,  &
!& ekin, ekib, occ, cgjr, rhcr, nplwex, nplw, nspnmx,  &
!& nband, nbnod1, nbnod2, nbnod, iod, dkgnrm, nspin, 1 )
!end if

!      if(loutfile(1)) write(nfile(1),*) ' ekin=', ekin

if( ltimecnt ) then
    ctt = timecnt()
    t_cpu1 = t_cpu1 + ctt - ctt0
    ctt0 = ctt
end if


!if( lvand ) then
!    !-----clear hcsr
!    call clearhcsr( nfile, myid, nodes,  &
!& hcsr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod )
!end if

!--- nonlocal pp energy & HC products
!if( lgamma ) then
    call enonloc( nfile, myid, nodes,  &
& lspin, nspin, lpaw, encl, enclv, occ, glocal, cgjr, rhcr, hcsr,  &
& nplwex, nplw, nband, nbnod1, nbnod2, nbnod,  &
& iod, mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, rvol,  &
& fft3x, fft3y, fftwork, ijkg, nplw2,  &
& eigr, eigi, apk, apki, scwr, scwi,  &
& ntype, nhk1, nhk2, iatoit,  &
& lkbpp_r, lvand_r, lkbpp_g, lvand_g, lkbppi, lvandi, lking,  &
& rdelv_c, 1, lmax, mxl, lchk, lclno,  &
& nylmmx, natom, ycos, ysin, nplwcs, tc, ltimecnt )
!else
!    if( lnlpp_g ) then
!        call set_ycosysin( nfile, myid, nodes,  &
!& natom, ycos, ysin, nplwcs )
!    end if
!    if( lkbpp_r .or. lvand_r ) then
!        call set_eikr( nfile, myid, nodes,  &
!& mshglb_c, kfft1, kfft2, kfft3 )
!        call set_eikl( nfile, myid, nodes,  &
!& ratm, natom, lkbpp_r, lvand_r )
!    end if
!    call enonloc_k( nfile, myid, nodes,  &
!& lspin, nspin, lpaw, encl, enclv, occ, glocal, cgjr, rhcr, hcsr,  &
!& nplwex, nplw, nband, nbnod1, nbnod2, nbnod,  &
!& iod, mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, rvol,  &
!& fft3x, fft3y, fftwork, ijkg, nplw2,  &
!& eigr, eigi, apk, apki, scwr, scwi,  &
!& ntype, nhk1, nhk2, iatoit,  &
!& lkbpp_r, lvand_r, lkbpp_g, lvand_g, lkbppi, lvandi, lking,  &
!& rdelv_c, 1, lmax, mxl, lchk, lclno,  &
!& nylmmx, natom, tc, ltimecnt )
!end if


!---HSE hybrid functional: short-ranged HF energy
    lhfstart = nspin == 1 .and. kvec == 1
    lhfend   = nspin == nspnmx .and. kvec == nknod .and. nknod == nknodx
!    call HF_SR( nfile, myid, nodes, ltimecnt, exhfsr, nspin, .true.,  &
!& lhfstart, lhfend, hfmix )


!---uniform electric field for insulator
!    call periodic_efield( nfile, myid, nodes, ltimecnt, nspin )


!---scissors corrections to eigenvalues at donor/acceptor interface
!if( lscissors ) then
!    !---determine weight_donor for scissors corrections
!    call cal_weight_donor( nfile, myid, nodes,  &
!& nspnmx, nspin, cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod,  &
!& iod, iodg, mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, rvol,  &
!& fft3x, fft3y, fftwork, ijkg, nplw2, rba, zboundary )
!    !--- unify weight_donor
!    call unifyweight_donor( nfile, myid, nodes,  &
!& nspin, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, iod, iodg,  &
!& prod, prodr )
!end if

if( ltimecnt ) then
    ctt = timecnt()
    t_cpu3 = t_cpu3 + ctt - ctt0
    ctt0 = ctt
end if


!    --- copy rhcr to cgjr ---
call cprhcg( nfile, myid, nodes, cgjr, rhcr, nplwex, nbnod )
!    --- to convert band decomposition to G decomposition ---
if( .not.lvand .and. jhybrid == 0 ) call set_all2all_out( leig_end )
call bdtogd( nfile, myid, nodes, ctt0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& rhcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iodg, idstnd, .true., ltimecnt )

!if( lvand ) then
!!    --- copy hcsr to cgjr ---
!    call cprhcg( nfile, myid, nodes, cgjr, hcsr, nplwex, nbnod )
!!    --- to convert band decomposition to G decomposition ---
!    if( jhybrid == 0 ) call set_all2all_out( leig_end )
!    call bdtogd( nfile, myid, nodes, ctt0,  &
!& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& hcsr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
!& iodg, idstnd, .true., ltimecnt )
!end if

!if( jhybrid /= 0 ) then
!!    --- copy hfrhcr(1,nspin) to cgjr ---
!    call cprhcg( nfile, myid, nodes, cgjr, hfrhcr(1,nspin), nplwex, nbnod )
!!    --- to convert band decomposition to G decomposition ---
!    call set_all2all_out( leig_end )
!    call bdtogd( nfile, myid, nodes, ctt0,  &
!& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& hfrhcr(1,nspin), bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
!& iodg, idstnd, .true., ltimecnt )
!end if

!if( lefield_islts ) then
!    call add_efrhcr_to_hfrhcr( nfile, myid, nodes,  &
!& hfrhcr(1,nspin), npnod1, npnod2, npnod, nband, jhybrid )
!end if


!--- Hamiltonian matrix
if( lgamma ) then
!    if( jhybrid == 0 .and. .not.lefield_islts ) then
        call hmatrix( nfile, myid, nodes,  &
& eigtmp(1,nspin), gdcr, rhcr, npnod1, npnod2, npnod,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& iods, bx0r, dmtrxr, dmtrxi, nbxxxx, prod, prodr )
!    else
!        call hfhmatrix( nfile, myid, nodes,  &
!& eigtmp(1,nspin), gdcr, rhcr, hfrhcr(1,nspin), npnod1, npnod2, npnod,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& iods, bx0r, dmtrxr, nbxxxx, prod, prodr )
!    end if
else
!    if( jhybrid == 0 ) then
!        call hmatrix_k( nfile, myid, nodes,  &
!& eigtmp(1,nspin), gdcr, rhcr, npnod1, npnod2, npnod,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& iods, bx0r, dmtrxr, dmtrxi, nbxxxx, prod, prodr )
!    else
!        call hfhmatrix_k( nfile, myid, nodes,  &
!& eigtmp(1,nspin), gdcr, rhcr, hfrhcr(1,nspin), npnod1, npnod2, npnod,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& iods, bx0r, dmtrxr, dmtrxi, nbxxxx, prod, prodr )
!    end if
end if

if( ltimecnt ) then
    ctt = timecnt()
    t_cpu4 = t_cpu4 + ctt - ctt0
    ctt0 = ctt
end if

!--- gradients and residuals
!if( lgamma ) then
!    if( jhybrid == 0 .and. .not.lefield_islts ) then
        call gradres( nfile, myid, nodes,  &
& lvand, eigtmp(1,nspin), occ, gdcr, rhcr, hcsr, npnod1, npnod2, npnod,  &
& nband, zansa1, zansa2, prod, prodr )
!    else
!        call hfgradres( nfile, myid, nodes,  &
!& lvand, eigtmp(1,nspin), occ, gdcr, rhcr, hfrhcr(1,nspin), hcsr, npnod1, npnod2, npnod,  &
!& nband, zansa1, zansa2, prod, prodr )
!    end if
!else
!    if( jhybrid == 0 ) then
!        call gradres_k( nfile, myid, nodes,  &
!& lvand, eigtmp(1,nspin), occ, gdcr, rhcr, hcsr, npnod1, npnod2, npnod,  &
!& nband, zansa1, zansa2, prod, prodr )
!    else
!        call hfgradres_k( nfile, myid, nodes,  &
!& lvand, eigtmp(1,nspin), occ, gdcr, rhcr, hfrhcr(1,nspin), hcsr, npnod1, npnod2, npnod,  &
!& nband, zansa1, zansa2, prod, prodr )
!    end if
!end if

if( ltimecnt ) then
    ctt = timecnt()
    t_cpu5 = t_cpu5 + ctt - ctt0
    ctt0 = ctt
end if

!---scissors corrections to eigenvalues at donor/acceptor interface
!if( lscissors ) then
!    call cal_scissors_corr( nfile, myid, nodes,  &
!& eigtmp(1,nspin), edc_scissors, nspin, nband, ecorr_acceptor, ecorr_donor, nel, occ )
!end if


!if( .not.lgamma ) then
!    call set_rhcr( rhcr, nspin )
!    if( lvand ) call set_hcsr( hcsr, nspin )
!    if( jhybrid /= 0 ) call set_hfrhcr( hfrhcr(1,nspin), nspin )
!end if

end do kvecdo3

!   if( .not.lgamma .and. nknod < nknodx ) then
!       !---HSE hybrid functional: short-ranged HF energy
!       call HF_SR( nfile, myid, nodes, ltimecnt, exhfsr, nspin, .false.,  &
!& .false., .true., hfmix )
!   end if

   !---unify eig across k-point decomposition domains
!   call unify_eig( nfile, myid, nodes, eigtmp(1,nspin), nband )

   if( lspin ) then
!       call sveig_k( nspin, eig, egvlud, nband, nkpnt )
       call dexchange( dmtrxr, svtrxr, nbxxxx*nband*nknod )
       if( .not.lgamma ) then
           call dexchange( dmtrxi, svtrxi, nbxxxx*nband*nknod )
       end if

       if( lgamma ) then
           call exrhcr( rhcr, rhcrsv, nspnod )
!           if( lvand ) call exrhcr( hcsr, hcsrsv, nspnod )
       end if
#ifndef VECTOR
       call exrhcr( glocal, hlocal, ntotfd_c )
#else
       call exrhcr( glocal, hlocal, kfft0 )
#endif
   end if
end do spindo

!---unify some quantities
unwork(1) = ekin
unwork(2) = encl
unwork(3) = enclv
unwork(4) = zansa2
call unify_sumn( unwork, 4, unbuff )
ekin  = unwork(1)
encl  = unwork(2)
enclv = unwork(3)
zansa2= unwork(4)
zansa2= abs(zansa2)  !  zansa2 could possibly be negative, if occupancies are negative.
call unify_max1( zansa1 )

!---fully relativistic SO energy
!call fullrela_unify( soenergy )
soenergy = 0d0

enclp  = encl + soenergy
enclpv = enclv
!----- onsite energies in the PAW method
!if( lpaw ) then
!    call onsite_ene( nfile, myid, nodes,  &
!& EKIN1, EHAR1, EEXT1, EEXCE1, ECORE1, edc1,  &
!& EHAR1ns, EEXCE1ns, ECORE1ns, EplusU, EplusUdc )
!    encl  = 0.d0
!    enclv = 0.d0
!end if


eext = elcl + encl
zansa2 = zansa2/dble(nel)
!---convergence stabilizer
!nstabi = 10
!xstabi = 10.d0
lchkzansa = zansa2 < zansa2min * 100.d0
if( jstab >= nstabi ) then
    lchkzansa = lchkzansa .and. zansa2 < prevzansa2(jstab-nstabi) / xstabi
end if
if( .not.lchkzansa .and. .not.lfssh_updtcg ) then
    lreset = .true.
    lresetud = .true.   ! <-- required?
    zansa2min = zansa2
    jstab = 0
    if( lmixhold ) then
        lmixhold = .false.
        aslh = aslh0    !* 0.5d0
        if(loutfile(1)) write(nfile(1),'(a,f7.4)') &
& ' *** lmixhold = .false. / lreset = .true. / aslh =', aslh
        if(loutfile(2)) write(nfile(2),'(a,f7.4)') &
& ' *** lmixhold = .false. / lreset = .true. / aslh =', aslh
      else
        aslhmin = 0.2d0
        if( aslh-aslhmin > 0.0d0 ) then
            aslh = max( aslh * 0.5d0, aslhmin )
            if(loutfile(1)) write(nfile(1),'(a,f7.4)') ' *** lreset = .true. / aslh =', aslh
            if(loutfile(2)) write(nfile(2),'(a,f7.4)') ' *** lreset = .true. / aslh =', aslh
          else
            if(loutfile(1)) write(nfile(1),'(a,f7.4)') ' *** lreset = .true.'
            if(loutfile(2)) write(nfile(2),'(a,f7.4)') ' *** lreset = .true.'
        end if
    end if
end if
zansa2min = min( zansa2, zansa2min )
prevzansa2(jstab) = zansa2
jstab = jstab + 1

ct = timecnt()
if( ltimecnt ) then
    do ii = 1, 2
    if( loutfile(ii) ) then
       write(nfile(ii),*) ' Local/nonlocal potential: ', elcl, encl
       write(nfile(ii),*) '       kinetic energy & HC products ',  &
&                         ' :          :', t_cpu1
       write(nfile(ii),*) '              unify local potential ',  &
&                         ' :          :', t_cpu2
       write(nfile(ii),*) '   nonlocal pp energy & HC products ',  &
&                         ' :    wvg2r :', tc(1)
       write(nfile(ii),*) '   nonlocal pp energy & HC products ',  &
&                         ' :   calnlc :', tc(2)
       write(nfile(ii),*) '   nonlocal pp energy & HC products ',  &
&                         ' :    rfft3 :', tc(3)
       write(nfile(ii),*) '   nonlocal pp energy & HC products ',  &
&                         ' :    total :', t_cpu3
       write(nfile(ii),*) '                 Hamiltonian matrix ',  &
&                         ' :          :', t_cpu4
       write(nfile(ii),*) '            gradients and residuals ',  &
&                         ' :          :', t_cpu5
       write(nfile(ii),*) '                              total ',  &
&                         ' : cpu-time :', ct-ct0
!            else
!             write(nfile(ii),*) ' HC products                   ',
!     &                          ' : cpu-time :', ct-ct0
    end if
    end do
end if
!-----------------------------------------------------------------------


!--- total energy ------------------------------------------------------
elvlt = 0.136058D+02

!-----calculate sum of eigenvalues
!if( lspin ) then
!    call dsum_a_by_b( sume, wegud, egvlud, nband*nkpnt*2 )
!  else
!    call dsum_a_by_b( sume, occ, eig, nband*nkpnt )
!end if
if( lspin ) then
    if( nkpnt == 1 ) then
        call dsum_a_by_b( sume, wegud, eigtmp, nband*nkpnt*2 )
    else
        call sum_eig_spin( sume, wegud, eigtmp, nband, nkpnt )
    end if
else
    call dsum_a_by_b( sume, occ, eigtmp, nband*nkpnt )
end if
sume0 = sume
!-----sume: Harris-Foulkes (HF) energy functional
sume = sume - ehar + ( eexc - evexc ) + ( ecor - evcor ) - entrpy  &
&    + edc1 + EplusUdc + edc_scissors - exhfsr + eefield

!---unify sume by broadcasting
call unify_by_dbcast1( sume )

!-----etot: Kohn-Sham (KS) energy functional
etot = ekin + ehar + eext + eexc + ecor - entrpy + exhfsr + eefield + soenergy &
&    + EKIN1 + EHAR1 + EEXT1 + EEXCE1 + ECORE1 + EplusU

!---define format
fmtHF = '("  * Total energy (HF) = ",f15.6," ( Ryd. )",f13.4," (eV)")'
fmtKS = '("  * Total energy (KS) = ",F15.6," ( Ryd. )",F13.4," (eV)")'
fmtprts = '("  * Energy parts ( Ryd. / eV )"/  &
& "       Kinetic        External       Hartree        Exchange    Correlation"/  &
& 5f15.6/5(f13.4,2x) )'

if( ltimecnt ) then
    do i = 1, 2
    if( loutfile(i) ) then
        write(nfile(i),trim(fmtHF)) sume, sume*elvlt
        write(nfile(i),trim(fmtKS)) etot, etot*elvlt
!            else
!              write(nfile(1),1310) sume, etot
!              write(nfile(2),1310) sume, etot
    end if
    end do
end if
!1110 format('  * Total energy (HF) = ',F15.6,' ( Ryd. )',F13.4,' (eV)')
!1210 format('  * Total energy (KS) = ',F15.6,' ( Ryd. )',F13.4,' (eV)')
!1310 format('  * Total E. (HF,KS) =',2F15.6,' ( Ryd. )')

!-----output energy parts
if( ltimecnt ) then
    do i = 1, 2
    if( loutfile(i) ) then
       if( .not.lvdw .or. lvdw_pre ) then
           write(nfile(i),trim(fmtprts))  &
&       ekin, eext, ehar, eexc, ecor,  &
&       ekin*elvlt, eext*elvlt,  &
&       ehar*elvlt, eexc*elvlt, ecor*elvlt
       else
           write(nfile(i),1103)  &
&       ekin, eext, ehar, eexc, ecor-ecnl, &
&       ekin*elvlt, eext*elvlt,  &
&       ehar*elvlt, eexc*elvlt, (ecor-ecnl)*elvlt, &
        ecnl, ecnl*elvlt
       end if

!       if( lpaw ) then
!           write(nfile(i),1101)  &
!&       EKIN1, EEXT1, EHAR1, EEXCE1, ECORE1,  &
!&       EKIN1*elvlt, EEXT1*elvlt, EHAR1*elvlt, EEXCE1*elvlt,  &
!&       ECORE1*elvlt
!       if( lpaw_sym ) then
!       !-----full symmetry calculations---------------------------
!           write(nfile(i),1102)  &
!&       0.d0, 0.d0, EHAR1ns, EEXCE1ns, ECORE1ns,  &
!&       0.d0, 0.d0, EHAR1ns*elvlt, EEXCE1ns*elvlt,  &
!&       ECORE1ns*elvlt
!       end if
!       end if

!       if( lplusU ) then
!       !-----for the mean-field Hubbard model (DFT+U)-------------
!           write(nfile(i),'(a/F15.6/F13.4)')  &
!& '  * DFT+U energy ( Ryd. / eV )', EplusU, EplusU*elvlt
!       end if

!       if( jhybrid /= 0 ) then
!       !-----HSE hybrid functional: short-ranged HF energy-------------
!           write(nfile(i),'(a/F15.6/F13.4)')  &
!& '  * short-ranged HF energy ( Ryd. / eV )', exhfsr, exhfsr*elvlt
!       end if

!       if( lefield_islts ) then
!       !---uniform electric field for insulator
!           write(nfile(i),'(a/F15.6/F13.4)')  &
!& '  * electric field energy ( Ryd. / eV )', eefield, eefield*elvlt
!       end if

!       if( lrela .and. lrela_full ) then
!       !---fully relativistic SO energy
!           write(nfile(i),'(a/F15.6/F13.4)')  &
!& '  * fully relativistic SO energy ( Ryd. / eV )', soenergy, soenergy*elvlt
!       end if

!       if( lsoc .and. lsoc_full ) then
!       !---non relativistic SCF SO energy
!           write(nfile(i),'(a/F15.6/F13.4)')  &
!& '  * non-relativistic SCF SO energy ( Ryd. / eV )', soenergy, soenergy*elvlt
!       end if
    end if
    end do
end if
!1100 format(  &
!&      '  * Energy parts ( Ryd. / eV )'/  &
!&      2X,'     Kinetic        External       Hartree        ',  &
!&      'Exchange    Correlation'/  &
!&      5F15.6/5(F13.4,2X) )
!!     &           ' ',78('-') )
1103 format(  &
&      '  * Energy parts ( Ryd. / eV )'/  &
&      2X,'     Kinetic        External       Hartree        ',  &
&      'Exchange   EcLDA / Ecnl'/  &
&      5F15.6/5(F13.4,2X) / 60x,F15.6/60x,(F13.4,2X) )
1101 format(  &
&      '  * On site energy parts ( Ryd. / eV )',  &
&      '   total / non-spherical E.'/  &
&       5F15.6/5(F13.4,2X) )
1102 format( 5F15.6/5(F13.4,2X) )


if( loutfile(1) ) then
if( ifmd.eq.0 .or. loutenergy ) then
!    if( lpaw ) then
!        eonsite = EKIN1 + EHAR1 + EEXT1 + EEXCE1 + ECORE1
!      else
        eonsite = 0.d0
!    end if
    write(nfile(5),1150) nstep, jgcycl, sume, etot,  &
&   ekin, eext, ehar, eexc, ecor-ecnl, ecnl, entrpy, eonsite, exhfsr, EplusU, eefield, soenergy
end if
end if
1150 format(2i6,2es18.10,20es14.6)
!-----------------------------------------------------------------------
!    convergence criterion
difene = abs( sume - eneold )/dble(nel)
difene2 = abs( sume - etot )/dble(nel)
ct = timecnt()
do i = 1, 2
if( loutfile(i) ) then
if( ltimecnt ) then
    if( jcycl.ge.1 ) then
        write(nfile(i),'(a22,2es12.4)')  &
&                            '  * Diff. total E/el.=', difene, difene2
        write(nfile(i),'(a22,2es12.4,2x,2es12.4)')  &
&                     '  * max.& av.residual=', zansa1, zansa2, bfzansa1, bfzansa2
    else
        write(nfile(i),'(a22,2es12.4)')  &
&                     '  * max.& av.residual=', zansa1, zansa2
    end if
    write(nfile(i),'(1x,a20,f10.4)')  &
&                     ' cpu-time ( sec )  :', ct - ct00
else
    if( jcycl == 0 ) then
        write(nfile(i),*) ' '
        write(nfile(i),*) '    E(HF) (Ryd.)   E(KS) (Ryd.)  ',  &
!&          'diff /el.  zansa1     zansa2        sec'
&          'diff /el.  diff2/el.  zansa2        sec'
    end if
    if( jcycl < 100) then
        write(nfile(i),'(i2,a1,2es15.7,es11.4,2es11.4,f9.3)')  &
!&        jcycl, ':', sume, etot, difene, zansa1, zansa2, ct - ct00
&        jcycl, ':', sume, etot, difene, difene2, zansa2, ct - ct00
    else
        write(nfile(i),'(i3,a1,es14.7,es15.7,es11.4,2es11.4,f9.3)')  &
!&        jcycl, ':', sume, etot, difene, zansa1, zansa2, ct - ct00
&        jcycl, ':', sume, etot, difene, difene2, zansa2, ct - ct00
    end if
end if
end if
end do
eneold = sume

if( jcycl < 1 ) then
    bfzansa1 = zansa1
    bfzansa2 = zansa2
end if

!if( lefield .and. lsawtooth .and. loutpolarization ) then
!    !--- electronic polarization with sawtooth electric field
!    call sawtooth_efield_pol( nfile, myid, nodes,  &
!& lsawtooth_shape, lsawtooth_xyz, rho,  &
!& lclust, rdel, rdelg, hcell, lorthrhmbc, nd1v, nd1vks, rdelv,  &
!& mshxyz(mulpit(1)),  &
!& mulx1(1), mulx2(1), muly1(1), muly2(1), mulz1(1), mulz2(1),  &
!& mshnx(mulpms(1)), mshny(mulpms(1)), mshnz(mulpms(1)), mshnod(1),  &
!& mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1) )
!end if
if( loutfile(1) ) then
    if( ifmd.eq.0 .or. loutzansa ) then
        write(nfile(29),'(2i6,6es12.4)') nstep, jgcycl, difene, difene2, zansa1, zansa2,  &
& bfzansa1, bfzansa2
    end if

!    if( lefield ) then
!        !---uniform electric field for insulator
!        call out_polarization( nfile, myid, nodes, loutpolarization, nstep, jgcycl )
!    end if
end if

      !==== stop here =========================
!      call fnormalstop( nfile, myid, nodes, 'pwlda' )
      !========================================

ct00 = ct
ct0  = ct
if( lvdw .and. zansa2 < max( 6.2d-07, tolres*10.d0 ) ) lvdw_pre = .false.

!if( laspc_exec ) then
!    lconverge = laspc_conv .or. zansa2 <= tolres_aspc
!else
    lconverge = difene.le.tolpot .or. zansa2.le.tolres .or. jcycl.ge.iscfmx  &
&          .or. jcycl > 1 .and. difene2.le.tolHF_KS
!end if
lconverge = lconverge .and. .not.lncinit
lncinit = .false.
if( lconverge .and. .not.lfssh_gsscf .or. lfssh_updtcg ) exit outer
scfif: if( .not.lconverge .or. .not.lfssh_gsscf ) then
!-----------------------------------------------------------------------

jcycl  = jcycl + 1
jgcycl = jgcycl + 1
if( ltimecnt ) then
    if(loutfile(1)) write(nfile(1),*) ('-',i=1,20),' iteration:',jcycl,' ',  &
&                         ('-',i=1,20)
    if(loutfile(2)) write(nfile(2),*) ('-',i=1,20),' iteration:',jcycl,' ',  &
&                         ('-',i=1,20)
end if


lhcexe = .false.
ct0 = timecnt()
csh0 = ct0
!-----------------------------------------------------------------------
!   subspace rotation (Unitary transformation)
!cc      if( lstart .or. jcycl.gt.1 .or. nstep.ge.nstep_ini+1 ) then
if( ltimecnt ) then
    if(loutfile(1)) write(nfile(1),*) ' Subspace rotation ( Unitary Tr. )  '
    if(loutfile(2)) write(nfile(2),*) ' Subspace rotation ( Unitary Tr. )  '
end if
do nspin = 1, nspnmx
   if( lspin ) then
       if( lgamma ) then
           call stspud( nspin, gdcr, gdcrsv, nspnod, lcgjsv )
           if( nspin == 1 .or. nspin == 2 .and. .not.lduplicate ) then
!               if( lvand ) call ldslmi( nfile, myid, nodes, nspin, lvandi, .false. )
           end if
       end if
   end if

   dupif: if( nspin == 1 .or. nspin == 2 .and. .not.lduplicate ) then

   !---set nspin in tddft-fssh.f90
   call set_nspin_in_tddft_fssh( nspin )

kvecdo4: do kvec = 1, nknod
!   if( .not.lgamma ) then
!       call set_kvec( kvec )
!       call get_wfk( gdcr, nspin )
!       call get_rhcr( rhcr, nspin )
!       if( lvand ) call get_hcsr( hcsr, nspin )
!       if( jhybrid /= 0 ) call get_hfrhcr( hfrhcr(1,nspin), nspin )
!       call get_plwk( nplw, nplwex, nspnod )
!       call get_indexk(  &
!& npnod1, npnod2, npnod, nplcnt, npldsp, nodes )
!       if( lvand ) then
!           call ldslmi_k( nfile, myid, nodes, lvandi, nspin, .false. )
!       end if
!   end if

 leig_start = kvec == 1
 leig_end   = kvec == nknod

!--- unitary transformation --------------------------------------------
if( lgamma ) then

!    if( .not.lefield_islts ) then
        call untryt( nfile, myid, nodes, csh0,  &
& gdcr, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& prod, prodr, prods, betk, w1r, w1i, bzan, ew1, tmpr, iblock,  &
& ltimecnt, leig_start, leig_end,  &
& bijr, biji, dmtrxr, dmtrxi, nbxxxx,  &
& node_r, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, bm1r, bm1i, bx0r )
!    else
!        !---uniform electric field for insulator
!        !---set the number of the occupied states
!        call set_occupied_states_in_efield( nspin )
!        !---set the Unitary matrix
!        call efuntryt( nfile, myid, nodes, csh0,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& prod, prodr, prods, betk, w1r, w1i, bzan, ew1,  &
!& ltimecnt, leig_start, leig_end,  &
!& bijr, biji, dmtrxr, dmtrxi, nbxxxx,  &
!& node_r, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, bm1r, bm1i, bx0r )
!        call untryhc( nfile, myid, nodes, csh0,  &
!& gdcr, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& prod, prodr, w1r, w1i, tmpr, iblock,  &
!& ltimecnt, leig_start, leig_end .and. .not.lvand .and. .not.laspc_rot,  &
!& bijr, biji, nbxxxx )
!    end if

    if( lhcunt ) then
        lhcexe = .true.
        !------ unitary transformation of HC products
        call untryhc( nfile, myid, nodes, csh0,  &
& rhcr, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& prod, prodr, w1r, w1i, tmpr, iblock,  &
& ltimecnt, leig_start, leig_end .and. .not.lvand .and. .not.laspc_rot,  &
& bijr, biji, nbxxxx )
        !------ unitary transformation of SC products
!        if( lvand ) then
!            call untryhc( nfile, myid, nodes, csh0,  &
!& hcsr, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& prod, prodr, w1r, w1i, tmpr, iblock,  &
!& ltimecnt, .false., leig_end .and. .not.laspc_rot,  &
!& bijr, biji, nbxxxx )
!            call trans_slmi_t( nfile, myid, nodes,  &
!& bijr, biji, nband, nbnod1, nbnod2, nbnod, nbxxxx,  &
!& ntype, nhk1_nat, nhk2_nat, natom, lvandi, lking, prod, prodr )
!        end if
        !------ unitary transformation of predicted w.f. in the ASPC method
        if( laspc_rot ) then
            call untryhc( nfile, myid, nodes, csh0,  &
& pgdcr, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& prod, prodr, w1r, w1i, tmpr, iblock,  &
& ltimecnt, .false., leig_end,  &
& bijr, biji, nbxxxx )
        end if
!        if( jhybrid /= 0 .or. lefield_islts ) then
!        !------ unitary transformation of short-range HF * C products
!            call untryhc( nfile, myid, nodes, csh0,  &
!& hfrhcr(1,nspin), npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& prod, prodr, w1r, w1i, tmpr, iblock,  &
!& ltimecnt, .false., leig_end,  &
!& bijr, biji, nbxxxx )
!        end if
!        if( lefield_islts ) then
!        !------ unitary transformation of the inverse of the S matrix
!            call untrysinv( nfile, myid, nodes,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& prod, prodr, tmpr, iblock, bijr, biji, nbxxxx, nspin )
!        end if
    end if

else

!    call untryt_k( nfile, myid, nodes, csh0,  &
!& gdcr, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& prod, prodr, prods, betk, w1r, w1i, bzan, ew1, tmpr, iblock,  &
!& ltimecnt, leig_start, leig_end,  &
!& bijr, biji, dmtrxr, dmtrxi, nbxxxx,  &
!& node_r, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, bm1r, bm1i, bx0r )
!    if( lhcunt ) then
!        lhcexe = .true.
!        !------ unitary transformation of HC products
!        call untryhc_k( nfile, myid, nodes, csh0,  &
!& rhcr, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& prod, prodr, w1r, w1i, tmpr, iblock,  &
!& ltimecnt, leig_start, leig_end .and. .not.lvand,  &
!& bijr, biji, nbxxxx )
!        !------ unitary transformation of SC products
!        if( lvand ) then
!            call untryhc_k( nfile, myid, nodes, csh0,  &
!& hcsr, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& prod, prodr, w1r, w1i, tmpr, iblock,  &
!& ltimecnt, .false., leig_end,  &
!& bijr, biji, nbxxxx )
!            call trans_slmi_t_k( nfile, myid, nodes,  &
!& bijr, biji, nband, nbnod1, nbnod2, nbnod, nbxxxx,  &
!& ntype, nhk1_nat, nhk2_nat, natom, lvandi, lking, prod, prodr )
!        end if
!        if( jhybrid /= 0 ) then
!            call untryhc_k( nfile, myid, nodes, csh0,  &
!& hfrhcr(1,nspin), npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& prod, prodr, w1r, w1i, tmpr, iblock,  &
!& ltimecnt, .false., leig_end,  &
!& bijr, biji, nbxxxx )
!        end if
!    end if

end if


!if( .not.lgamma ) then
!    call set_wfk( gdcr, nspin )
!    call set_rhcr( rhcr, nspin )
!    if( lvand ) then
!        call set_hcsr( hcsr, nspin )
!        call svslmi_k( nfile, myid, nodes, lvandi, nspin, .false. )
!    end if
!    if( jhybrid /= 0 ) call set_hfrhcr( hfrhcr(1,nspin), nspin )
!    if( lspin .and. lduplicate ) then
!        call set_wfk( gdcr, 2 )
!        call set_rhcr( rhcr, 2 )
!        if( lvand ) then
!            call set_hcsr( hcsr, 2 )
!            call svslmi_k( nfile, myid, nodes, lvandi, 2, .false. )
!        end if
!        if( jhybrid /= 0 ) call set_hfrhcr( hfrhcr(1,1), 2 )
!    end if
!end if

end do kvecdo4

   if( lspin ) then
       if( .not.lduplicate ) then
           call dexchange( dmtrxr, svtrxr, nbxxxx*nband*nknod )
           if( .not.lgamma ) then
               call dexchange( dmtrxi, svtrxi, nbxxxx*nband*nknod )
           end if

           if( lgamma ) then
               call exrhcr( rhcr, rhcrsv, nspnod )
!               if( lvand ) then
!                   call exrhcr( hcsr, hcsrsv, nspnod )
!                   call svslmi( nfile, myid, nodes, nspin, lvandi, .false. )
!               end if
               if( laspc_rot ) call exrhcr( pgdcr, pgdcrsv, nspnod )
           end if
       else
           call dcopy_a_to_b( dmtrxr, svtrxr, nbxxxx*nband*nknod )
           if( .not.lgamma ) then
               call dcopy_a_to_b( dmtrxi, svtrxi, nbxxxx*nband*nknod )
           end if

           if( lgamma ) then
               call dcopy_a_to_b( gdcr, gdcrsv, nspnod )
               call dcopy_a_to_b( rhcr, rhcrsv, nspnod )
!               if( lvand ) then
!                   call dcopy_a_to_b( hcsr, hcsrsv, nspnod )
!                   call svslmi( nfile, myid, nodes, nspin, lvandi, .false. )
!                   call svslmi( nfile, myid, nodes, 2, lvandi, .false. )
!               end if
               if( laspc_rot ) call dcopy_a_to_b( pgdcr, pgdcrsv, nspnod )
           end if
       end if
   end if

   end if dupif
end do
ct = timecnt()
if( ltimecnt ) then
    if(loutfile(1)) write(nfile(1),*) '                               total',  &
&                         ' : cpu-time :', ct-ct0
    if(loutfile(2)) write(nfile(2),*) '                               total',  &
&                         ' : cpu-time :', ct-ct0
!              write(nfile(1),*) '                                    ',
!     &                          ' : com-time :', t_comm
!              write(nfile(2),*) '                                    ',
!     &                          ' : com-time :', t_comm
!            else
!              write(nfile(1),*) ' Subspace rotation             ',
!     &                          ' : cpu-time :', ct-ct0
!              write(nfile(2),*) ' Subspace rotation             ',
!     &                          ' : cpu-time :', ct-ct0
end if
ct0 = ct
!cc      end if
!-----------------------------------------------------------------------


if( ltimecnt ) then
    if(loutfile(1)) write(nfile(1),*) ' Eigenvalue problem'
    if(loutfile(2)) write(nfile(2),*) ' Eigenvalue problem'
end if
!-----------------------------------------------------------------------
!    solving eigenvalue problem by iteration method
!-----------------------------------------------------------------------
!---set lstop in schmidt.f90
call setlstop( nfile, myid, nodes, .false. )

bfzansa1 = 0.d0
bfzansa2 = 0.d0
t_comm = 0.d0
lortho = .false.
if( nodes.ge.2 ) lortho = .false.
do nspin = 1, nspnmx
   if( lspin ) then
       if( lgamma ) then
           call stspud( nspin, gdcr, gdcrsv, nspnod, lcgjsv )
!           if( lvand ) then
!               call ldslmi( nfile, myid, nodes, nspin, lvandi, .false. )
!           end if
       end if
!       if( lvand ) then
!           call lddij( nfile, myid, nodes,  &
!& nspin, ntype, nhk1, nhk2, natom, lvandi, rvol )
!       end if
       if( ltimecnt ) then
           do i = 1, 2
           if( loutfile(i) ) then
               if( nspin.eq.1 ) then
                   write(nfile(i),*) '  ( for up-spin )'
               else
                   write(nfile(i),*) '  ( for down-spin )'
               end if
           end if
           end do
       end if
   end if

   dupif2: if( nspin == 1 .or. nspin == 2 .and. .not.lduplicate ) then

   !---set nspin in tddft-fssh.f90
   call set_nspin_in_tddft_fssh( nspin )
   !---set nspin in nlvand.f90
!   if( jhybrid /= 0 .or. lefield_islts ) call set_nspin_in_nlvand( nspin )


call set_all2all_out( .false. )
kvecdo5: do kvec = 1, nknod
!   if( .not.lgamma ) then
!       call set_kvec( kvec )
!       call get_wfk( gdcr, nspin )
!       call get_rhcr( rhcr, nspin )
!       if( lvand ) call get_hcsr( hcsr, nspin )
!       if( jhybrid /= 0 ) call get_hfrhcr( hfrhcr(1,nspin), nspin )
!       call get_plwk( nplw, nplwex, nspnod )
!       call get_indexk(  &
!& npnod1, npnod2, npnod, nplcnt, npldsp, nodes )
!       call get_plw7k( nplw7, nplw7ex, nspnod7 )
!       call get_index7k(  &
!& npnod71, npnod72, npnod7, nplcnt7, npldsp7, nodes )
!#if CCREAL4
!       call get_indexncg7k( ncgcnt7, ncgdsp7, nodes )
!#else
!       call get_indexncgk( ncgcnt, ncgdsp, nodes )
!#endif
!       if( lvand ) then
!           call ldslmi_k( nfile, myid, nodes, lvandi, nspin, .false. )
!       end if
!
!       if( ltimecnt ) then
!           k = kvec + nknod1 - 1
!           if(loutfile(1)) write(nfile(1),'(2x,a3,i5,a2,3f8.4,a17,l2)')  &
!& 'k =', k, ' (', (bzk(j,k),j=1,3), ' ) with lgammak =', lgammak(k)
!           if(loutfile(2)) write(nfile(2),'(2x,a3,i5,a2,3f8.4,a17,l2)')  &
!& 'k =', k, ' (', (bzk(j,k),j=1,3), ' ) with lgammak =', lgammak(k)
!        end if
!   end if

!   if( lefield_islts ) call set_occupied_states_in_efield( nspin )

 leig_start = kvec == 1
 leig_end   = kvec == nknod

!      if( ihldam.eq.1 ) then
    if( ltimecnt ) csh0 = timecnt()
!--- convert G decomposition to band decomposition ---------------------
    if( lhcexe .and. .not.lortho ) then
        call gdtobd( nfile, myid, nodes, csh0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& rhcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, idstnd, ltimecnt )
!        if( lvand ) then
!            call gdtobd( nfile, myid, nodes, csh0,  &
!& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& hcsr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
!& iod, iodg, idstnd, ltimecnt )
!        end if
    end if
!    if( jhybrid /= 0 .or. lefield_islts ) then
!        call gdtobd( nfile, myid, nodes, csh0,  &
!& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& hfrhcr(1,nspin), bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
!& iod, iodg, idstnd, ltimecnt )
!!        if( .not.lgamma ) call set_hfrhcr( hfrhcr(1,nspin), nspin )
!    end if
!--- unify wavefunctions -----------------------------------------------
    call set_all2all_out( leig_end )
!=========
orthdo: do
!=========
#if CCREAL4
!--- copy gdcr (real*8) to gdcr4 (real*4)
!          call cpgdg4( nfile, myid, nodes, gdcr4, gdcr, npnod, nband )
!          call unifywv( nfile, myid, nodes, csh0,
!     & cgjr4,nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,
!     & gdcr4,bufcr4,npnod1, npnod2, npnod, nplcnt, npldsp,
!     & iods, iodsg, idstnd, ncgcnt, ncgdsp, nspnod, npnnbnx, ltimecnt )
    if( npnod7 == npnod ) then
        call cpgdg4( nfile, myid, nodes,  &
& gdcr4, gdcr, npnod7, nband )
      else
        call cpgdg47( nfile, myid, nodes,  &
& gdcr4, gdcr, npnod, npnod7, nband )
    end if
    call unifywv( nfile, myid, nodes, csh0,  &
& cgjr4, nplw7ex, nplw7, nband, nbnod1, nbnod2,nbnod,nbncnt,nbndsp,  &
& gdcr4, bufcr4, npnod71, npnod72, npnod7, nplcnt7, npldsp7,  &
& iods, iodsg, idstnd, ncgcnt7, ncgdsp7, nspnod7, npnnbnx,  &
& ltimecnt )
!--- convert G decomposition to band decomposition 
    call gdtobd( nfile, myid, nodes, csh0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, idstnd, ltimecnt )
#else
    call unifywv( nfile, myid, nodes, csh0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iods, iodsg, idstnd, ncgcnt, ncgdsp, nspnod, npnnbnx, ltimecnt )
!--- copy cgjr to gdcr in correct order
    call cpcgrd( nfile, myid, nodes,  &
& cgjr, gdcr, nplwex, nbnod, iod )
#endif

!    if( lvand ) then
        !--- convert atom decomposition to band decomposition 
!        call slm_gdtobd( nfile, myid, nodes, csh0,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& natom, natnod1, natnod2, natnod, natcnt, natdsp,  &
!& iod, iodg, ioag, idstnd, ltimecnt )

!        if( lgamma ) then
!            !--- unify slmir, slmii
!            call unifyslm( nfile, myid, nodes,  &
!& nband, nbnod1, nbnod2, nbnod, iod, iodg,  &
!& ntype, nhk1, nhk2, natom, lvandi, bijr, nbxxxx, .true. )
!        end if
!    end if

    !---HSE hybrid functional: short-ranged HF energy
!    call HF_SR_self( nfile, myid, nodes, ltimecnt, nspin )

    !---uniform electric field for insulator
!    call periodic_efield_self( nfile, myid, nodes, ltimecnt, nspin )

!    if( lvand ) then
!        if( .not.lgamma ) then
!#if CCREAL4
!#else
            !--- unify slmir, slmii
!            call unifyslm( nfile, myid, nodes,  &
!& nband, nbnod1, nbnod2, nbnod, iod, iodg,  &
!& ntype, nhk1, nhk2, natom, lvandi, bijr, nbxxxx, .true. )
!!#endif
!        end if
!    end if

!--- conjugate-gradient minimization -----------------------------------
if( jgcycl == kstrial-1 ) then
    lreset = .true.
    lresetud = .true.   ! <-- required?
    if(loutfile(1)) write(nfile(1),'(a,f7.4)') ' *** lreset = .true.'
    if(loutfile(2)) write(nfile(2),'(a,f7.4)') ' *** lreset = .true.'
end if
if( jgcycl >= kstrial ) then
    itermx_  = itermx
    iteremx_ = iteremx
else
    itermx_  = 1
    iteremx_ = 1
end if
if( lgamma ) then
    call eigen( nfile, myid, nodes,  &
#if CCREAL4
& eig, gdcr, cgjr4, rhcr, hcsr, nplwex, nplw, nplw7ex,  &
#else
& eig, gdcr, cgjr, rhcr, hcsr, nplwex, nplw, nplwex,  &
#endif
& dkgnrm, prcd, ekib,  &
& gnk, hnk, sck, nband, nbnod1, nbnod2, nbnod, iod,  &
& glocal, mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, hfv, hfv_c,  &
& fft3x, fft3y, fftwork, ijkg, nplw2,  &
& eigr, eigi, apk, apki, scwr, scwi, rdelv_c, rvol,  &
& ntype, nhk1, nhk2, natom, iatoit, lvand,  &
& lkbpp_r, lvand_r, lkbpp_g, lvand_g, lkbppi, lvandi, lking,  &
& prod, prodr, w1r, w1i, iter, istdot, toleig, threinn, wdinner, lortho,  &
& bzan, nriter, ihyoji, itermx_, noband, iteremx_,  &
& lhcexe, ltimecnt, leig_start, leig_end,  &
& lmax, mxl, lchk, lclno, nylmmx, ycos, ysin, nplwcs, methodcg,  &
& hfrhcr(1,nspin), jhybrid, lefield_islts )
else
!    if( lnlpp_g ) then
!        call set_ycosysin( nfile, myid, nodes,  &
!& natom, ycos, ysin, nplwcs )
!    end if
!    if( lkbpp_r .or. lvand_r ) then
!        call set_eikr( nfile, myid, nodes,  &
!& mshglb_c, kfft1, kfft2, kfft3 )
!        call set_eikl( nfile, myid, nodes,  &
!& ratm, natom, lkbpp_r, lvand_r )
!    end if
!    call eigen_k( nfile, myid, nodes,  &
!#if CCREAL4
!& eig, gdcr, cgjr4, rhcr, hcsr, nplwex, nplw, nplw7ex,  &
!#else
!& eig, gdcr, cgjr, rhcr, hcsr, nplwex, nplw, nplwex,  &
!#endif
!& dkgnrm, prcd, ekib, nspin,  &
!& gnk, hnk, sck, nband, nbnod1, nbnod2, nbnod, iod, nspnmx,  &
!& glocal, mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, hfv, hfvi, hfv_c, hfvi_c,  &
!& fft3x, fft3y, fftwork, ijkg, nplw2,  &
!& eigr, eigi, apk, apki, scwr, scwi,  &
!& rdelv_c, rvol,  &
!& ntype, nhk1, nhk2, natom, iatoit, lvand,  &
!& lkbpp_r, lvand_r, lkbpp_g, lvand_g, lkbppi, lvandi, lking, ratm,  &
!& prod, prodr, w1r, w1i, iter, istdot, toleig, threinn, wdinner, lortho,  &
!& bzan, nriter, ihyoji, itermx_, noband, iteremx_,  &
!& lhcexe, ltimecnt, leig_start, leig_end,  &
!& lmax, mxl, lchk, lclno, nylmmx, methodcg,  &
!& hfrhcr(1,nspin), jhybrid, mxkpdupnm, cgjr_m, mshglb_c, kfft1, kfft2, kfft3, lefield_islts )
end if
!        else
!--- RMM-DIIS
!              call rmmdis( nband, hdiag, weigrd, iter, istdot, toleig,
!     &     rdelv, prod, prodr, rk, apk, betk, bzan, nriter, vk, x, xx,
!     &     dbuf, dbufr, noddatx, muldatx, mulext, nodexe, ihyoji,
!     &     itermx_, noband, iteremx_, ltimecnt, t_comm )
!      end if

!---calculation of bfzansa
k = kvec + nknod1 - 1
if( lspin ) then
    call cal_bfzansa( nfile, myid, nodes,  &
& bfzansa1, bfzansa2, bzan, wegud, nband, nbnod1, nbnod2, nbnod, iod, nkpnt, nspin, 2, k )
else
    call cal_bfzansa( nfile, myid, nodes,  &
& bfzansa1, bfzansa2, bzan, occ, nband, nbnod1, nbnod2, nbnod, iod, nkpnt, nspin, 1, k )
end if

!if( lnoncollinear ) then
!    !--- calculate wvratio : ratio of each component, if necessary.
!    call cal_wvratio_k( nfile, myid, nodes,  &
!& wvratio, wvmagxy, gdcr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, iod,  &
!& rvol, ntype, nhk1, nhk2, natom, lvand, lvandi )
!end if

!      !==== stop here =========================
!      call fstop( nfile, myid, nodes, 'pwlda2' )
!      !========================================

if( ltimecnt ) csh0 = timecnt()
!    --- copy gdcr to cgjr ---
   call cprhcg( nfile, myid, nodes, cgjr, gdcr, nplwex, nbnod )
!--- all to all communication ------------------------------------------
!--- to convert band decomposition to G decomposition
call bdtogd( nfile, myid, nodes, csh0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iodg, idstnd, .true., ltimecnt )

!if( lvand ) call slm_bdtogd( nfile, myid, nodes, csh0,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& natom, natnod1, natnod2, natnod, natcnt, natdsp,  &
!& iodg, ioag, idstnd, .true., ltimecnt )

if( .not.lortho ) then
!if( lgamma ) then

!--- Gram-Schmidt orthonormalizaion ------------------------------------
!    call schmidt( nfile, myid, nodes, csh0,  &
!& lvand, gdcr, npnod1, npnod2, npnod, nband, nbnod1, nbnod2, nbnod,  &
!& nbncnt, nbndsp, iod, iodg, prod, prodr, pdbuf,  &
!& ltimecnt, .false., leig_start, leig_end,  &
!& dmtrxr, nbxxxx,  &
!& node_c, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd, bijr, bm1r,  &
!& rvol, ntype, nhk1_nat, nhk2_nat, natom, lvandi )
    call schmidt( nfile, myid, nodes, csh0,  &
& gdcr, iod, iodg, ltimecnt, .false., leig_start, leig_end )
!--- unify eigenvalues -------------------------------------------------
    call unifyeg( nfile, myid, nodes, csh0,  &
& eig, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, iod, iodg,  &
& prod, prodr, ltimecnt, leig_start, leig_end )

!else

!--- Gram-Schmidt orthonormalizaion ------------------------------------
!   if( lgamma .or. lgammak(kvec+nknod1-1) ) then
!       call schmidt( nfile, myid, nodes, csh0,  &
!& gdcr, iod, iodg, ltimecnt, .false., leig_start, leig_end )
!   else
!       call schmidt_k( nfile, myid, nodes, csh0,  &
!& gdcr, iod, iodg, ltimecnt, .false., leig_start, leig_end )
!   end if
!!--- unify eigenvalues -------------------------------------------------
!    call unifyeg_k( nfile, myid, nodes, csh0,  &
!& eig, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, iod, iodg,  &
!& prod, prodr, ltimecnt, leig_start, leig_end )
!
!end if
end if

!if( lnoncollinear ) then
!    !--- calculate wvratio : ratio of each component, if necessary.
!    call unify_wvratio( nfile, myid, nodes,  &
!& wvratio, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, iod, iodg,  &
!& prod, prodr )
!    call unify_wvratio( nfile, myid, nodes,  &
!& wvmagxy, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, iod, iodg,  &
!& prod, prodr )
!end if

!---check whether Gram-Schmidt orthonormalizaion is successfull or not.
call getlcholesky( nfile, myid, nodes, lcholesky )
if( lcholesky .or. lortho ) exit orthdo
!============
end do orthdo
!============


!if( .not.lgamma ) then
!    call set_wfk( gdcr, nspin )
!    call set_rhcr( rhcr, nspin )
!    if( lvand ) then
!        call set_hcsr( hcsr, nspin )
!        call svslmi_k( nfile, myid, nodes, lvandi, nspin, .false. )
!    end if
!    if( lduplicate ) then
!        call set_wfk( gdcr, 2 )
!        call set_rhcr( rhcr, 2 )
!        if( lvand ) then
!            call set_hcsr( hcsr, 2 )
!            call svslmi_k( nfile, myid, nodes, lvandi, 2, .false. )
!        end if
!    end if
!end if

end do kvecdo5

   !---scissors corrections to eigenvalues at donor/acceptor interface
!   if( lscissors ) then
!       call cal_scissors_corr( nfile, myid, nodes,  &
!& eig, edc_scissors, nspin, nband, ecorr_acceptor, ecorr_donor, nel, occ )
!   end if


   !---unify eig across k-point decomposition domains
!   call unify_eig( nfile, myid, nodes, eig, nband )

!   if( lnoncollinear ) then
!       !--- unify wvratio across k-point decomposition domains, if necessary.
!       call unify_wvratio2( nfile, myid, nodes, wvratio, nband )
!       call unify_wvratio2( nfile, myid, nodes, wvmagxy, nband )
!    end if

   if( lspin ) then
       call sveig_k( nspin, eig, egvlud, nband, nkpnt )
       if( .not.lduplicate ) then
           if( lgamma ) then
               call exrhcr( rhcr, rhcrsv, nspnod )
!               if( lvand ) then
!                   call exrhcr( hcsr, hcsrsv, nspnod )
!                   call svslmi( nfile, myid, nodes, nspin, lvandi, .false. )
!               end if
           end if
       else
           call sveig_k( 2, eig, egvlud, nband, nkpnt )
           if( lgamma ) then
               call dcopy_a_to_b( gdcr, gdcrsv, nspnod )
               call dcopy_a_to_b( rhcr, rhcrsv, nspnod )
!               if( lvand ) then
!                   call dcopy_a_to_b( hcsr, hcsrsv, nspnod )
!                   call svslmi( nfile, myid, nodes, nspin, lvandi, .false. )
!                   call svslmi( nfile, myid, nodes, 2, lvandi, .false. )
!               end if
           end if
       end if
       if( .not.lduplicate ) then
#ifndef VECTOR
           call exrhcr( glocal, hlocal, ntotfd_c )
#else
           call exrhcr( glocal, hlocal, kfft0 )
#endif
       else
#ifndef VECTOR
           call dcopy_a_to_b( glocal, hlocal, ntotfd_c )
#else
           call dcopy_a_to_b( glocal, hlocal, kfft0 )
#endif
       end if
   end if
   end if dupif2
end do
ct = timecnt()
if( ltimecnt ) then
    if(loutfile(1)) write(nfile(1),*) '                               total',  &
&                         ' : cpu-time :', ct-ct0
    if(loutfile(2)) write(nfile(2),*) '                               total',  &
&                         ' : cpu-time :', ct-ct0
!              write(nfile(1),*) '                                    ',
!     &                          ' : com-time :', t_comm
!              write(nfile(2),*) '                                    ',
!     &                          ' : com-time :', t_comm
!            else
!              write(nfile(1),*) ' Eigenvalue problem            ',
!     &                          ' : cpu-time :', ct-ct0
!              write(nfile(2),*) ' Eigenvalue problem            ',
!     &                          ' : cpu-time :', ct-ct0
end if
ct0 = ct

!---unify bfzansa1, bfzansa2 across k-point decomposition domains, if necessary.
call unify_sum1( bfzansa2 )
call unify_max1( bfzansa1 )
bfzansa2= abs(bfzansa2)  ! bfzansa2 could possibly be negative, if occupancies are negative.
bfzansa2 = bfzansa2/dble(nel)


!---set lstop in schmidt.f90
call setlstop( nfile, myid, nodes, .true. )
!-----------------------------------------------------------------------
      !==== stop here =========================
!      call fnormalstop( nfile, myid, nodes, 'pwlda3' )
      !========================================


!--- set occupancies
if( .not.ltddft_fssh .or. lfssh_gsscf ) then

    if( lspin ) then
        call fermi_k_spin( nfile, myid, nodes,  &
& nband, nkpnt, noband, nel, egvlud, nsordr, wegud, wbzk, nocc,  &
& eig, norder, occ,  &
& entrpy, feneud, fermie, lfermi, tfermi, diffud, lfixud )
    else
        call fermi_k( nfile, myid, nodes,  &
& nband, nkpnt, noband, nel, eig, norder, occ, wbzk, nocc,  &
& entrpy, fermie, lfermi, tfermi )
    end if

    !---tetrahedron method for BZ-zone sampling
!    if( ltetra ) then
!        if( lspin ) then
!            call BZtetrahedron_spin( nfile, myid, nodes,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, nkpnt, noband,  &
!& nel, egvlud, nsordr, wegud, wbzk, eig, norder, occ,  &
!& entrpy, feneud, fermie, diffud, lfixud, ekib, prod, prodr )
!        else
!            call BZtetrahedron( nfile, myid, nodes,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, nkpnt, noband,  &
!& nel, eig, norder, occ, wbzk, entrpy, fermie, ekib, prod, prodr, lnoncollinear )
!        end if
!    end if

end if

if( ltddft_fssh .and. .not.lfssh_gsscf ) then
    !---TDDFT-FSSH  fractional occupations if necessary.
    if( lspin ) then
        call broadening_in_tddft_fssh( nfile, myid, nodes,  &
nband, nbnod, nspnmx, egvlud, nsordr, wegud, entrpy )
      else
        call broadening_in_tddft_fssh( nfile, myid, nodes,  &
nband, nbnod, nspnmx, eig, norder, occ, entrpy )
    end if
end if


!--- output eigenvalues and occupancies
!if( ifmd.eq.0 ) then
if( louteig ) then
    if( lspin ) then
        call out_eigocc_spin( nfile, myid, nodes,  &
& nband, nkpnt, egvlud, nsordr, wegud, bzk, wbzk,  &
& nstep, jgcycl, feneud, 3 )
    else
        call out_eigocc( nfile, myid, nodes,  &
& nband, nkpnt, eig, norder, occ, bzk, wbzk, nstep, jgcycl, fermie, 3 )
    end if
end if
!end if
!-----------------------------------------------------------------------


call gsync
ct0 = timecnt()
csh0 = ct0
!if( laspc_exec ) then
!    !---No iteration is enforced.  if( zansa2 > tolres .and. jcycl < iscfmx ) cycle outer
!    laspc_conv = jcycl >= iscfmx_aspc
!
!    !---if( jcycl == 1 .and. aspc_corr < 0.99999d0 ) then
!!    if( aspc_corr < 0.99999d0 ) then
!    if( aspc_corr < 0.99999d0 .and. .false. ) then !---not used
!    !-----Wavefunction correction
!      if( ltimecnt ) then
!          if(loutfile(1)) write(nfile(1),*) ' Corrector (w.f.) in the ASPC method'
!          if(loutfile(2)) write(nfile(2),*) ' Corrector (w.f.) in the ASPC method'
!      end if
!      if( lgamma ) then
!          call corr_wf_aspc( nfile, myid, nodes,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& gdcr, pgdcr, nplw, nplwex, npnod1, npnod2, npnod, nplcnt, npldsp,  &
!& rhcr, cgjr, iods, iodsg, ioa, ioag, idstnd, bufcr, prod, prodr, pdbuf, &
!& dmtrxr, bijr, bm1r, bx0r, nbxxxx, w1r, w1i, ew1, tmpr, iblock, &
!& nprvmx, lspin, gdcrsv, pgdcrsv, nspnod, lcgjsv, nspnmx, &
!& ntype, nhk1, nhk2, natom, nhk1_nat, nhk2_nat,  &
!& natnod1, natnod2, natnod, natcnt, natdsp,  &
!& lvand, lvandi, lking, lvand_r, lvand_g,  &
!& iatoit, mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, rvol, rdelv_c,  &
!& fft3x, fft3y, fftwork, ijkg, nplw2, eigr, eigi, ycos, ysin, nplwcs,  &
!& node_c, node_r, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd,  &
!& .false. )
!        else
!              call sbalin_k( nfile, myid, nodes,  &
!     & nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!     & gdcr, nplw, nplwex, npnod1, npnod2, npnod, nplcnt, npldsp,  &
!     & iods, iodsg, idstnd, bufcr, keypwv, prveig, nprvmx, nspnmx, prcd, gnk,  &
!     & ihest, nstep, nstep_ini, prod, prodr, pdbuf, facc1, facc2, facc3,  &
!     & dmtrxr, bijr, bm1r, bx0r, nbxxxx, w1r, w1i, ew1, tmpr, iblock,  &
!     & node_c, node_r, lbncnt, lbndsp, mbncnt, mbndsp, jdstnd,  &
!     & t_comm, lspin, gdcrsv, nspnod, lcgjsv, rhcr, cgjr,  &
!     & lvand, ntype, nhk1, nhk2, natom, lvandi, lking, lvand_r, lvand_g,  &
!     & iatoit, mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, rvol, rdelv_c,  &
!     & fft3x, fft3y, fftwork, ijkg, nplw2, eigr, ycos, ysin, nplwcs,  &
!     & nknod, mshglb_c, kfft1, kfft2, kfft3, ratm, eigi,  &
!     & prods, dmtrxi, biji, bm1i, bx0i, .false. )
!      end if
!
!      ct = timecnt()
!!          if( ltimecnt ) then
!      if(loutfile(1)) write(nfile(1),*) ' Corrector (w.f.) in the ASPC method', &
!&                          ' : cpu-time :', ct-ct0
!      if(loutfile(2)) write(nfile(2),*) ' Corrector (w.f.) in the ASPC method', &
!&                          ' : cpu-time :', ct-ct0
!!          end if
!    end if
!end if


end if scfif
lfssh_updtcg = lconverge
if( lfssh_updtcg ) then
    !--- output eigenvalues and occupancies
    if( lspin ) then
        call out_eigocc_spin( nfile, myid, nodes, &
& nband, nkpnt, egvlud, nsordr, wegud, bzk, wbzk, &
& nstep, jgcycl, feneud, 3 )
      else
        call out_eigocc( nfile, myid, nodes, &
& nband, nkpnt, eig, norder, occ, bzk, wbzk, nstep, jgcycl, fermie, 3 )
    end if

    !----SCF with the ground state in TDDFT-FSSH
    if( lspin ) then
        call tddft_fssh_exchg_occ( nfile, myid, nodes,  &
& nband, nbnod, nspnmx, egvlud, nsordr, wegud )
      else
        call tddft_fssh_exchg_occ( nfile, myid, nodes,  &
& nband, nbnod, nspnmx, eig, norder, occ )
    end if

    !---TDDFT-FSSH  fractional occupations if necessary.
    if( lspin ) then
        call broadening_in_tddft_fssh( nfile, myid, nodes,  &
nband, nbnod, nspnmx, egvlud, nsordr, wegud, entrpy )
      else
        call broadening_in_tddft_fssh( nfile, myid, nodes,  &
nband, nbnod, nspnmx, eig, norder, occ, entrpy )
    end if

    !--- output eigenvalues and excited-occupancies
    if( lspin ) then
        call out_eigocc_spin( nfile, myid, nodes, &
& nband, nkpnt, egvlud, nsordr, wegud, bzk, wbzk, &
& nstep, jgcycl, feneud, 20 )
      else
        call out_eigocc( nfile, myid, nodes, &
& nband, nkpnt, eig, norder, occ, bzk, wbzk, nstep, jgcycl, fermie, 20 )
    end if

    !---ground state forces are used in FSSH
    if( lfssh_gsfrc ) then
        !----SCF with the ground state in TDDFT-FSSH
        if( lspin ) then
            call tddft_fssh_restore_occ( nfile, myid, nodes,  &
& nband, nbnod, nspnmx, egvlud, nsordr, wegud )
        else
            call tddft_fssh_restore_occ( nfile, myid, nodes,  &
& nband, nbnod, nspnmx, eig, norder, occ )
        end if
        exit outer
    end if

    lreset   = .true.
    lresetud = .true.
    call tddft_fssh_store_rhgr( nfile, myid, nodes,  &
& rhgr, nplw5ex, nspnmx )

    !----change charge mixing method
    call imxchgg_change( nfile, myid, nodes, imxchg_fssh )
!    if( lspin ) call imxchggud_change( nfile, myid, nodes, 0 ) ! do not mix

!    call tddft_fssh_store_rhoij( nfile, myid, nodes )

!    if( lrtddft ) then
!        call lrtddft_store_rho( nfile, myid, nodes,  &
!& rho, mshnod(1), 1 )
!        if( lspin ) then
!            call lrtddft_store_rho( nfile, myid, nodes,  &
!& rhoud, mshnod(1), 2 )
!        end if
!        call lrtddft_fssh_store_sg1( nfile, myid, nodes )
!    end if

    !---store Hartree potential in vhar_out
    vhar_out(1:mshnod(1)) = vhar(1:mshnod(1))
      delrho(1:mshnod(1)) =  rho(1:mshnod(1))

end if

if( ltimecnt ) then
    if(loutfile(1)) write(nfile(1),*) ' New electron density '
    if(loutfile(2)) write(nfile(2),*) ' New electron density '
end if
!--- output charge density ---------------------------------------------
call pwpchg( nfile, myid, nodes, ltimecnt, lfssh_updtcg )

!if( lnoncollinear .and. ltimecnt ) then
!    call outmagnemom( nfile, myid, nodes,  &
!& nstep, jgcycl, rdelv, mshnod(1), ntype, nhk, natom )
!end if

ct = timecnt()
if( ltimecnt ) then
    if(loutfile(1)) write(nfile(1),*) '                              total ',  &
&                         ' : cpu-time :', ct-ct0
    if(loutfile(2)) write(nfile(2),*) '                              total ',  &
&                         ' : cpu-time :', ct-ct0
!            else
!              write(nfile(1),*) ' New electron density          ',
!     &                          ' : cpu-time :', ct-ct0
!              write(nfile(2),*) ' New electron density          ',
!     &                          ' : cpu-time :', ct-ct0
end if
ct0 = ct


!--- charge density for next iteration ---------------------------------
!---if( laspc_exec ) then
call set_mixprm_aspc( nfile, myid, nodes, &
& ncprvmx, aslh_aspc, bslh_aspc, aslh__, bslh__, laspc_exec )
!---end if
if( jgcycl >= itrial .or. imxchg_trial /= 0 .or. imxspin_trial /= 0 ) then
    !-----Check symmetry of charge density
!    if( .not.lnoncollinear ) then
!        if( myid == 0 ) then
!            call csymcd( nfile, myid, nodes, ltimecnt,  &
!& rhgr, nd1vks, nplw5, nplw5ex, kfft0d, nga, ngb, ngc, ijkgd, 1 )
!        end if
!    else
!        !-----noncollinear magnetism
!        call csymcd( nfile, myid, nodes, ltimecnt,  &
!& rhgr, nd1vks, nplw5, nplw5ex, kfft0d, nga, ngb, ngc, ijkgd, 3 )
!    end if

    if( .not.lfssh_updtcg ) then
        if( .not.laspc_exec ) then
            aslh_ = aslh
            bslh_ = bslh
            amxspin_ = amxspin
            bmxspin_ = bmxspin
          else
            aslh_ = aslh__
            bslh_ = bslh__
            amxspin_ = amxspin
            bmxspin_ = bmxspin
        end if
      else
        aslh_ = aslh_fssh
        bslh_ = bslh_fssh
        amxspin_ = amxspin
        bmxspin_ = bmxspin
    end if
    if( jgcycl < itrial ) then
        if( imxchg_trial /= 0 ) then
            aslh_ = aslh_trial
            bslh_ = bslh_trial
            call imxchgg_change( nfile, myid, nodes, imxchg_trial )
        end if
        if( imxspin_trial /= 0 ) then
            amxspin_ = amxspin_trial
            bmxspin_ = bmxspin_trial
            if( lspin ) then
                call imxchggud_change( nfile, myid, nodes, imxspin_trial )
            else
                !-----noncollinear magnetism
!                call imxchggmag_change( nfile, myid, nodes, imxspin_trial )
            end if
        end if
    end if
    if( jgcycl == itrial .and. itrial > 0 ) then
        if( imxchg_trial /= 0 ) then
            lreset = .true.
            if(loutfile(1)) write(nfile(1),'(a,f7.4)') ' *** end of trial step, lreset = .true.'
            if(loutfile(2)) write(nfile(2),'(a,f7.4)') ' *** end of trial step, lreset = .true.'
        else
            if(loutfile(1)) write(nfile(1),'(a,f7.4)') ' *** end of trial step'
            if(loutfile(2)) write(nfile(2),'(a,f7.4)') ' *** end of trial step'
        end if
    end if

    !-----Charge mixing
    call pwpmxchgg( nfile, myid, nodes, ct0,  &
& lreset, aslh_, bslh_, amxspin_, bmxspin_, lmixhold, ltimecnt, lresetud )

      paw_mix =   paw_mix0
    plusU_mix = plusU_mix0
    if( jgcycl < itrial ) then
        if( imxchg_trial /= 0 ) call imxchgg_change( nfile, myid, nodes, imxchg )
        if( imxspin_trial /= 0 ) then
            if( lspin ) then
                call imxchggud_change( nfile, myid, nodes, imxspin )
            else
                !-----noncollinear magnetism
!                call imxchggmag_change( nfile, myid, nodes, imxspin )
            end if
        end if
    end if
  else
!--- return rhgr to input charge density
    call pwpinichg( nfile, myid, nodes )
    !---to keep rhoij for the PAW
    !---lonsite = .true.
      paw_mix = 0.d0
    plusU_mix = 0.d0
end if

!--- check No. of electrons
call chknel( nfile, myid, nodes,  &
& rho, mshnod(1), rdelv, nel, noddatx, ltimecnt )
!      if( lspin ) then
!          call chknud( nfile, myid, nodes,  &
!     & rhoud, mshnod(1), rdelv, diffud, lfixud, noddatx, ltimecnt )
!      end if

if( lfssh_updtcg ) then
    !----restore charge mixing method
    call imxchgg_change( nfile, myid, nodes, imxchg )
    if( lspin ) call imxchggud_change( nfile, myid, nodes, imxspin )
end if

end do outer
!---just for the safe
  paw_mix =   paw_mix0
plusU_mix = plusU_mix0
!=======================================================================
!----------- loop ended ------------------------------------------------
!=======================================================================
!--- output eigenvalues and occupancies
!          write(nfile(3),'(10i7)') nstep, jgcycl, nband
if( .not.lfssh_gsscf ) then
if( lspin ) then
    call out_eigocc_spin( nfile, myid, nodes,  &
& nband, nkpnt, egvlud, nsordr, wegud, bzk, wbzk,  &
& nstep, jgcycl, feneud, 3 )
  else
    call out_eigocc( nfile, myid, nodes,  &
& nband, nkpnt, eig, norder, occ, bzk, wbzk, nstep, jgcycl, fermie, 3 )
end if
end if


!--- output E-k dispersion relation
!if( lspin ) then
!    call out_dispersion( nfile, myid, nodes,  &
!& egvlud, nsordr, 2, eig, norder, nstep, fermie )
!else
!    call out_dispersion( nfile, myid, nodes,  &
!& eig, norder, 1, eig, norder, nstep, fermie )
!end if


!--- first-order SO effects by purturbation calculation
!call SOp( nfile, myid, nodes, ltimecnt, sume0 )


!if( lnoncollinear ) then
!    call outmagnemom( nfile, myid, nodes,  &
!& nstep, jgcycl, rdelv, mshnod(1), ntype, nhk, natom )
!end if

!      !==== stop here =========================
!      call fstop( nfile, myid, nodes, 'pwlda4' )
!      !========================================

call gsync
ct0 = timecnt()
do nspin = 1, nspnmx
   if( lspin ) then
       if( lgamma ) then
           call stspud( nspin, gdcr, gdcrsv, nspnod, lcgjsv )
       end if
   end if

kvecdo8: do kvec = 1, nknod
!   if( .not.lgamma ) then
!       call set_kvec( kvec )
!       call get_wfk( gdcr, nspin )
!       call get_plwk( nplw, nplwex, nspnod )
!       call get_indexk(  &
!& npnod1, npnod2, npnod, nplcnt, npldsp, nodes )
!   end if

leig_end   = kvec == nknod .and. nspin == nspnmx

!   --- copy gdcr to rhcr
call cpgdrh( nfile, myid, nodes, rhcr, gdcr, npnod, nband )
!   --- to convert G decomposition to band decomposition ---
call set_all2all_out( leig_end )
if( ltimecnt ) ctt0 = timecnt()
call gdtobd( nfile, myid, nodes, ctt0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& rhcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iods, iodsg, idstnd, ltimecnt )

!if( .not.lgamma ) then
!    !-----store w.f. decomposed by bands
!    call set_wfk( rhcr, nspin )
!end if

end do kvecdo8

   if( lspin ) then
       if( lgamma ) then
           call stpsud( nspin, rhcr, gdcrsv, nspnod, lcgjsv )
           if( nspin.eq.1 ) then
!           --- copy rhcr to gdcr
               call cpgdrh( nfile, myid, nodes,  &
& gdcr, rhcr, npnod, nband )
           end if
       end if
   end if
end do


!------- linear-response TDDFT --------------------------------------
!if( lrtddft ) then
!    if( lfssh_gsscf ) then
!        !---restore variables, when lfssh_gsscf = .true.
!        if( lspin ) then
!            call tddft_fssh_exchange_occ( nfile, myid, nodes,  &
!& nband, nbnod, nspnmx, egvlud, wegud )
!          else
!            call tddft_fssh_exchange_occ( nfile, myid, nodes,  &
!& nband, nbnod, nspnmx, eig, occ )
!        end if
!        if( .not.ltddft_nscforce ) then
!            call lrtddft_exchange_rho( nfile, myid, nodes,  &
!& rho, mshnod(1), 1 )
!            if( lspin ) then
!                call lrtddft_exchange_rho( nfile, myid, nodes,  &
!& rhoud, mshnod(1), 2 )
!            end if
!            call lrtddft_exchange_sg1( nfile, myid, nodes )
!        end if
!    end if
!
!    call lr_tddft( nfile, myid, nodes, ltimecnt )
!
!    if( lfssh_gsscf ) then
!        !---restore variables, when lfssh_gsscf = .true.
!        if( lspin ) then
!            call tddft_fssh_exchange_occ( nfile, myid, nodes,  &
!& nband, nbnod, nspnmx, egvlud, wegud )
!          else
!            call tddft_fssh_exchange_occ( nfile, myid, nodes,  &
!& nband, nbnod, nspnmx, eig, occ )
!        end if
!        if( .not.ltddft_nscforce ) then
!            call lrtddft_exchange_rho( nfile, myid, nodes,  &
!& rho, mshnod(1), 1 )
!            if( lspin ) then
!                call lrtddft_exchange_rho( nfile, myid, nodes,  &
!& rhoud, mshnod(1), 2 )
!            end if
!            call lrtddft_exchange_sg1( nfile, myid, nodes )
!        end if
!    end if
!end if

!if( lhoteh ) then
!if( mod(nstep,nskip_hoteh).eq.0 .or. ifmd.eq.0 ) then
!    !--- hot electron/hole relaxation rate calculation
!    call hot_eh_relaxation_rate( nfile, myid, nodes, ltimecnt )
!end if
!end if


if(loutfile(1)) write(nfile(1),*) ' '
if(loutfile(2)) write(nfile(2),*) ' '
if( ltimecnt ) then
    if(loutfile(1)) write(nfile(1),*) ' Force calculation'
    if(loutfile(2)) write(nfile(2),*) ' Force calculation'
end if
!--- force and internal stress tensor ----------------------------------
t_comm = 0.d0
call force( nfile, myid, nodes,  &
& t_comm, lcstress, ltimecnt, enclpv, elcl )

!---HSE hybrid functional: short-ranged HF force
!call HF_SR_force( nfile, myid, nodes, lcstress, ltimecnt )

!---uniform electric field for insulator
!call periodic_efield_force( nfile, myid, nodes, lcstress, ltimecnt )


!      if( lshdep ) then
!          ecshdp = 0.d0
!          call shdfrc( nfile, myid, nodes,  &
!     & lvacuum, cshdep, nel, ntype, natom, nhk1, nhk2, zv,  &
!     & dpion, dpion2, dprho, ecshdp, frc )
!          ecolmb = ecolmb + ecshdp
!!--- etshdp is shape dependent energy calculated from sum of potential
!!--- etshd2 is shape dependent energy calculated from total dipole
!          etshdp = elshdp + ehshdp + ecshdp
!          etshd2 = cshdep*( ( dpion(1) - dprho(1) )**2
!     &                    + ( dpion(2) - dprho(2) )**2
!     &                    + ( dpion(3) - dprho(3) )**2 )
!      end if

sume = sume + ecolmb
etot = etot + ecolmb

!-----energy from an empirical correction for vdW : DFT-D
edisp = 0.d0
if( ldftd ) call add_energy_dftd( nfile, myid, nodes, sume, etot, edisp )

!-----ionic energy from uniform electric field
eefield_ion = 0.d0
!if( lefield ) call add_energy_efield( nfile, myid, nodes,  &
!& sume, etot, eefield_ion, frc, fefd, natom, strefd )

!-----set potential energy in TDDFT-FSSH
if( ltddft_fssh ) call set_pote_in_tddft_fssh( sume )

!--- correction to force
!if( .not.lefield ) call corfrc( ntype, nhk1, nhk2, frc, natom )
call corfrc( ntype, nhk1, nhk2, frc, natom )


if( .not.ltddft_fssh .or. .not.ltransition ) exit fsshcycle
!---TDDFT-FSSH: check acceptance and scale velocities if necessary
if( lspin ) then
    call tddft_fssh_judge_accept( nfile, myid, nodes,  &
& ltransition, nband, nspnmx, wegud, sume, &
& vatm, frc, watom, ntype, natom, nhk1, nhk2, dtmd )
  else
    call tddft_fssh_judge_accept( nfile, myid, nodes,  &
& ltransition, nband, nspnmx, occ, sume, &
& vatm, frc, watom, ntype, natom, nhk1, nhk2, dtmd )
end if
if( ltransition ) exit fsshcycle

ct = timecnt()
ct0 = ct
end do fsshcycle
!-----check if velocities need to be sent back
call set_lvsendback( nfile, myid, nodes, ltransition )


if( ltimecnt ) then

    do ii = 1, 2
    if( loutfile(ii) ) then
    write(nfile(ii),*) '                           total ',  &
&                      ': cpu-time :', ct- ct0
!          write(nfile(ii),*) '                                 ',
!     &                       ': com-time :', t_comm
    write(nfile(ii),*) ' '
    write(nfile(ii),1120) ecolmb, ecolmb*elvlt
!-----energy from an empirical correction for vdW : DFT-D
if( ldftd ) call out_energy_dftd( nfile, myid, nodes, elvlt, ii )
    write(nfile(ii),trim(fmtHF)) sume, sume*elvlt
    if( lfermi >= 3 .and. .not.ltetra ) then
        !---the oder of approximant
        npolyorder = lfermi - 3
        !---extrapolated energy to tfermi -> 0
        eextrpltd = ( dble(npolyorder+1)*sume + sume + entrpy )/dble(npolyorder+2)
        write(nfile(ii),1111) eextrpltd, eextrpltd*elvlt
    end if
!          if( lshdep ) then
!              write(nfile(ii),1130) etshdp, etshdp*elvlt
!              write(nfile(ii),1131) etshd2, etshd2*elvlt
!              write(nfile(ii),*) ' * Dipole moment'
!              write(nfile(ii),'((6x,a10,3e19.10))')
!     &         'atomic    ',(dpion(i), i = 1, 3),
!     &         'electronic',(-dprho(i), i = 1, 3),
!     &         'total     ',(dpion(i)-dprho(i), i = 1, 3)
!          end if
1120 format('  * Coulomb E.        = ',F15.6,' ( Ryd. )',F13.4,' (eV)')
1130 format('  * Shape dep. E.     = ',F15.6,' ( Ryd. )',F13.4,' (eV)')
1131 format('             ( check  = ',F15.6,' ( Ryd. )',F13.4,' (eV)',  &
& ' )')
1111 format('  * Extrapolated E.   = ',F15.6,' ( Ryd. )',F13.4,' (eV)')
        write(nfile(ii),*) ' * Force'
        do i = 1, nhk2(ntype)
!                 write(nfile(ii),'(i6,a12,3e19.10)') &
!     &                      i, ' (direct)  :', ( fclm(ix,i), ix = 1, 3)
!                 write(nfile(ii),'(6x,a12,3e19.10)') &
!     &                         ' (local)   :', ( floc(ix,i), ix = 1, 3)
!                 write(nfile(ii),'(6x,a12,3e19.10)') &
!     &                         ' (nonlocal):', ( fnlc(ix,i), ix = 1, 3)
!                 write(nfile(ii),'(6x,a12,3e19.10)') &
!     &                         ' (pcc)     :', ( fpcc(ix,i), ix = 1, 3)
!                 !-----energy from an empirical correction for vdW : DFT-D
!                 if( ldftd ) call out_force_dftd( nfile, myid, nodes, ii, i )
!                 if( jhybrid /= 0 ) write(nfile(ii),'(6x,a12,3e19.10)') &
!     &                         ' (HF-SR)   :', ( fhfs(ix,i), ix = 1, 3)
!                 if( lefield ) write(nfile(ii),'(6x,a12,3e19.10)') &
!     &                         ' (E field) :', ( fefd(ix,i), ix = 1, 3)
!                 write(nfile(ii),'(6x,a12,3e19.10)') &
!     &                         ' (total)   :', ( frc(ix,i), ix = 1, 3)
           write(nfile(ii),'(i6,3e19.10)')  &
&                                     i, ( frc(ix,i), ix = 1, 3)
        end do
    end if
    end do

  else

    do ii = 1, 2
    if( loutfile(ii) ) then
       write(nfile(ii),'(a27,f10.4)')  &
&                   ' Force calculation (sec) : ', ct- ct0
!-----energy from an empirical correction for vdW : DFT-D
if( ldftd ) call out_energy_dftd( nfile, myid, nodes, elvlt, ii )
       write(nfile(ii),trim(fmtHF)) sume, sume*elvlt
!             if( lshdep ) then
!                 write(nfile(ii),1130) etshdp, etshdp*elvlt
!                 write(nfile(ii),*) ' * Dipole moment'
!                 write(nfile(ii),'((6x,a10,3es15.6))')
!     &         'atomic    ',(dpion(i), i = 1, 3),
!     &         'electronic',(-dprho(i), i = 1, 3),
!     &         'total     ',(dpion(i)-dprho(i), i = 1, 3)
!             end if
    end if
    end do

end if

if( loutfile(1) ) then
!    if( lpaw ) then
!        eonsite = EKIN1 + EHAR1 + EEXT1 + EEXCE1 + ECORE1
!      else
        eonsite = 0.d0
!    end if
    write(nfile(5),1150) nstep, jgcycl, sume, etot,  &
&   ekin, eext, ehar, eexc, ecor-ecnl, ecnl, entrpy, eonsite, exhfsr, EplusU, eefield, soenergy,  &
&   ecolmb, edisp, eefield_ion
end if

!if( lefield .and. lsawtooth ) then
!    !--- electronic polarization with sawtooth electric field
!    call sawtooth_efield_pol( nfile, myid, nodes,  &
!& lsawtooth_shape, lsawtooth_xyz, rho,  &
!& lclust, rdel, rdelg, hcell, lorthrhmbc, nd1v, nd1vks, rdelv,  &
!& mshxyz(mulpit(1)),  &
!& mulx1(1), mulx2(1), muly1(1), muly2(1), mulz1(1), mulz2(1),  &
!& mshnx(mulpms(1)), mshny(mulpms(1)), mshnz(mulpms(1)), mshnod(1),  &
!& mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1) )
!end if
if( loutfile(1) ) then
    write(nfile(29),'(2i6,6es12.4)') nstep, jgcycl, difene, difene2, zansa1, zansa2,  &
& bfzansa1, bfzansa2

!    if( lefield ) then
!        !---uniform electric field for insulator
!        call out_polarization( nfile, myid, nodes, .true., nstep, jgcycl )
!    end if
end if

!      !==== stop here =========================
!      call fstop( nfile, myid, nodes, 'pwlda5' )
!      !========================================


ct0 = ct
if( lcstress ) then
!--- internal stress tensor from kinetic E. & Hartree E. & Exc ---
strexc = ( eexc - evexc ) + ( ecor + ecnl - evcor )
call stress( nfile, myid, nodes,  &
& strtot, strloc, strnlc, strclm, strgga, strpcc,  &
& ltimecnt, nstep,  &
& rhcr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, occ, rvol,  &
& iods,  &
& rhgr, nplw5ex, nplw5, recnrm, gx, gy, gz,  &
& ehar, vexc, strexc,  &
& rdel, rdelg, rdelv, hcell, hci, lorthrhmbc, lvacuum,  &
& ntype, natom, nhk1, nhk2, ratm, mx1, enclp, t_comm,  &
& lpcc, lpcc_r, lpcc_g, lpcci, rpcc, lpking,  &
&   mshxyz(mulpit(1)),  &
&   mulx1(1), mulx2(1), muly1(1), muly2(1), mulz1(1), mulz2(1),  &
&   mshnx(mulpms(1)), mshny(mulpms(1)), mshnz(mulpms(1)),  &
&   mshnod(1), mshnew(1),  &
&   mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1),  &
& lspin, nspnmx, gdcrsv, nspnod, lcgjsv, wegud, ldftd,  &
& strhfs, exhfsr, jhybrid, ncscale, strefd, lefield )
ct = timecnt()
do ii = 1, 2
   if( loutfile(ii) ) then
!             write(nfile(ii),*) ' '
       write(nfile(ii),'(a27,f10.4)')  &
&                   ' Stress calcu.     (sec) : ', ct- ct0
!             write(nfile(ii),*) '                           ',
!     &                          ' : com-time :', t_comm
   end if
end do
ct0 = ct
!-----------------------------------------------------------------------
end if

!      !==== stop here =========================
!      call fstop( nfile, myid, nodes, 'pwlda6' )
!      !========================================

!if( lintchg .and. mod(nstep,nskip_intchg) == 0 ) then
!!-----------------------------------------------------------------------
!!    calculate atomic charge by integrating charge density
!call calachg( nfile, myid, nodes,  &
!& nstep, lvacuum, lclust, ntype, natom, nhk1, nhk2, rintchg, ratm,  &
!& iatoit, rdel, rdelg, rdelv, hcell, hci, rho, bufatm, bufatmr,  &
!& mshxyz(mulpit(1)),  &
!& mulx1(1), mulx2(1), muly1(1), muly2(1), mulz1(1), mulz2(1),  &
!& mshnx(mulpms(1)), mshny(mulpms(1)), mshnz(mulpms(1)), mshnod(1),  &
!& mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1) )
!ct = timecnt()
!do i = 1, 2
!if( loutfile(i) ) then
!    write(nfile(i),*) ' '
!    write(nfile(i),*) ' Atomic charge      : cpu-time :',ct-ct0
!end if
!end do
!ct0 = ct
!!-----------------------------------------------------------------------
!end if


!-----------------------------------------------------------------------
!    Spherical harmonics expansion
!if( lsphexp ) then
!if( mod(nstep,nskpsphexp).eq.0 .or. ifmd.eq.0 ) then
!    call sphexp( nfile, myid, nodes, &
!& rhcr, nplwex, nplw, nband, nbnod, iods, &
!& lspin, nspnmx, gdcrsv, nspnod, lcgjsv, &
!& recnrm, gx,  gy,  gz, ycos, ysin, nplwcs, nknod, bzk, wbzk, nstep, &
!& ntype, natom, nhk1, nhk2, iatoit, lvand, lvandi, lmax, lchk, mxl, rvol, rsphexp )
!end if
!end if
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!    Mulliken population analysis
!if( lmulken ) then
!if( mod(nstep,nskpmulk).eq.0 .or. ifmd.eq.0 ) then
!!--- structure factors
!    call mlstruc( nfile, myid, nodes,  &
!& lclust, ntype, nhk1, nhk2, ratm, hcell, nd1v, nd1vks,  &
!& nplw, npnod1, npnod2, npnod,  &
!& prcd(1), gnk(1), hnk(1), prcd(nplw+2), gnk(nplw+2), hnk(nplw+2),  &
!& nplw/2, ltimecnt )
!if( lvand_r ) then
!    !-----initial set for real-space calculation
!    call setvnlc( nfile, myid, nodes,  &
!& lvacuum, lclust, ntype, nhk1, nhk2, natom, ratm, lvandi, lking,  &
!& iatoit, lmax, rdel_c, rdelg_c, hcell, hci, lorthrhmbc,  &
!& mshglb_c, kfft1, kfft2, kfft3, ylmr, ylmi, mxl, alloc_mem )
!end if

!if( leda ) call eda_clear( nfile, myid, nodes,  &
!& ntype, nhk1, nhk2 )

!do nspin = 1, nspnmx
!   if( lspin ) then
!       call stspud( nspin, rhcr, gdcrsv, nspnod, lcgjsv )
!       call ldocc( nspin, occ, wegud, nband )
!!       if( lvand ) then
!!           call ldslmi( nfile, myid, nodes, nspin, lvandi, .false. )
!!       end if
!!    --- copy rhcr to cgjr ---
!       call cprhcg( nfile, myid, nodes, cgjr, rhcr, nplwex, nbnod)
!       if( ltimecnt ) csh0 = timecnt()
!!    --- convert band decomposition to G decomposition ---
!       call bdtogd( nfile, myid, nodes, csh0,  &
!& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
!& iodg, idstnd, .false., ltimecnt )
!       do i = 1, 2
!       if( loutfile(i) ) then
!           if( nspin.eq.1 ) then
!               write(nfile(i),*) '  ( for up-spin )'
!           else
!               write(nfile(i),*) '  ( for down-spin )'
!           end if
!       end if
!       end do
!       do ib = 1, nband
!          norder(ib) = nsordr(ib+nband*(nspin-1))
!       end do
!   end if
!
!   !--- copy gdcr to cgjr
!   call cpgdrh( nfile, myid, nodes, cgjr, gdcr, npnod, nband )
!
!   call mulikf( nfile, myid, nodes, nspin, ltimecnt )
!
   !-----if lvand_r == .true., rhcr is destroyed.
!   if( lvand_r ) then
!       !--- copy gdcr to rhcr
!       call cpgdrh( nfile, myid, nodes, rhcr, gdcr, npnod, nband )
!       if( ltimecnt ) csh0 = timecnt()
!       !--- to convert G decomposition to band decomposition
!       call gdtobd( nfile, myid, nodes, csh0,  &
!& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& rhcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
!& iods, iodsg, idstnd, ltimecnt )
!   end if

!end do
!
!    ct = timecnt()
!    do i = 1, 2
!    if( loutfile(i) ) then
!       write(nfile(i),*) ' '
!       write(nfile(i),*) ' Mulliken analysis  : cpu-time :',ct-ct0
!    end if
!    end do
!    ct0 = ct
!end if
!end if
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!    Energy Density Analysis (EDA)
!if( leda ) then
!if( mod(nstep,nskpeda).eq.0 .or. ifmd.eq.0 ) then
!
!    !--- Atomic onsite energy
!    call onsite_eda( nfile, myid, nodes,  &
!& natom )

!    !--- Atomic direct coulomb energy
!    call eda_coulomb( nfile, myid, nodes,  &
!& atom_ecolmb, natom )

    !--- EDA for local potential (2)
!    call eda_local2( nfile, myid, nodes,  &
!& llclpp_r, llclpp_g, ldouble_grid_recip, lclust, lorthrhmbc,  &
!& nd1v, rdelv, ntype, nhk, nhk1, nhk2, nhk1_nod, nhk2_nod,  &
!& zv, nel, ratm, iatoit, iatmpt_nod, nion_nod, natom,  &
!& llking, mx1loc, rho, noddatx, rhgr,  &
!& hcell, rdel, lhfull, tablc, tablca, dltlc, rmxlc,  &
!& mshxyz(mulpit(1)),  &
!& mulx1(1), mulx2(1), muly1(1), muly2(1), mulz1(1), mulz2(1),  &
!& mshnx(mulpms(1)), mshny(mulpms(1)), mshnz(mulpms(1)), mshnod(1),  &
!& mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1),  &
!& nplw5, nplw5ex, nplw, nga, ngb, ngc, rvol, & !& ycos, ysin, nplwcs,
!& DCQ1, DCQ2, DCQ3, DSQ1, DSQ2, DSQ3, kmax1cs, kmax2cs, kmax3cs,  &
!& tmpk, bufatm, bufatmr, ltimecnt )

!          !--- Delaunay tetrahedralization
!          call delaunay( nfile, myid, nodes,
!     & lclust, lvacuum, ntype, nhk1, nhk2, ratm, hcell, rdelaunay )

    !--- Grid-based EDA for local potential (1), Hartree,
    !--- and exchange-correlation energies
!    call gdeda( nfile, myid, nodes,  &
!& lclust, rdel, rdelg, hcell, hci, lorthrhmbc, lvacuum,  &
!& nd1v, nd1vks, rdelv,  &
!& ntype, nhk1, nhk2, iatoit, ratm, radgeda, radeda,  &
!& rho, rhocore, noddatx, eeexc, vhar, vhshdp, vext, ltimecnt,  &
!& mshxyz(mulpit(1)),  &
!& mulx1(1), mulx2(1), muly1(1), muly2(1), mulz1(1), mulz2(1),  &
!& mshnx(mulpms(1)), mshny(mulpms(1)), mshnz(mulpms(1)), mshnod(1),  &
!& mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1) )

    !--- output results
!    call eda_out( nfile, myid, nodes,  &
!& nstep, ntype, nhk1, nhk2, ltimecnt )

!    ct = timecnt()
!    do i = 1, 2
!    if( loutfile(i) ) then
!        write(nfile(i),*) ' '
!        write(nfile(i),*) ' EDA  : cpu-time :',ct-ct0
!    end if
!    end do
!    ct0 = ct
!end if
!end if
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!    Maximally localized Wannier functions
! rhcr : band-decomposed variables, will be distroyed
! cgjr :    G-decomposed variables, will be distroyed
!
!if( lwannier ) then
!if( mod(nstep,nskpwann).eq.0 .or. ifmd.eq.0 ) then
!
!if( lspin ) then
!    tel = 0.d0
!    do m1 = 1, mshnod(1)
!       tel = tel + rhoud( m1 )
!    end do
!    call gdsum(tel,1,dbuf1r)
!    tel = tel * rdelv
!end if
!do nspin = 1, nspnmx
!   if( ltimecnt ) csh0 = timecnt()
!   if( lspin ) then
!   !--- spin-polarized case
!       call stspud( nspin, rhcr, gdcrsv, nspnod, lcgjsv )
!       if( lvand ) then
!           call ldslmi( nfile, myid, nodes, nspin, lvandi, .false. )
!           !--- convert atom decomposition to band decomposition 
!           call slm_gdtobd( nfile, myid, nodes, csh0,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& natom, natnod1, natnod2, natnod, natcnt, natdsp,  &
!& iods, iodsg, ioag, idstnd, ltimecnt )
!       end if
!!    --- copy rhcr to cgjr ---
!       call cprhcg( nfile, myid, nodes, cgjr, rhcr, nplwex, nbnod)
!!    --- convert band decomposition to G decomposition ---
!       call bdtogd( nfile, myid, nodes, csh0,  &
!& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
!& iodg, idstnd, .false., ltimecnt )
!       do i = 1, 2
!       if( loutfile(i) ) then
!           if( nspin.eq.1 ) then
!               write(nfile(i),*) '  ( for up-spin )'
!           else
!               write(nfile(i),*) '  ( for down-spin )'
!           end if
!       end if
!       end do
!       if( nspin.eq.1 ) then
!           nbwann = nint( ( dble(nel) + tel )*0.5d0 )
!         else
!           nbwann = nint( ( dble(nel) - tel )*0.5d0 )
!       end if
!   else
!   !--- spin-unpolarized case
!       nbwann = nel/2
!       if( mod(nel,2).ne.0 ) then
!           nbwann = nbwann + 1
!           do i = 1, 2
!           if( loutfile(i) ) then
!               write(nfile(i),*) ' '
!               write(nfile(i),*)  &
!&                 ' *** Warning : No. of electrons is not even.'
!           end if
!           end do
!       end if
!       if( lvand ) then
!           !--- convert atom decomposition to band decomposition 
!           call slm_gdtobd( nfile, myid, nodes, csh0,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& natom, natnod1, natnod2, natnod, natcnt, natdsp,  &
!& iods, iodsg, ioag, idstnd, ltimecnt )
!       end if
!   end if
!
!   !--- unify slmir, slmii
!   if( lvand ) call unifyslm( nfile, myid, nodes,  &
!& nband, nbnod1, nbnod2, nbnod, iods, iodsg,  &
!& ntype, nhk1, nhk2, natom, lvandi, bijr, nbxxxx, .false. )
!
!!      --- nbwann is No. of occupied states
!   iwbndtmp1 = iwbnd1
!   iwbndtmp2 = iwbnd2
!   if( iwbndtmp1.eq.0 .and. iwbndtmp2.eq.0 ) then
!       iwbndtmp1 = 1
!       iwbndtmp2 = nbwann
!   end if
!   if( iwbndtmp1.le.0      ) iwbndtmp1 = 1
!   if( iwbndtmp2.gt.nbwann ) iwbndtmp2 = nbwann
!   nbwann = iwbndtmp2 - iwbndtmp1 + 1
!
!   iwstttmp1 = iwstt1
!   iwstttmp2 = iwstt2
!   if( iwstttmp1.eq.0 .and. iwstttmp2.eq.0 ) then
!       iwstttmp1 = 1
!       iwstttmp2 = nbwann
!   end if
!   if( iwstttmp1.le.0      ) iwstttmp1 = 1
!   if( iwstttmp2.gt.nbwann ) iwstttmp2 = nbwann
!
!!      --- copy gdcr to cgjr
!   call cpgdrh( nfile, myid, nodes, cgjr, gdcr, npnod, nband )
!
!   if( nbwann > 0 ) then
!
!       if( iwbndtmp1 /= 1 ) then
!!               --- compress cgjr
!           call compcg( nfile, myid, nodes, cgjr, npnod, nband,  &
!& iwbndtmp1, iwbndtmp2 )
!!               --- copy cgjr to rhcr
!          call cpgdrh( nfile, myid, nodes, rhcr, cgjr,  &
!& npnod, nband )
!          if( ltimecnt ) csh0 = timecnt()
!!               --- to convert G decomposition to band decomposition ---
!          call gdtobd( nfile, myid, nodes, csh0,  &
!& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& rhcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
!& iods, iodsg, idstnd, ltimecnt )
!
!          !--- set start band index in nlvand.f
!          call setiwbnd1( nfile, myid, nodes, iwbndtmp1 )
!       end if
!
!       call mwannier( nfile, myid, nodes,  &
!& nstep, lclust, lvand, lvacuum, ltimecnt, hcell, nd1v,  &
!& nbwann, tolwan, iterwan, dname_wann, iwstttmp1, iwstttmp2,  &
!& natwan, natwno, nwanaa, loutuni,  &
!& nspin, lspin, bijr, bm1r, dmtrxr, bx0r, ew1, w1r, w1i,  &
!& cgjr, rhcr, bufcr, prod, prodr, prods, apk, eigr, glocal,  &
!& nplwex, nplw, npnod1, npnod2, npnod, nplcnt, npldsp,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, iodsg, idstnd,  &
!& iods, mfd2ft, ntotfd, nd1vks, kfft0d, ijkgd, nplw5, rvol,  &
!& mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, ijkg, nplw2,  &
!& kfft1b, kfft2b, kfft3b, kfft0b, ijkgb, nplw3,  &
!& ngboxa, ngboxb, ngboxc,  &
!& ntype, nhk1, nhk2, natom, natnod, iatoit, ioa, lvandi, lking, ratm, wrka,  &
!& fft3x, fft3y, fftwork, tmpr, iblock, mshglb )
!
!     else
!
!        if(loutfile(1)) write(nfile(1),*) ' stop Wannier functions,',  &
!& ' because of no bands: ', iwbndtmp1, iwbndtmp2
!        if(loutfile(2)) write(nfile(2),*) ' stop Wannier functions,',  &
!& ' because of no bands: ', iwbndtmp1, iwbndtmp2
!
!   end if
!
!!   --- copy gdcr to rhcr
!   call cpgdrh( nfile, myid, nodes, rhcr, gdcr, npnod, nband )
!   if( ltimecnt ) csh0 = timecnt()
!!   --- to convert G decomposition to band decomposition ---
!   call gdtobd( nfile, myid, nodes, csh0,  &
!& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& rhcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
!& iods, iodsg, idstnd, ltimecnt )
!
!end do
!
!    ct = timecnt()
!    do i = 1, 2
!    if( loutfile(i) ) then
!       write(nfile(i),*) ' '
!       write(nfile(i),*) ' Wannier functions  : cpu-time :',ct-ct0
!    end if
!    end do
!    ct0 = ct
!end if
!end if
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!    Conductivity calculation
!if( lconduct ) then
!if( mod(nstep,nskpconduct).eq.0 .or. ifmd.eq.0 ) then
!
!if(loutfile(1)) write(nfile(1),*) ' '
!if(loutfile(2)) write(nfile(2),*) ' '
!if( ltimecnt ) then
!    if(loutfile(1)) write(nfile(1),*) ' Conductivity calculation'
!    if(loutfile(2)) write(nfile(2),*) ' Conductivity calculation'
!end if
!
!--- initialization for conductivity calculation
!call conductivity_clear( nfile, myid, nodes,  &
!& ldcconduct, kfft0d )
!
!do nspin = 1, nspnmx
!   if( lspin ) then
!       call ldocck( nspin, occ, wegud, nband, nkpnt )
!       call ldeig_k( nspin, eig, egvlud, nband, nkpnt )
!       fermie = feneud(nspin)
!       if( lgamma ) then
!           call stspud( nspin, rhcr, gdcrsv, nspnod, lcgjsv )
!           if( lvand ) then
!               call ldslmi( nfile, myid, nodes, nspin, lvandi, .false. )
!           end if
!
!!    --- copy rhcr to cgjr ---
!           call cprhcg( nfile, myid, nodes, cgjr, rhcr, nplwex, nbnod)
!           if( ltimecnt ) csh0 = timecnt()
!!    --- convert band decomposition to G decomposition ---
!           call bdtogd( nfile, myid, nodes, csh0,  &
!& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
!& iodg, idstnd, .false., ltimecnt )
!       end if
!       do i = 1, 2
!       if( loutfile(i) ) then
!           if( nspin.eq.1 ) then
!               write(nfile(i),*) '  ( for up-spin )'
!           else
!               write(nfile(i),*) '  ( for down-spin )'
!           end if
!       end if
!       end do
!   end if
!
!   kvecdo9: do kvec = 1, nknod
!      if( .not.lgamma ) then
!          call set_kvec( kvec )
!          call get_wfk( rhcr, nspin )
!          call get_plwk( nplw, nplwex, nspnod )
!          call get_indexk(  &
!& npnod1, npnod2, npnod, nplcnt, npldsp, nodes )
!          if( lvand ) then
!              call ldslmi_k( nfile, myid, nodes, lvandi, nspin, .false. )
!          end if
!
!!    --- copy rhcr to cgjr ---
!          call cprhcg( nfile, myid, nodes, cgjr, rhcr, nplwex, nbnod)
!          leig_end   = kvec == nknod .and. nspin == nspnmx
!!    --- convert band decomposition to G decomposition ---
!          call set_all2all_out( leig_end )
!          if( ltimecnt ) csh0 = timecnt()
!          call bdtogd( nfile, myid, nodes, csh0,  &
!& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
!& iodg, idstnd, .false., ltimecnt )
!      end if
!
!   !--- copy gdcr to cgjr
!   call cpgdrh( nfile, myid, nodes, cgjr, gdcr, npnod, nband )
!
!   call conductivity( nfile, myid, nodes,  &
!& ldcconduct, lacconduct, wgconduct, tempconduct, efconduct, freqacmx,  &
!& icband1, icband2, ichband1, ichband2, iceband1, iceband2, immtband1, immtband2,  &
!& lspin, nspin, kvec, cgjr, rhcr, nplwex, nplw, npnod1, npnod2, npnod,  &
!& gx, gy, gz, wbzk(kvec),  &
!& nband, nbnod1, nbnod2, nbnod, iods, eig, fermie,  &
!& mfd2ft, ntotfd, nd1vks, kfft0d, nd1vks_c, rvol,  &
!& fft3x, fft3y, fftwork, ijkgd, nplw5,  &
!& apk, eigr, eigi, prod, prodr,  &
!& kfft1b, kfft2b, kfft3b, kfft0b, ijkgb, nplw3,  &
!& ngboxa, ngboxb, ngboxc,  &
!& ntype, natom, iatoit, ioa, lvand, lpaw, lvandi, lchk, mxl )
!
!   end do kvecdo9
!
!end do
!
!!--- output data
!call conductivity_out( nfile, myid, nodes,  &
!& ldcconduct, lacconduct, nstep, nband, nd1vks, kfft0d, rvol, nspnmx,  &
!& nkpnt, nknod, wbzk, fft3x, fft3y )
!
!    ct = timecnt()
!    do i = 1, 2
!    if( loutfile(i) ) then
!        write(nfile(i),*) ' '
!        write(nfile(i),*) ' Conductivity calculation : cpu-time :',ct-ct0
!    end if
!    end do
!    ct0 = ct
!end if
!end if
!-----------------------------------------------------------------------


!--- output supercell, atomic coordinates and atomic forces
call wrgeom( nfile, myid, nodes,  &
& nstep, lclust, ratm, frc, wrka, nhk, ntype, natom,  &
& nd1v, nd1vks, hcell )


!-----------------------------------------------------------------------
!--- dump electron density
!if( ldpchg ) then
!if( mod(nstep,nskip_dpchg).eq.0 .or. ifmd.eq.0 ) then
!
!if( lpaw .and. ldpchg_paw ) then
!
!    do nspin = 1, nspnmx
!       if( lspin ) then
!           call stspud( nspin, rhcr, gdcrsv, nspnod, lcgjsv )
!           call ldocc( nspin, occ, wegud, nband )
!       end if
!
!       call outsoftchg( nfile, myid, nodes, ct0,  &
!& glocal, occ, rhcr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod,  &
!& iods, rvol, kfft0d,  &
!& mfd2ft_c, ntotfd_c, nd1vks_c, kfft0, ijkg, nplw2ex, nplw2,  &
!& fft3x, fft3y, fftwork, rdelv, rhgrsoft(nplw5ex*(nspin-1)+1),  &
!& nplw5ex, ltimecnt )
!
!    end do
!
!    if( lspin ) then
!        do ig = 1, nplw5ex
!           densup = rhgrsoft(ig)
!           densdn = rhgrsoft(nplw5ex+ig)
!           rhgrsoft(ig)         = densup + densdn
!           rhgrsoft(nplw5ex+ig) = densup - densdn
!        end do
!    end if
!
!    call cddump( nfile, myid, nodes,  &
!& dname_cds, nstep, rhgrsoft, nplw5ex, nplw5, glocal,  &
!& mfd2ft, ntotfd, nd1vks, kfft0d, rvol,  &
!& fft3x, fft3y, fftwork, ijkgd, nplw5, .true.,  &
!& dc_xmin, dc_xmax, dc_ymin, dc_ymax, dc_zmin, dc_zmax,  &
!& lcompress_dpchg, ndigit_dpchg )
!
!    if( lspin ) then
!        call cddump( nfile, myid, nodes,  &
!& dname_cds, nstep, rhgrsoft(nplw5ex+1), nplw5ex, nplw5, glocal,  &
!& mfd2ft, ntotfd, nd1vks, kfft0d, rvol,  &
!& fft3x, fft3y, fftwork, ijkgd, nplw5, .false.,  &
!& dc_xmin, dc_xmax, dc_ymin, dc_ymax, dc_zmin, dc_zmax,  &
!& lcompress_dpchg, ndigit_dpchg )
!    end if
!
!  else
!
!    call cddump( nfile, myid, nodes,  &
!& dname_cds, nstep, rhgr, nplw5ex, nplw5, glocal,  &
!& mfd2ft, ntotfd, nd1vks, kfft0d, rvol,  &
!& fft3x, fft3y, fftwork, ijkgd, nplw5, .true.,  &
!& dc_xmin, dc_xmax, dc_ymin, dc_ymax, dc_zmin, dc_zmax,  &
!& lcompress_dpchg, ndigit_dpchg )
!
!    if( lspin ) then
!        call cddump( nfile, myid, nodes,  &
!& dname_cds, nstep, rhgr(nplw5ex+1), nplw5ex, nplw5, glocal,  &
!& mfd2ft, ntotfd, nd1vks, kfft0d, rvol,  &
!& fft3x, fft3y, fftwork, ijkgd, nplw5, .false.,  &
!& dc_xmin, dc_xmax, dc_ymin, dc_ymax, dc_zmin, dc_zmax,  &
!& lcompress_dpchg, ndigit_dpchg )
!    end if
!
!end if
!
!    ct = timecnt()
!    if(loutfile(1)) write(nfile(1),*) ' dump charge density : cpu-time :', ct- ct0
!    if(loutfile(2)) write(nfile(2),*) ' dump charge density : cpu-time :', ct- ct0
!    ct0 = ct
!
!end if
!end if
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!--- dump wavefunctions
if( ldpwav ) then
if( mod(nstep,nskip_dpwav).eq.0 .or. ifmd.eq.0 ) then
    do nspin = 1, nspnmx
       if( lspin ) then
           call stspud( nspin, rhcr, gdcrsv, nspnod, lcgjsv )
       end if

       call wvdump( nfile, myid, nodes,  &
& dname_eig, nstep, ibstt1, ibstt2,  &
& rhcr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, iods,  &
& mfd2ft, ntotfd, nd1vks, kfft0d, rvol,  &
& fft3x, fft3y, fftwork, ijkgd, nplw5, apk, eigr, nspin, lspin,  &
& kfft1b, kfft2b, kfft3b, kfft0b, ijkgb, nplw3, ngboxa, ngboxb, ngboxc,  &
& ntype, natom, iatoit, ioa, lvand, lvandi,  &
& wv_xmin, wv_xmax, wv_ymin, wv_ymax, wv_zmin, wv_zmax,  &
& lcompress_dpwav, ndigit_dpwav )

    end do
    ct = timecnt()
    if(loutfile(1)) write(nfile(1),*) ' dump wavefunctions  : cpu-time :', ct- ct0
    if(loutfile(2)) write(nfile(2),*) ' dump wavefunctions  : cpu-time :', ct- ct0
    ct0 = ct
end if
end if
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!--- dump local potential
!if( ldppot ) then
!if( mod(nstep,nskip_dppot).eq.0 .or. ifmd.eq.0 ) then
!    do nspin = 1, nspnmx
!       if( lspin ) then
!           call ldvloc( nspin, vlocud, vexc, mshnod(1) )
!       end if
!
!       call potdump( nfile, myid, nodes,  &
!& dname_pot, dname_potav, nstep, vexc, vhar, vhshdp, vext, nspin, lspin,  &
!& hdiag, mshnod(1), glocal, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d, fft3x,  &
!& pt_xmin, pt_xmax, pt_ymin, pt_ymax, pt_zmin, pt_zmax, nav_dppot,  &
!& lcompress_dppot, ndigit_dppot )
!
!    end do
!    ct = timecnt()
!    if(loutfile(1)) write(nfile(1),*) ' dump local potential  : cpu-time :', ct- ct0
!    if(loutfile(2)) write(nfile(2),*) ' dump local potential  : cpu-time :', ct- ct0
!    ct0 = ct
!end if
!end if
!-----------------------------------------------------------------------



if( ifmd.ge.1 ) then
    ct = timecnt()
    do i = 1, 2
    if( loutfile(i) ) then
        write(nfile(i),*) ' '
        write(nfile(i),'(a20,f10.4)')  &
&                         ' cpu-time ( sec )  :', ct - ctmd0
    end if
    end do
    ctmd0 = ct
    if( nstep < nstop ) then
!    if( lgamma ) then

        if( lspin ) then
            nspin = 1
            call stspud( nspin, rhcr, gdcrsv, nspnod, lcgjsv )
            do nspin = 1, 2
               if( ltimecnt ) csh0 = timecnt()
!    --- convert band decomposition to G decomposition ---
               call bdtogd( nfile, myid, nodes, csh0,  &
& rhcr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iodg, idstnd, .false., ltimecnt )
!--- save wavefunctions decomposed by components
               call stpsud( nspin, gdcr, gdcrsv, nspnod, lcgjsv )
               if( nspin.eq.1 ) then
!     --- copy gdcr to rhcr ---
                   call cprhcg( nfile, myid, nodes,  &
&                   rhcr, gdcr, nplwex, nbnod )
               end if
            end do
        end if

!    else
!
!        do nspin = 1, nspnmx
!        kvecdoa: do kvec = 1, nknod
!           call set_kvec( kvec )
!           call get_wfk( rhcr, nspin )
!           call get_plwk( nplw, nplwex, nspnod )
!           call get_indexk(  &
!& npnod1, npnod2, npnod, nplcnt, npldsp, nodes )
!
!           leig_end   = kvec == nknod .and. nspin == nspnmx
!!    --- convert band decomposition to G decomposition ---
!           call set_all2all_out( leig_end )
!           if( ltimecnt ) csh0 = timecnt()
!           call bdtogd( nfile, myid, nodes, csh0,  &
!& rhcr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& gdcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
!& iodg, idstnd, .false., ltimecnt )
!
!    !-----store w.f. decomposed by components
!           call set_wfk( gdcr, nspin )
!        end do kvecdoa
!        end do
!
!    end if
    end if
end if


if( lfssh_gsscf .and. .not.lfssh_gsfrc ) then
    !----SCF with the ground state in TDDFT-FSSH
    if( lspin ) then
        call tddft_fssh_restore_occ( nfile, myid, nodes,  &
& nband, nbnod, nspnmx, egvlud, nsordr, wegud )
      else
        call tddft_fssh_restore_occ( nfile, myid, nodes,  &
& nband, nbnod, nspnmx, eig, norder, occ )
    end if

    call tddft_fssh_restore_rhgr( nfile, myid, nodes,  &
& rhgr, nplw5ex, nspnmx )

!    call tddft_fssh_restore_rhoij( nfile, myid, nodes )
end if


!      !==== stop here =========================
!      call fstop( nfile, myid, nodes, 'pwlda' )
!      !========================================


return
end




subroutine pwlda_save( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!     save data for next MS step
!     LDA or GGA calculation with plane-wave method
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
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes
integer :: ierror


!      !==== stop here =========================
!      call fstop( nfile, myid, nodes, 'pwlda_save' )
!      !========================================

ct0 = timecnt()

if( lsave ) then
!=== save data =========================================================
    if(loutfile(1)) write(nfile(1),*) ' '
    if(loutfile(2)) write(nfile(2),*) ' '

!--- store eigenvalues and eigenvectors
!    if( lgamma ) then
        call svwfns( nfile, iogpsz,  &
& fname_eig, cgjr, nplwex, nplw,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, eig, egvlud, occ, wegud,  &
& rhcr, lsreal8, lspin, nspnmx, gdcrsv, nspnod, lcgjsv, ncscale, wvratio )
!    else
!        call sv_k_wfns( nfile,  &
!& fname_eigk, bzk,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, eig, egvlud,  &
!& rhcr, cgjr, lsreal8, lspin, nspnmx, ncscale, wvratio )
!    end if


    !--- store electron density
    call svcdty( nfile,  &
& fname_cds, rhgr, nplw5ex, nplw5, nspnmx, lnoncollinear, mshnod(1),  &
& glocal, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d, fft3x )

    !--- store hartree potential
    call svvhar( nfile,  &
& fname_hrt, vhar, mshnod(1), ldouble_grid,  &
& glocal, mftnod, mftdsp, mfd2ft, ntotfd, kfft0d, fft3x )


!--- store atomic configuration
    if( ifmd.eq.1 ) then
        nstepCG = nstep
    else if( ifmd.ge.2 ) then
        nstepMD = nstep
    end if
    call svacon( nfile, iogpsz,  &
& fname_ion, prevr, ifmd, nstepCG, nstepMD,  &
& ratm, frc, ntype, nhk1, nhk2, lclust, natom, hcell )

!--- store hcell, if lvshape = .true.
    call svhcell( nfile, hcell, lvshape )

!--- store data for tddft
    call save_tddft_fssh( nfile, iogpsz,  &
& cgjr, nplwex, nplw, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, nband,  &
& rhcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iods, iodsg, ioag, idstnd, bufcr, nspnmx )

!--- store data for uniform electric field calculation
!    call save_efield( nfile )

!-----------------------------------------------------------------------
!    store data at previous steps
!-----------------------------------------------------------------------
!--- save charge densities at previous time steps
    call savpcd( nfile,  &
& fname_pcds, ifmd, keypcd, prvrho, nplw5ex, nplw5, ncprvmx, nspnmx2,  &
& lnoncollinear, keypmg,  &
& glocal, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d, fft3x )
!--- save wavefunctions at previous time steps
!    if( lgamma ) then
        call savpwv( nfile, iogpsz,  &
& fname_peig, ifmd, lsreal8, keypwv,  &
& cgjr, nplwex, nplw, nbnod1, nbnod2,  &
& nbnod, nbncnt, nbndsp, nband, prveig,  &
& rhcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iods, iodsg, ioag, idstnd, bufcr, nprvmx, nspnmx, ncscale )
!    else
!        call savpwv_k( nfile,  &
!& fname_peigk, ifmd, lsreal8, keypwv, bzk,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
!& rhcr, cgjr, iods, iodsg, ioag, idstnd, bufcr, nprvmx, nspnmx, ncscale )
!    end if

!=== end of save data ==================================================
    ct = timecnt()
    if(loutfile(1)) write(nfile(1),*) ' dump data ... done  : cpu-time :', ct- ct0
    if(loutfile(2)) write(nfile(2),*) ' dump data ... done  : cpu-time :', ct- ct0
    ct0 = ct
end if


return
end




subroutine xc_nscforce( nfile, myid, nodes )
!-----------------------------------------------------------------------
!   xc contribution to NSC force
!-----------------------------------------------------------------------
use constants
use param
use param_atom
use pwlda_pp
use pwlda_atom
use pwlda_variables
use pwlda_proc
use pwlda_pw
use pwlda_grid
implicit real*8 ( a-h, o-z )
integer :: nfile(*), myid, nodes

!-----declare local variables


if( jgga.eq.1 ) then

    !--- LDA ( perdew-zunger + ceperley-alder )
!    if( lspin ) then
!
!        call xclsda( rho, rhoud, rhocore, mshnod(1), rdelv,  &
!& vlocud, eeexc, noddatx, eexc, evexc, ecor, evcor )
!
!      else
!
!        call xclda( rho, rhocore, mshnod(1), rdelv,  &
!&               vexc, eeexc, eexc, evexc, ecor, evcor )
!
!    end if

 else

    !--- GGA
!    if( lspin ) then
!
!        if( jgga.eq.2 ) then
!            !--- GGA ( PBE )
!            call xcpbe_spin( nfile, myid, nodes,  &
!& vlocud, eeexc, eexc, evexc, ecor, evcor,  &
!& rho, rhoud, rhocore, mshnod(1), rdelv,  &
!& x, rk, xx, vk, tmpk, tmpl, tmpm,  &
!& tmpn, hdiag, vexc, thrhgr, nplw5ex, nplw5,  &
!& glocal,apk,eigr, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d,  &
!& fft3x, fft3y, fftwork, ijkgd, nplw5, gx, gy, gz, t_comm,  &
!& lcstress, strgga )
!          else
!            !--- GGA ( RPBE )
!            call xcrpbe_spin( nfile, myid, nodes,  &
!& vlocud, eeexc, eexc, evexc, ecor, evcor,  &
!& rho, rhoud, rhocore, mshnod(1), rdelv,  &
!& x, rk, xx, vk, tmpk, tmpl, tmpm,  &
!& tmpn, hdiag, vexc, thrhgr, nplw5ex, nplw5,  &
!& glocal,apk,eigr, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d,  &
!& fft3x, fft3y, fftwork, ijkgd, nplw5, gx, gy, gz, t_comm,  &
!& lcstress, strgga )
!        end if
!
!      else
!
!        if( jgga.eq.2 ) then
!
!            if( .not.lvdw .or. lvdw_pre ) then
                !--- GGA ( PBE )
                call xcpbe_nscforce( nfile, myid, nodes,  &
& xexc, rho, rhocore, mshnod(1), rdelv,  &
& x, rk, xx, vk, delrho, ddelrho, thrhgr, nplw5ex, nplw5,  &
& glocal,apk,eigr, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d,  &
& fft3x, fft3y, fftwork, ijkgd, nplw5, gx, gy, gz )
!              else
!                !--- van der Waals DFT
!                call xcvdw( nfile, myid, nodes,  &
!& vexc, eeexc, eexc, evexc, ecor, ecnl, evcor, rho, rhocore,  &
!& rdelv, x, rk, xx, vk, thrhgr, nplw5ex, nplw5,  &
!& glocal,apk,eigr, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d,  &
!& fft3x, fft3y, fftwork, ijkgd, nplw5, gx, gy, gz, t_comm,  &
!& lcstress, strgga, lorthrhmbc, &
!& hcell, nd1v, mshglb, nd1vks(1), nd1vks(2), nd1vks(3), lvacuum,  &
!& mshxyz(mulpit(1)),  &
!& mulx1(1), mulx2(1), muly1(1), muly2(1), mulz1(1), mulz2(1),  &
!& mshnx(mulpms(1)), mshny(mulpms(1)), mshnz(mulpms(1)), mshnod(1),  &
!& mshx1(1), mshy1(1), mshz1(1), mshx(1), mshy(1), mshz(1), ltimecnt )
!            end if
!
!          else
!
!            !--- GGA ( RPBE )
!            call xcrpbe( nfile, myid, nodes,  &
!& vexc, eeexc, eexc, evexc, ecor, evcor, rho, rhocore,  &
!& mshnod(1), rdelv, x, rk, xx, vk, thrhgr, nplw5ex, nplw5,  &
!& glocal,apk,eigr, mftnod, mftdsp, mfd2ft, ntotfd, nd1vks, kfft0d,  &
!& fft3x, fft3y, fftwork, ijkgd, nplw5, gx, gy, gz, t_comm,  &
!& lcstress, strgga )
!        end if
!
!    end if

end if


return
end




subroutine copyQMforces( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!     just copy data, if root nodes are the same in both QM & MD regions
!-----------------------------------------------------------------------
use param
use pwlda_atom
use pwlda_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: ierror


call copyQM2MD( nfile, myid, nodes,  &
& frc, natom, sume, strtot, ierror )


return
end subroutine




subroutine copyMD2QM( nfile, myid, nodes,  &
& realatm_, ndata, hcell_, ierror )
!-----------------------------------------------------------------------
!     just copy data
!-----------------------------------------------------------------------
use outfile
use param
use param_atom
implicit none
integer :: nfile(*), myid, nodes
integer :: ndata
real*8, dimension(ndata) :: realatm_
real*8, dimension(3,3) :: hcell_
integer :: ierror

!-----declare local variables
integer :: i, j, ix


!-----error trap
if( ndata /= 3*natom ) then
    if(loutfile(1)) write(nfile(1),*) 'error in copyMD2QM : ndata /= 3*natom'
    if(loutfile(2)) write(nfile(2),*) 'error in copyMD2QM : ndata /= 3*natom'
    ierror = 1
    return
end if

!-----store atomic positions
do i = 1, natom
   realatm(1,i) = realatm_(3*i-2)
   realatm(2,i) = realatm_(3*i-1)
   realatm(3,i) = realatm_(3*i-0)
end do


!-----store supercell vectors
do i = 1, 3
do ix = 1, 3
   h_bkup(ix,i) = hcell_(ix,i)
end do
end do


return
end subroutine




subroutine sendQMforces( nfile, myid, nodes, id_md1, ierror )
!-----------------------------------------------------------------------
!     send QM forces to MD-root node
!-----------------------------------------------------------------------
use param
use pwlda_atom
use pwlda_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: id_md1
integer :: ierror

!-----declare local variables
integer :: nsd


nsd = 3*natom
call cisend(100,nsd,1,id_md1,0)
call cdsend(110,frc,nsd,id_md1,0)
call cdsend(120,sume,1,id_md1,0)
call cdsend(130,strtot,9,id_md1,0)


return
end subroutine




subroutine recvQMpositions( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!     recieve data from MD nodes
!-----------------------------------------------------------------------
use outfile
use param
use param_atom
implicit none
integer :: nfile(*), myid, nodes
integer :: ierror

!-----declare local variables
integer :: nrc, i, ix


call cirecv(100,nrc,1,0)

!-----error trap
if( nrc /= 3*natom ) then
    if(loutfile(1)) write(nfile(1),*) 'error in recvQMpositions :',  &
&                         ' nrc /= 3*natom'
    if(loutfile(2)) write(nfile(2),*) 'error in recvQMpositions :',  &
&                         ' nrc /= 3*natom'
    ierror = 2
    return
end if

call cdrecv(110,realatm,nrc,0)


!-----store supercell vectors
call cdrecv(120,h_bkup,9,0)


return
end subroutine




subroutine unifQMpositions( nfile, myid, nodes, lupdate, ierror )
!-----------------------------------------------------------------------
!     unify transferred coordinates in QM nodes
!-----------------------------------------------------------------------
use outfile
use param
use param_atom
use pwlda_variables
implicit none
integer :: nfile(*), myid, nodes
logical :: lupdate               ! lupdate = .true. : update coordinates
integer :: ierror

!-----declare local variables
integer :: i, ix, ndq, it, j
real*8  :: rvol_
real*8, dimension(3,3) :: b
real*8  :: q1, dq
real*8  :: small


call dbcast(realatm,3*natom,0)
call dbcast(h_bkup,9,0)

!-----set new supercell vectors hcell, while keep old ones in h_bkup
do i = 1, 3
do ix = 1, 3
   dq           =  hcell(ix,i)
    hcell(ix,i) = h_bkup(ix,i)
   h_bkup(ix,i) = dq
end do
end do


updateif: if( .not.lupdate ) then

!-----set displacement vector
   if( lclust ) then

       i = 1
       do ix = 1, 3
          disp_md(ix) = ratm(ix,i) - realatm(ix,i)
       end do

       !-----error trap
       do ix = 1, 3
       do i = 1, natom
          if( abs(disp_md(ix)-(ratm(ix,i)-realatm(ix,i)))>1.d-01 )  &
&                                                             then
              ierror = 3
              if(loutfile(1)) write(nfile(1),*) 'error in unifQMpositions(1) :',  &
&                          'inconsistent disp_md', i, ix,  &
&                          disp_md(ix)-(ratm(ix,i)-realatm(ix,i))
              if(loutfile(2)) write(nfile(2),*) 'error in unifQMpositions(1) :',  &
&                          'inconsistent disp_md', i, ix,  &
&                          disp_md(ix)-(ratm(ix,i)-realatm(ix,i))
              return
          end if
          if( abs(disp_md(ix)-(ratm(ix,i)-realatm(ix,i)))>1.d-05 )  &
&                                                             then
              do j = 1, 2
              if( loutfile(j) ) then
              write(nfile(j),*) '--------------------------------'
              write(nfile(j),*) 'warning in unifQMpositions(1) :',  &
&                          'inconsistent disp_md', i, ix,  &
&                          disp_md(ix)-(ratm(ix,i)-realatm(ix,i))
              write(nfile(j),*) ' may related to HSQM atoms',  &
&                               ' or symmetry operations'
             ! write(nfile(j),*) '--------------------------------'
              end if
              end do
          end if
       end do
       end do

     else

       !----- transpose matrix of b = inverse of hcell
       CALL RCIPRL( hcell, b, rvol_ )
       i = 1
       do ix = 1, 3
          q1 = b(1,ix)*realatm(1,i) + b(2,ix)*realatm(2,i)  &
&            + b(3,ix)*realatm(3,i)
          disp_md(ix) = ratm(ix,i) - q1
       end do

       !-----error trap
       do ix = 1, 3
       do i = 1, natom
          q1 = b(1,ix)*realatm(1,i) + b(2,ix)*realatm(2,i)  &
&            + b(3,ix)*realatm(3,i)
          dq = disp_md(ix) - ( ratm(ix,i) - q1 )
          ndq = nint(dq)
          if( abs(dq-dble(ndq)) > 1.d-02 ) then
              ierror = 3
              if(loutfile(1)) write(nfile(1),*) 'error in unifQMpositions(2) :',  &
&                 'inconsistent disp_md', i, ix, dq-dble(ndq)
              if(loutfile(2)) write(nfile(2),*) 'error in unifQMpositions(2) :',  &
&                 'inconsistent disp_md', i, ix, dq-dble(ndq)
              return
          end if
          if( abs(dq-dble(ndq)) > 1.d-06 ) then
              do j = 1, 2
              if( loutfile(j) ) then
              write(nfile(j),*) '--------------------------------'
              write(nfile(j),*) 'warning in unifQMpositions(2) :',  &
&                 'inconsistent disp_md', i, ix, dq-dble(ndq)
              write(nfile(j),*) ' may related to HSQM atoms',  &
&                               ' or symmetry operations'
             ! write(nfile(j),*) '--------------------------------'
              end if
              end do
          end if
       end do
       end do

   end if

end if updateif


!-----set coordinates by adding displacement vector
if( lclust ) then
    do ix = 1, 3
       do i = 1, natom
          ratm(ix,i) = realatm(ix,i) + disp_md(ix)
       end do
    end do
  else
    !----- transpose matrix of b = inverse of hcell
    CALL RCIPRL( hcell, b, rvol_ )
    small = 1.d-14
    do ix = 1, 3
       do i = 1, natom
          q1 = b(1,ix)*realatm(1,i) + b(2,ix)*realatm(2,i)  &
&            + b(3,ix)*realatm(3,i)
          q1 = q1 + disp_md(ix)
          if( abs(q1) < small ) q1 = 0.d0
          if( q1  < 0.d0 ) q1 = q1 - int(q1) + 1.d0
          if( q1 >= 1.d0 ) q1 = q1 - int(q1)
          if( abs(q1-1.d0) < small ) q1 = 0.d0
          ratm(ix,i) = q1
       end do
    end do
end if

!-----check
!      if( myid == 0 ) then
!         do i = 1, 2
!            write(nfile(i),*) 'check in unifQMpositions '
!            do it = 1, ntype
!               write(nfile(i),*) ' it=', it
!               write(nfile(i),*) ' nhk1, nhk2=', nhk1(it), nhk2(it)
!               write(nfile(i),*) ' displacement'
!               write(nfile(i),*) ( disp_md(ix), ix = 1,3 )
!               write(nfile(i),*) ' coordinates'
!               do j = nhk1(it), nhk2(it)
!                  write(nfile(i),'(i5,3f12.8)') j, (ratm(ix,j), ix =1,3)
!               end do
!            end do
!         end do
!      end if


return
end subroutine




subroutine copy_velocities_MD2QM( nfile, myid, nodes,  &
& vatm_, ndata, ierror )
!-----------------------------------------------------------------------
!     just copy data
!-----------------------------------------------------------------------
use outfile
use param
use param_atom
implicit none
integer :: nfile(*), myid, nodes
integer :: ndata
real*8, dimension(ndata) :: vatm_
integer :: ierror

!-----declare local variables
integer :: i, j, ix


!-----error trap
if( ndata /= 3*natom ) then
    if(loutfile(1)) write(nfile(1),*) 'error in copy_velocities_MD2QM : ndata /= 3*natom'
    if(loutfile(2)) write(nfile(2),*) 'error in copy_velocities_MD2QM : ndata /= 3*natom'
    ierror = 1
    return
end if

!-----store atomic positions
do i = 1, natom
   vatm(1,i) = vatm_(3*i-2)
   vatm(2,i) = vatm_(3*i-1)
   vatm(3,i) = vatm_(3*i-0)
end do


return
end subroutine




subroutine recvQMvelocities( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!     recieve data from MD nodes
!-----------------------------------------------------------------------
use outfile
use param
use param_atom
implicit none
integer :: nfile(*), myid, nodes
integer :: ierror

!-----declare local variables
integer :: nrc, i, ix


call cirecv(100,nrc,1,0)

!-----error trap
if( nrc /= 3*natom ) then
    if(loutfile(1)) write(nfile(1),*) 'error in recvQMvelocities :',  &
&                         ' nrc /= 3*natom'
    if(loutfile(2)) write(nfile(2),*) 'error in recvQMvelocities :',  &
&                         ' nrc /= 3*natom'
    ierror = 2
    return
end if

call cdrecv(110,vatm,nrc,0)


return
end subroutine




subroutine unifQMvelocities( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!     unify transferred velocities in QM nodes
!-----------------------------------------------------------------------
use outfile
use constants
use param
use param_atom
use pwlda_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: ierror

!-----declare local variables
integer :: i, it
real*8  :: tempchk


call dbcast(vatm,3*natom,0)

do i = 1, natom
   vatm(1:3,i) = vatm(1:3,i) * 2.d0    ! in [Rydberg units]
end do

!check---calculate temperature
tempchk = 0.d0
do it = 1, ntype
   do i = nhk1(it), nhk2(it)
      tempchk = tempchk + watom(it)*( &
& vatm(1,i)*vatm(1,i) + vatm(2,i)*vatm(2,i) + vatm(3,i)*vatm(3,i) )
   end do
end do
tempchk = 0.5d0 * tempchk
if(loutfile(1)) write(nfile(1),*) '*** check temperature in unifQMvelocities:', &
& tempchk / natom * (tempau/2.d0) * 2.d0/3.d0
if(loutfile(2)) write(nfile(2),*) '*** check temperature in unifQMvelocities:', &
& tempchk / natom * (tempau/2.d0) * 2.d0/3.d0

!-----set kinetic energy in TDDFT-FSSH
call set_ekin_in_tddft_fssh( tempchk )


return
end subroutine




module vsendback
!-----------------------------------------------------------------------
! type declaration and initialization for copying back to MD nodes
!-----------------------------------------------------------------------
implicit none

logical :: lvsendback

save
end module




subroutine set_lvsendback( nfile, myid, nodes, ltransition )
!-----------------------------------------------------------------------
!    set lvsendback
!-----------------------------------------------------------------------
use param
use vsendback
implicit none
integer :: nfile(*), myid, nodes
logical :: ltransition

lvsendback = ltddft_fssh .and. lfssh_vscale .and. ltransition

return
end subroutine




subroutine copy_velocities_toMD( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!     just copy data, if root nodes are the same in both QM & MD regions
!-----------------------------------------------------------------------
use param
use param_atom
use pwlda_atom
use pwlda_variables
use vsendback
implicit none
integer :: nfile(*), myid, nodes
integer :: ierror


call copylvsendback( nfile, myid, nodes, lvsendback )

if( .not.lvsendback ) return

call copyQM2MDv( nfile, myid, nodes,  &
& vatm, natom, ierror )


return
end subroutine




subroutine send_velocities_toMD( nfile, myid, nodes, id_md1, ierror )
!-----------------------------------------------------------------------
!     send QM velocities to MD-root node
!-----------------------------------------------------------------------
use param
use param_atom
use pwlda_atom
use pwlda_variables
use vsendback
implicit none
integer :: nfile(*), myid, nodes
integer :: id_md1
integer :: ierror

!-----declare local variables
integer :: nsd


call clsend(90,lvsendback,1,id_md1,0)

if( .not.lvsendback ) return

nsd = 3*natom
call cisend(100,nsd,1,id_md1,0)
call cdsend(110,vatm,nsd,id_md1,0)


return
end subroutine




subroutine statom( nfile, myid, nodes,  &
& ntype, nhk1, nhk2, nhk1_nod, nhk2_nod, iatmpt, nion, natom )
!-----------------------------------------------------------------------
!     set No. of atoms in each node
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: ntype
integer :: natom
integer, dimension(ntype) :: nhk1, nhk2
integer, dimension(ntype) :: nhk1_nod, nhk2_nod
integer, dimension(natom) :: iatmpt
integer :: nion

!-----declare local variables
integer :: ntot, i, ncur, it, ia


ntot = nhk2(ntype)

nion = ntot/nodes
do i = 1, nion
   if( mod(i,2).eq.1 ) then
       iatmpt(i) = nodes*( i - 1 ) + myid + 1
     else
       iatmpt(i) = nodes*i - myid
   end if
end do

ncur = nion*nodes
do i = 1, ntot - ncur
   if( myid.eq.nodes-i ) then
       nion = nion + 1
       iatmpt(nion) = ncur + i
   end if
end do


!-----set nhk1_nod, nhk2_nod
do it = 1, ntype
   nhk2_nod(it) = 0
   do i = 1, nion
      ia = iatmpt(i)
      if( ia >= nhk1(it) .and. ia <= nhk2(it) )  &
&         nhk2_nod(it) = nhk2_nod(it) + 1
   end do
end do
nhk1_nod(1) = 1
nhk2_nod(1) = nhk1_nod(1) + nhk2_nod(1) - 1
do it = 2, ntype
   nhk1_nod(it) = nhk2_nod(it-1) + 1
   nhk2_nod(it) = nhk1_nod(it) + nhk2_nod(it) - 1
end do


return
end subroutine




subroutine struc_atoms_alloc( nfile, myid, nodes,  &
& alloc_mem, natom, nplwcs, kmax1, kmax2, kmax3,  &
& kmax1_, kmax2_, kmax3_, kmax1d, kmax2d, kmax3d,  &
& lnlpp_g, llclpp_g, lpcc_g, lsphexp, pwscale )
!-----------------------------------------------------------------------
!     allocate memory for shared variables to calculate structure factor
!-----------------------------------------------------------------------
use pwlda_pw
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: natom
integer :: nplwcs
integer :: kmax1,  kmax2,  kmax3
integer :: kmax1_, kmax2_, kmax3_
integer :: kmax1d, kmax2d, kmax3d
logical :: lnlpp_g, llclpp_g, lpcc_g, lsphexp
real*8  :: pwscale

!-----declare local variables
integer :: status
real*8  :: the_mem
integer :: mplnplw = 0
integer :: mplkmax1 = 0
integer :: mplkmax2 = 0
integer :: mplkmax3 = 0
logical :: licall = .true.
save mplnplw, mplkmax1, mplkmax2, mplkmax3, licall


if( llclpp_g .or. lpcc_g .or. lsphexp ) then
    kmax1 = kmax1d
    kmax2 = kmax2d
    kmax3 = kmax3d
  else
    kmax1 = kmax1_
    kmax2 = kmax2_
    kmax3 = kmax3_
end if


if( lnlpp_g .or. llclpp_g .or. lpcc_g .or. lsphexp ) then

    if( lnlpp_g .or. lpcc_g .or. lsphexp ) then

        !-----check array sizes
        if( (nplwcs+1)*natom > mplnplw ) then

            !-----if already allocated, deallocate arrays
            if( allocated(ycos) ) then
                the_mem = 8.d0 * ( size(ycos) + size(ysin) )

                deallocate( ycos, ysin, stat=status )

                !------error trap
                call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'struc_atoms_alloc', .true. )
            end if

            !------allocate memory
            mplnplw = (nplwcs * pwscale + 1)*natom
            allocate( ycos(mplnplw), ysin(mplnplw), stat=status )

            the_mem = 8.d0 * ( size(ycos) + size(ysin) )

            !------error trap
            call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'struc_atoms_alloc', .true. )

        end if

      else

        nplwcs = 0
        if( licall ) then
            !------allocate memory, just for safe
            allocate( ycos(1), ysin(1), stat=status )
        end if

    end if


    !-----check array sizes
    if( (2*kmax1+1)*natom > mplkmax1 .or.  &
&       (2*kmax2+1)*natom > mplkmax2 .or.  &
&       (2*kmax3+1)*natom > mplkmax3     ) then

        !-----if already allocated, deallocate arrays
        if( allocated(DCQ1) ) then
            the_mem = 8.d0 * ( size(DCQ1) + size(DCQ2) + size(DCQ3)  &
&                            + size(DSQ1) + size(DSQ2) + size(DSQ3) )

            deallocate( DCQ1, DCQ2, DCQ3, DSQ1, DSQ2, DSQ3,  &
& stat=status )

            !------error trap
            call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'struc_atoms_alloc(2)', .true. )
        end if

        !------allocate memory
        mplkmax1 = (2*kmax1 * pwscale + 1)*natom
        mplkmax2 = (2*kmax2 * pwscale + 1)*natom
        mplkmax3 = (2*kmax3 * pwscale + 1)*natom
        allocate( DCQ1(mplkmax1), DCQ2(mplkmax2), DCQ3(mplkmax3),  &
&                 DSQ1(mplkmax1), DSQ2(mplkmax2), DSQ3(mplkmax3),  &
& stat=status )

        the_mem = 8.d0 * ( size(DCQ1) + size(DCQ2) + size(DCQ3)  &
&                        + size(DSQ1) + size(DSQ2) + size(DSQ3) )

        !------error trap
        call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'struc_atoms_alloc(2)', .true. )

    end if

    !---for k-point sampling
!    call struc_atoms_alloc_k( nfile, myid, nodes,  &
!& alloc_mem, natom, nplwcs, pwscale )

else

    if( licall ) then
        !------allocate memory, just for safe
        allocate( ycos(1), ysin(1),  &
& DCQ1(1), DCQ2(1), DCQ3(1), DSQ1(1), DSQ2(1), DSQ3(1),  &
& stat=status )
    end if

end if


licall = .false.

return
end subroutine




subroutine stscfn_pw( nfile, myid, nodes,  &
& ycos, ysin, DCQ1, DCQ2, DCQ3, DSQ1, DSQ2, DSQ3,  &
& kmax1, kmax2, kmax3, natom, nplwcs, nplw5, nga, ngb, ngc,  &
& lclust, ratm, hcell, nd1v, nd1vks, lnlpp_g, llclpp_g, lpcc_g, lsphexp )
!-----------------------------------------------------------------------
!     sine, cosine functions for atoms
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: natom
integer :: nplwcs
integer :: nplw5
integer :: kmax1, kmax2, kmax3
real*8,  dimension(0:nplwcs,*) :: ycos, ysin
real*8,  dimension(-kmax1:kmax1,natom) :: DCQ1, DSQ1
real*8,  dimension(-kmax2:kmax2,natom) :: DCQ2, DSQ2
real*8,  dimension(-kmax3:kmax3,natom) :: DCQ3, DSQ3
integer, dimension(0:nplw5) :: nga, ngb, ngc
logical :: lclust
real*8,  dimension(3,natom) :: ratm
real*8,  dimension(3,3) :: hcell  ! supercell vectors
integer, dimension(3) :: nd1v     ! the number of FD  meshes
integer, dimension(3) :: nd1vks   ! the number of FFT meshes
logical :: lnlpp_g, llclpp_g, lpcc_g, lsphexp

!-----declare local variables
real*8, parameter :: dpi = 3.141592653589793d0 * 2.d0
real*8  :: fact1, fact2, fact3
real*8  :: D2P1, DCSN1, DSGN1, D2P2, DCSN2, DSGN2,  &
&          D2P3, DCSN3, DSGN3, DCC12, DSC12
integer :: K, KM, KMM, ig, K1M, K2M, K3M


if( lclust ) then
    !-----if cluster, ratm are real coordinates,
    fact1 = 1.d0/hcell(1,1) * dble(nd1v(1))/dble(nd1vks(1))
    fact2 = 1.d0/hcell(2,2) * dble(nd1v(2))/dble(nd1vks(2))
    fact3 = 1.d0/hcell(3,3) * dble(nd1v(3))/dble(nd1vks(3))
  else
    !-----if cluster, ratm are scaled coordinates for nd1v
    fact1 = dble(nd1v(1))/dble(nd1vks(1))
    fact2 = dble(nd1v(2))/dble(nd1vks(2))
    fact3 = dble(nd1v(3))/dble(nd1vks(3))
end if

fact1 = fact1 * dpi
fact2 = fact2 * dpi
fact3 = fact3 * dpi
do K  = 1, natom

   !-----for ratm(1)
   D2P1  = ratm(1,K)*fact1
   DCSN1 = COS(D2P1)
   DSGN1 = SIN(D2P1)

   DCQ1(0,K)  = 1.0D0
   DSQ1(0,K)  = 0.0
   DCQ1(1,K)  = DCSN1
   DSQ1(1,K)  = DSGN1
   DCQ1(-1,K) = DCSN1 
   DSQ1(-1,K) = -DSGN1

   !-----for ratm(1)
   D2P2  = ratm(2,K)*fact2
   DCSN2 = COS(D2P2)
   DSGN2 = SIN(D2P2)

   DCQ2(0,K)  = 1.0D0
   DSQ2(0,K)  = 0.0
   DCQ2(1,K)  = DCSN2
   DSQ2(1,K)  = DSGN2
   DCQ2(-1,K) = DCSN2 
   DSQ2(-1,K) = -DSGN2

   !-----for ratm(1)
   D2P3  = ratm(3,K)*fact3
   DCSN3 = COS(D2P3)
   DSGN3 = SIN(D2P3)

   DCQ3(0,K)  = 1.0D0
   DSQ3(0,K)  = 0.0
   DCQ3(1,K)  = DCSN3
   DSQ3(1,K)  = DSGN3
   DCQ3(-1,K) = DCSN3 
   DSQ3(-1,K) = -DSGN3

end do


do KM = 2, kmax1
   KMM = KM - 1
   do K  = 1, natom
      DCQ1(KM,K)  = DCQ1(1,K)*DCQ1(KMM,K) - DSQ1(1,K)*DSQ1(KMM,K)
      DSQ1(KM,K)  = DSQ1(1,K)*DCQ1(KMM,K) + DCQ1(1,K)*DSQ1(KMM,K)
      DCQ1(-KM,K) = DCQ1(KM,K)
      DSQ1(-KM,K) = -DSQ1(KM,K)
   end do
end do

do KM = 2, kmax2
   KMM = KM - 1
   do K  = 1, natom
      DCQ2(KM,K)  = DCQ2(1,K)*DCQ2(KMM,K) - DSQ2(1,K)*DSQ2(KMM,K)
      DSQ2(KM,K)  = DSQ2(1,K)*DCQ2(KMM,K) + DCQ2(1,K)*DSQ2(KMM,K)
      DCQ2(-KM,K) = DCQ2(KM,K)
      DSQ2(-KM,K) = -DSQ2(KM,K)
   end do
end do

do KM = 2, kmax3
   KMM = KM - 1
   do K  = 1, natom
      DCQ3(KM,K)  = DCQ3(1,K)*DCQ3(KMM,K) - DSQ3(1,K)*DSQ3(KMM,K)
      DSQ3(KM,K)  = DSQ3(1,K)*DCQ3(KMM,K) + DCQ3(1,K)*DSQ3(KMM,K)
      DCQ3(-KM,K) = DCQ3(KM,K)
      DSQ3(-KM,K) = -DSQ3(KM,K)
   end do
end do


if( lnlpp_g .or. lpcc_g .or. lsphexp ) then
    do K = 1, natom
       ycos(0,K) = 1.d0
       ysin(0,K) = 0.d0
       do ig = 1, min( nplwcs, nplw5 )
          K1M = nga(ig)
          K2M = ngb(ig)
          K3M = ngc(ig)
          DCC12 = DCQ1(K1M,K)*DCQ2(K2M,K) -  &
&                 DSQ1(K1M,K)*DSQ2(K2M,K)
          DSC12 = DSQ1(K1M,K)*DCQ2(K2M,K) +  &
&                 DCQ1(K1M,K)*DSQ2(K2M,K)

          ycos(ig,K) = DCC12*DCQ3(K3M,K) - DSC12*DSQ3(K3M,K)
          ysin(ig,K) = DSC12*DCQ3(K3M,K) + DCC12*DSQ3(K3M,K)
!            ycos(-ig,K) =  ycos(ig,K)
!            ysin(-ig,K) = -ysin(ig,K)
       end do
    end do
end if


return
end subroutine




subroutine stpsud( nspin, gdcr, gdcrsv, nspnod, lcgjsv )
!-----------------------------------------------------------------------
!     store wavefunctions for nspin
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension gdcr(*)
dimension gdcrsv(*)
logical   lcgjsv

if( nspin.eq.1 .and. lcgjsv .or.  &
&   nspin.eq.2 .and. .not.lcgjsv ) then
    lcgjsv = .not.lcgjsv
    do ig = 1, nspnod
       buffr      = gdcrsv(ig)
       gdcrsv(ig) = gdcr(ig)
       gdcr(ig)   = buffr
    end do
end if

return
end


subroutine stspud( nspin, gdcr, gdcrsv, nspnod, lcgjsv )
!-----------------------------------------------------------------------
!     restore wavefunctions for nspin
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension gdcr(*)
dimension gdcrsv(*)
logical   lcgjsv

if( nspin.eq.1 .and. .not.lcgjsv .or.  &
&   nspin.eq.2 .and. lcgjsv ) then
    lcgjsv = .not.lcgjsv
    do ig = 1, nspnod
       buffr      = gdcrsv(ig)
       gdcrsv(ig) = gdcr(ig)
       gdcr(ig)   = buffr
    end do
end if

return
end


subroutine exrhcr( rhcr, rhcrsv, nspnod )
!-----------------------------------------------------------------------
!     exchange HC products for spin
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension rhcr(*)
dimension rhcrsv(*)

do ig = 1, nspnod
   buffr      = rhcrsv(ig)
   rhcrsv(ig) = rhcr(ig)
   rhcr(ig)   = buffr
end do

return
end


subroutine ldocc( nspin, occ, wegud, nband )
!-----------------------------------------------------------------------
!     restore occupations for nspin
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension occ(*)
dimension wegud(nband,2)

do ib = 1, nband
   occ(ib) = wegud(ib,nspin)
end do

return
end


subroutine ldocck( nspin, occ, wegud, nband, nkpnt )
!-----------------------------------------------------------------------
!     restore occupations for nspin
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension occ(nband,nkpnt)
dimension wegud(nband,2,nkpnt)

do k = 1, nkpnt
do ib = 1, nband
   occ(ib,k) = wegud(ib,nspin,k)
end do
end do

return
end


subroutine sveig( nspin, eig, egvlud, nband )
!-----------------------------------------------------------------------
!     store eigenvalues for nspin
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension eig(*)
dimension egvlud(nband,2)

do ib = 1, nband
   egvlud(ib,nspin) = eig(ib)
end do

return
end


subroutine sveig_k( nspin, eig, egvlud, nband, nkpnt )
!-----------------------------------------------------------------------
!     store eigenvalues for nspin
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension eig(nband,nkpnt)
dimension egvlud(nband,2,nkpnt)

do k = 1, nkpnt
do ib = 1, nband
   egvlud(ib,nspin,k) = eig(ib,k)
end do
end do

return
end


subroutine ldeig( nspin, eig, egvlud, nband )
!-----------------------------------------------------------------------
!     restore eigenvalues for nspin
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension eig(*)
dimension egvlud(nband,2)

do ib = 1, nband
   eig(ib) = egvlud(ib,nspin)
end do

return
end



subroutine ldvloc( nspin, vlocud, vexc, mshnod )
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension vlocud(mshnod,*)
dimension vexc(mshnod)

do i = 1, mshnod
   vexc(i) = vlocud(i,nspin)
end do

return
end


subroutine sum_eig_spin( sume, wegud, eigtmp, nband, nkpnt )
!-----------------------------------------------------------------------
implicit none
integer :: nband, nkpnt
real*8  :: sume
real*8,  dimension(nband,2,nkpnt) :: wegud
real*8,  dimension(nband,nkpnt,2) :: eigtmp

!-----declare local variables
integer :: nspin, ik, ib

sume = 0.d0
do nspin = 1, 2
do ik = 1, nkpnt
do ib = 1, nband
   sume = sume + wegud(ib,nspin,ik)*eigtmp(ib,ik,nspin)
end do
end do
end do

return
end




subroutine cal_bfzansa( nfile, myid, nodes,  &
& bfzansa1, bfzansa2, bzan, occ, nband, nbnod1, nbnod2, nbnod, iod,  &
& nkpnt, nspin, nnspin, k )
!-----------------------------------------------------------------------
!     restore occupations for nspin
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
real*8  :: bfzansa1, bfzansa2
integer :: nband, nbnod1, nbnod2, nbnod, iod(*), nkpnt, nspin, nnspin, k
real*8  :: bzan(nband), occ(nband,nnspin,nkpnt)

!------declare local variables
integer :: ib, ibb
real*8  :: thebfzansa1, thebfzansa2, buf1d


thebfzansa1 = 0.d0
thebfzansa2 = 0.d0
do ibb = 1, nbnod
   ib = iod(ibb)
   thebfzansa1 = max( thebfzansa1, bzan(ibb) )
   thebfzansa2 = thebfzansa2 + bzan(ibb)*occ(ib,nspin,k)
end do

call gdmax(thebfzansa1)
call gdsum(thebfzansa2,1,buf1d)

bfzansa1 = max( bfzansa1, thebfzansa1 )
bfzansa2 = bfzansa2 + thebfzansa2


return
end
