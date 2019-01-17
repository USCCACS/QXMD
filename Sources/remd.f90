



subroutine remd_set( nfile, myid, nodes, npx_, npy_, npz_, iogpsz_,  &
& lpureMD_, ierror )
!-----------------------------------------------------------------------
!     initial setting for MD calculation
!-----------------------------------------------------------------------
use constants
use remd_param
use remd_param_atom
use remd_variables
use remd_atom
use remd_qmmd
implicit none
integer :: nfile(*), myid, nodes
integer :: npx_, npy_, npz_, iogpsz_
logical :: lpureMD_
integer :: ierror

!-----declare local variables
real*8  :: ct, ct0, timecnt
logical :: lbath_clear = .true.


ct0 = timecnt()

lpureMD = lpureMD_
iogpsz  = iogpsz_

lstat   = lpureMD     ! Default set for statistical calculation
!-----Start parallel environment and keep my node ID --
npx = npx_
npy = npy_
npz = npz_
myx=myid/(npy*npz)
myy=mod(myid/npz,npy)
myz=mod(myid,npz)
!-----Reduced node origin
anxi=1d0/npx
anyi=1d0/npy
anzi=1d0/npz
sxog=anxi*myx
syog=anyi*myy
szog=anzi*myz
!-----Prepares a neighbor-node-ID table
!call ntset_md( nfile, myid, nodes,  &
!& myx, myy, myz, npx, npy, npz, anxi, anyi, anzi )


!----- read and check input data
call remd_input( nfile, myid, nodes, ierror )
if( ierror /= 0 ) return

ct = timecnt()
if( myid.eq.0 ) then
    write(nfile(1),*) ' read data (remd)   : cpu-time :', ct - ct0
    write(nfile(2),*) ' read data (remd)   : cpu-time :', ct - ct0
end if
ct0 = ct
!-----------------------------------------------------------------------


!-----allocate memory
call remd_atom_alloc( nfile, myid, nodes,  &
& alloc_mem, nmd, ntmd, nmdmax, nmd_alloc, natom,  &
& ifmd, lQMMDrun, lpureMD, volume, rc_buffer, ioptmze, lstress, latomic,  &
& lheat, lvmd, lidealref )

call remd_qmmd_alloc( nfile, myid, nodes,  &
& alloc_mem, natom, nHSQM, nmd_cluster, nHSMD, nprmax )


!if( lQMMDrun .or. lpureMD .or. .not.lidealref ) then
!    call remd_ba_var_alloc( nfile, myid, nodes,  &
!& alloc_mem, nmd_alloc, nmdmax, nmd_buffer, nmd_buffer_1d, rc_buffer, h )
!end if


!if( lQMMDrun .or. lpureMD ) then
!    call md_potential_alloc( nfile, myid, nodes,  &
!& alloc_mem, ntmd_ist, nprmax, nmd_buffer, nmdmax, nmd_alloc )
!
!    !-----Prepare potential tables, V 
!    call potset( nfile, myid, nodes )
!end if


!if( .not.lidealref ) then
!    call vmd_potential_alloc( nfile, myid, nodes,  &
!& alloc_mem, ntmd_ist, nprmax, nmd_buffer, nmdmax, nmd_alloc )
!
!    !-----Prepare potential tables, V 
!    call vmd_potset( nfile, myid, nodes )
!end if


!call remd_symm_alloc( nfile, myid, nodes,  &
!& alloc_mem, nmdmax )


!----- Allocate I/0 files
call remd_open_files( nfile, myid, nodes, ierror )

call remd_title_in_files( nfile, myid, nodes, ierror )


!-----read atomic configuration from files, if lstart == .true.
call remd_rdatoms( nfile, myid, nodes, npx, npy, npz, iogpsz,  &
& fname_toc, fname_mts, lstart, ifmd, nstepCG, nstepMD, dtmd,  &
& lQMMDrun, lpureMD,  &
& natom, ntype, nhk1, nhk2, ratm, lterminator, natomx, ntypex,  &
& nmd_cluster, ntmd_cluster, nhk1_cluster, nhk2_cluster, x_cluster,  &
& lMDterminator, nmdx_cluster, ntmdx_cluster,  &
& nmd, is, x, h, frc, pixcg, fnorm, nmd_alloc, nmdmax, nmdmax_,  &
& nnos, xlogs, vlogs, lbath_clear, pintlr, pintsr, pint, ierror,  &
& ioptmze, ioptmze_cell, lqninitial, xmsst, vxmsst, h0,  &
& nnos_atom, xlogs_atom, vlogs_atom )

!-----for virtual MD for thermodynamic integration
!if( lvmd .and. .not.lstart ) then
!    call select_lattice_vmd( nfile, myid, nodes )
!end if


!-----prepare basic parameters
call setup( nfile, myid, nodes,  &
& is, nmd, ntmd, ntot, fack, acon, watom_md, dtmd )


!-----check symmetry of atomic configuration
!call symcod( nfile, myid, nodes,  &
!& x(1,0), nmd, ntmd, ntot, h )


!----- Statistical calculation
!if( lstat ) call remd_init_stat( nfile, myid, nodes,  &
!& alloc_mem, ntmd, ntot, nmdmax, nmd_buffer, h, dtmd )


!-----initial velocities
call ivelct( nfile, myid, nodes,  &
& x(1,0), x(1,1), is, nmd, ntmd, ntot, watom_md, fack, h,  &
& sxog, syog, szog, treq, dtmd, ifmd, nstepMD, tempau, vrandom )


!-----Initialize the extended dynamics
if( ifmd == 3 ) then
    call init_nhc( nfile, myid, nodes,  &
& nnos, nresn, nyosh, tomega, gnkt, gkt,  &
& xlogs, vlogs, glogs, qmass, wdti, wdti2, wdti4, wdti8,  &
& ntot(0), treq, dtmd, lbath_clear, ierror )

else if( ifmd == 4 ) then
    call init_nhcp( nfile, myid, nodes,  &
& nnos, nresn, nyosh, tomega, gnkt, gkt,  &
& xlogs, vlogs, glogs, qmass, wdti, wdti2, wdti4, wdti8,  &
& tbomega, blkmod, volum0, gnd2kt, onf, bmass, fack2, fack, ntmd,  &
& h, ntot(0), treq, dtmd, lbath_clear, ierror )

!else if( ifmd == 5 ) then
!    call init_nhc_atom( nfile, myid, nodes,  &
!& lathermo, nnos, nresn, nyosh, tomega, gnkt, gkt,  &
!& xlogs, vlogs, glogs, qmass, wdti, wdti2, wdti4, wdti8,  &
!& ntmd, ntot, treq, dtmd, lbath_clear, ierror )
!
!    call init_nhc_atom2( nfile, myid, nodes,  &
!& lathermo, nnos_atom, nresn_atom, nyosh_atom, tomega_atom, g3kt, gkt,  &
!& xlogs_atom, vlogs_atom, qmass_atom,  &
!& wdti_atom, wdti2_atom, wdti4_atom, wdti8_atom,  &
!& is, nmd, ntmd, ntot, treq, dtmd, lbath_clear, ierror )

!else if( ifmd == 6 ) then
!    call init_nhc_atom( nfile, myid, nodes,  &
!& lathermo, nnos, nresn, nyosh, tomega, gnkt, gkt,  &
!& xlogs, vlogs, glogs, qmass, wdti, wdti2, wdti4, wdti8,  &
!& ntmd, ntot, treq, dtmd, lbath_clear, ierror )

!    call set_glethemo( nfile, myid, nodes, dtmd, nmd, lbath_clear )


else if( ifmd == 10 ) then
    call init_msst( nfile, myid, nodes,  &
& nnos, nresn, nyosh, tomega, gnkt, gkt,  &
& xlogs, vlogs, glogs, qmass, wdti, wdti2, wdti4, wdti8,  &
& tbomega, blkmod, h0, volum0, alpmatrix, betmatrix,  &
& gnd2kt, onf, bmass, fack2, fack, ntmd,  &
& h, totmass, watom_md, ntot, treq, dtmd, lbath_clear, ierror )
end if

!-----set external stress tensor
if( ifmd == 4 .or. ifmd == 10 .or. ioptmze_cell >= 0 ) then
    call set_pext( nfile, myid, nodes, pext, hpext )
end if


!-----looking for atoms related to QM terminators (HSQM)
!if( .not.lpureMD ) then
!    call init_HSQM( nfile, myid, nodes,  &
!& natom, ntype, nhk1, nhk2, ratm, lterminator,  &
!& x, nmd, h, sxog, syog, szog,  &
!& nHSQM, x_HSQM_qm, x_HSQM_md, alpha_HSQM, nHSQMx, ierror )
!end if

!-----looking for atoms related to MD terminators (HSMD)
!if( lQMMDrun ) then
!    call init_HSQM( nfile, myid, nodes,  &
!& nmd_cluster, ntmd_cluster, nhk1_cluster, nhk2_cluster, x_cluster,  &
!& lMDterminator,  &
!& x, nmd, h, sxog, syog, szog,  &
!& nHSMD, x_HSMD_qm, x_HSMD_md, alpha_HSMD, nHSMDx, ierror )
!end if


!-----set time step for MD
if( ifmd.eq.0 ) then
    nstep = max( nstepCG, nstepMD )
else if( ifmd.eq.1 ) then
    if( .not.lstart .or. lstart .and. nstepCG == 0 ) then
        nstepCG = nstep_ini
      else
        nstep_ini = 0
    end if
    nstep = nstepCG
else
    if( .not.lstart .or. lstart .and. nstepMD == 0 ) then
        nstepMD = nstep_ini
      else
        nstep_ini = 0
    end if
    nstep = nstepMD
end if


return
end subroutine




subroutine set_param_md(  nfile, myid, nodes,  &
& nstep_, nstop_, ifmd_, nstep_ini_, lmdstop_ )
!-----------------------------------------------------------------------
!     set basic parameters for MD calculation
!-----------------------------------------------------------------------
use constants
use remd_param
use remd_param_atom
use remd_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: nstep_, nstop_, ifmd_, nstep_ini_
logical :: lmdstop_

nstep_ = nstep
nstop_ = nstop
ifmd_  = ifmd
nstep_ini_ = nstep_ini
lmdstop_ = lmdstop

return
end subroutine




subroutine check_param_in_md(  nfile, myid, nodes,  &
& nstep_, nstop_, ifmd_, ierror  )
!-----------------------------------------------------------------------
!     set basic parameters for MD calculation
!-----------------------------------------------------------------------
use constants
use remd_param
use remd_param_atom
use remd_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: nstep_, nstop_, ifmd_
integer :: ierror

!-----error trap
if( nstep_ /= nstep ) ierror = 1
if( nstop_ /= nstop ) ierror = ierror + 10
if( ifmd_  /= ifmd  ) ierror = ierror + 100
if( ierror /= 0 ) then
    write(*,*) 'myid=',myid,' error in check_param_in_md',  &
&   nstep_, nstop_, ifmd_, nstep, nstop, ifmd
end if

return
end subroutine




subroutine set_nstep_in_md( nfile, myid, nodes, nstep_ )
!-----------------------------------------------------------------------
!     set time step
!-----------------------------------------------------------------------
use remd_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: nstep_

nstep = nstep_

return
end subroutine




subroutine get_origin( nfile, myid, nodes, sxog_, syog_, szog_ )
!-----------------------------------------------------------------------
!     initial setting for MD calculation
!-----------------------------------------------------------------------
use remd_param
implicit none
integer :: nfile(*), myid, nodes
real*8 :: sxog_, syog_, szog_

sxog_ = sxog
syog_ = syog
szog_ = szog

return
end subroutine




subroutine get_allocmem( nfile, myid, nodes, alloc_mem_ )
!-----------------------------------------------------------------------
!     get alloc_mem
!-----------------------------------------------------------------------
use remd_param
implicit none
integer :: nfile(*), myid, nodes
real*8 :: alloc_mem_

alloc_mem_ = alloc_mem

return
end subroutine




subroutine set_allocmem( nfile, myid, nodes, alloc_mem_ )
!-----------------------------------------------------------------------
!     set alloc_mem
!-----------------------------------------------------------------------
use remd_param
implicit none
integer :: nfile(*), myid, nodes
real*8 :: alloc_mem_

alloc_mem = alloc_mem_

return
end subroutine




subroutine remd_updtconfig( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!     update atomic configuration
!-----------------------------------------------------------------------
use constants
use remd_param
use remd_param_atom
use remd_variables
use remd_atom
implicit none
integer :: nfile(*), myid, nodes
integer :: ierror

!-----declare local variables
logical :: lclear, lclearcellh
logical :: lrunopt, lrunopt_cell
integer :: nstep0, nstepini = -1, nstep1, nstep2
save nstepini


if( nstep >= nstep_ini+1 ) then
!--- update atomic configuration ---------------------------------------

    if( ifmd == 1 ) then

        !-----for hibrid structural and cell optimization
        if( ioptmze >= 0 .and. ioptmze <= 9 .and. ioptmze_cell >= 0 ) then
            if( nstepini < 0 ) nstepini = nstep
            if( lhybridopt ) then
                !---optimize structure first
                nstep1 = nstep_hybrid
                nstep2 = nstep_hybrid_cell
            else
                !---optimize cell first
                nstep2 = nstep_hybrid
                nstep1 = nstep_hybrid_cell
            end if
            if( nstep - nstepini == 0 .or. nstep - nstepini == nstep1 ) then
                lqninitial      = .true.
                lqninitial_cell = .true.
            end if
            if( nstep - nstepini < nstep1 ) then
                lrunopt      = lhybridopt
                lrunopt_cell = .not.lhybridopt
            else
                lrunopt      = .not.lhybridopt
                lrunopt_cell = lhybridopt
                if( nstep - nstepini == nstep1+nstep2-1 ) nstepini = nstep + 1
            end if
            nstep0 = nstep - 1
        else
            lrunopt      = .true.
            lrunopt_cell = .true.
            nstep0 = 0
        end if

        lrunoptif: if( lrunopt ) then
            !-----structural optimization
            if( ioptmze == 0 ) then
                !-------Conjugate-gradient minimization
                call updtcg( nfile, myid, nodes,  &
& x, is, nmd, ntmd, frc, watom_md, h, dtmd, pixcg, fnorm, acon,  &
& lfixion )
            else if( ioptmze == 2 ) then
                if( ibfgsclear > 0 ) then
                    lclear = mod(nstep-nstep0,ibfgsclear) == 0
                else
                    lclear = .false.
                end if
                !-----for quasi-Newton methods with BFGS formula
                call updtqn( nfile, myid, nodes,  &
& x, is, nmd, ntmd, frc, watom_md, h, dtmd, pixcg, fnorm, acon,  &
& lfixion,  &
& deltax, deltaf, hessian, heigenvec, heigenval, scalex, scaleg,  &
& lqninitial, dist_max, lclear, lhessian, qnstabi, gammamin )
!            else if( ioptmze == 3 ) then
                !-----for RFO saddle poins search
!                call updtrfo( nfile, myid, nodes,  &
!& x, is, nmd, ntmd, frc, watom_md, h, dtmd, pixcg, fnorm, acon,  &
!& lfixion,  &
!& deltax, deltaf, hessian, heigenvec, heigenval, scalex, scaleg,  &
!& lqninitial, dist_max )
!            else if( ioptmze == 10 ) then
                !-----for Harmonic-mode analysis
!                call updthma( nfile, myid, nodes,  &
!& x, is, nmd, ntmd, frc, watom_md, h, pixcg, fnorm, acon,  &
!& deltax, deltaf, hessian, heigenvec, heigenval, scalex, scaleg,  &
!& hmadisp, nhmaord, lcentdif, ierror )
!            else if( ioptmze == 11 ) then
!                !-----for phonon-dispersion calculation
!                call updtphonon( nfile, myid, nodes,  &
!& x, is, nmd, ntmd, frc, watom_md, h, pixcg, fnorm, acon,  &
!& deltax, deltaf, hessian, heigenvec, heigenval, scalex, scaleg,  &
!& hmadisp, nhmaord, lcentdif, ierror )
!            else if( ioptmze == 20 ) then
!                !-----calculation with given atomic displacement
!                call updtdisp( nfile, myid, nodes,  &
!& x, is, nmd, h )
            end if
        end if lrunoptif

        lrunoptcellif: if( lrunopt_cell ) then
            !-----supercell optimization
            if( ioptmze_cell == 0 ) then
                !-------Conjugate-gradient minimization
                call updtcg_cell( nfile, myid, nodes,  &
& h, volume, pint, pext, dtcellcg, cellfnorm, irstrct, irstrct_sub, lcell_rstrct )
            else if( ioptmze_cell == 2 ) then
                if( iclearcellh > 0 ) then
                    lclearcellh = mod(nstep-nstep0,iclearcellh) == 0
                else
                    lclearcellh = .false.
                end if
                !-----for quasi-Newton methods with BFGS formula
                call updtqn_cell( nfile, myid, nodes,  &
& h(1,1,0), h(1,1,1), volume, pint, pext, dtcellcg, cellfnorm,  &
& hessian_cell, savefrc_cell, lqninitial_cell, irstrct, irstrct_sub, lcell_rstrct,  &
& lclearcellh, lhessian_cell, qnstabi, gammamin )
            end if
        end if lrunoptcellif

    else if( ifmd >= 2 ) then

        !-------Heat bath dynamics if ifmd=3
        if( ifmd == 3 )  &
&           call nhcint( nfile, myid, nodes,  &
& nnos, nresn, nyosh, gnkt, gkt,  &
& xlogs, vlogs, glogs, qmass, wdti, wdti2, wdti4, wdti8,  &
& x(1,1), is, nmd, ntmd, fack, h )

        !-------NPT dynamics if ifmd=4
        if( ifmd == 4 )  &
&           call nhcpfullint( nfile, myid, nodes,  &
& nnos, nresn, nyosh, gnkt, gkt,  &
& xlogs, vlogs, glogs, qmass, wdti, wdti2, wdti4, wdti8,  &
& x(1,1), is, nmd, ntmd, fack, fack2, h, volume, pint, pext,  &
& vu, bmass, onf, gnd2kt, irstrct, irstrct_sub, lcell_rstrct )

        !-------Heat bath dynamics if ifmd=5
!        if( ifmd == 5 ) then
!            call nhcint_atom2( nfile, myid, nodes,  &
!& lathermo, nnos_atom, nresn_atom, nyosh_atom, g3kt, gkt,  &
!& xlogs_atom, vlogs_atom, glogs_atom, qmass_atom,  &
!& wdti_atom, wdti2_atom, wdti4_atom, wdti8_atom,  &
!& x(1,1), is, nmd, ntmd, fack, h )
!
!            if( lgthermo )  &
!&               call nhcint_atom( nfile, myid, nodes,  &
!& lathermo, nnos, nresn, nyosh, gnkt, gkt,  &
!& xlogs, vlogs, glogs, qmass, wdti, wdti2, wdti4, wdti8,  &
!& x(1,1), is, nmd, ntmd, fack, h )
!        end if

        !-------GLE thermostat dynamics if ifmd=6
!        if( ifmd == 6 ) then
!            call glethermo( nfile, myid, nodes,  &
!& lathermo, x(1,1), is, nmd, ntmd, ntot, fack, watom_md, dtmd, h )

!            if( lgthermo )  &
!&               call nhcint_atom( nfile, myid, nodes,  &
!& lathermo, nnos, nresn, nyosh, gnkt, gkt,  &
!& xlogs, vlogs, glogs, qmass, wdti, wdti2, wdti4, wdti8,  &
!& x(1,1), is, nmd, ntmd, fack, h )
!        end if

        !-------MSST dynamics if ifmd=10
        if( ifmd == 10 )  &
&           call msstint( nfile, myid, nodes,  &
& x(1,1), is, nmd, ntmd, fack2, h, volume, xmsst, vxmsst,  &
& wdti(1), wdti2(1), wdti4(1), wdti8(1), vu, pint, pext, bmass,  &
& volum0, alpmatrix, betmatrix, totmass, shockspeed )

        !-----First half-step kick, v(t+Dt/2)
        call vkick( nfile, myid, nodes,  &
& x(1,1), frc, is, nmd, ntmd, lfixion )
!        call fixdisp_v( nfile, myid, nodes,  &
!& x(1,0), x(1,1), is, nmd, ntmd, ldispion, dvector, h )

        !-----Coordinate update, x(t+Dt)
        if( ifmd == 4 ) then

            call nhc_position_updt( nfile, myid, nodes,  &
& x(1,0), x(1,1), is, nmd, ntmd, lfixion, h, volume, dtmd, xu, vu,  &
& sxog, syog, szog, irstrct, irstrct_sub, lcell_rstrct )

        else if( ifmd == 10 ) then

            call msst_position_updt( nfile, myid, nodes,  &
& x(1,0), x(1,1), is, nmd, ntmd, lfixion, xu, vu,  &
& sxog, syog, szog, h, volume, dtmd,  &
& xmsst, vxmsst, alpmatrix, betmatrix, totmass )

        else

            !---projected velocity Verlet for structural optimization
            if( lpvv ) call pvv( nfile, myid, nodes,  &
& x(1,1), frc, is, nmd, ntmd, lfixion )

!            if( ifmd == 5 .and. lmomzero ) then
!                !---make total momentum zero
!                call momzero2( nfile, myid, nodes,  &
!& x(1,1), is, nmd, ntmd, ntot, watom_md )
!            end if

!            if( ncbonds == 0 ) then
                call updtps( nfile, myid, nodes,  &
& x(1,0), x(1,1), is, nmd, ntmd, lfixion )
!              else
!                !---constrained MD
!                call updtps_constraint( nfile, myid, nodes,  &
!& x(1,0), x(1,1), is, nmd, ntmd, lfixion, h, acon, &
!& ncbatm1, ncbatm2, cblength, cblambda, ncbonds )
!            end if

        end if
!        call fixdisp( nfile, myid, nodes,  &
!& x(1,0), x(1,1), is, nmd, ntmd, ldispion, dvector, h )

    end if
!-----------------------------------------------------------------------

!--- Sends moved-out atoms to neighbor nodes and receives moved-in atoms
!--- from neighbor nodes.
    if( ifmd >= 1 ) then
 !       if( lQMMDrun .or. lpureMD ) then

            !-----check atomic displacements
            !-----set imts = 0, if neighbor list must be updated
!            call chkdis( nfile, myid, nodes,  &
!& x(1,0), is, nmd, ntmd, h, lfixion, xrec, imts )

!            if( imts == 0 ) then
!                if( .not.lstat ) then
!                    call bamove( nfile, myid, nodes,  &
!& x(1,0), x(1,1), is, nmd, nmd_alloc, nmdmax,  &
!& leinstein, x_vmd, nnos_atom, xlogs_atom, vlogs_atom, qmass_atom, ifmd )
!                else
!                    call bamove_v0( nfile, myid, nodes,  &
!& x(1,0), x(1,1), is, nmd, nmd_alloc, nmdmax,  &
!& leinstein, x_vmd, nnos_atom, xlogs_atom, vlogs_atom, qmass_atom, ifmd )
!                end if
!            end if

!        else

            !-----periodic boundary conditions
!            if( .not.lstat ) then
                call updtps_pbc( nfile, myid, nodes,  &
& x(1,0), is, nmd, ntmd )
!            else
!                call updtps_pbc_v0( nfile, myid, nodes,  &
!& x(1,0), is, nmd, ntmd )
!            end if

!        end if
    end if

end if


return
end subroutine




subroutine remd_semifinal( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!     update atomic configuration
!-----------------------------------------------------------------------
use constants
use remd_param
use remd_param_atom
use remd_variables
use remd_atom
use remd_qmmd
implicit none
integer :: nfile(*), myid, nodes
integer :: ierror

!-----declare local variables
real*8  :: etotal, epot_total, ekin = 0.d0, hmltan, temp
real*8  :: bathpe, bathke, boxpe, boxke, frcmax
real*8  :: bathpe_atom, bathke_atom, ekindif
real*8  :: qmmd_energy(3)
integer :: i, noutqmmd, digit
character(50) :: frmt


!-----total potential energy
epot_total = epot + epotQM - epotCL


!-----force QM correction
if( .not.lpureMD ) then
    call correct_force( nfile, myid, nodes,  &
& natom, ntype, nhk1, nhk2, frcQM, lterminator, nmd, frc,  &
& ind_qm1_md, ind_HSQM_qm, ind_HSQM_md, alpha_HSQM, nHSQM, nHSQMx )

    call correct_stress( nfile, myid, nodes,  &
& pintlr, pintsr, pint, pit1, strQM )
end if


!-----force MD-cluster correction
!if( lQMMDrun ) then
!    call correct_force( nfile, myid, nodes,  &
!& nmd_cluster, ntmd_cluster, nhk1_cluster, nhk2_cluster, frcCL,  &
!& lMDterminator, nmd, frc,  &
!& ind_md1_md, ind_HSMD_qm, ind_HSMD_md, alpha_HSMD, nHSMD, nHSMDx )
!end if


!-----Check force by crystal symmetry
!call symfrc( nfile, myid, nodes,  &
!& nmd, frc )


if( ifmd == 1 ) then
    !-----get maximum norm of atomic forces
    call max_force( nfile, myid, nodes,  &
& nmd, ntmd, frc, is, lfixion, frcmax )
end if


!if( lvmd ) then
!    !-----for virtual MD for thermodynamic integration
!      epot_total = dlambda_vmd*epot_total   + (1.d0-dlambda_vmd)*epotref
!    frc(1:3*nmd) = dlambda_vmd*frc(1:3*nmd) + (1.d0-dlambda_vmd)*frcref(1:3*nmd)
!
!     pintlr(1:3,1:3) = dlambda_vmd*pintlr(1:3,1:3) + (1.d0-dlambda_vmd)*pintlrref(1:3,1:3)
!     pintsr(1:3,1:3) = dlambda_vmd*pintsr(1:3,1:3) + (1.d0-dlambda_vmd)*pintsrref(1:3,1:3)
!       pint(1:3,1:3) = pintlr(1:3,1:3) + pintsr(1:3,1:3)
!       pit1(1:3,1:3) =   pint(1:3,1:3)
!end if


!-----Calculate normalized accelerations, A
call normalized_a( nfile, myid, nodes,  &
& nmd, ntmd, frc, is, h, acon )


if( ifmd == 1 ) then

    !-----Cell optimization = enthalpy optimazation
    if( ioptmze_cell >= 0 ) then
         boxpe=(pext(1,1)+pext(2,2)+pext(3,3))*volume/3d0
         epot_total = epot_total + boxpe
    end if
    if( ioptmze >= 0 .and. ioptmze <= 9 .or. ioptmze_cell >= 0 ) then
    if( lhessian .or. lhessian_cell ) then
        ierror = -1
        if( myid.eq.0 ) then
            write(nfile(1),*) ' *** Negative eigenvalues of Hessian !'
            write(nfile(1),*) '     energy convergence =',  &
&                             epot_total - prevene
            write(nfile(1),*) '     force  convergence =', frcmax
            write(nfile(2),*) ' *** Negative eigenvalues of Hessian !'
            write(nfile(2),*) '     energy convergence =',  &
&                             epot_total - prevene
            write(nfile(2),*) '     force  convergence =', frcmax
        end if
        return
    end if

    !-----check energy convergence
    if( abs(epot_total-prevene) < tol_energy .or.  &
&                        frcmax < tol_force .and. ioptmze_cell < 0 ) then
        ierror = -1
        if( myid.eq.0 ) then
            write(nfile(1),*) ' *** CG minimization complete'
            write(nfile(1),*) '     energy convergence =',  &
&                             epot_total - prevene
            write(nfile(1),*) '     force  convergence =', frcmax
            write(nfile(2),*) ' *** CG minimization complete'
            write(nfile(2),*) '     energy convergence =',  &
&                             epot_total - prevene
            write(nfile(2),*) '     force  convergence =', frcmax
        end if
        return
    end if

    if( epot_total.gt.prevene ) then
!              ierror = 99
        if( myid.eq.0 ) then
            write(nfile(1),*) ' *** CG minimization error',  &
&                             ', still continue ...'
            write(nfile(2),*) ' *** CG minimization error',  &
&                             ', still continue ...'
        end if
!              return
        fnorm(1:nmd) = 0.d0
    end if
    end if
    prevene = epot_total

else if( ifmd >= 2 ) then

    if( nstep >= nstep_ini+1 ) then
        !-------Second half-step kick, v(t+Dt)
!        if( ncbonds == 0 ) then
            call vkick( nfile, myid, nodes,  &
& x(1,1), frc, is, nmd, ntmd, lfixion )
!          else
!            !---constrained MD
!            call vkick_constraint( nfile, myid, nodes,  &
!& x(1,0), x(1,1), frc, is, nmd, ntmd, lfixion, h, acon, &
!& ncbatm1, ncbatm2, cblambdav, ncbonds )
!            end if
!        call fixdisp_v( nfile, myid, nodes,  &
!& x(1,0), x(1,1), is, nmd, ntmd, ldispion, dvector, h )

        !-------Heat bath dynamics if ifmd=3
        if( ifmd == 3 )  &
&           call nhcint( nfile, myid, nodes,  &
& nnos, nresn, nyosh, gnkt, gkt,  &
& xlogs, vlogs, glogs, qmass, wdti, wdti2, wdti4, wdti8,  &
& x(1,1), is, nmd, ntmd, fack, h )

        !-------NPT dynamics if ifmd=4
        if( ifmd == 4 )  &
&           call nhcpfullint( nfile, myid, nodes,  &
& nnos, nresn, nyosh, gnkt, gkt,  &
& xlogs, vlogs, glogs, qmass, wdti, wdti2, wdti4, wdti8,  &
& x(1,1), is, nmd, ntmd, fack, fack2, h, volume, pint, pext,  &
& vu, bmass, onf, gnd2kt, irstrct, irstrct_sub, lcell_rstrct )

        !-------Heat bath dynamics if ifmd=5
!        if( ifmd == 5 ) then
!            if( lgthermo )  &
!&               call nhcint_atom( nfile, myid, nodes,  &
!& lathermo, nnos, nresn, nyosh, gnkt, gkt,  &
!& xlogs, vlogs, glogs, qmass, wdti, wdti2, wdti4, wdti8,  &
!& x(1,1), is, nmd, ntmd, fack, h )

!            call nhcint_atom2( nfile, myid, nodes,  &
!& lathermo, nnos_atom, nresn_atom, nyosh_atom, g3kt, gkt,  &
!& xlogs_atom, vlogs_atom, glogs_atom, qmass_atom,  &
!& wdti_atom, wdti2_atom, wdti4_atom, wdti8_atom,  &
!& x(1,1), is, nmd, ntmd, fack, h )
!        end if

        !-------GLE thermostat dynamics if ifmd=6
!        if( ifmd == 6 ) then
!            if( lgthermo )  &
!&               call nhcint_atom( nfile, myid, nodes,  &
!& lathermo, nnos, nresn, nyosh, gnkt, gkt,  &
!& xlogs, vlogs, glogs, qmass, wdti, wdti2, wdti4, wdti8,  &
!& x(1,1), is, nmd, ntmd, fack, h )

!            call glethermo( nfile, myid, nodes,  &
!& lathermo, x(1,1), is, nmd, ntmd, ntot, fack, watom_md, dtmd, h )
!        end if

        !-------MSST dynamics if ifmd=10
        if( ifmd == 10 )  &
&           call msstint( nfile, myid, nodes,  &
& x(1,1), is, nmd, ntmd, fack2, h, volume, xmsst, vxmsst,  &
& wdti(1), wdti2(1), wdti4(1), wdti8(1), vu, pint, pext, bmass,  &
& volum0, alpmatrix, betmatrix, totmass, shockspeed )
    end if

    !--- Kinetic E.
    call getekin( nfile, myid, nodes,  &
& ekin, x(1,1), is, nmd, ntmd, fack, h )
    temp = 2d0*ekin/3d0/ntot(0)

    !--- Heat-bath potential & kinetic energies
    if( ifmd == 3 )  &
&       call getbathe( nfile, myid, nodes,  &
& bathpe, bathke, nnos, gnkt, gkt, xlogs, vlogs, glogs, qmass )

    !--- Heat-bath & barostat potential & kinetic energies
    if( ifmd == 4 )  &
&       call getbathboxe( nfile, myid, nodes,  &
& bathpe, bathke, boxpe, boxke,  &
& nnos, gnkt, gkt, xlogs, vlogs, glogs, qmass,  &
& bmass, gnd2kt, pext, h, volume )

    !--- Heat-bath potential & kinetic energies
!    if( ifmd == 5 ) then
!        if( lgthermo ) then
!            call getbathe( nfile, myid, nodes,  &
!& bathpe, bathke, nnos, gnkt, gkt, xlogs, vlogs, glogs, qmass )
!        else
!            bathpe = 0.d0
!            bathke = 0.d0
!        end if
!
!        call getbathe_atom( nfile, myid, nodes,  &
!& bathpe_atom, bathke_atom, lathermo, nnos_atom, g3kt, gkt,  &
!& xlogs_atom, vlogs_atom, qmass_atom, is, nmd, ntmd )
!    end if

    !--- Heat-bath potential & kinetic energies
!    if( ifmd == 6 ) then
!        if( lgthermo ) then
!            call getbathe( nfile, myid, nodes,  &
!& bathpe, bathke, nnos, gnkt, gkt, xlogs, vlogs, glogs, qmass )
!        else
!            bathpe = 0.d0
!            bathke = 0.d0
!        end if
!
!        ekindif = 0d0
!!        call get_gle_ekindif( ekindif )
!    end if

    !--- Barostat potential & kinetic energies
    if( ifmd == 10 )  &
&       call getmsstboxe( nfile, myid, nodes,  &
& boxpe, boxke, h, volume, bmass, pext, &
& vxmsst, volum0, alpmatrix, totmass, shockspeed )


    !--- Total E.
    etotal = ekin + epot_total

    !--- conserved quantity
    if( ifmd == 2 ) then
        hmltan = etotal
    else if( ifmd == 3 ) then
        hmltan = etotal + bathpe + bathke
    else if( ifmd == 4 ) then
        hmltan = etotal + bathpe + bathke + boxpe + boxke
!    else if( ifmd == 5 ) then
!        hmltan = etotal + bathpe + bathke + bathpe_atom + bathke_atom
!    else if( ifmd == 6 ) then
!        hmltan = etotal + bathpe + bathke - ekindif
    else if( ifmd == 10 ) then
        hmltan = etotal + boxpe + boxke
    end if

    !--- check temperature
    if( liscale .and. mod(nstep,iscstp).eq.0 .and.  &
&       nstep > nstep_ini ) then
        isccnt = isccnt + 1
        if( isccnt.ge.iscnum ) liscale = .false.

        !-----Make total momentum zero
        if( lmomzero ) then
            call momzero( nfile, myid, nodes,  &
& x(1,0), x(1,1), is, nmd, ntmd, ntot, watom_md, h,  &
& sxog, syog, szog )
        end if

        !-----velocity scaling
        call vscale( nfile, myid, nodes,  &
& treq, ekin, x(1,1), is, nmd, ntmd, ntot )

        !-----scaling the thermostats
        if( ifmd == 3 .or. ifmd == 4 )  &
&           call nhc_scale( nfile, myid, nodes,  &
& nnos, gnkt, gkt, xlogs, vlogs, glogs, qmass, h )

!        if( ifmd == 5 ) then
!            if( lgthermo )  &
!&               call nhc_scale( nfile, myid, nodes,  &
!& nnos, gnkt, gkt, xlogs, vlogs, glogs, qmass, h )
!
!            call nhc_scale_atom( nfile, myid, nodes,  &
!& lathermo, nnos_atom, g3kt, gkt, xlogs_atom, vlogs_atom, qmass_atom, is, nmd, ntmd )
!        end if

!        if( ifmd == 6 ) then
!            if( lgthermo )  &
!&               call nhc_scale( nfile, myid, nodes,  &
!& nnos, gnkt, gkt, xlogs, vlogs, glogs, qmass, h )
!
!!            call gle_scale( nfile, myid, nodes,  &
!!& lathermo, is, nmd, ntmd, ntot )
!        end if

        !-----scaling the barostats
        if( ifmd == 10 )  &
&           call msst_scale( nfile, myid, nodes, vxmsst, h )
    end if

    if( ifmd == 10 ) then
    !--- clear barostat velocity
    if( lmsstscale .and. mod(nstep,msstscstp).eq.0 .and.  &
&       nstep > nstep_ini ) then
        msstsccnt = msstsccnt + 1
        if( msstsccnt.ge.msstscnum ) lmsstscale = .false.
        call msst_scale( nfile, myid, nodes, vxmsst, h )
    end if
    end if

end if


!if( lQMMDrun ) then
!    noutqmmd = 3
!    qmmd_energy(1) = epot
!    qmmd_energy(2) = epotQM
!    qmmd_energy(3) = epotCL
!  else
    noutqmmd = 0
    qmmd_energy(1) = 0.d0
    qmmd_energy(2) = 0.d0
    qmmd_energy(3) = 0.d0
!end if

if( myid.eq.0 .and. mod(nstep,ioskip) == 0 ) then

    digit = log(dble(nstep)+0.5d0)/log(10.d0) + 1.d0

    if( ifmd <= 1 ) then
        if( ioptmze_cell < 0 ) then
            frmt = '(i7,es18.10,3es16.8)'
            if( digit > 7 ) call get_frmt( frmt, digit )
            write(nfile(1),frmt) nstep, epot_total
            write(nfile(2),frmt) nstep, epot_total
            write(nfile(3),frmt) nstep, epot_total,  &
&                  ( qmmd_energy(i), i = 1, noutqmmd )
          else
            !-----epot_total = enthalpy
            frmt = '(i7,3es18.10,3es16.8)'
            if( digit > 7 ) call get_frmt( frmt, digit )
            write(nfile(1),frmt) nstep, epot_total, epot_total-boxpe, boxpe
            write(nfile(2),frmt) nstep, epot_total, epot_total-boxpe, boxpe
            write(nfile(3),frmt) nstep, epot_total, epot_total-boxpe, boxpe, &
&                         ( qmmd_energy(i), i = 1, noutqmmd )
        end if
        if( ifmd == 1 ) then
            write(nfile(1),'(6x,es18.10)') frcmax
            write(nfile(2),'(6x,es18.10)') frcmax
        end if
    else
        frmt = '(i7,es18.10,2es16.8,1x,a1,f11.4,a2)'
        if( digit > 7 ) call get_frmt( frmt, digit )
        write(nfile(1),frmt) nstep, hmltan, epot_total, ekin, '(', temp*tempau, ' )'
        write(nfile(2),frmt) nstep, hmltan, epot_total, ekin, '(', temp*tempau, ' )'

        frmt = '(i7,es18.10,2es16.8,f11.4,8es16.8)'
        if( digit > 7 ) call get_frmt( frmt, digit )
        if( ifmd == 2 ) then
            write(nfile(3),frmt) nstep, hmltan, epot_total, ekin, temp*tempau,  &
&                  ( qmmd_energy(i), i = 1, noutqmmd )
        else if( ifmd == 3 ) then
            write(nfile(3),frmt) nstep, hmltan, epot_total, ekin, temp*tempau,  &
&                    etotal, bathpe, bathke,  &
&                  ( qmmd_energy(i), i = 1, noutqmmd )
        else if( ifmd == 4 ) then
            write(nfile(3),frmt) nstep, hmltan, epot_total, ekin, temp*tempau,  &
&                    etotal, bathpe, bathke, boxpe, boxke,  &
&                  ( qmmd_energy(i), i = 1, noutqmmd )
!        else if( ifmd == 5 ) then
!            write(nfile(3),frmt) nstep, hmltan, epot_total, ekin, temp*tempau,  &
!&                    etotal, bathpe, bathke, bathpe_atom, bathke_atom,  &
!&                  ( qmmd_energy(i), i = 1, noutqmmd )
!        else if( ifmd == 6 ) then
!            write(nfile(3),frmt) nstep, hmltan, epot_total, ekin, temp*tempau,  &
!&                    etotal, bathpe, bathke, ekindif,  &
!&                  ( qmmd_energy(i), i = 1, noutqmmd )
        else if( ifmd == 10 ) then
            write(nfile(3),frmt) nstep, hmltan, epot_total, ekin, temp*tempau,  &
&                    etotal, boxpe, boxke,  &
&                  ( qmmd_energy(i), i = 1, noutqmmd )
        end if
    end if
end if

!if( lvmd ) then
!    if( myid.eq.0 ) then !.and. mod(nstep,ioskip) == 0 ) then
!        !-----for virtual MD for thermodynamic integration
!        frmt = '(i7,2es16.8)'
!        if( digit > 7 ) call get_frmt( frmt, digit )
!        write(nfile(24),frmt) nstep, epot+epotQM, epotref
!    end if
!end if


!if( lheat .and. mod(nstep,nskip_heat) == 0 ) then
!    !---heat flux vectors
!    call remd_heat_flux( nfile, myid, nodes,  &
!& nstep, digit, x(1,1), is, nmd, dtmd, ntmd, h, fack, volume,  &
!& wepot, wstrs, nmdmax__, ntot )
!end if


!if( ltotmom .and. mod(nstep,nskip_totmom) == 0 ) then
!    !---reduced total momentum
!    call totmom( nfile, myid, nodes,  &
!& nstep, digit, x(1,1), is, nmd, ntmd, ntot, watom_md )
!end if


!-----output atomic configurations
call remd_wrgeom( nfile, myid, nodes, iogpsz,  &
& nstep, ifmd, x(1,0), x(1,1), is, nmd, dtmd, ntmd, frc, h, acon, lfixion,  &
& sxog, syog, szog, pit1, fack, volume,  &
& lmdout, locoor, lovelo, loforc, ioskip, ioskipcoor, ioskipvelo, ioskipforc,  &
& ncbonds, lstress, nskip_stress, nskip_stress_out, latomic, nskip_atomic,  &
& wepot, wstrs, wstrsk, nmdmax__, ntot,  &
& lrestrict_area, xio_min, xio_max, yio_min, yio_max, zio_min, zio_max,  &
& ierror )


!-----output atomic configurations for QM & MD-cluster atoms
!if( lQMMDrun ) then
!
!    !-----for QM atoms
!    call remd_wrgeom_qmmd( nfile, myid, nodes,  &
!& nstep, ifmd, dtmd, ratm, vatm, frcQM, nhk, ntype, natom,  &
!& lterminator, h,  &
!& lmdout_qm, locoor_qm, lovelo_qm, loforc_qm, ioskip_qm,  &
!& nfile(10), nfile(11), nfile(12), ierror )
!
!    !-----for MD-cluster atoms
!    call remd_wrgeom_qmmd( nfile, myid, nodes,  &
!& nstep, ifmd, dtmd, x_cluster(1,0), x_cluster(1,1), frcCL,  &
!& ntot_cluster(1), ntmd_cluster, nmd_cluster, lMDterminator, h,  &
!& lmdout_cl, locoor_cl, lovelo_cl, loforc_cl, ioskip_cl,  &
!& nfile(13), nfile(14), nfile(15), ierror )
!
!end if


!----- Statistical calculation
!if( lstat ) call remd_cal_stat( nfile, myid, nodes,  &
!& nstep, x(1,0), x(1,1), is, nmd, npb, nmd_alloc, ntmd, h,  &
!& npx, npy, npz, anxi, anyi, anzi, myx, myy, myz, sxog, syog, szog )


return
end subroutine




subroutine remd_save( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!     save data in MD nodes
!-----------------------------------------------------------------------
use constants
use remd_param
use remd_param_atom
use remd_variables
use remd_atom
implicit none
integer :: nfile(*), myid, nodes
integer :: ierror


if( ifmd.eq.1 ) then
    nstepCG = nstep
else if( ifmd.ge.2 ) then
    nstepMD = nstep
end if

!if( ifmd == 1 ) then
    !-----Harmonic-mode analysis
!    if( ioptmze == 10 ) then
        !---Calculate Haromic mode
!        call hma_fin( nfile, myid, nodes,  &
!& x, is, nmd, ntmd, frc, watom_md, h, pixcg, fnorm, acon,  &
!& deltax, deltaf, hessian, heigenvec, heigenval, scalex, scaleg,  &
!& hmadisp, nhmaord, lcentdif )
!
!        !---Do not save, and return
!        return
!    end if

    !-----Phonon-dispersion calculation
!    if( ioptmze == 11 ) then
!        call phonon_fin( nfile, myid, nodes,  &
!& x, is, nmd, ntmd, frc, watom_md, h, pixcg, fnorm, acon,  &
!& deltax, deltaf, hessian, heigenvec, heigenval, scalex, scaleg,  &
!& hmadisp, nhmaord, lcentdif )
!
!        !---Do not save, and return
!        return
!    end if

    !-----calculation with given atomic displacement
!    if( ioptmze == 20 ) then
!
!        !---Do not save, and return
!        return
!    end if
!end if

!-----save atomic configuration to files, if lsave == .true.
call remd_svatoms( nfile, myid, nodes, npx, npy, npz, iogpsz,  &
& fname_toc, fname_mts, lsave, ifmd, nstepCG, nstepMD, dtmd,  &
& lQMMDrun, lpureMD,  &
& natom, ntype, nhk1, nhk2, ratm, lterminator, natomx, ntypex,  &
& nmd_cluster, ntmd_cluster, nhk1_cluster, nhk2_cluster, x_cluster,  &
& lMDterminator, nmdx_cluster, ntmdx_cluster,  &
& nmd, is, x, h, frc, pixcg, fnorm, nmd_alloc, nmdmax_,  &
& nnos, xlogs, vlogs, pintlr, pintsr, ierror,  &
& ioptmze, ioptmze_cell, xmsst, vxmsst, h0, nnos_atom, xlogs_atom, vlogs_atom )


call remd_close_files( nfile, myid, nodes, ierror )


!----- Statistical calculation
!if( lstat ) call remd_fin_stat( nfile, myid, nodes,  &
!& ntmd, ntot, dtmd )


return
end subroutine




subroutine remd_force( nfile, myid, nodes )
!-----------------------------------------------------------------------
!     calculate classical forces
!-----------------------------------------------------------------------
use constants
use remd_param
use remd_param_atom
use remd_variables
use remd_atom
use remd_qmmd
implicit none
integer :: nfile(*), myid, nodes

!-----declare local variables
logical :: lcstress, lcatomic


!-----zero clear classical energy and forces
call get_force_null( nfile, myid, nodes,  &
& epot, x, is, nmd, ntmd, frc, h, acon, pintlr, pintsr )

epotCL = 0.d0 


!-----classical energy and force for entire system
!if( lQMMDrun .or. lpureMD ) then
!    !-----check if stresses are calculated
!    lcstress = lstress .and. mod(nstep,nskip_stress) == 0
!    !-----check if atomic stresses & energies are calculated
!    lcatomic = latomic .and. mod(nstep,nskip_atomic) == 0  &
!&         .or. lheat   .and. mod(nstep,nskip_heat)   == 0
!
!    call get_force( nfile, myid, nodes,  &
!& epot, x, is, nmd, ntmd, frc, h, acon,  &
!& sxog, syog, szog, anxi, anyi, anzi, rc_buffer,  &
!& npb, nmd_buffer, nmd_buffer_1d, nmd_alloc, imts, lspr, nprmax,  &
!& is2ist, ist, xrec, pintlr, pintsr, pint, pit1, volume,  &
!& lcstress, lcatomic, wepot, wstrs, nmdmax__, ntot, lefield, efield )
!end if


!-----classical energy and force for cluster system
!if( lQMMDrun ) then
!
!    !-----get MD-cluster atom positions in cluster
!    call getMDCLpositions( nfile, myid, nodes,  &
!& natom, ntype, nhk1, nhk2, ratm, lterminator, ind_qm1_md,  &
!& nmd_cluster, ntmd_cluster, nhk1_cluster, nhk2_cluster, x_cluster,  &
!& lMDterminator, ind_md1_md )
!
!    !-----get MD terminator atom (HSQM) positions
!    call getHSQMpositions( nfile, myid, nodes,  &
!& nmd_cluster, ntmd_cluster, nhk1_cluster, nhk2_cluster, x_cluster,  &
!& lMDterminator, nmd, is, x, h,  &
!& ind_HSMD_qm, x_HSMD_qm,  &
!& ind_HSMD_md, x_HSMD_md, alpha_HSMD, nHSMD, nHSMDx,  &
!& sxog, syog, szog )
!
!    call get_CLforce( nfile, myid, nodes,  &
!& epotCL, x_cluster, is_cluster, nmd_cluster, ntmd_cluster, frcCL,  &
!& h, rc_buffer, imts, lspr_cluster, nprmax,  &
!& is2ist_cluster, ist_cluster )
!
!end if


!-----energies and forces from soft walls
!if( nwalls > 0 ) then
!    call get_force_walls( nfile, myid, nodes,  &
!& epot, x, is, nmd, ntmd, frc, h, nmd_alloc, sxog, syog, szog,  &
!& nwalls, wallp, wallv, wallf, nwallp )
!end if


!-----gravitational field
!if( lgravi ) then
!    call get_force_gravity( nfile, myid, nodes,  &
!& is, nmd, ntmd, frc, watom_md, nmd_alloc, gravdir )
!end if


!-----virtual MD for thermodynamic integration
!if( lvmd ) then
!    if( .not.lidealref ) then
!        !-----check if stresses are calculated
!        lcstress = lstress .and. mod(nstep,nskip_stress) == 0
!!        lcstress = .false.
!        lcatomic = .false.
!        call get_force_vmd( nfile, myid, nodes,  &
!& epotref, x, is, nmd, ntmd, frcref, h, acon,  &
!& sxog, syog, szog, anxi, anyi, anzi, rc_buffer,  &
!& npb, nmd_buffer, nmd_buffer_1d, nmd_alloc, lspr, nprmax,  &
!& is2ist, ist, xrec, pintlrref, pintsrref, pint, pit1, volume,  &
!& lcstress, lcatomic, wepot, wstrs, nmdmax__, ntot, lQMMDrun, lpureMD )
!    else if( leinstein ) then
!        !---Einstein solid
!        call get_force_vmd_einstein( nfile, myid, nodes,  &
!& epotref, x, is, nmd, ntmd, frcref, h, nmd_alloc, watom_md,  &
!& omega_vmd, x_vmd, nvmd )
!    end if
!end if


imts = 1


return
end subroutine




module heatbath_backup
!-----------------------------------------------------------------------
!     backup thermostat variables
!-----------------------------------------------------------------------
implicit none

integer :: nnosbk = 0
real*8,  allocatable, dimension(:) :: xlogsbk    ! Heat-bath coordinates
real*8,  allocatable, dimension(:) :: vlogsbk    ! Heat-bath velocities

save

end module




subroutine remd_rdatoms( nfile, myid, nodes, npx, npy, npz, iogpsz,  &
& fname_toc, fname_mts, lstart, ifmd, nstepCG, nstepMD, dt,  &
& lQMMDrun, lpureMD,  &
& natom, ntype, nhk1, nhk2, ratm, lterminator, natomx, ntypex,  &
& nmd_cluster, ntmd_cluster, nhk1_cluster, nhk2_cluster, x_cluster,  &
& lMDterminator, nmdx_cluster, ntmdx_cluster,  &
& n, is, x, h, frc, pixcg, fnorm, nmdx, nmdmax, nmdmax_,  &
& nnos, xlogs, vlogs, lbath_clear, pintlr, pintsr, pint, ierror,  &
& ioptmze, ioptmze_cell, lqninitial, xmsst, vxmsst, h0,  &
& nnos_atom, xlogs_atom, vlogs_atom )
!-----------------------------------------------------------------------
!     read atomic configuration from files, if lstart == .true.
!-----------------------------------------------------------------------
use heatbath_backup
implicit none
integer :: nfile(*), myid, nodes
integer :: npx, npy, npz
integer :: iogpsz
character(*):: fname_toc
character(*):: fname_mts
logical :: lstart
integer :: ifmd
integer :: nstepCG, nstepMD
real*8  :: dt
logical :: lQMMDrun
logical :: lpureMD

!-----for QM atoms
integer :: natom
integer :: ntype
integer :: natomx
integer :: ntypex
integer, dimension(ntypex) :: nhk1, nhk2
real*8,  dimension(3,natomx) :: ratm
logical, dimension(ntypex) :: lterminator

!-----for MD-cluster atoms
integer :: nmd_cluster
integer :: ntmd_cluster
integer :: nmdx_cluster, ntmdx_cluster
integer, dimension(ntmdx_cluster)  :: nhk1_cluster, nhk2_cluster
real*8,  dimension(3,nmdx_cluster) :: x_cluster
logical, dimension(ntmdx_cluster)  :: lMDterminator

!-----for MD atoms
integer :: n, nmdx, nmdmax
integer, dimension(nmdx) :: is
real*8,  dimension(3*nmdx,0:1) :: x
real*8,  dimension(3,3,0:1) :: h
real*8,  dimension(3*nmdmax) :: frc
integer :: nmdmax_
real*8,  dimension(3*nmdmax_) :: pixcg
real*8,  dimension(nmdmax_)   :: fnorm
integer :: nnos
real*8, dimension(nnos)  :: xlogs, vlogs
logical :: lbath_clear
integer :: nnos_atom
real*8, dimension(nnos_atom,*) :: xlogs_atom, vlogs_atom

!-----stress tensors
real*8, dimension(3,3) :: pintlr      ! stress by  long-range potentials
real*8, dimension(3,3) :: pintsr      ! stress by short-range potentials
real*8, dimension(3,3) :: pint        ! total internal stress

integer :: ierror
integer :: ioptmze
integer :: ioptmze_cell
logical :: lqninitial

real*8  :: xmsst, vxmsst, h0(3,3)

!-----declare local variables
integer, allocatable, dimension(:) :: nhk1r, nhk2r
logical, allocatable, dimension(:) :: terminator
character(50) :: filename
integer :: iunit
integer :: n3, i, it, i3, l, ia, ib, ix
integer :: ntyper, ifmdr = 0, ioptmzer = -2, ioptmze_cellr = -2
integer :: nr
integer :: nnosr
real*8  :: dtr = 0.d0
real*8  :: fact, fact2
real*8  :: xlogsr, vlogsr
integer :: istat, status
integer :: digit
real*8,  dimension(3,3) :: hred
character :: dummy
integer :: iQMMDversion = 0
integer :: npx_, npy_, npz_
logical :: lQMMDrun_, lpureMD_
!---for group I/O
integer :: nfiles, ii, n1, n1alloc, isnd
real(8), allocatable :: abuf(:)


!------allocate memory for local variables
ia = max( 1, ntype, ntmdx_cluster )
allocate( nhk1r(ia), nhk2r(ia), terminator(ia),  &
& stat=status )
if( status /= 0 ) then
    if( myid.eq.0 ) then
       write(nfile(1),*) 'memory allocation error in remd_rdatoms'
       write(nfile(2),*) 'memory allocation error in remd_rdatoms'
    end if
    ierror = 200
    return
end if

call allocate_unit_number( iunit )

ierror = 0
startif: if( lstart ) then

!-----read configuration for QM atoms
   if( myid == 0 ) then

       filename = fname_toc
       if( myid == 0 ) then
           write(nfile(1),*) 'open file(in remd_rdatoms): ',  &
&                            filename(1:len_trim(filename))
           write(nfile(2),*) 'open file(in remd_rdatoms): ',  &
&                            filename(1:len_trim(filename))
       end if
       open( iunit, file=filename, status='old', action='read',    &
&            iostat=ierror )

       blockif_b1: if( ierror == 0 ) then

          !-----load QM/MD version
          read(iunit,'(a1,i5)',iostat=istat) dummy, iQMMDversion
          if( istat > 0 .or. dummy /= '#' ) then
              !-----keep compatibility with old version
              istat = 0
              iQMMDversion = 0
              backspace iunit
          end if
          ierror = max( ierror, abs(istat) )

          if( iQMMDversion >= 1 ) then
              !-----load No. of parallel nodes
              read(iunit,'(3i5)',iostat=istat) npx_, npy_, npz_
              ierror = max( ierror, abs(istat) )

              !-----load logical variables for QM/MD run
              read(iunit,'(2l5)',iostat=istat) lQMMDrun_, lpureMD_
              ierror = max( ierror, abs(istat) )
            else
              !-----assume total node = 1
              npx_ = 1
              npy_ = 1
              npz_ = 1
              lQMMDrun_ = .false.
              lpureMD_  = .false.
          end if

          !-----error trap
          if( lQMMDrun .and. .not.lQMMDrun_ .or.  &
&             .not.lQMMDrun .and. lQMMDrun_      ) then
              ierror = 1000
              write(nfile(1),*) 'inconsistent with lQMMDrun ',  &
&                         ' in previous run:', lQMMDrun, lQMMDrun_
              write(nfile(2),*) 'inconsistent with lQMMDrun ',  &
&                         ' in previous run:', lQMMDrun, lQMMDrun_
          end if
          if( lpureMD .and. .not.lpureMD_ .or.  &
&             .not.lpureMD .and. lpureMD_      ) then
              ierror = 1000
              write(nfile(1),*) 'inconsistent with lpureMD ',  &
&                         ' in previous run:', lpureMD, lpureMD_
              write(nfile(2),*) 'inconsistent with lpureMD ',  &
&                         ' in previous run:', lpureMD, lpureMD_
          end if

          !-----load No. of QM atom for non-classical MD
          if( .not.lpureMD ) then
              read(iunit,'(i2)',iostat=istat) ntyper
              ierror = max( ierror, abs(istat) )
              if( ntyper /= ntype ) then
                 ierror = 10
                 write(nfile(1),*) 'ntyper.ne.ntype:',ntyper,ntype
                 write(nfile(2),*) 'ntyper.ne.ntype:',ntyper,ntype
              end if

              if( ierror == 0 ) then
                  read(iunit,'(2i7)',iostat=istat)  &
&                       ( nhk1r(it), nhk2r(it), it = 1, ntype )
                  ierror = max( ierror, abs(istat) )
                  do it = 1, ntype
                     if( nhk1r(it) /= nhk1(it) .or.  &
&                        nhk2r(it) /= nhk2(it)      ) then
                         ierror = 11
                         write(nfile(1),*) 'error in nhk1, nhk2'
                         write(nfile(2),*) 'error in nhk1, nhk2'
                     end if
                  end do
              end if
          end if

          if( ierror == 0 ) then
              if( iQMMDversion >= 4 ) then
                  if( iQMMDversion >= 5 ) then
                      read(iunit,'(3i5)',iostat=istat) ifmdr, ioptmzer, ioptmze_cellr
                    else
                      read(iunit,'(3i2)',iostat=istat) ifmdr, ioptmzer, ioptmze_cellr
                  end if
              else if( iQMMDversion >= 3 ) then
                  read(iunit,'(2i2)',iostat=istat) ifmdr, ioptmzer
              else
                  read(iunit,'(i2)',iostat=istat) ifmdr
              end if
              ierror = max( ierror, abs(istat) )
              if( ifmdr /= ifmd ) then
                  write(nfile(1),*) 'warning: ifmdr.ne.ifmd',  &
&                                             ifmdr,ifmd
                  write(nfile(2),*) 'warning: ifmdr.ne.ifmd',  &
&                                             ifmdr,ifmd
              end if
              if( iQMMDversion >= 3 .and. ioptmzer /= ioptmze ) then
                 write(nfile(1),*) 'warning: ioptmzer.ne.ioptmze',  &
&                                            ioptmzer,ioptmze
                 write(nfile(2),*) 'warning: ioptmzer.ne.ioptmze',  &
&                                            ioptmzer,ioptmze
              end if
              if( iQMMDversion >= 4 .and. ioptmze_cellr /= ioptmze_cell ) then
                 write(nfile(1),*) 'warning: ioptmze_cellr.ne.ioptmze_cell',  &
&                                            ioptmze_cellr,ioptmze_cell
                 write(nfile(2),*) 'warning: ioptmze_cellr.ne.ioptmze_cell',  &
&                                            ioptmze_cellr,ioptmze_cell
              end if
          end if

          if( ierror == 0 ) then
              if( iQMMDversion >= 6 ) then
                  read(iunit,'(2i11)',iostat=istat) nstepCG, nstepMD
              else
                  read(iunit,'(2i8)',iostat=istat) nstepCG, nstepMD
              end if
              ierror = max( ierror, abs(istat) )
          end if

          if( ierror == 0 ) then
              read(iunit,'(e23.15)',iostat=istat) dtr
              ierror = max( ierror, abs(istat) )
          end if

          !-----load QM-atom coordinates for non-classical MD
          if( .not.lpureMD ) then
              if( ierror == 0 ) then
                  read(iunit,'(3e23.15)',iostat=istat)  &
&                  ((ratm(ix,i), ix = 1, 3), i = 1, nhk2(ntype))
                  ierror = max( ierror, abs(istat) )
              end if
          end if


          !-----load parameters for QM/MD run
!          if( lQMMDrun .and. ierror == 0 ) then
!
!              if( ierror == 0 ) then
!
!                  read(iunit,*,iostat=istat)  &
!&                               ( terminator(it), it = 1, ntype )
!                  ierror = max( ierror, abs(istat) )
!
!                  do it = 1, ntype
!                  if( .not.lterminator(it).and.terminator(it) .or.  &
!&                  lterminator(it).and. .not.terminator(it) ) then
!                     ierror = 1000
!                     write(nfile(1),*) 'error in lterminator'
!                     write(nfile(2),*) 'error in lterminator'
!                  end if
!                  end do
!
!              end if
!
!
!           !-----load configuration for MD-cluster atoms
!
!              if( ierror == 0 ) then
!
!                  read(iunit,'(i2)',iostat=istat) ntyper
!                  ierror = max( ierror, abs(istat) )
!
!                  if( ntyper /= ntmd_cluster ) then
!                      ierror = 1010
!                      write(nfile(1),*) 'ntyper.ne.ntmd_cluster',  &
!&                                        ntyper, ntmd_cluster
!                      write(nfile(2),*) 'ntyper.ne.ntmd_cluster',  &
!&                                        ntyper, ntmd_cluster
!                  end if
!
!              end if
!
!
!              if( ierror == 0 ) then
!
!                  read(iunit,'(2i7)',iostat=istat)  &
!&                   ( nhk1r(it), nhk2r(it), it = 1, ntmd_cluster )
!                  ierror = max( ierror, abs(istat) )
!
!                  do it = 1, ntmd_cluster
!                     if( nhk1r(it) /= nhk1_cluster(it) .or.  &
!&                        nhk2r(it) /= nhk2_cluster(it)      ) then
!                         ierror = 1011
!                         write(nfile(1),*) 'error in nhk1_cluster'
!                         write(nfile(2),*) 'error in nhk1_cluster'
!                     end if
!                  end do
!
!              end if
!
!
!              if( ierror == 0 ) then
!                  read(iunit,'(3e23.15)',iostat=istat)  &
!&                    ((x_cluster(ix,i), ix = 1, 3),  &
!&                               i = 1, nhk2_cluster(ntmd_cluster))
!                  ierror = max( ierror, abs(istat) )
!              end if
!
!
!              if( ierror == 0 ) then
!
!                  read(iunit,*,iostat=istat)  &
!&                      ( terminator(it), it = 1, ntmd_cluster )
!                  ierror = max( ierror, abs(istat) )
!
!                do it = 1, ntmd_cluster
!                if( .not.lMDterminator(it).and.terminator(it) .or.  &
!&                lMDterminator(it).and. .not.terminator(it) ) then
!                     ierror = 1020
!                     write(nfile(1),*) 'error in lMDterminator'
!                     write(nfile(2),*) 'error in lMDterminator'
!                end if
!                end do
!
!              end if
!
!          end if


          if( iQMMDversion >= 2 ) then

              !-----read internal stress tensors
              if( ierror == 0 ) then
                  read(iunit,'(3e23.15)',iostat=istat)  &
&                  ( ( pintlr(ix,i), ix = 1, 3 ), i = 1, 3 ),  &
&                  ( ( pintsr(ix,i), ix = 1, 3 ), i = 1, 3 )
                  ierror = max( ierror, abs(istat) )
              end if

          end if


          if( ierror == 0 ) then
              read(iunit,'(i5)',iostat=istat) nnosr
              ierror = max( ierror, abs(istat) )
          end if


          if( ierror == 0 ) then

              if( ifmd >= 3 ) then
                  !-----check the number of thermostat
                  if( nnosr == nnos ) then
                      lbath_clear = .false.
                    else
                      write(nfile(1),*) 'warning: nnosr /= nnos',  &
&                                             nnosr, nnos
                      write(nfile(1),*) '*** info: ',  &
&                                  'Thermostats will be cleared.'
                      write(nfile(2),*) 'warning: nnosr /= nnos',  &
&                                             nnosr, nnos
                      write(nfile(2),*) '*** info: ',  &
&                                  'Thermostats will be cleared.'
                  end if
              end if

              if( .not.lbath_clear ) then

                  !-----if succeeded, read previous values
                  read(iunit,'(2e23.15)',iostat=istat)  &
&                        ( xlogs(i), vlogs(i), i = 1, nnosr )
                  ierror = max( ierror, abs(istat) )

              else if( ifmd == 0 .and. nnosr >=1 ) then

                  nnosbk = nnosr
                  !-----allocate memory
                  allocate( xlogsbk(nnosbk), vlogsbk(nnosbk),  &
&                           stat=istat )
                  ierror = max( ierror, abs(istat) )

                  if( ierror == 0 ) then
                      read(iunit,'(2e23.15)',iostat=istat)  &
&                        ( xlogsbk(i), vlogsbk(i), i = 1, nnosbk )
                      ierror = max( ierror, abs(istat) )
                  end if

              else

                  !-----if failed, read and throw away previous values
                  do i = 1, nnosr
                     if( ierror == 0 ) then
                         read(iunit,'(2e23.15)',iostat=istat)  &
&                                                xlogsr, vlogsr
                         ierror = max( ierror, abs(istat) )
                     end if
                  end do

              end if
          end if

          if( ierror == 0 ) then
              if( ifmdr == 10 ) then

                  read(iunit,'(3e23.15)',iostat=istat) xmsst, vxmsst
                  ierror = max( ierror, abs(istat) )

                  read(iunit,'(3e23.15)',iostat=istat) h0(3,3)
                  ierror = max( ierror, abs(istat) )

              end if
          end if

          if( ifmd >= 1 .and. ifmd /= ifmdr ) then
              if( ifmd == 1 ) then
                  nstepMD = 0
                else
                  nstepCG = 0
              end if
          end if

       else blockif_b1

          !----- error trap
          if( myid == 0 ) then
              write(nfile(1),*) 'cannot open file:',  &
&                               filename(1:len_trim(filename))
              write(nfile(2),*) 'cannot open file:',  &
&                               filename(1:len_trim(filename))
          end if

       end if blockif_b1
       close(iunit)

   end if
   !-----error trap
   call gimax(ierror)
   if( ierror /= 0 ) then
       if( myid.eq.0 ) then
           write(nfile(1),*) 'error : in file:',  &
&                            filename(1:len_trim(filename))
           write(nfile(2),*) 'error : in file:',  &
&                            filename(1:len_trim(filename))
       end if
       return
   end if
   if( .not.lpureMD ) call dbcast(ratm,3*natom,0)
   call ibcast(ifmdr,1,0)
   call ibcast(ioptmzer,1,0)
   call gimax(nstepCG)
   call gimax(nstepMD)
   call gdmax(dtr)
!   if( lQMMDrun ) call dbcast(x_cluster,3*nmd_cluster,0)

   call dbcast(pintlr,9,0)
   call dbcast(pintsr,9,0)
   do i = 1, 3
   do ix = 1, 3
      pint(ix,i) = pintlr(ix,i) + pintsr(ix,i)
   end do
   end do

   if( .not.lbath_clear ) ierror = 1
   call gimax(ierror)
   if( ierror == 1 ) then
       ierror = 0
       lbath_clear = .false.
       call dbcast(xlogs,nnos,0)
       call dbcast(vlogs,nnos,0)
   end if

   if( ifmdr == 10 .and. ifmd == 10 ) then
       call dbcast(xmsst,1,0)
       call dbcast(vxmsst,1,0)
       call dbcast(h0,9,0)
   end if

   if( myid.eq.0 ) then
       write(nfile(1),*) ' nstepCG, nstepMD =', nstepCG, nstepMD
       write(nfile(2),*) ' nstepCG, nstepMD =', nstepCG, nstepMD
   end if
   if( ifmd == 0 .and. nstepMD > 0 ) dt = dtr


!-----read configuration for MD atoms
   ierror = 0
   ioif: if( mod(myid,iogpsz) == 0 ) then
       nfiles=min(iogpsz,nodes-myid)

       digit = log(dble(nodes)+0.5d0)/log(10.d0) + 1.d0
       filename = fname_mts
       call get_dumpfname( filename, myid, digit )

       if( myid == 0 ) then
           write(nfile(1),*) 'open file(in remd_rdatoms): ',  &
&                            filename(1:len_trim(filename))
           write(nfile(2),*) 'open file(in remd_rdatoms): ',  &
&                            filename(1:len_trim(filename))
       end if
       open( iunit, file=filename, status='old', action='read',        &
&            form='unformatted', iostat=ierror )
       if( ierror /= 0 ) then
           !----- error trap
           write(nfile(1),*) 'myid=',myid,' cannot open file:',filename
           write(nfile(2),*) 'myid=',myid,' cannot open file:',filename
       end if
   else ioif
       isnd=(myid/iogpsz)*iogpsz
   end if ioif

   !----- error trap
   call gimax(ierror)

   blockif_b2: if( ierror == 0 ) then

   ioif2: if( mod(myid,iogpsz) == 0 ) then

       !--- read my data
       read(iunit,iostat=istat) nr
       ierror = max( ierror, abs(istat) )
!           if( n /= nr ) then
!               ierror = 100
!               write(nfile(1),*) 'myid=',myid,' n /= nr', n, nr
!               write(nfile(2),*) 'myid=',myid,' n /= nr', n, nr
!           end if
       n = nr

       n3=3*nr
       if( ierror == 0 ) then
           read(iunit,iostat=istat) (is(i),i=1,nr)
           ierror = max( ierror, abs(istat) )
       end if

       if( ierror == 0 ) then
           if( nstepMD > 0 ) then
               read(iunit,iostat=istat) ((x(i3,l),i3=1,n3),l=0,1)
             else
               !---if nstepMD == 0, do not read velocities
               read(iunit,iostat=istat) ((x(i3,l),i3=1,n3),l=0,0)
           end if
           ierror = max( ierror, abs(istat) )
       end if

       !-------Read an MD-box tensor
       if( ierror == 0 ) then
           read(iunit,iostat=istat)  &
&                      ((hred(ia,ib),  ia=1,3),ib=1,3),  &
&                     (((   h(ia,ib,l),ia=1,3),ib=1,3),l=1,1)
           ierror = max( ierror, abs(istat) )
!           if( ifmd == 4 .or. ifmd == 10 .or. ioptmze_cell > 0 ) then
!               !---for NPT or NPH ensemble, for future use.
!               do ib = 1, 3
!               do ia = 1, 3
!                  h(ia,ib,0) = hred(ia,ib)
!               end do
!               end do
!           end if
       end if

       if( ierror == 0 ) then
           read(iunit,iostat=istat) (frc(i3),i3=1,n3)
           ierror = max( ierror, abs(istat) )
       end if


       n1alloc = 0
       do ii=1, nfiles-1
          read(iunit,iostat=istat) n1
          ierror = max( ierror, abs(istat) )

          if( n1*10 > n1alloc ) then
              if( allocated(abuf) ) deallocate(abuf)
              n1alloc = n1*10
              allocate(abuf(n1alloc))
          end if

          read(iunit,iostat=istat) abuf(1:n1*10)
          ierror = max( ierror, abs(istat) )

          call cisend(myid+ii,n1,1,myid+ii,0)
          call cdsend(myid+ii,abuf,n1*10,myid+ii,0)
          call cdsend(myid+ii,hred,9,myid+ii,0)
          call cdsend(myid+ii,h(1,1,1),9,myid+ii,0)
       end do
!       if( allocated(abuf) ) deallocate(abuf)

   else ioif2

!       isnd=(myid/iogpsz)*iogpsz
       call cirecvs(myid,nr,1,isnd,0)
       n = nr
       n3=3*nr

       n1alloc = n*10
       allocate(abuf(n1alloc))

       call cdrecvs(myid,abuf,n*10,isnd,0)
       call cdrecvs(myid,hred,9,isnd,0)
       call cdrecvs(myid,h(1,1,1),9,isnd,0)

       is(1:n) = nint(abuf(1:n))
       x(1:n*3,0) = abuf(n+1:n*4)
       x(1:n*3,1) = abuf(n*4+1:n*7)
       frc(1:n*3) = abuf(n*7+1:n*10)

!       deallocate(abuf)

   end if ioif2

   if( ifmd == 4 .or. ifmd == 10 .or. ioptmze_cell > 0 ) then
       !---for NPT or NPH ensemble, for future use.
       do ib = 1, 3
          do ia = 1, 3
             h(ia,ib,0) = hred(ia,ib)
          end do
       end do
   end if

   !----- error trap
   call gimax(ierror)

   if( ierror == 0 ) then
   if( ifmdr == 1 .and. ifmd == 1 ) then

       ioif3: if( mod(myid,iogpsz) == 0 ) then

           !--- read my data
           read(iunit,iostat=istat) (pixcg(i3),i3=1,n3)
           ierror = max( ierror, abs(istat) )
           read(iunit,iostat=istat) (fnorm(i3),i3=1,nr)
           ierror = max( ierror, abs(istat) )

           do ii=1, nfiles-1
              read(iunit,iostat=istat) n1
              ierror = max( ierror, abs(istat) )

              if( n1*4 > n1alloc ) then
                  if( allocated(abuf) ) deallocate(abuf)
                  n1alloc = n1*4
                  allocate(abuf(n1alloc))
              end if

              read(iunit,iostat=istat) abuf(1:n1*4)
              ierror = max( ierror, abs(istat) )

!              call cisend(myid+ii,n1,1,myid+ii,0)
              call cdsend(myid+ii,abuf,n1*4,myid+ii,0)
           end do
!           if( allocated(abuf) ) deallocate(abuf)

       else ioif3

!           isnd=(myid/iogpsz)*iogpsz
!           call cirecvs(myid,nr,1,isnd,0)
!           n = nr
!           n3=3*nr

           if( n*4 > n1alloc ) then
               if( allocated(abuf) ) deallocate(abuf)
               n1alloc = n*4
               allocate(abuf(n1alloc))
           end if

           call cdrecvs(myid,abuf,n*4,isnd,0)

           pixcg(1:n3) = abuf(1:n3)
           fnorm(1:n)  = abuf(n3+1:n*4)

!           deallocate(abuf)

       end if ioif3


       if( ioptmzer >= 2 .and. ioptmzer <= 9 ) then
           if( ioptmze >= 2 .and. ioptmze <= 9 ) then
               call loadqnrfo( n, iunit, ierror, .true., iogpsz )
               lqninitial = .false.
           else
               call loadqnrfo( n, iunit, ierror, .false., iogpsz )
           end if
       end if

       if( ioptmze_cellr >= 0 ) then
           call loadqnrfo_cell( iunit, ierror, iogpsz )
       end if

       !----- error trap
       call gimax(ierror)

   end if
   end if

!   if( ierror == 0 ) then
!   if( ifmd == 5 ) then
!   if( ifmdr == 5 ) then
!
!       ioif4: if( mod(myid,iogpsz) == 0 ) then
!
!           !--- read my data
!           nnosr = 0
!           read(iunit,iostat=istat) nnosr
!           ierror = max( ierror, abs(istat) )
!           lbath_clear = nnosr /= nnos_atom
!
!           if( .not.lbath_clear ) then
!               read(iunit,iostat=istat) xlogs_atom(1:nnos_atom,1:n)
!               ierror = max( ierror, abs(istat) )
!
!               read(iunit,iostat=istat) vlogs_atom(1:nnos_atom,1:n)
!               ierror = max( ierror, abs(istat) )
!           else
!               read(iunit,iostat=istat) xlogsr
!               ierror = max( ierror, abs(istat) )
!
!               read(iunit,iostat=istat) vlogsr
!               ierror = max( ierror, abs(istat) )
!           end if
!
!           do ii=1, nfiles-1
!              call cisend(myid+ii,nnosr,1,myid+ii,0)
!
!              read(iunit,iostat=istat) n1
!              ierror = max( ierror, abs(istat) )
!
!              if( .not.lbath_clear ) then
!                  if( n1*nnos_atom*2 > n1alloc ) then
!                      if( allocated(abuf) ) deallocate(abuf)
!                      n1alloc = n1*nnos_atom*2
!                      allocate(abuf(n1alloc))
!                  end if
!
!                  read(iunit,iostat=istat) abuf(1:n1*nnos_atom*2)
!                  ierror = max( ierror, abs(istat) )
!
!!                  call cisend(myid+ii,n1,1,myid+ii,0)
!                  call cdsend(myid+ii,abuf,n1*nnos_atom*2,myid+ii,0)
!              else
!                  read(iunit,iostat=istat) xlogsr
!                  ierror = max( ierror, abs(istat) )
!              end if
!
!           end do
!!           if( allocated(abuf) ) deallocate(abuf)
!
!       else ioif4
!
!           call cirecvs(myid,nnosr,1,isnd,0)
!           lbath_clear = nnosr /= nnos_atom
!!           isnd=(myid/iogpsz)*iogpsz
!!           call cirecvs(myid,nr,1,isnd,0)
!!           n = nr
!!           n3=3*nr
!
!           if( .not.lbath_clear ) then
!               if( n*nnos_atom*2 > n1alloc ) then
!                   if( allocated(abuf) ) deallocate(abuf)
!                   n1alloc = n*nnos_atom*2
!                   allocate(abuf(n1alloc))
!               end if
!
!               call cdrecvs(myid,abuf,n*nnos_atom*2,isnd,0)
!
!               !-----copy abuf to xlogs_atom
!               call dcopy_a_to_b( abuf(            1), xlogs_atom, n*nnos_atom )
!               call dcopy_a_to_b( abuf(n*nnos_atom+1), vlogs_atom, n*nnos_atom )
!
!!               deallocate(abuf)
!           end if
!
!       end if ioif4
!
!       !----- error trap
!       call gimax(ierror)
!
!   else
!       lbath_clear = .true.
!   end if
!   end if
!   end if


!   if( ierror == 0 ) then
!   if( ifmd == 6 ) then
!   if( ifmdr == 6 ) then
!       call gle_read( nfile, myid, nodes, iunit, n, lbath_clear, ierror, iogpsz )
!
!       !----- error trap
!       call gimax(ierror)
!   else
!       lbath_clear = .true.
!   end if
!   end if
!   end if


!   if( ierror == 0 ) then
       !--- read variables for virtual MD for thermodynamic integration
!       call read_vmd( nfile, myid, nodes, iunit, n3, ierror )

       !----- error trap
!       call gimax(ierror)
!   end if

   if( mod(myid,iogpsz) == 0 ) close(iunit)

   if( allocated(abuf) ) deallocate(abuf)

   !---for all nodes
   !-----correct scaled velocities & forces, if dtr /= dt
   if( ifmd >= 2 .and. nstepMD > 0 .and.  &
&      abs(dt-dtr) > 1.d-13 ) then
       fact  = dt/dtr
       fact2 = fact*fact
       do i3 = 1, n3
          x(i3,1) = x(i3,1) * fact
          frc(i3) = frc(i3) * fact2
       end do
       if( .not.lbath_clear ) then
           do i = 1, nnos
              vlogs(i) = vlogs(i) * fact
           end do

           do ib = 1, 3
           do ia = 1, 3
              h(ia,ib,1) = h(ia,ib,1) * fact
           end do
           end do

!           if( ifmd == 5 .and. ifmdr == 5 ) then
!               vlogs_atom(1:nnos_atom,1:n) = vlogs_atom(1:nnos_atom,1:n) * fact
!           end if
       end if
   end if

!   else blockif_b2
!
!      !----- error trap
!      write(nfile(1),*) 'myid=',myid,' cannot open file:',filename
!      write(nfile(2),*) 'myid=',myid,' cannot open file:',filename

   end if blockif_b2

end if startif

!-----error trap
call gimax(ierror)
if( ierror /= 0 ) then
    if( myid == 0 ) then
        write(nfile(1),*) 'an error occurs in remd_rdatoms'
        write(nfile(2),*) 'an error occurs in remd_rdatoms'
    end if
end if


!------deallocate memory for local variables
deallocate( nhk1r, nhk2r, terminator, stat=status )
if( status /= 0 ) then
    if( myid.eq.0 ) then
     write(nfile(1),*) 'memory deallocation error in remd_rdatoms'
     write(nfile(2),*) 'memory deallocation error in remd_rdatoms'
    end if
    ierror = 201
    return
end if

call deallocate_unit_number( iunit )


return
end subroutine




subroutine remd_svatoms( nfile, myid, nodes, npx, npy, npz, iogpsz,  &
& fname_toc, fname_mts, lsave, ifmd, nstepCG, nstepMD, dt,  &
& lQMMDrun, lpureMD,  &
& natom, ntype, nhk1, nhk2, ratm, lterminator, natomx, ntypex,  &
& nmd_cluster, ntmd_cluster, nhk1_cluster, nhk2_cluster, x_cluster,  &
& lMDterminator, nmdx_cluster, ntmdx_cluster,  &
& n, is, x, h, frc, pixcg, fnorm, nmdx, nmdmax_,  &
& nnos, xlogs, vlogs, pintlr, pintsr, ierror,  &
& ioptmze, ioptmze_cell, xmsst, vxmsst, h0, nnos_atom, xlogs_atom, vlogs_atom )
!-----------------------------------------------------------------------
!     save atomic configuration to files, if lsave == .true.
!-----------------------------------------------------------------------
use heatbath_backup
implicit none
integer :: nfile(*), myid, nodes
integer :: npx, npy, npz
integer :: iogpsz
character(*):: fname_toc
character(*):: fname_mts
logical :: lsave
integer :: ifmd
integer :: nstepCG, nstepMD
real*8  :: dt
logical :: lQMMDrun
logical :: lpureMD

!-----for QM atoms
integer :: natom
integer :: ntype
integer :: natomx
integer :: ntypex
integer, dimension(ntypex) :: nhk1, nhk2
real*8,  dimension(3,natomx) :: ratm
logical, dimension(ntypex) :: lterminator

!-----for MD-cluster atoms
integer :: nmd_cluster
integer :: ntmd_cluster
integer :: nmdx_cluster, ntmdx_cluster
integer, dimension(ntmdx_cluster)  :: nhk1_cluster, nhk2_cluster
real*8,  dimension(3,nmdx_cluster) :: x_cluster
logical, dimension(ntmdx_cluster)  :: lMDterminator

!-----for MD atoms
integer :: n, nmdx
integer, dimension(nmdx) :: is
real*8,  dimension(3*nmdx,0:1) :: x
real*8,  dimension(3,3,0:1) :: h
real*8,  dimension(3*n) :: frc
integer :: nmdmax_
real*8,  dimension(3*nmdmax_) :: pixcg
real*8,  dimension(nmdmax_)   :: fnorm
integer :: nnos
real*8, dimension(nnos)  :: xlogs, vlogs
integer :: nnos_atom
real*8, dimension(nnos_atom,*) :: xlogs_atom, vlogs_atom

!-----stress tensors
real*8, dimension(3,3) :: pintlr      ! stress by  long-range potentials
real*8, dimension(3,3) :: pintsr      ! stress by short-range potentials

integer :: ierror
integer :: ioptmze
integer :: ioptmze_cell

real*8  :: xmsst, vxmsst, h0(3,3)

!-----declare local variables
character(50) :: filename
integer :: iunit
integer :: n3, i, it, i3, l, ia, ib, ix
integer :: istat
integer :: digit
integer :: iQMMDversion = 6
!---for group I/O
integer :: nfiles, ii, n1, n1alloc, isnd
real(8), allocatable :: abuf(:)


ierror = 0
startif: if( lsave ) then

    call allocate_unit_number( iunit )

   !-----save configuration for QM & MD-cluster atoms
   if( myid == 0 ) then

       filename = fname_toc
       if( myid == 0 ) then
           write(nfile(1),*) 'open file(in remd_svatoms): ',  &
&                            filename(1:len_trim(filename))
           write(nfile(2),*) 'open file(in remd_svatoms): ',  &
&                            filename(1:len_trim(filename))
       end if
       open( iunit, file=filename, status='unknown',  &
&            action='write', iostat=ierror )

       blockif_b1: if( ierror == 0 ) then

          !-----save QM/MD version
          write(iunit,'(a1,i5)',iostat=istat) '#', iQMMDversion
          ierror = max( ierror, abs(istat) )

          !-----save No. of parallel nodes
          write(iunit,'(3i5)',iostat=istat) npx, npy, npz
          ierror = max( ierror, abs(istat) )

          !-----save logical variables for QM/MD run
          write(iunit,'(2l5)',iostat=istat) lQMMDrun, lpureMD
          ierror = max( ierror, abs(istat) )

          !-----save No. of QM atoms for non-classical MD
          if( .not.lpureMD ) then
              write(iunit,'(i2)',iostat=istat) ntype
              ierror = max( ierror, abs(istat) )

              write(iunit,'(2i7)',iostat=istat)  &
&                       ( nhk1(it), nhk2(it), it = 1, ntype )
              ierror = max( ierror, abs(istat) )
          end if

          write(iunit,'(3i5)',iostat=istat) ifmd, ioptmze, ioptmze_cell
          ierror = max( ierror, abs(istat) )

          write(iunit,'(2i11)',iostat=istat) nstepCG, nstepMD
          ierror = max( ierror, abs(istat) )

          write(iunit,'(e23.15)',iostat=istat) dt
          ierror = max( ierror, abs(istat) )

          !-----save configuration for QM atoms for non-classical MD
          if( .not.lpureMD ) then
              write(iunit,'(3e23.15)',iostat=istat)  &
&              ((ratm(ix,i), ix = 1, 3), i = 1, nhk2(ntype))
              ierror = max( ierror, abs(istat) )
          end if


          !-----save parameters for QM/MD run
!          if( lQMMDrun ) then
!
!              write(iunit,*,iostat=istat)  &
!&                               ( lterminator(it), it = 1, ntype )
!              ierror = max( ierror, abs(istat) )
!
!           !-----save configuration for MD-cluster atoms
!
!              write(iunit,'(i2)',iostat=istat) ntmd_cluster
!              ierror = max( ierror, abs(istat) )
!
!              write(iunit,'(2i7)',iostat=istat)  &
!&                     ( nhk1_cluster(it), nhk2_cluster(it),  &
!&                                           it = 1, ntmd_cluster )
!              ierror = max( ierror, abs(istat) )
!
!              write(iunit,'(3e23.15)',iostat=istat)  &
!&                    ((x_cluster(ix,i), ix = 1, 3),  &
!&                               i = 1, nhk2_cluster(ntmd_cluster))
!              ierror = max( ierror, abs(istat) )
!
!              write(iunit,*,iostat=istat)  &
!&                      ( lMDterminator(it), it = 1, ntmd_cluster )
!              ierror = max( ierror, abs(istat) )
!
!          end if


          !-----save internal stress tensors
          write(iunit,'(3e23.15)',iostat=istat)  &
&              ( ( pintlr(ix,i), ix = 1, 3 ), i = 1, 3 ),  &
&              ( ( pintsr(ix,i), ix = 1, 3 ), i = 1, 3 )
          ierror = max( ierror, abs(istat) )


          if( ifmd >= 3 ) then

              write(iunit,'(i5)',iostat=istat) nnos
              ierror = max( ierror, abs(istat) )

              write(iunit,'(2e23.15)',iostat=istat)  &
&                        ( xlogs(i), vlogs(i), i = 1, nnos )
              ierror = max( ierror, abs(istat) )

          else if( ifmd == 0 .and. nnosbk >= 1 ) then

              write(iunit,'(i5)',iostat=istat) nnosbk
              ierror = max( ierror, abs(istat) )

              write(iunit,'(2e23.15)',iostat=istat)  &
&                        ( xlogsbk(i), vlogsbk(i), i = 1, nnosbk )
              ierror = max( ierror, abs(istat) )

          else

              write(iunit,'(i5)',iostat=istat) nnosbk
              ierror = max( ierror, abs(istat) )

          end if


          if( ifmd == 10 ) then

              write(iunit,'(3e23.15)',iostat=istat) xmsst, vxmsst
              ierror = max( ierror, abs(istat) )

              write(iunit,'(3e23.15)',iostat=istat) h0(3,3)
              ierror = max( ierror, abs(istat) )

          end if

       else blockif_b1

          !----- error trap
          if( myid == 0 ) then
              write(nfile(1),*) 'cannot open file:', filename
              write(nfile(2),*) 'cannot open file:', filename
          end if

       end if blockif_b1
       close(iunit)

   end if


!-----write configuration for MD atoms
   ierror = 0
   ioif: if( mod(myid,iogpsz) == 0 ) then
       nfiles=min(iogpsz,nodes-myid)

       digit = log(dble(nodes)+0.5d0)/log(10.d0) + 1.d0
       filename = fname_mts
       call get_dumpfname( filename, myid, digit )

       if( myid == 0 ) then
           write(nfile(1),*) 'open file(in remd_svatoms): ',  &
&                            filename(1:len_trim(filename))
           write(nfile(2),*) 'open file(in remd_svatoms): ',  &
&                            filename(1:len_trim(filename))
       end if
       open( iunit, file=filename, status='unknown',  &
&            action='write', form='unformatted', iostat=ierror )
   else ioif
       isnd=(myid/iogpsz)*iogpsz
   end if ioif

   !----- error trap
   call gimax(ierror)

   blockif_b2: if( ierror == 0 ) then

   ioif2: if( mod(myid,iogpsz) == 0 ) then

       !--- write my data
       write(iunit,iostat=istat) n
       ierror = max( ierror, abs(istat) )

       n3=3*n
       write(iunit,iostat=istat) (is(i),i=1,n)
       ierror = max( ierror, abs(istat) )

       write(iunit,iostat=istat) ((x(i3,l),i3=1,n3),l=0,1)
       ierror = max( ierror, abs(istat) )

       write(iunit,iostat=istat) (((h(ia,ib,l),ia=1,3),ib=1,3),l=0,1)
       ierror = max( ierror, abs(istat) )

       write(iunit,iostat=istat) (frc(i3),i3=1,n3)
       ierror = max( ierror, abs(istat) )


       n1alloc = 0
       do ii=1, nfiles-1
          call cirecvs(myid+ii,n1,1,myid+ii,0)

          if( n1*10 > n1alloc ) then
              if( allocated(abuf) ) deallocate(abuf)
              n1alloc = n1*10
              allocate(abuf(n1alloc))
          end if

          call cdrecvs(myid+ii,abuf,n1*10,myid+ii,0)

          write(iunit,iostat=istat) n1
          ierror = max( ierror, abs(istat) )

          write(iunit,iostat=istat) abuf(1:n1*10)
          ierror = max( ierror, abs(istat) )
       end do
!       if( allocated(abuf) ) deallocate(abuf)

   else ioif2

!       isnd=(myid/iogpsz)*iogpsz
       call cisend(myid,n,1,isnd,0)
       n3=3*n

       n1alloc = n*10
       allocate(abuf(n1alloc))

       abuf(1:n)=is(1:n)
       abuf(n+1:n*4)    = x(1:n*3,0)
       abuf(n*4+1:n*7)  = x(1:n*3,1)
       abuf(n*7+1:n*10) = frc(1:n*3)

       call cdsend(myid,abuf,n*10,isnd,0)
!       deallocate(abuf)

   end if ioif2

   !----- error trap
   call gimax(ierror)

   if( ifmd == 1 ) then

       ioif3: if( mod(myid,iogpsz) == 0 ) then

           !--- write my data
           write(iunit,iostat=istat) (pixcg(i3),i3=1,n3)
           ierror = max( ierror, abs(istat) )

           write(iunit,iostat=istat) (fnorm(i3),i3=1,n)
           ierror = max( ierror, abs(istat) )

           do ii=1, nfiles-1
              call cirecvs(myid+ii,n1,1,myid+ii,0)

              if( n1*4 > n1alloc ) then
                  if( allocated(abuf) ) deallocate(abuf)
                  n1alloc = n1*4
                  allocate(abuf(n1alloc))
              end if

              call cdrecvs(myid+ii,abuf,n1*4,myid+ii,0)

              write(iunit,iostat=istat) n1
              ierror = max( ierror, abs(istat) )

              write(iunit,iostat=istat) abuf(1:n1*4)
              ierror = max( ierror, abs(istat) )
           end do
!           if( allocated(abuf) ) deallocate(abuf)

       else ioif3

!           isnd=(myid/iogpsz)*iogpsz
           call cisend(myid,n,1,isnd,0)

           if( n*4 > n1alloc ) then
               if( allocated(abuf) ) deallocate(abuf)
               n1alloc = n*4
               allocate(abuf(n1alloc))
           end if

           abuf(1:n3)     = pixcg(1:n3)
           abuf(n3+1:n*4) = fnorm(1:n)

           call cdsend(myid,abuf,n*4,isnd,0)
!           deallocate(abuf)

       end if ioif3


       if( ioptmze >= 2 .and. ioptmze <= 9 ) then
           call outqnrfo( n, iunit, ierror, iogpsz )
       end if

       if( ioptmze_cell >= 0 ) then
           call outqnrfo_cell( iunit, ierror, iogpsz )
       end if

   end if

!   if( ifmd == 5 ) then
!
!       ioif4: if( mod(myid,iogpsz) == 0 ) then
!
!           !--- write my data
!           write(iunit,iostat=istat) nnos_atom
!           ierror = max( ierror, abs(istat) )
!
!           write(iunit,iostat=istat) xlogs_atom(1:nnos_atom,1:n)
!           ierror = max( ierror, abs(istat) )
!
!           write(iunit,iostat=istat) vlogs_atom(1:nnos_atom,1:n)
!           ierror = max( ierror, abs(istat) )
!
!           do ii=1, nfiles-1
!              call cirecvs(myid+ii,n1,1,myid+ii,0)
!
!              if( n1*nnos_atom*2 > n1alloc ) then
!                  if( allocated(abuf) ) deallocate(abuf)
!                  n1alloc = n1*nnos_atom*2
!                  allocate(abuf(n1alloc))
!              end if
!
!              call cdrecvs(myid+ii,abuf,n1*nnos_atom*2,myid+ii,0)
!
!              write(iunit,iostat=istat) n1
!              ierror = max( ierror, abs(istat) )
!
!              write(iunit,iostat=istat) abuf(1:n1*nnos_atom*2)
!              ierror = max( ierror, abs(istat) )
!           end do
!!           if( allocated(abuf) ) deallocate(abuf)
!
!       else ioif4
!
!!           isnd=(myid/iogpsz)*iogpsz
!           call cisend(myid,n,1,isnd,0)
!
!           if( n*nnos_atom*2 > n1alloc ) then
!               if( allocated(abuf) ) deallocate(abuf)
!               n1alloc = n*nnos_atom*2
!               allocate(abuf(n1alloc))
!           end if
!
!           !-----copy xlogs_atom to abuf
!           call dcopy_a_to_b( xlogs_atom, abuf(            1), n*nnos_atom )
!           call dcopy_a_to_b( vlogs_atom, abuf(n*nnos_atom+1), n*nnos_atom )
!
!           call cdsend(myid,abuf,n*nnos_atom*2,isnd,0)
!!           deallocate(abuf)
!
!       end if ioif4
!
!   end if

!   if( ifmd == 6 ) then
!       call gle_save( nfile, myid, nodes, iunit, n, ierror, iogpsz )
!   end if

   !--- save variables for virtual MD for thermodynamic integration
!   call save_vmd( nfile, myid, nodes, iunit, n3, ierror )

!   else blockif_b2
!
!      !----- error trap
!      write(nfile(1),*) 'myid=',myid,' cannot open file:',filename
!      write(nfile(2),*) 'myid=',myid,' cannot open file:',filename

   end if blockif_b2
   if( mod(myid,iogpsz) == 0 ) close(iunit)

   if( allocated(abuf) ) deallocate(abuf)

   call deallocate_unit_number( iunit )

end if startif

!-----error trap
call gimax(ierror)
if( ierror /= 0 ) then
    if( myid == 0 ) then
        write(nfile(1),*) 'an error occurs in remd_svatoms'
        write(nfile(2),*) 'an error occurs in remd_svatoms'
    end if
end if


return
end subroutine




subroutine remd_setQMpositions( nfile, myid, nodes )
!-----------------------------------------------------------------------
!     set coordinates for QM atoms in MD nodes to QM nodes
!-----------------------------------------------------------------------
use remd_param
use remd_param_atom
use remd_variables
use remd_qmmd
implicit none
integer :: nfile(*), myid, nodes


!-----get QM atom positions in cluster
call getQMpositions( nfile, myid, nodes,  &
& natom, ntype, nhk1, nhk2, ratm, lterminator, nmd, is, x, h,  &
& ind_qm1_md, sxog, syog, szog )


!-----get QM terminator atom (HSQM) positions
!call getHSQMpositions( nfile, myid, nodes,  &
!& natom, ntype, nhk1, nhk2, ratm, lterminator, nmd, is, x, h,  &
!& ind_HSQM_qm, x_HSQM_qm,  &
!& ind_HSQM_md, x_HSQM_md, alpha_HSQM, nHSQM, nHSQMx,  &
!& sxog, syog, szog )


!-----get real QM atom positions
call setrealQMpositions( nfile, myid, nodes,  &
& natom, ratm, realatm, h )


!-----get real QM atom positions
call setsupercellvectors( nfile, myid, nodes,  &
& lQMMDrun, h, hcell )


return
end subroutine




subroutine getQMpositions( nfile, myid, nodes,  &
& natom, ntype, nhk1, nhk2, Xqm1, lterminator, n, is, x, h,  &
& ind_qm1_md, sxog, syog, szog )
!-----------------------------------------------------------------------
!     get positions for QM atoms
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: natom
integer :: ntype
integer, dimension(ntype) :: nhk1, nhk2
real*8,  dimension(3,natom) :: Xqm1
logical, dimension(ntype) :: lterminator
integer :: n
integer, dimension(n) :: is
real*8,  dimension(3*n) :: x
real*8,  dimension(3,3,0:1) :: h
real*8  :: sxog, syog, szog
integer, dimension(2,natom) :: ind_qm1_md

!-----declare local variables
integer :: i, j, it
real*8 :: amin, amin_return
real*8 :: d1, d2, d3, dx, dy, dz, rr2
integer :: itemp
real*8,  dimension(3) :: dbuf
integer, dimension(2) :: ibuf


qmtypedo: do it = 1, ntype
hsqmif: if( .not.lterminator(it) ) then


    !-----for non-terminator QM atoms
    qmatomdo: do i = nhk1(it), nhk2(it)

       amin = 1d10
       do j = 1, n
          d1 = Xqm1(1,i)-(x(3*j-2)+sxog)
          d2 = Xqm1(2,i)-(x(3*j-1)+syog)
          d3 = Xqm1(3,i)-(x(3*j-0)+szog)
          if( abs(d1).gt.0.5d0 ) d1 = d1 - sign(1.d0,d1)
          if( abs(d2).gt.0.5d0 ) d2 = d2 - sign(1.d0,d2)
          if( abs(d3).gt.0.5d0 ) d3 = d3 - sign(1.d0,d3)
          dx = h(1,1,0)*d1+h(1,2,0)*d2+h(1,3,0)*d3
          dy = h(2,1,0)*d1+h(2,2,0)*d2+h(2,3,0)*d3
          dz = h(3,1,0)*d1+h(3,2,0)*d2+h(3,3,0)*d3
          rr2= dx**2+dy**2+dz**2
          if( rr2.lt.amin ) then
              itemp = j
              amin  = rr2
          end if
       end do
       amin_return = amin

       call gdmin( amin_return )

       if( abs(amin-amin_return).lt.1d-14 ) then
           Xqm1(1,i)=x(3*itemp-2)+sxog
           Xqm1(2,i)=x(3*itemp-1)+syog
           Xqm1(3,i)=x(3*itemp-0)+szog
           ind_qm1_md(1,i)=myid
           ind_qm1_md(2,i)=itemp
         else
           Xqm1(1,i)=-1d0
           Xqm1(2,i)=-1d0
           Xqm1(3,i)=-1d0
           ind_qm1_md(1,i)=-1
           ind_qm1_md(2,i)=-1
       end if
       call gdmaxpl( Xqm1(1,i), 3, dbuf )
       call gimaxpl( ind_qm1_md(1,i), 2, ibuf )

    end do qmatomdo


end if hsqmif
end do qmtypedo


return
end subroutine




subroutine setrealQMpositions( nfile, myid, nodes,  &
& natom, Xqm1, realatm, h )
!-----------------------------------------------------------------------
!     set real positions for QM atoms
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: natom
real*8,  dimension(3,natom) :: Xqm1
real*8,  dimension(3,natom) :: realatm
real*8,  dimension(3,3) :: h

!-----declare local variables
integer :: i
real*8 :: q1, q2, q3


!-----get real coordinates which will be transferred to QM nodes
do i = 1, natom
   q1 = Xqm1(1,i)
   q2 = Xqm1(2,i)
   q3 = Xqm1(3,i)
   realatm(1,i) = h(1,1)*q1+h(1,2)*q2+h(1,3)*q3
   realatm(2,i) = h(2,1)*q1+h(2,2)*q2+h(2,3)*q3
   realatm(3,i) = h(3,1)*q1+h(3,2)*q2+h(3,3)*q3
end do


return
end subroutine




subroutine setsupercellvectors( nfile, myid, nodes,  &
& lQMMDrun, h, hcell )
!-----------------------------------------------------------------------
!     set real positions for QM atoms
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
logical :: lQMMDrun
real*8, dimension(3,3,0:1) :: h
real*8, dimension(3,3)     :: hcell

!-----declare local variables
integer :: i, ix


if( .not.lQMMDrun ) then
    do i = 1, 3
    do ix = 1, 3
       hcell(ix,i) = h(ix,i,0)
    end do
    end do
end if


return
end subroutine




subroutine copyQMpositions( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!     just copy data, if root nodes are the same in both QM & MD regions
!-----------------------------------------------------------------------
use remd_param
use remd_param_atom
use remd_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: ierror


call copyMD2QM( nfile, myid, nodes,  &
& realatm, 3*natom, hcell, ierror )


return
end subroutine




subroutine sendQMpositions( nfile, myid, nodes, id_qm1, ierror )
!-----------------------------------------------------------------------
!     send QM positions to QM-root node
!-----------------------------------------------------------------------
use remd_param
use remd_param_atom
use remd_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: id_qm1
integer :: ierror

!-----declare local variables
integer :: nsd


nsd = 3*natom
call cisend(100,nsd,1,id_qm1,0)
call cdsend(110,realatm,nsd,id_qm1,0)
call cdsend(120,hcell,9,id_qm1,0)


return
end subroutine




subroutine copyQM2MD( nfile, myid, nodes,  &
& frcQM_, natom_, epotQM_, strQM_, ierror )
!-----------------------------------------------------------------------
!     just copy data
!-----------------------------------------------------------------------
use remd_param
use remd_param_atom
use remd_variables
use remd_atom
implicit none
integer :: nfile(*), myid, nodes
integer :: natom_
real*8, dimension(3*natom_) :: frcQM_
real*8  :: epotQM_
real*8, dimension(3,3) :: strQM_
integer :: ierror

!-----declare local variables
integer :: i, j


!-----error trap
if( natom /= natom_ ) ierror = 1

if( ierror /= 0 ) then
    if( myid == 0 ) then
        write(nfile(1),*) 'error in copyQM2MD : natom (or nmd) /= natom_'
        write(nfile(2),*) 'error in copyQM2MD : natom (or nmd) /= natom_'
    end if
    return
end if


!-----[Ryd.] unit -> [hartree] unit
epotQM = epotQM_ * 0.5d0
do i = 1, natom*3
   frcQM(i) = frcQM_(i) * 0.5d0
end do
do i = 1, 3
do j = 1, 3
   strQM(j,i) = strQM_(j,i) * 0.5d0
end do
end do


return
end subroutine




subroutine recvQMforces( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!     recieve data from QM nodes
!-----------------------------------------------------------------------
use remd_param
use remd_param_atom
use remd_variables
use remd_atom
implicit none
integer :: nfile(*), myid, nodes
integer :: ierror

!-----declare local variables
integer :: nrc, i, j


call cirecv(100,nrc,1,0)


!-----error trap
if( nrc /= 3*natom ) then
    if( myid == 0 ) then
        write(nfile(1),*) 'error in recvQMforces :',  &
&                         ' nrc /= 3*natom'
        write(nfile(2),*) 'error in recvQMforces :',  &
&                         ' nrc /= 3*natom'
    end if
    ierror = 2
    return
end if

call cdrecv(110,frcQM,nrc,0)
call cdrecv(120,epotQM,1,0)
call cdrecv(130,strQM,9,0)


!-----[Ryd.] unit -> [hartree] unit
epotQM = epotQM * 0.5d0
do i = 1, natom*3
   frcQM(i) = frcQM(i) * 0.5d0
end do
do i = 1, 3
do j = 1, 3
   strQM(j,i) = strQM(j,i) * 0.5d0
end do
end do


return
end subroutine




subroutine unifQMforces( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!     unify transferred forces in MD nodes
!-----------------------------------------------------------------------
use remd_param
use remd_param_atom
use remd_variables
use remd_atom
implicit none
integer :: nfile(*), myid, nodes
integer :: ierror

!-----declare local variables
integer :: i, it, i3


call dbcast(frcQM,3*natom,0)
call dbcast(strQM,9,0)


return
end subroutine




subroutine correct_force( nfile, myid, nodes,  &
& natom, ntype, nhk1, nhk2, frcQM, lterminator, n, a,  &
& ind_qm1_md, ind_HSQM_qm, ind_HSQM_md, alpha_HSQM, nHSQM, nHSQMx )
!-----------------------------------------------------------------------
!     correct force
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: natom
integer :: ntype
integer, dimension(ntype) :: nhk1, nhk2
real*8,  dimension(3*natom) :: frcQM
logical, dimension(ntype) :: lterminator
integer :: n
real*8,  dimension(3*n) :: a
integer, dimension(2,natom) :: ind_qm1_md
integer :: nHSQM, nHSQMx
integer, dimension(2,nHSQMx) :: ind_HSQM_qm
integer, dimension(2,nHSQMx) :: ind_HSQM_md
real*8,  dimension(nHSQMx) :: alpha_HSQM

!-----declare local variables
integer :: i, it, j1, j2, na


na = 0
qmtypedo: do it = 1, ntype
hsqmif: if( .not.lterminator(it) ) then


    !-----for non-terminator QM atoms
    qmatomdo: do i = nhk1(it), nhk2(it)

      !----    ind_qm1_md(1,*) -> node-id
      !----    ind_qm1_md(2,*) -> particle-id
      j1 = ind_qm1_md(1,i)
      j2 = ind_qm1_md(2,i)
      if( j1 == myid )then
          a(3*j2-2) = a(3*j2-2) + frcQM(3*i-2)
          a(3*j2-1) = a(3*j2-1) + frcQM(3*i-1)
          a(3*j2-0) = a(3*j2-0) + frcQM(3*i-0)
      end if

    end do qmatomdo


  else hsqmif


    !-----for terminator QM atoms
    qmHSatomdo: do i = nhk1(it), nhk2(it)

       na = na + 1
      !---- atoms in cluster
      !----    ind_HSQM_qm(1,*) -> node-id
      !----    ind_HSQM_qm(2,*) -> particle-id
      j1 = ind_HSQM_qm(1,na)
      j2 = ind_HSQM_qm(2,na)
      if( j1 == myid )then
        a(3*j2-2) = a(3*j2-2) + (1.d0-alpha_HSQM(na))*frcQM(3*i-2)
        a(3*j2-1) = a(3*j2-1) + (1.d0-alpha_HSQM(na))*frcQM(3*i-1)
        a(3*j2-0) = a(3*j2-0) + (1.d0-alpha_HSQM(na))*frcQM(3*i-0)
      end if

      !---- atoms not in cluster
      !----    ind_HSQM_md(1,*) -> node-id
      !----    ind_HSQM_md(2,*) -> particle-id
      j1 = ind_HSQM_md(1,na)
      j2 = ind_HSQM_md(2,na)
      if( j1 == myid )then
        a(3*j2-2) = a(3*j2-2) + alpha_HSQM(na)*frcQM(3*i-2)
        a(3*j2-1) = a(3*j2-1) + alpha_HSQM(na)*frcQM(3*i-1)
        a(3*j2-0) = a(3*j2-0) + alpha_HSQM(na)*frcQM(3*i-0)
      end if

    end do qmHSatomdo


end if hsqmif
end do qmtypedo


return
end subroutine




subroutine correct_stress( nfile, myid, nodes,  &
& pintlr, pintsr, pint, pit1, strQM )
!----------------------------------------------------------------------c
!     correct stress
!----------------------------------------------------------------------c
implicit none
integer :: nfile(*), myid, nodes
real*8, dimension(3,3) :: pintlr
real*8, dimension(3,3) :: pintsr
real*8, dimension(3,3) :: pint
real*8, dimension(3,3) :: pit1
real*8, dimension(3,3) :: strQM

!-----declare local variables
integer :: i, ix


do i = 1, 3
do ix = 1, 3
   pintlr(ix,i) = pintlr(ix,i) + strQM(ix,i)
     pint(ix,i) = pintlr(ix,i) + pintsr(ix,i)
     pit1(ix,i) = pint(ix,i)
end do
end do


return
end subroutine




subroutine remd_setQMvelocities( nfile, myid, nodes )
!-----------------------------------------------------------------------
!     set velocities of QM atoms in MD nodes to QM nodes
!-----------------------------------------------------------------------
use remd_param
use remd_param_atom
use remd_variables
use remd_qmmd
implicit none
integer :: nfile(*), myid, nodes


!-----get QM atom velocities
call getQMvelocities( nfile, myid, nodes,  &
& natom, ntype, nhk1, nhk2, vatm, lterminator, nmd, is, x(1,1), h,  &
& ind_qm1_md, dtmd )


return
end subroutine




subroutine getQMvelocities( nfile, myid, nodes,  &
& natom, ntype, nhk1, nhk2, Vqm1, lterminator, n, is, v, h,  &
& ind_qm1_md, dt )
!-----------------------------------------------------------------------
!     get real velocities of QM atoms
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: natom
integer :: ntype
integer, dimension(ntype) :: nhk1, nhk2
real*8,  dimension(3,natom) :: Vqm1
logical, dimension(ntype) :: lterminator
integer :: n
integer, dimension(n) :: is
real*8,  dimension(3*n) :: v
real*8,  dimension(3,3,0:1) :: h
integer, dimension(2,natom) :: ind_qm1_md
real*8  :: dt

!-----declare local variables
integer :: i, it, j1, j2
real*8 :: d1, d2, d3, dx, dy, dz
integer :: itemp
real*8,  dimension(3) :: dbuf


qmtypedo: do it = 1, ntype
hsqmif: if( .not.lterminator(it) ) then


    !-----for non-terminator QM atoms
    qmatomdo: do i = nhk1(it), nhk2(it)

      !----    ind_qm1_md(1,*) -> node-id
      !----    ind_qm1_md(2,*) -> particle-id
      j1 = ind_qm1_md(1,i)
      j2 = ind_qm1_md(2,i)
      if( j1 == myid )then
          d1 = v(3*j2-2)/dt
          d2 = v(3*j2-1)/dt
          d3 = v(3*j2-0)/dt
          dx = h(1,1,0)*d1+h(1,2,0)*d2+h(1,3,0)*d3
          dy = h(2,1,0)*d1+h(2,2,0)*d2+h(2,3,0)*d3
          dz = h(3,1,0)*d1+h(3,2,0)*d2+h(3,3,0)*d3
          Vqm1(1,i) = dx
          Vqm1(2,i) = dy
          Vqm1(3,i) = dz
        else
          Vqm1(1,i)=-1d10
          Vqm1(2,i)=-1d10
          Vqm1(3,i)=-1d10
      end if
      !---call gdmaxpl( Vqm1(1,i), 3, dbuf )
      call dbcast( Vqm1(1,i), 3, j1 )

    end do qmatomdo


end if hsqmif
end do qmtypedo


return
end subroutine




subroutine copyQMvelocities( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!     just copy data, if root nodes are the same in both QM & MD regions
!-----------------------------------------------------------------------
use remd_param
use remd_param_atom
use remd_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: ierror


call copy_velocities_MD2QM( nfile, myid, nodes,  &
& vatm, 3*natom, ierror )


return
end subroutine




subroutine sendQMvelocities( nfile, myid, nodes, id_qm1, ierror )
!-----------------------------------------------------------------------
!     send QM velocities to QM-root node
!-----------------------------------------------------------------------
use remd_param
use remd_param_atom
use remd_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: id_qm1
integer :: ierror

!-----declare local variables
integer :: nsd


nsd = 3*natom
call cisend(100,nsd,1,id_qm1,0)
call cdsend(110,vatm,nsd,id_qm1,0)


return
end subroutine




module vsendback_in_MD
!-----------------------------------------------------------------------
! type declaration and initialization for copying back to MD nodes
!-----------------------------------------------------------------------
implicit none

logical :: lvsendback

save
end module




subroutine copylvsendback( nfile, myid, nodes, lvsendback_ )
!-----------------------------------------------------------------------
!     set lvsendback
!-----------------------------------------------------------------------
use vsendback_in_MD
implicit none
integer :: nfile(*), myid, nodes
logical :: lvsendback_

lvsendback = lvsendback_

return
end subroutine




subroutine copyQM2MDv( nfile, myid, nodes,  &
& vatm_, ndata, ierror )
!-----------------------------------------------------------------------
!     just copy data
!-----------------------------------------------------------------------
use remd_param
use remd_param_atom
use remd_variables
use remd_atom
implicit none
integer :: nfile(*), myid, nodes
integer :: ndata
real*8, dimension(3,ndata) :: vatm_
integer :: ierror


!-----error trap
if( ndata /= natom ) then
    if( myid == 0 ) then
        write(nfile(1),*) 'error in copyQM2MDv : ndata /= natom'
        write(nfile(2),*) 'error in copyQM2MDv : ndata /= natom'
    end if
    ierror = 1
    return
end if


!-----[Ryd.] unit -> [hartree] unit
vatm(1:3,1:natom) = vatm_(1:3,1:natom) * 0.5d0


return
end subroutine




subroutine recv_velocities_fromQM( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!     recieve data from QM nodes
!-----------------------------------------------------------------------
use remd_param
use remd_param_atom
use remd_variables
use remd_atom
use vsendback_in_MD
implicit none
integer :: nfile(*), myid, nodes
integer :: ierror

!-----declare local variables
integer :: nrc


call clrecv(90,lvsendback,1,0)

if( .not.lvsendback ) return

call cirecv(100,nrc,1,0)


!-----error trap
if( nrc /= 3*natom ) then
    if( myid == 0 ) then
        write(nfile(1),*) 'error in recvQMvelocities :',  &
&                         ' nrc /= 3*natom'
        write(nfile(2),*) 'error in recvQMvelocities :',  &
&                         ' nrc /= 3*natom'
    end if
    ierror = 2
    return
end if

call cdrecv(110,vatm,nrc,0)


!-----[Ryd.] unit -> [hartree] unit
vatm(1:3,1:natom) = vatm(1:3,1:natom) * 0.5d0


return
end subroutine




subroutine unify_velocities_inMD( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!     unify transferred velocities in MD nodes
!-----------------------------------------------------------------------
use constants
use remd_param
use remd_param_atom
use remd_variables
use remd_qmmd
use remd_atom
use vsendback_in_MD
implicit none
integer :: nfile(*), myid, nodes
integer :: ierror

!-----declare local variables
real*8  :: ekin
integer :: i


call lbcast(lvsendback,1,0)

if( .not.lvsendback ) return

call dbcast(vatm,3*natom,0)


!-----get QM atom velocities
call correct_velocities( nfile, myid, nodes,  &
& natom, ntype, nhk1, nhk2, vatm, lterminator, nmd, is, x(1,1), h,  &
& ind_qm1_md, dtmd )

!---check temperature
call getekin( nfile, myid, nodes,  &
& ekin, x(1,1), is, nmd, ntype, fack, h )
if( myid == 0 ) then
    do i = 1, 2
       write(nfile(i),'(a,f10.4)')  &
&                  ' *** Temperature after sending back to MD nodes: ',  &
&                           2.d0*ekin*tempau/dble(3*ntot(0))
    end do
end if


return
end subroutine




subroutine correct_velocities( nfile, myid, nodes,  &
& natom, ntype, nhk1, nhk2, Vqm1, lterminator, n, is, v, h,  &
& ind_qm1_md, dt )
!-----------------------------------------------------------------------
!     get real velocities of QM atoms
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: natom
integer :: ntype
integer, dimension(ntype) :: nhk1, nhk2
real*8,  dimension(3,natom) :: Vqm1
logical, dimension(ntype) :: lterminator
integer :: n
integer, dimension(n) :: is
real*8,  dimension(3*n) :: v
real*8,  dimension(3,3,0:1) :: h
integer, dimension(2,natom) :: ind_qm1_md
real*8  :: dt

!-----declare local variables
integer :: i, it, j1, j2
real*8 :: d1, d2, d3, dx, dy, dz, vtemp
real*8,  dimension(3,3) :: hi
integer :: itemp
real*8,  dimension(3) :: dbuf


!----- transpose matrix of hi = inverse of h
CALL RCIPRL( h(1,1,0), hi, vtemp )

qmtypedo: do it = 1, ntype
hsqmif: if( .not.lterminator(it) ) then


    !-----for non-terminator QM atoms
    qmatomdo: do i = nhk1(it), nhk2(it)

      !----    ind_qm1_md(1,*) -> node-id
      !----    ind_qm1_md(2,*) -> particle-id
      j1 = ind_qm1_md(1,i)
      j2 = ind_qm1_md(2,i)
      if( j1 == myid )then
          dx = Vqm1(1,i)
          dy = Vqm1(2,i)
          dz = Vqm1(3,i)
          d1 = hi(1,1)*dx+hi(2,1)*dy+hi(3,1)*dz
          d2 = hi(1,2)*dx+hi(2,2)*dy+hi(3,2)*dz
          d3 = hi(1,3)*dx+hi(2,3)*dy+hi(3,3)*dz
          v(3*j2-2) = d1*dt
          v(3*j2-1) = d2*dt
          v(3*j2-0) = d3*dt
      end if

    end do qmatomdo


end if hsqmif
end do qmtypedo


return
end subroutine




subroutine setup( nfile, myid, nodes,  &
& is, n, ntype, ntot, fack, acon, am, dt )
!----------------------------------------------------------------------c
!     Sets up basic constants.
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: n
integer, dimension(n) :: is
integer :: ntype
integer, dimension(0:ntype) :: ntot
real*8,  dimension(ntype) :: fack
real*8,  dimension(ntype) :: acon
real*8,  dimension(ntype) :: am
real*8 :: dt

!-----declare local variables
integer :: i, it, ic
integer, allocatable, dimension(:)   :: ibuf
integer :: status


!------allocate memory
allocate( ibuf(ntype), stat=status )

!-----set the number of atoms
do it = 1, ntype
   ntot(it) = 0
   do i = 1, n
      if( is(i) == it ) ntot(it) = ntot(it) + 1
   end do
end do
call gisum(ntot(1),ntype,ibuf)
ntot(0) = 0
do it = 1, ntype
   ntot(0) = ntot(0) + ntot(it)
end do


!-----Prefactors for kinetic energy, FACK
do ic = 1, ntype
   fack(ic) = 0.5d0*am(ic)/dt**2
end do

!-----Prefactors for normalized acceleration, ACON
do ic = 1, ntype
   acon(ic) = 0.5d0*dt*dt/am(ic)
end do


!------deallocate memory
deallocate( ibuf, stat=status )


return
end subroutine




subroutine ivelct( nfile, myid, nodes,  &
& x, v, is, n, ntype, ntot, am, fack, h, sxog, syog, szog,  &
& treq, dt, ifmd, nstepMD, tempau, vrandom )
!-----------------------------------------------------------------------
!     initial velocities
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: n
real*8,  dimension(3*n) :: x
real*8,  dimension(3*n) :: v
integer, dimension(n) :: is
integer :: ntype
real*8, dimension(ntype) :: am
integer, dimension(0:ntype) :: ntot
real*8,  dimension(ntype) :: fack
real*8,  dimension(3,3,0:1) :: h
real*8 :: sxog, syog, szog
real*8  :: treq, dt
integer :: ifmd, nstepMD
real*8  :: tempau
logical, dimension(ntype) :: vrandom

!-----declare local variables
integer :: i, ic, i3, ia
real*8  :: seed, rnd1, rnd2, twopi, vtemp
real*8, dimension(3,3) :: hi
real*8  :: vx, vy, vz
real*8  :: ekin
integer :: ix, iy, iz
real*8, dimension(9) :: dbuf
real*8, allocatable, dimension(:)   :: facv
integer :: status
logical :: lclemom = .false.


if( nstepMD > 0 ) return

if( ifmd <= 1 .or. treq.lt.1.d-05/tempau ) then
    do i = 1, 3*n
       v(i) = 0.d0
    end do
    return
end if


!------allocate memory
allocate( facv(ntype), stat=status )

!-----Prepare Maxwellian velocities
do ic = 1, ntype
   facv(ic) = dsqrt(2d0*treq/am(ic)) * dt
end do

!      seed=13597d0
seed=25d0
!-----Use different random-number sequences for different nodes
if( myid > 0 ) then
    do i=1,myid*6*n
       CALL rnd00( rnd1, seed )
    end do
end if
!-----Physical velocity
twopi = 2.d0*acos(-1.d0)
do i3 = 1, 3*n
   i=(i3+2)/3
   CALL rnd00( rnd1, seed )
   CALL rnd00( rnd2, seed )
   if( vrandom(is(i)) ) then
       v(i3)=facv(is(i))*dsqrt(-dlog(rnd1))*dcos(twopi*rnd2)
     else
       v(i3) = v(i3) * dt
   end if
end do
!-----Reduced velocity
!      call matinv(h(1,1,0),hi)
!----- transpose matrix of hi = inverse of h
CALL RCIPRL( h(1,1,0), hi, vtemp )
do i = 1, n
   iz=3*i
   iy=iz-1
   ix=iy-1
   vx=v(ix)
   vy=v(iy)
   vz=v(iz)
   v(ix)=hi(1,1)*vx+hi(2,1)*vy+hi(3,1)*vz
   v(iy)=hi(1,2)*vx+hi(2,2)*vy+hi(3,2)*vz
   v(iz)=hi(1,3)*vx+hi(2,3)*vy+hi(3,3)*vz
end do

do i = 1, ntype
   lclemom = lclemom .or. .not.vrandom(i)
end do
lclemomif: if( .not.lclemom ) then

    !-----Make total momentum zero
    call momzero( nfile, myid, nodes,  &
& x, v, is, n, ntype, ntot, am, h, sxog, syog, szog )

end if lclemomif


!--- check temperature
!--- Kinetic E.
call getekin( nfile, myid, nodes,  &
& ekin, v, is, n, ntype, fack, h )

if( myid.eq.0 ) then
    do i = 1, 2
       write(nfile(i),'(a30,f10.4)')  &
&                  '  initial temp.             : ',  &
&                           2.d0*ekin*tempau/dble(3*ntot(0))
    end do
end if


!--- correct velocities
if( .not.lclemom ) then
    call vscale( nfile, myid, nodes,  &
& treq, ekin, v, is, n, ntype, ntot )
end if


!------deallocate memory
deallocate( facv, stat=status )


return
end




subroutine momzero( nfile, myid, nodes,  &
& x, v, is, n, ntype, ntot, am, h, sxog, syog, szog )
!-----------------------------------------------------------------------
!     Make the total momentum/angular momentum zero
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: n
real*8  :: x(3*n), v(3*n)
integer :: is(n)
integer :: ntype
real*8  :: am(ntype)
integer :: ntot(0:ntype)
real*8  :: h(3,3,0:1)
real*8  :: sxog, syog, szog

!-----declare local variables
integer :: i, ic, ia, ib
real*8  :: totmassl, totmass, vtemp
real*8  :: scom(3), pav(3)
real*8  :: ainer(3,3), aineri(3,3), hi(3,3)
real*8  :: prxl, pryl, przl
real*8  :: sx, sy, sz, rx, ry, rz, r2x, r2y, r2z, r2
real*8  :: svx, svy, svz, vx, vy, vz, px, py, pz, prx, pry, prz
real*8  :: svrx, svry, svrz, vrx, vry, vrz
real*8  :: omegax, omegay, omegaz
integer :: ix, iy, iz
real*8  :: dbuf(9)
!      real*8, allocatable, dimension(:,:) :: vav
integer :: status
logical :: langmom = .false.


!-----Make total momentum zero
totmassl = 0d0
do ia = 1, 3
   scom(ia) = 0d0
   pav(ia)  = 0d0
end do
do i = 1, n
   ic = is(i)
   totmassl = totmassl + am(ic)
   scom(1) = scom(1) + ( x(3*i-2) + sxog )*am(ic)
   scom(2) = scom(2) + ( x(3*i-1) + syog )*am(ic)
   scom(3) = scom(3) + ( x(3*i-0) + szog )*am(ic)
   do ia=1,3
      pav(ia) = pav(ia) + am(ic)*v(3*i-3+ia)
   end do
end do
totmass = totmassl
call gdsum(totmass,1,vtemp)
call gdsum(scom,3,dbuf)
do ia = 1, 3
   scom(ia) = scom(ia)/totmass
end do
call gdsum(pav,3,dbuf)
do ia=1,3
   pav(ia) = pav(ia)/ntot(0)
end do
if( myid ==0 ) then
    write(nfile(1),998) scom(1),scom(2),scom(3)
    write(nfile(2),998) scom(1),scom(2),scom(3)
    write(nfile(1),999) pav(1),pav(2),pav(3)
    write(nfile(2),999) pav(1),pav(2),pav(3)
end if
998   format(/'  Center-of-mass(reduced unit) = ',3e12.4)
999   format('  Momentum/atom(reduced unit)  = ',3e12.4)

!-----Shift velocities to make the total momentum zero
!      do ic = 1, ntype
!         do ia = 1, 3
!            vav(ia,ic) = pav(ia)/am(ic)
!         end do
!      end do
do ia = 1, 3
   pav(ia)  = pav(ia)/totmass*ntot(0)
end do
do i = 1, n
   do ia = 1, 3
!            v(3*i-3+ia) = v(3*i-3+ia) - vav(ia,is(i))
      v(3*i-3+ia) = v(3*i-3+ia) - pav(ia)
   end do
end do

!-----Make the angular momentum zero
angmomif: if( langmom ) then
!-----Calculate inertia tensor & angular momentum w.r.t the COM
do ib=1,3
   do ia=1,3
      ainer(ia,ib) = 0d0
   end do
end do
prxl=0d0
pryl=0d0
przl=0d0
do i = 1, n
  ic=is(i)
  iz=3*i
  iy=iz-1
  ix=iy-1
  sx=x(ix)+sxog-scom(1)
  sy=x(iy)+syog-scom(2)
  sz=x(iz)+szog-scom(3)
  rx=h(1,1,0)*sx+h(1,2,0)*sy+h(1,3,0)*sz
  ry=h(2,1,0)*sx+h(2,2,0)*sy+h(2,3,0)*sz
  rz=h(3,1,0)*sx+h(3,2,0)*sy+h(3,3,0)*sz
  r2x=rx**2
  r2y=ry**2
  r2z=rz**2
  r2=r2x+r2y+r2z
  ainer(1,1)=ainer(1,1)+am(ic)*(r2-r2x)
  ainer(2,2)=ainer(2,2)+am(ic)*(r2-r2y)
  ainer(3,3)=ainer(3,3)+am(ic)*(r2-r2z)
  ainer(2,3)=ainer(2,3)-am(ic)*ry*rz
  ainer(3,1)=ainer(3,1)-am(ic)*rz*rx
  ainer(1,2)=ainer(1,2)-am(ic)*rx*ry
  svx=v(ix)
  svy=v(iy)
  svz=v(iz)
  vx=h(1,1,0)*svx+h(1,2,0)*svy+h(1,3,0)*svz
  vy=h(2,1,0)*svx+h(2,2,0)*svy+h(2,3,0)*svz
  vz=h(3,1,0)*svx+h(3,2,0)*svy+h(3,3,0)*svz
  px=am(ic)*vx
  py=am(ic)*vy
  pz=am(ic)*vz
  prxl=prxl+ry*pz-rz*py
  pryl=pryl+rz*px-rx*pz
  przl=przl+rx*py-ry*px
end do
ainer(3,2)=ainer(2,3)
ainer(1,3)=ainer(3,1)
ainer(2,1)=ainer(1,2)
!      do ib=1,3
!        do ia=1,3
!          dbuf(ia+3*(ib-1))=ainer(ia,ib)
!        end do
!      end do
call gdsum(ainer,9,dbuf)
prx = prxl
call gdsum(prx,1,vtemp)
pry = pryl
call gdsum(pry,1,vtemp)
prz = przl
call gdsum(prz,1,vtemp)
if( myid == 0 ) then
    do i = 1, 2
       write(nfile(i),997) (ainer(1,ib),ib=1,3),prx
       write(nfile(i),997) (ainer(2,ib),ib=1,3),pry
       write(nfile(i),997) (ainer(3,ib),ib=1,3),prz
    end do
end if
997   format(' Inertia/anglm(reduced unit)  = ',3e12.4,5x,e12.4)

!-----Calculate angular velocity
!      call matinv(ainer,aineri)
!----- transpose matrix of aineri = inverse of ainer
CALL RCIPRL( ainer, aineri, vtemp )
omegax=aineri(1,1)*prx+aineri(2,1)*pry+aineri(3,1)*prz
omegay=aineri(1,2)*prx+aineri(2,2)*pry+aineri(3,2)*prz
omegaz=aineri(1,3)*prx+aineri(2,3)*pry+aineri(3,3)*prz

!-----Shift velocities to make the angular momentum zero
!----- transpose matrix of hi = inverse of h
CALL RCIPRL( h(1,1,0), hi, vtemp )
do i = 1, n
   iz=3*i
   iy=iz-1
   ix=iy-1
   sx=x(ix)+sxog-scom(1)
   sy=x(iy)+syog-scom(2)
   sz=x(iz)+szog-scom(3)
   rx=h(1,1,0)*sx+h(1,2,0)*sy+h(1,3,0)*sz
   ry=h(2,1,0)*sx+h(2,2,0)*sy+h(2,3,0)*sz
   rz=h(3,1,0)*sx+h(3,2,0)*sy+h(3,3,0)*sz
   vrx=omegay*rz-omegaz*ry
   vry=omegaz*rx-omegax*rz
   vrz=omegax*ry-omegay*rx
   svrx=hi(1,1)*vrx+hi(2,1)*vry+hi(3,1)*vrz
   svry=hi(1,2)*vrx+hi(2,2)*vry+hi(3,2)*vrz
   svrz=hi(1,3)*vrx+hi(2,3)*vry+hi(3,3)*vrz
   v(ix)=v(ix)-svrx
   v(iy)=v(iy)-svry
   v(iz)=v(iz)-svrz
end do
end if angmomif

!check
do ia = 1, 3
   pav(ia)  = 0d0
end do
do i = 1, n
   ic = is(i)
   do ia=1,3
      pav(ia) = pav(ia) + am(ic)*v(3*i-3+ia)
   end do
end do
call gdsum(pav,3,dbuf)
do ia=1,3
   pav(ia) = pav(ia)/ntot(0)
end do

if( myid.eq.0 ) then
    do i = 1, 2
       write(nfile(i),'(a30,3e15.6)')  &
& '  corrected total momentum  : ', ( pav(ia), ia = 1, 3 )
    end do
end if


return
end




subroutine init_nhc( nfile, myid, nodes,  &
& nnos, nresn, nyosh, tomega, gnkt, gkt,  &
& xlogs, vlogs, glogs, qmass, wdti, wdti2, wdti4, wdti8,  &
& ntot, treq, dmt, lbath_clear, ierror )
!----------------------------------------------------------------------c
!  Initilize the thermostats.
!----------------------------------------------------------------------c
implicit none
integer :: nfile(*), myid, nodes
integer :: nnos, nresn, nyosh
real*8  :: tomega
real*8  :: gnkt
real*8  :: gkt
real*8, dimension(nnos)  :: xlogs, vlogs, glogs, qmass
real*8, dimension(nyosh) :: wdti, wdti2, wdti4, wdti8
integer :: ntot
real*8  :: treq
real*8  :: dmt
logical :: lbath_clear
integer :: ierror

!-----declare local variables
real*8  :: qomega              ! frequency  for thermostat in [a.u.]
real*8  :: wdti_1, dseed, vsqr, fscale
integer :: i, iyosh


!-----Set thermostat masses
qomega = 2.d0*3.141592653589793d0/tomega
gnkt = 3d0*ntot*treq
gkt  =          treq
qmass(1)=gnkt/qomega**2
do i = 2, nnos
  qmass(i)=qmass(1)/dble(3*ntot)
end do

!-----Set higher-order parameters for Yoshida-Suzuki integration
wdti_1 = dmt/dble(nresn)
if( nyosh == 1 ) then
    wdti(1) = wdti_1
else if( nyosh == 3 ) then
  wdti(1) = (1d0/(2d0-2d0**(1d0/3d0)))*wdti_1
  wdti(2) = 1d0-2d0*wdti(1)
  wdti(3) = wdti(1)
else if( nyosh == 5 ) then
  wdti(1) = (1d0/(4d0-4d0**(1d0/3d0)))*wdti_1
  wdti(2) = wdti(1)
  wdti(3) = 1d0-4d0*wdti(1)
  wdti(4) = wdti(1)
  wdti(5) = wdti(1)
else
  write(nfile(1),*) 'unsupported nyosh selected--now quitting'
  ierror = 1000
  return
end if
do iyosh = 1, nyosh
  wdti2(iyosh) = 0.5d0*wdti(iyosh)
  wdti4(iyosh) = 0.5d0*wdti2(iyosh)
  wdti8(iyosh) = 0.5d0*wdti4(iyosh)
end do

!-----Initialize heat-bath state
if( lbath_clear ) then
    do i = 1, nnos
      xlogs(i) = 0d0
      vlogs(i) = 0d0
      glogs(i) = 0d0
    end do
    dseed = 23219d0
    vsqr  = 0d0
    do i = 1, nnos
       call rnd00( vlogs(i), dseed )
       vsqr = vsqr + qmass(i)*vlogs(i)*vlogs(i)
    end do
    !fscale = dsqrt(3d0*nnos*gkt/vsqr)
    fscale = dsqrt(nnos*gkt/vsqr)
    do i = 1, nnos
      vlogs(i) = vlogs(i)*fscale
    end do
end if


return
end subroutine




subroutine set_pext( nfile, myid, nodes, pext, hpext )
!-----------------------------------------------------------------------
!     set external stress tensor
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
real*8, dimension(3,3) :: pext
real*8  :: hpext

!-----declare local variables
integer :: i, j


do i = 1, 3
do j = 1, 3
   pext(j,i) = 0.d0
end do
end do
pext(1,1) = hpext
pext(2,2) = hpext
pext(3,3) = hpext


return
end subroutine



subroutine init_nhcp( nfile, myid, nodes,  &
& nnos, nresn, nyosh, tomega, gnkt, gkt,  &
& xlogs, vlogs, glogs, qmass, wdti, wdti2, wdti4, wdti8,  &
& tbomega, blkmod, volum0, gnd2kt, onf, bmass, fack2, fack, nc, h,  &
& ntot, treq, dmt, lbath_clear, ierror )
!-----------------------------------------------------------------------
!  Initilize the thermostat & the barostat.
!-----------------------------------------------------------------------
use constants
implicit none
integer :: nfile(*), myid, nodes
!-----for thermostat
integer :: nnos, nresn, nyosh
real*8  :: tomega
real*8  :: gnkt
real*8  :: gkt
real*8, dimension(nnos)  :: xlogs, vlogs, glogs, qmass
real*8, dimension(nyosh) :: wdti, wdti2, wdti4, wdti8
!-----for barostat
real*8  :: tbomega, blkmod, volum0
real*8  :: gnd2kt, onf, bmass
integer :: nc
real*8, dimension(nc) :: fack2, fack
real*8, dimension(3,3,0:1) :: h

integer :: ntot
real*8  :: treq
real*8  :: dmt
logical :: lbath_clear
integer :: ierror

!-----declare local variables
real*8  :: qomega              ! frequency  for thermostat in [a.u.]
real*8  :: bomega              ! frequency  for thermostat in [a.u.]
real*8  :: wdti_1, dseed, vsqr, fscale
integer :: i, iyosh, ic, ia, ib


!-----Set thermostat & barostat masses
qomega = 2.d0*3.141592653589793d0/tomega
bomega = 2.d0*3.141592653589793d0/tbomega
do ic=1,nc
   fack2(ic)=2d0*fack(ic)
end do
gnkt = 3d0*ntot*treq
gkt  =          treq
gnd2kt= gnkt+9d0*gkt
onf   = 1d0/dble(3*ntot)
!-----Masses for thermostats
qmass(1)=gnkt/qomega**2
do i = 2, nnos
  qmass(i)=qmass(1)*onf
end do
!-----Mass for barostat
!      bmass=gnd2kt/bomega**2
bmass=3.d0*blkmod*volum0/bomega**2
if( myid == 0 ) then
      write(nfile(1),'(a,2es16.8)') ' '
      write(nfile(1),'(a,2es16.8)')  &
& '  bmass [a.u. / kg m^-2] =', bmass, bmass*(welm*1.d-3*(audang*1.d-10)**2)
      write(nfile(1),'(a,2es16.8)') ' '
      write(nfile(2),'(a,2es16.8)') ' '
      write(nfile(2),'(a,2es16.8)')  &
& '  bmass [a.u. / kg m^-2] =', bmass, bmass*(welm*1.d-3*(audang*1.d-10)**2)
      write(nfile(2),'(a,2es16.8)') ' '
end if

!-----Set higher-order parameters for Yoshida-Suzuki integration
wdti_1 = dmt/dble(nresn)
if( nyosh == 1 ) then
  wdti(1) = wdti_1
else if( nyosh == 3 ) then
  wdti(1) = (1d0/(2d0-2d0**(1d0/3d0)))*wdti_1
  wdti(2) = 1d0-2d0*wdti(1)
  wdti(3) = wdti(1)
else if( nyosh == 5 ) then
  wdti(1) = (1d0/(4d0-4d0**(1d0/3d0)))*wdti_1
  wdti(2) = wdti(1)
  wdti(3) = 1d0-4d0*wdti(1)
  wdti(4) = wdti(1)
  wdti(5) = wdti(1)
else
  write(nfile(1),*) 'unsupported nyosh selected--now quitting'
  ierror = 1000
  return
end if
do iyosh = 1, nyosh
  wdti2(iyosh) = 0.5d0*wdti(iyosh)
  wdti4(iyosh) = 0.5d0*wdti2(iyosh)
  wdti8(iyosh) = 0.5d0*wdti4(iyosh)
end do

!-----Initialize thermostat & barostat states
if( lbath_clear ) then
    do i = 1, nnos
       xlogs(i) = 0d0
       vlogs(i) = 0d0
       glogs(i) = 0d0
    end do
    dseed = 23219d0
    vsqr  = 0d0
    do i = 1, nnos
       call rnd00( vlogs(i), dseed )
       vsqr = vsqr + qmass(i)*vlogs(i)*vlogs(i)
    end do
    fscale = dsqrt(3d0*nnos*gkt/vsqr)
    do i = 1, nnos
       vlogs(i) = vlogs(i)*fscale
    end do

    do ib=1,3
    do ia=1,3
       h(ia,ib,1) = 0.d0
    end do
    end do
end if

return
end



subroutine init_msst( nfile, myid, nodes,  &
& nnos, nresn, nyosh, tomega, gnkt, gkt,  &
& xlogs, vlogs, glogs, qmass, wdti, wdti2, wdti4, wdti8,  &
& tbomega, blkmod, h0, volum0, alpmatrix, betmatrix,  &
& gnd2kt, onf, bmass, fack2, fack, nc, h, totmass, am,  &
& ntot, treq, dmt, lbath_clear, ierror )
!-----------------------------------------------------------------------
!  Initilize the thermostat & the barostat.
!-----------------------------------------------------------------------
use constants
implicit none
integer :: nfile(*), myid, nodes
!-----for thermostat
integer :: nnos, nresn, nyosh
real*8  :: tomega
real*8  :: gnkt
real*8  :: gkt
real*8, dimension(nnos)  :: xlogs, vlogs, glogs, qmass
real*8, dimension(nyosh) :: wdti, wdti2, wdti4, wdti8
!-----for barostat
real*8  :: tbomega, blkmod, h0(3,3), volum0, alpmatrix(3,3), betmatrix(3,3)
real*8  :: gnd2kt, onf, bmass
integer :: nc
real*8, dimension(nc) :: fack2, fack
real*8, dimension(3,3,0:1) :: h
real*8  :: totmass
real*8, dimension(nc) :: am

integer, dimension(0:nc) :: ntot
real*8  :: treq
real*8  :: dmt
logical :: lbath_clear
integer :: ierror

!-----declare local variables
real*8  :: qomega              ! frequency  for thermostat in [a.u.]
real*8  :: bomega              ! frequency  for thermostat in [a.u.]
real*8  :: wdti_1, dseed, vsqr, fscale
integer :: i, iyosh, ic, ia, ib
real*8  :: b(3,3), val


!-----Set thermostat & barostat masses
qomega = 2.d0*3.141592653589793d0/tomega
bomega = 2.d0*3.141592653589793d0/tbomega
do ic=1,nc
   fack2(ic)=2d0*fack(ic)
end do
gnkt = 3d0*ntot(0)*treq
gkt  =          treq
gnd2kt= gnkt+9d0*gkt
onf   = 1d0/dble(3*ntot(0))
!-----Masses for thermostats
qmass(1)=gnkt/qomega**2
do i = 2, nnos
  qmass(i)=qmass(1)*onf
end do

!-----Total mass of system
totmass = 0.d0
do ic = 1, nc
   totmass = totmass + am(ic)*dble(ntot(ic))
end do
!-----Mass for barostat
!      bmass=gnd2kt/bomega**2
CALL RCIPRL( h0, b, volum0 )
bmass=blkmod/(volum0/totmass)/bomega**2
if( myid == 0 ) then
      write(nfile(1),'(a,2es16.8)') ' '
      write(nfile(1),'(a,2es16.8)')  &
& '  bmass [a.u. / kg^2 m^-4] =', bmass, bmass*(welm*welm*1.d-6/(audang*1.d-10)**4)
      write(nfile(1),'(a,2es16.8)') ' '
      write(nfile(2),'(a,2es16.8)') ' '
      write(nfile(2),'(a,2es16.8)')  &
& '  bmass [a.u. / kg^2 m^-4] =', bmass, bmass*(welm*welm*1.d-6/(audang*1.d-10)**4)
      write(nfile(2),'(a,2es16.8)') ' '
end if


!-----Set matrix
do ib = 1, 3
do ia = 1, 3
   val = 0.d0
   do ic = 1, 3
      val = val + h0(ia,ic)*alpmatrix(ic,ib)
   end do
   b(ia,ib) = val
end do
end do
alpmatrix(1:3,1:3) = b(1:3,1:3)

do ib = 1, 3
do ia = 1, 3
   val = 0.d0
   do ic = 1, 3
      val = val + h0(ia,ic)*betmatrix(ic,ib)
   end do
   b(ia,ib) = val
end do
end do
betmatrix(1:3,1:3) = b(1:3,1:3)


!-----Set higher-order parameters for Yoshida-Suzuki integration
!wdti_1 = dmt/dble(nresn)
wdti_1 = dmt
!if( nyosh == 1 ) then
  wdti(1) = wdti_1
!else if( nyosh == 3 ) then
!  wdti(1) = (1d0/(2d0-2d0**(1d0/3d0)))*wdti_1
!  wdti(2) = 1d0-2d0*wdti(1)
!  wdti(3) = wdti(1)
!else if( nyosh == 5 ) then
!  wdti(1) = (1d0/(4d0-4d0**(1d0/3d0)))*wdti_1
!  wdti(2) = wdti(1)
!  wdti(3) = 1d0-4d0*wdti(1)
!  wdti(4) = wdti(1)
!  wdti(5) = wdti(1)
!else
!  write(nfile(1),*) 'unsupported nyosh selected--now quitting'
!  ierror = 1000
!  return
!end if
!do iyosh = 1, nyosh
  iyosh = 1
  wdti2(iyosh) = 0.5d0*wdti(iyosh)
  wdti4(iyosh) = 0.5d0*wdti2(iyosh)
  wdti8(iyosh) = 0.5d0*wdti4(iyosh)
!end do

!-----Initialize thermostat & barostat states
if( lbath_clear ) then
    do i = 1, nnos
       xlogs(i) = 0d0
       vlogs(i) = 0d0
       glogs(i) = 0d0
    end do
    dseed = 23219d0
    vsqr  = 0d0
    do i = 1, nnos
       call rnd00( vlogs(i), dseed )
       vsqr = vsqr + qmass(i)*vlogs(i)*vlogs(i)
    end do
    fscale = dsqrt(3d0*nnos*gkt/vsqr)
    do i = 1, nnos
       vlogs(i) = vlogs(i)*fscale
    end do

    do ib=1,3
    do ia=1,3
       h(ia,ib,1) = 0.d0
    end do
    end do
end if


return
end




subroutine nhcint( nfile, myid, nodes,  &
& nnos, nresn, nyosh, gnkt, gkt,  &
& xlogs, vlogs, glogs, qmass, wdti, wdti2, wdti4, wdti8,  &
& x, is, n, ntype, fack, h )
!----------------------------------------------------------------------c
!   Nose-Hoover-chain part of time propagation from t=0 to t=dmt/2.
!----------------------------------------------------------------------c
implicit none
integer :: nfile(*), myid, nodes
integer :: nnos, nresn, nyosh
real*8  :: gnkt
real*8  :: gkt
real*8, dimension(nnos)  :: xlogs, vlogs, glogs, qmass
real*8, dimension(nyosh) :: wdti, wdti2, wdti4, wdti8
integer :: n, ntype
real*8  :: x(3*n)
integer :: is(n)
real*8  :: fack(ntype)
real*8  :: h(3,3,0:1)

!-----declare local variables
integer :: nnos1, iresn, iyosh, inos
real*8  :: scale, ekinv, akinv, aanhc
integer :: i, ix, iy, iz


nnos1 = nnos+1
scale = 1d0

!-----calculate the atom kinetic energy, ekin
call getekin( nfile, myid, nodes,  &
& ekinv, x, is, n, ntype, fack, h )
akinv = 2d0*ekinv

!-----Update heat-bath forces
glogs(1) = (akinv-gnkt)/qmass(1)
do i = 1, nnos-1
  glogs(i+1)=(qmass(i)*vlogs(i)*vlogs(i)-gkt)/qmass(i+1)
end do

!-----Start the multiple time step procedure
mtsdo: do iresn = 1, nresn

!-------Yoshida-Suzuki higher-order decomposition steps
  ysdo: do iyosh = 1, nyosh

!---------First half update of the thermostat velocities
    vlogs(nnos)=vlogs(nnos)+glogs(nnos)*wdti4(iyosh)
    do inos = 1, nnos-1
      aanhc = dexp(-wdti8(iyosh)*vlogs(nnos1-inos))
      vlogs(nnos-inos)=( vlogs(nnos-inos)*aanhc  &
&                      + wdti4(iyosh)*glogs(nnos-inos) )*aanhc
    end do

!---------Update the atom-velocity scaling factor
    aanhc = dexp(-wdti2(iyosh)*vlogs(1))
    scale = scale*aanhc

!---------Update the thermostat positions
    do inos = 1, nnos
      xlogs(inos)=xlogs(inos)+vlogs(inos)*wdti2(iyosh)
    end do

!---------Second half update of the thermostat velocities
    glogs(1) = (scale*scale*akinv-gnkt)/qmass(1)
    do inos = 1, nnos-1
      aanhc = dexp(-wdti8(iyosh)*vlogs(inos+1))
      vlogs(inos) = ( vlogs(inos)*aanhc  &
&                  + wdti4(iyosh)*glogs(inos) )*aanhc
      glogs(inos+1)=(qmass(inos)*vlogs(inos)*vlogs(inos)-gkt)  &
&                  /qmass(inos+1)
    end do
    vlogs(nnos) = vlogs(nnos)+glogs(nnos)*wdti4(iyosh)

!-------end do the YS decomposition steps
  end do ysdo

!-----end do the multiple time step procedure
end do mtsdo

!-----Update the atom velocities
do i = 1, n
  iz=3*i
  iy=iz-1
  ix=iy-1
  x(ix) = x(ix)*scale
  x(iy) = x(iy)*scale
  x(iz) = x(iz)*scale
end do


return
end subroutine




subroutine getbathe( nfile, myid, nodes,  &
& bathpe, bathke, nnos, gnkt, gkt, xlogs, vlogs, glogs, qmass )
!----------------------------------------------------------------------c
!     Heat-bath potential & kinetic energies
!----------------------------------------------------------------------c
implicit none
integer :: nfile(*), myid, nodes
real*8  :: bathpe, bathke
integer :: nnos
real*8  :: gnkt
real*8  :: gkt
real*8, dimension(nnos)  :: xlogs, vlogs, glogs, qmass

!-----declare local variables
integer :: inos


!-----Heat-bath potential energy
bathpe = gnkt*xlogs(1)
do inos = 2, nnos
  bathpe = bathpe + gkt*xlogs(inos)
end do

!-----Heat-bath kinetic energy
bathke = 0d0
do inos = 1, nnos
  bathke = bathke+0.5d0*qmass(inos)*vlogs(inos)*vlogs(inos)
end do


return
end subroutine




subroutine nhc_scale( nfile, myid, nodes,  &
& nnos, gnkt, gkt, xlogs, vlogs, glogs, qmass, h )
!----------------------------------------------------------------------c
!     scaling thermostats
!----------------------------------------------------------------------c
implicit none
integer :: nfile(*), myid, nodes
integer :: nnos
real*8  :: gnkt
real*8  :: gkt
real*8, dimension(nnos)  :: xlogs, vlogs, glogs, qmass
real*8  :: h(3,3,0:1)

!-----declare local variables
integer :: i, j
real*8  :: vsqr, fscale


do i = 1, nnos
   xlogs(i) = 0d0
end do
vsqr = 0d0
do i = 1, nnos
   vsqr = vsqr + qmass(i)*vlogs(i)*vlogs(i)
end do
!fscale=dsqrt(3d0*nnos*gkt/vsqr)
fscale=dsqrt(nnos*gkt/vsqr)
do i = 1, nnos
   vlogs(i) = vlogs(i)*fscale
end do

do i = 1, 3
do j = 1, 3
   h(j,i,1) = 0.d0
end do
end do


return
end subroutine




subroutine nhcpfullint( nfile, myid, nodes,  &
& nnos, nresn, nyosh, gnkt, gkt,  &
& xlogs, vlogs, glogs, qmass, wdti, wdti2, wdti4, wdti8,  &
& x, is, n, ntype, fack, fack2, h, volume, pint, pext,  &
& vu, bmass, onf, gnd2kt, irstrct, irstrct_sub, lcell_rstrct )
!----------------------------------------------------------------------c
!   NPT part of time integration from t=0 to t=DT/2.
!----------------------------------------------------------------------c
implicit none
integer :: nfile(*), myid, nodes
integer :: nnos, nresn, nyosh
real*8  :: gnkt
real*8  :: gkt
real*8, dimension(nnos)  :: xlogs, vlogs, glogs, qmass
real*8, dimension(nyosh) :: wdti, wdti2, wdti4, wdti8

integer :: n
real*8,  dimension(3*n) :: x
integer, dimension(n) :: is
integer :: ntype
real*8,  dimension(ntype) :: fack, fack2
real*8,  dimension(3,3,0:1) :: h
real*8  :: volume
real*8,  dimension(3,3) :: pint, pext

real*8,  dimension(3,n) :: vu
real*8  :: bmass, onf, gnd2kt, akinb
integer :: irstrct, irstrct_sub
logical :: lcell_rstrct(3)

!-----declare local variables
real*8,  dimension(3,3) :: vboxgt, vtemp, ubox, gboxg, hi
integer :: nnos1, iresn, iyosh, inos
real*8  :: aanhc
integer :: i, ia, ib
real*8,  dimension(3,3) :: akin, akintotm
real*8  :: akintot
real*8,  dimension(3)   :: veig
real*8,  dimension(3,3) :: veigv
real*8  :: trvg, sc1, sc2, sc3, vx, vy, vz, uv1, uv2, uv3


nnos1 = nnos+1

!----- transpose matrix of vtemp = inverse of h
CALL RCIPRL( h(1,1,0), vtemp, volume )
call mattrp( vtemp, hi )

!-----Transform scaled to unscaled velocities, VU
call vs2vu(n,x,vu,h(1,1,0))

!-----Get the atomic kinetic stress, AKIN, AKINTOT & AKINTOTM
call getkinp( akin, akintotm, akintot,  &
& is, n, ntype, fack2, vu, onf )

!-----Get the box kinetic energy, AKINB
call mattrp(h(1,1,1),vboxgt)
call submatmul(vboxgt,h(1,1,1),ubox)
akinb=bmass*(ubox(1,1)+ubox(2,2)+ubox(3,3))

!-----Update the forces on thermostats, GLOGS, & barostat, GBOXG
glogs(1)=(akintot+akinb-gnd2kt)/qmass(1)
do i=1,nnos-1
  glogs(i+1)=(qmass(i)*vlogs(i)*vlogs(i)-gkt)/qmass(i+1)
end do
do ib=1,3
  do ia=1,3
    gboxg(ia,ib)=( akintotm(ia,ib)+akin(ia,ib)  &
&                 +volume*(pint(ia,ib)-pext(ia,ib)) )/bmass
  end do
end do

!!-----restriction for MD cell
!call restrictMDcell2( irstrct, gboxg, h(1,1,0) )
!-----by symmetry operations
!if( irstrct >= 10 ) call restrictMDcell_sym( gboxg, h(1,1,0) )


!-----Start the multiple time step procedure---------------------------c

!-----MTS loop
mtsdo: do iresn = 1, nresn

!-------Yoshida-Nose factorization loop
  ysdo: do iyosh = 1, nyosh

!---------Update the thermostat velocities, VLOGS
    vlogs(nnos)=vlogs(nnos)+glogs(nnos)*wdti4(iyosh)
    do inos = 1, nnos-1
      aanhc = dexp(-wdti8(iyosh)*vlogs(nnos1-inos))
      vlogs(nnos-inos)=( vlogs(nnos-inos)*aanhc  &
&                      + wdti4(iyosh)*glogs(nnos-inos) )*aanhc
    end do

!---------Update the box velocities, H(,,1)
    aanhc=dexp(-wdti8(iyosh)*vlogs(1))
    do ib=1,3
      do ia=1,3
        h(ia,ib,1)=( h(ia,ib,1)*aanhc  &
&                  + wdti4(iyosh)*gboxg(ia,ib) )*aanhc
      end do
    end do

    !-----restriction for MD cell
!    call restrictMDcell( irstrct, irstrct_sub, h(1,1,0), lcell_rstrct )
!    call restrictMDcell_v( irstrct, irstrct_sub, h(1,1,0), h(1,1,1), lcell_rstrct )

!---------Update the thermostat positions, XLOGS
    do inos = 1, nnos
       xlogs(inos)=xlogs(inos)+vlogs(inos)*wdti2(iyosh)
    end do

!---------Update the atom velocities, VU
    trvg=onf*(h(1,1,1)+h(2,2,1)+h(3,3,1))
    do ib=1,3
      do ia=1,3
        vtemp(ia,ib)=h(ia,ib,1)
        if( ia == ib ) vtemp(ia,ib)=vtemp(ia,ib)+trvg+vlogs(1)
      end do
    end do
    call eigen3(vtemp,veig,veigv)
    sc1=dexp(-veig(1)*wdti2(iyosh))
    sc2=dexp(-veig(2)*wdti2(iyosh))
    sc3=dexp(-veig(3)*wdti2(iyosh))
    do i=1,n
      vx=vu(1,i)
      vy=vu(2,i)
      vz=vu(3,i)
!-----------Inverse orthogonal transformation
      uv1=vx*veigv(1,1)+vy*veigv(2,1)+vz*veigv(3,1)
      uv2=vx*veigv(1,2)+vy*veigv(2,2)+vz*veigv(3,2)
      uv3=vx*veigv(1,3)+vy*veigv(2,3)+vz*veigv(3,3)
!-----------Diagonal Liouville propagation
      uv1=uv1*sc1
      uv2=uv2*sc2
      uv3=uv3*sc3
!-----------Orthogonal transformation
      vu(1,i)=uv1*veigv(1,1)+uv2*veigv(1,2)+uv3*veigv(1,3)
      vu(2,i)=uv1*veigv(2,1)+uv2*veigv(2,2)+uv3*veigv(2,3)
      vu(3,i)=uv1*veigv(3,1)+uv2*veigv(3,2)+uv3*veigv(3,3)
    end do

!---------Get the atomic kinetic stress, AKIN, AKINTOT & AKINTOTM
    call getkinp( akin, akintotm, akintot,  &
& is, n, ntype, fack2, vu, onf )

!---------Update the box forces, GBOXG
    do ib=1,3
      do ia=1,3
        gboxg(ia,ib)=( akintotm(ia,ib)+akin(ia,ib)  &
&                     +volume*(pint(ia,ib)-pext(ia,ib)) )/bmass
      end do
    end do

!    !-----restriction for MD cell
!    call restrictMDcell2( irstrct, gboxg, h(1,1,0) )
    !-----by symmetry operations
!    if( irstrct >= 10 ) call restrictMDcell_sym( gboxg, h(1,1,0) )

!---------Update the box velocities, H(,,1)
    aanhc=dexp(-wdti8(iyosh)*vlogs(1))
    do ib=1,3
      do ia=1,3
        h(ia,ib,1)=( h(ia,ib,1)*aanhc  &
&                  + wdti4(iyosh)*gboxg(ia,ib) )*aanhc
      end do
    end do

    !-----restriction for MD cell
!    call restrictMDcell( irstrct, irstrct_sub, h(1,1,0), lcell_rstrct )
!    call restrictMDcell_v( irstrct, irstrct_sub, h(1,1,0), h(1,1,1), lcell_rstrct )

!---------Update the thermostat forces, GLOGS, & velocities, VLOGS
    call mattrp(h(1,1,1),vboxgt)
    call submatmul(vboxgt,h(1,1,1),ubox)
    akinb=bmass*(ubox(1,1)+ubox(2,2)+ubox(3,3))

    glogs(1)=(akintot+akinb-gnd2kt)/qmass(1)
    do inos=1,nnos-1
      aanhc=dexp(-wdti8(iyosh)*vlogs(inos+1))
      vlogs(inos)=( vlogs(inos)*aanhc  &
&                 + wdti4(iyosh)*glogs(inos) )*aanhc
      glogs(inos+1)=(qmass(inos)*vlogs(inos)*vlogs(inos)-gkt)  &
&                  /qmass(inos+1)
    end do
    vlogs(nnos)=vlogs(nnos)+glogs(nnos)*wdti4(iyosh)

!-------end do Yoshida-Suzuki higher-order loop
  end do ysdo

!-----end do MTS loop
end do mtsdo

!-----End the multiple time step procedure-----------------------------c

!-----Update scaled velocities
call vs2vu(n,vu,x,hi)


return
end




subroutine getbathboxe( nfile, myid, nodes,  &
& bathpe, bathke, boxpe, boxke,  &
& nnos, gnkt, gkt, xlogs, vlogs, glogs, qmass,  &
& bmass, gnd2kt, pext, h, volume )
!----------------------------------------------------------------------c
implicit none
integer :: nfile(*), myid, nodes
real*8  :: bathpe, bathke, boxpe, boxke
integer :: nnos
real*8  :: gnkt, gkt, gnd2kt
real*8, dimension(nnos) :: xlogs, vlogs, glogs, qmass
real*8, dimension(3,3)  :: pext
real*8, dimension(3,3,0:1)  :: h
real*8  :: bmass, volume, akinb

!-----declare local variables
integer :: inos
real*8, dimension(3,3)  :: vboxgt, ubox


!-----Calculate the extended potential energies, BATHPE & BOXPE
bathpe=gnd2kt*xlogs(1)
do inos=2,nnos
  bathpe=bathpe+gkt*xlogs(inos)
end do
boxpe=(pext(1,1)+pext(2,2)+pext(3,3))*volume/3d0

!-----Calculate the extended kinetic energies, BATHKE & BOXKE
bathke=0d0
do inos=1,nnos
  bathke=bathke+qmass(inos)*vlogs(inos)*vlogs(inos)
end do
bathke=0.5d0*bathke
!-----Get the box kinetic energy, AKINB
call mattrp(h(1,1,1),vboxgt)
call submatmul(vboxgt,h(1,1,1),ubox)
akinb=bmass*(ubox(1,1)+ubox(2,2)+ubox(3,3))
boxke=0.5d0*akinb


return
end




subroutine getkinp( akin, akintotm, akintot,  &
& is, n, ntype, fack2, vu, onf )
!-----------------------------------------------------------------------
!  Computes the atomic kinetic stress tensor, AKIN, twice the atomic 
!  kinetic energy, AKINTOT, & scaled, diagonal atomic-kinetic-stress 
!  tensor, AKINTOTM.
!-----------------------------------------------------------------------
implicit none
real*8,  dimension(3,3) :: akin, akintotm
real*8  :: akintot
integer :: n, ntype
integer, dimension(n) :: is
real*8,  dimension(ntype) :: fack2
real*8,  dimension(3,n) :: vu
real*8  :: onf

!-----declare local variables
real*8,  dimension(9) :: dbuf
integer :: ia, ib, i


!-----Get the atomic kinetic stress tensor, AKIN
do ib=1,3
  do ia=ib,3
    akin(ia,ib)=0d0
    do i=1,n
      akin(ia,ib)=akin(ia,ib)+fack2(is(i))*vu(ia,i)*vu(ib,i)
    end do
    akin(ib,ia)=akin(ia,ib)
  end do
end do
call gdsum(akin,9,dbuf)

!-----Total kinetic energy, AKINTOT
akintot=akin(1,1)+akin(2,2)+akin(3,3)

!-----Scaled, diagonal atomic-kinetic-stress tensor, AKINTOTM
do ib=1,3
  do ia=1,3
    if( ia == ib ) then
      akintotm(ia,ib)=onf*akintot
    else
      akintotm(ia,ib)=0d0
    end if
  end do
end do

return
end




subroutine msstint( nfile, myid, nodes,  &
& x, is, n, ntype, fack2, h, volume, xmsst, vxmsst,  &
& wdti, wdti2, wdti4, wdti8, vu, pint, pext, bmass,  &
& volum0, alpmatrix, betmatrix, totmass, shockspeed )
!----------------------------------------------------------------------c
!   MSST part of time integration from t=0 to t=DT/2.
!----------------------------------------------------------------------c
implicit none
integer :: nfile(*), myid, nodes
integer :: n
real*8,  dimension(3*n) :: x
integer, dimension(n) :: is
integer :: ntype
real*8,  dimension(ntype) :: fack2
real*8,  dimension(3,3,0:1) :: h
real*8  :: volume
real*8  :: xmsst, vxmsst
real*8  :: wdti, wdti2, wdti4, wdti8
real*8,  dimension(3,n) :: vu
real*8,  dimension(3,3) :: pint, pext
real*8  :: bmass, volum0, totmass, shockspeed
real*8  :: alpmatrix(3,3), betmatrix(3,3)

!-----declare local variables
integer :: ia, ib
real*8  :: hit(3,3)
real*8  :: fval, hx, gx, xdot


!----- transpose matrix of hit = inverse of h
CALL RCIPRL( h(1,1,0), hit, volume )
!call mattrp( hit, hi )  ! hi is the inverse of h

!-----Transform scaled to unscaled velocities, VU
call vs2vu(n,x,vu,h(1,1,0))

call msst_getgxhx( nfile, myid, nodes,  &
& is, n, ntype, fack2, h, volume, vu, pint, pext, bmass,  &
& volum0, alpmatrix, totmass, shockspeed, hit, fval, hx, gx )

!--- exp(i L4 * dt/2)
vxmsst = vxmsst + gx/bmass * wdti2

!--- time derivative of cell vectors
xdot = totmass/(fval*fval)*vxmsst
do ib = 1, 3
do ia = 1, 3
   h(ia,ib,1) = alpmatrix(ia,ib)*xdot
end do
end do

!!-----Update scaled velocities
!call vs2vu(n,vu,x,hi)


return
end




subroutine msst_getgxhx( nfile, myid, nodes,  &
& is, n, ntype, fack2, h, volume, vu, pint, pext, bmass,  &
& volum0, alpmatrix, totmass, shockspeed, hi, fval, hx, gx )
!----------------------------------------------------------------------c
!   MSST part of time integration from t=0 to t=DT/2.
!----------------------------------------------------------------------c
implicit none
integer :: nfile(*), myid, nodes
integer :: n
integer, dimension(n) :: is
integer :: ntype
real*8,  dimension(ntype) :: fack2
real*8,  dimension(3,3,0:1) :: h
real*8  :: volume
real*8,  dimension(3,n) :: vu
real*8,  dimension(3,3) :: pint, pext
real*8  :: bmass, volum0, totmass, shockspeed
real*8  :: alpmatrix(3,3)
real*8  :: hi(3,3), fval, hx, gx

!-----declare local variables
integer :: i, ia, ib
real*8  :: akin(3,3), akintotm(3,3), gboxg(3,3)
real*8  :: hmat(3,3), gmat(3,3)
real*8  :: akintot, vratio, val


!-----Get the atomic kinetic stress, AKIN, AKINTOT & AKINTOTM
call getkinp( akin, akintotm, akintot,  &
& is, n, ntype, fack2, vu, 0.d0 )

do ib = 1, 3
   do ia = 1, 3
      gboxg(ia,ib) = akin(ia,ib) + volume*(pint(ia,ib) - pext(ia,ib))
   end do
end do
vratio = volume/volum0
val    = totmass*shockspeed*shockspeed*(vratio-1.d0)*vratio
! 6/27/2014 to avoid volume expansion
val    = val * sign(1.d0,volum0-volume)
do ia = 1, 3
   gboxg(ia,ia) = gboxg(ia,ia) + val
end do

!---gmatrix
do ib = 1, 3
do ia = 1, 3
   val = 0.d0
   do i = 1, 3
      val = val + gboxg(ia,i)*hi(i,ib)
   end do
   gmat(ia,ib) = val
end do
end do

fval = 0.d0
do ib = 1, 3
do ia = 1, 3
   fval = fval + hi(ia,ib)*alpmatrix(ia,ib)
end do
end do
fval = fval*volume

!---hmatrix
do ib = 1, 3
do ia = 1, 3
   val = 0.d0
   do i = 1, 3
      val = val + alpmatrix(i,ia)*hi(i,ib)
   end do
   gboxg(ia,ib) = val
end do
end do
do ib = 1, 3
do ia = 1, 3
   val = 0.d0
   do i = 1, 3
      val = val + hi(ia,i)*gboxg(i,ib)
   end do
   hmat(ia,ib) = val*volume
end do
end do
val = bmass*totmass/(fval*fval*fval)
do ib = 1, 3
do ia = 1, 3
   hmat(ia,ib) = val*(fval*hi(ia,ib) - hmat(ia,ib))
end do
end do

!---g(x) and h(x)
gx = 0.d0
hx = 0.d0
do ib = 1, 3
do ia = 1, 3
   gx = gx + gmat(ia,ib)*alpmatrix(ia,ib)
   hx = hx + hmat(ia,ib)*alpmatrix(ia,ib)
end do
end do


return
end




subroutine getmsstboxe( nfile, myid, nodes,  &
& boxpe, boxke, h, volume, bmass, pext, &
& vxmsst, volum0, alpmatrix, totmass, shockspeed )
!----------------------------------------------------------------------c
implicit none
integer :: nfile(*), myid, nodes
real*8  :: boxpe, boxke
real*8, dimension(3,3,0:1)  :: h
real*8  :: volume, bmass
real*8, dimension(3,3)  :: pext
real*8  :: vxmsst, volum0, alpmatrix(3,3), totmass, shockspeed

!-----declare local variables
real*8  :: ubox(3,3), fval, vratio
integer :: ia, ib


!----- transpose matrix of ubox = inverse of h
CALL RCIPRL( h(1,1,0), ubox, volume )
!call mattrp( ubox, hi )

fval = 0.d0
do ib = 1, 3
do ia = 1, 3
   fval = fval + ubox(ia,ib)*alpmatrix(ia,ib)
end do
end do
fval = fval*volume

!-----Calculate the extended kinetic energies, BOXKE
boxke = 0.5d0*bmass*totmass/(fval*fval)*vxmsst*vxmsst

!-----Calculate the extended potential energies, BOXPE
vratio = shockspeed*(volume/volum0 - 1.d0)
!boxpe = (pext(1,1)+pext(2,2)+pext(3,3))*(volume-volum0)/3d0  &
!&     - 0.5d0*totmass*vratio*vratio
! 6/27/2014 to avoid volume expansion
boxpe = (pext(1,1)+pext(2,2)+pext(3,3))*(volume-volum0)/3d0  &
&     - 0.5d0*totmass*vratio*vratio * sign(1.d0,volum0-volume)


return
end




subroutine msst_scale( nfile, myid, nodes, vxmsst, h )
!----------------------------------------------------------------------c
!     scaling thermostats
!----------------------------------------------------------------------c
implicit none
integer :: nfile(*), myid, nodes
real*8  :: vxmsst
real*8, dimension(3,3,0:1) :: h

vxmsst = 0.d0

h(1:3,1:3,1) = 0.d0

return
end subroutine




subroutine xs2xu(n,xs,xu,h0,xog,yog,zog)
!----------------------------------------------------------------------c
!  Transforms from scaled to unscaled coordinates.
!  xog,yog,zog = origin of node in scaled coordinates.
!----------------------------------------------------------------------c
implicit none
integer :: n
real*8,  dimension(3,n) :: xs, xu
real*8,  dimension(3,3) :: h0
real*8  :: xog,yog,zog

!-----declare local variables
real*8  :: ax, bx, cx, ay, by, cy, az, bz, cz, sx, sy, sz
integer :: i

ax=h0(1,1)
bx=h0(1,2)
cx=h0(1,3)
ay=h0(2,1)
by=h0(2,2)
cy=h0(2,3)
az=h0(3,1)
bz=h0(3,2)
cz=h0(3,3)
do i=1,n
  sx=xs(1,i)+xog
  sy=xs(2,i)+yog
  sz=xs(3,i)+zog
  xu(1,i)=ax*sx+bx*sy+cx*sz
  xu(2,i)=ay*sx+by*sy+cy*sz
  xu(3,i)=az*sx+bz*sy+cz*sz
end do

return
end




subroutine xu2xs(n,xs,xu,h0,xog,yog,zog)
!----------------------------------------------------------------------c
!  Transforms from scaled to unscaled coordinates.
!  xog,yog,zog = origin of node in scaled coordinates.
!----------------------------------------------------------------------c
implicit none
integer :: n
real*8,  dimension(3,n) :: xs, xu
real*8,  dimension(3,3) :: h0
real*8  :: xog,yog,zog

!-----declare local variables
real*8  :: ax, bx, cx, ay, by, cy, az, bz, cz, sx, sy, sz
integer :: i

ax=h0(1,1)
bx=h0(1,2)
cx=h0(1,3)
ay=h0(2,1)
by=h0(2,2)
cy=h0(2,3)
az=h0(3,1)
bz=h0(3,2)
cz=h0(3,3)
do i=1,n
  sx=xs(1,i)
  sy=xs(2,i)
  sz=xs(3,i)
  xu(1,i)=ax*sx+bx*sy+cx*sz - xog
  xu(2,i)=ay*sx+by*sy+cy*sz - yog
  xu(3,i)=az*sx+bz*sy+cz*sz - zog
end do

return
end




subroutine vs2vu(n,vs,vu,h0)
!-----------------------------------------------------------------------
!  Transforms from scaled to unscaled velocities.
!-----------------------------------------------------------------------
implicit none
integer :: n
real*8,  dimension(3,n) :: vs, vu
real*8,  dimension(3,3) :: h0

!-----declare local variables
real*8  :: ax, bx, cx, ay, by, cy, az, bz, cz, sx, sy, sz
integer :: i


ax=h0(1,1)
bx=h0(1,2)
cx=h0(1,3)
ay=h0(2,1)
by=h0(2,2)
cy=h0(2,3)
az=h0(3,1)
bz=h0(3,2)
cz=h0(3,3)
do i = 1, n
   sx=vs(1,i)
   sy=vs(2,i)
   sz=vs(3,i)
   vu(1,i)=ax*sx+bx*sy+cx*sz
   vu(2,i)=ay*sx+by*sy+cy*sz
   vu(3,i)=az*sx+bz*sy+cz*sz
end do


return
end




subroutine submatmul(a,b,c)
!----------------------------------------------------------------------c
!  Multiplies 3-by-3 matrices a & b, and stores the result in c.
!----------------------------------------------------------------------c
implicit none
real*8,  dimension(3,3) :: a, b, c
integer :: i, j, k

do j=1,3
do i=1,3
   c(i,j)=0d0
   do k=1,3
      c(i,j) = c(i,j) + a(i,k)*b(k,j)
   end do
end do
end do

return
end 


subroutine mattrp(a,at)
!----------------------------------------------------------------------c
!  Transposes a 3-by-3 matrix, a, and stores the result in at.
!----------------------------------------------------------------------c
implicit none
real*8,  dimension(3,3) :: a, at
integer :: i, j

do j = 1, 3
do i = 1, 3
   at(i,j) = a(j,i)
end do
end do

return
end 




subroutine eigen3(a,d,v)
!----------------------------------------------------------------------c
!  Diagonalizes a real symmetric 3x3 matrix.
!     a(3,3): Input matrix
!     d(3):   Return eigenvalues
!     v(3,3): Return eigenvectors--v(*,i) = i-th eigenvector
!----------------------------------------------------------------------c
implicit none
real*8, dimension(3,3) :: a, v
real*8, dimension(3)   :: d
integer :: i,j

!-----Find eigenvalues
call jacobi(a,3,3,d,v,i)
!-----Restore the original matrix
do i=1,3
  do j=i+1,3
     a(i,j)=a(j,i)
  end do
end do

return
end


subroutine jacobi(a,n,np,d,v,nrot)
!----------------------------------------------------------------------c
!  Diagonalizes a real symmetric matrix by Jacobi transformation.
!  From "Numerical Recipes".
!----------------------------------------------------------------------c
implicit none
integer :: n,np,nrot
real*8, dimension(np,np) :: a, v
real*8, dimension(np)    :: d
integer, parameter :: nmax = 500
integer :: i,ip,iq,j
real*8  :: c,g,h,s,sm,t,tau,theta,tresh
real*8, dimension(nmax) :: b, z
real*8, parameter :: small = 1.d-30

do ip=1,n
  do iq=1,n
    v(ip,iq)=0.d0
  end do
  v(ip,ip)=1.d0
end do

do ip=1,n
  b(ip)=a(ip,ip)
  d(ip)=b(ip)
  z(ip)=0.d0
end do

nrot=0
do i=1,50
  sm=0.
  do ip=1,n-1
    do iq=ip+1,n
      sm=sm+abs(a(ip,iq))
    end do
  end do
  if( abs(sm) < small ) return
  if( i < 4 )then
    tresh=0.2d0*sm/n**2
  else
    tresh=0.d0
  end if
  do ip=1,n-1
    do iq=ip+1,n
      g=100.d0*abs(a(ip,iq))
      if( i > 4 .and. abs(abs(d(ip))+g-abs(d(ip))) < small .and.  &
&                     abs(abs(d(iq))+g-abs(d(iq))) < small ) then
        a(ip,iq)=0.d0
      else if( abs(a(ip,iq)) > tresh )then
        h=d(iq)-d(ip)
        if( abs(abs(h)+g-abs(h)) < small )then
          if( abs(h) < small ) cycle   ! added on 2015/05/09
          t=a(ip,iq)/h
        else
          theta=0.5d0*h/a(ip,iq)
          t=1.d0/(abs(theta)+sqrt(1.d0+theta**2))
          if( theta < 0.d0 ) t=-t
        end if
        c=1.d0/sqrt(1.d0+t**2)
        s=t*c
        tau=s/(1.d0+c)
        h=t*a(ip,iq)
        z(ip)=z(ip)-h
        z(iq)=z(iq)+h
        d(ip)=d(ip)-h
        d(iq)=d(iq)+h
        a(ip,iq)=0.
        do j=1,ip-1
          g=a(j,ip)
          h=a(j,iq)
          a(j,ip)=g-s*(h+g*tau)
          a(j,iq)=h+s*(g-h*tau)
        end do
        do j=ip+1,iq-1
          g=a(ip,j)
          h=a(j,iq)
          a(ip,j)=g-s*(h+g*tau)
          a(j,iq)=h+s*(g-h*tau)
        end do
        do j=iq+1,n
          g=a(ip,j)
          h=a(iq,j)
          a(ip,j)=g-s*(h+g*tau)
          a(iq,j)=h+s*(g-h*tau)
        end do
        do j=1,n
          g=v(j,ip)
          h=v(j,iq)
          v(j,ip)=g-s*(h+g*tau)
          v(j,iq)=h+s*(g-h*tau)
        end do
        nrot=nrot+1
      end if
    end do
  end do
  do ip=1,n
    b(ip)=b(ip)+z(ip)
    d(ip)=b(ip)
    z(ip)=0.d0
  end do
end do
stop 'too many iterations in jacobi'
return
end




subroutine getekin( nfile, myid, nodes,  &
& ekinv, v, is, n, ntype, fack, h )
!----------------------------------------------------------------------c
!  Calculates global kinetic energy.
!----------------------------------------------------------------------c
implicit none
integer :: nfile(*), myid, nodes
integer :: n, ntype
real*8  :: v(3,n)
integer :: is(n)
real*8  :: fack(ntype)
real*8  :: h(3,3,0:1)
real*8  :: ekinv

!-----declare local variables
integer :: i
real*8  :: ekinl, svx, svy, svz, vx, vy, vz, vtemp


ekinl = 0d0
do i = 1, n
   svx=v(1,i)
   svy=v(2,i)
   svz=v(3,i)
   vx=h(1,1,0)*svx+h(1,2,0)*svy+h(1,3,0)*svz
   vy=h(2,1,0)*svx+h(2,2,0)*svy+h(2,3,0)*svz
   vz=h(3,1,0)*svx+h(3,2,0)*svy+h(3,3,0)*svz
   ekinl=ekinl+fack(is(i))*(vx**2+vy**2+vz**2)
end do
ekinv = ekinl
call gdsum(ekinv,1,vtemp)


return
end




subroutine vscale( nfile, myid, nodes,  &
& treq, ekin, v, is, n, ntype, ntot )
!-----------------------------------------------------------------------
!     velocity scaling
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
real*8  :: treq, ekin
integer :: n
real*8,  dimension(3*n) :: v
integer, dimension(n)   :: is
integer :: ntype
integer, dimension(0:ntype) :: ntot

!-----declare local variables
integer :: i
real*8  :: tempc, fract


tempc = 2.d0*ekin/dble(3*ntot(0))
fract = sqrt(treq/tempc)
do i = 1, 3*n
   v(i) = fract*v(i)
end do


return
end




subroutine vkick( nfile, myid, nodes,  &
& v, a, is, n, ntype, lfixion )
!-----------------------------------------------------------------------
!     half-step kick, v(t+DT/2)
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: n
integer :: ntype
real*8,  dimension(3*n) :: v
real*8,  dimension(3*n) :: a
integer, dimension(n) :: is
logical, dimension(ntype) :: lfixion

!-----declare local variables
integer :: i, it

do i = 1, 3*n
it = is((i+2)/3)
if( .not.lfixion(it) ) then
   v(i) = v(i) + a(i)
end if
end do

return
end




subroutine pvv( nfile, myid, nodes,  &
& v, a, is, n, ntype, lfixion )
!-----------------------------------------------------------------------
!     projection of velocities on atomic forces
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: n
integer :: ntype
real*8,  dimension(3,n) :: v
real*8,  dimension(3,n) :: a
integer, dimension(n) :: is
logical, dimension(ntype) :: lfixion

!-----declare local variables
integer :: i, it
real*8 :: f, vf, v1, v2


do i = 1, n
   it = is(i)
   if( .not.lfixion(it) ) then
       vf = v(1,i)*a(1,i)+v(2,i)*a(2,i)+v(3,i)*a(3,i)
       f  = a(1,i)**2 + a(2,i)**2 + a(3,i)**2
       if( vf < 0.d0 .or. abs(f) < 1.d-30 ) then
           v(:,i) = 0.d0
         else
           v(:,i) = vf*a(:,i)/f
       end if
   end if
end do


return
end




subroutine nhc_position_updt( nfile, myid, nodes,  &
& x, v, is, n, ntype, lfixion, h, volume, dt, xu, vu,  &
& sxog, syog, szog, irstrct, irstrct_sub, lcell_rstrct )
!----------------------------------------------------------------------c
!  Atom & box position update in the modified Verlet integration in the
!  NPT dynamics.
!----------------------------------------------------------------------c
implicit none
integer :: nfile(*), myid, nodes
integer :: n
integer :: ntype
real*8,  dimension(3*n) :: x
real*8,  dimension(3*n) :: v
integer, dimension(n) :: is
logical, dimension(ntype) :: lfixion
real*8,  dimension(3,3,0:1) :: h
real*8  :: dt, volume
real*8,  dimension(3,n) :: xu, vu
real*8  :: sxog, syog, szog
integer :: irstrct, irstrct_sub
logical :: lcell_rstrct(3)

!-----declare local variables
real*8,  dimension(3)   :: adum2, bdum
real*8,  dimension(3,3) :: veigvt, veigv, ubox, hi
real*8,  dimension(3)   :: veig
real*8,  parameter :: e2=1d0/6d0,e4=e2/20d0,e6=e4/42d0,e8=e6/72d0
real*8  :: dt2, pp, adum, arg2, poly
integer :: ia, ib, i
real*8  :: rx, ry, rz, vx, vy, vz, u1, u2, u3, uv1, uv2, uv3


!-----Update diagonal box-velocity tensors
call eigen3(h(1,1,1),veig,veigv)
dt2=dt/2d0
do ia=1,3
  pp=dt2*veig(ia) 
  adum=dexp(pp)
  adum2(ia)=dexp(pp*2d0)
  arg2=pp*pp 
  poly=(((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1d0
  bdum(ia)=adum*poly
end do

!-----Transform scale to unscaled atomic coordinates & velocities
call xs2xu(n,x,xu,h(1,1,0),sxog,syog,szog)
call vs2vu(n,v,vu,h(1,1,0))

!-----Update atomic positions
do i=1,n
  rx=xu(1,i)
  ry=xu(2,i)
  rz=xu(3,i)
  vx=vu(1,i)
  vy=vu(2,i)
  vz=vu(3,i)
!-------Inverse orthogonal transformation of coordinate
  u1=rx*veigv(1,1)+ry*veigv(2,1)+rz*veigv(3,1)
  u2=rx*veigv(1,2)+ry*veigv(2,2)+rz*veigv(3,2)
  u3=rx*veigv(1,3)+ry*veigv(2,3)+rz*veigv(3,3)
!-------Inverse orthogonal transformation of velocity
  uv1=vx*veigv(1,1)+vy*veigv(2,1)+vz*veigv(3,1)
  uv2=vx*veigv(1,2)+vy*veigv(2,2)+vz*veigv(3,2)
  uv3=vx*veigv(1,3)+vy*veigv(2,3)+vz*veigv(3,3)
!-------Apply the box-velocity Liouville operator
  u1=u1*adum2(1)+uv1*bdum(1)
  u2=u2*adum2(2)+uv2*bdum(2)
  u3=u3*adum2(3)+uv3*bdum(3)
!-------Orthogonal transformation of coordinate
  xu(1,i)=u1*veigv(1,1)+u2*veigv(1,2)+u3*veigv(1,3)
  xu(2,i)=u1*veigv(2,1)+u2*veigv(2,2)+u3*veigv(2,3)
  xu(3,i)=u1*veigv(3,1)+u2*veigv(3,2)+u3*veigv(3,3)
end do

!-----Update box tensor
call mattrp(veigv,veigvt)
!-----Inverse orthogonal transformation of box tensor
call submatmul(veigvt,h(1,1,0),ubox)
!-----Apply the box-velocity Liouville operator
do ib=1,3
  do ia=1,3
    ubox(ia,ib)=adum2(ia)*ubox(ia,ib)
  end do
end do
!-----Orthogonal transformation of box tensor
call submatmul(veigv,ubox,h(1,1,0))

!-----restriction for MD cell
!call restrictMDcell( irstrct, irstrct_sub, h(1,1,0), lcell_rstrct )

!-----Update volume (box-related tensors, if necessary)
call RCIPRL( h(1,1,0), ubox, volume )
call mattrp(ubox,hi)

!-----Update scaled atomic coordinates & velocities
call xu2xs(n,xu,x,hi,sxog,syog,szog)
!-----Physical velocities unchanged, but box changed
call vs2vu(n,vu,v,hi)


return
end




subroutine msst_position_updt( nfile, myid, nodes,  &
& x, v, is, n, ntype, lfixion, xu, vu,  &
& sxog, syog, szog, h, volume, dt,  &
& xmsst, vxmsst, alpmatrix, betmatrix, totmass )
!----------------------------------------------------------------------c
!  Atom & box position update in the modified Verlet integration in the
!  MSST dynamics.
!----------------------------------------------------------------------c
implicit none
integer :: nfile(*), myid, nodes
integer :: n
integer :: ntype
real*8,  dimension(3*n) :: x
real*8,  dimension(3*n) :: v
integer, dimension(n) :: is
logical, dimension(ntype) :: lfixion
real*8,  dimension(3,n) :: xu, vu
real*8  :: sxog, syog, szog
real*8,  dimension(3,3,0:1) :: h
real*8  :: dt, volume
real*8  :: xmsst, vxmsst, alpmatrix(3,3), betmatrix(3,3), totmass

!-----declare local variables
real*8  :: fval, xdot
real*8  :: ubox(3,3), gmat(3,3), gimat(3,3)
real*8  :: px, py, pz
integer :: ia, ib, i, it, ix, iy, iz


!---gmat = h^t * h
gmat(:,:) = 0.d0
do ib = 1, 3
do ia = 1, 3
   fval = 0.d0
   do i = 1, 3
      fval = fval + h(i,ia,0)*h(i,ib,0)
   end do
   gmat(ia,ib) = fval
end do
end do


!----- transpose matrix of ubox = inverse of h
CALL RCIPRL( h(1,1,0), ubox, volume )

fval = 0.d0
do ib = 1, 3
do ia = 1, 3
   fval = fval + ubox(ia,ib)*alpmatrix(ia,ib)
end do
end do
fval = fval*volume

!--- exp(i L3 * dt/2)
xdot  = totmass/(fval*fval)*vxmsst
xmsst = xmsst + xdot * dt * 0.5d0

!--- new cell vectors and their time derivatives
do ib = 1, 3
do ia = 1, 3
   h(ia,ib,0) = alpmatrix(ia,ib)*xmsst + betmatrix(ia,ib)
   h(ia,ib,1) = alpmatrix(ia,ib)*xdot
end do
end do


!----- transpose matrix of ubox = inverse of h
CALL RCIPRL( h(1,1,0), ubox, volume )

!---gimat = h^-1 * h^-1^t
gimat(:,:) = 0.d0
do ib = 1, 3
do ia = 1, 3
   fval = 0.d0
   do i = 1, 3
      fval = fval + ubox(i,ia)*ubox(i,ib)
   end do
   gimat(ia,ib) = fval
end do
end do


do i = 1, n
   ix = 3*i-2
   iy = 3*i-1
   iz = 3*i
   px = gmat(1,1)*v(ix) + gmat(1,2)*v(iy) + gmat(1,3)*v(iz)
   py = gmat(2,1)*v(ix) + gmat(2,2)*v(iy) + gmat(2,3)*v(iz)
   pz = gmat(3,1)*v(ix) + gmat(3,2)*v(iy) + gmat(3,3)*v(iz)
   v(ix) = gimat(1,1)*px + gimat(1,2)*py + gimat(1,3)*pz
   v(iy) = gimat(2,1)*px + gimat(2,2)*py + gimat(2,3)*pz
   v(iz) = gimat(3,1)*px + gimat(3,2)*py + gimat(3,3)*pz
end do

do i = 1, 3*n
it = is((i+2)/3)
if( .not.lfixion(it) ) then
   x(i) = x(i) + v(i)
!         if( x(i) <  0.d0 ) x(i) = x(i) + 1.d0
!         if( x(i) >= 1.d0 ) x(i) = x(i) - 1.d0
end if
end do


!---gmat = h^t * h
gmat(:,:) = 0.d0
do ib = 1, 3
do ia = 1, 3
   fval = 0.d0
   do i = 1, 3
      fval = fval + h(i,ia,0)*h(i,ib,0)
   end do
   gmat(ia,ib) = fval
end do
end do

fval = 0.d0
do ib = 1, 3
do ia = 1, 3
   fval = fval + ubox(ia,ib)*alpmatrix(ia,ib)
end do
end do
fval = fval*volume

!--- exp(i L3 * dt/2)
xdot  = totmass/(fval*fval)*vxmsst
xmsst = xmsst + xdot * dt * 0.5d0


!--- new cell vectors and their time derivatives
do ib = 1, 3
do ia = 1, 3
   h(ia,ib,0) = alpmatrix(ia,ib)*xmsst + betmatrix(ia,ib)
   h(ia,ib,1) = alpmatrix(ia,ib)*xdot
end do
end do

!----- transpose matrix of ubox = inverse of h
CALL RCIPRL( h(1,1,0), ubox, volume )

!---gimat = h^-1 * h^-1^t
gimat(:,:) = 0.d0
do ib = 1, 3
do ia = 1, 3
   fval = 0.d0
   do i = 1, 3
      fval = fval + ubox(i,ia)*ubox(i,ib)
   end do
   gimat(ia,ib) = fval
end do
end do


do i = 1, n
   ix = 3*i-2
   iy = 3*i-1
   iz = 3*i
   px = gmat(1,1)*v(ix) + gmat(1,2)*v(iy) + gmat(1,3)*v(iz)
   py = gmat(2,1)*v(ix) + gmat(2,2)*v(iy) + gmat(2,3)*v(iz)
   pz = gmat(3,1)*v(ix) + gmat(3,2)*v(iy) + gmat(3,3)*v(iz)
   v(ix) = gimat(1,1)*px + gimat(1,2)*py + gimat(1,3)*pz
   v(iy) = gimat(2,1)*px + gimat(2,2)*py + gimat(2,3)*pz
   v(iz) = gimat(3,1)*px + gimat(3,2)*py + gimat(3,3)*pz
end do


return
end




subroutine updtps( nfile, myid, nodes,  &
& x, v, is, n, ntype, lfixion )
!-----------------------------------------------------------------------
!     Coordinate update, x(t+Dt)
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: n
integer :: ntype
real*8  :: x(3*n), v(3*n)
integer :: is(n)
logical :: lfixion(ntype)

!-----declare local variables
integer :: i, it

do i = 1, 3*n
it = is((i+2)/3)
if( .not.lfixion(it) ) then
   x(i) = x(i) + v(i)
!         if( x(i) <  0.d0 ) x(i) = x(i) + 1.d0
!         if( x(i) >= 1.d0 ) x(i) = x(i) - 1.d0
end if
end do

return
end




subroutine updtps_pbc( nfile, myid, nodes,  &
& x, is, n, ntype )
!-----------------------------------------------------------------------
!     Coordinate update, x(t+Dt)
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: n
integer :: ntype
real*8  :: x(3*n)
integer :: is(n)

!-----declare local variables
integer :: i

do i = 1, 3*n
   if( x(i) <  0.d0 ) x(i) = x(i) - int(x(i)) + 1.d0
   if( x(i) >= 1.d0 ) x(i) = x(i) - int(x(i))
end do

return
end




subroutine updtcg( nfile, myid, nodes,  &
& x, is, n, ntype, frc, watom, h, dtmd, pixcg, fnorm, acon,  &
& lfixion )
!-----------------------------------------------------------------------
!     Coordinate update, x(t+Dt) by CG
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: n
integer :: ntype
real*8,  dimension(3*n) :: x
integer, dimension(n) :: is
real*8,  dimension(3*n) :: frc
real*8,  dimension(ntype) :: watom
real*8,  dimension(3,3,0:1) :: h
real*8  :: dtmd
real*8,  dimension(3*n) :: pixcg
real*8,  dimension(n)   :: fnorm
real*8,  dimension(ntype) :: acon
logical, dimension(ntype) :: lfixion

!-----declare local variables
integer :: it, i, ix, iy, iz
real*8  :: twatom, frcm1, frcm2, frcm3, fnormc, betafc


twatom = 0.d0
do it = 1, ntype
   if( .not.lfixion(it) ) twatom = twatom + watom(it)
end do

do i = 1, n
   it = is(i)
   ix = 3*i-2
   iy = 3*i-1
   iz = 3*i
if( .not.lfixion(it) ) then
   frcm1 = frc(ix)/(watom(it)/twatom) / acon(it)
   frcm2 = frc(iy)/(watom(it)/twatom) / acon(it)
   frcm3 = frc(iz)/(watom(it)/twatom) / acon(it)
   fnormc = frcm1*frcm1 + frcm2*frcm2 + frcm3*frcm3
   if( fnormc < fnorm(i) ) then
       betafc = fnormc/fnorm(i)
     else
       betafc = 0.d0
   end if
   pixcg(ix) = frcm1 + betafc*pixcg(ix)
   pixcg(iy) = frcm2 + betafc*pixcg(iy)
   pixcg(iz) = frcm3 + betafc*pixcg(iz)
   fnorm(i) = fnormc
 else
   fnorm(i)  = 0.d0
   pixcg(ix) = 0.d0
   pixcg(iy) = 0.d0
   pixcg(iz) = 0.d0
end if
end do

do i = 1, 3*n
   x(i) = x(i) + pixcg(i)*dtmd
!         if( x(i) <  0.d0 ) x(i) = x(i) + 1.d0
!         if( x(i) >= 1.d0 ) x(i) = x(i) - 1.d0
end do


return
end




subroutine updtqn( nfile, myid, nodes,  &
& x, is, n, ntype, frc, watom, h, dtmd, pixcg, fnorm, acon,  &
& lfixion,  &
& deltax, deltaf, hessian, heigenvec, heigenval, scalex, scaleg,  &
& lqninitial, dist_max, lclear, lhessian, qnstabi, gammamin )
!-----------------------------------------------------------------------
!     Coordinate update, x(t+Dt) by 
!
!    Quasi-Newton methods using approximated Hessian with BFGS formula
!
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: n
integer :: ntype
real*8,  dimension(3*n) :: x
integer, dimension(n) :: is
real*8,  dimension(3*n) :: frc
real*8,  dimension(ntype) :: watom
real*8,  dimension(3,3) :: h
real*8  :: dtmd
real*8,  dimension(3*n) :: pixcg
real*8,  dimension(n)   :: fnorm
real*8,  dimension(ntype) :: acon
logical, dimension(ntype) :: lfixion
real*8,  dimension(3*n)     :: deltax, deltaf
real*8,  dimension(3*n,3*n) :: hessian, heigenvec
real*8,  dimension(3*n)     :: heigenval
real*8,  dimension(3*n)     :: scalex, scaleg
logical :: lqninitial
real*8  :: dist_max
logical :: lclear, lhessian
real*8  :: qnstabi, gammamin

!-----declare local variables
integer :: i, j, it, ix, iy, iz, k, l
real*8  :: frcm1, frcm2, frcm3, fmax, fnormc, dispmax
real*8  :: denm1, denm2, gamma
character(5) :: checkc


negativehdo: do
if( lqninitial ) then
    lqninitial = .false.
    !-----for first subroutine call

    !-----save force
    fmax = 0.d0
    do i = 1, n
       it = is(i)
       ix = 3*i-2
       iy = 3*i-1
       iz = 3*i
    if( .not.lfixion(it) ) then
       frcm1 = frc(ix) / acon(it)
       frcm2 = frc(iy) / acon(it)
       frcm3 = frc(iz) / acon(it)
       pixcg(ix) = frcm1
       pixcg(iy) = frcm2
       pixcg(iz) = frcm3
       fnormc = frcm1*frcm1 + frcm2*frcm2 + frcm3*frcm3
       fmax = max( fmax, fnormc )
      else
       pixcg(ix) = 0.d0
       pixcg(iy) = 0.d0
       pixcg(iz) = 0.d0
    end if
    end do
    fmax  = sqrt(fmax)

    !-----initial guess for displacements
!    dispmax = 0.01d0/max(h(1,1),h(2,2),h(3,3)) / fmax
!    deltax = dispmax*pixcg
    deltax = dtmd*pixcg

    !-----initial guess for Hessian
    hessian = 0.d0
    do i = 1, 3*n
       hessian(i,i) = 1.d0
    end do

    !-----coordinate update
    x = x + deltax

    !-----get real variables
    call get_realv( deltax, h, n )

    return
end if


if( lclear ) then
    !-----initial guess for Hessian
    hessian = 0.d0
    do i = 1, 3*n
       hessian(i,i) = 1.d0
    end do
end if


!-----force difference
do i = 1, n
   it = is(i)
   ix = 3*i-2
   iy = 3*i-1
   iz = 3*i
if( .not.lfixion(it) ) then
   frcm1 = frc(ix) / acon(it)
   frcm2 = frc(iy) / acon(it)
   frcm3 = frc(iz) / acon(it)
   deltaf(ix) = frcm1 - pixcg(ix)
   deltaf(iy) = frcm2 - pixcg(iy)
   deltaf(iz) = frcm3 - pixcg(iz)
   pixcg(ix) = frcm1
   pixcg(iy) = frcm2
   pixcg(iz) = frcm3
  else
   deltaf(ix) = 0.d0
   deltaf(iy) = 0.d0
   deltaf(iz) = 0.d0
   pixcg(ix) = 0.d0
   pixcg(iy) = 0.d0
   pixcg(iz) = 0.d0
end if
end do


!-----get real variables
call get_realv( deltaf, h, n )

!-----Hessian update: BFGS formula -------------------------------
denm1  = -1.d0/dot_product(deltaf,deltax)
scalex = matmul(hessian,deltax)
denm2  = -1.d0/dot_product(deltax,scalex)

do i = 1, 3*n
do j = 1, i
   hessian(j,i) = hessian(j,i)  &
& + denm1*deltaf(j)*deltaf(i) + denm2*scalex(j)*scalex(i)
end do
end do
do i = 1, 3*n
do j = 1, i-1
   hessian(i,j) = hessian(j,i)
end do
end do
!-----Hessian update: BFGS formula -------------------------------


!-----Unitary matrix 
heigenvec = hessian
checkc = 'OK   '
call RsEIGQR( heigenvec, heigenval, 3*n, 3*n,  &
& 1.0d-13, 1, 2, checkc, scalex, scaleg, deltax )

!---check negative eigenvalues
lhessian = heigenval(1) < 0.d0 .or. abs(heigenval(1)-1.d0) < 1.d-10
if( .not.lhessian ) exit negativehdo
    if( myid == 0 ) then
        write(nfile(1),*) 'heigenval =', heigenval
        write(nfile(2),*) 'heigenval =', heigenval
!        return
        write(nfile(1),*) 'negative or invalid Hessian, still continue ...'
        write(nfile(2),*) 'negative or invalid Hessian, still continue ...'
    end if
    lhessian = .false.
    lqninitial = .true.
!----check
!      do i = 1, 3*n
!         write(*,*) i, heigenval(i)
!      end do
!      do i = 1, 3*n
!      do j = 1, 3*n
!         denm1 = 0.d0
!         do k = 1, 3*n
!         do l = 1, 3*n
!            denm1 = denm1 + heigenvec(k,i)*hessian(k,l)*heigenvec(l,j)
!         end do
!         end do
!         write(*,*) i,j,denm1,heigenval(i)
!      end do
!      end do
end do negativehdo
if( myid == 0 ) then
    write(nfile(1),*) 'The status of eigenvalue problem: ', checkc
    write(nfile(2),*) 'The status of eigenvalue problem: ', checkc
    write(nfile(1),*) 'The lowest three eigenvalues of Hessian'
    write(nfile(2),*) 'The lowest three eigenvalues of Hessian'
    do i = 1, 3  !3*n
       write(nfile(1),*) i, heigenval(i)
       write(nfile(2),*) i, heigenval(i)
    end do
end if


!-----scaled gradients
deltaf = pixcg
!-----get real variables
call get_realv( deltaf, h, n )
do i = 1, 3*n
   denm1 = 0.d0
   do j = 1, 3*n
      denm1 = denm1 + heigenvec(j,i)*deltaf(j)
   end do
   scaleg(i) = denm1
end do

!-----Rational function
gamma = 0.5d0*( sqrt(heigenval(1)*heigenval(1)+4.d0*scaleg(1)*scaleg(1)) )
if( gamma >= gammamin ) then
    if( myid == 0 ) then
        write(nfile(1),*) ' gamma =', gamma
        write(nfile(2),*) ' gamma =', gamma
    end if
else
    if( myid == 0 ) then
        write(nfile(1),*) ' gamma =', gamma, ' ->', gammamin
        write(nfile(2),*) ' gamma =', gamma, ' ->', gammamin
    end if
    gamma = gammamin
end if
do i = 1, 3*n
   heigenval(i) = heigenval(i) + gamma
end do


!-----Newton-Raphson step
scalex = scaleg/heigenval

if( myid == 0 ) then
    write(nfile(1),*) 'The Newton-Raphson step for the lowest three eigenvalues'
    write(nfile(2),*) 'The Newton-Raphson step for the lowest three eigenvalues'
end if
do i = 1, 3*n
   if( i <= 3 ) then
       if( myid == 0 ) then
           write(nfile(1),*) i, scalex(i)
           write(nfile(2),*) i, scalex(i)
       end if
   end if
   if( heigenval(i) <= qnstabi ) then
       if( myid == 0 ) then
           write(nfile(1),*) 'Since heigenval <= qnstabi, Newton-Raphson step is cancelled for ', i
           write(nfile(2),*) 'Since heigenval <= qnstabi, Newton-Raphson step is cancelled for ', i
       end if
       scalex(i) = 0.d0
   end if
end do

!-----displacement vector
deltax = matmul(heigenvec,scalex)

!-----check displacement distance
call check_dist( deltax, dist_max, n )

!-----get scaled variables
scalex = deltax
call get_scalev( scalex, h, n )


!-----coordinate update
x = x + scalex


return
end




subroutine outqnrfo( n, iunit, ierror, iogpsz )
!-----------------------------------------------------------------------
!    save variables
!-----------------------------------------------------------------------
use remd_atom
implicit none
integer :: n, iunit, ierror, iogpsz

!-----declare local variables
integer :: i3, istat


write(iunit,iostat=istat) (hessian(i3),i3=1,3*n*3*n)
ierror = max( ierror, abs(istat) )

write(iunit,iostat=istat) (deltax(i3),i3=1,3*n)
ierror = max( ierror, abs(istat) )


return
end




subroutine loadqnrfo( n, iunit, ierror, lload, iogpsz )
!-----------------------------------------------------------------------
!    load variables
!-----------------------------------------------------------------------
use remd_atom
implicit none
integer :: n, iunit, ierror, iogpsz
logical :: lload

!-----declare local variables
integer :: i3, istat
real*8  :: dummy


if( lload ) then
    read(iunit,iostat=istat) (hessian(i3),i3=1,3*n*3*n)
    ierror = max( ierror, abs(istat) )

    read(iunit,iostat=istat) (deltax(i3),i3=1,3*n)
    ierror = max( ierror, abs(istat) )
else
    read(iunit,iostat=istat) (dummy,i3=1,3*n*3*n)
    ierror = max( ierror, abs(istat) )

    read(iunit,iostat=istat) (dummy,i3=1,3*n)
    ierror = max( ierror, abs(istat) )
end if

return
end




subroutine updtcg_cell( nfile, myid, nodes,  &
& h, volume, pint, pext, dtcellcg, cellfnorm, irstrct, irstrct_sub, lcell_rstrct )
!----------------------------------------------------------------------c
!  Supercell optimization by conjugate-gradient method
!----------------------------------------------------------------------c
implicit none
integer :: nfile(*), myid, nodes
real*8,  dimension(3,3,0:1) :: h
real*8  :: volume
real*8,  dimension(3,3) :: pint, pext
real*8  :: dtcellcg, cellfnorm
integer :: irstrct, irstrct_sub
logical :: lcell_rstrct(3)

!-----declare local variables
real*8,  dimension(3,3) :: vtemp, hi, cellfrc, pintext
integer :: i, j, k
real*8  :: fnormc, betafc


!-----restriction for stress tensor
pintext = pint - pext
!call restrictstr_for_opt( irstrct, pintext )


!----- transpose matrix of vtemp = inverse of h
CALL RCIPRL( h(1,1,0), vtemp, volume )
call mattrp( vtemp, hi )

do j = 1, 3
do i = 1, 3
   cellfrc(i,j) = 0.d0
   do k = 1, 3
      cellfrc(i,j) = cellfrc(i,j) + ( pintext(i,k) )*hi(j,k)
   end do
end do
end do
cellfrc = cellfrc * volume

fnormc = 0.d0
do j = 1, 3
do i = 1, 3
   fnormc = fnormc + cellfrc(i,j)*cellfrc(i,j)
end do
end do
if( fnormc < cellfnorm ) then
    betafc = fnormc/cellfnorm
  else
    betafc = 0.d0
end if
do j = 1, 3
do i = 1, 3
   h(i,j,1) = cellfrc(i,j) + betafc*h(i,j,1)
end do
end do
cellfnorm = fnormc

!-----by symmetry operations
!if( irstrct >= 10 ) call restrictMDcell_sym( h(1,1,1), h(1,1,0) )

do j = 1, 3
do i = 1, 3
   h(i,j,0) = h(i,j,0) + dtcellcg*h(i,j,1)
end do
end do

!-----restriction for MD cell
!call restrictMDcell_for_opt( irstrct, h )
!call restrictMDcell( irstrct, irstrct_sub, h, lcell_rstrct )


return
end subroutine




subroutine updtqn_cell( nfile, myid, nodes,  &
& h, dh, volume, pint, pext, dtcellcg, cellfnorm, hessian, savefrc,  &
& lqninitial, irstrct, irstrct_sub, lcell_rstrct,  &
& lclearcellh, lhessian, qnstabi, gammamin )
!----------------------------------------------------------------------c
!  Supercell optimization by
!    Quasi-Newton methods using approximated Hessian with BFGS formula
!----------------------------------------------------------------------c
implicit none
integer :: nfile(*), myid, nodes
real*8,  dimension(9) :: h, dh
real*8  :: volume
real*8,  dimension(3,3) :: pint, pext
real*8  :: dtcellcg, cellfnorm
real*8,  dimension(9,9) :: hessian
real*8,  dimension(9) :: savefrc
logical :: lqninitial
integer :: irstrct, irstrct_sub
logical :: lcell_rstrct(3)
logical :: lclearcellh, lhessian
real*8  :: qnstabi, gammamin

!-----declare local variables
real*8,  dimension(3,3) :: vtemp, hi, pintext
real*8,  dimension(9) :: cellfrc, deltaf, scalex, scaleg, heigenval
real*8,  dimension(9,9) :: heigenvec
real*8  :: fnormc, denm1, denm2, gamma
integer :: i, j, k, i1, i2, j1, j2
character(5) :: checkc
real*8  :: dttry = 1.d+10, dhnorm
save dttry


!-----restriction for stress tensor
pintext(1:3,1:3) = pint(1:3,1:3) - pext(1:3,1:3)
!call restrictstr_for_opt( irstrct, pintext )


!----- transpose matrix of vtemp = inverse of h
CALL RCIPRL( h, vtemp, volume )
call mattrp( vtemp, hi )

negativehdo: do
do j = 1, 3
do i = 1, 3
   cellfrc(i+3*(j-1)) = 0.d0
   do k = 1, 3
      cellfrc(i+3*(j-1)) = cellfrc(i+3*(j-1)) + ( pintext(i,k) )*hi(j,k)
   end do
end do
end do
cellfrc = cellfrc * volume

fnormc = 0.d0
do i = 1, 9
   fnormc = fnormc + cellfrc(i)*cellfrc(i)
end do

if( lqninitial ) then
    lqninitial = .false.
    !-----for first subroutine call

    !-----save force
    savefrc(1:9) = cellfrc(1:9)
    cellfnorm = fnormc

    dttry = min( dttry, dtcellcg )
    !-----initial guess for displacements
    dh(1:9) = dttry*cellfrc(1:9)

    !-----by symmetry operations
!    if( irstrct >= 10 ) call restrictMDcell_sym( dh, h )

    !-----coordinate update
    h(1:9) = h(1:9) + dh(1:9)

    !-----restriction for MD cell
!    call restrictMDcell_for_opt( irstrct, h )
!    call restrictMDcell( irstrct, irstrct_sub, h, lcell_rstrct )

    !-----initial guess for Hessian
    hessian = 0.d0
    do i = 1, 9
       hessian(i,i) = 1.d0
    end do

    return
end if

!if( fnormc > cellfnorm .or. lclearcellh ) then
if( lclearcellh ) then
    !-----initial guess for Hessian
    hessian = 0.d0
    do i = 1, 9
       hessian(i,i) = 1.d0
    end do
end if
cellfnorm = fnormc


!-----force difference
deltaf(1:9) = cellfrc(1:9) - savefrc(1:9)

!-----save force
savefrc(1:9) = cellfrc(1:9)

!-----Hessian update: BFGS formula -------------------------------
denm1  = -1.d0/dot_product(deltaf,dh)
scalex = matmul(hessian,dh)
denm2  = -1.d0/dot_product(dh,scalex)

do i = 1, 9
do j = 1, 9
   hessian(j,i) = hessian(j,i)  &
& + denm1*deltaf(j)*deltaf(i) + denm2*scalex(j)*scalex(i)
end do
end do
!-----Hessian update: BFGS formula -------------------------------


!-----Unitary matrix 
heigenvec = hessian
checkc = 'OK   '
call RsEIGQR( heigenvec, heigenval, 9, 9,  &
& 1.0d-13, 1, 2, checkc, scalex, scaleg, dh )

!---check negative eigenvalues
lhessian = heigenval(1) < 0.d0 .or. abs(heigenval(1)-1.d0) < 1.d-10
if( .not.lhessian ) exit negativehdo
    if( myid == 0 ) then
        write(nfile(1),*) 'updtqn_cell: heigenval =', heigenval
        write(nfile(2),*) 'updtqn_cell: heigenval =', heigenval
!        return
        write(nfile(1),*) 'negative or invalid Hessian, still continue ...'
        write(nfile(2),*) 'negative or invalid Hessian, still continue ...'
    end if
    lhessian = .false.
    lqninitial = .true.
end do negativehdo
if( myid == 0 ) then
    write(nfile(1),*) 'Quasi-Newton minimization for supercell'
    write(nfile(2),*) 'Quasi-Newton minimization for supercell'
    write(nfile(1),*) '  The status of eigenvalue problem: ', checkc
    write(nfile(2),*) '  The status of eigenvalue problem: ', checkc
    write(nfile(1),*) '  The lowest three eigenvalues of Hessian'
    write(nfile(2),*) '  The lowest three eigenvalues of Hessian'
    do i = 1, 3  !9
       write(nfile(1),*) i, heigenval(i)
       write(nfile(2),*) i, heigenval(i)
    end do
end if

!-----gradients
deltaf(1:9) = savefrc(1:9)

do i = 1, 9
   denm1 = 0.d0
   do j = 1, 9
      denm1 = denm1 + heigenvec(j,i)*deltaf(j)
   end do
   scaleg(i) = denm1
end do

!-----Rational function
gamma = 0.5d0*( sqrt(heigenval(1)*heigenval(1)+4.d0*scaleg(1)*scaleg(1)) )
if( gamma >= gammamin ) then
    if( myid == 0 ) then
        write(nfile(1),*) '   gamma =', gamma
        write(nfile(2),*) '   gamma =', gamma
    end if
else
    if( myid == 0 ) then
        write(nfile(1),*) '   gamma =', gamma, ' ->', gammamin
        write(nfile(2),*) '   gamma =', gamma, ' ->', gammamin
    end if
    gamma = gammamin
end if
do i = 1, 9
   heigenval(i) = heigenval(i) + gamma
end do


!-----Newton-Raphson step
scalex(1:9) = scaleg(1:9)/heigenval(1:9)

if( myid == 0 ) then
    write(nfile(1),*) '  The Newton-Raphson step for the lowest three eigenvalues'
    write(nfile(2),*) '  The Newton-Raphson step for the lowest three eigenvalues'
end if
do i = 1, 9
   if( i <= 3 ) then
       if( myid == 0 ) then
           write(nfile(1),*) i, scalex(i)
           write(nfile(2),*) i, scalex(i)
       end if
   end if
   if( heigenval(i) <= qnstabi ) then
       if( myid == 0 ) then
           write(nfile(1),*) 'Since heigenval <= qnstabi, Newton-Raphson step is cancelled for ', i
           write(nfile(2),*) 'Since heigenval <= qnstabi, Newton-Raphson step is cancelled for ', i
       end if
       scalex(i) = 0.d0
   end if
end do

!-----displacement vector
dh = matmul(heigenvec,scalex)

!-----by symmetry operations
!if( irstrct >= 10 ) call restrictMDcell_sym( dh, h )

!-----coordinate update
h(1:9) = h(1:9) + dh(1:9)


!-----restriction for MD cell
!call restrictMDcell_for_opt( irstrct, h )
!call restrictMDcell( irstrct, irstrct_sub, h, lcell_rstrct )


dhnorm = 0.d0
do i = 1, 9
   dhnorm = dhnorm + dh(i)*dh(i)
end do
dttry = min( dttry, sqrt(fnormc/dhnorm) )


return
end subroutine




subroutine outqnrfo_cell( iunit, ierror, iogpsz )
!-----------------------------------------------------------------------
!    save variables
!-----------------------------------------------------------------------
use remd_variables
implicit none
integer :: iunit, ierror, iogpsz

!-----declare local variables
integer :: istat


write(iunit,iostat=istat) cellfnorm
ierror = max( ierror, abs(istat) )

write(iunit,iostat=istat) hessian_cell(1:9,1:9)
ierror = max( ierror, abs(istat) )

write(iunit,iostat=istat) savefrc_cell(1:3,1:3)
ierror = max( ierror, abs(istat) )


return
end




subroutine loadqnrfo_cell( iunit, ierror, iogpsz )
!-----------------------------------------------------------------------
!    load variables
!-----------------------------------------------------------------------
use remd_variables
implicit none
integer :: iunit, ierror, iogpsz

!-----declare local variables
integer :: i, j
integer :: istat
real*8, parameter :: small = 1.d-30


read(iunit,iostat=istat) cellfnorm
if( istat == 0 )  &
& read(iunit,iostat=istat) hessian_cell(1:9,1:9)
if( istat == 0 )  &
& read(iunit,iostat=istat) savefrc_cell(1:3,1:3)

if( istat == 0 ) then
    outer: do j = 1, 3
       do i = 1, 3
          if( abs(savefrc_cell(i,j)) < small ) cycle
          lqninitial_cell = .false.
          exit outer
       end do
    end do outer
end if


return
end




subroutine check_dist( a, d, n )
!-----------------------------------------------------------------------
!    get real variables
!-----------------------------------------------------------------------
implicit none
integer :: n
real*8,  dimension(3,n) :: a
real*8  :: d

!-----declare local variables
integer :: i
real*8  :: t, d2


d2 = d*d
do i = 1, n
   t = a(1,i)*a(1,i) + a(2,i)*a(2,i) + a(3,i)*a(3,i)
   if( t > d2 ) then
       t = d/sqrt(t)
       a(:,i) = a(:,i)*t
   end if
end do


return
end




subroutine get_realv( a, h, n )
!-----------------------------------------------------------------------
!    get real variables
!-----------------------------------------------------------------------
implicit none
integer :: n
real*8,  dimension(3,n) :: a
real*8,  dimension(3,3) :: h

!-----declare local variables
integer :: i
real*8,  dimension(3) :: t


do i = 1, n
   t = h(:,1)*a(1,i) + h(:,2)*a(2,i) + h(:,3)*a(3,i)
   a(:,i) = t
end do


return
end




subroutine get_scalev( a, h, n )
!-----------------------------------------------------------------------
!    get scaled variables
!-----------------------------------------------------------------------
implicit none
integer :: n
real*8,  dimension(3,n) :: a
real*8,  dimension(3,3) :: h

!-----declare local variables
integer :: i
real*8,  dimension(3) :: t
real*8,  dimension(3,3) :: hi
real*8  :: vtemp


!-----transpose matrix of hi = inverse of h
CALL RCIPRL( h, hi, vtemp )
do i = 1, n
   t = hi(1,:)*a(1,i) + hi(2,:)*a(2,i) + hi(3,:)*a(3,i)
   a(:,i) = t
end do


return
end




subroutine get_ekinp( nfile, myid, nodes,  &
& pit1, v, is, n, ntype, h, fack, volume,  &
& lcatomic, wstrsk, nmdmax__, ntot )
!-----------------------------------------------------------------------
!    Kinetic energy contribution in MD
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
real*8,  dimension(3,3) :: pit1
integer :: n
integer :: ntype
real*8,  dimension(3,n) :: v
integer, dimension(n) :: is
real*8,  dimension(3,3,0:1) :: h
real*8,  dimension(ntype)   :: fack
real*8  :: volume
logical :: lcatomic
integer :: nmdmax__
real*8  :: wstrsk(6,nmdmax__)
integer, dimension(0:ntype) :: ntot

!-----declare local variables
integer :: i, ic
real*8  :: svx, svy, svz, vxi, vyi, vzi
real*8,  dimension(6) :: pitv, dbuf


pitv(1:6) = 0.d0
wstrsk(:,:) = 0.d0
do i = 1, n
  ic=is(i)
  svx=v(1,i)
  svy=v(2,i)
  svz=v(3,i)
  vxi=h(1,1,0)*svx+h(1,2,0)*svy+h(1,3,0)*svz
  vyi=h(2,1,0)*svx+h(2,2,0)*svy+h(2,3,0)*svz
  vzi=h(3,1,0)*svx+h(3,2,0)*svy+h(3,3,0)*svz
  dbuf(1) = 2d0*fack(ic)*vxi*vxi
  dbuf(2) = 2d0*fack(ic)*vyi*vyi
  dbuf(3) = 2d0*fack(ic)*vzi*vzi
  dbuf(4) = 2d0*fack(ic)*vyi*vzi
  dbuf(5) = 2d0*fack(ic)*vzi*vxi
  dbuf(6) = 2d0*fack(ic)*vxi*vyi
  pitv(1:6) = pitv(1:6) + dbuf(1:6)
  if( lcatomic ) wstrsk(1:6,i) = wstrsk(1:6,i) + dbuf(1:6)
end do
!-----global sum
call gdsum(pitv,6,dbuf)
do i = 1, 6
   pitv(i) = pitv(i)/volume
end do
pit1(1,1) = pit1(1,1) + pitv(1)
pit1(2,2) = pit1(2,2) + pitv(2)
pit1(3,3) = pit1(3,3) + pitv(3)
pit1(2,3) = pit1(2,3) + pitv(4)
pit1(3,1) = pit1(3,1) + pitv(5)
pit1(1,2) = pit1(1,2) + pitv(6)
pit1(3,2) = pit1(2,3)
pit1(1,3) = pit1(3,1)
pit1(2,1) = pit1(1,2)

if( lcatomic ) wstrsk(1:6,1:n) = wstrsk(1:6,1:n)*ntot(0)/volume


return
end




subroutine get_force_null( nfile, myid, nodes,  &
& epot, x, is, n, ntype, a, h, acon, pintlr, pintsr )
!-----------------------------------------------------------------------
!     classical energy and forces
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
real*8  :: epot
integer :: n
integer :: ntype
real*8,  dimension(3*n,0:1) :: x
integer, dimension(n)       :: is
real*8,  dimension(3*n)     :: a
real*8,  dimension(3,3,0:1) :: h
real*8,  dimension(ntype)   :: acon
real*8,  dimension(3,3) :: pintlr, pintsr

!-----declare local variables
integer :: i, ix


epot = 0.d0

!-----calculate physical forces, A
do i = 1, 3*n
   a(i) = 0.d0
end do

do i = 1, 3
do ix = 1, 3
   pintlr(ix,i) = 0.d0
   pintsr(ix,i) = 0.d0
end do
end do


!      !-----Calculate normalized accelerations, A
!      call normalized_a( nfile, myid, nodes, n, ntype, a, is, h, acon )


return
end subroutine




subroutine normalized_a( nfile, myid, nodes,  &
& n, ntype, a, is, h, acon )
!-----------------------------------------------------------------------
!     Calculate normalized accelerations, A
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: n
integer :: ntype
real*8,  dimension(3*n)     :: a
integer, dimension(n)       :: is
real*8,  dimension(3,3,0:1) :: h
real*8,  dimension(ntype)   :: acon

!-----declare local variables
integer :: i, i3, ix, iy, iz
real*8,  dimension(3,3) :: hi
real*8  :: vtemp, arx, ary, arz, ax, ay, az


!-----Inverse MD-cell transform
!-----transpose matrix of hi = inverse of h
CALL RCIPRL( h(1,1,0), hi, vtemp )
do i = 1, n
  iz=3*i
  iy=iz-1
  ix=iy-1
  !-----Physical force, AR
  arx=a(ix)
  ary=a(iy)
  arz=a(iz)
  !-----Reduced force, A=HI*AR
  ax=hi(1,1)*arx+hi(2,1)*ary+hi(3,1)*arz
  ay=hi(1,2)*arx+hi(2,2)*ary+hi(3,2)*arz
  az=hi(1,3)*arx+hi(2,3)*ary+hi(3,3)*arz
  !-----Multiply the normalization factor
  a(ix) = ax*acon(is(i))
  a(iy) = ay*acon(is(i))
  a(iz) = az*acon(is(i))
end do


return
end subroutine




subroutine set_pratm( nfile, myid, nodes, ratm, pratm, natom )
!-----------------------------------------------------------------------
!     set pratm
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: natom
real*8,  dimension(3*natom) :: ratm
real*8,  dimension(3*natom) :: pratm

!-----declare local variables
integer :: i

do i = 1, 3*natom
   pratm(i) = ratm(i)
end do

return
end subroutine




subroutine updt_atom_disp( nfile, myid, nodes,                     &
& lclust, ratm, pratm, prevr, h,  &
& ntype, natom, nhk1, nhk2 )
!-----------------------------------------------------------------------
!     update atomic displacements in MD
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
logical :: lclust
integer :: ntype
integer :: natom
real*8,  dimension(3,natom) :: ratm
real*8,  dimension(3,natom) :: pratm
real*8,  dimension(3,natom,3) :: prevr
real*8,  dimension(3,3) :: h
integer, dimension(ntype) :: nhk1, nhk2

!-----declare local variables
integer :: i, j
real*8  :: q1, q2, q3


do i = 1, nhk2(ntype)
do j = 1, 3
   prevr(j,i,1) = prevr(j,i,2)
   prevr(j,i,2) = prevr(j,i,3)
   prevr(j,i,3) = ratm(j,i) - pratm(j,i)
end do
end do

if( .not.lclust ) then
    do j = 1, 3
    do i = 1, nhk2(ntype)
       if( abs(prevr(j,i,3)).gt.0.5d0 )  &
&          prevr(j,i,3) = prevr(j,i,3) - sign(1.d0,prevr(j,i,3))
    end do
    end do
!          do i = 1, nhk2(ntype)
!             q1 = prevr(1,i,3)
!             q2 = prevr(2,i,3)
!             q3 = prevr(3,i,3)
!             prevr(1,i,3) = h(1,1)*q1 + h(1,2)*q2 + h(1,3)*q3
!             prevr(2,i,3) = h(2,1)*q1 + h(2,2)*q2 + h(2,3)*q3
!             prevr(3,i,3) = h(3,1)*q1 + h(3,2)*q2 + h(3,3)*q3
!          end do
end if


return
end subroutine




subroutine max_force( nfile, myid, nodes,  &
& n, ntype, a, is, lfixion, frcmax )
!-----------------------------------------------------------------------
!     get maximum norm of atomic forces
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: n
integer :: ntype
real*8,  dimension(3*n) :: a
integer, dimension(n) :: is
logical, dimension(ntype) :: lfixion
real*8 :: frcmax

!-----declare local variables
integer :: i, ix, iy, iz, it


frcmax = 0.d0
do i = 1, n
   it = is(i)
if( .not.lfixion(it) ) then
   iz=3*i
   iy=iz-1
   ix=iy-1
   frcmax = max( frcmax, a(ix)*a(ix) + a(iy)*a(iy) + a(iz)*a(iz) )
end if
end do
frcmax = sqrt(frcmax)
call gdmax(frcmax)


return
end subroutine




subroutine setacc(  &
& prevr, facc1, facc2, facc3, ntype, natom, nhk1, nhk2 )
!-----------------------------------------------------------------------
!     get coefficients for extraporation
!-----------------------------------------------------------------------
implicit real*8 ( a-h, o-z )
dimension   prevr(3,natom,3)
dimension   nhk1(ntype), nhk2(ntype)
dimension   amt(3,3)


do n1 = 1, 3
do n2 = 1, n1
   amt(n1,n2) = 0.d0
   do i  = 1, nhk2(ntype)
   do j = 1, 3
      amt(n1,n2) = amt(n1,n2) + prevr(j,i,n1)*prevr(j,i,n2)
   end do
   end do
   amt(n2,n1) = amt(n1,n2)
end do
end do
abrmmx = max( amt(1,1), amt(2,2), amt(3,3) )
if( abs(abrmmx).lt.1.d-14 ) then
    facc1 =  2.d0
    facc2 = -1.d0
    facc3 =  1.d0
    return
end if
do n1 = 1, 3
do n2 = 1, 3
   amt(n1,n2) = amt(n1,n2)/abrmmx
end do
end do


!--- for first order extraporation
if( amt(2,2).ge.1.d-05 ) then
    facc3 = amt(3,2)/amt(2,2)
else
    facc3 = 1.d0
end if

!--- for second order extraporation
bbb = amt(2,2)*amt(1,1) - amt(1,2)*amt(2,1)
if( abs(bbb).gt.1.d-05 ) then
    facc1 = (   amt(1,1)*amt(3,2) - amt(1,2)*amt(3,1) )/bbb
    facc2 = ( - amt(1,2)*amt(3,2) + amt(2,2)*amt(3,1) )/bbb
  else
    facc1 =  2.d0
    facc2 = -1.d0
end if


return
end




subroutine msst_setmatrix( nfile, myid, nodes,  &
& nshockv, h, alpmatrix, betmatrix, trsmatrix, ierror )
!-----------------------------------------------------------------------
!     set alpmatrix and betmatrix for MSST
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: nshockv(1:3)
real*8  :: h(3,3), alpmatrix(3,3), betmatrix(3,3), trsmatrix(3,3)
integer :: ierror

!-----declare local variables
integer :: i, j, k, ix, ip(3), ierr
real*8  :: amat(3,3), gmat(3,3), hmat(3,3), sum, det
real*8  :: newlmat(3,3), lenl3, cos23, sin23
real*8  :: v(3), hi(3,3)


!--- v = h^t h n1
hmat(:,:) = 0.d0
do i = 1, 3
do j = 1, 3
   do k = 1, 3
      hmat(j,i) = hmat(j,i) + h(k,j)*h(k,i)
   end do
end do
end do
v(:) = 0.d0
do i = 1, 3
   do j = 1, 3
      v(i) = v(i) + hmat(i,j)*nshockv(j)
   end do
end do


amat(:,:) = 0.d0
if( nshockv(1) /= 0 .and. nshockv(2) == 0 .and. nshockv(3) == 0 ) then
    amat(1,1) = 1.d0
    amat(2,2) = 1.d0
    amat(3,3) = 1.d0
else if( nshockv(1) == 0 .and. nshockv(2) /= 0 .and. nshockv(3) == 0 ) then
    amat(2,1) = 1.d0
    amat(3,2) = 1.d0
    amat(1,3) = 1.d0
else if( nshockv(1) == 0 .and. nshockv(2) == 0 .and. nshockv(3) /= 0 ) then
    amat(3,1) = 1.d0
    amat(1,2) = 1.d0
    amat(2,3) = 1.d0
else if( nshockv(3) == 0 ) then
    amat(1,1) = nshockv(1)
    amat(2,1) = nshockv(2)
    amat(3,1) = 0.d0
    amat(1,2) = -nshockv(2)
    amat(2,2) =  nshockv(1)
    amat(3,2) = 0.d0
    amat(1,3) = 0.d0
    amat(2,3) = 0.d0
    amat(3,3) = 1.d0
else if( nshockv(2) == 0 ) then
    amat(1,1) = nshockv(1)
    amat(2,1) = 0.d0
    amat(3,1) = nshockv(3)
    amat(1,2) = -nshockv(3)
    amat(2,2) = 0.d0
    amat(3,2) = nshockv(1)
    amat(1,3) = 0.d0
    amat(2,3) = 1.d0
    amat(3,3) = 0.d0
else if( nshockv(1) == 0 ) then
    amat(1,1) = 0.d0
    amat(2,1) = nshockv(2)
    amat(3,1) = nshockv(3)
    amat(1,2) = 0.d0
    amat(2,2) = -nshockv(3)
    amat(3,2) =  nshockv(2)
    amat(1,3) = 1.d0
    amat(2,3) = 0.d0
    amat(3,3) = 0.d0
else
    amat(1,1) = nshockv(1)
    amat(2,1) = nshockv(2)
    amat(3,1) = nshockv(3)
    amat(1,2) = -nshockv(2)
    amat(2,2) =  nshockv(1)
    amat(3,2) = 0.d0
    amat(1,3) = -nshockv(1)*nshockv(3)
    amat(2,3) = -nshockv(2)*nshockv(3)
    amat(3,3) = nshockv(1)*nshockv(1)+nshockv(2)*nshockv(2)
end if

!---preset vectors l1, l2, and l3
hmat(:,:) = 0.d0
do i = 1, 3
   do ix = 1, 3
      do j = 1, 3
         hmat(ix,i) = hmat(ix,i) + h(ix,j)*amat(j,i)
      end do
   end do
end do

!---normalization of l1
det = 0.d0
do k = 1, 3
   det = det + hmat(k,1)**2
end do
det = sqrt(det)
hmat(1:3,1) = hmat(1:3,1)/det

!---orthogonalization, if necessary
call matinv33( h, hi )
do j = 2, 3
   sum = 0.d0
   do ix = 1, 3
      sum = sum + hmat(ix,1)*hmat(ix,j)
   end do
   if( abs(sum) > 1.d-10 ) then
       do k = 1, 3
          hmat(k,j) = hmat(k,j) - sum*hmat(k,1)
       end do
       do i = 1, 3
          amat(i,j) = 0.d0
          do k = 1, 3
             amat(i,j) = amat(i,j) + hi(i,k)*hmat(k,j)
          end do
       end do
   end if
end do

!---check orthogonalization of cell vector
hmat(:,:) = 0.d0
do i = 1, 3
   do ix = 1, 3
      do j = 1, 3
         hmat(ix,i) = hmat(ix,i) + h(ix,j)*amat(j,i)
      end do
   end do
end do
!do i = 1, 3
do i = 1, 1
do j = i+1, 3
   sum = 0.d0
   do ix = 1, 3
      sum = sum + hmat(ix,i)*hmat(ix,j)
   end do
   if( abs(sum) > 1.d-10 ) then
       if( myid == 0 ) then
           do k = 1, 2
           write(nfile(k),*) 'MSST error : not orthogonalized cell vectors'
           write(nfile(k),*) 'original vectors'
           write(nfile(k),*) 'L1 =', h(1:3,1)
           write(nfile(k),*) 'L2 =', h(1:3,2)
           write(nfile(k),*) 'L3 =', h(1:3,3)
           write(nfile(k),*) 'extended vectors'
           write(nfile(k),*) 'L1 =', hmat(1:3,1)
           write(nfile(k),*) 'L2 =', hmat(1:3,2)
           write(nfile(k),*) 'L3 =', hmat(1:3,3)
           end do
       end if
       ierror = 100
       return
   end if
end do
end do


!---coordinates transformation so that hmat(1:3,1) // x-axis
newlmat(:,:) = 0.d0
newlmat(1,1) = sqrt(hmat(1,1)*hmat(1,1)+hmat(2,1)*hmat(2,1)+hmat(3,1)*hmat(3,1))
newlmat(2,2) = sqrt(hmat(1,2)*hmat(1,2)+hmat(2,2)*hmat(2,2)+hmat(3,2)*hmat(3,2))
lenl3        = sqrt(hmat(1,3)*hmat(1,3)+hmat(2,3)*hmat(2,3)+hmat(3,3)*hmat(3,3))
cos23 = (hmat(1,2)*hmat(1,3)+hmat(2,2)*hmat(2,3)+hmat(3,2)*hmat(3,3))/newlmat(2,2)/lenl3
sin23 = sqrt(1.d0-cos23*cos23)
newlmat(2,3) = lenl3*cos23
newlmat(3,3) = lenl3*sin23
call matinv33( hmat, gmat )

!---transform matrix
do i = 1, 3
do j = 1, 3
   sum = 0.d0
   do k = 1, 3
      sum = sum + newlmat(i,k)*gmat(k,j)
   end do
   trsmatrix(i,j) = sum
end do
end do

!---transform cell vectors
gmat(:,:) = h(:,:)
do i = 1, 3
do j = 1, 3
   sum = 0.d0
   do k = 1, 3
      sum = sum + trsmatrix(i,k)*gmat(k,j)
   end do
   h(i,j) = sum
end do
end do


!---set matrix
call matinv33( amat, gmat )
do i = 1, 3
do j = 1, 3
   alpmatrix(i,j) = amat(i,1)*gmat(1,j)
   betmatrix(i,j) = amat(i,2)*gmat(2,j) + amat(i,3)*gmat(3,j)
end do
end do


!---check determinant
amat(1:3,1:3) = alpmatrix(1:3,1:3) + betmatrix(1:3,1:3)

call DETMAT( det, amat, 3, IP, 3, ierr )

if( det < 0.d0 ) then
    if( myid == 0 ) then
        do k = 1, 2
           write(nfile(k),*) 'MSST error : negative determinant'
           write(nfile(k),*) 'original vectors'
           write(nfile(k),*) 'L1 =', h(1:3,1)
           write(nfile(k),*) 'L2 =', h(1:3,2)
           write(nfile(k),*) 'L3 =', h(1:3,3)
           write(nfile(k),*) 'extended vectors'
           write(nfile(k),*) 'L1 =', hmat(1:3,1)
           write(nfile(k),*) 'L2 =', hmat(1:3,2)
           write(nfile(k),*) 'L3 =', hmat(1:3,3)
        end do
    end if
    ierror = 101
    return
end if


return
end




subroutine out_hugoniot( nfile, myid, nodes, nstep, digit )
!-----------------------------------------------------------------------
!   Output Hugoniot relations
!-----------------------------------------------------------------------
use constants
use remd_param
implicit none
integer :: nfile(*), myid, nodes
integer :: nstep, digit

!-----declare local variables
real*8  :: partv, press, enegy, vratio, rho0
character(50) :: frmt


vratio = 1.d0 - volume/volum0
rho0   = totmass/volum0

partv = shockspeed*vratio
partv = partv / sqrt(welm*1.d-3/(hrdev*evdj))
!              [a.u.] -> [m/s]
!  T = L*sqrt(M/E) -> T/L = sqrt(M/E)

press = hpext + rho0*shockspeed*shockspeed*vratio

enegy = 0.5d0*(press + hpext)*vratio/rho0

press = press*prau
!              [a.u.] -> [GPa]

frmt = '(i7,f15.4,f17.8,es22.8)'
if( digit > 7 ) call get_frmt( frmt, digit )
!write(nfile(23),'(i7,f15.4,f17.8,es22.8)') nstep,  &
write(nfile(23),frmt) nstep,  &
partv, press, enegy

return
end
