



module input_variables
!-----------------------------------------------------------------------
! type declaration and initialization of variables in input.f90
!-----------------------------------------------------------------------
implicit none

logical :: losymop = .false.
save

end module




subroutine input( nfile, myid, nodes, lclust )
!-----------------------------------------------------------------------
!    read and check input data
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
logical :: lclust


!----- read input variables from an input file
call read_data( nfile, myid, nodes, lclust )


!----- read input variables for pseudopotential
call read_pp( nfile, myid, nodes )


!----- check input variables
call check_data( nfile, myid, nodes, lclust )


!----- write input variables
call write_data( nfile, myid, nodes, lclust )


return
end subroutine




subroutine check_data( nfile, myid, nodes, lclust )
!-----------------------------------------------------------------------
!    check input data
!-----------------------------------------------------------------------
use outfile
use param
use param_inispin
use param_atom
use constants
!use symmop
implicit none
integer :: nfile(*), myid, nodes
logical :: lclust

!-----declare local variables
integer :: ix, i, it, j, nzn, neband
real*8  :: x, y, z, q1, q2, q3
real*8  :: cmpos, small
real*8  :: ttmaxx
real*8  :: rdelr, rdelx, rdely, rdelz
real*8  :: rvol
real*8  :: qmax, qmin, qdis, volsuper, volMD
real*8, dimension(3,3) :: b
real*8, dimension(3)   :: sc
real*8  :: hinv(3,3)
integer :: ierror = 0


!-----the number of parallel nodes
nproc = npx*npy*npz


if( lnoncollinear ) then
    !--- for non-collinear magnetism, noband is doubled.
    neband = nband - noband
    noband = noband*2
    nband  = noband + neband

    !--- set the restriction for the direction of the magnetization
    if( latmref ) then
        !---if the restriction for each atom is on, disable the global restriction
        refmx = 0d0
        refmy = 0d0
        refmz = 0d0
!        lclearmagne_in_x = .false.
!        lclearmagne_in_y = .false.
!        lclearmagne_in_y = .false.
        if( inispin == 1 ) then
            !---set the global reference direction
            refmx = umx(1)
            refmy = umx(2)
            refmz = umx(3)
            latmref = .false.
        end if
        if( inispin < 1 .or. inispin > 4 ) then
            latmref = .false.
        end if
    end if

    if( inispin < 2 .or. inispin > 4 ) then
        lclearamag = .false.
    end if
end if

if( .not.lspin .and. .not.lnoncollinear ) then
    imxspin_trial = 0
end if


if( .not.lspin ) then
    lduplicate = .false.
end if

!----- check input variables & error trap ------------------------------
if( nd2v <= 0 ) then
    if(loutfile(1)) write(nfile(1),*) 'nd2v =', nd2v
    if(loutfile(2)) write(nfile(2),*) 'nd2v =', nd2v
    call fstop( nfile, myid, nodes, 'error : nd2v <= 0' )
end if

if( nband <= 0 ) then
    if(loutfile(1)) write(nfile(1),*) 'nband =', nband
    if(loutfile(2)) write(nfile(2),*) 'nband =', nband
    call fstop( nfile, myid, nodes, 'nband <= 0' )
end if

!--- not supported for nband < the number of processors
if( nband < nproc ) then
    nband = nproc
    if(loutfile(1)) write(nfile(1),*) '*** Since nband < nproc, ',  &
&                         'reset nband =', nband
    if(loutfile(2)) write(nfile(2),*) '*** Since nband < nproc, ',  &
&                         'reset nband =', nband
end if

!--- ecutsoft must be greater than ecut, and smaller than ecut*4.
!--- ecutdens must be greater than ecutsoft.
!--- if .not.lvand, ecutdens must be equal to ecutsoft.
if( .not.lvand ) ecutsoft = max( ecutdens, ecutsoft )
if( ecut > ecutsoft ) then
    ecutsoft = ecut
    if(loutfile(1)) write(nfile(1),*) '*** Since ecut > ecutsoft, ',  &
&                         'reset ecutsoft =', ecutsoft
    if(loutfile(2)) write(nfile(2),*) '*** Since ecut > ecutsoft, ',  &
&                         'reset ecutsoft =', ecutsoft
end if
if( ecut*4.d0 < ecutsoft ) then
    ecutsoft = ecut*4.d0
    if(loutfile(1)) write(nfile(1),*) '*** Since ecut*4 < ecutsoft, ',  &
&                         'reset ecutsoft =', ecutsoft
    if(loutfile(2)) write(nfile(2),*) '*** Since ecut*4 < ecutsoft, ',  &
&                         'reset ecutsoft =', ecutsoft
end if
if( ecutsoft > ecutdens ) then
    ecutdens = ecutsoft
    if(loutfile(1)) write(nfile(1),*) '*** Since ecutsoft > ecutdens, ',  &
&                         'reset ecutdens =', ecutdens
    if(loutfile(2)) write(nfile(2),*) '*** Since ecutsoft > ecutdens, ',  &
&                         'reset ecutdens =', ecutdens
end if
if( .not.lvand ) ecutdens = ecutsoft

!--- ecutorth must be smaller than or equal to ecut.
if( ecutorth > ecut .or. ecutorth < 1.d-14 ) then
    ecutorth = ecut
    if(loutfile(1)) write(nfile(1),*) '*** Since ecutorth > ecut or zero, ',  &
&                         'reset ecutorth =', ecutorth
    if(loutfile(2)) write(nfile(2),*) '*** Since ecutorth > ecut or zero, ',  &
&                         'reset ecutorth =', ecutorth
end if

!--- the default value of ecutc is ecut.
if( ecutc < 1.d-14 ) then
    ecutc = ecut
    if(loutfile(1)) write(nfile(1),*) '*** Since ecutc is zero, ',  &
&                         'reset ecutc =', ecutc
    if(loutfile(2)) write(nfile(2),*) '*** Since ecutc is zero, ',  &
&                         'reset ecutc =', ecutc
end if
if( ecutc > ecutdens ) then
    ecutc = ecutdens
    if(loutfile(1)) write(nfile(1),*) '*** Since ecutc > ecutdens, ',  &
&                         'reset ecutc =', ecutc
    if(loutfile(2)) write(nfile(2),*) '*** Since ecutc > ecutdens, ',  &
&                         'reset ecutc =', ecutc
end if


if( lfermi <= 0 ) lfermi = 1
if( lfermi >= 2 .and. tfermi < 1.d-15 ) then
    if(loutfile(1)) write(nfile(1),*) 'tfermi =', tfermi
    if(loutfile(2)) write(nfile(2),*) 'tfermi =', tfermi
    call fstop( nfile, myid, nodes, 'tfermi is too small' )
end if

lpcc = .false.
do it = 1, ntype
   if( rpcc(it) < -1.d-10 ) lpcci(it) = .false.
   lpcc = lpcc .or. lpcci(it)
end do

!---if lplusU = .false., no DFT+U
if( .not.lpaw   ) lplusU = .false.
if( .not.lplusU ) lplusU_at(1:ntype) = .false.
lplusU = .false.
do it = 1, ntype
   lplusU = lplusU .or. lplusU_at(it)
end do
if( lplusU ) then
    !----- frequencies [eV] -> [Ryd.]
    plusU_U(1:ntype) = plusU_U(1:ntype) / hrdev * 2.d0
    plusU_J(1:ntype) = plusU_J(1:ntype) / hrdev * 2.d0
end if

!---if lplusC = .false., no DFT+C
if( .not.lplusC ) lplusC_at(1:ntype) = .false.
lplusC = .false.
do it = 1, ntype
   lplusC = lplusC .or. lplusC_at(it)
end do

!---not supported yet for jhybrid >= 1 .and. jgga /= 2
if( jhybrid >= 1 .and. jgga_org /= 2 ) then
    call fstop( nfile, myid, nodes, 'not supported yet for jhybrid >= 1 .and. jgga /= 2' )
end if


if( multg <= 0 ) multg = 1
!-----for bulk calculations, the multigrid method is not used.
if( .not.lclust ) multg = 1


if( .not.lstart .and. ltddft ) then
    if( .not.ltddft_fssh .or. .not.lfssh_gsscf ) then
        ltddft = .false.
        ifmd   = 0
    end if
end if
if( .not.ltddft ) then
    ltddft_fssh = .false.
end if
if( ltddft .and. ltddft_fssh ) then
    ltddft = .false.
end if
if( ltddft_fssh .and. lfssh_gsscf ) then
    !tdbroad = 0.d0
    if( aslh_fssh < 1.d-05 ) aslh_fssh = aslh
    if( bslh_fssh < 1.d-05 ) bslh_fssh = bslh
end if
if( .not.ltddft_fssh ) then
    lfssh_gsscf = .false.
    lfssh_gsfrc = .false.
end if
if( ltddft_fssh ) then
    ldiagonal = .true.
    if( lrtddft ) then
        lfssh_gsscf = .true.
!        lexiton_dynamics = .true.
    end if
    if( .not.lfssh_gsscf ) then
        lfssh_gsfrc = .false.
    end if
    if( lfssh_gsfrc ) then
        ltddft_nscforce = .false.
    end if
end if
if( lrtddft ) then
    if( .not.ldiagonal .or. .not.lspin ) then
        lscissors_lrtddft = .false.
            lfission_rate = .false.
           lemission_rate = .false.
    end if
    if( .not.ltda ) lcic_tda = .false.
    !---for non-adiabatic FSSH simulations, CIC-TDA is not supported yet.
    if( ltddft_fssh ) lcic_tda = .false.
end if


if( ifmd < 0 ) ifmd = 0
if( ifmd > 6 .and. ifmd /= 10 ) ifmd = 0
!---if( lstart .or. ifmd == 0 ) nstep_ini = 0
if( ifmd == 0 ) nstep_ini = 0
if( liscale .and. iscstp <= 0 ) iscstp = 20

if( ifmd < 2 .and. ioptmze /= 1 ) laspc = .false.
if( .not.laspc ) lxlbomd = .false.

if( ichest < 0 ) ichest = 0
if( ihest  < 0 ) ihest  = 0
if( .not.laspc_chg .and. ichest > 2 ) ichest = 2
if( .not.laspc_wv  .and. ihest  > 2 ) ihest  = 2

call set_nprvmx( ichest, ihest, lxlbomd, kxlbomd )
!-----Set coefficients for the ASPC method
call set_coef_ASPC( nfile, myid, nodes, &
& ichest, ihest, laspc_chg, laspc_wv, nband, aspc_corr, lxlbomd, kxlbomd )

if( aslh_aspc < 1.d-05 ) aslh_aspc = aslh
if( bslh_aspc < 1.d-05 ) bslh_aspc = bslh


if( ibstt1 <= 0 ) ibstt1 = 1
if( ibstt2 <= 0 ) ibstt2 = nband

if( nskip_stress <= 0 ) nskip_stress = 1

if( nskpmulk <= 0 ) nskpmulk = 1

if( nskpeda /= nskpmulk ) nskpeda = nskpmulk

if( .not.lvshape ) pwscale = 1.d0
if( pwscale < 1.d0 ) pwscale = 1.2d0

!---noncollinear magnetism
if( .not.lnoncollinear ) then
    ncscale = 1
else
    ncscale = 2
end if

if( lvshape ) then
    lstress = .true.
    nskip_stress = 1
end if

!-----hold the original mixing parameters
aslh0 = aslh

!-----check structural optimization
if( ifmd == 1 .and. ioptmze == 1 ) then
    ifmd = 2
    treq = 0.d0
!          ioptmze_cell = -1
end if
!if( ifmd == 1 .and. ioptmze == 10  .or.  &
!&   ifmd == 1 .and. ioptmze == 11 ) then
!    lsave = .false.
!!          ioptmze_cell = -1
!    lhma  = .true.
!    if( lstart ) then
!        iscfmx = 0        ! no SCF iteration is needed
!        lstress = .false.
!        lmulken = .false.
!        lwannier = .false.
!        lstart = .false.
!    end if
!end if
!if( ifmd == 1 .and. ioptmze == 20 ) then
!    lsave = .false.
!!          ioptmze_cell = -1
!    lhma  = .true.
!!    if( lstart ) then
!!        iscfmx = 0        ! no SCF iteration is needed
!!        lstress = .false.
!!        lmulken = .false.
!!        lwannier = .false.
!!        lstart = .false.
!!    end if
!end if

treq = treq/(tempau/2.d0)
!            [K] -> [Ryd.]


!----- node_c and node_r must be 1 or even,
!----- and also mod(nodes,node_c) must be 0.
if( node_c <= 1 ) node_c = 1
if( node_c /= 1 .and. mod(node_c,2) == 1 ) node_c = node_c - 1
if( node_r <= 1 ) node_r = 1
if( node_r /= 1 .and. mod(node_r,2) == 1 ) node_r = node_r - 1

do
   if( mod(nodes,node_c) == 0 ) exit
   node_c = node_c - 2
   if( node_c <= 0 ) node_c = 1
end do

do
   if( mod(nodes,node_r) == 0 ) exit
   node_r = node_r - 2
   if( node_r <= 0 ) node_r = 1
end do


!-----set size of work array
iblock = dblock/8/nband
if( iblock < 100 ) iblock = 100
!-----------------------------------------------------------------------



!----- tfermi in [Ryd.]
tfermi = tfermi * akb / (evdj*hrdev)
tfermi = tfermi * 2.d0             !  [hartree] -> [Ryd.]


!--- total No. of electrons
nel = 0
do i = 1, ntype
   nel = nel + nint(zv(i))*nhk(i)
end do


!----- set No. of atoms
nhk1(1) = 1
nhk2(1) = nhk1(1) + nhk(1) - 1
do i = 2, ntype
   nhk1(i) = nhk2(i-1) + 1
   nhk2(i) = nhk1(i) + nhk(i) - 1
enddo


!----- set atomic mass & covalent radius
do it = 1, ntype
   nzn = nint(zatom(it))
   !----- error trap
   if( nzn <= 0 .or. nzn > 103 ) call fstop( nfile, myid, nodes,  &
&                                'nzn <= 0 .or. nzn > 103' )
   watom(it)  = dmassn(nzn)/avogad/(welm*2.d0)    ! in [Rydberg units]
   covrad(it) = covalent_r(nzn)/audang

   if( rsphexp(it) < 0.001d0 ) then
       if( radsphexp < 0.001d0 ) then
           rsphexp(it) = slater_r(nzn)/audang
         else
           rsphexp(it) = radsphexp
       end if
   end if

   if( radeda(it) < 0.001d0 )  &
&      radeda(it) = slater_r(nzn)/audang
end do


!---well potential
!if( lwell ) then
!    call set_well_radius( nfile, myid, nodes, covrad, ntype )
!end if

!----- check MD/super-cell vectors
!      CALL RCIPRL( h_MD, b, volMD )
!      CALL RCIPRL( hcell, b, volsuper )
!      if( volsuper > volMD )
!     &    call fstop( nfile, myid, nodes,  &
!     &                'error: MD cell must be larger than super cell' )


! if symmetry operation is taken into account, do not shift the center of mass
!if( lsymop ) lplcom = .false.

!clustif: if( lclust ) then
!!----- set and check variables for cluster calculations
!
!   !----- check supercell
!   call check_angle( nfile, myid, nodes, hcell, lorthrhmbc )
!   !----- error trap
!   if( .not.lorthrhmbc )  &
!&      call fstop( nfile, myid, nodes, 'incorrect supercell' )
!
!   do i = 1, 3
!      if( .not.lvacuum(i) ) vacuum(i) = hcell(i,i)
!   end do
!   ldouble_grid = .false.
!   ldouble_grid_recip = .false.
!   do i = 1, 3
!      ldouble_grid = ldouble_grid .or. lvacuum(i)
!   end do
!
!   do i = 1, 3
!      lvacuum(i) = .true.       ! for clusters, lvacuum must be .true.
!   enddo
!   lshdep  = .false.
!   lstress = .false.
!   nskip_stress = 1
!   do i = 1, 3
!      rmax(i)   = hcell(i,i) * 0.5d0
!   enddo
!
!
!   gridif: if( .not.ldouble_grid ) then
!
!      !----- get FD & FFT meshes
!      call get_mesh( nfile, myid, nodes,  &
!& ecutdens, hcell, nd1v, nd1vks )
!
!   else gridif
!
!      !----- get FD & FFT meshes with vacuum
!      call get_meshvac( nfile, myid, nodes,  &
!& lclust, ecutdens, hcell, vacuum, nd1v, nd1vks )
!
!   end if gridif
!
!
!   do i = 1, 3
!      ndisp_org(i) = ( nd1v(i) - nd1vks(i) )/2
!      disp_org(i) = dble(ndisp_org(i)) * hcell(i,i)/dble(nd1v(i))
!   enddo
!
!
!   !----- set the center of supercell
!   do ix = 1, 3
!      rdelr = vacuum(ix)/dble(nd1vks(ix))
!      cofmas(ix) = 0.5d0*( vacuum(ix) - rdelr )
!   end do
!
!
!   !----- set real coordinates for cluster calculations
!    do it = 1, ntype
!       if( icscale(it) == 1 ) then
!          do i = nhk1(it), nhk2(it)
!             q1 = ratm(1,i)
!             q2 = ratm(2,i)
!             q3 = ratm(3,i)
!             ratm(1,i) = h_MD(1,1)*q1+h_MD(1,2)*q2+h_MD(1,3)*q3
!             ratm(2,i) = h_MD(2,1)*q1+h_MD(2,2)*q2+h_MD(2,3)*q3
!             ratm(3,i) = h_MD(3,1)*q1+h_MD(3,2)*q2+h_MD(3,3)*q3
!          end do
!       endif
!    end do
!
!    ccenterif: if( lplcom ) then
!    !----- place the center of mass on the center of cell
!       ttmaxx = 0.d0
!       do it = 1, ntype
!          ttmaxx = ttmaxx + watom(it)*dble(nhk(it))
!       end do
!       do j = 1, 3
!          cmpos = 0.d0
!          do it = 1, ntype
!             do i = nhk1(it), nhk2(it)
!                cmpos = cmpos + ratm(j,i)*watom(it)
!             end do
!          end do
!          cmpos = cmpos / ttmaxx
!
!          do i = 1, nhk2(ntype)
!             ratm(j,i) = ratm(j,i) - cmpos + cofmas(j)
!          end do
!
!       end do
!
!    end if ccenterif
!
!
!    !--- error trap
!    do ix = 1, 3
!       qmax = -1.d+10
!       qmin =  1.d+10
!       do i = 1, nhk2(ntype)
!          qmax = max( qmax, ratm(ix,i) )
!          qmin = min( qmin, ratm(ix,i) )
!       end do
!       rdelr = vacuum(ix)/dble(nd1vks(ix))
!       qdis  = vacuum(ix) - rdelr
!       if( qmax - qmin > qdis )  &
!&          call fstop( nfile, myid, nodes,  &
!&                      'error in real coordinates' )
!
!       if( qmax >= qdis .or. qmin < 0.d0 ) then
!           cmpos = 0.5d0*(qmax + qmin)
!           do i = 1, nhk2(ntype)
!              ratm(ix,i) = ratm(ix,i) - cmpos + cofmas(ix)
!           end do
!       end if
!    end do
!
!else clustif
!----- set and check variables for bulk calculations

   do i = 1, 3
      if( .not.lvacuum(i) ) vacuum(i) = hcell(i,i)
   end do
   ldouble_grid = .false.
   ldouble_grid_recip = .false.
   do i = 1, 3
      ldouble_grid_recip = ldouble_grid_recip .or. lvacuum(i)
      if( lvacuum(i) )  &
&         idouble_grid_method = idouble_grid_method + 1
   end do
   lshdep = .false.
   ndisp_org(1) = 0
   ndisp_org(2) = 0
   ndisp_org(3) = 0
   disp_org(1)  = 0.d0
   disp_org(2)  = 0.d0
   disp_org(3)  = 0.d0
   if( ldouble_grid_recip ) lstress = .false.
   lsphere = .false.
   rmax(1) = 0.d0
   rmax(2) = 0.d0
   rmax(3) = 0.d0
!---Do the calculation for charged-periodic calculations with the uniform background
!---         if( idouble_grid_method /= 3 ) ncion = 0  ! if not cluster, set ion = 0

   !----- set ldoublegr
!   call set_ldoublegr( nfile, myid, nodes,  &
!& ldouble_grid_recip, lrtddft, llcexchange )

   !----- check supercell
!   call check_supercell( nfile, myid, nodes,  &
!& ldouble_grid_recip, idouble_grid_method, lvacuum, hcell,  &
!& lorthrhmbc, ierror )
   lorthrhmbc = .false.

   !----- error trap
   if( ierror /= 0 )  &
&      call fstop( nfile, myid, nodes, 'incorrect supercell' )


!   gridif2: if( .not.ldouble_grid_recip ) then

      !----- get FD & FFT meshes
      call get_mesh( nfile, myid, nodes,  &
& ecutdens, hcell, nd1v, nd1vks )

      !----- get double grid if lrtddft .and. llcexchange
!      call get_mesh_llcex( nfile, myid, nodes,  &
!& ecutdens, ecutlong, idouble_grid_method, large_cell_llcexchange,  &
!& hcell, nd1v, nd1vks, iosp )

!   else gridif2
!
!      !----- get FD & FFT meshes with vacuum
!      call get_mesh_bulkvac( nfile, myid, nodes,  &
!& ecutdens, ecutlong, ecutlong_small, idouble_grid_method,  &
!& hcell, lvacuum, vacuum, nd1v, nd1vks, disp_org, iosp )
!
!   end if gridif2


   !----- set the center of supercell
   if( ldouble_grid_recip ) then
      do ix = 1, 3
         cofmas(ix) = 0.d0
      end do
      do i = 1, 3
         if( lvacuum(i) ) then
             rdelx = hcell(1,i)/dble(nd1v(i))
             rdely = hcell(2,i)/dble(nd1v(i))
             rdelz = hcell(3,i)/dble(nd1v(i))
!                   cofmas(1) = cofmas(1)
!     &                       + 0.5d0*( rdelx*dble(nd1vks(i)-1) )
!                   cofmas(2) = cofmas(2)
!     &                       + 0.5d0*( rdely*dble(nd1vks(i)-1) )
!                   cofmas(3) = cofmas(3)
!     &                       + 0.5d0*( rdelz*dble(nd1vks(i)-1) )
             cofmas(1) = cofmas(1) + 0.5d0*rdelx*dble(nd1vks(i))
             cofmas(2) = cofmas(2) + 0.5d0*rdely*dble(nd1vks(i))
             cofmas(3) = cofmas(3) + 0.5d0*rdelz*dble(nd1vks(i))
         else
             do ix = 1, 3
                cofmas(ix) = cofmas(ix) + 0.5d0*hcell(ix,i)
             end do
         end if
      end do
   else
      do ix = 1, 3
        cofmas(ix) = 0.5d0*( hcell(ix,1)+hcell(ix,2)+hcell(ix,3) )
      end do
   end if

   !----- set scaled coordinates for bulk calculations
    !----- transpose matrix of b = inverse of hcell
    CALL RCIPRL( hcell, b, rvol )
    !----- scaled position of the center of supercell
    sc(1) = b(1,1)*cofmas(1) + b(2,1)*cofmas(2) + b(3,1)*cofmas(3)
    sc(2) = b(1,2)*cofmas(1) + b(2,2)*cofmas(2) + b(3,2)*cofmas(3)
    sc(3) = b(1,3)*cofmas(1) + b(2,3)*cofmas(3) + b(3,3)*cofmas(3)


    !----- get real coordinates
    do it = 1, ntype
       if( icscale(it) == 1 ) then
          do i = nhk1(it), nhk2(it)
             q1 = ratm(1,i)
             q2 = ratm(2,i)
             q3 = ratm(3,i)
             ratm(1,i) = h_MD(1,1)*q1+h_MD(1,2)*q2+h_MD(1,3)*q3
             ratm(2,i) = h_MD(2,1)*q1+h_MD(2,2)*q2+h_MD(2,3)*q3
             ratm(3,i) = h_MD(3,1)*q1+h_MD(3,2)*q2+h_MD(3,3)*q3
          end do
       else if( ifmd == 10 ) then
          !---for MSST, 
          call matinv33( hcell_rec, hinv )
          do i = nhk1(it), nhk2(it)
             x = ratm(1,i)
             y = ratm(2,i)
             z = ratm(3,i)
             q1 = hinv(1,1)*x + hinv(1,2)*y + hinv(1,3)*z
             q2 = hinv(2,1)*x + hinv(2,2)*y + hinv(2,3)*z
             q3 = hinv(3,1)*x + hinv(3,2)*y + hinv(3,3)*z
             ratm(1,i) = h_MD(1,1)*q1+h_MD(1,2)*q2+h_MD(1,3)*q3
             ratm(2,i) = h_MD(2,1)*q1+h_MD(2,2)*q2+h_MD(2,3)*q3
             ratm(3,i) = h_MD(3,1)*q1+h_MD(3,2)*q2+h_MD(3,3)*q3
          end do
       endif
    end do

    !----- get scaled coordinates
    small = 1.d-14
    do i = 1, nhk2(ntype)
       x = ratm(1,i)
       y = ratm(2,i)
       z = ratm(3,i)
       ratm(1,i) = b(1,1)*x + b(2,1)*y + b(3,1)*z
       ratm(2,i) = b(1,2)*x + b(2,2)*y + b(3,2)*z
       ratm(3,i) = b(1,3)*x + b(2,3)*y + b(3,3)*z
       if( abs(ratm(1,i)) < small ) ratm(1,i) = 0.d0
       if( abs(ratm(2,i)) < small ) ratm(2,i) = 0.d0
       if( abs(ratm(3,i)) < small ) ratm(3,i) = 0.d0
       if( abs(ratm(1,i)-1.d0) < small ) ratm(1,i) = 0.d0
       if( abs(ratm(2,i)-1.d0) < small ) ratm(2,i) = 0.d0
       if( abs(ratm(3,i)-1.d0) < small ) ratm(3,i) = 0.d0
    end do

    centerif: if( lplcom ) then
    !----- place the center of mass on the center of cell
       ttmaxx = 0.d0
       do it = 1, ntype
          ttmaxx = ttmaxx + watom(it)*dble(nhk(it))
       end do
       do j = 1, 3
          cmpos = 0.d0
          do it = 1, ntype
             do i = nhk1(it), nhk2(it)
                cmpos = cmpos + ratm(j,i)*watom(it)
             end do
          end do
          cmpos = cmpos / ttmaxx

       !-----   if( .not.lvacuum(j) ) then
             do i = 1, nhk2(ntype)
                ratm(j,i) = ratm(j,i) - cmpos + sc(j)
               !--- if( ratm(j,i) <0.d0 ) ratm(j,i) = ratm(j,i) + 1.d0
               !--- if( ratm(j,i)>=1.d0 ) ratm(j,i) = ratm(j,i) - 1.d0
             end do
       !-----   end if

       end do

    end if centerif


    !--- error trap
    do ix = 1, 3
       qmax = -1.d+10
       qmin =  1.d+10
       do i = 1, nhk2(ntype)
          qmax = max( qmax, ratm(ix,i) )
          qmin = min( qmin, ratm(ix,i) )
       end do
       if( lvacuum(ix) ) then
           qdis = dble(nd1vks(ix))/dble(nd1v(ix))
         else
           qdis = 1.d0
       end if
       if( qmax - qmin > qdis )  &
&          call fstop( nfile, myid, nodes,  &
&                      'error in scaled coordinates' )

       if( qmax >= qdis .or. qmin < 0.d0 ) then
           cmpos = 0.5d0*(qmax + qmin)
           do i = 1, nhk2(ntype)
              ratm(ix,i) = ratm(ix,i) - cmpos + sc(ix)
           end do
       end if
    end do

!end if clustif


!----- the grid spacing : rdel
do i = 1, 3
   rdel(i)  = hcell(i,i)/dble(nd1v(i))
enddo

do i = 1, 3
do j = 1, 3
   rdelg(j,i) = hcell(j,i)/dble(nd1v(i))
enddo
enddo


!----- the volume element for dense grid : rdelv
CALL RCIPRL( hcell, hci, volume )
rdelv = volume/dble( nd1v(1)*nd1v(2)*nd1v(3) )


!----- total No. of electrons, subtract the charge number of ion
nel = nel - ncion


!----- set logical variables releted to the calculation of pseudopotential

!-----if lpaw = .true. (the PAW method is used),
!-----lkbpp & llocl must be .false.
if( lpaw.and.lkbpp .or. lpaw.and.llocl )  &
&    call fstop( nfile, myid, nodes,  &
&                'NCPP cannot be used in the PAW method.' )

!      !-----real space calculatoin of ultrasoft pp. is not supported yet
!      do it = 1, ntype
!         if( lvandi(it) ) lking(it) = .false.
!      end do

!-----for cluster calculations, local pp must be calculated in real space
if( lclust .or. lshdep ) then
    do it = 1, ntype
    if( .not.llking(it) ) then
       llking(it)   = .true.
       rlking(it)   = 1.5d0
       glkgmax(it)  = 1.15d0
       glkgexct(it) = 0.8d0
       do i = 1, 2
          if( loutfile(i) ) then
              write(nfile(i),*)  &
&             '*** Since .not.llking(it) for it=', it, ','          
              write(nfile(i),*)  &
&                   '    reset llking,rlking,glkgmax,glkgexct =',  &
&             llking(it), rlking(it), glkgmax(it), glkgexct(it)
          end if
       end do
    end if
    end do
end if

!-----for bulk calculations, local pp must be calculated in reciprocal space
if( ldouble_grid_recip ) then
    do it = 1, ntype
       llking(it) = .false.
    end do
end if

!-----reciprocal space calculatoin of PCC is not supported yet, -> now supported 5/7/2011
!do it = 1, ntype
!   if( .not.lpking(it) ) then
!       lpking(it)   = .true.
!       rpking(it)   = 1.5d0
!       gpkgmax(it)  = 1.15d0
!       gpkgexct(it) = 0.8d0
!       if( myid == 0 ) then
!           do i = 1, 2
!              write(nfile(i),*)  &
!&             '*** Since .not.lpking(it) for it=', it, ','          
!              write(nfile(i),*)  &
!&                   '    reset lpking,rpking,gpkgmax,gpkgexct =',  &
!&             lpking(it), rpking(it), gpkgmax(it), gpkgexct(it)
!           end do
!       end if
!   end if
!end do


lkbpp_r  = .false.
lkbpp_g  = .false.
lvand_r  = .false.
lvand_g  = .false.
llclpp_r = .false.
llclpp_g = .false.
lpcc_r   = .false.
lpcc_g   = .false.
do it = 1, ntype
   lkbpp_r  = lkbpp_r  .or. ( lkbppi(it) .and.      lking(it) )
   lkbpp_g  = lkbpp_g  .or. ( lkbppi(it) .and. .not.lking(it) )
   lvand_r  = lvand_r  .or. ( lvandi(it) .and.      lking(it) )
   lvand_g  = lvand_g  .or. ( lvandi(it) .and. .not.lking(it) )
   llclpp_r = llclpp_r .or.      llking(it)
   llclpp_g = llclpp_g .or. .not.llking(it)
   lpcc_r   = lpcc_r   .or.      lpking(it)
   lpcc_g   = lpcc_g   .or. .not.lpking(it)
end do
lnlpp_r  = lkbpp_r .or. lvand_r
lnlpp_g  = lkbpp_g .or. lvand_g


!-----lconduct = .true., if at least one of ldcconduct and lacconduct is .true.
lconduct = lconduct .and. ldcconduct .or. lconduct .and. lacconduct


!-----for multiple k-point calculation, the following analyses are not supported yet.
if( .not.lgamma ) then
    lmulken  = .false.
    lwannier = .false.
    ldpwav   = .false.
    lscissors= .false.
    if( lefield .and. .not.lsawtooth ) lefield = .false.
end if


!-----for noncollinear magnetism, the following analyses are not supported yet.
if( lnoncollinear ) then
    if( .not.lgamma ) then
        lmulken  = .false.
    end if
    lwannier = .false.
    ldpwav   = .false.
    lscissors= .false.
    lconduct = .false.
    if( lefield .and. .not.lsawtooth ) lefield = .false.
end if


!-----for tddft, the following analyses are not supported yet.
if( ltddft ) then
    lmulken  = .false.
    lwannier = .false.
    ldpwav   = .false.
    lconduct = .false.
    if( ifmd ==4 .or. ifmd == 10 ) call fstop( nfile, myid, nodes,  &
&                      'TDDFT with NPT ensemble is not supported yet' )
end if
if( ltddft .or. ltddft_fssh ) then
    if( .not.lgamma ) call fstop( nfile, myid, nodes,  &
&                      'TDDFT with multiple k-point is not supported yet' )
#if PCHOLESKY
    call fstop( nfile, myid, nodes,  &
&                    'TDDFT with MACRO PCHOLESKY is not supported yet' )
#endif
end if


!-----uniform electric field
!if( lefield ) then
!    efield(1:3) = efield(1:3) * 2.d0 ! [hartree] -> [Ryd.]
!
!    if( lsawtooth ) lconstraintD = .false.
!
!    !---store efield as the displacement vector
!    if( lconstraintD ) call set_dfield_in_efield( efield )
!
!    if( ifmd == 0 ) loutpolarization = .true.
!
!    if( lsawtooth .and. lsawtooth_shape ) then
!!        !----- check supercell
!!        call check_angle( nfile, myid, nodes, hcell, lorthrhmbc )
!!        !----- error trap
!!        if( .not.lorthrhmbc )  &
!!&           call fstop( nfile, myid, nodes, 'supercell should be rectrangular' )
!
!        !-----periodic sawtooth potential
!        lsawtooth_xyz(1) = abs(efield(1)) > 1.d-30
!        lsawtooth_xyz(2) = abs(efield(2)) > 1.d-30
!        lsawtooth_xyz(3) = abs(efield(3)) > 1.d-30
!    end if
!    if( .not.lsawtooth ) lsawtooth_shape = .false.
!    if( .not.lstart ) lefield_start = .false.
!end if
lefield_islts = lefield .and. .not.lsawtooth


!-----for ultrasoft pp, wannier calculation has been supported.
!-----if( lvand ) lwannier = .false.


!-----check the existence of files for atomic orbitals
!call check_mulliken( nfile, myid, nodes,  &
!& jgga, ntype, zatom, lkbppi, lvandi, lpaw,  &
!& rdecmp, rxmulk, lmulao, lmulbsis, nmulao, mulaox,  &
!& lchk, lmax, mxl, ms, aname, ierror )
!if( ierror /= 0 )  &
!&      call fstop( nfile, myid, nodes, 'error: check_mulliken' )


!-----set lsphexp in sphexp.f90
!!if( .not.lnlpp_g .and. .not.llclpp_g .and. .not.lpcc_g ) then
!if( .not.lnlpp_g .and. .not.lpcc_g ) then
!    !--- Stop spherical harmonics expansion, because (ycos, ysin) are not allocated.
!    if( lsphexp .and. myid == 0 ) then
!        do i = 1, 2
!           write(nfile(i),*)  &
!              '*** Since both lnlpp_g and lpcc_g are .false.,', &
!              ' the spherical harmonics expansion is stopped.'
!        end do
!    end if
!    lsphexp = .false.
!end if
!call set_lsphexp( lsphexp )


!-----if lmulken = .false., set leda = .false.
if( .not.lmulken ) leda = .false.
if( .not.lpaw )    leda = .false.
!call setleda( leda )


!-----set ltddft in tddft_fssh.f90
!----- tdbroad in [Ryd.]
tdbroad = tdbroad * akb / (evdj*hrdev)
tdbroad = tdbroad * 2.d0             !  [hartree] -> [Ryd.]
call set_ltddft( nfile, myid, nodes, &
& ltddft, ltddft_fssh, ltddft_start, nel, dttddft, lfssh_switch, &
& lfssh_gsscf, tdbroad, lfssh_random, rseed_fssh, lfssh_boltzmn, treq, &
& lfssh_vscale, tminimum, lfssh_parallel, ltddft_nscforce, lrtddft )


!-----conductivity calculation
if( lconduct ) then
    !-----unit change
    !----- tempconduct [K] -> [Ryd.]
    tempconduct = tempconduct * akb / (evdj*hrdev)
    tempconduct = tempconduct * 2.d0             !  [hartree] -> [Ryd.]

    !----- electric-field
!    efconduct(1:3) = efconduct(1:3) * ??

!    !----- frequencies [eV] -> [Ryd.]
!    freqacmx  = freqacmx / hrdev * 2.d0

    if( .not.lacconduct ) freqacmx = 0.d0
else
    ldcconduct = .false.
    lacconduct = .false.
end if
!-----set lconduct in pppaw.f90
!call set_lconduct_in_pppaw( nfile, myid, nodes, lconduct )


!-----set lsreal8 in module savedata
call set_lsreal8( nfile, myid, nodes, lsreal8 )


!-----set lnoncollinear in commwv.f90
call set_lnoncollinear_in_commwv( nfile, myid, nodes, lnoncollinear )

!-----set lnoncollinear in engrad.f90
call set_lnoncollinear_in_engrad( nfile, myid, nodes, lnoncollinear )

!-----set lnoncollinear in eigen.f90
call set_parameters_in_eigen( nfile, myid, nodes,  &
& lnoncollinear, roughfeig, roughfzan )

!-----set lnoncollinear in fermi.f90
call set_lnoncollinear_in_fermi( nfile, myid, nodes, lnoncollinear )

!-----set lnoncollinear in mulliken.f90
!call set_lnoncollinear_in_mulliken( nfile, myid, nodes, lnoncollinear )

!-----set lmixreal in module ncmagne_mixing
!call set_lmixreal( nfile, myid, nodes, lmixreal )

if( lrela ) then
    if( .not.lnoncollinear ) lrela_full = .false.
    if( lrela_full ) lrela_sop = .false.
else
    lrela_full = .false.
    lrela_sop  = .false.
end if


if( lrela .or. .not.lpaw ) then
    lsoc = .false.
end if
if( lsoc ) then
    if( .not.lnoncollinear ) lsoc_full = .false.
end if

!-----set lrela in module SOp_variables
!call set_lrela_in_so( nfile, myid, nodes,  &
!& lrela, lrela_full, lrela_sop, lsoc, lsoc_full )


!---if single k-point calculation or MD, stop dispersion relation calculation
!call check_dispersion( nfile, myid, nodes, nkpnt, ifmd )


return
end subroutine




subroutine get_mesh( nfile, myid, nodes,  &
& ecutdens, hcell, nd1v, nd1vks )
!-----------------------------------------------------------------------
!    get No. of FD & FFT meshes
!        in the case that FD & FFT meshes are the same.
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
real*8  :: ecutdens
real*8,  dimension(3,3) :: hcell
integer, dimension(3) :: nd1v, nd1vks

real*8,  dimension(3,3) :: B
real*8  :: rvol, dpi, GXFFT
real*8  :: A1X, A1Y, A1Z, tri, AAA, PK3
integer, dimension(3) :: IIP1 = (/ 2, 3, 1 /)
integer, dimension(3) :: IIP2 = (/ 3, 1, 2 /)
integer, dimension(3) :: IIP3 = (/ 1, 2, 3 /)
integer :: i, IP1, IP2, IP3, KMAXM


!   ... reciprocal lattice ...
CALL RCIPRL( hcell, B, rvol )

dpi = 2.d0*acos(-1.d0)
GXFFT  = sqrt( ecutdens )/dpi


do i = 1, 3
   IP1 = IIP1(i)
   IP2 = IIP2(i)
   IP3 = IIP3(i)
   A1X = B(2,IP1)*B(3,IP2) - B(3,IP1)*B(2,IP2)
   A1Y = B(3,IP1)*B(1,IP2) - B(1,IP1)*B(3,IP2)
   A1Z = B(1,IP1)*B(2,IP2) - B(2,IP1)*B(1,IP2)
   tri = B(1,IP3)*A1X + B(2,IP3)*A1Y + B(3,IP3)*A1Z
   AAA = A1X*A1X + A1Y*A1Y + A1Z*A1Z
   PK3 = tri/SQRT(AAA)
   KMAXM = GXFFT/PK3
   nd1vks(i) = 2*KMAXM + 1
end do


!--- check FFT mesh
call chkftm( nfile, myid, nd1vks(1), nd1vks(2), nd1vks(3), .true. )

nd1v(1) = nd1vks(1)
nd1v(2) = nd1vks(2)
nd1v(3) = nd1vks(3)


return
end subroutine




subroutine write_data( nfile, myid, nodes, lclust )
!-----------------------------------------------------------------------
!    write input data
!-----------------------------------------------------------------------
use outfile
use param
use param_inispin
use param_atom
use constants
!use symmop
use input_variables
implicit none
integer :: nfile(*), myid, nodes
logical :: lclust
integer :: i, it, l, j, ix, k
real*8  :: rvdwmax
integer :: imod1, imod2, iosp_vdw
real*8  :: q_cut, dlambda, d_soft, phi_at_0
integer :: nqmesh


do i = 1, 2

   myidif: if( loutfile(i) ) then

      write(nfile(i),*) ' '
!      if( lclust ) then
!         write(nfile(i),*) ' ---- Ab initio calculation',  &
!&                          ' for atomic cluster by PW method ----'
!      else
         write(nfile(i),*) ' ---- Ab initio calculation',  &
&                          ' for bulk by PW method ----'
!      end if
      write(nfile(i),*) ' '
!      if( lpaw ) then
!         write(nfile(i),*) '         <<< The projector-augmented wave method >>> '
!      else
         write(nfile(i),*) '         <<< The standard pseudopotential method >>> '
!      end if
!      if( lrela ) then
!          if( lrela_full ) then
!             write(nfile(i),*) '             using fully relativistic two-component spinors'
!          else
!             write(nfile(i),*) '             within scalar-relativistic approach'
!          end if
!      end if
!      if( lsoc ) then
!          if( lsoc_full ) then
!             write(nfile(i),*) '             with self-consistent spin-orbit coupling interaction'
!          else
!             write(nfile(i),*) '             with perturbative spin-orbit coupling interaction'
!          end if
!      end if
      if( ldftd ) then
         write(nfile(i),*) '             with an empirical correction for vdW (DFT-D)'
      end if
!      if( lplusU ) then
!         write(nfile(i),*) '             with the mean-field Hubbard model (DFT+U)'
!      end if
!      if( lplusC ) then
!         write(nfile(i),*) '             with an empirical correction to non-local pp. (DFT+C)'
!      end if
!      if( jhybrid >= 1 ) then
!         write(nfile(i),*) '             with range-separated hybrid functionals'
!      end if
!      if( lnoncollinear ) then
!         write(nfile(i),*) '             with the noncollinear magnetism'
!      end if
      if( ltddft ) then
         write(nfile(i),*) '             based on time-dependent', &
&                          ' density-functional theory (TDDFT) '
      end if
      if( ltddft_fssh ) then
         write(nfile(i),*) '             based on time-dependent', &
&                          ' density-functional theory (TDDFT) '
         write(nfile(i),*) '             with fewest switching surface hopping (FSSH) '
         if( lrtddft ) then
             write(nfile(i),*) '             for the EXCITON diffusion'
         end if
      end if
!      if( lefield ) then
!         write(nfile(i),*) '             in the presence of uniform electric field'
!      end if

      write(nfile(i),*) ' '
!      if( lpaw ) then
!          write(nfile(i),'(1x,a12,l6,6x,a12,l6,6x,a12,l6)')  &
!&                              '      lpaw :', lpaw,  &
!&                              '  lpaw_sym :', lpaw_sym
!      end if
!      if( lrela ) then
!          write(nfile(i),'(1x,a12,l6,6x,a12,l6,6x,a12,l6)')  &
!&                              '     lrela :', lrela,  &
!&                              'lrela_full :', lrela_full,  &
!&                              ' lrela_sop :', lrela_sop
!      end if
!      if( lsoc ) then
!          write(nfile(i),'(1x,a12,l6,6x,a12,l6,6x,a12,l6)')  &
!&                              '      lsoc :', lsoc,  &
!&                              ' lsoc_full :', lsoc_full
!      end if
      write(nfile(i),2606) '    lclust :', lclust,  &
&                          '      jgga :', jgga_org,  &
&                          '    lstart :', lstart
!      if( jhybrid >= 1 ) then
!          write(nfile(i),'(1x, a12,es11.4,1x, a12,es11.4,1x, a12,i6)')  &
!&                              '   hmixing :', hmixing,  &
!&                              '    hrange :', hrange,  &
!&                              '  mgridred :', mgridred
!      end if
      if( ltddft .or. ltddft_fssh ) then
          write(nfile(i),2646) 'ltddft_star:', ltddft_start, &
&                              '   dttddft :', dttddft, &
&                              'lfssh_switc:', lfssh_switch
          write(nfile(i),'(1x, a12,l6,6x,    a12,l6,6x, a12,f10.2)')  &
&                              'lfssh_gsscf:', lfssh_gsscf,  &
&                              'lfssh_gsfrc:', lfssh_gsfrc
          write(nfile(i),2044) 'imxchg_fssh:', imxchg_fssh,  &
&                              ' aslh_fssh :', aslh_fssh, &
&                              ' bslh_fssh :', bslh_fssh
          write(nfile(i),'(1x, a12,l6,6x,    a12,f10.2,2x,a12,l6)') &
&                              'lfssh_rando:', lfssh_random, &
&                              '   tdbroad :', tdbroad/(akb/(evdj*hrdev))/2.d0, &
&                              'lfssh_paral:', lfssh_parallel
          write(nfile(i),*)    'rseed_fssh :', rseed_fssh
          write(nfile(i),'(1x, a12,l6,6x,    a12,l6,6x, a12,f10.2)') &
&                              'lfssh_boltz:', lfssh_boltzmn, &
&                              'lfssh_vscal:', lfssh_vscale, &
&                              '  tminimum :', tminimum
          write(nfile(i),'(1x, a12,l6,6x,    a12,l6,6x, a12,f10.2)') &
&                              'ltddft_nscf:', ltddft_nscforce
      end if
      if( lrtddft ) then
          write(nfile(i),'(1x, a12,l6,6x,    a12,l6,6x, a12,es11.4)') &
&                              '   lrtddft :', lrtddft,  &
&                              'lrspecific :', lrspecific,  &
&                              '  fdmeshxc :', fdmeshxc
          write(nfile(i),'(1x, a12,es11.4,1x, a12,es11.4,1x, a12,es11.4)') &
&                              ' threshrpa :', threshrpa,  &
&                              'threshdipol:', threshdipole,  &
&                              'threshexcit:', threshexcite
          write(nfile(i),'(1x, a12,l6,6x,    a12,l6,6x, a12,es11.4)') &
&                              '      ltda :', ltda,  &
&                              '  lcic_tda :', lcic_tda,  &
&                              'cic_tda_cou:', cic_tda_coupling
          write(nfile(i),'(1x, a12,l6,6x, a12,es11.4,1x, a12,es11.4)') &
&                              'llcexchange:', llcexchange,  &
&                              'alpha_llcex:', alpha_llcexchange,  &
&                              'large_cell_:', large_cell_llcexchange
          write(nfile(i),'(1x, a12,es11.4,1x, a12,es11.4,1x,a12,l6)') &
&                              'hybrid_mixi:', hybrid_mixing,  &
&                              ' foersigma :', foersigma,  &
&                              ' ldiagonal :', ldiagonal
          write(nfile(i),'(1x, a12,l6,6x, a12,l6,6x, a12,es11.4)') &
&                              'lfission_ra:', lfission_rate,  &
&                              'lemission_r:', lemission_rate,  &
&                              'refractive :', refractive
          write(nfile(i),'(1x, a12,l6,6x,    a12,es11.4,1x, a12,es11.4)') &
&                              'lscissors_l:', lscissors_lrtddft,  &
&                              'ecorr_tirpl:', ecorr_tirplet,  &
&                              'ecorr_singl:', ecorr_singlet
          if( .not.lrspecific ) then
              write(nfile(i),'(1x, a12,es11.4,1x)') &
&                              '   enediff :', enediff
          end if
      end if
!      if( lvdw ) then
!          write(nfile(i),'(1x, a12,l6,6x,    a12,l6,6x, a12,es11.4)') &
!&                              '      lvdw :', lvdw,  &
!&                              'lvdw_factor:', lvdw_factorize
!          if( lvdw_factorize ) then
!              call get_qrange( q_cut, dlambda, nqmesh )
!              call get_d_soft( d_soft, phi_at_0 )
!              write(nfile(i),'(1x, a12,es11.4,1x, a12,es11.4,1x, a12,i6)') &
!&                              '     q_cut :', q_cut,  &
!&                              '   dlambda :', dlambda,  &
!&                              '    nqmesh :', nqmesh
!              write(nfile(i),'(1x, a12,es11.4,1x, a12,es11.4,1x, a12,i6)') &
!&                              '    d_soft :', d_soft,  &
!&                              '  phi_at_0 :', phi_at_0
!          else
!              call get_rvdwmax( rvdwmax )
!              call get_imod( imod1, imod2 )
!              call get_iosp( iosp_vdw )
!              write(nfile(i),'(1x, a12,es11.4,1x, a12,es11.4,1x, a12,es11.4)') &
!&                                  '   rvdwmax :', rvdwmax
!              write(nfile(i),2000) '     imod1 :', imod1,  &
!&                                  '     imod2 :', imod2,  &
!&                                  '      iosp :', iosp_vdw
!          end if
!      end if
      if( lscissors ) then
          write(nfile(i),'(1x, a12,l6,6x,    a12,es11.4,1x, a12,es11.4)') &
&                              ' lscissors :', lscissors,  &
&                              'ecorr_accep:', ecorr_acceptor,  &
&                              'ecorr_donor:', ecorr_donor
          write(nfile(i),'(1x, a12,f10.5,2x, a12,f10.5)')  &
&                              ' zboundary :', zboundary
      end if
      write(nfile(i),2002) ' electrons :', nel,  &
&                          '     ntype :', ntype
      write(nfile(i),2002) '    iscfmx :', iscfmx,  &
&                          '    itrial :', itrial
      write(nfile(i),'(1x, a12,es11.4,1x, a12,es11.4,1x, a12,es11.4)')  &
&                          '    tolpot :', tolpot,  &
&                          '    tolres :', tolres,  &
&                          '  tolHF_KS :', tolHF_KS
      write(nfile(i),2602) '    lhcunt :', lhcunt,  &
&                          '    ihldam :', ihldam,  &
&                          '    toleig :', toleig
      write(nfile(i),'(1x, a12,f10.5,2x, a12,f10.5,2x, a12,f10.5)') &
&                          ' roughfeig :', roughfeig,  &
&                          ' roughfzan :', roughfzan
      write(nfile(i),2000) '    itermx :', itermx,  &
&                          '   iteremx :', iteremx,  &
&                          '   kstrial :', kstrial
      write(nfile(i),'(1x, a12,i6,6x, a12,es11.4,1x, a12,es11.4)')  &
&                          '  methodcg :', methodcg,  &
&                          '   threinn :', threinn,  &
&                          '   wdinner :', wdinner
      write(nfile(i),2044) '    imxchg :', imxchg,  &
&                          '      aslh :', aslh,  &
&                          '      bslh :', bslh
      write(nfile(i),2044) 'imxchg_tria:', imxchg_trial,  &
&                          'aslh_trial :', aslh_trial,  &
&                          'bslh_trial :', bslh_trial
!      if( lpaw ) then
!          if( lplusU ) then
!              write(nfile(i),'(1x, a12,f10.5,2x, a12,f10.5)')  &
!&                          '   paw_mix :', paw_mix,  &
!&                          ' plusU_mix :', plusU_mix
!          else
!              write(nfile(i),2410) '   paw_mix :', paw_mix
!          end if
!      end if
      write(nfile(i),'(1x, a12,i6,6x,    a12,f10.5,2x,    a12,l6)')  &
&                          '    itratn :', itratn,  &
&                          '   cmetric :', cmetric
      write(nfile(i),'(1x, a12,l6,6x,    a12,l6,6x,    a12,l6)')  &
&                          '   louteig :', louteig,  &
&                          'loutenergy :', loutenergy,  &
&                          ' loutzansa :', loutzansa
      write(nfile(i),'(1x, a12,l6,6x,    a12,l6,6x,    a12,l6)')  &
&                          '   loutlog :', loutlog
      write(nfile(i),'(1x, a12,i6,6x,    a12,es11.4,1x,    a12,i6)')  &
&                          '    nstabi :', nstabi,  &
&                          '    xstabi :', xstabi
      if( lspin .or. lnoncollinear ) then
!          if( lnoncollinear ) then
!              write(nfile(i),'(1x, a12,l6,6x,    a12,l6,6x,    a12,l6)')  &
!&                              '  lmixreal :', lmixreal
!          end if
          write(nfile(i),2044) '   imxspin :', imxspin,  &
&                              '   amxspin :', amxspin,  &
&                              '   bmxspin :', bmxspin
          write(nfile(i),'(1x, a12,i6,6x,    a12,f10.5,2x,    a12,f10.5)')  &
&                              '   nmxspin :', nmxspin,  &
&                              'spinmetric :', spinmetric,  &
&                              '   wkerker :', wkerker
          write(nfile(i),2044) 'imxspin_tri:', imxspin_trial,  &
&                              'amxspin_tri:', amxspin_trial,  &
&                              'bmxspin_tri:', bmxspin_trial
      end if
      write(nfile(i),2002) '      nd2v :', nd2v,  &
&                          '     multg :', multg,  &
&                          '     tolcg :', tolcg
      write(nfile(i),2410) '    weigrd :', weigrd,  &
&                          '    iblock :', iblock
      if( ifmd >= 1 ) then
         write(nfile(i),2050) '      ifmd :', ifmd,  &
&                             '      dtmd :', dtmd
!---     &                              '     nstop :', nstop
         if( .not.lstart ) then
             write(nfile(i),2050) ' nstep_ini :', nstep_ini
         end if
         write(nfile(i),2600) '   liscale :', liscale,  &
&                             '    iscnum :', iscnum,  &
&                             '    iscstp :', iscstp
         write(nfile(i),2600) '  lmomzero :', lmomzero
         write(nfile(i),2500) '      treq :', treq*tempau/2.d0,  &
&                             '    ichest :', ichest,  &
&                             '     ihest :', ihest
         write(nfile(i),2666) '     laspc :', laspc,  &
&                             ' laspc_chg :', laspc_chg,  &
&                             '  laspc_wv :', laspc_wv
         write(nfile(i),2604) '   lxlbomd :', lxlbomd,  &
&                             '   kxlbomd :', kxlbomd
         write(nfile(i),2042) 'iscfmx_aspc:', iscfmx_aspc,  &
&                             ' aspc_corr :', aspc_corr,  &
&                             'tolres_aspc:', tolres_aspc
         write(nfile(i),2044) ' nskp_aspc :', nskp_aspc,  &
&                             ' aslh_aspc :', aslh_aspc,  &
&                             ' bslh_aspc :', bslh_aspc
      endif
      write(nfile(i),2660) '     lsave :', lsave,  &
&                          '   lsreal8 :', lsreal8
      write(nfile(i),2606) '   lmulken :', lmulken,  &
&                          '  nskpmulk :', nskpmulk,  &
&                          ' ldecmpovp :', ldecmpovp
      write(nfile(i),2666) '  lspdatom :', lspdatom,  &
&                          '  lpmchgat :', lpmchgat,  &
&                          'lmdecmpdos :', lmdecmpdos
      write(nfile(i),2604) '   lsphexp :', lsphexp,  &
&                          'nskpsphexp :', nskpsphexp,  &
&                          ' radsphexp :', radsphexp
      write(nfile(i),2604) '      leda :', leda,  &
&                          '   nskpeda :', nskpeda,  &
&                          '   radgeda :', radgeda
!            write(nfile(i),2410) ' rdelaunay :', rdelaunay
      write(nfile(i),2600) '  lwannier :', lwannier,  &
&                          '    iwbnd1 :', iwbnd1,  &
&                          '    iwbnd2 :', iwbnd2
      write(nfile(i),2020) '   iterwan :', iterwan,  &
&                          '    tolwan :', tolwan
      write(nfile(i),2000) '  nskpwann :', nskpwann,  &
&                          '    iwstt1 :', iwstt1,  &
&                          '    iwstt2 :', iwstt2
      write(nfile(i),2602) '   loutuni :', loutuni,  &
&                          '    natwan :', natwan
!      if( lconduct ) then
!          write(nfile(i),'(1x, a12,l6,6x,    a12,l6,6x,    a12,l6)')  &
!&                          '  lconduct :', lconduct,  &
!&                          'ldcconduct :', ldcconduct,  &
!&                          'lacconduct :', lacconduct
!          write(nfile(i),'(1x, a12,i6,6x,    a12,f10.5,2x, a12,f10.2)')  &
!&                          'nskpconduct:', nskpconduct,  &
!&                          ' wgconduct :', wgconduct,  &
!&                          'tempconduct:', tempconduct/(akb/(evdj*hrdev))/2.d0
!          write(nfile(i),'(1x, a12,3es12.4)')  &
!&                          ' efconduct :', efconduct(1:3)
!      end if
      write(nfile(i),'(1x, a12,f10.5,2x, a12,f10.5,2x, a12,f10.2)')  &
&                          '  freqacmx :', freqacmx !*hrdev/2.d0
      write(nfile(i),'(1x, a12,i6,6x,    a12,i6,6x,    a12,i6)')  &
&                          '   icband1 :', icband1,  &
&                          '   icband2 :', icband2
      write(nfile(i),'(1x, a12,i6,6x,    a12,i6,6x,    a12,i6)')  &
&                          '  ichband1 :', ichband1,  &
&                          '  ichband2 :', ichband2
      write(nfile(i),'(1x, a12,i6,6x,    a12,i6,6x,    a12,i6)')  &
&                          '  iceband1 :', iceband1,  &
&                          '  iceband2 :', iceband2
      write(nfile(i),'(1x, a12,i6,6x,    a12,i6,6x,    a12,i6)')  &
&                          ' immtband1 :', immtband1,  &
&                          ' immtband2 :', immtband2
      write(nfile(i),2606) '   lstress :', lstress,  &
&                          ' nskip_stre:', nskip_stress,  &
&                          '    lshdep :', lshdep
      write(nfile(i),2600) '   lintchg :', lintchg,  &
&                          ' nskip_intc:', nskip_intchg
      write(nfile(i),2600) '    ldpchg :', ldpchg,  &
&                          ' nskip_dpch:', nskip_dpchg
      if( ldpchg ) then
        write(nfile(i),2555) '   dc_xmin :', dc_xmin,  &
&                            '   dc_xmax :', dc_xmax,  &
&                            '   dc_ymin :', dc_ymin
        write(nfile(i),2555) '   dc_ymax :', dc_ymax,  &
&                            '   dc_zmin :', dc_zmin,  &
&                            '   dc_zmax :', dc_zmax
        write(nfile(i),2600) 'lcompress_d:', lcompress_dpchg,  &
&                            'ndigit_dpch:', ndigit_dpchg
      end if
      write(nfile(i),2600) '    ldpwav :', ldpwav,  &
&                          ' nskip_dpwa:', nskip_dpwav
      if( ldpwav ) then
        write(nfile(i),2555) '   wv_xmin :', wv_xmin,  &
&                            '   wv_xmax :', wv_xmax,  &
&                            '   wv_ymin :', wv_ymin
        write(nfile(i),2555) '   wv_ymax :', wv_ymax,  &
&                            '   wv_zmin :', wv_zmin,  &
&                            '   wv_zmax :', wv_zmax
        write(nfile(i),2600) 'lcompress_d:', lcompress_dpwav,  &
&                            'ndigit_dpwa:', ndigit_dpwav
      end if
      write(nfile(i),2000) '    ibstt1 :', ibstt1,  &
&                          '    ibstt2 :', ibstt2
      write(nfile(i),2600) '    ldppot :', ldppot,  &
&                          ' nskip_dppo:', nskip_dppot,  &
&                          ' nav_dppot :', nav_dppot
      if( ldppot ) then
        write(nfile(i),2555) '   pt_xmin :', pt_xmin,  &
&                            '   pt_xmax :', pt_xmax,  &
&                            '   pt_ymin :', pt_ymin
        write(nfile(i),2555) '   pt_ymax :', pt_ymax,  &
&                            '   pt_zmin :', pt_zmin,  &
&                            '   pt_zmax :', pt_zmax
        write(nfile(i),2600) 'lcompress_d:', lcompress_dppot,  &
&                            'ndigit_dppo:', ndigit_dppot
      end if
      write(nfile(i),2600) '    lhoteh :', lhoteh,  &
&                          'nskip_hoteh:', nskip_hoteh
      if( lhoteh ) then
          write(nfile(i),'(1x,a12,i6,6x,a12,i6,6x,a12,l6)')  &
&                          '    ibhoth :', ibhoth,  &
&                          '    ibhote :', ibhote
      end if
      write(nfile(i),2555) '      ecut :', ecut,  &
&                          '  ecutdens :', ecutdens,  &
&                          '  ecutsoft :', ecutsoft
      write(nfile(i),2555) '  ecutorth :', ecutorth,  &
&                          '     ecutc :', ecutc
      write(nfile(i),2640) '   lvshape :', lvshape,  &
&                          '   pwscale :', pwscale
!      if( ldouble_grid_recip ) then
!        write(nfile(i),2550) '  ecutlong :', ecutlong,  &
!&                            'ecutlong_s :', ecutlong_small,  &
!&                            '      iosp :', iosp
!      end if
      write(nfile(i),'(1x,a12,i6,6x,a12,l6,6x,a12,l6)')  &
&                          '     nkpnt :', nkpnt,  &
&                          '    lgamma :', lgamma,  &
&                          '    ltetra :', ltetra
      write(nfile(i),'(1x,a12,l6,6x,a12,l6,6x,a12,l6)')  &
&                          '  lksgamma :', lksgamma
      if( .not.lnoncollinear ) then
          write(nfile(i),2000) '     nband :', nband,  &
&                              '    noband :', noband,  &
&                              ' cf. nel/2 :', nel/2
      else
          write(nfile(i),2000) '     nband :', nband,  &
&                              '    noband :', noband,  &
&                              '   cf. nel :', nel
      end if
      write(nfile(i),2050) '    lfermi :', lfermi,  &
&                  '    tfermi :', tfermi/(akb/(evdj*hrdev))/2.d0, &
&                          '     ncion :', ncion
      write(nfile(i),2600) '    lplcom :', lplcom
!&                          '    node_c :', node_c,  &
!&                          '    node_r :', node_r
!      write(nfile(i),2000) '       mx1 :', mx1,  &
!&                          '    mx1loc :', mx1loc
      write(nfile(i),2666) '     lkbpp :', lkbpp,  &
&                          '     lvand :', lvand,  &
&                          '     lvfull:', lvfull
      write(nfile(i),2666) '     llocl :', llocl,  &
&                          '     lpcc  :', lpcc
      write(nfile(i),2666) '   lnlpp_r :', lnlpp_r,  &
&                          '  llclpp_r :', llclpp_r,  &
&                          '    lpcc_r :', lpcc_r
      write(nfile(i),2666) '   lnlpp_g :', lnlpp_g,  &
&                          '  llclpp_g :', llclpp_g,  &
&                          '    lpcc_g :', lpcc_g
      write(nfile(i),2666) '   lkbpp_r :', lkbpp_r,  &
&                          '   lkbpp_g :', lkbpp_g
      write(nfile(i),2666) '   lvand_r :', lvand_r,  &
&                          '   lvand_g :', lvand_g
      write(nfile(i),2646) '     lspin :', lspin
      if( lspin ) then
          write(nfile(i),2646) '    lfixud :', lfixud,  &
&                          '    diffud :', diffud,  &
&                          '   lwfrand :', lwfrand
          write(nfile(i),2006) '   inispin :', inispin,  &
&                              '   nspinat :', nspinat,  &
&                              'lduplicate :', lduplicate
      end if
!      if( lnoncollinear ) then
!          write(nfile(i),2006) '   inispin :', inispin,  &
!&                              '   nspinat :', nspinat,  &
!&                              '   latmref :', latmref
!          write(nfile(i),2646) 'lclearamag :', lclearamag
!          write(nfile(i),'(1x,a12,3f8.4)')  &
!&                              ' refmx,y,z :', refmx, refmy, refmz
!          write(nfile(i),'(1x,a,3l5)')  &
!&                              ' lclearmagne_in_x,y,z :',  &
!&              lclearmagne_in_x, lclearmagne_in_y, lclearmagne_in_z
!          write(nfile(i),'(1x,a,3l5)')  &
!&                              '   lfixmagne_in_x,y,z :',  &
!&              lfixmagne_in_x, lfixmagne_in_y, lfixmagne_in_z
!          write(nfile(i),'(1x,a,3f8.4)')  &
!&                              '    fixmagne_in_x,y,z :',  &
!&              fixmagne_in_x, fixmagne_in_y, fixmagne_in_z
!      end if
!      if( lwell ) then
!          write(nfile(i),2646) '     lwell :', lwell,  &
!&                              'wellheight :', wellheight
!      else
!          write(nfile(i),2646) '     lwell :', lwell
!      end if
!      write(nfile(i),'(1x,a12,i6,6x,a12,i6,5x,a13,3l4)')  &
!&                          '   irstrct :', irstrct,  &
!&                          'irstrct_sub:', irstrct_sub,  &
!&                         'lcell_rstrct:', lcell_rstrct(1:3)
!      if( lefield ) then
!          write(nfile(i),'(1x,a12,l6,6x,a12,l6,6x,a12,l6)')  &
!&                          '   lefield :', lefield,  &
!&                          'lefield_sta:', lefield_start
!          write(nfile(i),'(1x,a12,l6,6x,a12,l6,6x,a12,l6)')  &
!&                          'constraintD:', lconstraintD,  &
!&                          ' lsawtooth :', lsawtooth,  &
!&                          'lsawtooth_s:', lsawtooth_shape
!          write(nfile(i),'(1x,a12,l6,6x,a12,l6,6x,a12,l6)')  &
!&                          'loutpolariz:', loutpolarization
!          if( .not.lconstraintD ) then
!              write(nfile(i),'(1x,a12,3f10.6)')  &
!&                          'efield [Ry]:', efield(1:3)
!          else
!              call out_dfield( nfile, myid, nodes, i )
!          end if
!          call out_efield_mix( nfile, myid, nodes, i )
!!          write(nfile(i),'(1x,a12,3l6)')  &
!!&                          'lsawtooth_x:', lsawtooth_xyz(1:3)
!      else
!          write(nfile(i),'(1x,a12,l6,6x,a12,f11.2,1x,a12,es12.6)')  &
!&                          '   lefield :', lefield
!      end if

      write(nfile(i),*) ' '
!      write(nfile(i),3000) '   h(MD) (L_1) :', ( h_MD(j,1),j=1,3)
!      write(nfile(i),3000) '         (L_2) :', ( h_MD(j,2),j=1,3)
!      write(nfile(i),3000) '         (L_3) :', ( h_MD(j,3),j=1,3)
      write(nfile(i),3000) '   hcell (L_1) :', ( hcell(j,1),j=1,3)
      write(nfile(i),3000) '         (L_2) :', ( hcell(j,2),j=1,3)
      write(nfile(i),3000) '         (L_3) :', ( hcell(j,3),j=1,3)
!      write(nfile(i),3010) '    lorthrhmbc :', lorthrhmbc
!      write(nfile(i),3010) '       lvacuum :', ( lvacuum(j),j=1,3)
!      write(nfile(i),3000) '        vacuum :', ( vacuum(j), j=1,3)
!      write(nfile(i),3000) '        cofmas :', ( cofmas(j), j=1,3)
      write(nfile(i),3020) '          nd1v :', ( nd1v(j),   j=1,3)
      write(nfile(i),3020) '        nd1vks :', ( nd1vks(j), j=1,3)
!      write(nfile(i),3020) '    msrhmx,y,z :',  &
!&                                           msrhmx, msrhmy, msrhmz
!      write(nfile(i),3010) '       lsphere :', lsphere
!      write(nfile(i),3010) '  ldouble_grid :', ldouble_grid
!      write(nfile(i),3010) 'ldouble_grid_r :', ldouble_grid_recip
!      write(nfile(i),3020) 'idouble_grid_m :', idouble_grid_method
!      write(nfile(i),3020) '     ndisp_org :',(ndisp_org(j),j=1,3)
!      write(nfile(i),3000) '      disp_org :',( disp_org(j),j=1,3)
!      write(nfile(i),3000) 'alpha_ldouble_g:', alpha_ldouble_grid_recip

!      if( nkpnt > 0 ) then
!          write(nfile(i),*) ' '
!          write(nfile(i),3100)
!          do k = 1, nkpnt
!          if( k <= 10 .or. k >= nkpnt-10 ) then
!             write(nfile(i),3110) k,(bzk(j,k),j=1,3),wbzk(k),lgammak(k),kpsymop(0,k)
!             if( losymop .and. allocated(kpdupnm) ) then
!                 do j = 1, kpdupnm(k)
!                    write(nfile(i),'(8x,i2,a,i3,a,14i4)')  &
!& j, '. (', kpsymop(j,k),' ) ', ( kpdupsymop(l,j,k),l=1,min(14,kpdupsymop(0,j,k)) )
!                    if( kpdupsymop(0,j,k) >= 15 ) write(nfile(i),'((3x,18i4))')  &
!&                               ( kpdupsymop(l,j,k),l=15,kpdupsymop(0,j,k) )
!                 end do
!                 do j = 1, kpdupnm(k)
!                    write(nfile(i),'(8x,i2,a,3i5,a)')  &
!& j, '. (', ibzk(1:3,j,k),' ) '
!                 end do
!             end if
!             if( k == 10 .and. k /= nkpnt ) write(nfile(i),'(a)') '     ......'
!          end if
!          end do
!          if( losymop ) then
!              write(nfile(i),'(2x,a)') 'symmetry operations'
!              do k = 1, nkpnt
!                write(nfile(i),'(i7,a,(16i4))') k, ':', ( kpsymop(j,k), j = 1, kpsymop(0,k) )
!              end do
!          end if
!
!3100 format(2x,'k point     bzk                        wbzk   lgammak kpsymop(0,k)')
!3110 format(i7,':   (',3f8.4,' )',f9.4,l5,i9)
!      end if

!      call out_dispersion_info( nfile(i) )

      if( ltddft .or. ltddft_fssh ) then
          call out_occupations( nfile(i), lspin )
!          call out_dish( nfile(i) )
      end if


      write(nfile(i),*) ' '
      write(nfile(i),3200)
      do it = 1, ntype
         write(nfile(i),3210) it, aname(nint(zatom(it))),  &
&         zatom(it), zv(it), watom(it)*avogad*welm*2.d0, covrad(it),  &
&         nhk(it), lmax(it), lclno(it),  &
&         lkbppi(it), lvandi(it), lvflag(it), llocli(it)
      end do
!      call out_lsomax( nfile(i), aname )
      write(nfile(i),3300)
      do it = 1, ntype
         write(nfile(i),3310) it, aname(nint(zatom(it))),  &
&         lking(it), rking(it), gkgmax(it), gkgexct(it),  &
&   llking(it), rlking(it), glkgmax(it), glkgexct(it), rctflc,  &
&         lpking(it), rpking(it), gpkgmax(it), gpkgexct(it),  &
&         lpcci(it), rpcc(it), rintchg(it),  &
&         icscale(it), vrandom(it)
      end do
!      if( lvand ) then
!          write(nfile(i),3500) mxref, rctmax
!          do it = 1, ntype
!             if( lvandi(it) )  &
!&                write(nfile(i),3510) it, aname(nint(zatom(it))),  &
!&                   ( nrefe(l,it), vdrcut(l,it), l = 0, lmax(it) )
!          end do
!      end if
!      if( lpaw ) then
!          write(nfile(i),'(/6x,a)')  &
!& 'frac_rcomp_it (compensation-charge cutoff = max(r_cl)/frac_rcomp_it)'
!          do it = 1, ntype
!             if( lvandi(it) )  &
!&                write(nfile(i),'(i2,a,a2,4f8.4)') it, ':', aname(nint(zatom(it))),  &
!&                   frac_rcomp_it(it)
!          end do
!      end if
!      if( lmulken ) then
!          write(nfile(i),'(/6x,a)') 'lmulao lmulbsis rxmulk rdecmp(1,l) rdecmp(2,l)'
!          do it = 1, ntype
!             write(nfile(i),3410) it, aname(nint(zatom(it))),  &
!&    ( lmulao(l,it), lmulbsis(l,it), rxmulk(l,it), rdecmp(1,l,it), rdecmp(2,l,it),  &
!&      l = 1, nmulao(it) )
!          end do
!          !---output parameters for bond charges
!          call out_bondc( nfile(i) )
!          !---ouput variables for band-decomposed overlap population
!          call out_nbdovp( nfile(i) )
!      end if
!      if( lsphexp ) then
!          write(nfile(i),'(/6x,a)') 'rsphexp'
!          do it = 1, ntype
!             write(nfile(i),3710) it, aname(nint(zatom(it))),  &
!&                                 rsphexp(it)
!          end do
!      end if
!      if( leda ) then
!          write(nfile(i),'(/6x,a)') 'radeda'
!          do it = 1, ntype
!             write(nfile(i),3710) it, aname(nint(zatom(it))),  &
!&                                 radeda(it)
!          end do
!      end if
!      if( lwannier .and. natwan > 0 ) then
!          write(nfile(i),'(/6x,a)') 'natwno, nwanaa'
!          do j = 1, natwan
!             write(nfile(i),3610) j, natwno(j), nwanaa(j)
!          end do
!      end if
!      if( lplusU ) then
!          write(nfile(i),'(/6x,a)') 'lplusU_at, plusU_l, plusU_U [eV], plusU_J [eV]'
!          do it = 1, ntype
!             write(nfile(i),'(i2,a,a2,l6,i9,f14.4,f14.4)')  &
!& it, ':', aname(nint(zatom(it))),  &
!& lplusU_at(it), plusU_l(it), plusU_U(it)/2.d0*hrdev, plusU_J(it)/2.d0*hrdev
!          end do
!      end if
!      if( lplusC ) then
!          write(nfile(i),'(/6x,a)') 'lplusC_at, ( plusC_r, plusC_e, l=0,lmax)'
!          do it = 1, ntype
!             write(nfile(i),'(i2,a,a2,l6,4(f8.2,f8.4))')  &
!& it, ':', aname(nint(zatom(it))),  &
!& lplusC_at(it), ( plusC_r(l,it), plusC_e(l,it), l = 0, lmax(it) )
!          end do
!      end if

      if( lspin ) then
          if( inispin > 1 ) then
              write(nfile(i),'(/6x,a)') ' iatspin    atmmagne'
              do j = 1, nspinat
                 write(nfile(i),'(i4,a,i6,f14.4)') j, ':', iatspin(j), atmmagne(j)
              end do
          end if
      end if

!      if( lnoncollinear ) then
!          if( inispin == 1 ) then
!              write(nfile(i),'(/6x,a)') ' atmmagne    umx     umy     umz'
!              write(nfile(i),'(f14.4,2x,3f8.4)') atmmagne(1), umx(1), umy(1), umz(1)
!          else if( inispin == 2 .or. inispin == 3 ) then
!              write(nfile(i),'(/6x,a)') ' iatspin         atmmagne    umx     umy     umz'
!              do j = 1, nspinat
!                 write(nfile(i),'(i4,a,i6,a,a2,a,f14.4,2x,3f8.4)')  &
!& j, ':', iatspin(j), ' (', aname(nint(zatom(iatspin(j)))), ')', atmmagne(j), umx(j), umy(j), umz(j)
!              end do
!          else if( inispin == 4 ) then
!              write(nfile(i),'(/6x,a)') ' # of atom  atmmagne    umx     umy     umz'
!              do j = 1, nspinat
!                 write(nfile(i),'(i4,a,i6,f14.4,2x,3f8.4)')  &
!& j, ':', iatspin(j), atmmagne(j), umx(j), umy(j), umz(j)
!              end do
!          end if
!      end if

      !---well potential
!      if( lwell ) then
!          call out_well_radius( nfile, myid, nodes, i,  &
!& ntype, aname, zatom )
!      end if

      write(nfile(i),*) ' '
      do it = 1, ntype
         write(nfile(i),*) ' it=', it
         write(nfile(i),*) ' nhk1, nhk2=', nhk1(it), nhk2(it)
         write(nfile(i),*) ' coordinates'
         do j = nhk1(it), nhk2(it)
            write(nfile(i),'(i5,3f12.8)') j, (ratm(ix,j), ix =1,3)
         end do
         write(nfile(i),*) ' velocities'
         do j = nhk1(it), nhk2(it)
            write(nfile(i),'(i5,3f12.8)') j, (vatm(ix,j), ix =1,3)
         end do
      end do

   end if myidif

end do

!
! 0:i6, 1:i10, 2:e11.4, 4:f10.5, 5:f10.2, 6:l6 descriptors
!
2000 format(1x, a12,i6,6x,    a12,i6,6x,    a12,i6)
2006 format(1x, a12,i6,6x,    a12,i6,6x,    a12,l6)
2020 format(1x, a12,i6,6x,    a12,es11.4,1x,a12,i6)
2602 format(1x, a12,l6,6x,    a12,i6,6x,    a12,es11.4)
2002 format(1x, a12,i6,6x,    a12,i6,6x,    a12,es11.4)
2044 format(1x, a12,i6,6x,    a12,f10.5,2x, a12,f10.5)
2042 format(1x, a12,i6,6x,    a12,f10.5,2x, a12,es11.4)
2410 format(1x, a12,f10.5,2x, a12,i10,2x,   a12,i6)
2050 format(1x, a12,i6,6x,    a12,f10.2,2x, a12,i6)
2600 format(1x, a12,l6,6x,    a12,i6,6x,    a12,i6)
2604 format(1x, a12,l6,6x,    a12,i6,6x,    a12,f10.5)
2606 format(1x, a12,l6,6x,    a12,i6,6x,    a12,l6)
2500 format(1x, a12,f10.2,2x, a12,i6,6x,    a12,i6)
2660 format(1x, a12,l6,6x,    a12,l6,6x,    a12,i6)
2666 format(1x, a12,l6,6x,    a12,l6,6x,    a12,l6)
2550 format(1x, a12,f10.2,2x, a12,f10.2,2x, a12,i6)
2555 format(1x, a12,f10.2,2x, a12,f10.2,2x, a12,f10.2)
2640 format(1x, a12,l6,6x,    a12,f10.5,2x, a12,i6)
2644 format(1x, a12,l6,6x,    a12,f10.5,2x, a12,f10.5)
2646 format(1x, a12,l6,6x,    a12,f10.5,2x, a12,l6)
2656 format(1x, a12,l6,6x,    a12,f10.2,2x, a12,l6)

3000 format(1x,a16,3es16.8)
3010 format(1x,a16,3l5) 
3020 format(1x,a16,3i5) 

3200 format(6x,'zatom   zv   watom  covrad  nhk lmax local',  &
&         ' lkbppi lvandi lvflag llocli' )
3300 format(/6x,'lking  rking  gkgmax  gkgexct'/  &
&       6x,'llking rlking glkgmax glkgexct     rctflc'/  &
&       6x,'lpking rpking gpkgmax gpkgexct lpcci rpcc  rint',  &
&         ' icscale vrandom' )
3500 format(/6x,'# of ref.E. & cutoff length for USPP for l=0,lmax :',  &
&          ' max. =', i2, f10.6)
3210 format(i2,':',a2,f6.1,f5.1,f9.4,f7.3,i5,i4,i5,4l7)
3310 format(i2,':',a2,l4,f9.3,2f8.3/  &
&             5x,l4,f9.3,2f8.3,   f13.2/  &
&             5x,l4,f9.3,2f8.3, l6,f7.2,f6.2,i5,l8 )
3410 format(i2,':',a2,i6,1x,l6,2x,f7.3,2x,f10.3,f12.3/  &
&           (5x,i6,1x,l6,2x,f7.3,(2x,f10.3,f12.3)) )
3510 format(i2,':',a2,4(i4,f10.6))
3610 format(i4,':',i6,i8)
3710 format(i2,':',a2,f7.3)


return
end subroutine




subroutine read_data( nfile, myid, nodes, lclust )
!-----------------------------------------------------------------------
!     read input variables from an input file
!-----------------------------------------------------------------------
use outfile
use param
use param_atom
!use symmop
implicit none
integer :: nfile(*), myid, nodes
logical :: lclust

!-----declare local variables
character(10) :: word
character(50) :: fname1
character(50) :: fname = 'control/filename'
integer :: iunit
integer :: ierror, istat
integer :: i, ix, j
logical :: lnotread_QMatoms = .true.
logical :: lnotread_supercell = .true.
logical :: lnotread_MDcell    = .true.
real*8,  dimension(3,3) :: hcell_
real*8  :: alpmatrix(3,3), betmatrix(3,3), sum


call allocate_unit_number( iunit )

!--- get file name : fname1
if(loutfile(1)) write(nfile(1),*) 'open file: ', fname(1:len_trim(fname))
if(loutfile(2)) write(nfile(2),*) 'open file: ', fname(1:len_trim(fname))
open( iunit, file=fname, status='old', action='read',  &
&     iostat=ierror )

!--- global maximum
call gimax(ierror)

blockif_a1: if( ierror == 0 ) then

   read(iunit,*,iostat=istat) fname1
   !----- error trap
   if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                   'error in file :'//fname )

else blockif_a1

   !----- error trap
   call fstop( nfile, myid, nodes, 'cannot open file :'//fname )

end if blockif_a1
close(iunit)


!--- open input file
if(loutfile(1)) write(nfile(1),*) 'open file: ', fname1(1:len_trim(fname1))
if(loutfile(2)) write(nfile(2),*) 'open file: ', fname1(1:len_trim(fname1))
open( iunit, file=fname1, status='old', action='read',  &
&     iostat=ierror )

!--- global maximum
call gimax(ierror)


blockif_b1: if( ierror == 0 ) then


   !-----get the numbers of species & atoms : ntype & natom
   prereaddo: do
      read(iunit,'(a10)',iostat=istat) word
      preblockif_b2: if( istat == 0 ) then

         prereadif: if( word == '*atoms    ' ) then

            !-----get the numbers of species & atoms : ntype & natom
            call read_atom_number( nfile, myid, nodes, iunit,  &
& ntype, natom )

            lnotread_QMatoms = .false.

            !----- exit
            exit

         end if prereadif

      else if( istat < 0 ) then preblockif_b2

         !----- exit, if at EOF
         exit

      else preblockif_b2

         !----- error trap
         call fstop( nfile, myid, nodes,'error in file :'//fname1)

      end if preblockif_b2

   end do prereaddo


   !----- error trap
   if( lnotread_QMatoms ) call fstop( nfile, myid, nodes,  &
&                         'no QM atoms in file :'//fname1)

   natom_alloc = natom
   natnod_alloc = (natom+nodes)/nodes


   !-----allocate memory for atoms
   call param_atom_alloc( nfile, myid, nodes,  &
& alloc_mem, ntype, natom, natom_alloc )


   !----- rewind file
   rewind(iunit)


!   kreaddo: do
!      read(iunit,'(a10)',iostat=istat) word
!      kblockif_b2: if( istat == 0 ) then
!
!         kreadif: if( word == '*k-points ' ) then
!
!            !-----get the numbers of k points  : nkpnt
!            call read_k_number( nfile, myid, nodes, iunit,  &
!& nkpnt )
!
!            !----- exit
!            exit
!
!         end if kreadif
!
!      else if( istat < 0 ) then kblockif_b2
!
!         !----- exit, if at EOF
!         exit
!
!      else kblockif_b2
!
!         !----- error trap
!         call fstop( nfile, myid, nodes,'error in file :'//fname1)
!
!      end if kblockif_b2
!
!   end do kreaddo

   if( nkpnt > 0 ) then
       !----- allocate memory for k points
       allocate( bzk(3,nkpnt), wbzk(nkpnt), lgammak(nkpnt),  &
!& kpsymop(0:nofsym,nkpnt), ibzk(3,1,nkpnt),  &
& stat=istat )

       !----- set default values
       bzk(1,1) = 0.d0
       bzk(2,1) = 0.d0
       bzk(3,1) = 0.d0
       wbzk(1) = 1.d0
       lgammak(1) = .true.
!       kpsymop(0,1) = 1
!       kpsymop(1,1) = 0
!       ibzk(:,:,:) = 0
   end if

   !----- rewind file
   rewind(iunit)


   !----- read input data
   readdo: do
      read(iunit,'(a10)',iostat=istat) word
      blockif_b2: if( istat == 0 ) then

         readif: if( word == '*start(on/' .or. word == '*start    ' ) then

            call read_start( nfile, myid, nodes, iunit, lstart )

         else if( word == '*TDDFT-MD ' ) then readif

            call read_tddft( nfile, myid, nodes, iunit,  &
& ltddft, ltddft_fssh, ltddft_start, dttddft, lfssh_switch, &
& lfssh_gsscf, lfssh_gsfrc, imxchg_fssh, aslh_fssh, bslh_fssh, tdbroad,  &
& lfssh_random, rseed_fssh, &
& lfssh_boltzmn, lfssh_vscale, tminimum, lfssh_parallel, ltddft_nscforce )

!         else if( word == '*linear-re' ) then readif
!
!            call read_lrtddft( nfile, myid, nodes, iunit,  &
!& lrtddft, fdmeshxc, lrspecific, enediff, threshrpa, threshdipole,  &
!& threshexcite, llcexchange, alpha_llcexchange, large_cell_llcexchange,  &
!& hybrid_mixing, foersigma, ldiagonal,  &
!& lscissors_lrtddft, ecorr_tirplet, ecorr_singlet, lfission_rate,  &
!& lemission_rate, refractive, ltda, lcic_tda, cic_tda_coupling )

!         else if( word == '*PAW      ' ) then readif
!
!            call read_paw( nfile, myid, nodes, iunit,  &
!& lpaw, lpaw_sym, paw_mix )

!         else if( word == '*Relativis' ) then readif
!
!            call read_rela( nfile, myid, nodes, iunit,  &
!& lrela, lrela_full, lrela_sop )

!         else if( word == '*spin-orbi' ) then readif
!
!            call read_soc( nfile, myid, nodes, iunit,  &
!& lsoc, lsoc_full )

!         else if( word == '*cluster  ' ) then readif

!            call read_clust( nfile, myid, nodes, iunit, lclust )

         else if( word == '*approxima' ) then readif

            call read_approx( nfile, myid, nodes, iunit,  &
& jgga, jgga_org, lvdw, lvdw_factorize, ldftd, lplusU, plusU_mix, lplusC,  &
& jhybrid, hmixing, hrange, mgridred )

!         else if( word == '*scissors ' ) then readif

!            call read_scissors( nfile, myid, nodes, iunit,  &
!& lscissors, ecorr_acceptor, ecorr_donor, zboundary )

!         else if( word == '*k-points ' ) then readif

!             if( nkpnt > 0 ) then
!                call read_kp( nfile, myid, nodes, iunit,  &
!& nkpnt, bzk, wbzk, lgammak, lgamma )
!             else
!                call read_kp2( nfile, myid, nodes, iunit,  &
!& npkx, ltetra, lksgamma, lgamma )
!             end if

         else if( word == '*SCF itera' ) then readif

            call read_scf( nfile, myid, nodes, iunit,  &
& iscfmx, itrial, tolpot, tolres, tolHF_KS, imxchg, aslh, bslh, itratn,  &
& louteig, loutenergy, loutzansa, loutlog,  &
& imxspin, amxspin, bmxspin, nmxspin, lhcunt, nstabi, xstabi,  &
& imxchg_trial, aslh_trial, bslh_trial, imxspin_trial, amxspin_trial, bmxspin_trial,  &
& cmetric, spinmetric, wkerker, lmixreal )

         else if( word == '*Kohn-Sham' ) then readif

            call read_ks( nfile, myid, nodes, iunit,  &
& ihldam, toleig, roughfeig, roughfzan, threinn, itermx, iteremx, kstrial, methodcg )

!         else if( word == '*Poisson e' ) then readif

!            call read_pois( nfile, myid, nodes, iunit,  &
!& multg, tolcg, weigrd, nd2v, msrhmx, msrhmy, msrhmz )

         else if( word == '*molecular' ) then readif

            call read_md( nfile, myid, nodes, iunit,  &
& ifmd, dtmd, nstep_ini, treq, liscale, iscnum, iscstp,  &
& ichest, ihest, laspc_chg, laspc_wv, laspc, nskp_aspc, aspc_corr, aslh_aspc, bslh_aspc,  &
& iscfmx_aspc, tolres_aspc, pwscale, lmomzero, ioptmze, ioptmze_cell, nshockv,  &
& lxlbomd, kxlbomd, irstrct, irstrct_sub, lcell_rstrct )

         else if( word == '*save data' ) then readif

            call read_save( nfile, myid, nodes, iunit,  &
& lsave, lsreal8 )

!         else if( word == '*Mulliken ' ) then readif
!
!            call read_mul( nfile, myid, nodes, iunit,  &
!& lmulken, nskpmulk, ldecmpovp, lspdatom, lpmchgat, lmdecmpdos )

!         else if( word == '*Spherical' ) then readif
!
!            call read_sphexp( nfile, myid, nodes, iunit,  &
!& lsphexp, nskpsphexp, radsphexp )

!         else if( word == '*EDA      ' ) then readif
!
!            call read_eda( nfile, myid, nodes, iunit,  &
!& leda, nskpeda, radgeda )!, rdelaunay )

!         else if( word == '*Wannier f' ) then readif
!
!            call read_wan( nfile, myid, nodes, iunit,  &
!& lwannier, iwbnd1, iwbnd2, iterwan, tolwan, iwstt1, iwstt2,  &
!& nskpwann, loutuni )

!         else if( word == '*Conductiv' ) then readif
!
!            call read_conduct( nfile, myid, nodes, iunit,  &
!& lconduct, ldcconduct, lacconduct, nskpconduct, wgconduct, tempconduct,  &
!& efconduct, freqacmx, icband1, icband2, ichband1, ichband2,  &
!& iceband1, iceband2, immtband1, immtband2 )

         else if( word == '*stress ca' ) then readif

            call read_str( nfile, myid, nodes, iunit,  &
& lstress, nskip_stress )

!         else if( word == '*atomic ch' ) then readif
!
!            call read_achg( nfile, myid, nodes, iunit,  &
!& lintchg, nskip_intchg )

!         else if( word == '*dump char' ) then readif
!
!            call read_dpchg( nfile, myid, nodes, iunit,  &
!& ldpchg, nskip_dpchg, dc_xmin, dc_xmax, dc_ymin, dc_ymax, dc_zmin, dc_zmax,  &
!& lcompress_dpchg, ndigit_dpchg )

         else if( word == '*dump wave' ) then readif

            call read_dpwav( nfile, myid, nodes, iunit,  &
& ldpwav, nskip_dpwav, ibstt1, ibstt2, wv_xmin, wv_xmax, wv_ymin, wv_ymax, wv_zmin, wv_zmax,  &
& lcompress_dpwav, ndigit_dpwav )

!         else if( word == '*dump loca' ) then readif
!
!            call read_dppot( nfile, myid, nodes, iunit,  &
!& ldppot, nskip_dppot, nav_dppot, pt_xmin, pt_xmax, pt_ymin, pt_ymax, pt_zmin, pt_zmax,  &
!& lcompress_dppot, ndigit_dppot )

!         else if( word == '*hot elect' ) then readif
!
!            call read_hoteh( nfile, myid, nodes, iunit,  &
!& lhoteh, nskip_hoteh, ibhoth,ibhote )

         else if( word == '*MD cell  ' ) then readif

            call read_cell( nfile, myid, nodes, iunit, h_MD )

            lnotread_MDcell = .false.

         else if( word == '*supercell' ) then readif

            call read_cell( nfile, myid, nodes, iunit, hcell )

            lnotread_supercell = .false.

!         else if( word == '*vacuum   ' ) then readif
!
!            call read_vac( nfile, myid, nodes, iunit,  &
!& lvacuum, vacuum, alpha_ldouble_grid_recip )

!         else if( word == '*spherical' ) then readif

!            call read_sphere( nfile, myid, nodes, iunit, lsphere )

         else if( word == '*planewave' ) then readif

            call read_pw( nfile, myid, nodes, iunit,  &
& ecut, ecutdens, ecutsoft, ecutorth, ecutc )

!         else if( word == '*double-gr' ) then readif

!            call read_dg( nfile, myid, nodes, iunit,  &
!& iosp, ecutlong, ecutlong_small )

         else if( word == '*electroni' ) then readif

            call read_band( nfile, myid, nodes, iunit,  &
& noband, nband, lfermi, tfermi, ncion )

         else if( word == '*spin pola' ) then readif

            call read_lsda( nfile, myid, nodes, iunit,  &
& lspin, lnoncollinear, lfixud, diffud, lwfrand, lduplicate,  &
& latmref, lclearamag, refmx, refmy, refmz,  &
& lclearmagne_in_x, lclearmagne_in_y, lclearmagne_in_z,  &
& lfixmagne_in_x, lfixmagne_in_y, lfixmagne_in_z,  &
& fixmagne_in_x, fixmagne_in_y, fixmagne_in_z )

         else if( word == '*atoms    ' ) then readif

            call read_atom( nfile, myid, nodes, iunit,  &
& ntype, natom, lplcom, zatom, lkbppi, lvandi, lvflag, llocli,  &
& lkbpp, lvand, lvfull, llocl, lking, rking, gkgmax, gkgexct,  &
& llking, rlking, glkgmax, glkgexct,  &
& lpcci, rpcc, lpking, rpking, gpkgmax, gpkgexct, nhk,  &
& icscale, vrandom, ratm, vatm,  &
& rdecmp, rxmulk, lmulao, lmulbsis, nmulao, mulaox, rsphexp, radeda,  &
& rintchg, frac_rcomp_it, lplusU_at, plusU_l, plusU_U, plusU_J,  &
& lplusC_at, plusC_r, plusC_e, mxlx )
!
!               else if( word == '*local pot' ) then readif
!
!                  call read_local( nfile, myid, nodes, iunit, rctflc )

!         else if( word == '*Cholesky ' ) then readif

!            call read_chol( nfile, myid, nodes, iunit, node_c )

!         else if( word == '*eigenvalu' ) then readif

!            call read_eig( nfile, myid, nodes, iunit, node_r )

!         else if( word == '*work arra' ) then readif

!            call read_wrk( nfile, myid, nodes, iunit, dblock )

!         else if( word == '*table dim' ) then readif

!            call read_tabd( nfile, myid, nodes, iunit,  &
!& mx1, mx1loc )

!         else if( word == '*symmetry ' ) then readif

!            if( .not.lreadsym ) then
!                call read_symm( nfile, myid, nodes, iunit,  &
!& lsymop, nofsym, nmsyop, symrtn, symtrs, linvsn,  &
!& symm, nsymm, trsfrm, ltrsfm, trsfri, ltrsfi, radsym, lsymcd, lsymspin,  &
!& ddmin, ddmax )
!
!                call chksymm( nfile, myid, nodes, ierror )
!                !----- error trap
!                if( ierror /= 0 )  &
!&          call fstop( nfile, myid, nodes,'error in symmetry op.')
!
!                lreadsym = .true.
!            end if

!         else if( word == '*well pote' ) then readif
!
!            call read_well( nfile, myid, nodes, iunit,  &
!& ntype, lwell, wellheight )

!         else if( word == '*electric ' ) then readif
!
!            call read_efield( nfile, myid, nodes, iunit,  &
!& lefield, efield, lsawtooth, lsawtooth_shape, loutpolarization,  &
!& lefield_start, lconstraintD )

         end if readif

      else if( istat < 0 ) then blockif_b2

         !----- exit, if at EOF
         exit

      else blockif_b2

         !----- error trap
         call fstop( nfile, myid, nodes,'error in file :'//fname1)

      end if blockif_b2

   end do readdo

else blockif_b1

   !----- error trap
   call fstop( nfile, myid, nodes, 'cannot open file :'//fname1 )

end if blockif_b1
close(iunit)

call deallocate_unit_number( iunit )



!-----set MD/super-cell vectors
!----- error trap
if( lnotread_supercell .and. lnotread_MDcell ) then
    call fstop( nfile, myid, nodes,  &
& 'error in setting MD/super-cell vectors' )
end if

!-----check cell vectors for MSST
if( ifmd == 10 ) then
    irstrct = 0
    !---store original supercell
    if( lnotread_MDcell ) then
        hcell_rec(:,:) = hcell(:,:)
    else
        hcell_rec(:,:) = h_MD(:,:)
    end if
end if

lvshape = ifmd == 4 .or. ifmd == 10 .or. ioptmze_cell >= 0
!-----read previous hcell, if lvshape = .true.
if( lstart .and. lvshape ) then
    call rdhcell( nfile, hcell_, lvshape )
    !-----check hcell_
    if( abs(hcell_(1,1)) > 1.d-13 .or.  &
&       abs(hcell_(1,2)) > 1.d-13 .or.  &
&       abs(hcell_(1,3)) > 1.d-13 .or.  &
&       abs(hcell_(2,1)) > 1.d-13 .or.  &
&       abs(hcell_(2,2)) > 1.d-13 .or.  &
&       abs(hcell_(2,3)) > 1.d-13 .or.  &
&       abs(hcell_(3,1)) > 1.d-13 .or.  &
&       abs(hcell_(3,2)) > 1.d-13 .or.  &
&       abs(hcell_(3,3)) > 1.d-13       ) then
        do i = 1, 3
        do ix = 1, 3
           hcell(ix,i) = hcell_(ix,i)
        end do
        end do
        lnotread_supercell = .false.
    end if
end if

!if( irstrct /= 0 .and. irstrct /= 10 ) then
!    if( lnotread_MDcell ) then
!        !-----restriction for MD cell
!        call restrictMDcell( irstrct, irstrct_sub, hcell, lcell_rstrct )
!    else
!        !-----restriction for MD cell
!        call restrictMDcell( irstrct, irstrct_sub, h_MD, lcell_rstrct )
!    end if
!end if

if( lnotread_MDcell ) then
    do i = 1, 3
    do ix = 1, 3
       h_MD(ix,i) = hcell(ix,i)
    end do
    end do
end if

if( lnotread_supercell ) then
    do i = 1, 3
    do ix = 1, 3
       hcell(ix,i) = h_MD(ix,i)
    end do
    end do
end if


!-----check cell vectors for MSST
if( ifmd == 10 ) then
    call msst_setmatrix( nfile, myid, nodes,  &
& nshockv, hcell, alpmatrix, betmatrix, trsmatrix, ierror )

    if( ierror /= 0 )  &
& call fstop( nfile, myid, nodes,'error in msst_setmatrix')

    alpmatrix(:,:) = h_MD(:,:)
    do i = 1, 3
    do j = 1, 3
       sum = 0.d0
       do ix = 1, 3
          sum = sum + trsmatrix(i,ix)*alpmatrix(ix,j)
       end do
       h_MD(i,j) = sum
    end do
    end do
end if


!if( .not.lgamma ) then
!    !---set k vectors
!    if( nkpnt == 0 ) then
!
!        do i = 1, 2
!        if( loutfile(i) ) then
!            write(nfile(i),*) ' '
!            if( lksgamma ) then
!                write(nfile(i),*) '*** k-point generation centered at the Gamma point'
!            else
!                write(nfile(i),*) '*** k-point generation by MonkhostPack method'
!            end if
!        end if
!        end do
!
!        call MonkhostPack( nfile, myid, nodes, npkx, hcell, ltetra, lksgamma )
!
!        !--- check if GAMMA point
!        do i = 1, nkpnt
!           lgammak(i) = abs(bzk(1,i)) < 1.d-08  &
!&                 .and. abs(bzk(2,i)) < 1.d-08  &
!&                 .and. abs(bzk(3,i)) < 1.d-08
!        end do
!
!        lgamma = nkpnt == 1 .and. lgammak(1)
!
!    else
!
!        do i = 1, nkpnt
!           kpsymop(0,i) = 1
!           kpsymop(1,i) = 0
!        end do
!
!    end if
!end if


return
end subroutine




subroutine read_start( nfile, myid, nodes, iunit, lstart )
!-----------------------------------------------------------------------
! read variables for restarting the program
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
logical :: lstart
character(10) :: word
integer :: istat


readdo: do
   read(iunit,'(a10)',iostat=istat) word

   readif: if( istat == 0 ) then

      selectif: if( word == '(on/off)  ' .or. word == '(how of it' ) then

         read(iunit,*,iostat=istat) lstart
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0001 in read_start' )

      else if( word == '*end      ' ) then selectif

         exit

      endif selectif

   else readif

      !----- error trap
      call fstop( nfile, myid, nodes, 'error-0000 in read_start' )

   end if readif
end do readdo


return
end subroutine




subroutine read_tddft( nfile, myid, nodes, &
& iunit, ltddft, ltddft_fssh, ltddft_start, dttddft, lfssh_switch, &
& lfssh_gsscf, lfssh_gsfrc, imxchg_fssh, aslh_fssh, bslh_fssh, tdbroad,  &
& lfssh_random, rseed_fssh, &
& lfssh_boltzmn, lfssh_vscale, tminimum, lfssh_parallel, ltddft_nscforce )
!-----------------------------------------------------------------------
!     set ltddft : .true. = execute MD based on TDDFT
!                  .false.= execute MD based on conventional DFT
!   ltddft_start : .true. = restart
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
logical :: ltddft, ltddft_fssh, ltddft_start, lfssh_switch, lfssh_gsscf, lfssh_gsfrc
logical :: lfssh_random
integer :: imxchg_fssh
real*8  :: dttddft, aslh_fssh, bslh_fssh, tdbroad, rseed_fssh
logical :: lfssh_boltzmn, lfssh_vscale, lfssh_parallel, ltddft_nscforce
real*8  :: tminimum
character(10) :: word
integer :: istat, nocc_change, nexciton


readdo: do
   read(iunit,'(a10)',iostat=istat) word

   readif: if( istat == 0 ) then

      selectif: if( word == '(on/off)  ' .or. word == '(how of it' ) then

         read(iunit,*,iostat=istat) ltddft
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0001 in read_tddft' )

      else if( word == '(restart) ' ) then selectif

         read(iunit,*,iostat=istat) ltddft_start
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0002 in read_tddft' )

      else if( word == '(occupatio' ) then selectif

         read(iunit,*,iostat=istat) nocc_change
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0003 in read_tddft' )

         call set_occupations( nfile, myid, nodes, iunit, nocc_change )

!      else if( word == '(initial e' ) then selectif

!         read(iunit,*,iostat=istat) nexciton
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-1003 in read_tddft' )
!
!         call set_exciton( nfile, myid, nodes, iunit, nexciton )

      else if( word == '(FSSH)    ' ) then selectif

         read(iunit,*,iostat=istat) ltddft_fssh
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0004 in read_tddft' )

      else if( word == '(time step' ) then selectif

         read(iunit,*,iostat=istat) dttddft
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0005 in read_tddft' )

         dttddft = dttddft / 2.d0   ! in [Rydberg units]

      else if( word == '(FSSH-swit' ) then selectif

         read(iunit,*,iostat=istat) lfssh_switch
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0006 in read_tddft' )

      else if( word == '(FSSH-grou' ) then selectif

         read(iunit,*,iostat=istat) lfssh_gsscf
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0007 in read_tddft' )

      else if( word == '(FSSH-mixi' ) then selectif

         read(iunit,*,iostat=istat) aslh_fssh, bslh_fssh
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0008 in read_tddft' )

      else if( word == '(FSSH-char' ) then selectif

         read(iunit,*,iostat=istat) imxchg_fssh
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0018 in read_tddft' )

         read(iunit,*,iostat=istat) aslh_fssh, bslh_fssh
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0028 in read_tddft' )

      else if( word == '(broadenin' ) then selectif

         read(iunit,*,iostat=istat) tdbroad
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0009 in read_tddft' )

      else if( word == '(FSSH-rand' ) then selectif

         read(iunit,*,iostat=istat) lfssh_random
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0010 in read_tddft' )

         read(iunit,*,iostat=istat) rseed_fssh
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0011 in read_tddft' )

      else if( word == '(Boltzmann' ) then selectif

         read(iunit,*,iostat=istat) lfssh_boltzmn
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0012 in read_tddft' )

      else if( word == '(velocity ' ) then selectif

         read(iunit,*,iostat=istat) lfssh_vscale
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0013 in read_tddft' )

         read(iunit,*,iostat=istat) tminimum
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0014 in read_tddft' )

      else if( word == '(parallel ' ) then selectif

         read(iunit,*,iostat=istat) lfssh_parallel
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0015 in read_tddft' )

      else if( word == '(NSC force' ) then selectif

         read(iunit,*,iostat=istat) ltddft_nscforce
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0016 in read_tddft' )

!      else if( word == '(DISH)    ' ) then selectif

!         call set_dish( nfile, myid, nodes, iunit )

      else if( word == '(ground st' ) then selectif

         read(iunit,*,iostat=istat) lfssh_gsfrc
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0017 in read_tddft' )

      else if( word == '*end      ' ) then selectif

         exit

      endif selectif

   else readif

      !----- error trap
      call fstop( nfile, myid, nodes, 'error-0000 in read_tddft' )

   end if readif
end do readdo


return
end subroutine




subroutine read_approx( nfile, myid, nodes, iunit,  &
& jgga, jgga_org, lvdw, lvdw_factorize, ldftd, lplusU, plusU_mix, lplusC,  &
& jhybrid, hmixing, hrange, mgridred )
!-----------------------------------------------------------------------
! read variables for approximation for Exc
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
integer :: jgga, jgga_org
logical :: lvdw, lvdw_factorize
logical :: ldftd, lplusU, lplusC
real*8  :: plusU_mix
integer :: jhybrid
real*8  :: hmixing, hrange
integer :: mgridred

!---Declare local variables
character(10) :: word
integer :: istat
integer :: imod1, imod2, iosp
real*8  :: rvdwmax
real*8  :: q_cut, dlambda, d_soft, phi_at_0
integer :: nqmesh


readdo: do
   read(iunit,'(a10)',iostat=istat) word

   readif: if( istat == 0 ) then

      selectif: if( word == '(approxima' ) then

         read(iunit,*,iostat=istat) jgga
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
                          'error-0001 in read_approx' )

      else if( word == '(cutoff le' ) then selectif

!         read(iunit,*,iostat=istat) rvdwmax
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!                          'error-0002 in read_approx' )
!
!         call set_rvdwmax( rvdwmax )

      else if( word == '(integrati' ) then selectif

!         read(iunit,*,iostat=istat) imod1
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!                          'error-0003 in read_approx' )
!
!         read(iunit,*,iostat=istat) imod2
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!                          'error-0004 in read_approx' )
!
!         call set_imod( imod1, imod2 )

      else if( word == '(order of ' ) then selectif

!         read(iunit,*,iostat=istat) iosp
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!                          'error-0005 in read_approx' )
!
!         call set_iosp( iosp )

      else if( word == '(DFT-D)   ' ) then selectif

         read(iunit,*,iostat=istat) ldftd
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
                          'error-0006 in read_approx' )

!      else if( word == '(DFT+U)   ' ) then selectif
!
!         read(iunit,*,iostat=istat) lplusU
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!                          'error-0007 in read_approx' )
!
!      else if( word == '(onsite ch' ) then selectif
!
!         read(iunit,*,iostat=istat) plusU_mix
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!                          'error-0007-2 in read_approx' )
!
!      else if( word == '(DFT+C)   ' ) then selectif
!
!         read(iunit,*,iostat=istat) lplusC
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!                          'error-0008 in read_approx' )

      else if( word == '(kernel fa' ) then selectif

!         read(iunit,*,iostat=istat) lvdw_factorize
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!                          'error-0009 in read_approx' )
!
!         call set_lvdw_factorize( lvdw_factorize )

      else if( word == '(q-mesh ra' ) then selectif

!         read(iunit,*,iostat=istat) q_cut, dlambda, nqmesh
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!                          'error-0010 in read_approx' )
!
!         call set_qrange( q_cut, dlambda, nqmesh )

      else if( word == '(kernel cu' ) then selectif

!         read(iunit,*,iostat=istat) d_soft, phi_at_0
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!                          'error-0011 in read_approx' )
!
!         call set_d_soft( d_soft, phi_at_0 )

!      else if( word == '(hybrid fu' ) then selectif
!
!         read(iunit,*,iostat=istat) jhybrid
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!                          'error-0020 in read_approx' )
!
!      else if( word == '(hybrid mi' ) then selectif
!
!         read(iunit,*,iostat=istat) hmixing
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!                          'error-0021 in read_approx' )
!
!      else if( word == '(range-sep' ) then selectif
!
!         read(iunit,*,iostat=istat) hrange
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!                          'error-0022 in read_approx' )
!
!      else if( word == '(grid-redu' ) then selectif
!
!         read(iunit,*,iostat=istat) mgridred
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!                          'error-0023 in read_approx' )

      else if( word == '*end      ' ) then selectif

         exit

      endif selectif

   else readif

      !----- error trap
      call fstop(nfile, myid, nodes, 'error-0000 in read_approx')

   end if readif
end do readdo


jgga_org = jgga
!if( jgga == 4 .or. jgga == 5 .or. jgga == 6 ) then
!    if( jgga == 4 .or. jgga == 5 ) then
!        call set_param_revPBE
!    end if
!    if( jgga == 6 ) then
!        call set_param_vdW_DF2
!    end if
!    lvdw = jgga == 5 .or. jgga == 6
!    jgga = 2
!end if
!if( lvdw .and. .not.lvdw_factorize ) call set_chkftm_even


return
end subroutine




subroutine read_scf( nfile, myid, nodes, iunit,  &
& iscfmx, itrial, tolpot, tolres, tolHF_KS, imxchg, aslh, bslh, itratn,  &
& louteig, loutenergy, loutzansa, loutlog,  &
& imxspin, amxspin, bmxspin, nmxspin, lhcunt, nstabi, xstabi,  &
& imxchg_trial, aslh_trial, bslh_trial, imxspin_trial, amxspin_trial, bmxspin_trial,  &
& cmetric, spinmetric, wkerker, lmixreal )
!-----------------------------------------------------------------------
! read variables for SCF iteration
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
integer :: iscfmx, itrial
real*8  :: tolpot, tolres, tolHF_KS
integer :: imxchg
real*8  :: aslh, bslh
integer :: itratn
logical :: louteig, loutenergy, loutzansa, loutlog
integer :: imxspin
real*8  :: amxspin, bmxspin
integer :: nmxspin
logical :: lhcunt
integer :: nstabi
real*8  :: xstabi
integer :: imxchg_trial
real*8  :: aslh_trial, bslh_trial
integer :: imxspin_trial
real*8  :: amxspin_trial, bmxspin_trial
real*8  :: cmetric, spinmetric, wkerker
logical :: lmixreal

character(10) :: word
integer :: istat


readdo: do
   read(iunit,'(a10)',iostat=istat) word

   readif: if( istat == 0 ) then

      selectif: if( word == '(global it' ) then

         read(iunit,*,iostat=istat) iscfmx
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0001 in read_scf' )

      else if( word.eq.'(trial glo' ) then selectif

         read(iunit,*,iostat=istat) itrial
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0002 in read_scf' )

      else if( word.eq.'(tolerance' ) then selectif

         read(iunit,*,iostat=istat) tolpot
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0003 in read_scf' )

         read(iunit,*,iostat=istat) tolres
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0003-2 in read_scf' )

      else if( word.eq.'(HF-KS tol' ) then selectif

         read(iunit,*,iostat=istat) tolHF_KS
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0004 in read_scf' )

      else if( word.eq.'(charge mi' ) then selectif

         read(iunit,*,iostat=istat) imxchg
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0005 in read_scf' )

         read(iunit,*,iostat=istat) aslh, bslh
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0006 in read_scf' )

      else if( word.eq.'(HC produc' ) then selectif

         read(iunit,*,iostat=istat) lhcunt
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0007 in read_scf' )

      else if( word.eq.'(number of' ) then selectif

         read(iunit,*,iostat=istat) itratn
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0008 in read_scf' )

      else if( word.eq.'(spin-dens' ) then selectif

         read(iunit,*,iostat=istat) imxspin
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0009 in read_scf' )

         read(iunit,*,iostat=istat) amxspin, bmxspin
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0010 in read_scf' )

      else if( word.eq.'(# of spin' ) then selectif

         read(iunit,*,iostat=istat) nmxspin
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0011 in read_scf' )

      else if( word.eq.'(output ei' ) then selectif

         read(iunit,*,iostat=istat) louteig
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0012 in read_scf' )

      else if( word.eq.'(output en' ) then selectif

         read(iunit,*,iostat=istat) loutenergy
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0013 in read_scf' )

      else if( word.eq.'(output re' ) then selectif

         read(iunit,*,iostat=istat) loutzansa
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0013-2 in read_scf' )

      else if( word.eq.'(output de' ) then selectif

         read(iunit,*,iostat=istat) loutlog
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0013-3 in read_scf' )

      else if( word.eq.'(convergen' ) then selectif

         read(iunit,*,iostat=istat) nstabi
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0014 in read_scf' )

         read(iunit,*,iostat=istat) xstabi
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0015 in read_scf' )

      else if( word.eq.'(trial cha' ) then selectif

         read(iunit,*,iostat=istat) imxchg_trial
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0016 in read_scf' )

         read(iunit,*,iostat=istat) aslh_trial, bslh_trial
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0017 in read_scf' )

      else if( word.eq.'(trial spi' ) then selectif

         read(iunit,*,iostat=istat) imxspin_trial
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0018 in read_scf' )

         read(iunit,*,iostat=istat) amxspin_trial, bmxspin_trial
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0019 in read_scf' )

      else if( word.eq.'(metric in' ) then selectif

         read(iunit,*,iostat=istat) cmetric
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0020 in read_scf' )

      else if( word.eq.'(spin metr' ) then selectif

         read(iunit,*,iostat=istat) spinmetric
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0021 in read_scf' )

      else if( word.eq.'(Kerker-mi' ) then selectif

         read(iunit,*,iostat=istat) wkerker
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0022 in read_scf' )

      else if( word.eq.'(magnetic ' ) then selectif

         read(iunit,*,iostat=istat) lmixreal
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0023 in read_scf' )

      else if( word == '*end      ' ) then selectif

         exit

      endif selectif

   else readif

      !----- error trap
      call fstop(nfile, myid, nodes, 'error-0000 in read_scf')

   end if readif
end do readdo


return
end subroutine




subroutine read_ks( nfile, myid, nodes, iunit,  &
& ihldam, toleig, roughfeig, roughfzan, threinn, itermx, iteremx, kstrial, methodcg )
!-----------------------------------------------------------------------
! read variables for Kohn-Sham equation
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
integer :: ihldam
real*8  :: toleig, roughfeig, roughfzan, threinn
integer :: itermx, iteremx, kstrial, methodcg
character(10) :: word
integer :: istat


readdo: do
   read(iunit,'(a10)',iostat=istat) word

   readif: if( istat == 0 ) then

      selectif: if( word == '(method)  ' .or. word == '(how of it' ) then

         read(iunit,*,iostat=istat) ihldam
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0001 in read_ks' )

      else if( word.eq.'(tolerance' ) then selectif

         read(iunit,*,iostat=istat) toleig
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0002 in read_ks' )

      else if( word.eq.'(iteration' ) then selectif

         read(iunit,*,iostat=istat) itermx
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0003 in read_ks' )

      else if( word.eq.'(empty-ban' ) then selectif

         read(iunit,*,iostat=istat) iteremx
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0004 in read_ks' )

      else if( word.eq.'(threshold' ) then selectif

         read(iunit,*,iostat=istat) threinn
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0005 in read_ks' )

      else if( word.eq.'(trial glo' ) then selectif

         read(iunit,*,iostat=istat) kstrial
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0006 in read_ks' )

      else if( word.eq.'(CG method' ) then selectif

         read(iunit,*,iostat=istat) methodcg
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0007 in read_ks' )

      else if( word.eq.'(rough tol' ) then selectif

         read(iunit,*,iostat=istat) roughfeig
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0008 in read_ks' )

         read(iunit,*,iostat=istat) roughfzan
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0009 in read_ks' )

      else if( word == '*end      ' ) then selectif

         exit

      endif selectif

   else readif

      !----- error trap
      call fstop( nfile, myid, nodes, 'error-0000 in read_ks' )

   end if readif
end do readdo


return
end subroutine




subroutine read_md( nfile, myid, nodes, iunit,  &
& ifmd, dtmd, nstep_ini, treq, liscale, iscnum, iscstp,  &
& ichest, ihest, laspc_chg, laspc_wv, laspc, nskp_aspc, aspc_corr, aslh_aspc, bslh_aspc,  &
& iscfmx_aspc, tolres_aspc, pwscale, lmomzero, ioptmze, ioptmze_cell, nshockv,  &
& lxlbomd, kxlbomd, irstrct, irstrct_sub, lcell_rstrct )
!-----------------------------------------------------------------------
! read variables for molecular dynamics
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
integer :: ifmd
real*8  :: dtmd
!---      integer :: nstop
integer :: nstep_ini
real*8  :: treq
logical :: liscale, lmomzero
integer :: iscnum, iscstp
integer :: ichest, ihest
logical :: laspc_chg, laspc_wv, laspc
integer :: nskp_aspc, iscfmx_aspc
real*8  :: tolres_aspc, aspc_corr, aslh_aspc, bslh_aspc
real*8  :: pwscale
integer :: ioptmze, ioptmze_cell
integer :: nshockv(3)
logical :: lxlbomd
integer :: kxlbomd
integer :: irstrct, irstrct_sub
logical :: lcell_rstrct(3)

!-----declare local variables
real*8  :: dummy
character(10) :: word
integer :: i, istat


!-----set default value
call set_ioptmze( ioptmze, ioptmze_cell )

readdo: do
   read(iunit,'(a10)',iostat=istat) word

   readif: if( istat == 0 ) then

      selectif: if( word == '(method)  ' .or. word == '(how of it' ) then

         read(iunit,*,iostat=istat) ifmd
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0001 in read_md' )

      else if( word == '(time step' ) then selectif

         read(iunit,*,iostat=istat) dtmd   !---, nstop
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0002 in read_md' )

         dtmd = dtmd / 2.d0   ! in [Rydberg units]

      else if( word == '(temperatu' ) then selectif

         read(iunit,*,iostat=istat) treq
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0003 in read_md' )

      else if( word == '(check tem' ) then selectif

         read(iunit,*,iostat=istat) liscale
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0004 in read_md' )

         read(iunit,*,iostat=istat) iscnum
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0005 in read_md' )

         read(iunit,*,iostat=istat) iscstp
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0006 in read_md' )

      else if( word == '(charge es' ) then selectif

         read(iunit,*,iostat=istat) ichest
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0007 in read_md' )

      else if( word == '(wavefunct' ) then selectif

         read(iunit,*,iostat=istat) ihest
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0008 in read_md' )

      else if( word == '(planewave' ) then selectif

         read(iunit,*,iostat=istat) pwscale
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0009 in read_md' )

      else if( word == '(make tota' ) then selectif

         read(iunit,*,iostat=istat) lmomzero
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0010 in read_md' )

      else if( word == '(initial s' ) then selectif

         read(iunit,*,iostat=istat) nstep_ini
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0011 in read_md' )

      else if( word == '(optimizat' ) then selectif

         read(iunit,*,iostat=istat) ioptmze
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0012 in read_md' )

      else if( word == '(cell opti' ) then selectif

         read(iunit,*,iostat=istat) ioptmze_cell
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0013 in read_md' )

      else if( word == '(ASPC char' ) then selectif

         read(iunit,*,iostat=istat) laspc_chg
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0014 in read_md' )

      else if( word == '(ASPC wave' ) then selectif

         read(iunit,*,iostat=istat) laspc_wv
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0015 in read_md' )

      else if( word == '(ASPC acce' ) then selectif

         read(iunit,*,iostat=istat) laspc
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0016 in read_md' )

      else if( word == '(ASPC para' ) then selectif

         read(iunit,*,iostat=istat) aspc_corr
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0017 in read_md' )

      else if( word == '(ASPC mixi' ) then selectif

         read(iunit,*,iostat=istat) aslh_aspc, bslh_aspc
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0018 in read_md' )

      else if( word == '(ASPC iter' ) then selectif

         read(iunit,*,iostat=istat) iscfmx_aspc
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0019 in read_md' )

      else if( word == '(ASPC tole' ) then selectif

         read(iunit,*,iostat=istat) tolres_aspc
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0020 in read_md' )

      else if( word == '(ASPC corr' ) then selectif

         read(iunit,*,iostat=istat) nskp_aspc
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0021 in read_md' )

      else if( word == '(shock wav' ) then selectif

         read(iunit,*,iostat=istat) dummy
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0022 in read_md' )

         read(iunit,*,iostat=istat) nshockv(1:3)
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0023 in read_md' )

      else if( word == '(Extended ' ) then selectif

         read(iunit,*,iostat=istat) lxlbomd
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0024 in read_md' )

      else if( word == '(previous ' ) then selectif

         read(iunit,*,iostat=istat) kxlbomd
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0025 in read_md' )

!      else if( word == '(restricti' ) then selectif

!         read(iunit,*,iostat=istat) irstrct
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0027 in read_md' )

!      else if( word == '(sub-restr' ) then selectif

!         read(iunit,*,iostat=istat) irstrct_sub
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0028 in read_md' )

!      else if( word == '(MD cell e' ) then selectif

!         do i = 1, 3
!            read(iunit,*,iostat=istat) lcell_rstrct(i)
!            !----- error trap
!            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0029 in read_md' )
!         end do

      else if( word == '*end      ' ) then selectif

         exit

      endif selectif

   else readif

      !----- error trap
      call fstop( nfile, myid, nodes, 'error-0000 in read_md' )

   end if readif
end do readdo


!      !-----check structural optimization
!      if( ifmd == 1 .and. ioptmze == 1 ) then
!          ifmd = 2
!          treq = 0.d0
!      end if

!----- check cell optimization
if(    ifmd /= 1  &
& .or. ifmd == 1 .and. ioptmze == 1  &
& .or. ifmd == 1 .and. ioptmze == 10  &
& .or. ifmd == 1 .and. ioptmze == 11  &
& .or. ifmd == 1 .and. ioptmze == 20 ) ioptmze_cell = -1

!----- check MD cell edge restriction
if( irstrct == 0 ) then
    if( lcell_rstrct(1) .or. lcell_rstrct(2) .or. lcell_rstrct(3) ) irstrct = 7
end if


return
end subroutine




subroutine read_save( nfile, myid, nodes, iunit,  &
& lsave, lsreal8 )
!-----------------------------------------------------------------------
! read variables for save data
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
logical :: lsave, lsreal8

!-----declare local variables
character(10) :: word
integer :: istat


readdo: do
   read(iunit,'(a10)',iostat=istat) word

   readif: if( istat == 0 ) then

      selectif: if( word == '(on/off)  ' .or. word == '(how of it' ) then

         read(iunit,*,iostat=istat) lsave
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0001 in read_save' )

      else if( word == '(data type' ) then selectif

         read(iunit,*,iostat=istat) lsreal8
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0002 in read_save' )

      else if( word == '*end      ' ) then selectif

         exit

      endif selectif

   else readif

      !----- error trap
      call fstop( nfile, myid, nodes, 'error-0000 in read_save' )

   end if readif
end do readdo


return
end subroutine




subroutine read_str( nfile, myid, nodes, iunit,  &
& lstress, nskip_stress )
!-----------------------------------------------------------------------
! read variables for stress calculation
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
logical :: lstress
integer :: nskip_stress
character(10) :: word
integer :: istat


readdo: do
   read(iunit,'(a10)',iostat=istat) word

   readif: if( istat == 0 ) then

      selectif: if( word == '(on/off)  ' .or. word == '(how of it' ) then

         read(iunit,*,iostat=istat) lstress
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0001 in read_str' )

      else if( word == '(skip step' ) then selectif

         read(iunit,*,iostat=istat) nskip_stress
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0002 in read_str' )

      else if( word == '*end      ' ) then selectif

         exit

      endif selectif

   else readif

      !----- error trap
      call fstop( nfile, myid, nodes, 'error-0000 in read_str' )

   end if readif
end do readdo


return
end subroutine




subroutine read_dpwav( nfile, myid, nodes, iunit,  &
& ldpwav, nskip_dpwav, ibstt1, ibstt2, wv_xmin, wv_xmax, wv_ymin, wv_ymax, wv_zmin, wv_zmax,  &
& lcompress_dpwav, ndigit_dpwav )
!-----------------------------------------------------------------------
! read variables for dumping wavefunctions
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
logical :: ldpwav
integer :: nskip_dpwav, ibstt1, ibstt2
real*8  :: wv_xmin, wv_xmax, wv_ymin, wv_ymax, wv_zmin, wv_zmax
logical :: lcompress_dpwav
integer :: ndigit_dpwav

character(10) :: word
integer :: istat


readdo: do
   read(iunit,'(a10)',iostat=istat) word

   readif: if( istat == 0 ) then

      selectif: if( word == '(on/off)  ' .or. word == '(how of it' ) then

         read(iunit,*,iostat=istat) ldpwav
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0001 in read_dpwav' )

      else if( word == '(bands)   ' ) then selectif

         read(iunit,*,iostat=istat) ibstt1, ibstt2
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0002 in read_dpwav' )

      else if( word == '(skip step' ) then selectif

         read(iunit,*,iostat=istat) nskip_dpwav
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0003 in read_dpwav' )

      else if( word == '(output ar' ) then selectif

         read(iunit,*,iostat=istat) wv_xmin, wv_xmax
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0004 in read_dpwav' )

         read(iunit,*,iostat=istat) wv_ymin, wv_ymax
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0005 in read_dpwav' )

         read(iunit,*,iostat=istat) wv_zmin, wv_zmax
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0006 in read_dpwav' )

      else if( word == '(compresse' ) then selectif

         read(iunit,*,iostat=istat) lcompress_dpwav
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0007 in read_dpwav' )

      else if( word == '(significa' ) then selectif

         read(iunit,*,iostat=istat) ndigit_dpwav
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0008 in read_dpwav' )

      else if( word == '*end      ' ) then selectif

         exit

      endif selectif

   else readif

      !----- error trap
      call fstop( nfile, myid, nodes, 'error-0000 in read_dpwav' )

   end if readif
end do readdo


return
end subroutine




subroutine read_cell( nfile, myid, nodes, iunit, hcell )
!-----------------------------------------------------------------------
! read variables for supercell
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
real*8, dimension(3,3) :: hcell

!-----declare local variables
character(10) :: word
character(10) :: unit = '(bohr)    '
integer :: istat, i, ix
real*8  :: dltca, dltcb, dltcc, angalf, angbet, anggam
real*8  :: angcon, E1X, E1Y, E1Z, E2X, E2Y, E2Z, E3X, E3Y, E3Z
real*8,  parameter :: audang = 0.529177249d0


readdo: do
   read(iunit,'(a10)',iostat=istat) word

   readif: if( istat == 0 ) then

      selectif: if( word == '(cell vect' ) then

         do i = 1, 3
            read(iunit,*,iostat=istat) ( hcell(ix,i), ix = 1, 3 )
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0001 in read_cell' )
         end do

      else if( word == '(lengths &' ) then selectif

         read(iunit,*,iostat=istat) dltca, dltcb, dltcc
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0002 in read_cell' )

         read(iunit,*,iostat=istat) angalf, angbet, anggam
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0003 in read_cell' )

         angcon = acos(-1.d0)/180.d0
         angalf = angalf * angcon
         angbet = angbet * angcon
         anggam = anggam * angcon
      !--- unit vectors parallel to cell vectors ---
         E1X = 1.0
         E1Y = 0.0
         E1Z = 0.0
         E2X = DCOS(anggam)
         E2Y = DSIN(anggam)
         E2Z = 0.0
         E3X = DCOS(angbet)
         E3Y = DCOS(angalf) - E3X*E2X
         E3Y = E3Y/E2Y
         E3Z = 1.0 - E3X*E3X - E3Y*E3Y
         E3Z = DSQRT(E3Z)

         hcell(1,1) = dltca*E1X
         hcell(2,1) = dltca*E1Y
         hcell(3,1) = dltca*E1Z
         hcell(1,2) = dltcb*E2X
         hcell(2,2) = dltcb*E2Y
         hcell(3,2) = dltcb*E2Z
         hcell(1,3) = dltcc*E3X
         hcell(2,3) = dltcc*E3Y
         hcell(3,3) = dltcc*E3Z

      else if( word == '(unit of l' ) then selectif

         read(iunit,'(a10)',iostat=istat) unit
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0004 in read_cell' )

      else if( word == '*end      ' ) then selectif

         !----- length unit in [bohr]
         if( unit == '(ang)     ' ) then
             do i = 1, 3
             do ix = 1, 3
                hcell(ix,i) = hcell(ix,i) / audang
             end do
             end do
         end if

         exit

      end if selectif

   else readif

      !----- error trap
      call fstop( nfile, myid, nodes, 'error-0000 in read_cell' )

   end if readif
end do readdo


return
end subroutine




subroutine read_pw( nfile, myid, nodes, iunit,  &
& ecut, ecutdens, ecutsoft, ecutorth, ecutc )
!-----------------------------------------------------------------------
! read variables for plane-wave cutoff energy
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
real*8  :: ecut, ecutdens, ecutsoft, ecutorth, ecutc
character(10) :: word
character(10) :: unit = '(ry)      '
integer :: istat
real*8  :: rydev = 27.2116d0/2.d0


readdo: do
   read(iunit,'(a10)',iostat=istat) word

   readif: if( istat == 0 ) then

      selectif: if( word == '(for wavef' ) then

         read(iunit,*,iostat=istat) ecut
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0001 in read_pw' )

      else if( word == '(for elect' ) then selectif

         read(iunit,*,iostat=istat) ecutdens
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0002 in read_pw' )

      else if( word == '(for soft ' ) then selectif

         read(iunit,*,iostat=istat) ecutsoft
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0003 in read_pw' )

      else if( word == '(for ortho' ) then selectif

         read(iunit,*,iostat=istat) ecutorth
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0004 in read_pw' )

      else if( word == '(for charg' ) then selectif

         read(iunit,*,iostat=istat) ecutc
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0005 in read_pw' )

      else if( word == '(unit of c' ) then selectif

         read(iunit,'(a10)',iostat=istat) unit
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0006 in read_pw' )

      else if( word == '*end      ' ) then selectif

         !----- energy unit in [Ryd.]
         unitif: if( unit == '(hr)      ' ) then
            ecut     = ecut     * 2.d0
            ecutdens = ecutdens * 2.d0
            ecutsoft = ecutsoft * 2.d0
            ecutorth = ecutorth * 2.d0
            ecutc    = ecutc    * 2.d0
         else if( unit == '(ev)      ' ) then unitif
            ecut     = ecut     / rydev
            ecutdens = ecutdens / rydev
            ecutsoft = ecutsoft / rydev
            ecutorth = ecutorth / rydev
            ecutc    = ecutc    / rydev
         end if unitif
         exit

      endif selectif

   else readif

      !----- error trap
      call fstop( nfile, myid, nodes, 'error-0000 in read_pw' )

   end if readif
end do readdo


return
end subroutine




subroutine read_band( nfile, myid, nodes, iunit,  &
& noband, nband, lfermi, tfermi, ncion )
!-----------------------------------------------------------------------
! read variables for electronic bands
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
integer :: noband, nband
integer :: neband = 0             ! the number of empty bands
integer :: lfermi
real*8  :: tfermi
integer :: ncion
character(10) :: word
integer :: istat


readdo: do
   read(iunit,'(a10)',iostat=istat) word

   readif: if( istat == 0 ) then

      selectif: if( word == '(occupied ' ) then

         read(iunit,*,iostat=istat) noband
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0001 in read_band' )

      else if( word == '(empty ban' ) then selectif

         read(iunit,*,iostat=istat) neband
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0002 in read_band' )

      else if( word == '(broadenin' ) then selectif

         read(iunit,*,iostat=istat) lfermi, tfermi
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0003 in read_band' )

      else if( word == '(charge nu' ) then selectif

         read(iunit,*,iostat=istat) ncion
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0004 in read_band' )

      else if( word == '*end      ' ) then selectif

         !----- the number of total bands
         nband = noband + neband
         exit

      endif selectif

   else readif

      !----- error trap
      call fstop( nfile, myid, nodes, 'error-0000 in read_band' )

   end if readif
end do readdo


return
end subroutine




subroutine read_lsda( nfile, myid, nodes,  &
& iunit, lspin, lnoncollinear, lfixud, diffud, lwfrand, lduplicate,  &
& latmref, lclearamag, refmx, refmy, refmz,  &
& lclearmagne_in_x, lclearmagne_in_y, lclearmagne_in_z,  &
& lfixmagne_in_x, lfixmagne_in_y, lfixmagne_in_z,  &
& fixmagne_in_x, fixmagne_in_y, fixmagne_in_z )
!-----------------------------------------------------------------------
! read variables for spin-polarized calculations
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
logical :: lspin, lnoncollinear
logical :: lfixud
real*8  :: diffud
logical :: lwfrand, lduplicate
logical :: latmref, lclearamag
real*8  :: refmx, refmy, refmz
logical :: lclearmagne_in_x, lclearmagne_in_y, lclearmagne_in_z
logical :: lfixmagne_in_x, lfixmagne_in_y, lfixmagne_in_z
real(8) :: fixmagne_in_x, fixmagne_in_y, fixmagne_in_z

!---declare local variables
character(10) :: word
integer :: inispin =  1
integer :: nspinat = -1
real(8) :: refabs
integer :: istat


readdo: do
   read(iunit,'(a10)',iostat=istat) word

   readif: if( istat == 0 ) then

      selectif: if( word == '(on/off)  ' .or. word == '(how of it' ) then

         read(iunit,*,iostat=istat) lspin
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0001 in read_lsda' )

      else if( word == '(fix spin ' ) then selectif

         read(iunit,*,iostat=istat) lfixud
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0002 in read_lsda' )

      else if( word == '(spin pola' ) then selectif

         read(iunit,*,iostat=istat) diffud
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0003 in read_lsda' )

      else if( word == '(initial w' ) then selectif

         read(iunit,*,iostat=istat) lwfrand
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0004 in read_lsda' )

      else if( word == '(duplicate' ) then selectif

         read(iunit,*,iostat=istat) lduplicate
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0005 in read_lsda' )

      else if( word == '(initial s' ) then selectif

         read(iunit,*,iostat=istat) inispin
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0006 in read_lsda' )

      else if( word == '(ferromagn' ) then selectif

         if( inispin == 2 ) then

             read(iunit,*,iostat=istat) nspinat
             !----- error trap
             if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0007 in read_lsda' )

             call read_lsda2( nfile, myid, nodes, iunit, inispin, nspinat, lnoncollinear )

         end if

      else if( word == '(antiferro' ) then selectif

         if( inispin == 3 ) then

             read(iunit,*,iostat=istat) nspinat
             !----- error trap
             if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0008 in read_lsda' )

             call read_lsda2( nfile, myid, nodes, iunit, inispin, nspinat, lnoncollinear )

         end if

      else if( word == '(atomic sp' ) then selectif

         if( inispin == 4 ) then

             read(iunit,*,iostat=istat) nspinat
             !----- error trap
             if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0009 in read_lsda' )

             call read_lsda2( nfile, myid, nodes, iunit, inispin, nspinat, lnoncollinear )

         end if

!      else if( word == '(uniform s' ) then selectif
!
!         if( lnoncollinear .and. inispin == 1 ) then
!
!             nspinat = 1
!             call read_lsda3( nfile, myid, nodes, iunit )
!
!         end if
!
!      else if( word == '(noncollin' ) then selectif
!
!         read(iunit,*,iostat=istat) lnoncollinear
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0010 in read_lsda' )

!      else if( word == '(reference' ) then selectif
!
!         read(iunit,*,iostat=istat) refmx, refmy, refmz
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0011 in read_lsda' )

!      else if( word == '(clear mag' ) then selectif
!
!         read(iunit,*,iostat=istat) lclearmagne_in_x
!         read(iunit,*,iostat=istat) lclearmagne_in_y
!         read(iunit,*,iostat=istat) lclearmagne_in_z
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0012 in read_lsda' )

!      else if( word == '(atomic re' ) then selectif
!
!         read(iunit,*,iostat=istat) latmref
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0013 in read_lsda' )

!      else if( word == '(clear ato' ) then selectif
!
!         read(iunit,*,iostat=istat) lclearamag
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0014 in read_lsda' )

!      else if( word == '(fix magne' ) then selectif
!
!         read(iunit,*,iostat=istat) lfixmagne_in_x, fixmagne_in_x
!         read(iunit,*,iostat=istat) lfixmagne_in_y, fixmagne_in_y
!         read(iunit,*,iostat=istat) lfixmagne_in_z, fixmagne_in_z
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0015 in read_lsda' )

      else if( word == '*end      ' ) then selectif

         if( nspinat < 0 )  &
&            call read_lsda2( nfile, myid, nodes, iunit, inispin, 0, lnoncollinear )

         exit

      endif selectif

   else readif

      !----- error trap
      call fstop( nfile, myid, nodes, 'error-0000 in read_lsda' )

   end if readif
end do readdo


!---make sure that lnoncollinear = .false. when lspin = .false.
if( .not.lspin ) lnoncollinear = .false.

!if( lnoncollinear ) then
!    !---set lspin = .false. when lnoncollinear = .false.
!    lspin = .false.
!    !---normalize the reference vector if non-zero vector
!    refabs = max(sqrt(refmx*refmx + refmy*refmy + refmz*refmz), 1d-05)
!    refmx = refmx/refabs
!    refmy = refmy/refabs
!    refmz = refmz/refabs
!end if


return
end subroutine




subroutine read_lsda2( nfile, myid, nodes, iunit, inispin_, nspinat_, lnoncollinear )
!-----------------------------------------------------------------------
! read variables for spin-polarized calculations
!-----------------------------------------------------------------------
!use param
use param_inispin
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
integer :: inispin_, nspinat_
logical :: lnoncollinear

!-----declare local variables
integer :: status, istat
integer :: i
real*8  :: um


inispin = inispin_
nspinat = nspinat_

!------allocate memory
allocate( iatspin(nspinat+1), atmmagne(nspinat+1), stat=status )
if( lnoncollinear .and. status == 0 ) then
    allocate( umx(nspinat+1), umy(nspinat+1), umz(nspinat+1), stat=status )
end if

!----- error trap
if( status /= 0 ) call fstop( nfile, myid, nodes,  &
&                 'allocate memory error in read_lsda2' )

do i = 1, nspinat

   if( .not.lnoncollinear ) then
       read(iunit,*,iostat=istat) iatspin(i), atmmagne(i)
   else
       read(iunit,*,iostat=istat) iatspin(i), atmmagne(i), umx(i), umy(i), umz(i)
       !---normalization
       um = sqrt(umx(i)*umx(i) + umy(i)*umy(i) + umz(i)*umz(i))
       umx(i) = umx(i)/um
       umy(i) = umy(i)/um
       umz(i) = umz(i)/um
   end if
   !----- error trap
   if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                               'error-0001 in read_lsda2' )

end do


return
end subroutine




subroutine read_atom_number( nfile, myid, nodes, iunit,  &
& ntype, natom )
!-----------------------------------------------------------------------
!     get the numbers of species & atoms
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
integer :: ntype
integer :: natom

!-----declare local variables
character(10) :: word
integer :: istat
integer :: it


readdo: do
   read(iunit,'(a10)',iostat=istat) word

   readif: if( istat == 0 ) then

      selectif: if( word == '(species) ' ) then

         read(iunit,*,iostat=istat) ntype

         atomif: if( istat == 0 ) then

            natom = 0
            do it = 1, ntype
               call read_natom( nfile, myid, nodes, iunit, natom )
            end do

         else atomif
            !----- error trap
            call fstop( nfile, myid, nodes,  &
&                         'error-0001 in read_atom_number' )
         end if atomif

      else if( word == '*end      ' ) then selectif

         exit

      endif selectif

   else readif

      !----- error trap
      call fstop( nfile, myid, nodes,  &
&                 'error-0000 in read_atom_number' )

   end if readif
end do readdo


!      !-----rewind file
!      rewind(iunit)
!      
!      !-----return correct position
!      rereaddo: do
!         read(iunit,'(a10)',iostat=istat) word
!
!         rereadif: if( istat == 0 ) then
!
!            reselectif: if( word == '*atoms    ' ) then
!
!               exit
!
!            endif reselectif
!
!         else rereadif
!
!            !----- error trap
!            call fstop( nfile, myid, nodes,
!     &                  'error-0010 in read_atom_number' )
!
!         end if rereadif
!      end do rereaddo



return
end subroutine




subroutine read_natom( nfile, myid, nodes, iunit, natom )
!-----------------------------------------------------------------------
! read the number of atoms
!-----------------------------------------------------------------------
use outfile
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
integer :: natom

!-----declare local variables
character(10) :: word
character(50) :: fname
integer :: istat
integer :: nhk
integer :: icscale
integer :: keyword
integer :: ierror, ntotal, i, n_i
real*8  :: x_i, y_i, z_i
integer :: junit


call allocate_unit_number( junit )

readdo: do
   read(iunit,'(a10)',iostat=istat) word

   readif: if( istat == 0 ) then

      selectif: if( word == '(the numbe' ) then

         read(iunit,*,iostat=istat) nhk
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0001 in read_natom' )

         config_1: if( nhk > 0 ) then

            natom = natom + nhk

         else config_1

            rereaddo: do
               read(iunit,'(a10)',iostat=istat) word

               rereadif: if( istat == 0 ) then

                  reselectif: if( word == '(position ' ) then

         !----- read configuration from a file
            read(iunit,*,iostat=istat) fname
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0010 in read_natom' )

            read(iunit,*,iostat=istat) icscale
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0011 in read_natom' )

            !----- open configuration file
            if(loutfile(1)) write(nfile(1),*) 'open file: ',  &
&                                 fname(1:len_trim(fname))
            if(loutfile(2)) write(nfile(2),*) 'open file: ',  &
&                                 fname(1:len_trim(fname))
            open( junit, file=fname, status='old', action='read',  &
&                 iostat=ierror )
            blockif_1: if( ierror == 0 ) then

         readdo_config: do

               read(iunit,*,iostat=istat) keyword
               !----- error trap
               if( istat < 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0012 in read_natom' )

               if( istat > 0 ) exit

               read(junit,*,iostat=istat) ntotal
               !----- error trap
               if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                               'error-0013 in read_natom' )

               !----- read configuration
               do i = 1, ntotal
                  read(junit,*,iostat=istat) n_i, x_i, y_i, z_i
                  !----- error trap
                  if( istat /= 0 ) call fstop( nfile, myid, nodes, &
&                                  'error-0014 in read_natom' )
                  if( n_i == keyword ) then
                      natom = natom + 1
                  end if
               end do

               rewind(junit)

         end do readdo_config

               close(junit)

            else blockif_1

               !----- error trap
               call fstop( nfile, myid, nodes,  &
&                          'cannot open file :'//fname )

            endif blockif_1

                     exit

                  endif reselectif

               else rereadif

                  !----- error trap
                  call fstop( nfile, myid, nodes,  &
&                             'error-1000 in read_natom' )

               end if rereadif
            end do rereaddo

         end if config_1

         exit

      endif selectif

   else readif

      !----- error trap
      call fstop( nfile, myid, nodes, 'error-0000 in read_natom' )

   end if readif
end do readdo

call deallocate_unit_number( junit )


return
end subroutine




subroutine read_atom( nfile, myid, nodes, iunit,  &
& ntype, natom, lplcom, zatom, lkbppi, lvandi, lvflag, llocli,  &
& lkbpp, lvand, lvfull, llocl, lking, rking, gkgmax, gkgexct,  &
& llking, rlking, glkgmax, glkgexct,  &
& lpcci, rpcc, lpking, rpking, gpkgmax, gpkgexct, nhk,  &
& icscale, vrandom, ratm, vatm,  &
& rdecmp, rxmulk, lmulao, lmulbsis, nmulao, mulaox, rsphexp, radeda,  &
& rintchg, frac_rcomp_it, lplusU_at, plusU_l, plusU_U, plusU_J,  &
& lplusC_at, plusC_r, plusC_e, mxlx )
!-----------------------------------------------------------------------
! read variables for atoms
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
integer :: ntype
integer :: natom
logical :: lplcom
real*8,  dimension(ntype) :: zatom
logical ::                   lkbpp,  lvand,  lvfull, llocl
logical, dimension(ntype) :: lkbppi, lvandi, lvflag, llocli
logical, dimension(ntype) :: lking
real*8,  dimension(ntype) :: rking, gkgmax, gkgexct
logical, dimension(ntype) :: llking
real*8,  dimension(ntype) :: rlking, glkgmax, glkgexct
logical, dimension(ntype) :: lpcci
real*8,  dimension(ntype) :: rpcc
logical, dimension(ntype) :: lpking
real*8,  dimension(ntype) :: rpking, gpkgmax, gpkgexct
integer, dimension(ntype) :: nhk
integer, dimension(ntype) :: icscale
logical, dimension(ntype) :: vrandom
real*8,  dimension(3,natom) :: ratm
real*8,  dimension(3,natom) :: vatm
integer :: mulaox
real*8,  dimension(2,mulaox,ntype) :: rdecmp
real*8,  dimension(mulaox,ntype)   :: rxmulk
integer, dimension(mulaox,ntype)   :: lmulao
logical, dimension(mulaox,ntype)   :: lmulbsis
integer, dimension(ntype)          :: nmulao
real*8,  dimension(ntype) :: rsphexp
real*8,  dimension(ntype) :: radeda
real*8,  dimension(ntype) :: rintchg
real*8,  dimension(ntype) :: frac_rcomp_it
logical, dimension(ntype) :: lplusU_at
integer, dimension(ntype) :: plusU_l
real*8,  dimension(ntype) :: plusU_U, plusU_J
logical, dimension(ntype) :: lplusC_at
integer :: mxlx
real*8,  dimension(0:mxlx,ntype) :: plusC_r, plusC_e

!-----declare local variables
character(10) :: word
integer :: istat, it, nhk1


readdo: do
   read(iunit,'(a10)',iostat=istat) word

   readif: if( istat == 0 ) then

      selectif: if( word == '(species) ' ) then

         read(iunit,*,iostat=istat) ntype

         atomif: if( istat == 0 ) then

            nhk1 = 1
            lkbpp = .false.
            lvand = .false.
            lvfull= .false.
            llocl = .false.
            do it = 1, ntype

               call read_type( nfile, myid, nodes, iunit,  &
& zatom(it), lkbppi(it), lvandi(it), lvflag(it), llocli(it),  &
& lking(it), rking(it), gkgmax(it), gkgexct(it),  &
& llking(it), rlking(it), glkgmax(it), glkgexct(it),  &
& lpcci(it), rpcc(it), lpking(it), rpking(it), gpkgmax(it),  &
& gpkgexct(it), nhk(it), icscale(it), vrandom(it),  &
& ratm(1,nhk1), vatm(1,nhk1),  &
& rdecmp(1,1,it), rxmulk(1,it), lmulao(1,it), lmulbsis(1,it),  &
& nmulao(it), mulaox,  &
& rsphexp(it), radeda(it), rintchg(it), frac_rcomp_it(it),  &
& lplusU_at(it), plusU_l(it), plusU_U(it), plusU_J(it),  &
& lplusC_at(it), plusC_r(0,it), plusC_e(0,it), mxlx, natom-nhk1+1 )

               nhk1 = nhk1 + nhk(it)

               lkbpp = lkbpp .or. lkbppi(it)
               lvand = lvand .or. lvandi(it)
               lvfull= lvfull.or. lvflag(it)
               llocl = llocl .or. llocli(it)
            end do

         else atomif
            !----- error trap
            call fstop( nfile, myid, nodes,  &
&                         'error-0001 in read_atom' )
         end if atomif

      else if( word == '(place c.m' ) then selectif

         read(iunit,*,iostat=istat) lplcom
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0002 in read_atom' )

      else if( word == '*end      ' ) then selectif

         exit

      endif selectif

   else readif

      !----- error trap
      call fstop( nfile, myid, nodes, 'error-0000 in read_atom' )

   end if readif
end do readdo


return
end subroutine




subroutine read_type( nfile, myid, nodes, iunit,  &
& zatom, lkbppi, lvandi, lvflag, llocli,  &
& lking, rking, gkgmax, gkgexct, llking, rlking, glkgmax, glkgexct,  &
& lpcci, rpcc, lpking, rpking, gpkgmax, gpkgexct, nhk,  &
& icscale, vrandom, ratm, vatm,  &
& rdecmp, rxmulk, lmulao, lmulbsis, nmulao, mulaox, rsphexp, radeda,  &
& rintchg, frac_rcomp_it, lplusU_at, plusU_l, plusU_U, plusU_J,  &
& lplusC_at, plusC_r, plusC_e, mxlx, mxn )
!-----------------------------------------------------------------------
! read variables for atomic data
!-----------------------------------------------------------------------
use outfile
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
real*8  :: zatom
logical :: lkbppi, lvandi, lvflag, llocli
logical :: lking
real*8  :: rking, gkgmax, gkgexct
logical :: llking
real*8  :: rlking, glkgmax, glkgexct
logical :: lpcci
real*8  :: rpcc
logical :: lpking
real*8  :: rpking, gpkgmax, gpkgexct
integer :: nhk
integer :: icscale
logical :: vrandom
real*8,  dimension(3,*) :: ratm
real*8,  dimension(3,*) :: vatm
integer :: mulaox
real*8,  dimension(2,mulaox) :: rdecmp
real*8,  dimension(mulaox)   :: rxmulk
integer, dimension(mulaox)   :: lmulao
logical, dimension(mulaox)   :: lmulbsis
integer                      :: nmulao
real*8  :: rsphexp
real*8  :: radeda
real*8  :: rintchg
real*8  :: frac_rcomp_it
logical :: lplusU_at
integer :: plusU_l
real*8  :: plusU_U, plusU_J
logical :: lplusC_at
integer :: mxlx
real*8  :: plusC_r(0:mxlx), plusC_e(0:mxlx)
integer :: mxn

!-----declare local variables
character(10) :: word, pptype
character(50) :: fname
character(10) :: unit
real*8,  parameter :: audang = 0.529177249d0
integer :: istat, keyword, ierror, ntotal, i, l, l2
integer :: n_i, nion, nionv, ncl1, ncl2, ncl3, ix, iy, iz
real*8  :: x_i, y_i, z_i
real*8  :: vmax, rrxmulk
integer :: junit
logical :: lreadmul, lreadpp, lreadnlc, lreadlc, lreadpcc


call allocate_unit_number( junit )

nion  = 0
nionv = 0
lreadmul = .true.
lreadnlc = .false.
lreadlc  = .false.
lreadpcc = .false.
unit = '(bohr)    '

lreadpp= .false.
lkbppi = .false.
lvandi = .false.
lvflag = .false.
llocli = .false.

readdo: do
   read(iunit,'(a10)',iostat=istat) word

   readif: if( istat == 0 ) then

      selectif: if( word == '(atomic nu' ) then

         read(iunit,*,iostat=istat) zatom
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0001 in read_type' )

      else if( word == '(pseudopot' ) then selectif

         read(iunit,'(a10)',iostat=istat) pptype
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0002 in read_type' )

         !----- set pseudopotential type
         lkbppi = pptype == 'kbpp      '
         lvandi = pptype == 'uspp      '.or.pptype == 'vand      '
         lvflag = pptype == 'uspp      '
         llocli = pptype == 'local     '

         lreadpp  = .true.

      else if( word == '(nonlocal ' ) then selectif

         read(iunit,*,iostat=istat) lking, rking, gkgmax, gkgexct
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0003 in read_type' )

         lreadnlc = .true.

      else if( word == '(local pot' ) then selectif

         read(iunit,*,iostat=istat) llking,rlking,glkgmax,glkgexct
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0004 in read_type' )

         lreadlc  = .true.

      else if( word == '(partial c' ) then selectif

         read(iunit,*,iostat=istat) lpcci, rpcc
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0005 in read_type' )

         read(iunit,*,iostat=istat) lpking,rpking,gpkgmax,gpkgexct
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0006 in read_type' )

         lreadpcc = .true.

      else if( word == '(the numbe' ) then selectif

         read(iunit,*,iostat=istat) nhk
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0007 in read_type' )

         !----- error trap
         if( nhk > mxn ) call fstop( nfile, myid, nodes,  &
&                                    'nhk > mxn' )

      else if( word == '(compensat' ) then selectif

         read(iunit,*,iostat=istat) frac_rcomp_it
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0008 in read_type' )

      else if( word == '(position ' ) then selectif

         config_1: if( nhk == 0 ) then

         !----- read configuration from a file
            read(iunit,*,iostat=istat) fname
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0010 in read_type' )

            read(iunit,*,iostat=istat) icscale
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0011 in read_type' )

            !----- open configuration file
            if(loutfile(1)) write(nfile(1),*) 'open file: ',  &
&                                 fname(1:len_trim(fname))
            if(loutfile(2)) write(nfile(2),*) 'open file: ',  &
&                                 fname(1:len_trim(fname))
            open( junit, file=fname, status='old', action='read',  &
&                 iostat=ierror )
            blockif_1: if( ierror == 0 ) then

            nion = 0
            readdo_config: do

               read(iunit,*,iostat=istat) keyword
               !----- error trap
               if( istat < 0 ) call fstop( nfile, myid, nodes,  &
&                               'error-0012 in read_type' )

               if( istat > 0 ) exit

               read(junit,*,iostat=istat) ntotal
               !----- error trap
               if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                               'error-0013 in read_type' )

               !----- read configuration
               do i = 1, ntotal
                  read(junit,*,iostat=istat) n_i, x_i, y_i, z_i
                  !----- error trap
                  if( istat /= 0 ) call fstop( nfile, myid, nodes, &
&                                  'error-0014 in read_type' )
                  if( n_i == keyword ) then
                      nion = nion + 1
                      !----- error trap
                      if( nion > mxn ) call fstop( nfile, myid,  &
&                                      nodes, 'nion > mxn' )
                      ratm(1,nion) = x_i
                      ratm(2,nion) = y_i
                      ratm(3,nion) = z_i
                  end if
               end do

               rewind(junit)

            end do readdo_config

               close(junit)

            else blockif_1

               !----- error trap
               call fstop( nfile, myid, nodes,  &
&                          'cannot open file :'//fname )

            endif blockif_1

         end if config_1

      else if( word == '(positions' ) then selectif

         config_2: if( nhk > 0 ) then

            read(iunit,*,iostat=istat) ncl1, ncl2, ncl3
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0020 in read_type' )

            !----- error trap
            if( mod(nhk,(ncl1*ncl2*ncl3)) /= 0 )  &
&               call fstop( nfile, myid, nodes,  &
&                           'wrong # of ions in read_type' )

            read(iunit,*,iostat=istat) icscale
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0021 in read_type' )

            !----- error trap
            if( icscale /= 1 .and. ncl1*ncl2*ncl3 > 1 )  &
&               call fstop( nfile, myid, nodes,  &
&      'not supported for icscale /= 1 .and. ncl1*ncl2*ncl3 > 1' )

            !----- read configuration
            nion = 0
            do i = 1, nhk/(ncl1*ncl2*ncl3)
               read(iunit,*,iostat=istat) x_i, y_i, z_i
               !----- error trap
               if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                               'error-0022 in read_type' )

               scaleif: if( icscale == 1 ) then
                  !----- set scaled coordinates
                  do iz = 1, ncl3
                  do iy = 1, ncl2
                  do ix = 1, ncl1
                     nion = nion + 1
                     ratm(1,nion) = (x_i + dble(ix-1))/dble(ncl1)
                     ratm(2,nion) = (y_i + dble(iy-1))/dble(ncl2)
                     ratm(3,nion) = (z_i + dble(iz-1))/dble(ncl3)
                  end do
                  end do
                  end do
               else scaleif
                  !----- set real coordinates
                  nion = nion + 1
                  ratm(1,nion) = x_i
                  ratm(2,nion) = y_i
                  ratm(3,nion) = z_i
               end if scaleif

            end do

         end if config_2

      else if( word == '(velocity ' ) then selectif

         velocity_1: if( nhk == 0 ) then

            vrandom = .false.

         !----- read velocities from a file
            read(iunit,*,iostat=istat) fname
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0030 in read_type' )

            read(iunit,*,iostat=istat) keyword
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0031 in read_type' )

            !----- open velocity file
            if(loutfile(1)) write(nfile(1),*) 'open file: ',  &
&                                 fname(1:len_trim(fname))
            if(loutfile(2)) write(nfile(2),*) 'open file: ',  &
&                                 fname(1:len_trim(fname))
            open( junit, file=fname, status='old', action='read',  &
&                 iostat=ierror )
            blockif_2: if( ierror == 0 ) then

               read(junit,*,iostat=istat) ntotal
               !----- error trap
               if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                               'error-0032 in read_type' )

               read(junit,*,iostat=istat) vmax
               !----- error trap
               if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                               'error-0033 in read_type' )

               !----- read velocities
               nionv = 0
               do i = 1, ntotal
                  read(junit,*,iostat=istat) n_i, x_i, y_i, z_i
                  !----- error trap
                  if( istat /= 0 ) call fstop( nfile, myid, nodes, &
&                                  'error-0034 in read_type' )
                  if( n_i == keyword ) then
                      nionv = nionv + 1
                      !----- error trap
                      if( nionv > mxn ) call fstop( nfile, myid,  &
&                                       nodes, 'nionv > mxn' )
                      vatm(1,nionv) = x_i * vmax
                      vatm(2,nionv) = y_i * vmax
                      vatm(3,nionv) = z_i * vmax
                  end if
               end do
               close(junit)

            else blockif_2

               vrandom = .true.
               !----- error trap
!                     call fstop( nfile, myid, nodes,  &
!     &                           'cannot open file :'//fname)

            endif blockif_2

         end if velocity_1

      else if( word == '(velocitie' ) then selectif

         velocity_2: if( nhk > 0 ) then

            vrandom = .false.

            read(iunit,*,iostat=istat) ncl1, ncl2, ncl3
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0040 in read_type' )

            !----- error trap
            if( mod(nhk,(ncl1*ncl2*ncl3)) /= 0 )  &
&               call fstop( nfile, myid, nodes,  &
&                           'wrong # of ions (2) in read_type' )

            read(iunit,*,iostat=istat) vmax
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0041 in read_type' )

            !----- read velocities
            nionv = 0
            do i = 1, nhk/(ncl1*ncl2*ncl3)
               read(iunit,*,iostat=istat) x_i, y_i, z_i
               !----- error trap
               if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                               'error-0042 in read_type' )

               do iz = 1, ncl3
               do iy = 1, ncl2
               do ix = 1, ncl1
                  nionv = nionv + 1
                  vatm(1,nionv) = x_i * vmax
                  vatm(2,nionv) = y_i * vmax
                  vatm(3,nionv) = z_i * vmax
               end do
               end do
               end do

            end do

         end if velocity_2

      !----- Note: classical style ----------------------
      else if( word == '(Mulliken ' .and. lreadmul ) then selectif

         read(iunit,*,iostat=istat) rrxmulk
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0050 in read_type' )

         nmulao = 4
         do i = 1, nmulao
            lmulao(i) = (i-1)*10
            lmulbsis(i) = .true.
            rxmulk(i) = rrxmulk
            read(iunit,*,iostat=istat) rdecmp(1,i), rdecmp(2,i)
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0051 in read_type' )
         end do

      !----- Note: new style ----------------------------
      else if( word == '(Populatio' .and. lreadmul ) then selectif

         read(iunit,*,iostat=istat) nmulao
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0052 in read_type' )

         do i = 1, nmulao
            read(iunit,*,iostat=istat) lmulao(i), rxmulk(i),  &
&                                      rdecmp(1,i), rdecmp(2,i)
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0053 in read_type' )
         end do

         do i = 1, nmulao
            lmulbsis(i) = .true.
         end do

!         lreadmul = .false.

      !----- Note: new style ----------------------------
      else if( word == '(Rev. Popu' ) then selectif

         read(iunit,*,iostat=istat) nmulao
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0052 in read_type' )

         do i = 1, nmulao
            read(iunit,*,iostat=istat) lmulao(i), lmulbsis(i), rxmulk(i),  &
&                                      rdecmp(1,i), rdecmp(2,i)
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0053 in read_type' )
         end do

         lreadmul = .false.

      else if( word == '(spherical' ) then selectif

         read(iunit,*,iostat=istat) rsphexp
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0054 in read_type' )

      else if( word == '(EDA)     ' ) then selectif

         read(iunit,*,iostat=istat) radeda
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0055 in read_type' )

      else if( word == '(atomic ch' ) then selectif

         read(iunit,*,iostat=istat) rintchg
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0056 in read_type' )

      else if( word == '(unit of l' ) then selectif

         read(iunit,'(a10)',iostat=istat) unit
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0057 in read_type' )

!      else if( word == '(DFT+U)   ' ) then selectif
!
!         read(iunit,*,iostat=istat) lplusU_at
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0058 in read_type' )
!
!         read(iunit,*,iostat=istat) plusU_l
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0059 in read_type' )
!
!         read(iunit,*,iostat=istat) plusU_U
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0060 in read_type' )
!
!         read(iunit,*,iostat=istat) plusU_J
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0061 in read_type' )
!
!      else if( word == '(DFT+C)   ' ) then selectif
!
!         read(iunit,*,iostat=istat) lplusC_at
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0065 in read_type' )
!
!         do l = 0, mxlx
!            read(iunit,*,iostat=istat) plusC_r(l), plusC_e(l)
!            !----- error trap
!            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                            'error-0066 in read_type' )
!         end do

      else if( word == '(end)     ' ) then selectif

         exit

      endif selectif

   else readif

      !----- error trap
      call fstop( nfile, myid, nodes, 'error-0000 in read_type' )

   end if readif
end do readdo


!----- length unit in [bohr]
if( icscale /= 1 ) then
    if( unit == '(ang)     ' ) then
        ratm(1:3,1:nion) = ratm(1:3,1:nion) / audang
    end if
end if

!----- error trap
if( nion == 0 ) call fstop( nfile, myid, nodes, 'nion == 0' )
if( nionv /= 0 .and. nionv /= nion ) then
    if(loutfile(1)) write(nfile(1),*) 'nion, nionv =', nion, nionv
    if(loutfile(2)) write(nfile(2),*) 'nion, nionv =', nion, nionv
    call fstop( nfile, myid, nodes, 'nionv /= nion' )
endif

!----- set the number of ions in the supercell
nhk = nion

call deallocate_unit_number( junit )


if( lreadpp .and. lreadnlc .and. lreadlc .and. lreadpcc ) return

!---default settings
select case(nint(zatom))
case(1)
    !---H
    if( .not.lreadpp  ) lkbppi = .true.
    if( .not.lreadnlc )  lking = .true.;  rking = 1.5d0;  gkgmax = 1.25d0;  gkgexct = 0.8d0
    if( .not.lreadlc  ) llking =.false.; rlking = 1.5d0; glkgmax = 1.15d0; glkgexct = 0.8d0
    if( .not.lreadpcc ) then
        lpcci = .false.; rpcc = 1.4d0
        lpking=.true.; rpking = 1.1d0; gpkgmax = 1.15d0; gpkgexct = 0.8d0
    end if
case(8)
    !---O
    if( .not.lreadpp  ) lkbppi = .true.
    if( .not.lreadnlc )  lking = .true.;  rking = 1.5d0;  gkgmax = 1.25d0;  gkgexct = 0.5d0
    if( .not.lreadlc  ) llking =.false.; rlking = 1.5d0; glkgmax = 1.15d0; glkgexct = 0.8d0
    if( .not.lreadpcc ) then
        lpcci = .true.; rpcc = 1.4d0
        lpking= .true.; rpking = 1.1d0; gpkgmax = 1.15d0; gpkgexct = 0.8d0
    end if
case(14)
    !---Se
    if( .not.lreadpp  ) lkbppi = .true.
    if( .not.lreadnlc )  lking = .true.;  rking = 1.5d0;  gkgmax = 1.25d0;  gkgexct = 0.5d0
    if( .not.lreadlc  ) llking =.false.; rlking = 1.5d0; glkgmax = 1.15d0; glkgexct = 0.8d0
    if( .not.lreadpcc ) then
        lpcci = .true.; rpcc = 2.0d0
        lpking= .true.; rpking = 1.1d0; gpkgmax = 1.15d0; gpkgexct = 0.8d0
    end if
case(34)
    !---Se
    if( .not.lreadpp  ) lkbppi = .true.
    if( .not.lreadnlc )  lking = .true.;  rking = 1.5d0;  gkgmax = 1.25d0;  gkgexct = 0.5d0
    if( .not.lreadlc  ) llking =.false.; rlking = 1.5d0; glkgmax = 1.15d0; glkgexct = 0.8d0
    if( .not.lreadpcc ) then
        lpcci = .true.; rpcc = 1.8d0
        lpking= .true.; rpking = 1.1d0; gpkgmax = 1.15d0; gpkgexct = 0.8d0
    end if
case(42)
    !---Mo
    if( .not.lreadpp  ) lkbppi = .true.
    if( .not.lreadnlc )  lking = .true.;  rking = 1.5d0;  gkgmax = 1.25d0;  gkgexct = 0.8d0
    if( .not.lreadlc  ) llking =.false.; rlking = 1.5d0; glkgmax = 1.15d0; glkgexct = 0.8d0
    if( .not.lreadpcc ) then
        lpcci = .true.; rpcc = 2.0d0
        lpking= .true.; rpking = 1.5d0; gpkgmax = 1.15d0; gpkgexct = 0.8d0
    end if
case default
    call fstop( nfile, myid, nodes, 'no defaults setting for atomic data' )
end select


return
end subroutine




subroutine read_pp( nfile, myid, nodes )
!-----------------------------------------------------------------------
!     read input variables for pseudopotential
!-----------------------------------------------------------------------
use param
use param_atom
use constants
implicit none
integer :: nfile(*), myid, nodes

!-----declare local variables
logical, allocatable, dimension(:,:) ::   lchk_tmp
real*8,  allocatable, dimension(:,:) :: vdrcut_tmp
integer, allocatable, dimension(:,:) :: nrefe_tmp
integer :: status
integer :: it
integer :: l, ifpp


!-----allocate memory for work array
allocate( lchk_tmp(0:mxlx,ntype), vdrcut_tmp(0:mxlx,ntype),  &
& nrefe_tmp(0:mxlx,ntype),  &
& stat=status )
if( status /= 0 ) call fstop( nfile, myid, nodes,  &
&                 'memory allocation error in readpp' )

if( lkbpp )  &
&  call pptmpp( nfile, myid, nodes,  &
& jgga, lrela, ntype, lkbppi, zatom, zv, lmax, lsomax, lclno,  &
& lchk_tmp, mxlx, aname )

if( .not.lvand ) lpaw = .false.
!if( lvand ) then
!    if( lpaw ) then
!        ifpp = 3
!      else
!        ifpp = 2
!    end if
!    call ppvand( nfile, myid, nodes,  &
!& ifpp, jgga, lrela, ntype, lvandi, zatom, zv, lmax, lsomax, lclno,  &
!& lchk_tmp, vdrcut_tmp, nrefe_tmp, mxlx, ms, aname )
!end if


!-----get maximum angular momentum : mxl
mxl = 0
do it = 1, ntype
   mxl = max( mxl, lmax(it) )
end do
!nylmmx = mxl*(mxl+2) + 1


!-----allocate memory related to mxl
call mxl_alloc( nfile, myid, nodes )

do it = 1, ntype
   do l = 0, lmax(it)
      lchk(l,it) =   lchk_tmp(l,it)
   end do
end do

rctmax = 0.d0
mxref  = 0.d0
do it = 1, ntype
!   if( lvandi(it) ) then
!       do l = 0, lmax(it)
!          vdrcut(l,it) = vdrcut_tmp(l,it)
!          rctmax = max( rctmax, vdrcut(l,it) )
!          nrefe(l,it) = nrefe_tmp(l,it)
!          mxref = max( mxref, nrefe(l,it) )
!       end do
!     else
       do l = 0, lmax(it)
          nrefe(l,it) = 1
          mxref = max( mxref, nrefe(l,it) )
       end do
!   end if
end do


!-----deallocate memory for work array
deallocate( lchk_tmp, vdrcut_tmp, nrefe_tmp, stat=status )


!-----set variables for arrays
call set_nlpp_parameters( nfile, myid, nodes,  &
& lrela, ntype, lchk, lvandi, lvflag, lvand, mxl, nylmmx, mxref, lmax, nrefe )


!-----set variables for DFT+C in module psvariables
!if( lplusC ) call set_dftC_in_psvariables( nfile, myid, nodes,  &
!& lplusC, lplusC_at, plusC_r, plusC_e, mxlx, ntype, mxl )


return
end subroutine




subroutine fstop( nfile, myid, nodes, message )
!-----------------------------------------------------------------------
!    error trap: display message and stop program
!-----------------------------------------------------------------------
use outfile
implicit none
integer :: nfile(*), myid, nodes
character(*) :: message
integer :: ierr

if( myid == 0 ) then
    if(loutfile(1)) write(nfile(1),*) message(1:len_trim(message))
    if(loutfile(2)) write(nfile(2),*) message(1:len_trim(message))
end if
!-----Finalize the parallel environment
!call end_parallel(ierr)
call abort_parallel(ierr)
stop

end subroutine




subroutine fnormalstop( nfile, myid, nodes, message )
!-----------------------------------------------------------------------
!    error trap: display message and stop program
!-----------------------------------------------------------------------
use outfile
implicit none
integer :: nfile(*), myid, nodes
character(*) :: message
integer :: ierr

if( myid == 0 ) then
    if(loutfile(1)) write(nfile(1),*) message(1:len_trim(message))
    if(loutfile(2)) write(nfile(2),*) message(1:len_trim(message))
end if
!-----Finalize the parallel environment
call end_parallel(ierr)
!call abort_parallel(ierr)
stop

end subroutine




module check_alloc_lgimax
!-----------------------------------------------------------------------
! type declaration and initialization of variables in input.f90
!-----------------------------------------------------------------------
implicit none

logical :: lgimax = .true.
save

end module




subroutine set_lgimax_false
!-----------------------------------------------------------------------
use check_alloc_lgimax
implicit none

lgimax = .false.

end subroutine




subroutine check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, sub_name, lstop )
!-----------------------------------------------------------------------
!    error trap: check memory allocation
!
!    if lstop = .true.  & an error occurs, stop program
!    if lstop = .false. & an error occurs, display message
!-----------------------------------------------------------------------
use outfile
use check_alloc_lgimax
implicit none
integer :: nfile(*), myid, nodes
integer :: status
real*8  :: alloc_mem, the_mem
character(*) :: sub_name
logical :: lstop

!-----declare local variables
integer :: i
integer :: u1, u2
real*8, dimension(3) :: dev = (/ 1.d0, 1.d+03, 1.d+06 /)
character(2), dimension(3) :: unit = (/ 'B ', 'KB', 'MB' /)


if( lgimax ) then
    status = abs(status)
    call gimax(status)
end if

statif: if( status == 0 ) then

    !-----update counter for allocated memory
    alloc_mem = alloc_mem + the_mem

        !-----display messages
        u1 = log(max(the_mem,1.d0))/log(10.d0)
        u1 = u1/3 + 1
        if( u1.gt.3 ) u1 = 3
        u2 = log(max(alloc_mem,1.d0))/log(10.d0)
        u2 = u2/3 + 1
        if( u2.gt.3 ) u2 = 3
        do i = 1, 2
        if( loutfile(i) ) then
           write(nfile(i),*)  &
&      nint(  the_mem/dev(u1)), ' ', unit(u1),  &
& ' (', nint(alloc_mem/dev(u2)), ' ', unit(u2), ' ), allocated (',  &
& sub_name(1:len_trim(sub_name)), ')'
        end if
        end do

else statif

    if( lstop ) then
        call fstop( nfile, myid, nodes,  &
& ('memory allocation error in '//sub_name(1:len_trim(sub_name))) )
    else
       if(loutfile(1)) write(nfile(1),*) 'memory allocation error in ',  &
&                            sub_name(1:len_trim(sub_name))
       if(loutfile(2)) write(nfile(2),*) 'memory allocation error in ',  &
&                            sub_name(1:len_trim(sub_name))
    end if

end if statif


return
end subroutine




subroutine check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, sub_name, lstop )
!-----------------------------------------------------------------------
!    error trap: check memory deallocation
!
!    if lstop = .true.  & an error occurs, stop program
!    if lstop = .false. & an error occurs, display message
!-----------------------------------------------------------------------
use outfile
use check_alloc_lgimax
implicit none
integer :: nfile(*), myid, nodes
integer :: status
real*8  :: alloc_mem, the_mem
character(*) :: sub_name
logical :: lstop

!-----declare local variables
integer :: i
integer :: u1, u2
real*8, dimension(3) :: dev = (/ 1.d0, 1.d+03, 1.d+06 /)
character(2), dimension(3) :: unit = (/ 'B ', 'KB', 'MB' /)


if( lgimax ) then
    status = abs(status)
    call gimax(status)
end if

statif: if( status == 0 ) then

    !-----update counter for allocated memory
    alloc_mem = alloc_mem - the_mem

        !-----display messages
        u1 = log(max(the_mem,1.d0))/log(10.d0)
        u1 = u1/3 + 1
        if( u1.gt.3 ) u1 = 3
        u2 = log(max(alloc_mem,1.d0))/log(10.d0)
        u2 = u2/3 + 1
        if( u2.gt.3 ) u2 = 3
        do i = 1, 2
        if( loutfile(i) ) then
           write(nfile(i),*)  &
&      nint( -the_mem/dev(u1)), ' ', unit(u1),  &
& ' (', nint(alloc_mem/dev(u2)), ' ', unit(u2),  &
& ' ), deallocated (',  &
& sub_name(1:len_trim(sub_name)), ')'
        end if
        end do

else statif

    if( lstop ) then
        call fstop( nfile, myid, nodes,  &
&('memory deallocation error in '//sub_name(1:len_trim(sub_name))))
    else
       if(loutfile(1)) write(nfile(1),*) 'memory deallocation error in ',  &
&                            sub_name(1:len_trim(sub_name))
       if(loutfile(2)) write(nfile(2),*) 'memory deallocation error in ',  &
&                            sub_name(1:len_trim(sub_name))
    end if

end if statif


return
end subroutine




subroutine check_alloc_accum( nfile, myid, nodes,  &
& status, the_mem, sub_name, lstop )
!-----------------------------------------------------------------------
!    error trap: check memory allocation
!
!    if lstop = .true.  & an error occurs, stop program
!    if lstop = .false. & an error occurs, display message
!-----------------------------------------------------------------------
use param
implicit none
integer :: nfile(*), myid, nodes
integer :: status
real*8  :: the_mem
character(*) :: sub_name
logical :: lstop


call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, sub_name, lstop )


end subroutine




subroutine check_dealloc_accum( nfile, myid, nodes,  &
& status, the_mem, sub_name, lstop )
!-----------------------------------------------------------------------
!    error trap: check memory deallocation
!
!    if lstop = .true.  & an error occurs, stop program
!    if lstop = .false. & an error occurs, display message
!-----------------------------------------------------------------------
use param
implicit none
integer :: nfile(*), myid, nodes
integer :: status
real*8  :: the_mem
character(*) :: sub_name
logical :: lstop


call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, sub_name, lstop )


return
end subroutine




subroutine get_alloc_mem( alloc_mem_ )
!-----------------------------------------------------------------------
use param
implicit none
real*8  :: alloc_mem_

alloc_mem_ = alloc_mem

return
end subroutine




subroutine put_alloc_mem( alloc_mem_ )
!-----------------------------------------------------------------------
use param
implicit none
real*8  :: alloc_mem_

alloc_mem = alloc_mem_

return
end subroutine
