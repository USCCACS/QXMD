



subroutine remd_input( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!    read and check input data
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: ierror


!----- read input variables from an input file 
!----- & allocate memory for variables in module remd_param_atom
call remd_read_data( nfile, myid, nodes )


!----- check input variables
call remd_check_data( nfile, myid, nodes, ierror )

if( ierror /= 0 ) return

!----- write input variables
call remd_write_data( nfile, myid, nodes )


return
end subroutine




subroutine remd_check_data( nfile, myid, nodes, ierror )
!-----------------------------------------------------------------------
!    check input data
!-----------------------------------------------------------------------
use remd_param
use remd_param_atom
use constants
!use symmop
implicit none
integer :: nfile(*), myid, nodes
integer :: ierror

!-----declare local variables
integer :: i, j, it, it2, nzn, ix, ibuf, nonfixed
real*8  :: tx, ty, tz, volsuper, volMD, sum
real*8, dimension(3,3) :: b
logical :: lfrag


!-----the number of parallel nodes
nproc = npx*npy*npz

!----- check input variables & error trap ------------------------------


if( ifmd < 0 ) ifmd = 0
if( ifmd > 6 .and. ifmd /= 10 ) ifmd = 0
!---if( lstart .or. ifmd == 0 ) nstep_ini = 0
if( ifmd == 0 ) nstep_ini = 0
if( liscale .and. iscstp <= 0 ) iscstp = 20

if( ifmd == 0 .or. ifmd >= 10 ) lmdstop = .false.

!-----check temperature
!----- error trap
if( ifmd >= 3 .and. ifmd < 10 .and. treq <= 1.d-30 ) then
    if( myid == 0 ) then
        write(nfile(1),*) 'error for setting temperature:',  &
&                         ' ifmd, treq=', ifmd, treq
        write(nfile(2),*) 'error for setting temperature:',  &
&                         ' ifmd, treq=', ifmd, treq
    end if
    ierror = 10
    return
end if

treq = treq/tempau
!            [K] -> [a.u.]


!-----check pressure
if( ifmd == 4 .or. ifmd == 10 .or. ioptmze_cell >= 0 ) then
    hpext  = hpext/prau       ! external pressure
    blkmod = blkmod/prau      ! bulk modulus
!                [GPa] -> [a.u.]

    CALL RCIPRL( h(1,1,0), b, volume )
    volum0 = volume           ! set initial volume

    lstress = .true.
    nskip_stress = 1
    if( .not.lQMMDrun .and. .not.lpureMD ) nskip_stress_out = 1
end if


!-----check shock wave velocity
if( ifmd == 10 ) then
    shockspeed = shockspeed * sqrt(welm*1.d-3/(hrdev*evdj))
!                [m/s] -> [a.u.]
!  T = L*sqrt(M/E) -> T/L = sqrt(M/E)

    call msst_setmatrix( nfile, myid, nodes,  &
& nshockv, h(1,1,0), alpmatrix, betmatrix, trsmatrix, ierror )

    if( ierror /= 0 ) return

    !---set initial value of cell vectors
    h0(1:3,1:3) = h(1:3,1:3,0)

    b(:,:) = hcell(:,:)
    do i = 1, 3
    do j = 1, 3
       sum = 0.d0
       do ix = 1, 3
          sum = sum + trsmatrix(i,ix)*b(ix,j)
       end do
       hcell(i,j) = sum
    end do
    end do
end if

!------gravitational field
!if( lgravi ) then
!    !---magnitude of gravitational field
!    gravg0 = 9.8d0
!!                  in [m/s^2]
!    gravg0 = gravg0 * audang*1.d-10*welm*1.d-3/(hrdev*evdj)
!!                [m/s^2] -> [a.u.]
!!  T = L*sqrt(M/E) -> T^2/L = L*M/E
!
!    !---direction of gravitational field
!    gravdir(1:3) = gravdir(1:3)*gravmag*gravg0
!end if


!----- set No. of QM atoms
if( ntype > 0 ) then
    nhk1(1) = 1
    nhk2(1) = nhk1(1) + nhk(1) - 1
    do i = 2, ntype
       nhk1(i) = nhk2(i-1) + 1
       nhk2(i) = nhk1(i) + nhk(i) - 1
    end do
end if


!----- set No. of MD-cluster atoms
if( ntmd_cluster > 0 ) then
    nhk1_cluster(1) = 1
    nhk2_cluster(1) = nhk1_cluster(1) + ntot_cluster(1) - 1
    do i = 2, ntmd_cluster
       nhk1_cluster(i) = nhk2_cluster(i-1) + 1
       nhk2_cluster(i) = nhk1_cluster(i) + ntot_cluster(i) - 1
    end do
    ntot_cluster(0) = 0
    do it = 1, ntmd_cluster
       ntot_cluster(0) = ntot_cluster(0) + ntot_cluster(it)
    end do
end if


!-----set the number of MD atoms : <--- the same as subroutine setup
do it = 1, ntmd
   ntot(it) = 0
   do i = 1, nmd
      if( is(i) == it ) ntot(it) = ntot(it) + 1
   end do
   call gisum(ntot(it),1,ibuf)
end do
ntot(0) = 0
do it = 1, ntmd
   ntot(0) = ntot(0) + ntot(it)
end do


!----- set atomic mass

!----- for QM atoms
do it = 1, ntype
   nzn = nint(zatom(it))
   !----- error trap
   if( nzn <= 0 .or. nzn > 103 ) then
       if( myid == 0 ) then
           write(nfile(1),*) 'nzn <= 0 .or. nzn > 103 :', nzn
           write(nfile(2),*) 'nzn <= 0 .or. nzn > 103 :', nzn
       end if
       ierror = 11
       return
   end if
   watom(it) = dmassn(nzn)/avogad/welm
end do


!----- for MD atoms
do it = 1, ntmd
   nzn = nint(zatom_md(it))
   !----- error trap
   if( nzn <= 0 .or. nzn > 103 ) then
       if( myid == 0 ) then
           write(nfile(1),*) 'nzn <= 0 .or. nzn > 103 :', nzn
           write(nfile(2),*) 'nzn <= 0 .or. nzn > 103 :', nzn
       end if
       ierror = 11
       return
   end if
   watom_md(it) = dmassn(nzn)/avogad/welm * ficmass(it)
end do


!----- for MD cluster atoms
do it = 1, ntmd_cluster
   nzn = nint(zatom_cluster(it))
   !----- error trap
   if( nzn <= 0 .or. nzn > 103 ) then
       if( myid == 0 ) then
           write(nfile(1),*) 'nzn <= 0 .or. nzn > 103 :', nzn
           write(nfile(2),*) 'nzn <= 0 .or. nzn > 103 :', nzn
       end if
       ierror = 11
       return
   end if
   watom_cluster(it) = dmassn(nzn)/avogad/welm
end do


!-----set translation table from is to ist
ntmd_ist = 0
do it = 1, ntmd
   nzn = nint(zatom_md(it))
   lfrag = .false.
   do it2 = 1, ntmd_ist
      if( nzn == nint(zatom_md_ist(it2)) ) then
          lfrag = .true.
          is2ist(it) = it2
      end if
   end do
   if( .not.lfrag ) then
       ntmd_ist = ntmd_ist + 1
       is2ist(it) = ntmd_ist
       zatom_md_ist(ntmd_ist) = zatom_md(it)
   end if
end do

!----- for MD cluster atoms
do it = 1, ntmd_cluster
   nzn = nint(zatom_cluster(it))
   lfrag = .false.
   do it2 = 1, ntmd_ist
      if( nzn == nint(zatom_md_ist(it2)) ) then
          lfrag = .true.
          is2ist_cluster(it) = it2
      end if
   end do

   !----- error trap
   if( .not.lfrag ) then
       if( myid == 0 ) then
           write(nfile(1),*) 'incorrect species in MD cluster'
           write(nfile(2),*) 'incorrect species in MD cluster'
       end if
       ierror = 12
       return
   end if
end do


!-----check input variables and allocate memory for thermostat variables
if( ifmd >= 3 ) then
    !-----check variables
    if( nnos  < 1 ) nnos  = 1
    if( nresn < 1 ) nresn = 1
    if( nyosh < 1 ) nyosh = 1
    if(  tomega <= 0.d0 )  tomega = 100.d0*dtmd
    if( tbomega <= 0.d0 ) tbomega = 100.d0*dtmd
    if( blkmod  <= 0.d0 )  blkmod = 250.d0/prau
else
    nnos  = 1
    nresn = 1
    nyosh = 1
end if
call remd_thermostat_alloc( nfile, myid, nodes,  &
& alloc_mem, nnos, nresn, nyosh )

!if( ifmd == 5 ) then
!    !-----check variables
!    if( nnos_atom  < 1 ) nnos_atom  = 1
!    if( nresn_atom < 1 ) nresn_atom = 1
!    if( nyosh_atom < 1 ) nyosh_atom = 1
!    if( tomega_atom <= 0.d0 ) tomega_atom = 100.d0*dtmd
!    call remd_atom_thermostat_alloc( nfile, myid, nodes,  &
!& alloc_mem, nnos_atom, nresn_atom, nyosh_atom, nmd_alloc )
!else
    nnos_atom  = 1
    nresn_atom = 1
    nyosh_atom = 1
    call remd_atom_thermostat_alloc( nfile, myid, nodes,  &
& alloc_mem, nnos_atom, nresn_atom, nyosh_atom, 1 )
!end if

!---set GLE thermostat
!if( ifmd == 6 ) then
!    call set_gle( nfile, myid, nodes,  &
!& dtmd, treq, nmd, lmomzero, ierror )
!    if( ierror /= 0 ) return
!end if

!----- set nstop for Harmonic-mode analysis
!if( ifmd == 1 .and. ioptmze == 10 ) then
!    if( lcentdif ) then
!        nstop = 3*nmd * 2*nhmaord
!    else
!        nstop = 3*nmd * nhmaord
!    end if
!    call set_phnlstart( nfile, myid, nodes, lstart )
!    lhma  = .true.
!    lstart = .false.
!end if

!----- set data for phonon-dispersion calculation
!if( ifmd == 1 .and. ioptmze == 11 ) then
!    call set_phnstop( nfile, myid, nodes,  &
!& nstop, nmd, nhmaord, lcentdif, h(1,1,0) )
!    call set_phnlstart( nfile, myid, nodes, lstart )
!    lhma  = .true.
!    lstart = .false.
!end if

!----- set data for calculations with given atomic configurations
!if( ifmd == 1 .and. ioptmze == 20 ) then
!    call set_nstop_from_remd_disp( nfile, myid, nodes, nstop )
!    lhma  = .true.
!!    lstart = .false.
!end if


!----- set irstrct when lsymop = .true.
!if( lsymop .and. irstrct == 0 ) irstrct = 10


!----- set lmdout
lmdout = locoor .or. lovelo .or. loforc

if( .not.lQMMDrun ) then
    locoor_qm = .false.
    lovelo_qm = .false.
    loforc_qm = .false.
    locoor_cl = .false.
    lovelo_cl = .false.
    loforc_cl = .false.
end if
!-----velocity file: not supported yet
lovelo_qm = .false.
lovelo_cl = .false.

lmdout_qm = locoor_qm .or. lovelo_qm .or. loforc_qm
lmdout_cl = locoor_cl .or. lovelo_cl .or. loforc_cl


if( .not.lQMMDrun .and. .not.lpureMD ) then
    latomic = .false.
    lheat   = .false.
end if


!-----check IO parameters
nskip_atomic = ioskip
ioskip_cl = ioskip
ioskip_qm = ioskip
!if( lQMMDrun .or. lpureMD ) then
!    if( lstress ) then
!        if( mod(nskip_atomic,nskip_stress) /= 0 ) then
!            nskip_stress = nskip_atomic
!        end if
!    end if
!    if( lheat ) then
!        if( .not.lstress ) then
!            lstress = .true.
!            nskip_stress = nskip_heat
!        else
!            if( nskip_stress /= 1 ) then
!                do i = 1, min(nskip_atomic, nskip_heat)
!                   if( mod(nskip_atomic,i) == 0 .and. mod(nskip_heat,i) == 0 ) then
!                       nskip_stress = i
!                   end if
!                end do
!            end if
!        end if
!    end if
!    if( lstress ) then
!        if( nskip_stress_out == 0 ) then
!!            !---originally, lstress = .false.
!!            nskip_stress_out = nstop + 1
!        else
!            if( mod(nskip_stress_out,nskip_stress) /= 0 ) then
!                nskip_stress_out = max( (nskip_stress_out/nskip_stress)*nskip_stress, 1 )
!            end if
!        end if
!    end if
!end if
lrestrict_area = xio_min < xio_max .or. yio_min < yio_max .or. &
&                zio_min < zio_max
if( lrestrict_area ) then
    if( xio_min > xio_max ) then
        xio_min = 0.d0 - 1.d-06
        xio_max = 1.d0 + 1.d-06
    end if
    if( yio_min > yio_max ) then
        yio_min = 0.d0 - 1.d-06
        yio_max = 1.d0 + 1.d-06
    end if
    if( zio_min > zio_max ) then
        zio_min = 0.d0 - 1.d-06
        zio_max = 1.d0 + 1.d-06
    end if
end if


!----- check MD/super-cell vectors
!      CALL RCIPRL( h(1,1,0), b, volMD )
!      CALL RCIPRL( hcell, b, volsuper )
!      if( volsuper > volMD )
!&         call fstop( nfile, myid, nodes,  &
!&                     'error: MD cell must be larger than super cell' )


!-----<input data>
!-----ratm : variables for QM atoms,
!-----       must be real coordinates or scaled coordinates corresponding to 'h'
!-----x    : variables for MD atoms,
!-----       must be real coordinates or scaled coordinates corresponding to 'h'
!-----
!-----<in main routine>
!-----ratm : scaled coordinates for QM atoms corresponding to 'h'
!-----x    : scaled coordinates for MD atoms corresponding to 'h'

!----- transpose matrix of b = inverse of h
CALL RCIPRL( h(1,1,0), b, volume )

!-----if real coordinates, get scaled coordinates

!-----QM atoms
do it = 1, ntype
   if( icscale(it) /= 1 ) then
       do i = nhk1(it), nhk2(it)
          tx = ratm(1,i)
          ty = ratm(2,i)
          tz = ratm(3,i)
          ratm(1,i) = b(1,1)*tx + b(2,1)*ty + b(3,1)*tz
          ratm(2,i) = b(1,2)*tx + b(2,2)*ty + b(3,2)*tz
          ratm(3,i) = b(1,3)*tx + b(2,3)*ty + b(3,3)*tz
       end do
   end if
end do


!-----MD atoms
do it = 1, ntmd
   if( icscale_md(it) /= 1 ) then
       do i = 1, nmd
          if( is(i) == it ) then
              tx = x(3*i-2,0)
              ty = x(3*i-1,0)
              tz = x(3*i-0,0)
              x(3*i-2,0) = b(1,1)*tx + b(2,1)*ty + b(3,1)*tz
              x(3*i-1,0) = b(1,2)*tx + b(2,2)*ty + b(3,2)*tz
              x(3*i-0,0) = b(1,3)*tx + b(2,3)*ty + b(3,3)*tz
          end if
       end do
   end if

!         if( .not.vrandom(it) .and. .not.lfixion(it) ) then
!             do i = 1, nmd
!                if( is(i) == it ) then
!                    tx = x(3*i-2,1)
!                    ty = x(3*i-1,1)
!                    tz = x(3*i-0,1)
!                    x(3*i-2,1) = b(1,1)*tx + b(2,1)*ty + b(3,1)*tz
!                    x(3*i-1,1) = b(1,2)*tx + b(2,2)*ty + b(3,2)*tz
!                    x(3*i-0,1) = b(1,3)*tx + b(2,3)*ty + b(3,3)*tz
!                end if
!             end do
!         end if

   if( lfixion(it) ) then
       vrandom(it)  = .false.
       watom_md(it) = 1.d+12
       do i = 1, nmd
          if( is(i) == it ) then
              x(3*i-2,1) = 0.d0
              x(3*i-1,1) = 0.d0
              x(3*i-0,1) = 0.d0
          end if
       end do
   end if
   if( lfixion(it) .or. ifmd <= 1 ) then
       ldispion(it)  = .false.
   end if
end do

!-----check temperature if lfixion = .true.
nonfixed = 0
do it = 1, ntmd
   if( .not.lfixion(it) ) nonfixed = nonfixed + ntot(it)
end do
if( nonfixed /= ntot(0) ) then
    treq = treq * dble(nonfixed)/dble(ntot(0))
    if( ifmd >= 2 ) then
    if( myid == 0 ) then
        do i = 1, 2
           write(nfile(i),*) ( '*', it = 1, 50 )
           write(nfile(i),*) ' No. of fixed/non-fixed ions    :',  &
&                            ntot(0)-nonfixed, nonfixed
           write(nfile(i),*) ' Equilibrium T has been changed :',  &
&                            treq * tempau
           write(nfile(i),*) ( '*', it = 1, 50 )
        end do
    end if
    end if
end if


!-----MD cluster atoms
do it = 1, ntmd_cluster
   if( icscale_cluster(it) /= 1 ) then
       do i = 1, nmd_cluster
          if( is_cluster(i) == it ) then
              tx = x_cluster(3*i-2,0)
              ty = x_cluster(3*i-1,0)
              tz = x_cluster(3*i-0,0)
            x_cluster(3*i-2,0) = b(1,1)*tx + b(2,1)*ty + b(3,1)*tz
            x_cluster(3*i-1,0) = b(1,2)*tx + b(2,2)*ty + b(3,2)*tz
            x_cluster(3*i-0,0) = b(1,3)*tx + b(2,3)*ty + b(3,3)*tz
          end if
       end do
   end if
end do


!--- error trap
!--- check configuration
!-----for MD atoms
do i = 1, nmd
   if( x(3*i-2,0) < 0.d0 .or. x(3*i-2,0) >= anxi .or.  &
&      x(3*i-1,0) < 0.d0 .or. x(3*i-1,0) >= anyi .or.  &
&      x(3*i-0,0) < 0.d0 .or. x(3*i-0,0) >= anzi )  &
&      ierror = 901
end do

!-----for MD cluster atoms
do i = 1, nmd_cluster
   if( x_cluster(3*i-2,0) <  0.d0 .or.  &
&      x_cluster(3*i-2,0) >= 1.d0 .or.  &
&      x_cluster(3*i-1,0) <  0.d0 .or.  &
&      x_cluster(3*i-1,0) >= 1.d0 .or.  &
&      x_cluster(3*i-0,0) <  0.d0 .or.  &
&      x_cluster(3*i-0,0) >= 1.d0       )  &
&      ierror = 902
end do

!-----for QM atoms
do i = 1, natom
   if( ratm(1,i) < 0.d0 .or. ratm(1,i) >= 1.d0 .or.  &
&      ratm(2,i) < 0.d0 .or. ratm(2,i) >= 1.d0 .or.  &
&      ratm(3,i) < 0.d0 .or. ratm(3,i) >= 1.d0 )  &
&      ierror = 903
end do

call gimax(ierror)
if( myid == 0 ) then
    do i = 1, 2
       if( ierror == 901 ) then
           write(nfile(i),*) 'error in coordinates of MD atoms'
       end if
       if( ierror == 902 ) then
           write(nfile(i),*) 'error in coordinates of CL atoms'
       end if
       if( ierror == 903 ) then
           write(nfile(i),*) 'error in coordinates of QM atoms'
       end if
    end do
end if


!-----for virtual MD for thermodynamic integration
!if( lvmd ) then
!    leinstein = .false.
!    if( lidealref ) then
!        do it = 1, ntmd
!           leinstein = leinstein .or. nvmd(it) == 2
!        end do
!    end if
!else
    lidealref = .true.
    leinstein = .false.
    if( .not.allocated(x_vmd) ) allocate( x_vmd(3) )
!end if


return
end subroutine




subroutine remd_write_data( nfile, myid, nodes )
!-----------------------------------------------------------------------
!    write input data
!-----------------------------------------------------------------------
use remd_param
use remd_param_atom
use constants
implicit none
integer :: nfile(*), myid, nodes

!-----declare local variables
integer :: i, it, l, j, ix


myidif: if( myid == 0 ) then

   do i = 1, 2
      write(nfile(i),*) ' '
!   if( lQMMDrun ) then
!      write(nfile(i),*) ' ----      Hybrid QM/MD calcula',  &
!&                       'tion with PW method      ----'
!   else if( lpureMD ) then
!      write(nfile(i),*) ' ---- Classical MD calculation ',  &
!&                       'with empirical potential ----'
!   else
      write(nfile(i),*) ' ----      Ab initio MD calcula',  &
&                       'tion with PW method      ----'
!   end if
!   if( lvmd ) then
!      write(nfile(i),*) ' ----     virtual MD for thermo',  &
!&                       'dynamic integration      ----'
!   end if
!   if( lQMMDrun .or. lpureMD ) then
!       call out_potential( nfile(i) )
!   end if

      write(nfile(i),*) ' '
      write(nfile(i),2606) '    lstart :', lstart,  &
&                          '    iogpsz :', iogpsz
      write(nfile(i),2000) '      ntmd :', ntmd,  &
&                          '     ntype :', ntype,  &
&                          ' ntmd_clus :', ntmd_cluster
      write(nfile(i),2000) '  ntmd_ist :', ntmd_ist
      if( ntot(0) < 1000000 ) then
          write(nfile(i),'(1x,a12,i6,6x,a12,i6,6x,a12,i6)')  &
&                          '   ntot(0) :', ntot(0),  &
&                          '     natom :', natom,  &
&                          '  nmd_clus :', nmd_cluster
      else
          write(nfile(i),'(1x,a12,i12,  a12,i6,6x,a12,i6)')  &
&                          '   ntot(0) :', ntot(0),  &
&                          '     natom :', natom,  &
&                          '  nmd_clus :', nmd_cluster
      end if
      write(nfile(i),2606) '    lmdout :', lmdout,  &
&                          '    ioskip :', ioskip,  &
&                          'lrestrict_a:', lrestrict_area
!      if( lrestrict_area ) then
!          write(nfile(i),2555) '   xio_min :', xio_min,  &
!&                              '   xio_max :', xio_max,  &
!&                              '   yio_min :', yio_min
!          write(nfile(i),2555) '   yio_max :', yio_max,  &
!&                              '   zio_min :', zio_min,  &
!&                              '   zio_max :', zio_max
!      end if
      if( lmdout ) then
          write(nfile(i),2666) '    locoor :', locoor,  &
&                              '    lovelo :', lovelo,  &
&                              '    loforc :', loforc
          write(nfile(i),'(1x,a12,i6,6x,a12,i6,6x,a12,i6)')  &
&                              'ioskipcoor :', ioskipcoor,  &
&                              'ioskipvelo :', ioskipvelo,  &
&                              'ioskipforc :', ioskipforc
      end if
!      if( lQMMDrun ) then
!         write(nfile(i),2600) ' lmdout_qm :', lmdout_qm,  &
!&                             ' ioskip_qm :', ioskip_qm
!         if( lmdout_qm ) then
!            write(nfile(i),2666) ' locoor_qm :', locoor_qm,  &
!&                                ' lovelo_qm :', lovelo_qm,  &
!&                                ' loforc_qm :', loforc_qm
!         end if
!         write(nfile(i),2600) ' lmdout_cl :', lmdout_cl,  &
!&                             ' ioskip_cl :', ioskip_cl
!         if( lmdout_cl ) then
!            write(nfile(i),2666) ' locoor_cl :', locoor_cl,  &
!&                                ' lovelo_cl :', lovelo_cl,  &
!&                                ' loforc_cl :', loforc_cl
!         end if
!      end if
      if( ifmd >= 1 ) then
         if( nstop < 1000000 ) then
             write(nfile(i),2050) '      ifmd :', ifmd,  &
&                                 '      dtmd :', dtmd,  &
&                                 '     nstop :', nstop
         else
             write(nfile(i),'(1x,a12,i6,6x,a12,f10.2,2x,a12,i11)')  &
&                                 '      ifmd :', ifmd,  &
&                                 '      dtmd :', dtmd,  &
&                                 '     nstop :', nstop
         end if
         write(nfile(i),2666) '   lmdstop :', lmdstop
         if( .not.lstart ) then
             write(nfile(i),2050) ' nstep_ini :', nstep_ini
         end if
         write(nfile(i),2060) '   ioptmze :', ioptmze,  &
&                             '      lpvv :', lpvv
         if( ifmd == 1 ) then
             if( ioptmze >= 2 .and. ioptmze <= 9 ) then
                 write(nfile(i),2400) '  dist_max :', dist_max
!             else if( ioptmze == 10 ) then
!                 write(nfile(i),'(1x, a12,f10.5,2x, a12,i6,6x, a12,l6)')  &
!&                                     '   hmadisp :', hmadisp, &
!&                                     '   nhmaord :', nhmaord, &
!&                                     '  lcentdif :', lcentdif
             end if
             if( ioptmze == 2 ) then
                 write(nfile(i),'(1x, a12,i6,6x,a12,es11.4,1x,a12,es11.4)')  &
&                             'ibfgsclear :', ibfgsclear,  &
&                             '   qnstabi :', qnstabi,  &
&                             '  gammamin :', gammamin
             end if
         end if
         write(nfile(i),2060) 'ioptmze_cel:', ioptmze_cell
         if( ifmd == 1 ) then
             if( ioptmze_cell == 0 ) then
                 write(nfile(i),2400) '  dtcellcg :', dtcellcg
             end if
             if( ioptmze_cell == 2 ) then
                 write(nfile(i),2060) 'iclearcellh:', iclearcellh
             end if
         end if
         if( ifmd == 1 ) then
             if( ioptmze >= 0 .and. ioptmze <= 9 .and. ioptmze_cell >= 0 ) then
                 write(nfile(i),'(1x,a12,l6,5x,a13,i6,1x,a17,i6)')  &
&                                     'lhybridopt :', lhybridopt,  &
&                                    'nstep_hybrid:', nstep_hybrid,  &
&                                'nstep_hybrid_cel:', nstep_hybrid_cell
             end if
         end if
         write(nfile(i),2600) '   liscale :', liscale,  &
&                             '    iscnum :', iscnum,  &
&                             '    iscstp :', iscstp
         write(nfile(i),2600) '  lmomzero :', lmomzero
         write(nfile(i),2500) '      treq :', treq * tempau
         if( ifmd == 4 .or. ifmd == 10 .or. ioptmze_cell >= 0 ) then
              write(nfile(i),2500) '     hpext :', hpext * prau
         end if
!         write(nfile(i),'(1x,a12,i6,6x,a12,i6,5x,a13,3l4)')  &
!&                          '   irstrct :', irstrct,  &
!&                          'irstrct_sub:', irstrct_sub,  &
!&                         'lcell_rstrct:', lcell_rstrct(1:3)
         if( ifmd >= 3 ) then
             if( ifmd < 10 ) then
                write(nfile(i),2000) '      nnos :', nnos,  &
&                                    '     nresn :', nresn,  &
&                                    '     nyosh :', nyosh
             end if
             if( ifmd == 3 ) then
                 write(nfile(i),2500) '    tomega :', tomega
             else if( ifmd == 4 ) then
                 write(nfile(i),2550) '    tomega :', tomega,  &
&                                     '   tbomega :', tbomega
                 write(nfile(i),2550) '    blkmod :', blkmod * prau
!             else if( ifmd == 5 ) then
!                 write(nfile(i),2500) '    tomega :', tomega
!                 write(nfile(i),2000) ' nnos_atom :', nnos_atom,  &
!&                                     'nresn_atom :', nresn_atom,  &
!&                                     'nyosh_atom :', nyosh_atom
!                 write(nfile(i),2500) 'tomega_atom:', tomega_atom
             else if( ifmd == 10 ) then
                 write(nfile(i),'(1x, a12,f10.2,2x, a12,f10.2,2x, a12,es11.4)')  &
&                                     '   tbomega :', tbomega,  &
&                                     '    blkmod :', blkmod * prau,  &
&                                     '    vxmsst :', vxmsst
                 write(nfile(i),2600) 'lmsstscale :', lmsstscale,  &
&                                     ' msstscnum :', msstscnum,  &
&                                     ' msstscstp :', msstscstp
             end if
         end if
         write(nfile(i),2660) '     lstat :', lstat,  &
&                             '   latomic :', latomic,  &
&                             'nskip_atomi:', nskip_atomic
!         write(nfile(i),2600) '     lheat :', lheat,  &
!&                             'nskip_heat :', nskip_heat
         if( ifmd == 1 ) then
             write(nfile(i),'(1x, a12,es11.4,1x, a12,es11.4,1x, a12,i6)')  &
&                                 'tol_energy :', tol_energy,  &
&                                 ' tol_force :', tol_force
         end if
      end if
      write(nfile(i),2600) '   lstress :', lstress,  &
&                          ' nskip_stre:', nskip_stress,  &
&                          'nskp_strs_o:', nskip_stress_out
      write(nfile(i),2600) '   ltotmom :', ltotmom,  &
&                          ' nskip_totm:', nskip_totmom
      write(nfile(i),2660) '     lsave :', lsave,  &
&                          '   lsreal8 :', lsreal8
      write(nfile(i),'(1x, a12,f11.4,2x, a12,f11.4,2x, a12,i6)')  &
&                          ' rc_buffer :', rc_buffer
      write(nfile(i),2000) ' nmd_alloc :', nmd_alloc,  &
&                          'nmd_buffer :', nmd_buffer
!      write(nfile(i),2000) '    nwalls :', nwalls
!      if( lgravi ) then
!          write(nfile(i),'(1x,a12,l6,6x,a12,es12.5,a12,es12.5)')  &
!&                          '    lgravi :', lgravi,  &
!&                          '   gravmag :', gravmag,  &
!&                          ' in [a.u.] :', gravmag*gravg0
!          write(nfile(i),'(1x,a12,3f10.6)')  &
!&                          '    gravdir :', gravdir(1:3)/gravmag/gravg0
!      else
!          write(nfile(i),'(1x,a12,l6,6x,a12,f11.2,1x,a12,es12.6)')  &
!&                          '    lgravi :', lgravi
!      end if
!      if( lefield ) then
!          write(nfile(i),'(1x,a12,l6,6x,a12,es12.5,a12,es12.5)')  &
!&                          '   lefield :', lefield
!          write(nfile(i),'(1x,a12,3f10.6)')  &
!&                          'efield [hr]:', efield(1:3)
!      else
!          write(nfile(i),'(1x,a12,l6,6x,a12,f11.2,1x,a12,es12.6)')  &
!&                          '   lefield :', lefield
!      end if

      write(nfile(i),*) ' '
      write(nfile(i),3000) '   h(MD) (L_1) :', ( h(j,1,0),j=1,3)
      write(nfile(i),3000) '         (L_2) :', ( h(j,2,0),j=1,3)
      write(nfile(i),3000) '         (L_3) :', ( h(j,3,0),j=1,3)
      if( .not.lpureMD ) then
      write(nfile(i),3000) '   hcell (L_1) :', ( hcell(j,1),j=1,3)
      write(nfile(i),3000) '    (QM) (L_2) :', ( hcell(j,2),j=1,3)
      write(nfile(i),3000) '         (L_3) :', ( hcell(j,3),j=1,3)
      end if

      if( ifmd == 10 ) then
          write(nfile(i),*) ' '
          write(nfile(i),'(1x,a,es15.6,f12.2)')  &
& 'shock speed (a.u., m/s) :', shockspeed, shockspeed/sqrt(welm*1.d-3/(hrdev*evdj))
          write(nfile(i),'(1x,a,3i5)')  &
& '           nshockv(1:3) :', nshockv(1:3)
          write(nfile(i),'(1x,a,3es15.6)')  &
& '              alpmatrix :', alpmatrix(1,1:3),  &
& '                         ', alpmatrix(2,1:3),  &
& '                         ', alpmatrix(3,1:3)
          write(nfile(i),'(1x,a,3es15.6)')  &
& '              betmatrix :', betmatrix(1,1:3),  &
& '                         ', betmatrix(2,1:3),  &
& '                         ', betmatrix(3,1:3)
      end if
! call fstop( nfile, myid, nodes, 'write_data' )

!      if( ifmd == 1 .and. ioptmze == 11 ) then
!          write(nfile(i),*) ' '
!          write(nfile(i),*) '* Phonon-dispersion calculation'
!          write(nfile(i),'(1x, a12,f10.5,2x, a12,i6,6x, a12,l6)')  &
!&                              '   hmadisp :', hmadisp, &
!&                              '   nhmaord :', nhmaord, &
!&                              '  lcentdif :', lcentdif
!          call out_phdata( nfile(i) )
!      end if

!      if( ifmd == 1 .and. ioptmze == 20 ) then
!          write(nfile(i),*) ' '
!          write(nfile(i),*) '* Calculations with given atomic displacement'
!          call out_dispdata( nfile(i) )
!      end if

      !-----for MD atoms
      write(nfile(i),*) ' '
      write(nfile(i),3200)
      do it = 1, ntmd
         write(nfile(i),3210) it, aname(nint(zatom_md(it))),  &
&         zatom_md(it), watom_md(it)*avogad*welm, ntot(it),  &
&         icscale_md(it), vrandom(it), lfixion(it), is2ist(it)
      end do
      write(nfile(i),3220)
      do it = 1, ntmd
         if( ldispion(it) ) then
             write(nfile(i),3230) it, aname(nint(zatom_md(it))),  &
&         ficmass(it), ldispion(it), (dvector(j,it),j=1,3)
            else
             write(nfile(i),3230) it, aname(nint(zatom_md(it))),  &
&         ficmass(it), ldispion(it)
         end if
      end do
      write(nfile(i),*) ' (reduced species)'
      do it = 1, ntmd_ist
         write(nfile(i),3210) it, aname(nint(zatom_md_ist(it))),  &
&         zatom_md_ist(it)
      end do

      if( ntype > 0 ) then
          !-----for QM atoms
          write(nfile(i),*) ' '
          write(nfile(i),3300)
          do it = 1, ntype
             write(nfile(i),3210) it, aname(nint(zatom(it))),  &
&             zatom(it), watom(it)*avogad*welm, nhk(it),  &
&             icscale(it), lterminator(it)
          end do
      end if

      if( ntmd_cluster > 0 ) then
          !-----for MD-cluster atoms
          write(nfile(i),*) ' '
          write(nfile(i),3400)
          do it = 1, ntmd_cluster
         write(nfile(i),3210) it, aname(nint(zatom_cluster(it))),  &
&         zatom_cluster(it), watom_cluster(it)*avogad*welm,  &
&         ntot_cluster(it), icscale_cluster(it), lMDterminator(it)
          end do
      end if


!      if( ncbonds > 0 ) then
!          !-----bond-length constraint
!          write(nfile(i),*) ' '
!          write(nfile(i),'(a)') ' * Bond-length constraints'
!          write(nfile(i),'(a)') '        atom1 - atom2      length'
!          do j = 1, ncbonds
!             write(nfile(i),'(i6,a,i6,a,i6,f14.6)') &
!& j, ':', ncbatm1(j), ' -',ncbatm2(j), cblength(j)
!          end do
!      end if

!-----check
!            write(nfile(i),*) ' '
!            write(nfile(i),*) 'coordinates of MD atoms '
!            do it = 1, ntmd
!               write(nfile(i),*) ' it=', it
!               write(nfile(i),*) ' ntot=', ntot(it)
!               write(nfile(i),*) ' coordinates'
!               do j = 1, nmd
!               if( is(j) == it ) then
!               write(nfile(i),'(i5,3f12.8)') j, (x(3*j-3+ix,0), ix =1,3)
!               end if
!               end do
!               write(nfile(i),*) ' velocities'
!               do j = 1, nmd
!               if( is(j) == it ) then
!               write(nfile(i),'(i5,3f12.8)') j, (x(3*j-3+ix,1), ix =1,3)
!               end if
!               end do
!            end do
!
!            if( ntype > 0 ) then
!            write(nfile(i),*) ' '
!            write(nfile(i),*) 'coordinates of QM atoms '
!            do it = 1, ntype
!               write(nfile(i),*) ' it=', it
!               write(nfile(i),*) ' nhk1, nhk2=', nhk1(it), nhk2(it)
!               write(nfile(i),*) ' coordinates'
!               do j = nhk1(it), nhk2(it)
!                  write(nfile(i),'(i5,3f12.8)') j, (ratm(ix,j), ix =1,3)
!               end do
!!               write(nfile(i),*) ' velocities'
!!               do j = nhk1(it), nhk2(it)
!!                  write(nfile(i),'(i5,3f12.8)') j, (vatm(ix,j), ix =1,3)
!!               end do
!            end do
!            end if
!
!            if( ntmd_cluster > 0 ) then
!            write(nfile(i),*) ' '
!            write(nfile(i),*) 'coordinates of MD-cluster atoms '
!            do it = 1, ntmd_cluster
!               write(nfile(i),*) ' it=', it
!               write(nfile(i),*) ' ntot=', ntot_cluster(it)
!               write(nfile(i),*) ' coordinates'
!               do j = 1, nmd_cluster
!               if( is_cluster(j) == it ) then
!                   write(nfile(i),'(i5,3f12.8)')
!&                    j, (x_cluster(3*j-3+ix,0), ix =1,3)
!               end if
!               end do
!!               write(nfile(i),*) ' velocities'
!!               if( is(j) == it ) then
!!                   write(nfile(i),'(i5,3f12.8)')
!!&                    j, (x_cluster(3*j-3+ix,1), ix =1,3)
!!               end if
!            end do
!            end if
!-----check

!      if( nwalls > 0 ) then
!          write(nfile(i),3600)
!          do j = 1, nwalls
!             write(nfile(i),3610) j, ( wallp(ix,j), ix = 1, 3 ),  &
!& ( wallv(ix,j), ix = 1, 3 ), wallf(j), nwallp(j)
!          end do
!      end if


!      if( ifmd == 5 .or. ifmd == 6 ) then
!          !-----NVT for each atom
!          write(nfile(i),*) ' '
!          if( ifmd == 5 ) then
!              write(nfile(i),'(a)') ' * Nose-Hoover thermostats for each atom type'
!          else
!              write(nfile(i),'(a)') ' * GLE thermostats for each atom type'
!          end if
!          write(nfile(i),'(a,l2)') '    lgthermo =', lgthermo
!          write(nfile(i),'(a)') '         lathermo'
!          do it = 1, ntmd
!             if( lathermo(it) ) then
!                 write(nfile(i),'(i2,a,a2,4x,a)') it, ':',  &
!& aname(nint(zatom_md(it))), 'T:respective thermostat'
!             else
!                 write(nfile(i),'(i2,a,a2,4x,a)') it, ':',  &
!& aname(nint(zatom_md(it))), 'F:global     thermostat'
!             end if
!          end do
!!          write(nfile(i),'(a)') '      ( nvmd = 1:ideal gas / 2:Einstein solid )'
!      end if

!      if( ifmd == 6 ) then
!          !-----GLE thermostat
!          write(nfile(i),*) ' '
!
!          call out_GLE_thermostat( nfile(i) )
!      end if

!      if( lvmd ) then
!          !-----virtual MD for thermodynamic integration
!          write(nfile(i),*) ' '
!          write(nfile(i),'(a,a)') ' * Parameters for virtual MD ',  &
!&                                 'for thermodynamic integration'
!          if( lidealref ) then
!              write(nfile(i),'(a,f10.5,l5)') '    dlambda_vmd, leinstein =',  &
!& dlambda_vmd, leinstein
!              write(nfile(i),'(a)') '         nvmd               omega_vmd'
!              do it = 1, ntmd
!                if( nvmd(it) == 1 ) then
!                   write(nfile(i),'(i2,a,a2,4x,a,1x,es13.5)') it, ':',  &
!& aname(nint(zatom_md(it))), '1:ideal gas'
!                else
!                   write(nfile(i),'(i2,a,a2,4x,a,1x,es13.5)') it, ':',  &
!& aname(nint(zatom_md(it))), '2:Einstein solid', omega_vmd(it)
!                end if
!              end do
!!          write(nfile(i),'(a)') '      ( nvmd = 1:ideal gas / 2:Einstein solid )'
!          else
!              call out_vmd_potential( nfile(i) )
!          end if
!      end if

      write(nfile(i),*) ' '

   end do

end if myidif

!
! 0:i6, 2:e11.4, 4:f10.5, 5:f10.2, 6:l6 descriptors
!
2000 format(1x, a12,i6,6x,    a12,i6,6x,    a12,i6)
2060 format(1x, a12,i6,6x,    a12,l6,6x,    a12,i6)
2602 format(1x, a12,l6,6x,    a12,i6,6x,    a12,es11.4)
2002 format(1x, a12,i6,6x,    a12,i6,6x,    a12,es11.4)
2044 format(1x, a12,i6,6x,    a12,f10.5,2x, a12,f10.5)
2400 format(1x, a12,f10.5,2x, a12,i6,6x,    a12,i6)
2050 format(1x, a12,i6,6x,    a12,f10.2,2x, a12,i6)
2600 format(1x, a12,l6,6x,    a12,i6,6x,    a12,i6)
2606 format(1x, a12,l6,6x,    a12,i6,6x,    a12,l6)
2500 format(1x, a12,f10.2,2x, a12,i6,6x,    a12,i6)
2660 format(1x, a12,l6,6x,    a12,l6,6x,    a12,i6)
2666 format(1x, a12,l6,6x,    a12,l6,6x,    a12,l6)
2550 format(1x, a12,f10.2,2x, a12,f10.2,2x, a12,i6)
2555 format(1x, a12,f10.2,2x, a12,f10.2,2x, a12,f10.2)
2640 format(1x, a12,l6,6x,    a12,f10.5,2x, a12,i6)

3000 format(1x,a16,3es16.8)
3010 format(1x,a16,3l5) 
3020 format(1x,a16,3i5) 
3200 format(' * MD atoms'/  &
& 6x,'zatom    watom            ntot icscale vrandom lfixion is2ist' )
3210 format(i2,':',a2,f6.1,es13.5,i12,i5,2l8,i8 )
3220 format(6x,'ficmass ldispion dvector' )
3230 format(i2,':',a2,f8.2,l4,5x,3f7.4)
3300 format(' * QM atoms'/  &
& 6x,'zatom    watom            ntot icscale lterminator' )
3400 format(' * MD-cluster atoms'/  &
& 6x,'zatom    watom            ntot icscale lterminator' )
3600 format(/6x,'wallp                wallv             wallf',  &
& '         nwallp')
3610 format(i4,':',3f6.2,2x,3f6.2,es13.5,i5)


return
end subroutine




subroutine remd_read_data( nfile, myid, nodes )
!-----------------------------------------------------------------------
!     read input variables from an input file
!-----------------------------------------------------------------------
use remd_param
use remd_param_atom
!use symmop
implicit none
integer :: nfile(*), myid, nodes

!-----declare local variables
character(10) :: word
character(50) :: fname1
character(50) :: fname = 'control/filename'
integer :: iunit
integer :: ierror, istat
integer :: ichest, ihest
integer :: i, j, ix, it, na, na1
logical :: lnotread_QMatoms = .true.
logical :: lnotread_MDatoms = .true.
logical :: lnotread_supercell = .true.
logical :: lnotread_MDcell    = .true.
real*8, dimension(3,3) :: b
!real*8, dimension(3) :: disp = 0.d0
real*8, dimension(3) :: r1max, r1min
real*8  :: q1, q2, q3, vv !, vmd_buffer


call allocate_unit_number( iunit )

!--- get file name : fname1
if( myid == 0 ) then
    write(nfile(1),*) 'open file: ', fname(1:len_trim(fname))
    write(nfile(2),*) 'open file: ', fname(1:len_trim(fname))
endif
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
if( myid == 0 ) then
    write(nfile(1),*) 'open file: ', fname1(1:len_trim(fname1))
    write(nfile(2),*) 'open file: ', fname1(1:len_trim(fname1))
endif
open( iunit, file=fname1, status='old', action='read',  &
&     iostat=ierror )

!--- global maximum
call gimax(ierror)


blockif_b1: if( ierror == 0 ) then

   !-----get the MD cell vectors
   cellreaddo: do
      read(iunit,'(a10)',iostat=istat) word
      cellblockif: if( istat == 0 ) then

         cellreadif: if( word == '*MD cell  ' ) then

            call read_cell( nfile, myid, nodes, iunit, h )

            lnotread_MDcell = .false.

         else if( word == '*supercell' ) then cellreadif

            call read_cell( nfile, myid, nodes, iunit, hcell )

            lnotread_supercell = .false.

         end if cellreadif

      else if( istat < 0 ) then cellblockif

         !----- exit, if at EOF
         exit

      else cellblockif

         !----- error trap
         call fstop( nfile, myid, nodes,'error in file :'//fname1)

      end if cellblockif

   end do cellreaddo


   !-----set MD/super-cell vectors
   !----- error trap
   if( lnotread_supercell .and. lnotread_MDcell ) then
    call fstop( nfile, myid, nodes,  &
& 'error in setting MD/super-cell vectors' )
   end if

   if( lnotread_MDcell ) then
       do i = 1, 3
       do ix = 1, 3
          h(ix,i,0) = hcell(ix,i)
       end do
       end do
   end if

   if( lnotread_supercell ) then
       do i = 1, 3
       do ix = 1, 3
          hcell(ix,i) = h(ix,i,0)
       end do
       end do
   end if

   !----- transpose matrix of b = inverse of h
   CALL RCIPRL( h(1,1,0), b, volume )


   !----- rewind file
   rewind(iunit)

   !-----get the numbers of species & atoms
   !----- ntype & natom  for QM atoms
   !----- ntmd  & nmd    for MD atoms
   ntmd = 0
   prereaddo: do
      read(iunit,'(a10)',iostat=istat) word
      preblockif_b2: if( istat == 0 ) then

         prereadif: if( word == '*atoms    ' ) then

            !-----get the numbers of species & atoms : ntype & natom
            !-----for QM atoms
            call read_atom_number( nfile, myid, nodes, iunit,  &
& ntype, natom )

            lnotread_QMatoms = .false.

         else if( word == '*MD atoms ' ) then prereadif

            !-----get displacement vector : disp
            call remd_read_atom_number( nfile, myid, nodes, iunit, &
& ntmd, nmd, sxog, syog, szog, anxi, anyi, anzi, b,  &
& disp, r1max, r1min )

            do ix = 1, 3
               !--- error trap
               if( r1max(ix) - r1min(ix) > 1.d0 )  &
&                  call fstop( nfile, myid, nodes,  &
&               'error in scaled coordinates in remd_read_data' )

               disp(ix) = 0.5d0 *(1.d0 - (r1max(ix) + r1min(ix)) )
            end do

            !----- rewind file
            rewind(iunit)
            prereaddo3: do
                read(iunit,'(a10)',iostat=istat) word
                preblockif_b3: if( istat == 0 ) then

                    if( word == '*MD atoms ' ) exit

                else preblockif_b3

                    !----- error trap
                    call fstop( nfile, myid, nodes,  &
&                           'error in file :'//fname1)

                end if preblockif_b3
            end do prereaddo3


            !-----get the numbers of species & atoms : ntmd & nmd
            !-----for MD atoms
            call remd_read_atom_number( nfile, myid, nodes, iunit, &
& ntmd, nmd, sxog, syog, szog, anxi, anyi, anzi, b,  &
& disp, r1max, r1min )


            !--- error trap
            do ix = 1, 3
               if( r1max(ix) >= 1.d0 .or. r1min(ix) < 0.d0 )  &
&                  call fstop( nfile, myid, nodes,  &
&            'error in scaled coordinates (2) in remd_read_data' )
            end do

            lnotread_MDatoms = .false.

         else if( word == '*MD termin' ) then prereadif

            !-----get the numbers of species & atoms : ntHSMD & nHSMD
            !-----for MD terminator atoms
            call read_atom_number( nfile, myid, nodes, iunit,  &
& ntHSMD, nHSMD )

!         else if( word == '*constrain' ) then prereadif

            !-----for constraint conditions
!            call remd_read_constraint( nfile, myid, nodes, iunit )

!         else if( word == '*virtual m' ) then prereadif
!
!            call remd_preread_vmd( nfile, myid, nodes, iunit )

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
   if( .not.lpureMD .and. lnotread_QMatoms )  &
&      call fstop( nfile, myid, nodes,  &
&                         'no QM atoms in file :'//fname1)
   if( lnotread_MDatoms .and. lnotread_QMatoms )  &
&      call fstop( nfile, myid, nodes,  &
&                         'no QM & MD atoms in file :'//fname1)


   !----- if no MD atoms, this is not the QM/MD run
   !-----    copy   QM atoms as MD atoms
   !-----    ignore MD terminator atoms
   !-----    ignore QM terminator atoms
   if( lnotread_MDatoms ) then
       lQMMDrun = .false.
       ntmd      = ntype
       nmd       = natom
       ntHSMD    = 0
       nHSMD     = 0
   end if

   !----- if purely classical MD run, this is not the QM/MD run
   !-----    ignore MD terminator atoms
   !-----    ignore QM terminator atoms
!   if( lpureMD ) then
!       lQMMDrun = .false.
!       ntHSMD    = 0
!       nHSMD     = 0
!   end if

   rc_buffer = 0d0
   !----- if lvmd == .false. or. lidealref == .true., 
   !-----    set vmd_buffer = 0.d0
!   call vmd_potential_prealloc( nfile, myid, nodes, alloc_mem, ntmd )
!   call set_vmd_rc_buffer( nfile, myid, nodes, vmd_buffer )

!   if( lQMMDrun .or. lpureMD ) then
!       !----- pre-allocate memory for variables for MD potentials
!       call md_potential_prealloc( nfile, myid, nodes, alloc_mem, ntmd )

!       !----- set rc_buffer ( = 0, if not QM/MD run )
!       call set_rc_buffer( nfile, myid, nodes, rc_buffer )

!       rc_buffer = max( rc_buffer, vmd_buffer )
!     else
!       rc_buffer = vmd_buffer
!   end if

   !-----allocate memory for atoms
   call remd_param_atom_alloc( nfile, myid, nodes,  &
& alloc_mem, ntype, natom, ntmd, nmd, nmd_alloc, b,  &
& rc_buffer, anxi, anyi, anzi, ntHSMD, nHSMD )


   !----- rewind file
   rewind(iunit)

   !----- read input data
   readdo: do
      read(iunit,'(a10)',iostat=istat) word
      blockif_b2: if( istat == 0 ) then

         readif: if( word == '*start(on/' .or. word == '*start    ' ) then

            call read_start( nfile, myid, nodes, iunit, lstart )

         else if( word == '*molecular' ) then readif

            call remd_read_md( nfile, myid, nodes, iunit,  &
& ifmd, dtmd, nstop, lmdstop, nstep_ini, &
& ioptmze, lpvv, hmadisp, nhmaord, lcentdif, ibfgsclear, qnstabi, gammamin, &
& ioptmze_cell, dtcellcg, iclearcellh,  &
& lhybridopt, nstep_hybrid, nstep_hybrid_cell,  &
& treq, liscale, iscnum, iscstp, lmomzero, ichest, ihest,  &
& nnos, nresn, nyosh, tomega, hpext, tbomega, blkmod,  &
& nnos_atom, nresn_atom, nyosh_atom, tomega_atom,  &
& irstrct, irstrct_sub, lcell_rstrct,  &
& ioskip, locoor, lovelo, loforc, ioskipcoor, ioskipvelo, ioskipforc,  &
& ioskip_qm, locoor_qm, lovelo_qm, loforc_qm,  &
& ioskip_cl, locoor_cl, lovelo_cl, loforc_cl, lstat,  &
& latomic, nskip_atomic, lheat, nskip_heat,  &
& xio_min, xio_max, yio_min, yio_max, zio_min, zio_max,  &
& tol_energy, tol_force, shockspeed, nshockv,  &
& lmsstscale, msstscnum, msstscstp, vxmsst )

         else if( word == '*save data' ) then readif

            call read_save( nfile, myid, nodes, iunit,  &
& lsave, lsreal8 )

         else if( word == '*stress ca' ) then readif

            call read_str( nfile, myid, nodes, iunit,  &
& lstress, nskip_stress )

            !---store original skip step
            if( lstress ) nskip_stress_out = nskip_stress
!               else if( word == '*MD cell  ' ) then readif
!
!                  call read_cell( nfile, myid, nodes, iunit, h )
!
!                  lnotread_MDcell = .false.
!
!               else if( word == '*supercell' ) then readif
!
!                  call read_cell( nfile, myid, nodes, iunit, hcell )
!
!                  lnotread_supercell = .false.

         else if( word == '*total mom' ) then readif

!            call read_totmom( nfile, myid, nodes, iunit,  &
!& ltotmom, nskip_totmom )

         else if( word == '*atoms    ' ) then readif

            call remd_read_atom( nfile, myid, nodes, iunit,  &
& ntype, natom, zatom, nhk, icscale, vrandom, lfixion, ratm, vatm,  &
& lterminator, lQMMDrun, lnotread_MDatoms, b, disp,  &
& ldispion, dvector, ficmass )

            if( .not.lQMMDrun .and. .not.lpureMD ) ntot(1:ntype) = nhk(1:ntype)

         else if( word == '*MD atoms ' ) then readif

            call remd_read_MDatom( nfile, myid, nodes, iunit,  &
& ntmd, nmd, zatom_md, icscale_md, vrandom, lfixion,  &
& x(1,0), x(1,1), is(1), nmd_alloc,  &
& sxog, syog, szog, anxi, anyi, anzi, b, disp,  &
& ldispion, dvector, ficmass, ntot )

!         else if( word == '*MD termin' ) then readif
!
!            if( lQMMDrun ) then
!               call remd_read_HSMDatom( nfile, myid, nodes, iunit, &
!& ntHSMD, nHSMD, zatom_HSMD, ntot_HSMD, icscale_HSMD,  &
!& x_HSMD, x_HSMD, b, disp )
!            end if

!         else if( word == '*symmetry ' ) then readif
!
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

!         else if( word == '*soft wall' ) then readif
!
!            call remd_read_walls( nfile, myid, nodes, iunit )

!         else if( word == '*gravitati' ) then readif
!
!            call remd_read_gravity( nfile, myid, nodes, iunit,  &
!& lgravi, gravmag, gravdir )

!         else if( word == '*electric ' ) then readif
!
!            call read_efield( nfile, myid, nodes, iunit,  &
!& lefield, efield, lsawtooth, lsawtooth_shape, loutpolarization,  &
!& lefield_start, lconstraintD )

!         else if( word == '*virtual m' ) then readif
!
!            call remd_read_vmd( nfile, myid, nodes,  &
!& iunit, b )

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


!-----if no MD atoms in the file, copy QM atoms
if( lnotread_MDatoms ) then
    call copy_QM_to_MD( nfile, myid, nodes,  &
& ntype, natom, zatom, nhk, icscale, ratm, vatm,  &
& ntmd, nmd, zatom_md, icscale_md, x, is, nmd_alloc,  &
& sxog, syog, szog, anxi, anyi, anzi, b, disp, r1max, r1min )
    !---Einstein solid in virtual MD
!    call shift_lattice_vmd( nfile, myid, nodes )!, disp )
end if

!-----zero clear the number of QM atoms if pure classical MD
!if( lpureMD ) then
!    ntype = 0
!    natom = 0
!end if


!if( lQMMDrun ) then
!    !-----get the number of QM terminator (HSQM) atoms
!    do i = 1, ntype
!       if( lterminator(i) ) then
!           ntHSQM = ntHSQM + 1
!            nHSQM =  nHSQM + nhk(i)
!       end if
!    end do
!
!    !-----get the number of MD cluster atoms
!    ntmd_cluster = ntype - ntHSQM + ntHSMD
!     nmd_cluster = natom -  nHSQM +  nHSMD
!end if


!-----allocate memory for MD-cluster atoms
call remd_param_atom2_alloc( nfile, myid, nodes,  &
& alloc_mem, ntmd_cluster, nmd_cluster )


!if( lQMMDrun ) then
!    !-----set variables for MD-cluster atoms
!    it = 0
!    na = 0
!    na1 = 1
!    do i = 1, ntype
!       if( .not.lterminator(i) ) then
!           it = it + 1
!           zatom_cluster(it)   = zatom(i)
!           icscale_cluster(it) = icscale(i)
!           lMDterminator(it)   = .false.
!           ntot_cluster(it)    = nhk(i)
!           do j = na1, na1 + nhk(i) - 1
!              na = na + 1
!              x_cluster(3*na-2,0) = ratm(1,j)
!              x_cluster(3*na-1,0) = ratm(2,j)
!              x_cluster(3*na-0,0) = ratm(3,j)
!              is_cluster(na)      = it
!           end do
!       end if
!       na1 = na1 + nhk(i)
!    end do
!
!    na1 = 1
!    do i = 1, ntHSMD
!       it = it + 1
!       zatom_cluster(it)   = zatom_HSMD(i)
!       icscale_cluster(it) = icscale_HSMD(i)
!       lMDterminator(it)   = .true.
!       ntot_cluster(it)    = ntot_HSMD(i)
!       do j = na1, na1 + ntot_HSMD(i) - 1
!          na = na + 1
!          x_cluster(3*na-2,0) = x_HSMD(3*j-2)
!          x_cluster(3*na-1,0) = x_HSMD(3*j-1)
!          x_cluster(3*na-0,0) = x_HSMD(3*j-0)
!          is_cluster(na)      = it
!       end do
!       na1 = na1 + ntot_HSMD(i)
!    end do
!
!    !-----error trap
!    if( it /= ntmd_cluster ) call fstop( nfile, myid, nodes,  &
!&       'error: it /= ntmd_cluster' )
!    if( na /= nmd_cluster ) call fstop( nfile, myid, nodes,  &
!&       'error: na /= nmd_cluster' )
!
!end if


!-----set origin of symmetry operations for MD atoms
!call setsymmmd( nfile, myid, nodes, disp )

if( ifmd == 10 ) irstrct = 0

!if( irstrct /= 0 .and. irstrct /= 10 ) then
!    !-----restriction for MD cell
!    call restrictMDcell( irstrct, irstrct_sub, h(1,1,0), lcell_rstrct )
!    if( .not.lQMMDrun ) then
!        hcell(1:3,1:3) = h(1:3,1:3,0)
!    end if
!end if


!-----set soft walls
!do i = 1, nwalls
!   q1 = b(1,1)*wallp(1,i) + b(2,1)*wallp(2,i) + b(3,1)*wallp(3,i)
!   q2 = b(1,2)*wallp(1,i) + b(2,2)*wallp(2,i) + b(3,2)*wallp(3,i)
!   q3 = b(1,3)*wallp(1,i) + b(2,3)*wallp(2,i) + b(3,3)*wallp(3,i)
!   q1 = q1 + disp(1)
!   q2 = q2 + disp(2)
!   q3 = q3 + disp(3)
!   wallp(1,i) = h(1,1,0)*q1 + h(1,2,0)*q2 + h(1,3,0)*q3
!   wallp(2,i) = h(2,1,0)*q1 + h(2,2,0)*q2 + h(2,3,0)*q3
!   wallp(3,i) = h(3,1,0)*q1 + h(3,2,0)*q2 + h(3,3,0)*q3
!
!   vv = wallv(1,i)*wallv(1,i) + wallv(2,i)*wallv(2,i)  &
!&     + wallv(3,i)*wallv(3,i)
!   vv = sqrt(vv)
!   wallv(1,i) = wallv(1,i) / vv
!   wallv(2,i) = wallv(2,i) / vv
!   wallv(3,i) = wallv(3,i) / vv
!end do


if( irstrct /= 0 .and. irstrct /= 10 ) then
    !----- transpose matrix of b = inverse of h
    CALL RCIPRL( h(1,1,0), b, volume )
end if


return
end subroutine




subroutine remd_read_atom_number( nfile, myid, nodes, iunit,  &
& ntype, natom, sxog, syog, szog, anxi, anyi, anzi, b,  &
& disp, r1max, r1min )
!-----------------------------------------------------------------------
!     get the numbers of species & atoms in each node
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
integer :: ntype
integer :: natom
real*8 :: anxi, anyi, anzi
real*8 :: sxog, syog, szog            ! Reduced node origin
real*8,  dimension(3,3) :: b
real*8, dimension(3) :: disp
real*8, dimension(3) :: r1max, r1min

!-----declare local variables
character(10) :: word, cpot
integer :: istat
integer :: it, ix
logical :: lPME
integer :: np_PME


readdo: do
   read(iunit,'(a10)',iostat=istat) word

   readif: if( istat == 0 ) then

      selectif: if( word == '(species) ' ) then

         read(iunit,*,iostat=istat) ntype

         atomif: if( istat == 0 ) then

            do ix = 1, 3
               r1max(ix) = -1.d+10
               r1min(ix) =  1.d+10
            end do

            natom = 0
            do it = 1, ntype
               call remd_read_natom( nfile, myid, nodes,  &
& iunit, natom, sxog, syog, szog, anxi, anyi, anzi, b,  &
& disp, r1max, r1min, it )
            end do

         else atomif
            !----- error trap
            call fstop( nfile, myid, nodes,  &
&                         'error-0001 in remd_read_atom_number' )
         end if atomif

      else if( word == '(interatom' ) then selectif

!         read(iunit,'(a10)',iostat=istat) cpot
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0002 in remd_read_atom_number' )
!
!         call set_potential( nfile, myid, nodes, cpot, ntype )

      else if( word == '(PME metho' ) then selectif

!         read(iunit,*,iostat=istat) lPME
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0003 in remd_read_atom_number' )
!
!         call set_lPME( nfile, myid, nodes, lPME )

      else if( word == '(interpola' ) then selectif

!         read(iunit,*,iostat=istat) np_PME
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0004 in remd_read_atom_number' )
!
!         call set_np_PME( nfile, myid, nodes, np_PME )

      else if( word == '*end      ' ) then selectif

         exit

      endif selectif

   else readif

      !----- error trap
      call fstop( nfile, myid, nodes,  &
&                 'error-0000 in remd_read_atom_number' )

   end if readif
end do readdo


return
end subroutine




subroutine remd_read_natom( nfile, myid, nodes,  &
& iunit, natom, sxog, syog, szog, anxi, anyi, anzi, b,  &
& disp, r1max, r1min, it )
!-----------------------------------------------------------------------
! read the number of atoms
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
integer :: natom
real*8 :: anxi, anyi, anzi
real*8 :: sxog, syog, szog            ! Reduced node origin
real*8,  dimension(3,3) :: b
real*8, dimension(3) :: disp
real*8, dimension(3) :: r1max, r1min
integer :: it

!-----declare local variables
character(10) :: word
character(50) :: fname
character(10) :: unit
real*8,  parameter :: audang = 0.529177249d0
integer :: istat
integer :: nhk
integer :: icscale
integer :: keyword
integer :: ierror, ntotal, i, n_i, ix, iy, iz
integer :: ncl1, ncl2, ncl3
real*8  :: x_i, y_i, z_i
real*8  :: q1, q2, q3
integer :: junit
real(8) :: zatom
character(80) :: ANNfiles
character(10) :: ANNlunit, ANNeunit
logical :: lANNascii, lANNsc
real(8) :: ANNscepsilon, ANNscsigma, ANNscpower


call allocate_unit_number( junit )

unit = '(bohr)    '

readdo: do
   read(iunit,'(a10)',iostat=istat) word

   readif: if( istat == 0 ) then

      selectif: if( word == '(the numbe' ) then

         read(iunit,*,iostat=istat) nhk
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0001 in remd_read_natom' )

         config_1: if( nhk > 0 ) then

            rereaddo1: do
               read(iunit,'(a10)',iostat=istat) word

               rereadif1: if( istat == 0 ) then

                  reselectif1: if( word == '(positions' ) then

            read(iunit,*,iostat=istat) ncl1, ncl2, ncl3
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0020 in remd_read_natom' )

            !----- error trap
            if( mod(nhk,(ncl1*ncl2*ncl3)) /= 0 )  &
&               call fstop( nfile, myid, nodes,  &
&                           'wrong # of ions in remd_read_natom' )

            read(iunit,*,iostat=istat) icscale
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0021 in remd_read_natom' )

            !----- error trap
            if( icscale /= 1 .and. ncl1*ncl2*ncl3 > 1 )  &
&               call fstop( nfile, myid, nodes,  &
&      'not supported for icscale /= 1 .and. ncl1*ncl2*ncl3 > 1' )

            !----- read configuration
            do i = 1, nhk/(ncl1*ncl2*ncl3)
               read(iunit,*,iostat=istat) x_i, y_i, z_i
               !----- error trap
               if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                               'error-0022 in remd_read_natom' )

               !----- length unit in [bohr]
               if( icscale /= 1 ) then
                   if( unit == '(ang)     ' ) then
                       x_i = x_i / audang
                       y_i = y_i / audang
                       z_i = z_i / audang
                   end if
               end if

               do iz = 1, ncl3
               do iy = 1, ncl2
               do ix = 1, ncl1
                  !----- set scaled coordinates
                  scaleif: if( icscale == 1 ) then
                     q1 = (x_i + dble(ix-1))/dble(ncl1)
                     q2 = (y_i + dble(iy-1))/dble(ncl2)
                     q3 = (z_i + dble(iz-1))/dble(ncl3)
                  else scaleif
                     q1 = b(1,1)*x_i + b(2,1)*y_i + b(3,1)*z_i
                     q2 = b(1,2)*x_i + b(2,2)*y_i + b(3,2)*z_i
                     q3 = b(1,3)*x_i + b(2,3)*y_i + b(3,3)*z_i
                  end if scaleif
                  q1 = q1 + disp(1)
                  q2 = q2 + disp(2)
                  q3 = q3 + disp(3)
                  r1max(1) = max( r1max(1), q1 )
                  r1min(1) = min( r1min(1), q1 )
                  r1max(2) = max( r1max(2), q2 )
                  r1min(2) = min( r1min(2), q2 )
                  r1max(3) = max( r1max(3), q3 )
                  r1min(3) = min( r1min(3), q3 )
                  if( q1 >= sxog .and. q1 < sxog + anxi .and.  &
&                     q2 >= syog .and. q2 < syog + anyi .and.  &
&                     q3 >= szog .and. q3 < szog + anzi  ) then
                      natom = natom + 1
                  end if
               end do
               end do
               end do

            end do

                     exit

                  else if( word == '(unit of l' ) then reselectif1

                      read(iunit,'(a10)',iostat=istat) unit
                      if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                                     'error-0030 in remd_read_natom' )

                  endif reselectif1

               else rereadif1

                  !----- error trap
                  call fstop( nfile, myid, nodes,  &
&                             'error-0100 in remd_read_natom' )

               end if rereadif1
            end do rereaddo1

         else config_1

            rereaddo: do
               read(iunit,'(a10)',iostat=istat) word

               rereadif: if( istat == 0 ) then

                  reselectif: if( word == '(position ' ) then

         !----- read configuration from a file
            read(iunit,*,iostat=istat) fname
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0010 in remd_read_natom' )

            read(iunit,*,iostat=istat) icscale
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0011 in remd_read_natom' )

            !----- open configuration file
            if( myid == 0 ) then
                write(nfile(1),*) 'open file: ',  &
&                                 fname(1:len_trim(fname))
                write(nfile(2),*) 'open file: ',  &
&                                 fname(1:len_trim(fname))
            endif
            open( junit, file=fname, status='old', action='read',  &
&                 iostat=ierror )
            blockif_1: if( ierror == 0 ) then

         readdo_config: do

               read(iunit,*,iostat=istat) keyword
               !----- error trap
               if( istat < 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0012 in remd_read_natom' )

               if( istat > 0 ) exit

               read(junit,*,iostat=istat) ntotal
               !----- error trap
               if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                               'error-0013 in remd_read_natom' )

               !----- read configuration
               do i = 1, ntotal
                  read(junit,*,iostat=istat) n_i, x_i, y_i, z_i
                  !----- error trap
                  if( istat /= 0 ) call fstop( nfile, myid, nodes, &
&                                'error-0014 in remd_read_natom' )
                  if( n_i == keyword ) then
                      !----- set scaled coordinates
                      scaleif2: if( icscale == 1 ) then
                         q1 = x_i
                         q2 = y_i
                         q3 = z_i
                      else scaleif2

                         !----- length unit in [bohr]
                         if( unit == '(ang)     ' ) then
                             x_i = x_i / audang
                             y_i = y_i / audang
                             z_i = z_i / audang
                         end if

                         q1 = b(1,1)*x_i + b(2,1)*y_i + b(3,1)*z_i
                         q2 = b(1,2)*x_i + b(2,2)*y_i + b(3,2)*z_i
                         q3 = b(1,3)*x_i + b(2,3)*y_i + b(3,3)*z_i
                      end if scaleif2
                      q1 = q1 + disp(1)
                      q2 = q2 + disp(2)
                      q3 = q3 + disp(3)
                      r1max(1) = max( r1max(1), q1 )
                      r1min(1) = min( r1min(1), q1 )
                      r1max(2) = max( r1max(2), q2 )
                      r1min(2) = min( r1min(2), q2 )
                      r1max(3) = max( r1max(3), q3 )
                      r1min(3) = min( r1min(3), q3 )
                      if( q1 >= sxog .and. q1 < sxog + anxi .and.  &
&                         q2 >= syog .and. q2 < syog + anyi .and.  &
&                         q3 >= szog .and. q3 < szog + anzi ) then
                          natom = natom + 1
                      end if
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

                  else if( word == '(unit of l' ) then reselectif

                      read(iunit,'(a10)',iostat=istat) unit
                      if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                                     'error-0015 in remd_read_natom' )

                  endif reselectif

               else rereadif

                  !----- error trap
                  call fstop( nfile, myid, nodes,  &
&                             'error-1001 in remd_read_natom' )

               end if rereadif
            end do rereaddo

         end if config_1

!         exit

      else if( word == '(unit of l' ) then selectif

         read(iunit,'(a10)',iostat=istat) unit
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                                     'error-2000 in remd_read_natom' )

      else if( word == '(atomic nu' ) then selectif

         read(iunit,*,iostat=istat) zatom
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                                     'error-2001 in remd_read_natom' )

         !----- set atom type name
!         call potpreset_ANNatom( nfile, myid, nodes, it, zatom )

      else if( word == '(ANN poten' ) then selectif

!         read(iunit,*,iostat=istat) ANNfiles
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                                     'error-2002 in remd_read_natom' )
!
!         read(iunit,'(a10)',iostat=istat) ANNlunit
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                                     'error-2003 in remd_read_natom' )
!
!         read(iunit,'(a10)',iostat=istat) ANNeunit
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                                     'error-2004 in remd_read_natom' )
!
!         read(iunit,*,iostat=istat) lANNascii
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                                     'error-2005 in remd_read_natom' )
!
!         !----- set ANN file name
!         call potpreset_ANNfiles( nfile, myid, nodes,  &
!& it, ANNfiles, ANNlunit, ANNeunit, lANNascii )

      else if( word == '(ANN soft ' ) then selectif

!         read(iunit,*,iostat=istat) lANNsc
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                                     'error-2006 in remd_read_natom' )
!
!         read(iunit,*,iostat=istat) ANNscepsilon
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                                     'error-2007 in remd_read_natom' )
!
!         read(iunit,*,iostat=istat) ANNscsigma
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                                     'error-2008 in remd_read_natom' )
!
!         read(iunit,*,iostat=istat) ANNscpower
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                                     'error-2009 in remd_read_natom' )
!
!         !----- set ANN soft core
!         call potpreset_ANNsc( nfile, myid, nodes,  &
!& it, lANNsc, ANNscepsilon, ANNscsigma, ANNscpower )

      else if( word == '(end)     ' ) then selectif

         exit

      endif selectif

   else readif

      !----- error trap
      call fstop( nfile, myid, nodes,  &
&                 'error-0000 in remd_read_natom' )

   end if readif
end do readdo

call deallocate_unit_number( junit )


return
end subroutine




subroutine remd_read_atom( nfile, myid, nodes, iunit,  &
& ntype, natom, zatom, nhk, icscale, vrandom, lfixion, ratm, vatm,  &
& lterminator, lQMMDrun, lnotread_MDatoms, b, disp,  &
& ldispion, dvector, ficmass )
!-----------------------------------------------------------------------
! read variables for atoms
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
integer :: ntype
integer :: natom
real*8,  dimension(ntype) :: zatom
integer, dimension(ntype) :: nhk
integer, dimension(ntype) :: icscale
logical, dimension(ntype) :: vrandom
logical, dimension(ntype) :: lfixion
real*8,  dimension(3,natom) :: ratm
real*8,  dimension(3,natom) :: vatm
logical, dimension(ntype) :: lterminator
logical :: lQMMDrun
logical :: lnotread_MDatoms
real*8,  dimension(3,3) :: b
real*8, dimension(3) :: disp
logical, dimension(ntype) :: ldispion
real*8,  dimension(3,ntype) :: dvector
real*8,  dimension(ntype) :: ficmass

!-----declare local variables
character(10) :: word
integer :: istat, it, nhk1, i


readdo: do
   read(iunit,'(a10)',iostat=istat) word

   readif: if( istat == 0 ) then

      selectif: if( word == '(species) ' ) then

         read(iunit,*,iostat=istat) ntype

         atomif: if( istat == 0 ) then

            nhk1 = 1
            do it = 1, ntype
               call remd_read_type( nfile, myid, nodes, iunit,  &
& zatom(it), nhk(it), icscale(it), vrandom(it), lfixion(it),  &
& ratm(1,nhk1), vatm(1,nhk1), lterminator(it),  &
& lQMMDrun, lnotread_MDatoms, .true., natom-nhk1+1,  &
& 0.d0, 0.d0, 0.d0, 1.d0, 1.d0, 1.d0, b, disp,  &
& ldispion(it), dvector(1,it), ficmass(it), it )
               nhk1 = nhk1 + nhk(it)
            end do

         else atomif
            !----- error trap
            call fstop( nfile, myid, nodes,  &
&                         'error-0001 in remd_read_atom' )
         end if atomif

      else if( word == '*end      ' ) then selectif

         exit

      endif selectif

   else readif

      !----- error trap
      call fstop( nfile, myid, nodes,  &
&                 'error-0000 in remd_read_atom' )

   end if readif
end do readdo


return
end subroutine




subroutine remd_read_type( nfile, myid, nodes, iunit,  &
& zatom, nhk, icscale, vrandom, lfixion, ratm, vatm,  &
& lterminator, lQMMDrun, lnotread_MDatoms, lreadratm, mxn,  &
& sxog, syog, szog, anxi, anyi, anzi, b, disp,  &
& ldispion, dvector, ficmass, it )
!-----------------------------------------------------------------------
! read variables for atomic data
!-----------------------------------------------------------------------
use constants
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
real*8  :: zatom
integer :: nhk
integer :: icscale
logical :: vrandom
logical :: lfixion
real*8  :: ratm(3,*), vatm(3,*)
logical :: lterminator
logical :: lQMMDrun
logical :: lnotread_MDatoms
logical :: lreadratm
integer :: mxn
real*8  :: anxi, anyi, anzi
real*8  :: sxog, syog, szog            ! Reduced node origin
real*8  :: b(3,3)
real*8  :: disp(3)
logical :: ldispion
real*8  :: dvector(3)
real*8  :: ficmass
integer :: it

!-----declare local variables
character(10) :: word
character(50) :: fname_pos, fname_vel
character(10) :: unit_length, unit_energy
!real*8,  parameter :: audang = 0.529177249d0
integer :: istat, keyword, ierror, ntotal, ntotal2, i
integer :: n_i, nion, nionv, ncl1, ncl2, ncl3, ix, iy, iz, nv_i
real*8  :: x_i, y_i, z_i, q1, q2, q3, vx_i, vy_i, vz_i
real*8  :: vmax
real*8  :: scepsilon, scsigma, scpower
integer :: junit, junit2
integer :: nkeyword_pos, keyword_pos(100), nkeyword_vel, keyword_vel(100)
character(10) :: unit_length_pos
integer :: natm
logical, allocatable :: laccept(:)


nion  = 0
nionv = 0
unit_length = '(bohr)    '
unit_energy = '(hr)      '
fname_pos = ''
fname_vel = ''
nkeyword_pos = 0
nkeyword_vel = 0

readdo: do
   read(iunit,'(a10)',iostat=istat) word

   readif: if( istat == 0 ) then

      selectif: if( word == '(atomic nu' ) then

         read(iunit,*,iostat=istat) zatom
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0001 in remd_read_type' )

      else if( word == '(the numbe' ) then selectif

         read(iunit,*,iostat=istat) nhk
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0007 in remd_read_type' )

!               !----- error trap
!               if( nhk > mxn ) call fstop( nfile, myid, nodes,  &
!&                                          'nhk > mxn' )

      else if( word == '(position ' ) then selectif

         config_1: if( nhk == 0 ) then

            !----- read configuration from a file
            read(iunit,*,iostat=istat) fname_pos
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0010 in remd_read_type' )

            read(iunit,*,iostat=istat) icscale
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0011 in remd_read_type' )

            nkeyword_pos = 0
            do
               read(iunit,*,iostat=istat) keyword
               !----- error trap
               if( istat < 0 ) call fstop( nfile, myid, nodes,  &
&                               'error-0012 in remd_read_type' )

               if( istat > 0 ) exit
               nkeyword_pos = nkeyword_pos + 1
               keyword_pos(nkeyword_pos) = keyword
            end do

            !---save the current unit of length
            unit_length_pos = unit_length

         end if config_1

      else if( word == '(positions' ) then selectif

         config_2: if( nhk > 0 ) then

            read(iunit,*,iostat=istat) ncl1, ncl2, ncl3
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0020 in remd_read_type' )

            !----- error trap
            if( mod(nhk,(ncl1*ncl2*ncl3)) /= 0 )  &
&               call fstop( nfile, myid, nodes,  &
&                           'wrong # of ions in remd_read_type' )

            read(iunit,*,iostat=istat) icscale
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0021 in remd_read_type' )

            !----- error trap
            if( icscale /= 1 .and. ncl1*ncl2*ncl3 > 1 )  &
&               call fstop( nfile, myid, nodes,  &
&      'not supported for icscale /= 1 .and. ncl1*ncl2*ncl3 > 1' )

            !---allocate memory
            allocate( laccept(nhk) )
            laccept(1:nhk) = .false.

            !----- read configuration
            nion = 0
            natm = 0
            do i = 1, nhk/(ncl1*ncl2*ncl3)
               read(iunit,*,iostat=istat) x_i, y_i, z_i
               !----- error trap
               if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                               'error-0022 in remd_read_type' )

               !----- length unit in [bohr]
               if( icscale /= 1 ) then
                   if( unit_length == '(ang)     ' ) then
                       x_i = x_i / audang
                       y_i = y_i / audang
                       z_i = z_i / audang
                   end if
               end if

               do iz = 1, ncl3
               do iy = 1, ncl2
               do ix = 1, ncl1
                  natm = natm + 1
                  !----- set scaled coordinates
                  scaleif: if( icscale == 1 ) then
                     q1 = (x_i + dble(ix-1))/dble(ncl1)
                     q2 = (y_i + dble(iy-1))/dble(ncl2)
                     q3 = (z_i + dble(iz-1))/dble(ncl3)
                  else scaleif
                     q1 = b(1,1)*x_i + b(2,1)*y_i + b(3,1)*z_i
                     q2 = b(1,2)*x_i + b(2,2)*y_i + b(3,2)*z_i
                     q3 = b(1,3)*x_i + b(2,3)*y_i + b(3,3)*z_i
                  end if scaleif
                  q1 = q1 + disp(1)
                  q2 = q2 + disp(2)
                  q3 = q3 + disp(3)
                  if( lreadratm .or.  &
&                   ( q1 >= sxog .and. q1 < sxog + anxi .and.  &
&                     q2 >= syog .and. q2 < syog + anyi .and.  &
&                     q3 >= szog .and. q3 < szog + anzi ) ) then
                      nion = nion + 1
                      !----- error trap
                      if( nion > mxn ) call fstop( nfile, myid,  &
&                                      nodes, 'nion > mxn' )
                      laccept(natm) = .true.
                      ratm(1,nion) = q1 - sxog
                      ratm(2,nion) = q2 - syog
                      ratm(3,nion) = q3 - szog
                  end if
               end do
               end do
               end do

            end do
            icscale = 1

         end if config_2

      else if( word == '(velocity ' ) then selectif

         velocity_1: if( nhk == 0 .and. lnotread_MDatoms ) then

         !----- read velocities from a file
            read(iunit,*,iostat=istat) fname_vel
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0030 in remd_read_type' )

            nkeyword_vel = 0
            do
               read(iunit,*,iostat=istat) keyword
               !----- error trap
               if( istat < 0 ) call fstop( nfile, myid, nodes,  &
&                              'error-0031 in remd_read_type' )

               if( istat > 0 ) exit
               nkeyword_vel = nkeyword_vel + 1
               keyword_vel(nkeyword_vel) = keyword
            end do

         end if velocity_1

      else if( word == '(velocitie' ) then selectif

         velocity_2: if( nhk > 0 .and. lnotread_MDatoms ) then

            vrandom = .false.

            read(iunit,*,iostat=istat) ncl1, ncl2, ncl3
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0040 in remd_read_type' )

            !----- error trap
            if( mod(nhk,(ncl1*ncl2*ncl3)) /= 0 )  &
&               call fstop( nfile, myid, nodes,  &
&                        'wrong # of ions (2) in remd_read_type' )

            read(iunit,*,iostat=istat) vmax
            !----- error trap
            if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                            'error-0041 in remd_read_type' )

            !----- error trap
            if( .not.allocated(laccept) ) call fstop( nfile, myid, nodes,  &
&                                         'laccept is not allocated in remd_read_type' )

            !----- read velocities
            nionv = 0
            natm = 0
            do i = 1, nhk/(ncl1*ncl2*ncl3)
               read(iunit,*,iostat=istat) x_i, y_i, z_i
               !----- error trap
               if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                               'error-0042 in remd_read_type' )

               do iz = 1, ncl3
               do iy = 1, ncl2
               do ix = 1, ncl1
                  natm = natm + 1
                  if( laccept(natm) ) then
                      nionv = nionv + 1
                      vatm(1,nionv) = x_i * vmax
                      vatm(2,nionv) = y_i * vmax
                      vatm(3,nionv) = z_i * vmax
                  end if
               end do
               end do
               end do

            end do

         end if velocity_2

      else if( word == '(fix posit' ) then selectif

         if( lnotread_MDatoms ) then
             read(iunit,*,iostat=istat) lfixion
             !----- error trap
             if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                             'error-0043 in remd_read_type' )
         end if

!      else if( word == '(fix displ' ) then selectif
!
!         if( lnotread_MDatoms ) then
!             read(iunit,*,iostat=istat) ldispion
!             !----- error trap
!             if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                             'error-0045 in remd_read_type' )
!             read(iunit,*,iostat=istat) (dvector(ix),ix=1,3)
!             !----- error trap
!             if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                             'error-0046 in remd_read_type' )
!         end if

      else if( word == '(fictitiou' ) then selectif

         if( lnotread_MDatoms ) then
             read(iunit,*,iostat=istat) ficmass
             !----- error trap
             if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                             'error-0047 in remd_read_type' )
         end if

!      else if( word == '(terminato' ) then selectif
!
!         if( lQMMDrun ) then
!             read(iunit,*,iostat=istat) lterminator
!             !----- error trap
!             if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                             'error-0044 in remd_read_type' )
!         end if

      else if( word == '(unit of l' ) then selectif

         read(iunit,'(a10)',iostat=istat) unit_length
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                             'error-0045 in remd_read_type' )

      else if( word == '(unit of e' ) then selectif

         read(iunit,'(a10)',iostat=istat) unit_energy
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                             'error-0046 in remd_read_type' )

      else if( word == '(soft core' ) then selectif

!         read(iunit,*,iostat=istat) scepsilon
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                             'error-0047 in remd_read_type' )
!
!         read(iunit,*,iostat=istat) scsigma
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                             'error-0048 in remd_read_type' )
!
!         read(iunit,*,iostat=istat) scpower
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                             'error-0049 in remd_read_type' )
!
!         !----- energy unit in [hr.]
!         if( unit_energy == '(ry)      ' ) then
!             scepsilon = scepsilon / 2.d0
!         else if( unit_energy == '(ev)      ' ) then
!             scepsilon = scepsilon / hrdev
!         end if
!
!         !----- length unit in [bohr]
!         if( unit_length == '(ang)     ' ) then
!             scsigma = scsigma / audang
!         end if
!
!         !----- set potential parameters
!         call potpreset_soft_core( nfile, myid, nodes,  &
!& it, scepsilon, scsigma, scpower )


      else if( word == '(end)     ' ) then selectif

         exit

      endif selectif

   else readif

      !----- error trap
      call fstop( nfile, myid, nodes,  &
&                 'error-0000 in remd_read_type' )

   end if readif
end do readdo


!---read configurations and/or velocities from files
call allocate_unit_number( junit )
call allocate_unit_number( junit2 )

config_3: if( nhk == 0 ) then

    !----- error trap
    if( fname_pos == '' ) call fstop( nfile, myid, nodes,  &
&                                     'error-0100 in remd_read_type' )

    vrandom = fname_vel == ''

    !----- open configuration file
    if( myid == 0 ) then
        write(nfile(1),*) 'open file: ', trim(fname_pos)
        write(nfile(2),*) 'open file: ', trim(fname_pos)
    end if
    open( junit, file=trim(fname_pos), status='old', action='read', iostat=ierror )

    !----- error trap
    if( ierror /= 0 ) call fstop( nfile, myid, nodes,  &
&                                 'cannot open file :'//trim(fname_pos) )

    if( .not.vrandom ) then
        !----- open velocity file
        if( myid == 0 ) then
            write(nfile(1),*) 'open file: ', trim(fname_vel)
            write(nfile(2),*) 'open file: ', trim(fname_vel)
        end if
        open( junit2, file=trim(fname_vel), status='old', action='read', iostat=ierror )

        if( ierror /= 0 ) then
            if( myid == 0 ) then
                write(nfile(1),*) 'cannot open file: ', trim(fname_vel), ', set vrandom = .true.'
                write(nfile(2),*) 'cannot open file: ', trim(fname_vel), ', set vrandom = .true.'
            end if
            ierror = 0
            vrandom = .true.
        end if
        !----- error trap
        if( .not.vrandom .and. nkeyword_vel /= nkeyword_pos ) then
            if( myid == 0 ) then
                write(nfile(1),*) 'error in file(1): ', trim(fname_vel), ', set vrandom = .true.'
                write(nfile(2),*) 'error in file(1): ', trim(fname_vel), ', set vrandom = .true.'
            end if
            vrandom = .true.
        end if
    end if

    nion = 0
    nionv = 0
    readdo_config2: do keyword = 1, nkeyword_pos

       read(junit,*,iostat=istat) ntotal
       !----- error trap
       if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                                   'error-0101 in remd_read_type' )

       if( .not.vrandom ) then
           read(junit2,*,iostat=istat) ntotal2
           !----- error trap
           if( istat /= 0 ) then
               if( myid == 0 ) then
                   write(nfile(1),*) 'error in file(2): ', trim(fname_vel), ', set vrandom = .true.'
                   write(nfile(2),*) 'error in file(2): ', trim(fname_vel), ', set vrandom = .true.'
               end if
               istat = 0
               vrandom = .true.
           end if
           !----- error trap
           if( ntotal2 /= ntotal ) then
               if( myid == 0 ) then
                   write(nfile(1),*) 'error in file(3): ', trim(fname_vel), ', set vrandom = .true.'
                   write(nfile(2),*) 'error in file(3): ', trim(fname_vel), ', set vrandom = .true.'
               end if
               vrandom = .true.
           end if

           read(junit2,*,iostat=istat) vmax
           !----- error trap
           if( istat /= 0 ) then
               if( myid == 0 ) then
                   write(nfile(1),*) 'error in file(4): ', trim(fname_vel), ', set vrandom = .true.'
                   write(nfile(2),*) 'error in file(4): ', trim(fname_vel), ', set vrandom = .true.'
               end if
               istat = 0
               vrandom = .true.
           end if
       end if

       !----- read configuration
       do i = 1, ntotal
          read(junit,*,iostat=istat) n_i, x_i, y_i, z_i
          !----- error trap
          if( istat /= 0 ) call fstop( nfile, myid, nodes, &
&                                      'error-0102 in remd_read_type' )

          if( .not.vrandom ) then
              read(junit2,*,iostat=istat) nv_i, vx_i, vy_i, vz_i
              !----- error trap
              if( istat /= 0 ) then
                  if( myid == 0 ) then
                      write(nfile(1),*) 'error in file(5): ', trim(fname_vel), ', set vrandom = .true.'
                      write(nfile(2),*) 'error in file(5): ', trim(fname_vel), ', set vrandom = .true.'
                  end if
                  istat = 0
                  vrandom = .true.
              end if
              !----- error trap
              if( n_i /= nv_i ) then
                  if( myid == 0 ) then
                      write(nfile(1),*) 'error in file(6): ', trim(fname_vel), ', set vrandom = .true.'
                      write(nfile(2),*) 'error in file(6): ', trim(fname_vel), ', set vrandom = .true.'
                  end if
                  vrandom = .true.
              end if
          end if

          if( n_i == keyword_pos(keyword) ) then
              !----- set scaled coordinates
              scaleif2: if( icscale == 1 ) then
                  q1 = x_i
                  q2 = y_i
                  q3 = z_i
              else scaleif2
                  !----- length unit in [bohr]
                  if( unit_length_pos == '(ang)     ' ) then
                      x_i = x_i / audang
                      y_i = y_i / audang
                      z_i = z_i / audang
                  end if

                  q1 = b(1,1)*x_i + b(2,1)*y_i + b(3,1)*z_i
                  q2 = b(1,2)*x_i + b(2,2)*y_i + b(3,2)*z_i
                  q3 = b(1,3)*x_i + b(2,3)*y_i + b(3,3)*z_i
              end if scaleif2
              q1 = q1 + disp(1)
              q2 = q2 + disp(2)
              q3 = q3 + disp(3)
              if( lreadratm .or.  &
&                 (q1 >= sxog .and. q1 < sxog + anxi .and.  &
&                  q2 >= syog .and. q2 < syog + anyi .and.  &
&                  q3 >= szog .and. q3 < szog + anzi)) then
                  nion = nion + 1
                  !----- error trap
                  if( nion > mxn ) call fstop( nfile,myid, &
&                                             nodes, 'nion > mxn' )
                  ratm(1,nion) = q1 - sxog
                  ratm(2,nion) = q2 - syog
                  ratm(3,nion) = q3 - szog
                  if( .not.vrandom ) then
                      nionv = nionv + 1
                      !----- error trap
                      if( nionv > mxn ) call fstop( nfile, myid,  &
&                                       nodes, 'nionv > mxn' )
                      vatm(1,nionv) = vx_i * vmax
                      vatm(2,nionv) = vy_i * vmax
                      vatm(3,nionv) = vz_i * vmax
                  end if
              end if
          end if
       end do

       rewind(junit)
       if( .not.vrandom ) rewind(junit2)

     end do readdo_config2

     icscale = 1
     if( vrandom ) nionv = 0
     close(junit)
     close(junit2)

end if config_3

call deallocate_unit_number( junit )
call deallocate_unit_number( junit2 )


!----- error trap
!      if( nion == 0 ) call fstop( nfile, myid, nodes, 'nion == 0' )
if( nionv /= 0 .and. nionv /= nion ) then
    if( myid == 0 ) then
        write(nfile(1),*) 'nion, nionv =', nion, nionv
        write(nfile(2),*) 'nion, nionv =', nion, nionv
    endif
    call fstop( nfile, myid, nodes, 'nionv /= nion' )
endif

!----- set the number of ions in the supercell
nhk = nion

!---deallocate memory
if( allocated(laccept) ) deallocate( laccept )


return
end subroutine




subroutine copy_QM_to_MD( nfile, myid, nodes,  &
& ntype, natom, zatom, nhk, icscale, ratm, vatm,  &
& ntmd, nmd, zatom_md, icscale_md, x, is, nmd_alloc,  &
& sxog, syog, szog, anxi, anyi, anzi, b, disp, r1max, r1min )
!-----------------------------------------------------------------------
!     if this is not the QM/MD run, copy QM atoms
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: ntype
integer :: natom
real*8,  dimension(ntype) :: zatom
integer, dimension(ntype) :: nhk
integer, dimension(ntype) :: icscale
real*8,  dimension(3,natom) :: ratm
real*8,  dimension(3,natom) :: vatm

integer :: ntmd
integer :: nmd
integer :: nmd_alloc
real*8,  dimension(ntmd) :: zatom_md
integer, dimension(ntmd) :: icscale_md
real*8,  dimension(3*nmd_alloc,0:1) :: x
integer, dimension(nmd_alloc)  :: is
real*8 :: anxi, anyi, anzi
real*8 :: sxog, syog, szog            ! Reduced node origin
real*8,  dimension(3,3) :: b
real*8, dimension(3) :: disp
real*8, dimension(3) :: r1max, r1min

!---declare local variables
integer :: it, nhk1, i, ix
real*8  :: x_i, y_i, z_i
real*8  :: q1, q2, q3


!-----get displacement vector
do ix = 1, 3
   r1max(ix) = -1.d+10
   r1min(ix) =  1.d+10

   do i = 1, natom
      r1max(ix) = max( r1max(ix), ratm(ix,i) )
      r1min(ix) = min( r1min(ix), ratm(ix,i) )
   end do
   !--- error trap
   if( r1max(ix) - r1min(ix) > 1.d0 )  &
&      call fstop( nfile, myid, nodes,  &
&                 'error in scaled coordinates in copy_QM_to_MD' )

   disp(ix) = 0.5d0 *(1.d0 - (r1max(ix) + r1min(ix)) )
end do


!-----set zero for the number of MD atoms per node
nmd = 0

!-----set coordinates for MD atoms
nhk1 = 1
do it = 1, ntype
   zatom_md(it) = zatom(it)
   do i = nhk1, nhk1 + nhk(it) - 1
      x_i = ratm(1,i)
      y_i = ratm(2,i)
      z_i = ratm(3,i)
      !----- set scaled coordinates
      scaleif: if( icscale(it) == 1 ) then
          q1 = x_i
          q2 = y_i
          q3 = z_i
      else scaleif
          q1 = b(1,1)*x_i + b(2,1)*y_i + b(3,1)*z_i
          q2 = b(1,2)*x_i + b(2,2)*y_i + b(3,2)*z_i
          q3 = b(1,3)*x_i + b(2,3)*y_i + b(3,3)*z_i
      end if scaleif
      q1 = q1 + disp(1)
      q2 = q2 + disp(2)
      q3 = q3 + disp(3)
      ratm(1,i) = q1
      ratm(2,i) = q2
      ratm(3,i) = q3

      if( q1 >= sxog .and. q1 < sxog + anxi .and.  &
&         q2 >= syog .and. q2 < syog + anyi .and.  &
&         q3 >= szog .and. q3 < szog + anzi  ) then
          nmd = nmd + 1

          x(nmd*3-2,0) = ratm(1,i) - sxog
          x(nmd*3-1,0) = ratm(2,i) - syog
          x(nmd*3-0,0) = ratm(3,i) - szog
          x(nmd*3-2,1) = vatm(1,i)
          x(nmd*3-1,1) = vatm(2,i)
          x(nmd*3-0,1) = vatm(3,i)
          is(nmd) = it
      end if

   end do
   nhk1 = nhk1 + nhk(it)
   icscale_md(it) = 1
   icscale(it)    = 1
end do


return
end subroutine




subroutine remd_read_MDatom( nfile, myid, nodes, iunit,  &
& ntype, nmd, zatom, icscale, vrandom, lfixion,  &
& ratm, vatm, is, natom,  &
& sxog, syog, szog, anxi, anyi, anzi, b, disp,  &
& ldispion, dvector, ficmass, ntot )
!-----------------------------------------------------------------------
! read variables for MD atoms
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
integer :: ntype
integer :: nmd
integer :: natom
real*8,  dimension(ntype) :: zatom
integer, dimension(ntype) :: icscale
logical, dimension(ntype) :: vrandom
logical, dimension(ntype) :: lfixion
real*8,  dimension(3,natom) :: ratm
real*8,  dimension(3,natom) :: vatm
integer, dimension(natom) :: is
real*8 :: anxi, anyi, anzi
real*8 :: sxog, syog, szog            ! Reduced node origin
real*8,  dimension(3,3) :: b
real*8,  dimension(3) :: disp
logical, dimension(ntype) :: ldispion
real*8,  dimension(3,ntype) :: dvector
real*8,  dimension(ntype) :: ficmass
integer :: ntot(0:ntype)

!-----declare local variables
character(10) :: word
integer :: istat, it, nhk, i
logical :: ldummy


readdo: do
   read(iunit,'(a10)',iostat=istat) word

   readif: if( istat == 0 ) then

      selectif: if( word == '(species) ' ) then

         read(iunit,*,iostat=istat) ntype

         atomif: if( istat == 0 ) then

            nmd = 1
            do it = 1, ntype

               call remd_read_type( nfile, myid, nodes, iunit,  &
& zatom(it), nhk, icscale(it), vrandom(it), lfixion(it),  &
& ratm(1,nmd), vatm(1,nmd), ldummy, .false., .true., .false.,  &
& natom-nmd+1, sxog, syog, szog, anxi, anyi, anzi, b, disp,  &
& ldispion(it), dvector(1,it), ficmass(it), it )

               do i = nmd, nmd + nhk - 1
                  is(i) = it
               end do

               nmd = nmd + nhk
               ntot(it) = nhk

            end do
            nmd = nmd - 1

         else atomif
            !----- error trap
            call fstop( nfile, myid, nodes,  &
&                         'error-0001 in remd_read_MDatom' )
         end if atomif

      else if( word == '*end      ' ) then selectif

         exit

      endif selectif

   else readif

      !----- error trap
      call fstop( nfile, myid, nodes,  &
&                 'error-0000 in remd_read_MDatom' )

   end if readif
end do readdo


return
end subroutine




subroutine remd_read_md( nfile, myid, nodes, iunit,  &
& ifmd, dtmd, nstop, lmdstop, nstep_ini,  &
& ioptmze, lpvv, hmadisp, nhmaord, lcentdif, ibfgsclear, qnstabi, gammamin, &
& ioptmze_cell, dtcellcg, iclearcellh,  &
& lhybridopt, nstep_hybrid, nstep_hybrid_cell,  &
& treq, liscale, iscnum, iscstp, lmomzero, ichest, ihest,  &
& nnos, nresn, nyosh, tomega, hpext, tbomega, blkmod,  &
& nnos_atom, nresn_atom, nyosh_atom, tomega_atom,  &
& irstrct, irstrct_sub, lcell_rstrct,  &
& ioskip, locoor, lovelo, loforc, ioskipcoor, ioskipvelo, ioskipforc,  &
& ioskip_qm, locoor_qm, lovelo_qm, loforc_qm,  &
& ioskip_cl, locoor_cl, lovelo_cl, loforc_cl, lstat,  &
& latomic, nskip_atomic, lheat, nskip_heat,  &
& xio_min, xio_max, yio_min, yio_max, zio_min, zio_max,  &
& tol_energy, tol_force, shockspeed, nshockv,  &
& lmsstscale, msstscnum, msstscstp, vxmsst )
!-----------------------------------------------------------------------
! read variables for molecular dynamics for MD nodes
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
integer :: iunit
integer :: ifmd
real*8  :: dtmd
integer :: nstop
logical :: lmdstop
integer :: nstep_ini
integer :: ioptmze
logical :: lpvv
real*8  :: hmadisp
integer :: nhmaord
logical :: lcentdif
integer :: ibfgsclear
real*8  :: qnstabi, gammamin
integer :: ioptmze_cell
real*8  :: dtcellcg
integer :: iclearcellh
logical :: lhybridopt
integer :: nstep_hybrid, nstep_hybrid_cell
real*8  :: treq
logical :: liscale, lmomzero
integer :: iscnum, iscstp
integer :: ichest, ihest
integer :: nnos, nresn, nyosh
real*8  :: tomega
real*8  :: hpext
real*8  :: tbomega, blkmod
integer :: irstrct, irstrct_sub
logical :: lcell_rstrct(3)
integer :: nnos_atom, nresn_atom, nyosh_atom
real*8  :: tomega_atom
integer :: ioskip
logical :: locoor, lovelo, loforc
integer :: ioskipcoor, ioskipvelo, ioskipforc
integer :: ioskip_qm
logical :: locoor_qm, lovelo_qm, loforc_qm
integer :: ioskip_cl
logical :: locoor_cl, lovelo_cl, loforc_cl
logical :: lstat
logical :: latomic
integer :: nskip_atomic
logical :: lheat
integer :: nskip_heat
real*8  :: xio_min, xio_max, yio_min, yio_max, zio_min, zio_max
real*8  :: tol_energy, tol_force
real*8  :: shockspeed
integer :: nshockv(3)
logical :: lmsstscale
integer :: msstscnum
integer :: msstscstp
real*8  :: vxmsst

!-----declare local variables
character(10) :: word
integer :: i, istat
logical :: lnotread_nathermo, lnotread_nnos_atom


lnotread_nathermo  = .true.
lnotread_nnos_atom = .true.

readdo: do
   read(iunit,'(a10)',iostat=istat) word

   readif: if( istat == 0 ) then

      selectif: if( word == '(method)  ' .or. word == '(how of it' ) then

         read(iunit,*,iostat=istat) ifmd
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0001 in read_md' )

      else if( word == '(time step' ) then selectif

         read(iunit,*,iostat=istat) dtmd, nstop
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0002 in read_md' )

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

      else if( word == '(thermosta' ) then selectif

         read(iunit,*,iostat=istat) nnos
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0020 in read_md' )

         read(iunit,*,iostat=istat) nresn
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0021 in read_md' )

         read(iunit,*,iostat=istat) nyosh
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0022 in read_md' )

         read(iunit,*,iostat=istat) tomega
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0023 in read_md' )

!      else if( word == '(atom type' ) then selectif

!         lnotread_nathermo  = .false.
!
!         call read_atom_type_with_thermostat( nfile, myid, nodes, iunit )

!      else if( word == '(atom ther' ) then selectif

!         lnotread_nnos_atom = .false.
!
!         read(iunit,*,iostat=istat) nnos_atom
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0120 in read_md' )
!
!         read(iunit,*,iostat=istat) nresn_atom
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0121 in read_md' )
!
!         read(iunit,*,iostat=istat) nyosh_atom
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0122 in read_md' )
!
!         read(iunit,*,iostat=istat) tomega_atom
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0123 in read_md' )

      else if( word == '(pressure)' ) then selectif

         read(iunit,*,iostat=istat) hpext
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0024 in read_md' )

      else if( word == '(barostat ' ) then selectif

         read(iunit,*,iostat=istat) tbomega
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0025 in read_md' )

         read(iunit,*,iostat=istat) blkmod
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0026 in read_md' )

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
!&                            'error-0029 in read_md' )
!         end do

      else if( word == '(output da' ) then selectif

         read(iunit,*,iostat=istat) ioskip
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0030 in read_md' )

         read(iunit,*,iostat=istat) locoor, ioskipcoor
         if( istat /= 0 ) then
             backspace iunit
             read(iunit,*,iostat=istat) locoor
             !----- error trap
             if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                             'error-0031 in read_md' )

             ioskipcoor = ioskip
         end if

         read(iunit,*,iostat=istat) lovelo, ioskipvelo
         if( istat /= 0 ) then
             backspace iunit
             read(iunit,*,iostat=istat) lovelo
             !----- error trap
             if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                             'error-0032 in read_md' )

             ioskipvelo = ioskip
         end if

         read(iunit,*,iostat=istat) loforc, ioskipforc
         if( istat /= 0 ) then
             backspace iunit
             read(iunit,*,iostat=istat) loforc
             !----- error trap
             if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                             'error-0033 in read_md' )

             ioskipforc = ioskip
         end if

      else if( word == '(output QM' ) then selectif

         read(iunit,*,iostat=istat) ioskip_qm
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0040 in read_md' )

         read(iunit,*,iostat=istat) locoor_qm
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0041 in read_md' )

         read(iunit,*,iostat=istat) lovelo_qm
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0042 in read_md' )

         read(iunit,*,iostat=istat) loforc_qm
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0043 in read_md' )

      else if( word == '(output MD' ) then selectif

         read(iunit,*,iostat=istat) ioskip_cl
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0050 in read_md' )

         read(iunit,*,iostat=istat) locoor_cl
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0051 in read_md' )

         read(iunit,*,iostat=istat) lovelo_cl
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0052 in read_md' )

         read(iunit,*,iostat=istat) loforc_cl
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0053 in read_md' )

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

      else if( word == '(statistic' ) then selectif

         read(iunit,*,iostat=istat) lstat
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0060 in read_md' )

      else if( word == '(tolerance' ) then selectif

         read(iunit,*,iostat=istat) tol_energy
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0061 in read_md' )

         read(iunit,*,iostat=istat) tol_force
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0062 in read_md' )

      else if( word == '(optimizat' ) then selectif

         read(iunit,*,iostat=istat) ioptmze
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0063 in read_md' )

!      else if( word == '(displacem' ) then selectif

!         read(iunit,*,iostat=istat) hmadisp
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0064 in read_md' )

!      else if( word == '(order of ' ) then selectif

!         read(iunit,*,iostat=istat) nhmaord
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0065 in read_md' )

!      else if( word == '(central d' ) then selectif

!         read(iunit,*,iostat=istat) lcentdif
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0065-2 in read_md' )

!      else if( word == '(transform' ) then selectif

!         call read_phtrsfrm( nfile, myid, nodes, iunit, 0 )

!      else if( word == '(inverse t' ) then selectif

!         call read_phtrsfrm( nfile, myid, nodes, iunit, 1 )

!      else if( word == '(k points)' ) then selectif

!         call read_phkpoints( nfile, myid, nodes, iunit, .true. )

!      else if( word == '(phonon k ' ) then selectif

!         call read_phkpoints( nfile, myid, nodes, iunit, .false. )

!      else if( word == '(plot inte' ) then selectif

!         call read_phinterval( nfile, myid, nodes, iunit )

!      else if( word == '(maximum d' ) then selectif

!         call read_phrmax( nfile, myid, nodes, iunit )

!      else if( word == '(k points ' ) then selectif

!         call read_phdoskpoints( nfile, myid, nodes, iunit )

!      else if( word == '(filter fo' ) then selectif

!         call read_phdosfileter( nfile, myid, nodes, iunit )

      else if( word == '(clear Hes' ) then selectif

         read(iunit,*,iostat=istat) ibfgsclear
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0066 in read_md' )

      else if( word == '(stabilize' ) then selectif

!         read(iunit,*,iostat=istat) qnstabi
         read(iunit,*,iostat=istat) gammamin
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0066-2 in read_md' )

      else if( word == '(cell opti' ) then selectif

         read(iunit,*,iostat=istat) ioptmze_cell
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0067 in read_md' )

      else if( word == '(cell CG t' ) then selectif

         read(iunit,*,iostat=istat) dtcellcg
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0068 in read_md' )

      else if( word == '(clear cel' ) then selectif

         read(iunit,*,iostat=istat) iclearcellh
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0069 in read_md' )

      else if( word == '(hybrid op' ) then selectif

         read(iunit,*,iostat=istat) lhybridopt
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0070 in read_md' )

         read(iunit,*,iostat=istat) nstep_hybrid
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0071 in read_md' )

         read(iunit,*,iostat=istat) nstep_hybrid_cell
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0072 in read_md' )

!      else if( word == '(atomic st' ) then selectif

!         read(iunit,*,iostat=istat) latomic
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0073 in read_md' )

!         read(iunit,*,iostat=istat) nskip_atomic
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0074 in read_md' )

!      else if( word == '(output ar' ) then selectif
!
!         read(iunit,*,iostat=istat) xio_min, xio_max
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0074 in read_md' )
!
!         read(iunit,*,iostat=istat) yio_min, yio_max
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0075 in read_md' )
!
!         read(iunit,*,iostat=istat) zio_min, zio_max
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0076 in read_md' )

      else if( word == '(shock wav' ) then selectif

         read(iunit,*,iostat=istat) shockspeed
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0077 in read_md' )

         read(iunit,*,iostat=istat) nshockv(1:3)
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0078 in read_md' )

      else if( word == '(clear bar' ) then selectif

         read(iunit,*,iostat=istat) lmsstscale
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0079 in read_md' )

         read(iunit,*,iostat=istat) msstscnum
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0080 in read_md' )

         read(iunit,*,iostat=istat) msstscstp
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0081 in read_md' )

      else if( word == '(initial b' ) then selectif

         read(iunit,*,iostat=istat) vxmsst
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0082 in read_md' )

!      else if( word == '(GLE therm' ) then selectif

!         call read_GLEthermostat( nfile, myid, nodes, iunit )

!      else if( word == '(GLE frequ' ) then selectif

!         call read_GLEfrequency( nfile, myid, nodes, iunit )

!      else if( word == '(GLE param' ) then selectif

!         call read_GLEparameters( nfile, myid, nodes, iunit )

      else if( word == '(activate ' ) then selectif

         read(iunit,*,iostat=istat) lmdstop
         !----- error trap
         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
&                         'error-0083 in read_md' )

!      else if( word == '(heat flux' ) then selectif
!
!         read(iunit,*,iostat=istat) lheat
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0084 in read_md' )
!
!         read(iunit,*,iostat=istat) nskip_heat
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0085 in read_md' )

!      else if( word == '(atomic di' ) then selectif
!
!         call read_atomic_displacement( nfile, myid, nodes, iunit )

      else if( word == '*end      ' ) then selectif

         exit

      endif selectif

   else readif

      !----- error trap
      call fstop( nfile, myid, nodes, 'error-0000 in read_md' )

   end if readif
end do readdo


!----- check cell optimization
if( ifmd /= 1 ) ioptmze_cell = -1
if(    ifmd == 1 .and. ioptmze == 1  &
& .or. ifmd == 1 .and. ioptmze == 10  &
& .or. ifmd == 1 .and. ioptmze == 11  &
& .or. ifmd == 1 .and. ioptmze == 20 ) then
    if( ioptmze_cell /= -1 ) then
        if( myid == 0 ) then
            write(nfile(1),*) '*** Warning: stop cell optimization, because', &
&                             'ifmd, ioptmze=', ifmd, ioptmze
            write(nfile(2),*) '*** Warning: stop cell optimization, because', &
&                             'ifmd, ioptmze=', ifmd, ioptmze
        end if
        ioptmze_cell = -1
    end if
end if


!-----check structural optimization
lpvv = ifmd == 1 .and. ioptmze == 1
if( lpvv ) then
    ifmd = 2
    treq = 0.d0
end if


!-----check atom thermostats
!if( ifmd == 5 .or. ifmd == 6 ) then
!    if( lnotread_nathermo ) &
!&       call read_atom_type_with_thermostat2( nfile, myid, nodes )
!
!    if( lnotread_nnos_atom ) then
!        nnos_atom   = nnos
!        nresn_atom  = nresn
!        nyosh_atom  = nyosh
!        tomega_atom = tomega
!    end if
!end if


!----- check MD cell edge restriction
if( irstrct == 0 ) then
    if( lcell_rstrct(1) .or. lcell_rstrct(2) .or. lcell_rstrct(3) ) irstrct = 7
end if


return
end subroutine




subroutine set_ioptmze( ioptmze_, ioptmze_cell_ )
!-----------------------------------------------------------------------
! set default value of ioptmze
!-----------------------------------------------------------------------
use remd_param
implicit none
integer :: ioptmze_, ioptmze_cell_

ioptmze_ = ioptmze
ioptmze_cell_ = ioptmze_cell

return
end subroutine




!subroutine read_totmom( nfile, myid, nodes, iunit,  &
!& ltotmom, nskip_totmom )
!!-----------------------------------------------------------------------
!! read variables for total momentum calculation
!!-----------------------------------------------------------------------
!implicit none
!integer :: nfile(*), myid, nodes
!integer :: iunit
!logical :: ltotmom
!integer :: nskip_totmom
!character(10) :: word
!integer :: istat
!
!
!readdo: do
!   read(iunit,'(a10)',iostat=istat) word
!
!   readif: if( istat == 0 ) then
!
!      selectif: if( word == '(on/off)  ' .or. word == '(how of it' ) then
!
!         read(iunit,*,iostat=istat) ltotmom
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0001 in read_totmom' )
!
!      else if( word == '(skip step' ) then selectif
!
!         read(iunit,*,iostat=istat) nskip_totmom
!         !----- error trap
!         if( istat /= 0 ) call fstop( nfile, myid, nodes,  &
!&                         'error-0002 in read_totmom' )
!
!      else if( word == '*end      ' ) then selectif
!
!         exit
!
!      endif selectif
!
!   else readif
!
!      !----- error trap
!      call fstop( nfile, myid, nodes, 'error-0000 in read_totmom' )
!
!   end if readif
!end do readdo
!
!
!return
!end subroutine




