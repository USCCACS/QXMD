



module tddft_variables
!-----------------------------------------------------------------------
! type declaration of shared variables in tddft.f90
!-----------------------------------------------------------------------
implicit none

logical :: ltddft = .false.    ! .true. = execute MD based on TDDFT
                               ! .false.= execute MD based on conventional DFT
logical :: ltddft_fssh         ! .true. = FSSH, .false. = Ehrenfest
logical :: ltddft_start        ! .true. = restart
logical :: lfssh_switch        ! .true.  = switching aveilable
                               ! .false. = cccupations are fixed
logical :: lfssh_gsscf         ! .true.  = SCF with the ground  state
                               ! .false. = SCF with the excited state
logical :: ltddft_nscforce     ! .true. = Non-self consistent (NSC) force
logical :: lrtddft             ! .true. = calculate Casida coupling matrix

integer :: nel                 ! # of electrons

real*8  :: dttddft = 0.02d0    ! time step in  [a.u.]  Caution !! [Rydberg units]
logical :: lfssh_parallel

logical :: ltdbroad
real*8  :: tdbroad             ! width of Gaussian broadening in [Rydberg units]

logical :: lfssh_random

logical :: lfssh_boltzmn
real*8  :: treq

logical :: lfssh_vscale
real*8  :: tminimum

integer :: nbase, noccmx       !-----Number of wave packets : noccmx

integer :: nocc_change = 0
integer, allocatable, dimension(:) :: numband
real*8,  allocatable, dimension(:,:) :: occ_new

integer :: nexciton = 0
integer, allocatable, dimension(:) :: iband_hole, iband_electron
logical, allocatable, dimension(:) :: ldegenerate
integer, allocatable, dimension(:) :: iband_hole_, iband_electron_
logical, allocatable, dimension(:) :: ldegenerate_
integer, allocatable, dimension(:) :: numexciton
logical, allocatable, dimension(:) :: ldmcalc_ex

!---for DISH
logical :: lfssh_dish = .false.
integer :: ndishpair = 0
integer, allocatable, dimension(:) :: ndishi, ndishj
real*8,  allocatable, dimension(:) :: temp_decoherence_rate
real*8,  allocatable, dimension(:,:) :: decoherence_rate
integer, allocatable, dimension(:,:) :: dish_prev_step

save

end module




module tddft_fssh_variables
!-----------------------------------------------------------------------
! type declaration of shared variables in tddft.f90
!-----------------------------------------------------------------------
implicit none

logical :: ldensity_matrix_repre = .true.   
           ! .true.  = density matrix representation
           ! .false. = expansion coefficient representation

logical :: lcheck_unitary = .false.
logical :: lcheck_trimat  = .false.

logical :: lspin
logical :: lsetmatrix = .false.
integer :: noband
real*8,  allocatable, dimension(:,:) :: occ_rec
real*8,  allocatable, dimension(:,:,:) :: phisphi

real*8,  allocatable, dimension(:,:,:) :: accumprob
logical, allocatable, dimension(:,:) :: lsignc
real*8,  allocatable, dimension(:,:) :: probinteg
integer :: iudx

! for expansion coefficient representation
real*8,  allocatable, dimension(:,:,:) :: cr, ci
real*8,  allocatable, dimension(:,:) :: amtr, amti
real*8,  allocatable, dimension(:,:) :: bmthr, bmthi, bmtar, bmtai
real*8,  allocatable, dimension(:,:) :: lmtr, lmti
real*8,  allocatable, dimension(:,:) :: umtr, umti
real*8,  allocatable, dimension(:,:) :: cmtr, cmti
real*8,  allocatable, dimension(:,:) :: wrkr, wrki
integer, allocatable, dimension(:) :: nbcnt, nbdsp

! for density matrix representation
integer, parameter :: nupdn = 3
real*8,  allocatable, dimension(:,:,:,:) :: dnmxr, dnmxi
real*8,  allocatable, dimension(:,:,:) :: kr, ki
real*8,  allocatable, dimension(:,:) :: dmtr, dmti
real*8  :: rkc(4)

! for density matrix representation for exiton diffusion
real*8,  allocatable, dimension(:,:) :: accumprob_ex
integer, parameter :: nupdn_ex = 20
real*8,  allocatable, dimension(:,:,:) :: dnmxr_ex, dnmxi_ex
real*8,  allocatable, dimension(:,:,:) :: kr_ex, ki_ex
real*8,  allocatable, dimension(:,:) :: dmtr_ex, dmti_ex
!real*8  :: rkc(4)
integer :: nlrstates
integer :: nlrreduce = 0
real*8,  allocatable, dimension(:) :: excite
real*8,  allocatable, dimension(:,:) :: weight
integer, allocatable, dimension(:,:) :: ibstate
logical :: lset_ibstate_rec = .false.
integer, allocatable, dimension(:,:,:) :: ibstate_rec
integer, allocatable, dimension(:) :: state_ref
logical, allocatable, dimension(:) :: lexcite
real*8,  allocatable, dimension(:,:) :: dkj_ex
integer, allocatable, dimension(:,:) :: nnews
integer, allocatable, dimension(:) :: iud_ex
real*8,  allocatable, dimension(:,:) :: gamma_dipole

integer, allocatable, dimension(:,:,:) :: iunit_prob

real*8  :: rndx

real*8,  allocatable, dimension(:,:) :: occ_store
real*8,  allocatable, dimension(:,:) :: eig_store
real*8,  allocatable, dimension(:) :: rhgr_store
real*8,  allocatable, dimension(:,:,:) :: rhoij_store
real*8,  allocatable, dimension(:,:,:,:) :: rhomm_store
integer :: lfmax, mmx
integer :: nionx = 0
real*8,  allocatable, dimension(:,:) :: gdcr_rec
logical :: lset_gdcr_rec = .false.

!-----variables for Gaussian broadening
integer, parameter :: ntbfi = 10000
real*8,  allocatable, dimension(:) :: tbfi ,tbfia
real*8 :: tbxe
real*8 :: tbmx
real*8 :: piru

real*8  :: theekin = 0.d0 ! the current kinetic energy of ions
real*8  :: thepote = 0.d0 ! the current potential energy
logical :: lsetekin = .false.
integer, parameter :: nprevex = 3
real*8  :: totene(nprevex) = 0.d0 ! the previous total energies
integer :: npreve = 0
logical :: lreject = .false.
real*8,  allocatable, dimension(:,:) :: couplvec  ! non-adiabatic coupling vector

save

end module




subroutine get_ltddft( ltddft_, ltddft_fssh_ )
!-----------------------------------------------------------------------
use tddft_variables
implicit none
logical :: ltddft_, ltddft_fssh_

ltddft_      = ltddft
ltddft_fssh_ = ltddft_fssh

return
end subroutine





subroutine set_ltddft( nfile, myid, nodes, &
& ltddft_, ltddft_fssh_, ltddft_start_, nel_, dttddft_, lfssh_switch_, &
& lfssh_gsscf_, tdbroad_, lfssh_random_, rseed_fssh, lfssh_boltzmn_, treq_, &
& lfssh_vscale_, tminimum_, lfssh_parallel_, ltddft_nscforce_, lrtddft_ )
!-----------------------------------------------------------------------
!     set ltddft
!-----------------------------------------------------------------------
use constants
use tddft_variables
use tddft_fssh_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
logical :: ltddft_, ltddft_fssh_, ltddft_start_, lfssh_switch_, lfssh_gsscf_
logical :: lfssh_random_, lfssh_parallel_, ltddft_nscforce_
integer :: nel_
real*8  :: dttddft_, tdbroad_, rseed_fssh, treq_
logical :: lfssh_boltzmn_, lfssh_vscale_, lrtddft_
real*8  :: tminimum_

!-----declare local variables
integer :: i
integer :: date_time(8)
character(10) :: big_ben(3)
real*8  :: a1, a2, a3, rn


ltddft       = ltddft_
ltddft_fssh  = ltddft_fssh_
ltddft_start = ltddft_start_

nel = nel_

dttddft = dttddft_

if( ldensity_matrix_repre ) lfssh_parallel_ = .false.
lfssh_parallel = lfssh_parallel_

ltddft_nscforce = ltddft_nscforce_

lfssh_switch = lfssh_switch_

lfssh_gsscf = lfssh_gsscf_

tdbroad = tdbroad_
ltdbroad = tdbroad > 1.d-05

lfssh_random = lfssh_random_

lfssh_boltzmn = lfssh_boltzmn_
treq = treq_

lfssh_vscale = lfssh_vscale_
tminimum = tminimum_ * 3.d0/2.d0/(tempau/2.d0)   !  [K] -> [Ryd]

lrtddft = lrtddft_


!---initialize random number
if( lfssh_random ) then
    !---set manually
    rndx = rseed_fssh
  else
    !---set automatically
    call date_and_time( big_ben(1), big_ben(2), big_ben(3), date_time )
    !        print *,'date_time array values:'
    !        print *,'year=',date_time(1)
    !        print *,'month_of_year=',date_time(2)
    !        print *,'day_of_month=',date_time(3)
    !        print *,'time difference in minutes=',date_time(4)
    !        print *,'hour of day=',date_time(5)
    !        print *,'minutes of hour=',date_time(6)
    !        print *,'seconds of minute=',date_time(7)
    !        print *,'milliseconds of second=',date_time(8)
    !        print *, 'DATE=',big_ben(1)
    !        print *, 'TIME=',big_ben(2)
    !        print *, 'ZONE=',big_ben(3)

    a3 = date_time(8)/100
    a2 = date_time(8)/10 - a3*10
    a1 = mod(date_time(8),10)
    rndx = a2*10.d0**a1 + a3*6.d0**a2 + a1*3.d0**a3
    call rnd00( rn, rndx )
    do i = 1, int(rn*100000)
       call rnd00( rn, rndx )
       !---for DISH
       call rndexp( rn, 1.d0 )
    end do
    rseed_fssh = rndx
end if


return
end subroutine




subroutine set_ltddft_fssh( ltddft_fssh_, lfssh_gsscf_ )
!-----------------------------------------------------------------------
use tddft_variables
implicit none
logical :: ltddft_fssh_, lfssh_gsscf_

lfssh_gsscf = lfssh_gsscf_
ltddft_fssh = ltddft_fssh_

return
end subroutine




subroutine tddft_fssh_alloc( nfile, myid, nodes,  &
& alloc_mem, nband, nbnod, noband_, nspnmx, natom )
!-----------------------------------------------------------------------
!     allocate memory for TDDFT-FSSH
!-----------------------------------------------------------------------
use tddft_variables
use tddft_fssh_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
real*8  :: alloc_mem
integer :: nband, nbnod, noband_, nspnmx, natom, nb4

!-----declare local variables
real*8  :: the_mem
integer :: i
integer :: status

if( .not.ltddft_fssh ) return


lspin  = nspnmx == 2
noband = noband_

if( .not.lrtddft ) then
!    iudx = 2
    iudx = nupdn
else
    iudx = nupdn_ex
end if

!-----allocate arrays
allocate( occ_rec(nband,nspnmx), &
& phisphi(nband,nband,nspnmx), lsignc(nband,nspnmx),  &
& bmthr(nband,nband), bmthi(nband,nband),  &
& iunit_prob(-iudx:iudx,nband,nspnmx), &
& stat=status )

the_mem =  &
& 8.d0 * ( size(occ_rec) + size(phisphi) + size(bmthr) + size(bmthi) )  &
& + 4.d0 * ( size(iunit_prob) )  &
& + 1.d0 * ( size(lsignc) )

if( .not.lrtddft .and. .not.ldensity_matrix_repre .and. status == 0 ) then
    allocate( bmtar(nband,nband), bmtai(nband,nband), &
& lmtr(nband,nband), lmti(nband,nband), umtr(nband,nband), umti(nband,nband), &
& stat=status )
    the_mem = the_mem  &
& + 8.d0 * ( size(bmtar) + size(bmtai) &
& + size(lmtr) + size(lmti) + size(umtr) + size(umti) )
end if

if( .not.lrtddft .and. status == 0 ) then
    allocate( accumprob(nband,nband,nspnmx),  &
& stat=status )
    the_mem = the_mem  &
& + 8.d0 * ( size(accumprob) )
end if

if( lfssh_parallel .and. status == 0 ) then
    allocate(  &
& nbcnt(nodes), nbdsp(nodes), &
& stat=status )
    the_mem = the_mem  &
& + 4.d0 * ( size(nbcnt) + size(nbdsp) )
end if

if( status == 0 ) then
    if( lrtddft ) then
        !---exciton dynamics
        nb4 = 4
        allocate(  &
& dnmxr_ex(0:nupdn_ex,0:nupdn_ex,nexciton),  &
& dnmxi_ex(0:nupdn_ex,0:nupdn_ex,nexciton),  &
& kr_ex(0:nupdn_ex,0:nupdn_ex,nb4),  &
& ki_ex(0:nupdn_ex,0:nupdn_ex,nb4),  &
& dmtr_ex(0:nupdn_ex,0:nupdn_ex),  &
& dmti_ex(0:nupdn_ex,0:nupdn_ex),  &
& probinteg(0:nupdn_ex,nexciton), &
& nnews(0:nupdn_ex,nexciton), &
& iband_hole_(nexciton), iband_electron_(nexciton),  &
& ldegenerate_(nexciton), numexciton(nexciton), ldmcalc_ex(nexciton),  &
& iud_ex(nexciton),  &
& stat=status )
        the_mem = the_mem  &
& + 8.d0 * ( size(dnmxr_ex) + size(dnmxi_ex) + size(kr_ex) + size(ki_ex)  &
&          + size(dmtr_ex) + size(dmti_ex) + size(probinteg) )  &
& + 4.d0 * ( size(nnews) + size(iband_hole_)*4 )  &
& + 1.d0 * size(ldegenerate_)*2
    else

        if( ldensity_matrix_repre ) then
            nb4 = max( nband, 4 )
            allocate(  &
& dnmxr(-nupdn:nupdn,-nupdn:nupdn,nband,nspnmx),  &
& dnmxi(-nupdn:nupdn,-nupdn:nupdn,nband,nspnmx),  &
& kr(-nupdn:nupdn,-nupdn:nupdn,nb4), ki(-nupdn:nupdn,-nupdn:nupdn,nb4),  &
& dmtr(-nupdn:nupdn,-nupdn:nupdn), dmti(-nupdn:nupdn,-nupdn:nupdn),  &
& probinteg(-iudx:iudx,nband), &
& stat=status )
            the_mem = the_mem  &
& + 8.d0 * ( size(dnmxr) + size(dnmxi) + size(kr) + size(ki)  &
&          + size(dmtr) + size(dmti) + size(probinteg) )
        else
            allocate(  &
& cr(nband,nband,nspnmx), ci(nband,nband,nspnmx),  &
& amtr(nband,nband), amti(nband,nband),  &
& cmtr(nband,nband), cmti(nband,nband),  &
& wrkr(nband,nband), wrki(nband,nband), &
& probinteg(-iudx:iudx,nband), &
& stat=status )
            the_mem = the_mem  &
& + 8.d0 * ( size(cr) + size(ci) + size(amtr) + size(amti)  &
& + size(cmtr) + size(cmti) + size(wrkr) + size(wrki) + size(probinteg) )
        end if

    end if
end if

if( ltdbroad .and. status == 0 ) then
    allocate( tbfi(0:ntbfi) ,tbfia(0:ntbfi), &
& stat=status )
    the_mem = the_mem + 8.d0*( size(tbfi) + size(tbfia) )
end if

!if( lfssh_gsscf .and. status == 0 ) then
if( status == 0 ) then
    allocate( occ_store(nband,nspnmx), eig_store(nband,nspnmx), &
& stat=status )
    the_mem = the_mem + 8.d0*( size(occ_store) + size(eig_store) )
end if

if( lfssh_vscale .and. status == 0 ) then
    allocate( couplvec(3,natom), &
& stat=status )
    the_mem = the_mem + 8.d0*( size(couplvec) )
end if

if( lfssh_dish .and. status == 0 ) then
    allocate( decoherence_rate(nband,nband), dish_prev_step(nband,nspnmx), &
& stat=status )
    the_mem = the_mem + 8.d0*( size(decoherence_rate) )  &
& + 4.d0 * size(dish_prev_step)
end if


!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'tddft_fssh_alloc', .true. )

occ_rec(:,:) = 0.d0
iunit_prob(-iudx:iudx,1:nband,1:nspnmx) = 0
lsignc(1:nband,1:nspnmx) = .false.

if( lrtddft ) then
    !---exciton dynamics
    dnmxr_ex = 0.d0
    dnmxi_ex = 0.d0
    dnmxr_ex(0,0,1:nexciton) = 1.d0
    rkc(1) = 0.d0
    rkc(2) = dttddft/2.d0
    rkc(3) = dttddft/2.d0
    rkc(4) = dttddft
else
    accumprob = 0.d0
    if( ldensity_matrix_repre ) then
        !--- density matrix representation
        dnmxr = 0.d0
        dnmxi = 0.d0
        dnmxr(0,0,1:nband,1:nspnmx) = 1.d0
        rkc(1) = 0.d0
        rkc(2) = dttddft/2.d0
        rkc(3) = dttddft/2.d0
        rkc(4) = dttddft
    else
        cr = 0.d0
        ci = 0.d0
        do i = 1, nband
           cr(i,i,1:nspnmx) = 1.d0
        end do
        amtr(1:nband,1:nband) = 0.d0
        amti(1:nband,1:nband) = 0.d0
        cmtr(1:nband,1:nband) = 0.d0
        cmti(1:nband,1:nband) = 0.d0
        wrkr(1:nband,1:nband) = 0.d0
        wrki(1:nband,1:nband) = 0.d0
    end if
end if

if( lfssh_dish ) then
    decoherence_rate(:,:) = 0.d0
    do i = 1, ndishpair
       decoherence_rate(ndishi(i),ndishj(i)) = temp_decoherence_rate(i) * 2.d0
       decoherence_rate(ndishj(i),ndishi(i)) = temp_decoherence_rate(i) * 2.d0
                                            ! [Hartree a.u.} -> [Rydberg a.u.]
    end do
    dish_prev_step(:,:) = -1
    !--- deallocate memory
    deallocate( ndishi, ndishj, temp_decoherence_rate,  &
& stat=status )
end if


!----- initial set for Gaussian broadening -----
if( ltdbroad ) then
    call setgbr( tbfi ,tbfia, ntbfi, tbxe, tbmx, piru )
end if


return
end subroutine




subroutine tddft_fssh_pw_alloc( nfile, myid, nodes, &
& alloc_mem, nspnmx, nplw5ex, pwscale, lpnbndx )
!-----------------------------------------------------------------------
! allocate memory for the plane wave method in TDDFT-FSSH
!-----------------------------------------------------------------------
use tddft_variables
use tddft_fssh_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
real*8  :: alloc_mem
integer :: nspnmx, nplw5ex
real*8  :: pwscale
integer :: lpnbndx

!-----declare local variables
integer :: ngenhex    = 0    ! nplw5ex
integer :: lpnbndxx   = 0
real*8  :: the_mem
integer :: status
save ngenhex, lpnbndxx

!if( .not.ltddft_fssh .or. .not.lfssh_gsscf ) return
if( .not.ltddft_fssh ) return


if( lfssh_gsscf ) then
if( nplw5ex > ngenhex ) then

    !-----if already allocated, deallocate arrays
    if( allocated(rhgr_store) ) then
        the_mem =  8.d0 * size(rhgr_store)
        !-----deallocate arrays
        deallocate( rhgr_store, stat=status )

        !------error trap
        call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'tddft_fssh_pw_alloc', .true. )
     end if

     !-----allocate arrays
     ngenhex    = nplw5ex* pwscale
     allocate( rhgr_store(ngenhex*nspnmx), stat=status )

     the_mem =  8.d0 * size(rhgr_store)

     !------error trap
     call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'tddft_fssh_pw_alloc', .true. )

end if
end if


if( lpnbndx > lpnbndxx ) then

    !-----if already allocated, deallocate arrays
    if( allocated(gdcr_rec) ) then
        the_mem =  8.d0 * size(gdcr_rec)
        !-----deallocate arrays
        deallocate( gdcr_rec, stat=status )

        lset_gdcr_rec = .false.

        !------error trap
        call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'tddft_fssh_pw_alloc (2)', .true. )
     end if

     !-----allocate arrays
     lpnbndxx    = lpnbndx
     allocate( gdcr_rec(lpnbndx,1:nspnmx), stat=status )

     the_mem =  8.d0 * size(gdcr_rec)

     !------error trap
     call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'tddft_fssh_pw_alloc (2)', .true. )

end if


return
end subroutine




subroutine set_occupations( nfile, myid, nodes, iunit, nocc_change_ )
!-----------------------------------------------------------------------
!     set electronic occupations
!-----------------------------------------------------------------------
use tddft_variables
implicit none

integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: iunit
integer :: nocc_change_

!-----declare local variables
integer :: ib
integer :: status, istat


nocc_change = nocc_change_
if( nocc_change <= 0 ) return

!--- allocate memory
allocate( numband(nocc_change), occ_new(nocc_change,2), stat=status )

do ib = 1, nocc_change
   read(iunit,*,iostat=istat) numband(ib), occ_new(ib,1), occ_new(ib,2)
   !----- error trap
   if( istat /= 0 ) call fstop( nfile, myid, nodes,          &
                                'error-0001 in set_occupations' )
end do


return
end subroutine




subroutine out_occupations( nfile, lspin )
!-----------------------------------------------------------------------
!     set electronic occupations
!-----------------------------------------------------------------------
use tddft_variables
implicit none
integer :: nfile
logical :: lspin

!-----declare local variables
integer :: ib


write(nfile,*) ' '

if( lrtddft ) then

    write(nfile,*) ' Exciton in TDDFT-FSSH (only available when ltddft_start=.false.)'
    write(nfile,*) '  nexciton : ', nexciton
    write(nfile,*) '      iband_hole, iband_electron, ldegenerate'
    do ib = 1, nexciton
       write(nfile,'(i5,a1,2i8,l6)') ib, ':',  &
& iband_hole(ib), iband_electron(ib), ldegenerate(ib)
    end do

else

    write(nfile,*) ' Occupation changes in TDDFT (only available when ltddft_start=.false.)'
    write(nfile,*) '  nocc_change : ', nocc_change

    if( nocc_change > 0 ) then
        if( .not.lspin ) then
            write(nfile,3100)
            do ib = 1, nocc_change
               write(nfile,3110) numband(ib), occ_new(ib,1)
            end do
          else
            write(nfile,3200)
            do ib = 1, nocc_change
               write(nfile,3110) numband(ib), occ_new(ib,1), occ_new(ib,2)
            end do
        end if

 3100 format(2x,' numband        occ_new')
 3200 format(2x,' numband    occ_new(up)    occ_new(down)')
 3110 format(i8,2f17.4)
    end if

end if


return
end subroutine




subroutine tddft_fssh_ini( nfile, myid, nodes,  &
& nband, nbnod, nspnmx, eig, norder, occ, adde_fssh )
!-----------------------------------------------------------------------
!    Initial occupancies in TDDFT-FSSH
!-----------------------------------------------------------------------
use outfile
use tddft_variables
use tddft_fssh_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: nband, nbnod, nspnmx
real*8,  dimension(nband,nspnmx) :: eig, occ
integer, dimension(nband,nspnmx) :: norder
real*8  :: adde_fssh

!-----declare local variables
integer :: ib, nspin
real*8  :: sum


!-----Number of wave packets : noccmx
nbase = nel/2
nbase = nbase + nel - 2*nbase
!noccmx = nbase
!do ib = 1, nocc_change
!   noccmx = max( noccmx, numband(ib) )
!end do
noccmx = noband


if( ltddft_start ) then
!    if( ltdbroad ) then
!        call broadening_in_tddft_fssh2( nfile, myid, nodes,  &
!& nband, nbnod, nspnmx, eig, norder, occ, adde_fssh )
!    end if
    return
end if


!---set occupations of wave packets
occ_rec(:,:) = 0.d0
do nspin = 1, nspnmx
   do ib = 1, nbase
      occ_rec(ib,nspin) = 2.d0/dble(nspnmx)
   end do
   if( lrtddft ) then
       if( .not.lspin ) then
           do ib = 1, nexciton
           if( iband_hole(ib) > 0 ) then
               occ_rec(    iband_hole(ib),nspin) = occ_rec(    iband_hole(ib),nspin) - 1.d0
               occ_rec(iband_electron(ib),nspin) = occ_rec(iband_electron(ib),nspin) + 1.d0
           end if
           end do
       else
           do ib = 1, nexciton
           if( iband_hole(ib) > 0 ) then
           if( ldegenerate(ib) .and. nspin == 2 .or. .not.ldegenerate(ib) .and. nspin == 1 ) then
               occ_rec(    iband_hole(ib),nspin) = occ_rec(    iband_hole(ib),nspin) - 1.d0
           end if
           if( nspin == 1 ) then
               occ_rec(iband_electron(ib),nspin) = occ_rec(iband_electron(ib),nspin) + 1.d0
           end if
           end if
           end do
       end if
   else
       do ib = 1, nocc_change
          occ_rec(numband(ib),nspin) = occ_new(ib,nspin)
       end do
   end if
end do
!---check sum of occupations
sum = 0.d0
do nspin = 1, nspnmx
   do ib = 1, noccmx
      sum = sum + occ_rec(ib,nspin)
   end do
end do
if( abs(sum - dble(nel)) > 1.d-05 ) then
    if(loutfile(1)) write(nfile(1),*) 'sum, nel =', sum, nel
    if(loutfile(2)) write(nfile(2),*) 'sum, nel =', sum, nel
    call fstop( nfile, myid, nodes, 'Error: sum /= nel in tddft_fssh_ini' )
end if

!---save occupations
if( .not.lfssh_gsscf ) then
    occ = occ_rec
end if


return
end subroutine




subroutine tddft_fssh_exchg_occ( nfile, myid, nodes,  &
& nband, nbnod, nspnmx, eig, norder, occ )
!-----------------------------------------------------------------------
!    Exchange occupancies in TDDFT-FSSH
!-----------------------------------------------------------------------
use tddft_variables
use tddft_fssh_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: nband, nbnod, nspnmx
real*8,  dimension(nband,nspnmx) :: eig, occ
integer, dimension(nband,nspnmx) :: norder


occ_store(1:nband,1:nspnmx) = occ(1:nband,1:nspnmx)  ! ground -state occupations
occ(1:nband,1:nspnmx) = occ_rec(1:nband,1:nspnmx)    ! excited-state occupations

eig_store(1:nband,1:nspnmx) = eig(1:nband,1:nspnmx)  ! ground -state eigenvalues


return
end subroutine




subroutine tddft_fssh_restore_occ( nfile, myid, nodes,  &
& nband, nbnod, nspnmx, eig, norder, occ )
!-----------------------------------------------------------------------
!    Restore ground-state occupancies in TDDFT-FSSH
!-----------------------------------------------------------------------
use tddft_variables
use tddft_fssh_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: nband, nbnod, nspnmx
real*8,  dimension(nband,nspnmx) :: eig, occ
integer, dimension(nband,nspnmx) :: norder


occ(1:nband,1:nspnmx) = occ_store(1:nband,1:nspnmx)  ! ground -state occupations
eig(1:nband,1:nspnmx) = eig_store(1:nband,1:nspnmx)  ! ground -state eigenvalues


return
end subroutine




subroutine tddft_fssh_store_rhgr( nfile, myid, nodes,  &
& rhgr, nplw5ex, nspnmx )
!-----------------------------------------------------------------------
!    Store charge density in TDDFT-FSSH
!-----------------------------------------------------------------------
use tddft_variables
use tddft_fssh_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: nplw5ex, nspnmx
real*8,  dimension(nplw5ex*nspnmx) :: rhgr


rhgr_store(1:nplw5ex*nspnmx) = rhgr(1:nplw5ex*nspnmx)


return
end subroutine




subroutine tddft_fssh_restore_rhgr( nfile, myid, nodes,  &
& rhgr, nplw5ex, nspnmx )
!-----------------------------------------------------------------------
!    Restore charge density in TDDFT-FSSH
!-----------------------------------------------------------------------
use tddft_variables
use tddft_fssh_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: nplw5ex, nspnmx
real*8,  dimension(nplw5ex*nspnmx) :: rhgr


rhgr(1:nplw5ex*nspnmx) = rhgr_store(1:nplw5ex*nspnmx)


return
end subroutine




subroutine broadening_in_tddft_fssh( nfile, myid, nodes,  &
& nband, nbnod, nspnmx, eig, norder, occ, adde_fssh )
!-----------------------------------------------------------------------
! Gaussian broadening of ccupancies in TDDFT-FSSH
!-----------------------------------------------------------------------
use tddft_variables
use tddft_fssh_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: nband, nbnod, nspnmx
real*8,  dimension(nband,nspnmx) :: eig, occ
integer, dimension(nband,nspnmx) :: norder
real*8  :: adde_fssh


if( ltdbroad ) then
    call broadening_in_tddft_fssh2( nfile, myid, nodes,  &
& nband, nbnod, nspnmx, eig, norder, occ, adde_fssh )
end if


return
end subroutine




subroutine broadening_in_tddft_fssh2( nfile, myid, nodes,  &
& nband, nbnod, nspnmx, eig, norder, occ, adde_fssh )
!-----------------------------------------------------------------------
! Gaussian broadening
!-----------------------------------------------------------------------
use outfile
use tddft_variables
use tddft_fssh_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: nband, nbnod, nspnmx
real*8,  dimension(nband,nspnmx) :: eig, occ
integer, dimension(nband,nspnmx) :: norder
real*8  :: adde_fssh

!-----declare local variables
integer :: ib, nspin, m
real*8  :: a, di, dim, dip, dm, d, di1, di2, dimo, dipo
real*8  :: sum, pairu, aexp


pairu = -tdbroad/sqrt(acos(-1.d0))/2.d0
adde_fssh = 0.d0
do nspin = 1, nspnmx
   do ib = 1, noccmx
      if( ib == 1 ) then
          di1 = 1.d0
        else
          a = (eig(ib-1,nspin)-eig(ib,nspin))/tdbroad
          if( a > tbmx) then
              di1 = 0.d0
          else if( a < -tbmx ) then
              di1 = 1.d0
          else
              dm = (a+tbmx)/tbxe
              m  = 0.5d0*dm
              m  = 2*m
              d  = 0.5d0*( dm - dble(m) )
              di1 = d*( (d-1.d0)*tbfia(m) + tbfi(m+2) - tbfi(m) ) + tbfi(m)
          end if
      end if
      if( ib == noccmx ) then
          di2 = 0.d0
        else
          a = (eig(ib+1,nspin)-eig(ib,nspin))/tdbroad
          if( a > tbmx) then
              di2 = 0.d0
          else if( a < -tbmx ) then
              di2 = 1.d0
          else
              dm = (a+tbmx)/tbxe
              m  = 0.5d0*dm
              m  = 2*m
              d  = 0.5d0*( dm - dble(m) )
              di2 = d*( (d-1.d0)*tbfia(m) + tbfi(m+2) - tbfi(m) ) + tbfi(m)
          end if
      end if
      di = (di1 - di2)*occ_rec(ib,nspin)

      dimo = 0.d0
      if( ib > 1 ) then
          a = (eig(ib,nspin)-eig(ib-1,nspin))/tdbroad
          if( a > tbmx) then
              dim = 0.d0
          else if( a < -tbmx ) then
              dim = 1.d0
          else
              dm = (a+tbmx)/tbxe
              m  = 0.5d0*dm
              m  = 2*m
              d  = 0.5d0*( dm - dble(m) )
              dim = d*( (d-1.d0)*tbfia(m) + tbfi(m+2) - tbfi(m) ) + tbfi(m)
          end if

          dimo = dim*occ_rec(ib-1,nspin)

          !---additional energy
!          a = a*a
!          if( a.lt.3.3d+01 ) then
!              aexp = exp(-a)
!            else
!              aexp = 0.d0
!          end if
!          adde_fssh = adde_fssh + (occ_rec(ib,nspin)-occ_rec(ib-1,nspin)) &
!&                   * ( pairu*aexp + eig(ib-1,nspin)*dim )

      end if

      dipo = 0.d0
      if( ib < noccmx ) then
          a = (eig(ib+1,nspin)-eig(ib,nspin))/tdbroad
          if( a > tbmx) then
              dip = 0.d0
          else if( a < -tbmx ) then
              dip = 1.d0
          else
              dm = (a+tbmx)/tbxe
              m  = 0.5d0*dm
              m  = 2*m
              d  = 0.5d0*( dm - dble(m) )
              dip = d*( (d-1.d0)*tbfia(m) + tbfi(m+2) - tbfi(m) ) + tbfi(m)
          end if

          dipo = dip*occ_rec(ib+1,nspin)

          !---additional energy
!          a = a*a
!          if( a.lt.3.3d+01 ) then
!              aexp = exp(-a)
!            else
!              aexp = 0.d0
!          end if
!          adde_fssh = adde_fssh + (occ_rec(ib+1,nspin)-occ_rec(ib,nspin)) &
!&                   * ( pairu*aexp + eig(ib+1,nspin)*dip )

      end if

      occ(ib,nspin) = di + dimo + dipo

   end do
end do
adde_fssh = -adde_fssh


!---check sum of occupations
sum = 0.d0
do nspin = 1, nspnmx
   do ib = 1, noccmx
      sum = sum + occ(ib,nspin)
   end do
end do
if( abs(sum - dble(nel)) > 1.d-12 ) then
    if(loutfile(1)) write(nfile(1),*) 'sum, nel =', sum, nel
    if(loutfile(2)) write(nfile(2),*) 'sum, nel =', sum, nel
    call fstop( nfile, myid, nodes, 'Error: sum /= nel in tddft_fssh_ini' )
end if


return
end subroutine




module nspin_in_tddft_fssh
!-----------------------------------------------------------------------
! type declaration of shared variables in tddft.f90
!-----------------------------------------------------------------------
implicit none

integer :: nspin

save

end module




subroutine set_nspin_in_tddft_fssh( nspin_ )
use nspin_in_tddft_fssh
implicit none
integer :: nspin_

nspin = nspin_

return
end subroutine




subroutine set_ekin_in_tddft_fssh( theekin_ )
use tddft_fssh_variables
implicit none
real*8 :: theekin_

theekin = theekin_
lsetekin = .true.

return
end subroutine




subroutine set_pote_in_tddft_fssh( thepote_ )
use tddft_fssh_variables
implicit none
real*8 :: thepote_

thepote = thepote_

return
end subroutine




subroutine tddft_fssh_check_unitary( nfile, myid, nodes,  &
& bijr, nband, nbnod, nbxxxx, prod, prodr )
!-----------------------------------------------------------------------
! Check wavefunction exchange by Unitary transformation
!-----------------------------------------------------------------------
use outfile
use tddft_variables
use tddft_fssh_variables
use nspin_in_tddft_fssh
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: nband, nbnod, nbxxxx
real*8  :: bijr(nbxxxx,*)
real*8,  dimension(nband) :: prod, prodr

!-----declare local variables
integer :: ib, jb, ibb, jbb
real*8  :: umatmax, tmpv


if( .not.ltddft_fssh ) return


do ib = 1, nband
   umatmax = 0.d0
   do jb = 1, nband
      if( umatmax < abs(bijr(jb,ib)) ) then
          umatmax = abs(bijr(jb,ib))
          prod(ib) = jb
      end if
   end do

   !---Do not change signs to calculate transfer matrix by finite differencing
   jbb = nint(prod(ib))
   if( bijr(jbb,ib) < 0.d0 ) then
       do jb = 1, nband
          bijr(jb,ib) = - bijr(jb,ib)
       end do
    end if
end do


if( .not.lcheck_unitary ) return

prodr(1:nband) = 0.d0
do ib = 1, nband
   if( prodr(ib) < 0.5d0 ) then
      do
         ibb = nint(prod(ib))
         if( ibb == ib ) exit
         if( abs(occ_rec(ibb,nspin)-occ_rec(ib,nspin)) < 1.d-05 ) exit
         if( prodr(ibb) > 0.5d0 ) then
             call fstop( nfile, myid, nodes, &
     &                  'TDDFT-FSSH: strange Unitary matrix in untryt' )
         end if
         do jb = 1, nband
            tmpv = bijr(jb,ib)
            bijr(jb,ib)  = bijr(jb,ibb)
            bijr(jb,ibb) = tmpv
         end do
if(loutfile(1)) write(nfile(1),*) '*** information exchange column of Unitrary matrix', ib, ibb
if(loutfile(2)) write(nfile(2),*) '*** information exchange column of Unitrary matrix', ib, ibb
         prod(ib)  = prod(ibb) 
         prod(ibb) = ibb
         prodr(ibb) = 1.d0
      end do
      prodr(ib) = 1.d0
   end if
end do


return
end subroutine




subroutine tddft_fssh_check_trimat( nfile, myid, nodes,  &
& dmtrxr, nband, nbnod, nbxxxx, prod, prod2 )
!-----------------------------------------------------------------------
! Check wavefunction exchange by triangular matrix in Schmitd
!-----------------------------------------------------------------------
use outfile
use tddft_variables
use tddft_fssh_variables
use nspin_in_tddft_fssh
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: nband, nbnod, nbxxxx
real*8  :: dmtrxr(nbxxxx,*)
real*8,  dimension(nband) :: prod, prod2

!-----declare local variables
integer :: ib, jb, ibb, jbb
real*8  :: umatmax, tmpv


if( .not.ltddft_fssh ) return
if( .not.lcheck_trimat ) return


do ib = 1, nband
   umatmax = 0.d0
   do jb = 1, nband
      if( umatmax < abs(dmtrxr(ib,jb)) ) then
          umatmax = abs(dmtrxr(ib,jb))
          prod(ib) = jb
      end if
   end do
end do


prod2(1:nband) = 0.d0
do ib = 1, nband
   if( prod2(ib) < 0.5d0 ) then
      do
         ibb = nint(prod(ib))
         if( ibb == ib ) exit
         if( abs(occ_rec(ibb,nspin)-occ_rec(ib,nspin)) < 1.d-05 ) exit
         if( prod2(ibb) > 0.5d0 ) then
             call fstop( nfile, myid, nodes, &
     &                  'TDDFT-FSSH: strange triangular matrix in schmidt' )
         end if
         do jb = 1, nband
            tmpv = dmtrxr(ib,jb)
            dmtrxr(ib,jb)  = dmtrxr(ibb,jb)
            dmtrxr(ibb,jb) = tmpv
         end do
if(loutfile(1)) write(nfile(1),*) '*** information exchange low of triangular matrix', ib, ibb
if(loutfile(2)) write(nfile(2),*) '*** information exchange low of triangular matrix', ib, ibb
         prod(ib)  = prod(ibb) 
         prod(ibb) = ibb
         prod2(ibb) = 1.d0
      end do
      prod2(ib) = 1.d0
   end if
end do


return
end subroutine




subroutine tddft_fssh_set_matrix( nfile, myid, nodes,  &
& gdcr, nplw, nplwex, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& lspin_, gdcrsv, nspnod, lcgjsv, lvand, lvandi, ntype, rvol,  &
& aijr, nbxxxx, prod, prodr, prodx )
!-----------------------------------------------------------------------
! Set wavefunction matrix
!-----------------------------------------------------------------------
use outfile
use tddft_variables
use tddft_fssh_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: nplw, nplwex, npnod1, npnod2, npnod, nplcnt(*), npldsp(*)
integer :: nband, nbnod1, nbnod2, nbnod, nbncnt(*), nbndsp(*)
real*8  :: gdcr(2*npnod,*)
logical :: lspin_, lcgjsv
real*8  :: gdcrsv(*)
integer :: nspnod
integer :: ntype
logical :: lvand
logical, dimension(ntype) :: lvandi
real*8  :: rvol
integer :: nbxxxx
real*8  :: aijr(nbxxxx,*)
real*8  :: prod(nband), prodr(nband), prodx(nband,nband)

!---declare local variables
integer :: nnspin, nspin


if( .not.ltddft_fssh ) return


if( lspin ) then
    nnspin = 2
  else
    nnspin = 1
end if

if( lset_gdcr_rec ) then

do nspin = 1, nnspin
   if( lspin ) then
       call stspud( nspin, gdcr, gdcrsv, nspnod, lcgjsv )
!       if( lvand ) call ldslmi( nfile, myid, nodes, nspin, lvandi, .false. )
   end if

   call tddft_fssh_cal_matrix( nfile, myid, nodes,  &
& gdcr, gdcr_rec(1,nspin), nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& npnod1, npnod2, npnod, prod, prodr, lvand, rvol, nspin,  &
& aijr, nbxxxx )

    call tddft_fssh_set_matrix2( nfile, myid, nodes,  &
& aijr, nband, nbnod, nbxxxx, nspin, prod, prodr, prodx )
end do

end if


do nspin = 1, nnspin
   if( lspin ) then
       call stspud( nspin, gdcr, gdcrsv, nspnod, lcgjsv )
!       if( lvand ) call ldslmi( nfile, myid, nodes, nspin, lvandi, .false. )
   end if

   call dcopy_a_to_b( gdcr, gdcr_rec(1,nspin), nspnod )

!   if( lvand ) call sv_slmir_atm( nfile, myid, nodes, nband, nspin )

end do

lset_gdcr_rec = .true.


return
end subroutine




subroutine tddft_fssh_cal_matrix( nfile, myid, nodes,  &
& gdcr, gdcr_rec, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& npnod1, npnod2, npnod, prod, prodr, lvand, rvol, nspin,  &
& aijr, nbxxxx )
!-----------------------------------------------------------------------
!     subspace alignment by  cgjr and rhcr
!-----------------------------------------------------------------------
use outfile
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: npnod1, npnod2, npnod
real*8  :: gdcr(2*npnod,*)
real*8  :: gdcr_rec(2*npnod,*)
integer :: nband, nbnod1, nbnod2, nbnod, nbncnt(*), nbndsp(*)
real*8  :: prod(*), prodr(*)
logical :: lvand
real*8  :: rvol
integer :: nspin
integer :: nbxxxx
real*8  :: aijr(nbxxxx,*)

!---declare local variables
integer :: root, i, j, m, mshnod
real*8  :: cscr


mshnod = 2*npnod

root = 0
do i = 1, nband
   if( i.gt.nbndsp(root+1)+nbncnt(root+1) ) root = root + 1
   do j = 1, nband
      prod(j) = 0.d0
      do m = 1, mshnod
         prod(j) = prod(j) + gdcr_rec(m,i)*gdcr(m,j)
      enddo
      if( myid.eq.0 ) then
          prod(j) = prod(j)*2.d0 - gdcr_rec(1,i)*gdcr(1,j)
        else
          prod(j) = prod(j)*2.d0
      endif
   enddo
!   if( lvand ) then
!       do j = 1, nband
!          call calpcscc_tddft( nfile, myid, nodes,  &
!& cscr, i, j, rvol, nspin )
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


return
end subroutine




subroutine tddft_fssh_set_matrix2( nfile, myid, nodes,  &
& aijr, nband, nbnod, nbxxxx, nspin, prod, prodr, prodx )
!-----------------------------------------------------------------------
! Set wavefunction matrix
!-----------------------------------------------------------------------
use outfile
use tddft_variables
use tddft_fssh_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: nband, nbnod, nbxxxx, nspin
real*8  :: aijr(nbxxxx,*)
real*8  :: prod(nband), prodr(nband), prodx(nband,nband)

!-----declare local variables
integer :: ib, jb, ibb, jbb, is, iss
integer :: ismin, ismax, icmin, icmax, num, i, j, l
real*8  :: umatmax, wrkary(100), wrkaryI
integer :: iwrkary(100), isort(100)
logical :: lexchange
real*8  :: amx
integer :: ic, icx, ib1, ib2, jbb2, jbx, ihband, ipband


if( .not.ltddft_fssh ) return

!if( lrtddft ) then
!    call lrtddft_get_ihpband( ihband, ipband )
!else
    ihband = 1
    ipband = noband
!end if

lsetmatrix = .true.

!---Search the lowest-unoccupied and highest-occupied states
do is = 1, noccmx
   if( .not.lspin .and. occ_rec(is,nspin) < 1.5d0 .or. &
&           lspin .and. occ_rec(is,nspin) < 0.5d0 ) then
       ismin = is
       exit
   end if
end do
do is = noccmx, 1, -1
   if( occ_rec(is,nspin) > 0.5d0 ) then
       ismax = is
       exit
   end if
end do
if( .not.lrtddft ) then
    ihband = max( ismin - 2, 1 )
    ipband = min( ismax + 2, noccmx )
end if

!!---check matrix
!if( myid == 0 ) then
!    ib1 = 290
!    ib2 = 330
!    write(*,'(3x,100i4)') (jb,jb=ib1, ib2)
!    do ib = ib1, ib2
!    write(*,'(i3,100f4.1)') ib, (abs(aijr(ib,jb)),jb=ib1, ib2)
!    end do
!end if

!---Search the largest component in each row
!---Note: aijr(ib,jb) = <phi_ib(t)|S|phi_jb(t+dt)>
prod(1:nband)  = 0.d0
prodr(1:nband) = 0.d0
do ib = 1, noband
   umatmax = 0.d0
   do jb = 1, noband
      if( umatmax < abs(aijr(ib,jb)) ) then
          umatmax = abs(aijr(ib,jb))
          jbb = jb
      end if
   end do
   prod(ib) = jbb
   prodr(jbb) = prodr(jbb) + 1.d0
   prodx(nint(prodr(jbb)),jbb) = ib
end do
!!---check
!if( myid == 0 ) then
!    do ib = ihband, ipband
!        write(*,*) ib, prod(ib), prodr(ib)
!    end do
!end if

!---guarantee one-by-one correspondence
do jbb = 1, noband
   if( nint(prodr(jbb)) > 1 ) then
       amx = 0.d0
       do ic = 1, nint(prodr(jbb))
          ib1 = nint(prodx(ic,jbb))
          if( abs(aijr(ib1,jbb)) > amx ) then
              amx = abs(aijr(ib1,jbb))
              ib2 = ib1
          end if
       end do
       !---error trap
       if( nint(prod(ib2)) /= jbb ) call fstop( nfile, myid, nodes, &
     &                  'error in tddft_fssh_set_matrix2' )
       icx = nint(prodr(jbb))
       do ic = 1, icx
          ib1 = nint(prodx(ic,jbb))
          if( ib1 == ib2 ) cycle
          amx = 0.d0
          do jbb2 = 1, noband
             if( nint(prodr(jbb2)) == 0 ) then
             if( abs(aijr(ib1,jbb2)) > amx ) then
                 amx = abs(aijr(ib1,jbb2))
                 jbx = jbb2
             end if
             end if
          end do
          prodr(jbx) = 1.d0
          prodr(jbb) = prodr(jbb) - 1.d0
          prod(ib1) = jbx
       end do
   end if
end do
do ib = 1, noband
   jbb = prod(ib)
   prodx(1,jbb) = ib
end do

!---avoid excitation exceeding the limit of LRTDDFT
do ib = 1, ihband - 1
   prod(ib) = ib
   prodr(ib) = 1.d0
end do
do ib = ipband+1, noband
   prod(ib) = ib
   prodr(ib) = 1.d0
end do
do ib = ihband, ipband
   if( nint(prod(ib)) < ihband .or. nint(prod(ib)) > ipband ) then
       ib1 = ib
       prod(ib1) = ib1
       do
          ib1 = prodx(1,ib1)
          if( ib1 == ib ) exit
          prod(ib1) = ib1
       end do
   end if
end do
!!---check
!if( myid == 0 ) then
!    do ib = ihband, ipband
!        write(*,*) ib, prod(ib), prodr(ib)
!    end do
!end if

!---Check consistency
lsetmatrix = .true.
jbb = 0
do jb = 1, noband
   num = nint(prodr(jb))
   if( num == 0 ) cycle
   i = 0
   do ib = 1, noband
      if( jb == nint(prod(ib)) ) then
          i = i + 1
          iwrkary(i) = ib
           wrkary(i) = abs(aijr(ib,jb))
      end if
   end do
   !---error trap
   if( i /= num ) then
       lsetmatrix = .false.
       exit
   end if
   if( num == 1 ) then
       jbb = jbb + 1
       bmthi(iwrkary(1),1) = jbb
   else
       !---sorting
       isort(1) = 1
       ido: do i = 2, num
          jdo: do j = 1, i - 1
            IF( wrkary(i) < wrkary(j) ) cycle jdo
            wrkaryi = wrkary(i)
            do l = i, j + 1, -1
               wrkary(l) = wrkary(l-1)
               isort(l) = isort(l-1)
            end do
            wrkary(j) = wrkaryi
            isort(j) = i
            cycle ido
          end do jdo
          isort(i) = i
       end do ido
       do i = 1, num
          jbb = jbb + 1
          bmthi(iwrkary(isort(i)),1) = jbb
       end do
    end if
end do
prod(1:noband) = bmthi(1:noband,1)
!---check
prodr(1:nband) = 0.d0
do ib = 1, noband
   jbb = nint(prod(ib))
   prodr(jbb) = prodr(jbb) + 1.d0
end do
do jb = 1, noband
   num = nint(prodr(jb))
   if( num /= 1 ) then
       lsetmatrix = .false.
       exit
   end if
end do

icmin = noband + 1
icmax = 0
do ib = 1, noband
   if( ib /= nint(prod(ib)) ) then
       icmin = ib
       exit
   end if
end do
do ib = noband, 1, -1
   if( ib /= nint(prod(ib)) ) then
       icmax = ib
       exit
   end if
end do
!do jbb = 1, noband
!   if( 0.5d0 < prodr(jbb) .and. prodr(jbb) < 1.5d0 ) cycle
!   if( jbb < ismin ) then
!       icmin = max( icmin, jbb )
!   else if( jbb > ismax ) then
!       icmax = min( icmax, jbb )
!   else
!       lsetmatrix = .false.
!   end if
!end do
!do ib = 1, noband
!   jbb = nint(prod(ib))
!   if( 0.5d0 < prodr(jbb) .and. prodr(jbb) < 1.5d0 ) cycle
!   if( ib < ismin ) then
!       icmin = max( icmin, ib )
!   else if( ib > ismax ) then
!       icmax = min( icmax, ib )
!   else
!       lsetmatrix = .false.
!   end if
!end do


if( lsetmatrix ) then
    lexchange = .false.
    do ib = 1, noband
       prodr(ib) = ib
    end do
    occ_store(1:nband,nspin) = occ_rec(1:nband,nspin)
    if( .not.lrtddft ) then
        if( ldensity_matrix_repre ) then
               kr(-nupdn:nupdn,-nupdn:nupdn,1:noccmx)  &
&         = dnmxr(-nupdn:nupdn,-nupdn:nupdn,1:noccmx,nspin)
               ki(-nupdn:nupdn,-nupdn:nupdn,1:noccmx)  &
&         = dnmxi(-nupdn:nupdn,-nupdn:nupdn,1:noccmx,nspin)
        else
            amtr(1:noband,1:noccmx)  = cr(1:noband,1:noccmx,nspin)
            amti(1:noband,1:noccmx)  = ci(1:noband,1:noccmx,nspin)
        end if
        bmthr(1:noccmx,1:noccmx) = accumprob(1:noccmx,1:noccmx,nspin)
    end if
    do ib = 1, noccmx
       if( lsignc(ib,nspin) ) then
           bmthi(ib,1) = -1.d0
       else
           bmthi(ib,1) =  1.d0
       end if
    end do

    do ib = 1, icmin - 1
       phisphi(ib,1:noband,nspin) = aijr(ib,1:noband)
    end do
    if( icmax /= 0 ) then
        do ib = icmax + 1, noband
           phisphi(ib,1:noband,nspin) = aijr(ib,1:noband)
        end do
    end if
    do ib = icmin, icmax
       if( prod(ib) < 0.5d0 ) cycle
       ibb = ib
       do
          jbb = nint(prod(ibb))
          phisphi(jbb,1:noband,nspin) = aijr(ibb,1:noband)
          if( jbb /= ibb ) then
              !---Exchange occupations and wave packets
              lexchange = .true.
              !if( abs(occ_store(jbb,nspin)-occ_store(ibb,nspin)) > 0.5d0 ) then
              if(loutfile(1)) write(nfile(1),'(a,i4,a,i4,a)') &
& '*** exchange wave packets bet. the', ibb, 'th and', jbb, 'th states (not transition)'
              if(loutfile(2)) write(nfile(2),'(a,i4,a,i4,a)') &
& '*** exchange wave packets bet. the', ibb, 'th and', jbb, 'th states (not transition)'
              !end if
              prodr(ibb) = jbb
              occ_rec(jbb,nspin) = occ_store(ibb,nspin)
              if( .not.lrtddft ) then
                  if( ldensity_matrix_repre ) then
                      dnmxr(-nupdn:nupdn,-nupdn:nupdn,jbb,nspin)  &
&                      = kr(-nupdn:nupdn,-nupdn:nupdn,ibb)
                      dnmxi(-nupdn:nupdn,-nupdn:nupdn,jbb,nspin)  &
&                      = ki(-nupdn:nupdn,-nupdn:nupdn,ibb)
                  else
                      cr(1:noband,jbb,nspin) = amtr(1:noband,ibb)
                      ci(1:noband,jbb,nspin) = amti(1:noband,ibb)
                  end if
              end if
          end if
          prod(ibb) = 0
          if( ib == jbb ) exit
          ibb = jbb
       end do
    end do


    if( lexchange ) then

        if( lrtddft ) then
            do jb = 1, nexciton
               if( .not.lspin ) then
                   do ib = 1, noband
                      if( ib == iband_hole(jb) ) then
                          ibb = nint(prodr(ib))
                          if( ibb <= nbase .or. ibb > nbase .and. iband_electron(jb) == ibb ) then
                              iband_hole(jb) = ibb
                          end if
                          exit
                      end if
                   end do
                   do ib = 1, noband
                      if( ib == iband_electron(jb) ) then
                          ibb = nint(prodr(ib))
                          if( ibb > nbase .or. ibb <= nbase .and. iband_hole(jb) > nbase ) then
                              iband_electron(jb) = ibb
                          end if
                          exit
                      end if
                   end do
                   if( iband_hole(jb) > iband_electron(jb) ) then
                       !---this is ground state
                           iband_hole(jb) = 0
                       iband_electron(jb) = 0
                   end if
               else
                   if( ldegenerate(jb) ) then
                       !---triplet
                       if( nspin == 2 ) then
                           do ib = 1, noband
                              if( ib == iband_hole(jb) ) then
                                  ibb = nint(prodr(ib))
                                  if( ibb <= nbase ) then
                                      iband_hole(jb) = ibb
                                  end if
                                  exit
                              end if
                           end do
                       else
                           do ib = 1, noband
                              if( ib == iband_electron(jb) ) then
                                  ibb = nint(prodr(ib))
                                  if( ibb > nbase ) then
                                      iband_electron(jb) = ibb
                                  end if
                                  exit
                              end if
                           end do
                       end if
!                       if( iband_hole(jb) > iband_electron(jb) ) then
!                           !---this is ground state
!                               iband_hole(jb) = 0
!                           iband_electron(jb) = 0
!                       end if
                   else
                       !---singlet
                       if( nspin == 1 ) then
                           do ib = 1, noband
                              if( ib == iband_hole(jb) ) then
                                  ibb = nint(prodr(ib))
                                  if( ibb <= nbase .or. ibb > nbase .and. iband_electron(jb) == ibb ) then
                                      iband_hole(jb) = ibb
                                  end if
                                  exit
                              end if
                           end do
                           do ib = 1, noband
                              if( ib == iband_electron(jb) ) then
                                  ibb = nint(prodr(ib))
                                  if( ibb > nbase .or. ibb <= nbase .and. iband_hole(jb) > nbase ) then
                                      iband_electron(jb) = ibb
                                  end if
                                  exit
                              end if
                           end do
                       else
                           if( iband_hole(jb) > iband_electron(jb) ) then
                               if( nbase+1 == nint(prodr(nbase)) .and. nbase == nint(prodr(nbase+1)) ) then
                                       iband_hole(jb) = nbase
                                   iband_electron(jb) = nbase + 1
                               end if
                           end if
                           if( iband_hole(jb) > iband_electron(jb) ) then
                               !---this is ground state
                                   iband_hole(jb) = 0
                               iband_electron(jb) = 0
                           end if
                       end if
                   end if
               end if
            end do
            do jb = 0, nlrreduce
               do ib = 1, noband
                  if( ib == ibstate(1,jb) ) then
                      ibb = nint(prodr(ib))
                      ibstate(1,jb) = ibb
                      exit
                  end if
               end do
               do ib = 1, noband
                  if( ib == ibstate(2,jb) ) then
                      ibb = nint(prodr(ib))
                      ibstate(2,jb) = ibb
                      exit
                  end if
               end do
            end do
            if( lset_ibstate_rec ) then
            do iss = 1, nexciton
            do jb = 0, nlrreduce
               do ib = 1, noband
                  if( ib == ibstate_rec(1,jb,iss) ) then
                      ibb = nint(prodr(ib))
                      ibstate_rec(1,jb,iss) = ibb
                      exit
                  end if
               end do
               do ib = 1, noband
                  if( ib == ibstate_rec(2,jb,iss) ) then
                      ibb = nint(prodr(ib))
                      ibstate_rec(2,jb,iss) = ibb
                      exit
                  end if
               end do
            end do
            end do
            end if
            !---check the band index of hole and electron
!            call lrtddft_chkorder( nfile, myid, nodes, prodr, noband, nspin )

        else

            do ib = 1, noband
               ibb = nint(prodr(ib))
            do jb = 1, noband
               jbb = nint(prodr(jb))
               accumprob(ibb,jbb,nspin) = bmthr(ib,jb)
            end do
            end do
        end if

            do ib = 1, noband
               ibb = nint(prodr(ib))
               lsignc(ibb,nspin) = bmthi(ib,1) < 0.d0
            end do

    end if
    !---Consider sign change
    do ib = 1, noband
        bmthi(ib,1) = phisphi(ib,ib,nspin)
    end do
    do ib = 1, noband
       if( lsignc(ib,nspin) ) then
           phisphi(ib,1:noband,nspin) = - phisphi(ib,1:noband,nspin)
       end if
!       if( phisphi(ib,ib,nspin) < 0.d0 ) lsignc(ib,nspin) = .not.lsignc(ib,nspin)
       if( bmthi(ib,1) < 0.d0 ) lsignc(ib,nspin) = .not.lsignc(ib,nspin)
       if( lsignc(ib,nspin) ) then
           phisphi(1:noband,ib,nspin) = - phisphi(1:noband,ib,nspin)
       end if
    end do

  else

    do i = 1, 2
    if( loutfile(i) ) then
        write(nfile(i),*) '*** TDDFT-FSSH: error in tddft_fssh_set_matrix'
        write(nfile(i),*) '*** stop integration of TDDFT equations'
        do ib = ismin, ismax
!        do ib = 1, nband
            write(nfile(i),'(i5,100f8.4)') ib, aijr(ib,ismin:ismax)
        end do
    end if
    end do
!       call fstop( nfile, myid, nodes, &
!&                  'TDDFT-FSSH: 1-by-1 correspondence error' )
end if

!---check
!if( myid == 0 ) then
!    do ib = 1, noband
!       write(*,'(100f8.4)') ( phisphi(ib,jb,nspin), jb = 1, noband )
!    end do
!end if


return
end subroutine




subroutine tddft_fssh( nfile, myid, nodes,  &
& nband, nbnod, nspnmx, eig, norder, occ, dtmd, &
& dkj, cmr, cmi, cdotr, cdoti, egv, w1r, w1i, w2r, bm1r, bm1i, nstep_md, &
& keypwv, keypcd, laspc_exec, ltransition, lexcitedstates )
!-----------------------------------------------------------------------
! TDDFT-FSSH main routine
!-----------------------------------------------------------------------
use tddft_variables
use tddft_fssh_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: nband, nbnod, nspnmx
real*8,  dimension(nband,nspnmx) :: eig, occ
integer, dimension(nband,nspnmx) :: norder
real*8  :: dtmd
real*8,  dimension(nband,nband) :: dkj, cmr, cmi
real*8,  dimension(nband) :: cdotr, cdoti, egv, w1r, w1i, w2r, bm1r, bm1i
integer :: nstep_md
integer :: keypwv, keypcd
logical :: laspc_exec
logical :: ltransition
logical :: lexcitedstates


!if( lrtddft ) then
!
!    !---exciton dynamics
!    call tddft_fssh_exciton( nfile, myid, nodes,  &
!& nband, nbnod, nspnmx, eig, norder, occ, dtmd, &
!& dkj, cmr, cmi, cdotr, cdoti, egv, w1r, w1i, w2r, bm1r, bm1i, nstep_md, &
!& keypwv, keypcd, laspc_exec, ltransition, lexcitedstates )
!
!else

    lexcitedstates = .true.
    if( ldensity_matrix_repre ) then
        !--- density matrix representation
        call tddft_fssh_dm( nfile, myid, nodes,  &
& nband, nbnod, nspnmx, eig, norder, occ, dtmd, &
& dkj, cmr, cmi, cdotr, cdoti, egv, w1r, w1i, w2r, bm1r, bm1i, nstep_md, &
& keypwv, keypcd, laspc_exec, ltransition )
    else
!        call tddft_fssh_coef( nfile, myid, nodes,  &
!& nband, nbnod, nspnmx, eig, norder, occ, dtmd, &
!& dkj, cmr, cmi, cdotr, cdoti, egv, w1r, w1i, w2r, bm1r, bm1i, nstep_md, &
!& keypwv, keypcd, laspc_exec, ltransition )
    end if

!end if


return
end subroutine




subroutine tddft_fssh_dm( nfile, myid, nodes,  &
& nband, nbnod, nspnmx, eig, norder, occ, dtmd, &
& dkj, cmr, cmi, cdotr, cdoti, egv, w1r, w1i, w2r, bm1r, bm1i, nstep_md, &
& keypwv, keypcd, laspc_exec, ltransition )
!-----------------------------------------------------------------------
! TDDFT-FSSH main routine
!-----------------------------------------------------------------------
use outfile
use tddft_variables
use tddft_fssh_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
integer :: nband, nbnod, nspnmx
real*8,  dimension(nband,nspnmx) :: eig, occ
integer, dimension(nband,nspnmx) :: norder
real*8  :: dtmd
real*8,  dimension(nband,nband) :: dkj, cmr, cmi
real*8,  dimension(nband) :: cdotr, cdoti, egv, w1r, w1i, w2r, bm1r, bm1i
integer :: nstep_md
integer :: keypwv, keypcd
logical :: laspc_exec
logical :: ltransition

!-----declare local variables
integer :: nstep, nspin, jb, kb, is, i, iud, news
integer :: nupdn1, nupdn2, jbb, kbb, lbb, lb, irk
real*8  :: rn, probsum, anr, acurr, boltzmn, bkl, bkl0
logical :: ldmcalc
integer :: digit, ierror
character(50) :: fname
character(1) :: cud(2)
real*8  :: decotime, rnexp, rn2, dlambda
real*8  :: ct, ct00, ct0, timecnt
logical :: ltimecnt = .false.
save ltimecnt


if( .not.ltddft_fssh ) return

!---store the previous total energies
if( lfssh_vscale ) then
if( lsetekin ) then
    if( npreve == nprevex ) then
        do i = 1, npreve - 1
           totene(i) = totene(i+1)
        end do
        npreve = npreve - 1
    end if
    npreve = npreve + 1
    totene(npreve) = theekin + thepote
  else
    return
end if
end if

if( .not.lsetmatrix ) then
    !---Initialization of coefficients
    dnmxr = 0.d0
    dnmxi = 0.d0
    dnmxr(0,0,1:nband,1:nspnmx) = 1.d0
    accumprob = 0.d0
    if( lfssh_dish ) then
        dish_prev_step(:,:) = -1
    end if
    lsignc = .false.
    return
end if

if( lreject ) then
    lreject = .false.
    if( .not.lfssh_gsscf ) then
        keypwv = 0
        keypcd = 0
        laspc_exec = .false.
    end if
    return
end if

if( lfssh_dish ) then
    do nspin = 1, nspnmx
    do jb = 1, noband
       if( dish_prev_step(jb,nspin) == -1 ) dish_prev_step(jb,nspin) = nstep_md - 1
    end do
    end do
end if


!---generate random number
if( myid == 0 ) then
    call rnd00( rn, rndx )
end if
call dbcast(rn,1,0)


ct = timecnt()
ct00 = ct
ct0  = ct

occ_store = occ_rec

probsum = 0.d0
cud(1) = 'u'
cud(2) = 'd'
nstep = dtmd/dttddft
spin: do nspin = 1, nspnmx

   do jb = 1, noband
   do kb = 1, noband
      dkj(kb,jb) = -(phisphi(kb,jb,nspin)-phisphi(jb,kb,nspin))/2.d0/dtmd
   end do
   end do

   probinteg(-iudx:iudx,1:noccmx) = 0.d0

   !---Integration
   integdo: do i = 1, nstep

      !---Runge-Kutta integration
      isdo: do is = 1, noccmx
         ldmcalc = .false.
         if( occ_store(is,nspin) > 0.5d0 ) then
!             do iud = -1, 1, 2
             do iud = -iudx, iudx
                news = is + iud
                if( news <= 0 .or. news > noband .or. news == is ) cycle
                if( nspnmx == 1 .and. occ_store(news,nspin) > 1.5d0 .or. &
                  & nspnmx == 2 .and. occ_store(news,nspin) > 0.5d0 ) cycle
                    ldmcalc = .true.
             end do
         end if
         if( .not.ldmcalc ) cycle isdo

         nupdn1 = max( -nupdn, 1 - is )
         nupdn2 = min(  nupdn, noband - is )

         do kbb = nupdn1, nupdn2
            kb = is + kbb
         do jbb = nupdn1, nupdn2
            jb = is + jbb
            kr(jbb,kbb,1) =  ( eig(jb,nspin)-eig(kb,nspin) )  &
&                           *dnmxi(jbb,kbb,is,nspin)
            ki(jbb,kbb,1) = -( eig(jb,nspin)-eig(kb,nspin) )  &
&                           *dnmxr(jbb,kbb,is,nspin)
            do lbb = nupdn1, nupdn2
               lb = is + lbb
               kr(jbb,kbb,1) = kr(jbb,kbb,1)  &
&    - dnmxr(jbb,lbb,is,nspin)*dkj(lb,kb) + dnmxr(lbb,kbb,is,nspin)*dkj(jb,lb)
               ki(jbb,kbb,1) = ki(jbb,kbb,1)  &
&    - dnmxi(jbb,lbb,is,nspin)*dkj(lb,kb) + dnmxi(lbb,kbb,is,nspin)*dkj(jb,lb)
            end do
         end do
         end do

         do irk = 2, 4
            do kbb = nupdn1, nupdn2
            do jbb = nupdn1, nupdn2
               dmtr(jbb,kbb) = dnmxr(jbb,kbb,is,nspin) + rkc(irk)*kr(jbb,kbb,irk-1)
               dmti(jbb,kbb) = dnmxi(jbb,kbb,is,nspin) + rkc(irk)*ki(jbb,kbb,irk-1)
            end do
            end do

            do kbb = nupdn1, nupdn2
               kb = is + kbb
            do jbb = nupdn1, nupdn2
               jb = is + jbb
               kr(jbb,kbb,irk) =  ( eig(jb,nspin)-eig(kb,nspin) )*dmti(jbb,kbb)
               ki(jbb,kbb,irk) = -( eig(jb,nspin)-eig(kb,nspin) )*dmtr(jbb,kbb)
               do lbb = nupdn1, nupdn2
                  lb = is + lbb
                  kr(jbb,kbb,irk) = kr(jbb,kbb,irk)  &
&    - dmtr(jbb,lbb)*dkj(lb,kb) + dmtr(lbb,kbb)*dkj(jb,lb)
                  ki(jbb,kbb,irk) = ki(jbb,kbb,irk)  &
&    - dmti(jbb,lbb)*dkj(lb,kb) + dmti(lbb,kbb)*dkj(jb,lb)
               end do
            end do
            end do
         end do

         do kbb = nupdn1, nupdn2
         do jbb = nupdn1, nupdn2
            dnmxr(jbb,kbb,is,nspin) = dnmxr(jbb,kbb,is,nspin)  &
&            + (      kr(jbb,kbb,1) + 2.d0*kr(jbb,kbb,2)  &
&               + 2d0*kr(jbb,kbb,3) +      kr(jbb,kbb,4) )*dttddft/6.d0
            dnmxi(jbb,kbb,is,nspin) = dnmxi(jbb,kbb,is,nspin)  &
&            + (      ki(jbb,kbb,1) + 2.d0*ki(jbb,kbb,2)  &
&               + 2d0*ki(jbb,kbb,3) +      ki(jbb,kbb,4) )*dttddft/6.d0
         end do
         end do

         !---Integration of transition probability
!         do iud = -1, 1, 2
!         do iud = -2, 2
         do iud = -iudx, iudx
            news = is + iud
            if( news <= 0 .or. news > noband .or. news == is ) cycle
            anr   = dnmxr(iud,0,is,nspin)
            acurr = dnmxr(0,0,is,nspin)
            probinteg(iud,is) = probinteg(iud,is) + anr/acurr
         end do

      end do isdo

   end do integdo

   if( ltimecnt ) then
       ct = timecnt()
       if(loutfile(1)) write(nfile(1),*) ' Integration of DM : cpu-time :', ct-ct0
       if(loutfile(2)) write(nfile(2),*) ' Integration of DM : cpu-time :', ct-ct0
       ct0 = ct
   end if

!   !---orthonormalization check
!   do is = 1, noccmx
!   do js = 1, noccmx
!      tr = 0.d0
!      do kb = 1, noband
!         tr = tr + cr(kb,is,nspin)*cr(kb,js,nspin) + ci(kb,is,nspin)*ci(kb,js,nspin)
!      end do
!      write(*,*) is, js, tr
!   end do
!   end do

!   if( nspin == 2 ) stop

   ltransition = .false.
   !---Probability
   isdo2: do is = 1, noccmx

      if( lfssh_dish ) then
          !---Decoherence time
          dlambda = 0.d0
          do iud = -iudx, iudx
             news = is + iud
             if( news <= 0 .or. news > noband .or. news == is ) cycle
             dlambda = dlambda + decoherence_rate(is,news)  &
&                     * ( dnmxr(iud,iud,is,nspin) + dnmxi(iud,iud,is,nspin) )
          end do
          if( dlambda > 1.d-30 ) then
              decotime = 1.d0/max(dlambda, 1.d-30)
              if( myid == 0 ) then
              if( decotime/dtmd < 1.d+09 ) then
                  write(nfile(1),'(a,i4,a,2x,a,f10.0,i8)') &
& 'State:', is, cud(nspin), 'Decoherence/Elapsed time (in MD step) =', decotime/dtmd,  &
& nstep_md - dish_prev_step(is,nspin)
                  write(nfile(2),'(a,i4,a,2x,a,f10.0,i8)') &
& 'State:', is, cud(nspin), 'Decoherence/Elapsed time (in MD step) =', decotime/dtmd,  &
& nstep_md - dish_prev_step(is,nspin)
              else
                  write(nfile(1),'(a,i4,a,2x,a,es14.6,i8)') &
& 'State:', is, cud(nspin), 'Decoherence/Elapsed time (in MD step) =', decotime/dtmd,  &
& nstep_md - dish_prev_step(is,nspin)
                  write(nfile(2),'(a,i4,a,2x,a,es14.6,i8)') &
& 'State:', is, cud(nspin), 'Decoherence/Elapsed time (in MD step) =', decotime/dtmd,  &
& nstep_md - dish_prev_step(is,nspin)
              end if
              end if

              if( myid == 0 ) then
                  call rndexp( rnexp, dlambda )
              end if
              call dbcast(rnexp,1,0)

              if( myid == 0 ) then
                  write(nfile(1),'(17x,a,f19.0)')  &
& 'vs. random number (Poisson dist.) =', rnexp/dtmd
                  write(nfile(2),'(17x,a,f19.0)')  &
& 'vs. random number (Poisson dist.) =', rnexp/dtmd
              end if
              if( rnexp < (nstep_md - dish_prev_step(is,nspin))*dtmd ) then
                  !---the state vector is reduced.
                  if( myid == 0 ) then
                      call rnd00( rn2, rndx )
                  end if
                  call dbcast(rn2,1,0)
                  if( rn2 < dnmxr(0,0,is,nspin) + dnmxi(0,0,is,nspin) ) then
                      if( myid == 0 ) then
                          write(nfile(1),'(a)') &
& '---> DISH!!! w.f. is projected onto the decohering state.'
                          write(nfile(2),'(a)') &
& '---> DISH!!! w.f. is projected onto the decohering state.'
                      end if
                      dnmxr(:,:,is,nspin) = 0.d0
                      dnmxi(:,:,is,nspin) = 0.d0
                      dnmxr(0,0,is,nspin) = 1.d0
                      dish_prev_step(is,nspin) = nstep_md
                      cycle isdo2
                  end if
              end if
          end if
      end if

      if( occ_store(is,nspin) > 0.5d0 ) then
!          do iud = -1, 1, 2
!          do iud = -2, 2
          do iud = -iudx, iudx
             news = is + iud
             if( news <= 0 .or. news > noband .or. news == is ) cycle
             if( nspnmx == 1 .and. occ_store(news,nspin) > 1.5d0 .or. &
               & nspnmx == 2 .and. occ_store(news,nspin) > 0.5d0 ) cycle
!             if( nspnmx == 1 .and. iud == 1 .and. occ_store(is,nspin) < 1.5d0 &
!               & .and. occ_store(news,nspin) > 0.5d0 .and. occ_store(news,nspin) < 1.5d0 ) cycle
!             anr = cr(is,is,nspin)*cr(news,is,nspin) + ci(is,is,nspin)*ci(news,is,nspin)
!             acurr = cr(is,is,nspin)**2 + ci(is,is,nspin)**2
!             bkl = 2.d0*dkj(news,is)*dtmd*anr/acurr
             if( lfssh_boltzmn .and. iud == 1 ) then
                 !---Boltzmann factor for upward transitions
                 boltzmn = exp( -(eig(news,nspin)-eig(is,nspin))/treq)
               else
                 boltzmn = 1.d0
             end if
             bkl0 = 2.d0*dkj(news,is)*dttddft*probinteg(iud,is)
             bkl  = bkl0*boltzmn

             if( bkl > 0.d0 ) &
&               accumprob(news,is,nspin) = accumprob(news,is,nspin) + bkl

             !---output probability
             if( loutfile(1) ) then
                 !---if file is not opened, open it.
                 if( iunit_prob(iud,is,nspin) == 0 ) then
                     call allocate_unit_number( iunit_prob(iud,is,nspin) )
                     digit = log(dble(noband)+0.5d0)/log(10.d0) + 1.d0
                     fname = 'data/qm_fsshprob_'
                     call get_dumpfname( fname, is, digit )
                     fname = ( fname(1:len_trim(fname)) // 'to' )
                     call get_dumpfname( fname, news, digit )
                     if( nspin == 1 ) then
                         fname = ( fname(1:len_trim(fname)) // '-u.d' )
                       else
                         fname = ( fname(1:len_trim(fname)) // '-d.d' )
                     end if
                     open( iunit_prob(iud,is,nspin), file=fname, status='unknown',  &
     &                     action='write', iostat=ierror )
                     write(iunit_prob(iud,is,nspin),'(a)') "# step  probability  accumulation"
                 end if

                 write(iunit_prob(iud,is,nspin),'(i6,2es13.5)') &
&                             nstep_md, bkl, accumprob(news,is,nspin)

                 if( boltzmn > 0.99999d0 ) then
                     write(nfile(1),'(a,f9.6,11x,a,f9.6,a,i4,a,a,i4,a)') &
& '*** Probability=', bkl, '(', accumprob(news,is,nspin), &
&' ) from state', is, cud(nspin), ' to', news, cud(nspin)
                     write(nfile(2),'(a,f9.6,11x,a,f9.6,a,i4,a,a,i4,a)') &
& '*** Probability=', bkl, '(', accumprob(news,is,nspin), &
&' ) from state', is, cud(nspin), ' to', news, cud(nspin)
                   else
                     write(nfile(1),'(a,f9.6,a,f9.6,a,f9.6,a,i4,a,a,i4,a)') &
& '*** Probability=', bkl0, '->', bkl, '(', accumprob(news,is,nspin), &
&' ) from state', is, cud(nspin), ' to', news, cud(nspin)
                     write(nfile(2),'(a,f9.6,a,f9.6,a,f9.6,a,i4,a,a,i4,a)') &
& '*** Probability=', bkl0, '->', bkl, '(', accumprob(news,is,nspin), &
&' ) from state', is, cud(nspin), ' to', news, cud(nspin)
                 end if
             end if


             !---compare probsum with random number
             if( lfssh_switch ) then
             if( probsum <= rn .and. rn < probsum + max( bkl, 0.d0 ) ) then
                 ltransition = .true.
                 if( lfssh_vscale ) then
                     if(loutfile(1)) write(nfile(1),'(a)') 'Electronic transition is temporarily accepted.'
                     if(loutfile(2)) write(nfile(2),'(a)') 'Electronic transition is temporarily accepted.'
                 else
                     if(loutfile(1)) write(nfile(1),'(a)') 'Electronic transition is accepted.'
                     if(loutfile(2)) write(nfile(2),'(a)') 'Electronic transition is accepted.'
                 end if
                 occ_rec(is,nspin)   = occ_rec(is,nspin)   - 1.d0
                 occ_rec(news,nspin) = occ_rec(news,nspin) + 1.d0
             end if
             end if
             probsum = probsum + max( bkl, 0.d0 )

          end do
      end if
   end do isdo2

end do spin
if( .not.lfssh_gsscf ) then
    occ = occ_rec
end if
if( ltransition ) then
    !---Initialization of coefficients
!    cr = 0.d0
!    ci = 0.d0
!    do is = 1, noccmx
!       cr(is,is,1:nspnmx) = 1.d0
!    end do
!    accumprob = 0.d0
!    lsignc = .false.
    if( .not.lfssh_gsscf ) then
        keypwv = 0
        keypcd = 0
        laspc_exec = .false.
    end if
end if

ct = timecnt()
if( ltimecnt ) then
    if(loutfile(1)) write(nfile(1),*) 'tddft_fssh-2:', ct - ct0
    if(loutfile(2)) write(nfile(2),*) 'tddft_fssh-2:', ct - ct0
end if
if(loutfile(1)) write(nfile(1),*) 'tddft_fssh : cpu-time (sec)  :', ct - ct00
if(loutfile(2)) write(nfile(2),*) 'tddft_fssh : cpu-time (sec)  :', ct - ct00


return
end subroutine




subroutine tddft_fssh_judge_accept( nfile, myid, nodes,  &
& ltransition, nband, nspnmx, occ, sume, &
& vatm, frc, watom, ntype, natom, nhk1, nhk2, dtmd )
!-----------------------------------------------------------------------
! TDDFT-FSSH: check acceptance and scale velocities if necessary
!-----------------------------------------------------------------------
use outfile
use constants
use tddft_variables
use tddft_fssh_variables
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
logical :: ltransition
integer :: nband, nspnmx
real*8,  dimension(nband,nspnmx) :: occ
real*8  :: sume
integer :: ntype, natom
real*8,  dimension(3,natom) :: vatm, frc
real*8,  dimension(ntype)   :: watom
integer, dimension(ntype)   :: nhk1, nhk2
real*8  :: dtmd

!-----declare local variables
integer :: is, icase, it, i
real*8  :: eest, de1, de2, ex, en, a, b, c
real*8  :: curekin, targetekin


!---scale velocity if necessary
if( lfssh_vscale ) then
    if( npreve == 1 ) then
        eest = totene(1)
    else if( npreve == 2 ) then
        eest = totene(2)      ! + totene(2) - totene(1)
    else if( npreve >= 3 ) then
        de1 = totene(npreve)   - totene(npreve-1)
        de2 = totene(npreve-1) - totene(npreve-2)
        de1 = de1 / natom * 13.6d04
        de2 = de2 / natom * 13.6d04
        ex  = max( abs(de1), abs(de2) )
        en  = min( abs(de1), abs(de2) )
        if(          en >= 100.d0 .and. abs(de1-de2)/en < 0.2d0 .or. &
& en >=  20.d0 .and. en <  100.d0 .and. abs(de1-de2)/en < 0.7d0-en/200d0 .or. &
& en <   20.d0 .and. ex >=  20.d0 .and. abs(de1-de2)/en < 1.d0 .or. &
& en <   20.d0 .and. ex <   20.d0 ) then
            a = 0.5d0*( totene(npreve) - 2.d0*totene(npreve-1) + totene(npreve-2) )
            b = 0.5d0*( totene(npreve) - totene(npreve-2) )
            c = totene(npreve-1)
            eest = 4.d0*a + 2.d0*b + c
            if(loutfile(1)) write(nfile(1),'(a,3es15.7)') ' *** Prev. total E.:', &
& totene(npreve-2), totene(npreve-1), totene(npreve)
            if(loutfile(2)) write(nfile(2),'(a,3es15.7)') ' *** Prev. total E.:', &
& totene(npreve-2), totene(npreve-1), totene(npreve)
          else
            eest = totene(npreve)
        end if
    end if
    targetekin = eest - sume    ! Target kinetic energy
    ltransition = targetekin/natom > tminimum
    do i = 1, 2
    if( loutfile(i) ) then
        write(nfile(i),*) '*** Estimated total E. :', eest
        write(nfile(i),*) '*** The current pot. E.:', sume
        write(nfile(i),*) '*** Target & minimum temperature [K]:', &
& targetekin/natom* (tempau/2.d0) * 2.d0/3.d0, tminimum* (tempau/2.d0) * 2.d0/3.d0
    end if
    end do
    if( ltransition ) then
        !---calculate the current kinetic energy
        curekin = 0.d0
        do it = 1, ntype
           do i = nhk1(it), nhk2(it)
              curekin = curekin + watom(it)*( &
&             vatm(1,i)*vatm(1,i) + vatm(2,i)*vatm(2,i) + vatm(3,i)*vatm(3,i) )
           end do
        end do
        curekin = 0.5d0 * curekin

        !---set non-adiabatic coupling vector <- currently is approximated by velocity
        do it = 1, ntype
           do i = nhk1(it), nhk2(it)
              couplvec(1:3,i) = watom(it)*vatm(1:3,i)
           end do
        end do

        !---velocity scaling
        call tddft_fssh_vscale( nfile, myid, nodes,  &
& ltransition, couplvec, curekin, targetekin, vatm, frc, watom, &
& ntype, natom, nhk1, nhk2, dtmd )
    end if
end if
lreject = .not.ltransition
if( lreject ) then
    if(loutfile(1)) write(nfile(1),*) '*** The electronic transition has been rejected. ***'
    if(loutfile(2)) write(nfile(2),*) '*** The electronic transition has been rejected. ***'
  else
    if(loutfile(1)) write(nfile(1),*) '*** The electronic transition has been accepted. ***'
    if(loutfile(2)) write(nfile(2),*) '*** The electronic transition has been accepted. ***'
end if

if( ltransition ) then
    !---Initialization of coefficients
    if( lrtddft ) then
        !---exciton dynamics
        dnmxr_ex = 0.d0
        dnmxi_ex = 0.d0
        dnmxr_ex(0,0,1:nexciton) = 1.d0
        accumprob_ex(:,:) = 0.d0
    else
        if( ldensity_matrix_repre ) then
            !--- density matrix representation
            dnmxr = 0.d0
            dnmxi = 0.d0
            dnmxr(0,0,1:nband,1:nspnmx) = 1.d0
            if( lfssh_dish ) then
                dish_prev_step(:,:) = -1
            end if
        else
            cr = 0.d0
            ci = 0.d0
            do is = 1, noccmx
               cr(is,is,1:nspnmx) = 1.d0
            end do
        end if
        accumprob = 0.d0
    end if
    lsignc = .false.
  else
    !---restore occupancies, if rejected.
    if( lrtddft ) then
        iband_hole(1:nexciton) =     iband_hole_(1:nexciton)
    iband_electron(1:nexciton) = iband_electron_(1:nexciton)
       ldegenerate(1:nexciton) =    ldegenerate_(1:nexciton)
    end if
    occ_rec = occ_store
    if( .not.lfssh_gsscf ) then
        occ = occ_rec
    end if
end if


return
end subroutine




subroutine tddft_fssh_vscale( nfile, myid, nodes,  &
& ltransition, couplvec, curekin, targetekin, vatm, frc, watom, &
& ntype, natom, nhk1, nhk2, dtmd )
!-----------------------------------------------------------------------
! TDDFT-FSSH: scale velocities
!-----------------------------------------------------------------------
use outfile
use constants
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
logical :: ltransition
integer :: ntype, natom
real*8,  dimension(3,natom) :: vatm, frc, couplvec
real*8  :: curekin, targetekin
real*8,  dimension(ntype)   :: watom
integer, dimension(ntype)   :: nhk1, nhk2
real*8  :: dtmd

!-----declare local variables
integer :: it, i
real*8  :: a, b, gamma, kernel, chkt


a = 0.d0
b = 0.d0
do it = 1, ntype
   do i = nhk1(it), nhk2(it)
      a = a + ( couplvec(1,i)*couplvec(1,i) + couplvec(2,i)*couplvec(2,i) &
&             + couplvec(3,i)*couplvec(3,i) )/watom(it)
      b = b + couplvec(1,i)*vatm(1,i) + couplvec(2,i)*vatm(2,i) &
&           + couplvec(3,i)*vatm(3,i)
   end do
end do
a = 0.5d0 * a

kernel = b*b + 4.d0*a*(targetekin - curekin)
ltransition = kernel >= 0.d0
if( .not.ltransition ) then
    if(loutfile(1)) write(nfile(1),*) '*** The kernel is negative:', kernel
    if(loutfile(2)) write(nfile(2),*) '*** The kernel is negative:', kernel
end if
if( ltransition ) then
    if( b < 0.d0 ) then
        gamma = (b + sqrt(kernel))/(2.d0*a)
      else
        gamma = (b - sqrt(kernel))/(2.d0*a)
    end if
    do it = 1, ntype
       do i = nhk1(it), nhk2(it)
          vatm(1:3,i) = vatm(1:3,i) - gamma*couplvec(1:3,i)/watom(it)
       end do
    end do

    !---check the current kinetic energy
    chkt = 0.d0
    do it = 1, ntype
       do i = nhk1(it), nhk2(it)
          chkt = chkt + watom(it)*( &
&         vatm(1,i)*vatm(1,i) + vatm(2,i)*vatm(2,i) + vatm(3,i)*vatm(3,i) )
       end do
    end do
    chkt = 0.5d0 * chkt
    if(loutfile(1)) write(nfile(1),*) '*** Check temperature [K] aftre scaling:', &
& chkt/natom* (tempau/2.d0) * 2.d0/3.d0
    if(loutfile(2)) write(nfile(2),*) '*** Check temperature [K] aftre scaling:', &
& chkt/natom* (tempau/2.d0) * 2.d0/3.d0
 
    !---correction for the velocity Verlet algorithm
    do it = 1, ntype
       do i = nhk1(it), nhk2(it)
          vatm(1:3,i) = vatm(1:3,i) - 0.5d0*frc(1:3,i)/watom(it)*dtmd
       end do
    end do
end if


return
end subroutine




subroutine read_tddft_fssh( nfile, iogpsz, occ, nstepMD,  &
& cgjr, rhcr, bufcr, iodg, ioag, idstnd,  &
& nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& npnod1, npnod2, npnod, nplcnt, npldsp, nspnmx )
!-----------------------------------------------------------------------
! read data in TDDFT-FSSH
!-----------------------------------------------------------------------
use outfile
use tddft_variables
use tddft_fssh_variables
implicit none
integer :: nfile(*)
integer :: iogpsz
integer :: nplwex, nplw, nspnmx
integer :: nband, nbnod1, nbnod2, nbnod, nbncnt(*), nbndsp(*)
real*8,  dimension(nband,nspnmx) :: occ
integer :: nstepMD
real*8  :: cgjr(nplwex,*)
real*8  :: rhcr(nplwex,*)
real*8  :: bufcr(*)
integer :: iodg(*), ioag(*), idstnd(*)
integer :: npnod1, npnod2, npnod, nplcnt(*), npldsp(*)

!-----declare local variables
integer :: iunit
integer :: myid, nodes, nkd
integer :: myid_qm_un, nodes_qm_un, ircv
integer :: ierror, istat
integer :: noband_, nspnmx_, nprevex_
real*8  :: dttddft_
logical :: lfssh_vscale_  = .false.
logical :: lfssh_dish_    = .false.
logical :: lset_gdcr_rec_ = .false.
logical :: ldensity_matrix_repre_, lrtddft_
integer :: nupdn_, nexciton_, nlrstates_, nlrreduce_
integer :: myid_lr, nodes_lr, myid_pw, nodes_pw
integer :: nversion, nbando, nplwo, nplwoex, nspin
real*8  :: dummy
real*8  :: ct, ct0, ct00, timecnt
logical :: lrdslmir
integer :: nfiles, ii, isnd, ibuf(4)


if( .not.ltddft_fssh ) return

ltddft_start = ltddft_start .and. nstepMD > 0
if( .not.ltddft_start ) return

ct00 = timecnt()

!-----set communicator
call get_worldqm( myid, nodes )

!-----set communicator
call get_worldqmun( myid_qm_un, nodes_qm_un )


ierror = 0

if( myid == 0 ) then

  ioif: if( mod(myid_qm_un,iogpsz) == 0 ) then ! always .ture. for no divide-and-conquer MD
     nfiles=min(iogpsz,nodes_qm_un-myid_qm_un) ! always   1    for no divide-and-conquer MD

    call allocate_unit_number( iunit )

    call open_tddftfsshfiles( nfile, myid, nodes,  &
& myid_qm_un, nodes_qm_un, iunit, ierror )

    read(iunit,iostat=istat) lrtddft_
    if(istat==0) then
    if( lrtddft .and. .not.lrtddft_ .or. .not.lrtddft .and. lrtddft_ ) then
        istat = 5
        write(nfile(1),*) 'ldensity_matrix_repre /= ldensity_matrix_repre_ in read_tddft_fssh',  &
& ldensity_matrix_repre, ldensity_matrix_repre_
        write(nfile(2),*) 'ldensity_matrix_repre /= ldensity_matrix_repre_ in read_tddft_fssh',  &
& ldensity_matrix_repre, ldensity_matrix_repre_
    end if
    end if

    if(istat==0) then
    read(iunit,iostat=istat) ldensity_matrix_repre_
    if( ldensity_matrix_repre .and. .not.ldensity_matrix_repre_ .or.  &
&  .not.ldensity_matrix_repre .and.      ldensity_matrix_repre_ ) then
        istat = 3
        write(nfile(1),*) 'ldensity_matrix_repre /= ldensity_matrix_repre_ in read_tddft_fssh',  &
& ldensity_matrix_repre, ldensity_matrix_repre_
        write(nfile(2),*) 'ldensity_matrix_repre /= ldensity_matrix_repre_ in read_tddft_fssh',  &
& ldensity_matrix_repre, ldensity_matrix_repre_
    end if
    end if

    if(istat==0) then
    read(iunit,iostat=istat) noband_, nspnmx_, dttddft_
    if( noband /= noband_ ) then
        istat = 1
        write(nfile(1),*) 'noband /= noband_ in read_tddft_fssh', noband, noband_
        write(nfile(2),*) 'noband /= noband_ in read_tddft_fssh', noband, noband_
    end if
    if( nspnmx /= nspnmx_ ) then
        istat = 2
        write(nfile(1),*) 'nspnmx /= nspnmx_ in read_tddft_fssh', nspnmx, nspnmx_
        write(nfile(2),*) 'nspnmx /= nspnmx_ in read_tddft_fssh', nspnmx, nspnmx_
    end if
    end if

    if(istat==0) read(iunit,iostat=istat) occ_rec(1:noband,1:nspnmx)
    if( lrtddft ) then
        if(istat==0) read(iunit,iostat=istat) nupdn_, nexciton_, nlrstates_, nlrreduce_
        if( nupdn_ > nupdn_ex ) then
            istat = 6
            write(nfile(1),*) 'nupdn_ > nupdn_ex in read_tddft_fssh', nupdn_ , nupdn_ex
            write(nfile(2),*) 'nupdn_ > nupdn_ex in read_tddft_fssh', nupdn_ , nupdn_ex
        end if
        if( nexciton_ /= nexciton ) then
            istat = 6
            write(nfile(1),*) 'nexciton_ /= nexciton in read_tddft_fssh', nexciton_ , nexciton
            write(nfile(2),*) 'nexciton_ /= nexciton in read_tddft_fssh', nexciton_ , nexciton
        end if
        if( nlrstates_ /= nlrstates ) then
            istat = 6
            write(nfile(1),*) 'nlrstates_ /= nlrstates in read_tddft_fssh', nlrstates_ , nlrstates
            write(nfile(2),*) 'nlrstates_ /= nlrstates in read_tddft_fssh', nlrstates_ , nlrstates
        end if
        if(istat==0) read(iunit,iostat=istat)  &
& iband_hole(1:nexciton), iband_electron(1:nexciton), ldegenerate(1:nexciton)
        if(istat==0) read(iunit,iostat=istat)  &
&                    dnmxr_ex(0:nupdn_,0:nupdn_,1:nexciton)
        if(istat==0) read(iunit,iostat=istat)  &
&                    dnmxi_ex(0:nupdn_,0:nupdn_,1:nexciton)
        if(istat==0) read(iunit,iostat=istat) accumprob_ex(0:nlrstates,1:nexciton)
        if(istat==0) read(iunit,iostat=istat) excite(0:nlrstates)
        if(istat==0) read(iunit,iostat=istat) weight(1:2,1:nlrstates)
        if(istat==0) read(iunit,iostat=istat) ibstate(1:2,0:nlrstates)
        if(istat==0) read(iunit,iostat=istat) lexcite(0:nlrstates)
    else
        if( ldensity_matrix_repre ) then
            if(istat==0) read(iunit,iostat=istat) nupdn_
            if( nupdn_ > nupdn ) then
                istat = 4
                write(nfile(1),*) 'nupdn_ > nupdn in read_tddft_fssh', nupdn_ , nupdn
                write(nfile(2),*) 'nupdn_ > nupdn in read_tddft_fssh', nupdn_ , nupdn
            end if
            if(istat==0) read(iunit,iostat=istat)  &
&                        dnmxr(-nupdn_:nupdn_,-nupdn_:nupdn_,1:noband,1:nspnmx)
            if(istat==0) read(iunit,iostat=istat)  &
&                        dnmxi(-nupdn_:nupdn_,-nupdn_:nupdn_,1:noband,1:nspnmx)
        else
            if(istat==0) read(iunit,iostat=istat) cr(1:noband,1:noband,1:nspnmx)
            if(istat==0) read(iunit,iostat=istat) ci(1:noband,1:noband,1:nspnmx)
        end if
        if(istat==0) read(iunit,iostat=istat) accumprob(1:noband,1:noband,1:nspnmx)
    end if
    if(istat==0) read(iunit,iostat=istat) lsignc(1:noband,1:nspnmx)

    if(istat==0) then
       read(iunit,iostat=istat) lfssh_vscale_
       if( istat/=0 ) then
           lfssh_vscale_ = .false.
           istat = 0
       end if
       if( lfssh_vscale_ ) then
        if(istat==0) read(iunit,iostat=istat) nprevex_, npreve
        if(istat==0) then
           if( nprevex /= nprevex_ ) then
               write(nfile(1),*) 'nprevex /= nprevex_ in read_tddft_fssh', nprevex, nprevex_
               write(nfile(2),*) 'nprevex /= nprevex_ in read_tddft_fssh', nprevex, nprevex_
           end if
           npreve = min(npreve,nprevex)
        end if
        if(istat==0) read(iunit,iostat=istat) totene(1:min(nprevex,nprevex_))
        if(istat==0) read(iunit,iostat=istat) theekin, thepote
       end if
    end if

    if(istat==0) then
       read(iunit,iostat=istat) lfssh_dish_
       if( istat/=0 ) then
           lfssh_dish_ = .false.
           istat = 0
       end if
       if( lfssh_dish_ ) then
           read(iunit,iostat=istat) dish_prev_step(1:noband,1:nspnmx)
       end if
    end if

    if(istat==0) then
       read(iunit,iostat=istat) lset_gdcr_rec_
       if( istat/=0 ) then
           lset_gdcr_rec_ = .false.
           istat = 0
       end if
    end if

    !-----for divide-and-conquer MD
!    if( nfiles > 1 ) then
!        call read_tddft_fssh_send( myid_qm_un, nodes_qm_un, istat, iunit, nfiles )
!    end if

  else ioif

    !-----for divide-and-conquer MD
!    ircv=(myid_qm_un/iogpsz)*iogpsz
!    call read_tddft_fssh_recv( myid_qm_un, nodes_qm_un, istat, ircv, nspnmx,  &
!& ldensity_matrix_repre_, lrtddft_, noband_, nspnmx_, dttddft_ ,  &
!& nupdn_, nexciton_, nlrstates_, nlrreduce_, lfssh_vscale_, nprevex_, npreve,  &
!& lfssh_dish_, lset_gdcr_rec_ )

  end if ioif

end if


!-----set communicator
call get_worldqm( myid, nodes )

call ibcast(istat,1,0)
if( istat /= 0 ) then
    ltddft_start = .false.
    if( myid == 0 ) then
  ioif2: if( mod(myid_qm_un,iogpsz) == 0 ) then ! always .ture. for no divide-and-conquer MD
        close(iunit)
        call deallocate_unit_number( iunit )
  end if ioif2
    end if

    !-----set communicator
    call get_worldkd( myid, nodes, nkd )

    return
end if

call dbcast(occ_rec,nband*nspnmx,0)
if( myid == 0 ) then
    if(loutfile(1)) write(nfile(1),*) 'read_tddft_fssh : successful'
    if(loutfile(2)) write(nfile(2),*) 'read_tddft_fssh : successful'
end if
if( lrtddft ) then
    if( myid == 0 ) then
        nlrreduce = nlrreduce_
    end if
    call ibcast(nlrreduce,1,0)
    call ibcast(iband_hole,nexciton,0)
    call ibcast(iband_electron,nexciton,0)
    call lbcast(ldegenerate,nexciton,0)
    call dbcast(dnmxr_ex,(nupdn_ex+1)*(nupdn_ex+1)*nexciton,0)
    call dbcast(dnmxi_ex,(nupdn_ex+1)*(nupdn_ex+1)*nexciton,0)
    call dbcast(accumprob_ex,(nlrstates+1)*nexciton,0)
    call dbcast(excite,nlrstates+1,0)
    call dbcast(weight,2*nlrstates,0)
    call ibcast(ibstate,2*(nlrstates+1),0)
    call lbcast(lexcite,nlrstates+1,0)
else
    if( ldensity_matrix_repre ) then
        call dbcast(dnmxr,(2*nupdn+1)*(2*nupdn+1)*nband*nspnmx,0)
        call dbcast(dnmxi,(2*nupdn+1)*(2*nupdn+1)*nband*nspnmx,0)
    else
        call dbcast(cr,nband*nband*nspnmx,0)
        call dbcast(ci,nband*nband*nspnmx,0)
    end if
    call dbcast(accumprob,nband*nband*nspnmx,0)
end if
call lbcast(lsignc,nband*nspnmx,0)
call lbcast(lfssh_vscale_,1,0)
if( lfssh_vscale_ ) then
    call ibcast(npreve,1,0)
    call dbcast(totene,npreve,0)
    call dbcast(theekin,1,0)
    call dbcast(thepote,1,0)
    lsetekin = .true.
end if
call lbcast(lfssh_dish_,1,0)
if( lfssh_dish_ ) then
    call ibcast(dish_prev_step,nband*nspnmx,0)
end if
noccmx = noband
if( .not.lfssh_gsscf ) then
    occ = occ_rec
end if


call lbcast(lset_gdcr_rec_,1,0)
lset_gdcr_rec = lset_gdcr_rec_

iflset_gdcr_rec: if( lset_gdcr_rec ) then

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

  ioif3: if( mod(myid_qm_un,iogpsz) == 0 ) then ! always .ture. for no divide-and-conquer MD
     nfiles=min(iogpsz,nodes_qm_un-myid_qm_un) ! always   1    for no divide-and-conquer MD

        read(iunit,iostat=istat) nversion    ! for backward compatibility
        if( istat == 0 ) read(iunit,iostat=istat) nbando
        if( istat == 0 ) read(iunit,iostat=istat) nplwo

        if( nbando /= nband ) then
            istat = 1
            write(nfile(1),*) 'nbando /= nband in read_tddft_fssh', nbando, nband
            write(nfile(2),*) 'nbando /= nband in read_tddft_fssh', nbando, nband
        end if
        if( nplwo /= nplw ) then
            istat = 1
            write(nfile(1),*) 'nplwo /= nplw in read_tddft_fssh', nplwo, nplw
            write(nfile(2),*) 'nplwo /= nplw in read_tddft_fssh', nplwo, nplw
        end if
        if( istat == 0 ) then
            nplwoex = 2*( nplwo + 1 )
        else
            ierror = 1
        end if

    !-----for divide-and-conquer MD
    do ii=1, nfiles-1
       ibuf(1) = nversion
       ibuf(2) = nbando
       ibuf(3) = nplwo
       ibuf(4) = ierror
       call cisend(myid_qm_un+ii,ibuf,4,myid_qm_un+ii,0)
    end do

  else ioif3

    !-----for divide-and-conquer MD
    isnd=(myid_qm_un/iogpsz)*iogpsz
    call cirecvs(myid_qm_un,ibuf,4,isnd,0)
    nversion = ibuf(1)
    nbando   = ibuf(2)
    nplwo    = ibuf(3)
    ierror   = ibuf(4)
    nplwoex = 2*( nplwo + 1 )

  end if ioif3

    end if

    !-----set communicator
    call get_worldkd( myid, nodes, nkd )

    call gimax(ierror)
    errif1: if( ierror == 0 ) then

        ct0 = timecnt()
        do nspin = 1, nspnmx
           call rdwfn22( nfile, myid, nodes, iogpsz, iunit, &
& rhcr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& cgjr, .true., nplwoex, nband, nband, .true., dummy, ierror )
           if( ierror.gt.0 ) exit

           ct = timecnt()
           if( myid.eq.0 ) then
               if(loutfile(1)) write(nfile(1),*) '                read & distrib. w.f.',  &
&                                ' : cpu-time :', ct-ct0
               if(loutfile(2)) write(nfile(2),*) '                read & distrib. w.f.',  &
&                                ' : cpu-time :', ct-ct0
           end if
           ct0 = ct

           !--- convert band decomposition to G decomposition ---
           call bdtogd( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& gdcr_rec(1,nspin), bufcr,  &
& npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iodg, idstnd, .false., .true. )

        end do

    end if errif1

    call gimax(ierror)

    errif2: if( ierror == 0 ) then

!        call rdslmir_tddft_fssh( nfile, myid, nodes, iogpsz, iunit,  &
!& lrdslmir, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, iodg, ioag, idstnd )
        lrdslmir = .false.

        if( lrdslmir ) then
            if( myid.eq.0 ) then
                if(loutfile(1)) write(nfile(1),*) 'read rdslmir_tddft_fssh : successful'
                if(loutfile(2)) write(nfile(2),*) 'read rdslmir_tddft_fssh : successful'
            end if
        else
            ierror = 1
        end if

    end if errif2

    errif3: if( ierror /= 0 ) then

        if( myid.eq.0 ) then
            write(nfile(1),*) 'error : in files QM_tddftfssh ...'
            write(nfile(2),*) 'error : in files QM_tddftfssh ...'
        end if

    end if errif3

    end if npwif


    if( nodes_pw > 1 ) then
        !-----set communicator
        call get_worldpw( myid, nodes )
        !-----internode synchronization
        call gsync

        call ibcast(ierror,1,0)
        if( ierror == 0 ) then
            do nspin = 1, nspnmx
               call dbcast(gdcr_rec(1,nspin),2*npnod*nband,0)
            end do
            call lbcast(lrdslmir,1,0)
!            if( lrdslmir ) then
!                call bcastslmir_tddft_fssh( nfile, myid, nodes )
!            end if
        end if

    end if

    end if nlrif

    if( nodes_lr > 1 ) then
        !-----set communicator
        call get_worldlr( myid, nodes )
        !-----internode synchronization
        call gsync

        call ibcast(ierror,1,0)
        if( ierror == 0 ) then
            do nspin = 1, nspnmx
               call dbcast(gdcr_rec(1,nspin),2*npnod*nband,0)
            end do
            call lbcast(lrdslmir,1,0)
!            if( lrdslmir ) then
!                call bcastslmir_tddft_fssh( nfile, myid, nodes )
!            end if
        end if

    end if

end if iflset_gdcr_rec


!-----set communicator
call get_worldqm( myid, nodes )

if( myid == 0 ) then
  ioif4: if( mod(myid_qm_un,iogpsz) == 0 ) then ! always .ture. for no divide-and-conquer MD
    close(iunit)
    call deallocate_unit_number( iunit )
  end if ioif4
end if


!-----set communicator
call get_worldkd( myid, nodes, nkd )

ct = timecnt()
if(loutfile(1)) write(nfile(1),*) '    read_tddft_fssh : cpu-time :', ct-ct00
if(loutfile(2)) write(nfile(2),*) '    read_tddft_fssh : cpu-time :', ct-ct00
ct0 = ct


return
end




subroutine save_tddft_fssh( nfile, iogpsz,  &
& cgjr, nplwex, nplw, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, nband,  &
& rhcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, ioag, idstnd, bufcr, nspnmx )
!-----------------------------------------------------------------------
! save data in TDDFT-FSSH
!-----------------------------------------------------------------------
use outfile
use tddft_variables
use tddft_fssh_variables
implicit none
integer :: nfile(*)
integer :: iogpsz
real*8  :: cgjr(*)
real*8  :: rhcr(*)
integer :: nplwex, nplw, nbnod1, nbnod2, nbnod, nbncnt(*), nbndsp(*), nband
integer :: npnod1, npnod2, npnod, nplcnt(*), npldsp(*)
integer :: iod(*), iodg(*), ioag(*), idstnd(*)
real*8  :: bufcr(*)
integer :: nspnmx

!-----declare local variables
integer :: iunit
integer :: myid, nodes, nkd
integer :: myid_qm_un, nodes_qm_un, isnd
integer :: myid_lr, nodes_lr, myid_pw, nodes_pw
integer :: nversion, nspin
real*8  :: ct0, ct, ct00, timecnt
integer :: ierror
integer :: nfiles


if( .not.ltddft_fssh ) return

ct00 = timecnt()

!-----set communicator
call get_worldqm( myid, nodes )

!-----set communicator
call get_worldqmun( myid_qm_un, nodes_qm_un )


if( myid == 0 ) then

  ioif: if( mod(myid_qm_un,iogpsz) == 0 ) then ! always .ture. for no divide-and-conquer MD
     nfiles=min(iogpsz,nodes_qm_un-myid_qm_un) ! always   1    for no divide-and-conquer MD

    call allocate_unit_number( iunit )

    call open_tddftfsshfiles( nfile, myid, nodes,  &
& myid_qm_un, nodes_qm_un, iunit, ierror )

    write(iunit) lrtddft
    write(iunit) ldensity_matrix_repre
    write(iunit) noband, nspnmx, dttddft
    write(iunit) occ_rec(1:noband,1:nspnmx)
    if( lrtddft ) then
        write(iunit) nupdn_ex, nexciton, nlrstates, nlrreduce
        write(iunit) iband_hole(1:nexciton), iband_electron(1:nexciton), ldegenerate(1:nexciton)
        write(iunit) dnmxr_ex(0:nupdn_ex,0:nupdn_ex,1:nexciton)
        write(iunit) dnmxi_ex(0:nupdn_ex,0:nupdn_ex,1:nexciton)
        write(iunit) accumprob_ex(0:nlrstates,1:nexciton)
        write(iunit) excite(0:nlrstates)
        write(iunit) weight(1:2,1:nlrstates)
        write(iunit) ibstate(1:2,0:nlrstates)
        write(iunit) lexcite(0:nlrstates)
    else
        if( ldensity_matrix_repre ) then
            write(iunit) nupdn
            write(iunit) dnmxr(-nupdn:nupdn,-nupdn:nupdn,1:noband,1:nspnmx)
            write(iunit) dnmxi(-nupdn:nupdn,-nupdn:nupdn,1:noband,1:nspnmx)
        else
            write(iunit) cr(1:noband,1:noband,1:nspnmx)
            write(iunit) ci(1:noband,1:noband,1:nspnmx)
        end if
        write(iunit) accumprob(1:noband,1:noband,1:nspnmx)
    end if
    write(iunit) lsignc(1:noband,1:nspnmx)
    write(iunit) lfssh_vscale
    if( lfssh_vscale ) then
        write(iunit) nprevex, npreve
        write(iunit) totene(1:nprevex)
        write(iunit) theekin, thepote
    end if
    write(iunit) lfssh_dish
    if( lfssh_dish ) then
        write(iunit) dish_prev_step(1:noband,1:nspnmx)
    end if

    write(iunit) lset_gdcr_rec

    !-----for divide-and-conquer MD
!    if( nfiles > 1 ) then
!        call save_tddft_fssh_recv( myid_qm_un, nodes_qm_un, iunit, nfiles )
!    end if

  else ioif

    !-----for divide-and-conquer MD
!    isnd=(myid_qm_un/iogpsz)*iogpsz
!    call save_tddft_fssh_send( myid_qm_un, nodes_qm_un, isnd, nspnmx, nband )

  end if ioif

end if


iflset_gdcr_rec: if( lset_gdcr_rec ) then

    !-----set communicator
    call get_worldlr( myid_lr, nodes_lr )

    !-----set communicator
    call get_worldpw( myid_pw, nodes_pw )

    !-----set communicator
    call get_worldkd( myid, nodes, nkd )

    nlrif: if( myid_lr == 0 ) then
    npwif: if( myid_pw == 0 ) then
    if( myid.eq.0 ) then

  ioif2: if( mod(myid_qm_un,iogpsz) == 0 ) then ! always .ture. for no divide-and-conquer MD

        nversion = -2
        write(iunit) nversion    ! for backward compatibility
        write(iunit) nband
        write(iunit) nplw

  end if ioif2

    end if

    do nspin = 1, nspnmx

       ct0 = timecnt()
       !--- copy prveig to rhcr
       call cpgdrh( nfile, myid, nodes, rhcr, gdcr_rec(1,nspin),  &
& npnod, nband )
       !--- to convert G decomposition to band decomposition ---
       call gdtobd( nfile, myid, nodes, ct0,  &
& cgjr, nplwex, nplw, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp,  &
& rhcr, bufcr, npnod1, npnod2, npnod, nplcnt, npldsp,  &
& iod, iodg, idstnd, .true. )

       !----- save previous wavefunctions
       call svwfn22( nfile, myid, nodes, iogpsz, iunit,  &
& rhcr, nplwex, nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, cgjr )

    end do

    !----- save previous slmir
!    call svslmir_tddft_fssh( nfile, myid, nodes, iogpsz, iunit,  &
!& nband, nbnod1, nbnod2, nbnod, nbncnt, nbndsp, iod, iodg, ioag, idstnd )

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

end if iflset_gdcr_rec


!-----set communicator
call get_worldqm( myid, nodes )

if( myid == 0 ) then
  ioif3: if( mod(myid_qm_un,iogpsz) == 0 ) then ! always .ture. for no divide-and-conquer MD
    close(iunit)
    call deallocate_unit_number( iunit )
  end if ioif3
end if


!-----set communicator
call get_worldkd( myid, nodes, nkd )

ct = timecnt()
if(loutfile(1)) write(nfile(1),*) '     save_tddft_fssh :          :', ct- ct00
if(loutfile(2)) write(nfile(2),*) '     save_tddft_fssh :          :', ct- ct00
ct0 = ct


return
end
