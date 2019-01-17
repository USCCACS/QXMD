



module remd_param
!-----------------------------------------------------------------------
! type declaration and initialization of input variables
!-----------------------------------------------------------------------
use allocated_memory
implicit none

integer :: npx = 1, npy = 1, npz = 1  ! No. of nodes in x,y,z directions
integer :: nproc                      ! No. of total nodes
integer :: iogpsz

logical :: lpureMD                    ! .true. = pure classical MD

real*8 :: anxi, anyi, anzi
real*8 :: sxog, syog, szog            ! Reduced node origin

real*8 :: rc_buffer = 0.d0            ! length of buffer region


!real*8  :: alloc_mem = 0.d0    ! counter for allocated memory

logical :: lstart = .false.    ! .true. = restart

integer :: ifmd = 0            !  0:non, 1:CG, 2:NVE-MD, 3:NVT-MD, 4:NPT-MD,
                               !  5:NVT for each atom, 6:NVT with GLE thermostat,
                               ! 10:MSST (multiscale shock technique)
real*8  :: dtmd = 100.d0       ! time step in [a.u.]
logical :: lmdstop = .false.   ! .true. = a STOP file is created
! If you wanna stop MD run, change .false. to .true. in the STOP file.
integer :: nstop= 10           ! total step
integer :: nstep_ini= 0        ! initial step number only for lstart == .false.
real*8  :: treq = 300.d0       ! temperature in [K]
real*8  :: hpext= 0.d0         ! pressure in [GPa]
logical :: liscale = .false.   ! .true. = check temperature
integer :: iscnum = 25         ! number of temperature check
integer :: iscstp = 20         ! skip step for temperature check
logical :: lmomzero = .false.  ! .true. = make momentum zero
integer :: ioptmze = 2         ! for structural optimization
                               !  -1: do not optimize atomic coordinates
                               !   0: Conjugate gradient
                               !   1: Projected velocity Verlet
                               !   2: Quasi-Newton method with BFGS formula
                               !   3: RFO saddle points search
                               !  10: Harmonic-mode analysis
                               !  11: Phonon-dispersion calculation
logical :: lpvv = .false.      ! .true. = Projected velocity Verlet
real*8  :: dist_max = 0.5d0    ! maximum displacement for ioptmze >= 2
!-----only for quasi-Newton method (ifmd==1 & ioptmze=2)
                               ! two stabilizer for quasi-Newton method
real*8  :: qnstabi    = 0.d0   ! if heigenval(i) <= qnstabi, set scalex(i) = 0.d0
real*8  :: gammamin   = 0.1d0  ! gamma <- min( gamma, gammamin )
                               !       in subroutine updtqn (remd.f90)
integer :: ibfgsclear =  0     ! clear Hessian every ibfgsclear step
                               ! if ibfgsclear == 0, Hessian is not cleared.
!-----only for Harmonic-mode analysis (ifmd==1 & ioptmze=10)
real*8  :: hmadisp  = 0.01d0   ! displacement [a.u.] 
integer :: nhmaord  = 1        ! order of numerical differentiation
logical :: lcentdif = .true.   ! .true. = central difference

integer :: ioptmze_cell = -1   ! for supercell optimization
                               !  -1: do not optimize supercell
                               !   0: Conjugate gradient
                               !   1: not used yet
                               !   2: Quasi-Newton method with BFGS formula
!-----only for Conjugate gradient method (ifmd==1 & ioptmze_cell==0)
real*8  :: dtcellcg = 0.1d0    ! CG time step
!-----only for quasi-Newton method (ifmd==1 & ioptmze_cell==2)
integer :: iclearcellh = 0     ! clear Hessian every iclearcellh step
                               ! if iclearcellh == 0, Hessian is not cleared.

!-----only for hybrid optimization (ifmd==1 & ioptmze>=0 & ioptmze_cell>=0 )
logical :: lhybridopt = .true. ! .true. = do structural optimization first
integer :: nstep_hybrid      = 10  ! time step for structural optimization
integer :: nstep_hybrid_cell = 10  ! time step for cell optimization

!-----thermostat parameters for NVT (ifmd == 3, 5, 6) & NPT (ifmd == 4)
integer :: nnos  = 1           ! # of thermostats
integer :: nresn = 1           ! MTS steps for heat bath
integer :: nyosh = 1           ! Yoshida-Suzuki decomposition step
real*8  :: tomega= 5500.d0     ! time scale for thermostat in [a.u.]
real*8  :: gnkt                ! 3*ntot(0)*treq
real*8  :: gkt                 ! treq

!-----thermostat parameters for NVT for each atom (ifmd == 5)
integer :: nnos_atom  = 1       ! # of thermostats
integer :: nresn_atom = 1       ! MTS steps for heat bath
integer :: nyosh_atom = 1       ! Yoshida-Suzuki decomposition step
real*8  :: tomega_atom= 5500.d0 ! time scale for thermostat in [a.u.]
real*8  :: g3kt                 ! 3*treq
logical :: lgthermo = .true.    ! global thermostat

!-----barostat parameters for NPT (ifmd == 4)
real*8  :: tbomega= 5500.d0    ! time scale for barostat in [a.u.]
real*8  :: blkmod =  250.d0    ! bulk modulus in [GPa]
real*8  :: volum0              ! volume to determine mass of barostat
real*8  :: gnd2kt, onf, bmass
integer :: irstrct = 0         ! restriction for MD cell in NPT-MD
!          0: no restriction (triclinic)
!          1: cubic          a  = b  = c, alpha = beta = gamma = 90
!          2: orthorhombic   a /= b /= c, alpha = beta = gamma = 90
!          3: tetragonal     a  = b /= c, alpha = beta = gamma = 90
!          4: monoclinic     a /= b /= c, alpha = beta = 90, gamma /= 90
!          5: hexagonal      a  = b /= c, alpha = beta = 90, gamma  = 120
!          6: trigonal       a  = b  = c, alpha = beta = gamma /= 90
!         10: by symmetry opearations (set automatically, when lsymop=.true. & irstrct==0)
! To enforce the lattice shape with symmetry opearations,
! add 10 to the above values. For example, 
!         12 for an orthorhombic cell with symmetry opearations.
integer :: irstrct_sub = 1     ! restriction for MD cell in NPT-MD
!      Sub-option for the restriction for MD cell
!      if irstrct == 3 (tetragonal),
!         irstrct_sub = 1 (default)...... a  = b /= c, alpha = beta = gamma = 90
!                     = 2          ...... b  = c /= a, alpha = beta = gamma = 90
!                    else          ...... c  = a /= b, alpha = beta = gamma = 90
!      if irstrct == 4 (monoclinic),
!         irstrct_sub = 1 (default)...... a /= b /= c, alpha = beta  = 90, gamma /= 90
!                     = 2          ...... a /= b /= c, beta  = gamma = 90, alpha /= 90
!                    else          ...... a /= b /= c, gamma = alpha = 90, beta  /= 90
!      if irstrct == 5 (hexagonal),
!         irstrct_sub = 1 (default)...... a  = b /= c, alpha = beta  = 90, gamma = 120
!                     = 2          ...... b  = c /= a, beta  = gamma = 90, alpha = 120
!                    else          ...... c  = a /= b, gamma = alpha = 90, beta  = 120
!
logical :: lcell_rstrct(3) = .false.    ! MD cell edge restriction
!      if lcell_rstrct(1) = .true., the length of L1 is fixed.
!      if lcell_rstrct(2) = .true., the length of L2 is fixed.
!      if lcell_rstrct(3) = .true., the length of L3 is fixed.

!-----barostat parameters for MSST (ifmd == 10)
real*8  :: shockspeed = 2000.d0 ! speed of shock wave in [m/s]
integer :: nshockv(1:3) = (/1,0,0/) ! direction with respect to super cell
                                    !  1 0 0 = Vs // L1 direction
                                    !  0 1 0 = Vs // L2 direction
                                    !  0 0 1 = Vs // L3 direction
                                    !  1 1 0 = Vs // L1 + L2 direction
                                    ! etc.
real*8  :: alpmatrix(3,3), betmatrix(3,3), h0(3,3)
real*8  :: xmsst = 1.d0, vxmsst = 0.d0
real*8  :: totmass
real*8  :: trsmatrix(3,3)           ! coordinate transformation natrix
logical :: lmsstscale = .false.     ! .true. = clear barostat velocity
integer :: msstscnum = 20           ! number of clear
integer :: msstscstp = 500          ! skip step


!-----parameter to output data
integer :: ioskip = 1          ! skip step
integer :: ioskipcoor = 1      ! skip step
integer :: ioskipvelo = 1      ! skip step
integer :: ioskipforc = 1      ! skip step
logical :: lmdout              !  = locoor .or. lovelo .or. loforc
logical :: locoor = .false.    ! .true. = output scaled coordinates
logical :: lovelo = .false.    ! .true. = output scaled velocities
logical :: loforc = .false.    ! .true. = output scaled forces

!-----parameters for output area
real*8  :: xio_min = 1.d0, xio_max = 0.d0 ! scaled x coordinates of output area
real*8  :: yio_min = 1.d0, yio_max = 0.d0 ! scaled y coordinates of output area
real*8  :: zio_min = 1.d0, zio_max = 0.d0 ! scaled z coordinates of output area
logical :: lrestrict_area

!-----parameter to output QM data
integer :: ioskip_qm = 1       ! skip step
logical :: lmdout_qm           !  = locoor_qm .or. lovelo_qm .or. loforc_qm
logical :: locoor_qm = .false. ! .true. = output scaled coordinates
logical :: lovelo_qm = .false. ! .true. = output scaled velocities
logical :: loforc_qm = .false. ! .true. = output scaled forces

!-----parameter to output MD-cluster data
integer :: ioskip_cl = 1       ! skip step
logical :: lmdout_cl           !  = locoor_cl .or. lovelo_cl .or. loforc_cl
logical :: locoor_cl = .false. ! .true. = output scaled coordinates
logical :: lovelo_cl = .false. ! .true. = output scaled velocities
logical :: loforc_cl = .false. ! .true. = output scaled forces


logical :: lsave   = .true.    ! .true. = save data
logical :: lsreal8 = .true.    ! .true. = in real*8 data

logical :: lstress = .false.   ! .true. = stress calculation
integer :: nskip_stress = 5    ! skip step
integer :: nskip_stress_out = 0 ! skip step

logical :: ltotmom = .false.   ! .true. = total-momentum calculation
integer :: nskip_totmom = 5    ! skip step


!-----parameters for soft walls
integer :: nwalls = 0           ! No. of soft walls
real*8,  allocatable, dimension(:,:) :: wallp, wallv
                               ! wallp(3,*) = any point on the wall
                               ! wallv(3,*) = the normal vector to the wall
real*8,  allocatable, dimension(:) :: wallf  ! strength of repulsive potential
integer, allocatable, dimension(:) :: nwallp ! power    of repulsive potential


!-----parameters for gravitational field
logical :: lgravi = .false.    ! .true. = apply gravitational field
real*8  :: gravmag = 0.d0      ! magnitude of gravitational field in unit of g = 9.8 [m/s^2]
real*8  :: gravg0  = 0.d0      ! g in [a.u.]
real*8  :: gravdir(3) = 0.d0   ! direction of gravitational field


logical :: lQMMDrun = .true.   ! .true. = QM/MD run

logical :: lstat = .true.      ! .true. = statistical calculation

logical :: latomic = .false.   ! .true. = output atomic stress & energy
integer :: nskip_atomic        ! skip step, should be equal to ioskip

logical :: lheat = .false.     ! .true. = calculate heat flux vector
integer :: nskip_heat = 1      ! skip step


!-----parameters for electric field
logical :: lefield = .false.   ! .true. = apply uniform electric field
logical :: lefield_start = .true. ! .true. = restart
real*8  ::  efield(3) = 0.d0   ! electric field vector
logical :: lsawtooth = .false. ! .true.  for isolated/periodic systems
                               ! .false. for periodic insulator systems
logical :: lsawtooth_shape = .false.  ! .true.  for periodic sawtooth potential
                                      ! .false. for non-periodic sawtooth potential
logical :: loutpolarization = .false. ! .true. = output polarization through SCF iterations
logical :: lconstraintD = .false. 


!-----tolerance for CG minimization (ifmd == 1)
real*8  :: tol_energy = 1.d-08 ! tolerance for energy/atom in [a.u.]
real*8  :: tol_force  = 5.d-04 ! tolerance for max. force  in [a.u.]

!-----constraint conditions
integer :: ncbonds = 0                                    ! No. of bonds constrained
integer, allocatable, dimension(:) :: ncbatm1, ncbatm2    ! atom No.
real*8,  allocatable, dimension(:) :: cblength            ! constraint length
real*8,  allocatable, dimension(:) :: cblambda, cblambdav ! Lagrange constant


integer :: natom  = 0                     ! the number of QM atoms
integer :: ntype  = 0                     ! the number of QM atomic species
real*8, dimension(3,3) :: hcell = 0.d0    ! supercell vectors for QM region
real*8 :: rvol                            ! volume of supercell
integer :: nHSQM  = 0                     ! the number of handshake QM atoms
integer :: ntHSQM = 0                     ! the number of handshake QM atomic species

integer :: nmd  = 0                       ! the number of MD atoms per node
integer :: nmd_alloc = 0                  ! allocated array size for MD atoms
integer :: ntmd = 0                       ! the number of MD atomic species
integer :: ntmd_ist = 0                   ! the number of reduced MD atomic species
real*8, dimension(3,3,0:1) :: h = 0.d0    ! supercell vectors for MD region
real*8 :: volume                          ! volume of MD cell
integer :: nHSMD  = 0                     ! the number of handshake MD atoms
integer :: ntHSMD = 0                     ! the number of handshake MD atomic species

integer :: nmd_cluster  = 0               ! the number of MD cluster atoms
integer :: ntmd_cluster = 0               ! the number of MD cluster atomic species


!-----parameters for virtual MD for thermodynamic integration
logical :: lvmd = .false.       ! .true. = virtual MD
real*8  :: dlambda_vmd = 1.d0   ! weight of real system
!    An MD simulation is performed with the virtual potential
!     U = dlambda_vmd*U(real) + (1-dlambda_vmd)*U(reference).
logical :: lidealref = .true.   ! .true.  = ideal gas/Einstein solid
                                ! .false. = classical-potential system specified
!integer :: ntmd_vmd = 0        ! the number of lattice-site species
integer :: nmd_vmd = 0          ! the number of lattice sites
logical :: leinstein = .false.  ! .true. = Einstein solid exists

save


end module




module remd_param_atom
!-----------------------------------------------------------------------
! type declaration and initialization of variables for atoms
!-----------------------------------------------------------------------
implicit none


!-----variables for QM atoms
integer :: natomx, ntypex                      ! array size
real*8,  allocatable, dimension(:) :: zatom    ! atomic number
integer, allocatable, dimension(:) :: nhk      ! the number of QM atoms
integer, allocatable, dimension(:) :: nhk1, nhk2
integer, allocatable, dimension(:) :: icscale     ! 1:scaled, 2:real coordinates
real*8,  allocatable, dimension(:,:) :: ratm      ! atomic scaled coordinates
real*8,  allocatable, dimension(:,:) :: realatm   ! atomic real coordinates
real*8,  allocatable, dimension(:,:) :: vatm      ! atomic scaled velocities
real*8,  allocatable, dimension(:) :: watom       ! mass number
logical, allocatable, dimension(:) :: lterminator ! .true. = QM terminator (HSQM)

!-----variables for MD atoms
integer :: nmdmax                                 ! array size for # of atoms in each node
integer :: npb = 0                                ! # of buffer atoms
integer :: nmd_buffer, nmd_buffer_1d              ! array size for # of buffer atoms
real*8,  allocatable, dimension(:) :: zatom_md    ! atomic number
integer, allocatable, dimension(:) :: icscale_md  ! 1:scaled, 2:real coordinates
real*8,  allocatable, dimension(:,:) :: x      ! atomic scaled coordinates & velocities
integer, allocatable, dimension(:) :: is       ! atomic species for each atom
integer, allocatable, dimension(:) :: ntot     ! the number of atoms
real*8,  allocatable, dimension(:) :: watom_md ! mass number
logical, allocatable, dimension(:) :: vrandom  ! .true.=random velocities
logical, allocatable, dimension(:) :: lfixion  ! .true.=fix atomic positions
real*8,  allocatable, dimension(:) :: zatom_md_ist ! atomic number
integer, allocatable, dimension(:) :: is2ist       ! translation table to true species
logical, allocatable, dimension(:) :: ldispion    ! .true.=update position with constant d
real*8,  allocatable, dimension(:,:) :: dvector   ! constant displacement vectors
real*8,  allocatable, dimension(:) :: ficmass     ! mass = watom_md * ficmass
real*8, dimension(3) :: disp = 0.d0               ! displacement vector

!-----variables for MD cluster atoms
integer :: nmdx_cluster, ntmdx_cluster              ! array size
real*8,  allocatable, dimension(:) :: zatom_cluster    ! atomic number
integer, allocatable, dimension(:) :: nhk1_cluster, nhk2_cluster
integer, allocatable, dimension(:) :: icscale_cluster  ! 1:scaled, 2:real coordinates
real*8,  allocatable, dimension(:,:) :: x_cluster   ! atomic scaled coordinates & velocities
integer, allocatable, dimension(:) :: is_cluster    ! atomic species for each atom
integer, allocatable, dimension(:) :: ntot_cluster  ! the number of atoms
real*8,  allocatable, dimension(:) :: watom_cluster ! mass number
logical, allocatable, dimension(:) :: lMDterminator ! .true. = MD terminator (HSMD)
integer, allocatable, dimension(:) :: is2ist_cluster ! translation table to true species

!-----variables for MD handshake atoms
real*8,  allocatable, dimension(:) :: zatom_HSMD    ! atomic number
integer, allocatable, dimension(:) :: icscale_HSMD  ! 1:scaled, 2:real coordinates
real*8,  allocatable, dimension(:) :: x_HSMD        ! atomic scaled coordinates
integer, allocatable, dimension(:) :: ntot_HSMD     ! the number of atoms


!-----vaeiables for thermostat only for NVT-MD (ifmd == 3 )
!-----related to nnos
real*8,  allocatable, dimension(:) :: xlogs    ! Heat-bath coordinates
real*8,  allocatable, dimension(:) :: vlogs    ! Heat-bath velocities
real*8,  allocatable, dimension(:) :: glogs    ! Heat-bath accelerations
real*8,  allocatable, dimension(:) :: qmass    ! Heat-bath masses
!-----related to nyosh
real*8,  allocatable, dimension(:) :: wdti     ! Time-steps for YS integration
real*8,  allocatable, dimension(:) :: wdti2
real*8,  allocatable, dimension(:) :: wdti4
real*8,  allocatable, dimension(:) :: wdti8

!-----vaeiables for thermostat only for NVT-MD for each atom (ifmd == 5 or 6)
logical, allocatable, dimension(:) :: lathermo
                     ! lathermo =  .true. : atom type with resepective thermostat
                     !          = .false. : atom type with global thermostat
!-----related to nnos_atom
real*8,  allocatable, dimension(:,:) :: xlogs_atom ! Heat-bath coordinates
real*8,  allocatable, dimension(:,:) :: vlogs_atom ! Heat-bath velocities
real*8,  allocatable, dimension(:)   :: glogs_atom    ! Heat-bath accelerations
real*8,  allocatable, dimension(:,:) :: qmass_atom ! Heat-bath masses
!-----related to nyosh_atom
real*8,  allocatable, dimension(:) :: wdti_atom    ! Time-steps for YS integration
real*8,  allocatable, dimension(:) :: wdti2_atom
real*8,  allocatable, dimension(:) :: wdti4_atom
real*8,  allocatable, dimension(:) :: wdti8_atom


!-----vaeiables for virtual MD for thermodynamic integration
integer, allocatable, dimension(:) :: nvmd      ! reference system
                                                ! 1:ideal gas, 2:Einstein solid
real*8,  allocatable, dimension(:) :: omega_vmd ! reference system
real*8,  allocatable, dimension(:) :: x_vmd     ! atomic scaled coordinates
                                                ! of lattice points
integer, allocatable, dimension(:) :: ntot_vmd  ! the number of lattice sites
!integer, allocatable, dimension(:) :: icscale_vmd  ! 1:scaled, 2:real coordinates

save


end module




!subroutine remd_get_lstress( nfile, myid, nodes, lstress_ )
!!-----------------------------------------------------------------------
!!    check input data
!!-----------------------------------------------------------------------
!use remd_param
!implicit none
!integer :: nfile(*), myid, nodes
!logical :: lstress_
!
!lstress_ = lstress
!
!return
!end subroutine




subroutine remd_param_atom_alloc( nfile, myid, nodes,  &
& alloc_mem, ntype, natom, ntmd, nmd, nmd_alloc, b,  &
& rc_buffer, anxi, anyi, anzi, ntHSMD, nHSMD )
!-----------------------------------------------------------------------
!     allocate memory for atoms
!-----------------------------------------------------------------------
use remd_param_atom
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: ntype, natom
integer :: ntmd, nmd
integer :: nmd_alloc
real*8, dimension(3,3) :: b
real*8  :: rc_buffer
real*8  :: anxi, anyi, anzi
integer :: ntHSMD, nHSMD

!-----declare local variables
integer :: ntmax
integer :: ntHSMDx, nHSMDx
integer :: status
real*8  :: the_mem
integer :: it, i, j
real*8  :: vala, valb, valc, rca, rcb, rcc, vbuffer


!-----determine size of array for the atomic coordinates
nmdmax = nmd
call gimax(nmdmax)
nmdmax = dble(nmdmax) * 1.5d0 + 100          ! <---- should be optimized

vala = dsqrt( b(1,1)**2 + b(2,1)**2 + b(3,1)**2 )
valb = dsqrt( b(1,2)**2 + b(2,2)**2 + b(3,2)**2 )
valc = dsqrt( b(1,3)**2 + b(2,3)**2 + b(3,3)**2 )
rca  = rc_buffer*vala / anxi
rcb  = rc_buffer*valb / anyi
rcc  = rc_buffer*valc / anzi

!-----scaled volume of buffer region in 1 direction (node-cell volume = 1)
vbuffer = max( rca, (1.d0 + 2.d0*rca )*rcb,  &
&                   (1.d0 + 2.d0*rca )*( 1.d0 + 2.d0*rcb )*rcc )
nmd_buffer_1d = nmdmax * vbuffer          ! <---- should be optimized
nmd_buffer_1d = max( 1, nmd_buffer_1d )

!-----scaled volume of buffer region (node-cell volume = 1)
vbuffer = ( 1.d0 + 2.d0*rca )*( 1.d0 + 2.d0*rcb )  &
&        *( 1.d0 + 2.d0*rcc ) - 1.d0
nmd_buffer = nmdmax * vbuffer             ! <---- should be optimized
nmd_buffer = max( 1, nmd_buffer )

nmd_alloc  = nmdmax + nmd_buffer


ntmax  = max( ntype, ntmd )
ntypex = max( 1, ntype )
natomx = max( 1, natom )


!------allocate memory for QM atoms
allocate( zatom(ntypex),  &
& nhk(ntypex), nhk1(ntypex), nhk2(ntypex), icscale(ntypex),  &
& ratm(3,natomx), realatm(3,natomx), vatm(3,natomx),  &
& watom(ntypex), lterminator(ntypex),  &
& stat=status )

!------allocate memory for MD atoms
if( status == 0 ) allocate( zatom_md(ntmd),  &
& x(3*nmd_alloc,0:1), is(nmd_alloc),  &
& ntot(0:ntmd), icscale_md(ntmd), watom_md(ntmd),  &
& vrandom(ntmax), lfixion(ntmax), is2ist(ntmd), zatom_md_ist(ntmd),  &
& ldispion(ntmax), dvector(3,ntmax), ficmass(ntmax),  &
& stat=status )

ntHSMDx = max( 1, ntHSMD )
 nHSMDx = max( 1,  nHSMD )

!------allocate memory for MD atoms
if( status == 0 ) allocate( zatom_HSMD(ntHSMDx),  &
& x_HSMD(3*nHSMDx), ntot_HSMD(ntHSMDx),  &
& icscale_HSMD(ntHSMDx),  &
& stat=status )

the_mem =  &
&  8.d0 * ntypex  &
& + 4.d0 * ntypex*6  &
& + 8.d0 * ( 3*natomx*3 + ntypex ) + 1.d0 * ntypex  &
& + 8.d0 * ntmd  &
& + 8.d0 * ( 3*nmd_alloc*2 )  &
& + 4.d0 * ( nmd_alloc + 2*ntmd + 1 )  &
& + 8.d0 * ntmd  &
& + 1.d0 * ( 2*ntmax )  &
& + 4.d0 * ( ntmd )  &
& + 8.d0 * ( ntmd )  &
& + 8.d0 *  ntHSMDx  &
& + 8.d0 * 3*nHSMDx  &
& + 4.d0 *   nHSMDx  &
& + 4.d0 *  ntHSMDx * 2

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'remd_param_atom_alloc', .true. )


!-----set initial values
do it = 1, ntypex
   zatom(it)  = 0.d0     ! atomic number
   nhk(it)    = 0          ! No. of atoms
   nhk1(it)   = 0
   nhk2(it)   = 0
   icscale(it) = 1       ! 1:scaled, 2:real coordinates
   watom(it)  = 0.d0     ! mass number
   lterminator(it) = .false. ! .true.= QM terminator (HSQM)
end do

do j = 1, 3
do i = 1, natomx
   ratm(j,i) = 0.d0      ! atomic coordinates
   vatm(j,i) = 0.d0      ! atomic velocities
end do
end do

do it = 1, ntmd
   zatom_md(it)   = 0.d0     ! atomic number
   icscale_md(it) = 1
   watom_md(it)   = 0.d0
   ntot(it)   = 0
end do
ntot(0)   = 0
do it = 1, ntmax
   vrandom(it) = .true.   ! .true.=random velocities
   lfixion(it) = .false.  ! .true.=fix atomic positions
   ficmass(it) = 1.d0     ! fictitious mass
end do
do it = 1, ntmax
   ldispion(it) = .false. ! .true.=fix atomic displacements
   do i = 1, 3
      dvector(i,it) = 0.d0
   end do
end do

do i = 1, nmd_alloc*3
   x(i,0) = 0.d0
   x(i,1) = 0.d0
end do
do i = 1, nmd_alloc
   is(i)  = 0
end do

do i = 1, ntHSMDx
   zatom_HSMD(i) = 0.d0
   ntot_HSMD(i)  = 0
   icscale_HSMD(i) = 1
end do
do i = 1, nHSMDx*3
   x_HSMD(i) = 0.d0
end do


return
end subroutine




subroutine remd_param_atom2_alloc( nfile, myid, nodes,  &
& alloc_mem, ntype, natom )
!-----------------------------------------------------------------------
!     allocate memory for MD cluster atoms
!-----------------------------------------------------------------------
use remd_param_atom
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: ntype, natom

!-----declare local variables
integer :: ntmax, namax
integer :: status
real*8  :: the_mem
integer :: it, i, j


ntmax = max( 1, ntype )
namax = max( 1, natom )
ntmdx_cluster = ntmax
nmdx_cluster  = namax

!------allocate memory for MD cluster atoms
allocate( zatom_cluster(ntmax),  &
& nhk1_cluster(ntmax), nhk2_cluster(ntmax), icscale_cluster(ntmax),  &
& x_cluster(3*namax,0:1), is_cluster(namax),  &
& ntot_cluster(0:ntmax), watom_cluster(ntmax),  &
& lMDterminator(ntmax), is2ist_cluster(ntmax),  &
& stat=status )

the_mem =  &
&  8.d0 * ntmax  &
& + 4.d0 * ntmax * 3  &
& + 8.d0 * 3*namax*2  &
& + 4.d0 * namax  &
& + 4.d0 * ( ntmax + 1 )  &
& + 8.d0 * ntmax  &
& + 1.d0 * ntmax  &
& + 4.d0 * ntmax

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'remd_param_atom2_alloc', .true. )


!-----set initial values
do it = 1, ntmax
   zatom_cluster(it)  = 0.d0     ! atomic number
   nhk1_cluster(it)  = 0
   nhk2_cluster(it)  = 0
   watom_cluster(it)  = 0.d0     ! mass number
   lMDterminator(it) = .false.      ! .true.= MD terminator (HSQM)
   ntot_cluster(it)  = 0
end do
ntot_cluster(0)  = 0

do i = 1, namax*3
   x_cluster(i,0) = 0.d0
   x_cluster(i,1) = 0.d0
end do
do i = 1, namax
   is_cluster(i)  = 0
end do


return
end subroutine




subroutine remd_thermostat_alloc( nfile, myid, nodes,  &
& alloc_mem, nnos, nresn, nyosh )
!-----------------------------------------------------------------------
!     allocate memory for thermostat variables
!-----------------------------------------------------------------------
use remd_param_atom
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: nnos, nresn, nyosh

!-----declare local variables
integer :: status
real*8  :: the_mem


!------allocate memory
allocate( xlogs(nnos), vlogs(nnos), glogs(nnos), qmass(nnos),  &
& wdti(nyosh), wdti2(nyosh), wdti4(nyosh), wdti8(nyosh),  &
& stat=status )

the_mem =  &
& + 8.d0 * ( nnos*4 + nyosh*4 )

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'remd_thermostat_alloc', .true. )


return
end subroutine




subroutine remd_atom_thermostat_alloc( nfile, myid, nodes,  &
& alloc_mem, nnos, nresn, nyosh, nmd_alloc )
!-----------------------------------------------------------------------
!     allocate memory for thermostat variables
!-----------------------------------------------------------------------
use remd_param_atom
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: nnos, nresn, nyosh, nmd_alloc

!-----declare local variables
integer :: status
real*8  :: the_mem


!------allocate memory
allocate( xlogs_atom(nnos,nmd_alloc), vlogs_atom(nnos,nmd_alloc),  &
& glogs_atom(nnos), qmass_atom(nnos,nmd_alloc),  &
& wdti_atom(nyosh), wdti2_atom(nyosh), wdti4_atom(nyosh), wdti8_atom(nyosh),  &
& stat=status )

the_mem =  &
& + 8.d0 * ( size(xlogs_atom) + size(vlogs_atom) + size(glogs_atom) &
&          + size(qmass_atom) + size(wdti_atom) + size(wdti2_atom)  &
&          + size(wdti4_atom) + size(wdti8_atom) )

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'remd_atom_thermostat_alloc', .true. )


return
end subroutine




module remd_variables
!-----------------------------------------------------------------------
! type declaration of shared variables in remd.f
!-----------------------------------------------------------------------
implicit none

integer :: myx, myy, myz

integer :: nstepCG = 0
integer :: nstepMD = 0
integer :: nstep = 0

real*8  :: epot   = 0.d0  ! classical energy for entire  system
real*8  :: epotCL = 0.d0  ! classical energy for cluster system
real*8  :: epotQM = 0.d0  ! QM        energy for cluster system
real*8  :: epotref = 0.d0 ! potential energy for reference ideal system in virtual MD

real*8  :: prevene = 1.d+10
integer :: isccnt = 0
integer :: msstsccnt = 0

character(50):: fname_toc   ! data-file name for the next run
character(50):: fname_mts

integer :: imts = 0

real*8, dimension(3,3) :: pintlr=0.d0 ! stress by  long-range potentials
real*8, dimension(3,3) :: pintsr=0.d0 ! stress by short-range potentials
real*8, dimension(3,3) :: pint  =0.d0 ! stress tensor
real*8, dimension(3,3) :: pit1  =0.d0 ! stress tensor for output
real*8, dimension(3,3) :: pext  =0.d0 ! external stress tensor
real*8, dimension(3,3) :: pintlrref=0.d0 ! stress from reference system in virtual MD
real*8, dimension(3,3) :: pintsrref=0.d0 ! 

real*8  :: cellfnorm = 0.d0
real*8  :: hessian_cell(9,9) = 0.d0, savefrc_cell(3,3) = 0.d0
logical :: lqninitial_cell = .true.

save


end module




module remd_atom
!-----------------------------------------------------------------------
! type declaration of variables for atoms
!-----------------------------------------------------------------------
implicit none


!-----variables for MD atoms
real*8, allocatable, dimension(:) :: frc

!-----variables for optimization (ifmd == 1)
integer :: nmdmax_
real*8, allocatable, dimension(:) :: fnorm
real*8, allocatable, dimension(:) :: pixcg

real*8, allocatable, dimension(:) :: fack  ! Prefactors for kinetic energy
real*8, allocatable, dimension(:) :: acon  ! Prefactors for forces

!-----variables for (NPT)-MD
real*8, allocatable, dimension(:) :: fack2  ! Prefactors for kinetic energy
real*8, allocatable, dimension(:) :: xu, vu

!-----variables for atomic energy & stress
integer :: nmdmax__
real*8, allocatable, dimension(:) :: wepot
real*8, allocatable, dimension(:) :: wstrs, wstrsk

!-----variables for QM/MD or pure classical runs
real*8, allocatable, dimension(:) :: xrec   ! to calculate atomic displacements

integer, allocatable, dimension(:,:) :: lspr   ! pair list
integer :: nprmax                              ! array size for lspr

integer, allocatable, dimension(:) :: ist      ! true atomic species for each atom


!-----variables for QM atoms
real*8, allocatable, dimension(:) :: frcQM
!-----variables for QM stress
real*8, dimension(3,3) :: strQM


!-----variables for quasi-Newton/Rational function optimazation (RFO) methods
real*8, allocatable, dimension(:) :: deltax, deltaf
real*8, allocatable, dimension(:) :: hessian, heigenvec, heigenval
real*8, allocatable, dimension(:) :: scalex, scaleg
logical :: lqninitial = .true.
!--      integer :: ibfgsclear = 10 ! Hessian is clear every ibfgsclear step
logical :: lhessian      = .false. ! lhessian = .true. if Hessian has negative eigenvalues
logical :: lhessian_cell = .false.

!-----variables for virtual MD for thermodynamic integration
real*8, allocatable, dimension(:) :: frcref

save


end module




module remd_qmmd
!-----------------------------------------------------------------------
! type declaration of variables for QM/MD method
!-----------------------------------------------------------------------
implicit none


!----- ind_*(1,i) = MD-node id
!----- ind_*(2,i) = MD-atom id in the MD node
!----- x_*        = coordinates


!----- ind_qm1_md   for QM (not HSQM) atom
integer, allocatable, dimension(:,:) :: ind_qm1_md

!----- ind_HSQM_qm, x_HSQM_qm   for QM atom related to HSQM
!----- ind_HSQM_md, x_HSQM_md   for MD atom related to HSQM
integer :: nHSQMx
integer, allocatable, dimension(:,:) :: ind_HSQM_qm
real*8,  allocatable, dimension(:)   :: x_HSQM_qm
integer, allocatable, dimension(:,:) :: ind_HSQM_md
real*8,  allocatable, dimension(:)   :: x_HSQM_md
real*8,  allocatable, dimension(:)   :: alpha_HSQM


!----- ind_md1_md   for MD-cluster (not HSMD) atom
integer, allocatable, dimension(:,:) :: ind_md1_md

!----- ind_HSMD_qm, x_HSMD_qm   for QM atom related to HSMD
!----- ind_HSMD_md, x_HSMD_md   for MD atom related to HSMD
integer :: nHSMDx
integer, allocatable, dimension(:,:) :: ind_HSMD_qm
real*8,  allocatable, dimension(:)   :: x_HSMD_qm
integer, allocatable, dimension(:,:) :: ind_HSMD_md
real*8,  allocatable, dimension(:)   :: x_HSMD_md
real*8,  allocatable, dimension(:)   :: alpha_HSMD

!-----forces for MD-cluster atoms
real*8,  allocatable, dimension(:) :: frcCL
integer, allocatable, dimension(:) :: ist_cluster ! true atomic species for each atom
integer, allocatable, dimension(:,:) :: lspr_cluster   ! pair list

save


end module




subroutine remd_atom_alloc( nfile, myid, nodes,  &
& alloc_mem, nmd, ntype, nmdmax, nmd_alloc, natom,  &
& ifmd, lQMMDrun, lpureMD, volume, rc_buffer, ioptmze, lstress, latomic,  &
& lheat, lvmd, lidealref )
!-----------------------------------------------------------------------
!     allocate memory for variables for atoms
!-----------------------------------------------------------------------
use remd_atom
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: nmd
integer :: ntype
integer :: nmdmax
integer :: nmd_alloc
integer :: natom
integer :: ifmd
logical :: lQMMDrun, lpureMD
real*8  :: volume, rc_buffer
integer :: ioptmze
logical :: lstress, latomic, lheat, lvmd, lidealref

!-----declare local variables
integer :: status
real*8  :: the_mem
integer :: natomx
integer :: nmdtmp


natomx = max( 1, natom )

!------allocate memory
allocate( frc(3*nmdmax),  &
& fack(ntype), acon(ntype),  &
& frcQM(3*natomx),  &
& stat=status )

the_mem =  &
&  8.d0 * ( size(frc) + size(fack) + size(acon) + size(frcQM) )


if( lvmd ) then
    !------allocate memory
    allocate( frcref(3*nmdmax), stat=status )
    the_mem = the_mem + 8.d0 * size(frcref) 
    frcref(:) = 0.d0
end if


if( status == 0 ) then
    if( ifmd == 1 ) then
        nmdmax_ = nmdmax
    else
        nmdmax_ = 1
    end if
    allocate( fnorm(nmdmax_), pixcg(3*nmdmax_),  &
& stat=status )

    the_mem = the_mem  &
& + 8.d0 * ( size(fnorm) + size(pixcg) )
end if


!-----initial guess for max # of atoms around each atom
nprmax = nmdmax * nodes / volume * 4.5d0*rc_buffer**3
nprmax = nprmax * 2
nprmax = max( 1, nprmax )

if( status == 0 ) then
    if( lQMMDrun .or. lpureMD .or. .not.lidealref ) then
        allocate( xrec(3*nmdmax), lspr(0:nprmax,nmdmax), ist(nmd_alloc),  &
& stat=status )

    the_mem = the_mem  &
& + 8.d0 * ( nmdmax*3 )  &
& + 4.d0 * ( (nprmax+1)*nmdmax + nmd_alloc )
    end if
end if


if( status == 0 ) then
if( ifmd == 4 .or. ifmd == 10 ) then
    !-----for (NPT)-MD
    allocate( fack2(ntype), vu(3*nmdmax), xu(3*nmdmax),  &
& stat=status )

    the_mem = the_mem  &
& + 8.d0 * ( ntype )  &
& + 8.d0 * ( 3*nmdmax*2 )
end if
end if


if( status == 0 ) then
    if( latomic .or. lheat ) then
        nmdmax__ = nmdmax
    else
        nmdmax__ = 1
    end if
    allocate( wepot(nmdmax__), wstrs(6*nmdmax__), wstrsk(6*nmdmax__),  &
& stat=status )

    the_mem = the_mem  &
& + 8.d0 * ( size(wepot) + size(wstrs) + size(wstrsk) )
end if


if( status == 0 .and. ifmd == 1 .and. ioptmze >= 2 ) then
    !-----error trap
    if( nodes > 1 ) then
        call fstop( nfile, myid, nodes,  &
& 'not supported yet for parallel calculations for RFO methods' )
    end if

    !-----for quasi-Newton/RFO methods
    allocate( deltax(3*nmd), deltaf(3*nmd),  &
& hessian(3*nmd*3*nmd), heigenvec(3*nmd*3*nmd),  &
& heigenval(3*nmd),  &
& scalex(3*nmd), scaleg(3*nmd),  &
& stat=status )

    the_mem = the_mem  &
& + 8.d0 * ( size(deltax) + size(deltaf)  &
&          + size(hessian) + size(heigenvec) + size(heigenval)  &
&          + size(scalex) + size(scaleg) )
end if


!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'remd_atom_alloc', .true. )


!if( ifmd == 6 ) then
!    !-----GLE thermostat
!    call glethemo_alloc( nfile, myid, nodes, alloc_mem, nmdmax )
!end if


!-----clear variables
pixcg(:) = 0.d0
fnorm(:) = 0.d0
if( lQMMDrun .or. lpureMD .or. .not.lidealref ) then
    xrec(1:3*nmdmax) = 0.d0
end if


return
end subroutine




subroutine remd_qmmd_alloc( nfile, myid, nodes,  &
& alloc_mem, natom, nHSQM, nmd_cluster, nHSMD, nprmax )
!-----------------------------------------------------------------------
!     allocate memory for variables for QM/MD method
!-----------------------------------------------------------------------
use remd_qmmd
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: natom
integer :: nHSQM
integer :: nmd_cluster
integer :: nHSMD
integer :: nprmax

!-----declare local variables
integer :: nmdx_cluster
integer :: status
real*8  :: the_mem
integer :: natomx


nHSQMx = max( 1, nHSQM )
nHSMDx = max( 1, nHSMD )
nmdx_cluster = max( 1, nmd_cluster )
natomx = max( 1, natom )

!------allocate memory
allocate( ind_qm1_md(2,natomx),  &
& ind_HSQM_qm(2,nHSQMx), x_HSQM_qm(3*nHSQMx),  &
& ind_HSQM_md(2,nHSQMx), x_HSQM_md(3*nHSQMx),  &
& alpha_HSQM(nHSQMx),  &
& ind_md1_md(2,nmdx_cluster),  &
& ind_HSMD_qm(2,nHSMDx), x_HSMD_qm(3*nHSMDx),  &
& ind_HSMD_md(2,nHSMDx), x_HSMD_md(3*nHSMDx),  &
& alpha_HSMD(nHSMDx),  &
& frcCL(3*nmdx_cluster), ist_cluster(nmdx_cluster),  &
& lspr_cluster(0:nprmax,nmdx_cluster),  &
& stat=status )

the_mem =  &
&  4.d0 * ( natomx*2 )  &
& + 4.d0 * ( nHSQMx*2 )  &
& + 8.d0 * ( nHSQMx*3 )  &
& + 4.d0 * ( nHSQMx*2 )  &
& + 8.d0 * ( nHSQMx*3 )  &
& + 8.d0 * ( nHSQMx )  &
& + 4.d0 * ( nmdx_cluster*2 )  &
& + 4.d0 * ( nHSMDx*2 )  &
& + 8.d0 * ( nHSMDx*3 )  &
& + 4.d0 * ( nHSMDx*2 )  &
& + 8.d0 * ( nHSMDx*3 )  &
& + 8.d0 * ( nHSMDx )  &
& + 8.d0 * 3*nmdx_cluster  &
& + 4.d0 *   nmdx_cluster

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'remd_qmmd_alloc', .true. )


frcCL(1:3*nmdx_cluster) = 0.d0


return
end subroutine
