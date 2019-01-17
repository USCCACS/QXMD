



module param
!-----------------------------------------------------------------------
! type declaration and initialization of input variables
!-----------------------------------------------------------------------
use allocated_memory
implicit none

integer :: npx = 1, npy = 1, npz = 1  ! No. of nodes in x,y,z directions
integer :: nproc                      ! No. of total nodes
integer :: iogpsz

!real*8  :: alloc_mem = 0.d0    ! counter for allocated memory

logical :: lstart = .false.    ! .true. = restart

logical :: ltddft = .false.    ! .true. = execute MD based on TDDFT
                               ! .false.= execute MD based on conventional DFT
logical :: ltddft_start = .false.  ! .true. = restart
logical :: ltddft_fssh = .true.    ! .true. = FSSH, .false. = Ehrenfest
logical :: lfssh_switch = .true.   ! .true.  = switching aveilable
                                   ! .false. = cccupations are fixed
logical :: lfssh_gsscf = .true.    ! .true.  = SCF with the ground  state
                                   ! .false. = SCF with the excited state
logical :: lfssh_gsfrc = .false.   ! .true.  = the ground  state forces are used
                                   ! .false. = the excited state forces are used.
                                   ! only for lfssh_gsscf = .true.
integer :: imxchg_fssh = 1         ! charge mixing = 0:non, 1:Pulay, 2:Anderson, 3:simple
real*8  :: aslh_fssh = 0.9d0       ! mixing parameters, only for lfssh_gsscf = .true.
real*8  :: bslh_fssh = 0.6d0       ! = aslh, bslh, when 0.d0.
logical :: lfssh_random = .false.  ! .true. = manual, .false. = automatic
real*8  :: rseed_fssh = 25d0       ! seed for random numbers, only for lfssh_random = .true.
logical :: lfssh_boltzmn = .true.  ! .true. = multiply the upward probability
                                   !          by the Boltzmann factor
logical :: lfssh_vscale = .false.  ! .true. = scale velocity when switching
real*8  :: tminimum = 0d0          ! minimum temperature in [Kelvin] to accept switching
real*8  :: tdbroad = 0d0           ! width of Gaussian broadening in [Kelvin]
                                   ! 0.0: no broadening
real*8  :: dttddft = 0.02d0        ! time step in  [a.u.]  Caution !! [Rydberg units]
logical :: lfssh_parallel = .true. ! .true. = do parallel calculation
logical :: ltddft_nscforce = .true.! .true. = Non-self consistent (NSC) force

logical :: lrtddft = .false.       ! .true. = calculate Casida coupling matrix
real*8  :: fdmeshxc = 1.d-06       ! FD mesh for the 2nd derivative of Exc in [a.u.]
logical :: lrspecific = .false.    !  .true.  = for specified states
real*8  :: enediff = 0.3d0         ! max energy diff between occ & unocc states in [Ryd.]
real*8  :: threshrpa = 0.d0        ! threshold value for RPA1 (accelerate calculations)
real*8  :: threshdipole = 0.d0     ! threshold value for the transition dipole
real*8  :: threshexcite = 1.d-04   ! threshold value for the difference
                                   !           between the energies of hole and particle
logical :: llcexchange = .false.   ! long-range exchange correction scheme
real*8  :: alpha_llcexchange = 0.17d0   ! paramter to divide 1/r in to short/long parts
real*8  :: large_cell_llcexchange = 2d0 ! paramter to give a CUBIC large cell
                                        !  The cell-vector length is given by 
                                        !  large_cell_llcexchange * max(L1,L2,L3).
real*8  :: hybrid_mixing = 0d0     ! mixing paramter for Hybrid method
real*8  :: foersigma = 0.1d0       ! width of Gaussian to calculate Foerster transfer rate
logical :: ldiagonal = .true.      ! 
!      if ltddft == .true. .and. ltddft_fssh == .true., set ldiagonal = .true.
!      if ldiagonal == .true., the matrix elements are calculated only for
!      ihole_1 == ihole_2 .and. iparticle_1 == iparticle_2
!logical :: lexiton_dynamics = .false.   !  = ltddft_fssh .and. lrtddft
logical :: lscissors_lrtddft = .false. ! scissors correction to the excitation energy
real*8  :: ecorr_tirplet = 0.d0    ! correction to the triplet state
real*8  :: ecorr_singlet = 0.d0    ! correction to the singlet state
logical :: lfission_rate = .false. ! singlet-fission-rate calculation
logical :: lemission_rate = .true. ! spontaneous-emission-rate calculation
real*8  :: refractive = 1.5d0      ! the refractive index
logical :: ltda = .false.            ! .true. = Tamm-Dancoff approximation (TDA)
logical :: lcic_tda = .false.        ! .true. = configuration-interaction-corrected (CIC) TDA
real*8  :: cic_tda_coupling = 0.01d0 ! Coupling between the S_0 and S_1 states in the CIC-TDA

logical :: lpaw = .true.       ! .true. = the projector-augmented wave method
                               ! .false.= the pseudopotential method
logical :: lpaw_sym = .false.  ! .true. = full      symmetry calculation
                               ! .false.= spherical symmetry calculation
                               !          for on site energy

logical :: lrela = .false.     ! .true. =     relativistic calculation
                               ! .false.= non relativistic calculation
! Note:     relativistic date sets = control/RPAW or RUSPP or RNCPP
!       non relativistic date sets = control/ PAW or  USPP or  NCPP
logical :: lrela_full = .false. ! .true. = fully relativistic calculation
                                !          with two-component spinors
                                ! .false.= scalar relativistic calculation
logical :: lrela_sop = .true.   ! only for lrela_full = .false.
! Hamiltonian matrix including the spin-orbit interaction is diagonalized,
! and the results are recorded to data/qm_eig_SO.d

logical :: lsoc = .false.      ! only for lpaw = .true. & lrela = .false.
                               ! .true. = with    spin-orbit coupling
                               ! .false.= without spin-orbit coupling
logical :: lsoc_full = .false. ! only for lnoncollinear = .true.
! For no spin or collinear spin, always lsoc_full = .false.,
! and Hamiltonian matrix including the spin-orbit interaction is diagonalized,
!     and the results are recorded to data/qm_eig_SO.d


integer :: jgga = 2            ! 1:LDA, 2:GGA(PBE), 3:GGA(RPBE)
                               ! 4:revPBE, 5:van der Waals DFT(vdW-DF), 6:vdW-DF2
integer :: jgga_org = 2
logical :: lvdw = .false.
logical :: lvdw_factorize = .true.
logical :: ldftd = .false.     ! on/off DFT-D empirical correction for vdW
logical :: lplusU = .false.    ! on/off the mean-field Hubbard model to take into account 
                               !        the strong on-site Coulomb repulsion (DFT+U)
logical :: lplusC = .false.    ! on/off an empirical correction to non-local pp.(DFT+C)
integer :: jhybrid = 0         ! 0:none, 1:HSE with PBE
real*8  :: hmixing = 0.25d0    ! EX = hmixing*Ex^HF + (1 - hmixing)*Ex^PBE
real*8  :: hrange = 0.11d0     !  1/r = erfc(hrange*r)/r + erf(hrange*r)/r
                               !          Short Range        Long Range
integer :: mgridred  = 1       ! grid-reduction factor
                               ! The number of k-point grids should be the multiple of mgridred.

logical :: lscissors = .false.  ! on/off scissors correction at donor/acceptor interface
                                ! to eigenvalues
real*8  :: ecorr_acceptor = 0d0 ! corrections for acceptor
real*8  :: ecorr_donor = 0d0    ! corrections for donor
real*8  :: zboundary = 0d0      ! boundary in the z direction

integer :: iscfmx = 100        ! maximum No. of global iterations
integer :: itrial = 0          ! No. of trial global iterations
real*8  :: tolpot   = 5.0d-09  ! tolerance for total energy
real*8  :: tolres   = 5.0d-09  ! tolerance for average residual
real*8  :: tolHF_KS = 5.0d-09  ! tolerance for HF-KS functional difference
integer :: imxchg_trial = 0    ! trial charge mixing = 0:non, 1:Pulay, 2:Anderson, 3:simple
real*8  :: aslh_trial = 0.1d0  ! trial mixing parameters
real*8  :: bslh_trial = 0.64d0
integer :: imxspin_trial = 0   ! trial charge mixing = 0:non, 1:Pulay, 2:Anderson, 3:simple
real*8  :: amxspin_trial = 0.1d0  ! trial mixing parameters
real*8  :: bmxspin_trial = 0.64d0

integer :: imxchg = 1          ! charge mixing = 0:non, 1:Pulay, 2:Anderson, 3:simple,
                               !                 4:Srivastava, 5:Johnson
real*8  :: aslh = 0.8d0        ! mixing parameters
real*8  :: bslh = 0.64d0
real*8  :: aslh0               ! hold the original mixing parameters
integer :: itratn = 10         ! # of mixing
real*8  :: cmetric = 20.d0     ! the metric in charge mixing
!  The weighting factor w(g) in the calculation of residual is given by
!     w(g) = (cmetric-1)*G_min^2*G_max^2/(G_max^2 - cmetric*G_min^2)
logical :: louteig = .false.   ! output eigenvalues  through SCF iterations
logical :: loutenergy = .false.! output energy parts through SCF iterations
logical :: loutzansa = .false. ! output residuals    through SCF iterations
logical :: loutlog = .false.   ! output detailed log through SCF iterations
real*8  :: paw_mix = 0.5d0     ! onsite charge mixing parameters in the PAW method
real*8  :: plusU_mix = 0.5d0   ! onsite charge mixing parameters in the PAW method & DFT+U
integer :: nstabi = 10         ! convergence stabilizer
real*8  :: xstabi = 10.d0      ! If zansa2 < zansa2/xstabi at nstabi iterations before,
                               ! set lreset = .false.

logical :: lmixreal = .true.   ! only for lnoncollinear = .true.
                               !  .true.  =       real-space mixing
                               !  .false. = recoprocal-space mixing
integer :: imxspin = 0         ! spin-density mixing = 0:non, 1:Pulay, 2:Anderson, 3:simple,
                               !                       4:Srivastava, 5:Johnson
real*8  :: amxspin = 0.8d0     ! mixing parameters
real*8  :: bmxspin = 0.64d0
integer :: nmxspin = 10        ! # of spin-density mixing
real*8  :: spinmetric = 1.d0   ! the metric in spin-density mixing
real*8  :: wkerker = 0.d0      ! Kerker-mixing weight in spin-density mixing
!  Caution!!  sigma_new(G) = sigma_in(G) + K(G)*(sigma_out(G) - sigma_in(G)),
!             where sigma is the spin density, and 
!             K(G) = amxspin*(1 - wkerker + wkerker*G~2/(G^2 + bmxspin)).
!  Note that sigma(G=0) is not updated when wkerker = 1.d0,
!  i.e., the difference between the numbers of up- and down-spin is fixed.

logical :: lhcunt = .true.     ! .true. = HC products by Unitary tr.
integer :: ihldam = 1          ! 1:CG, 2:RMM-DIIS
real*8  :: toleig = 1.d-10     ! tolerance for Kohn-Sham equation
real*8  :: roughfeig = 0.3d0   ! rough tolerance factor for eigenvalues
real*8  :: roughfzan = 0.3d0   ! rough tolerance factor for residuals
real*8  :: threinn = 0.0d0     ! threshold for direction of new w.f.
real*8  :: wdinner = 0.9d0     ! window for direction of new w.f.
integer :: itermx  = 4         ! maximum No. of iterations
integer :: iteremx = 4         ! maximum No. of iterations for empty bands
integer :: ihyoji  = 100       ! display every 100 iterations
integer :: kstrial = 0         ! No. of trial global iterations with itermx = iteremx = 1
integer :: methodcg = 2        ! 1:line minimization, 2:BKL

integer :: multg = 2           ! multigrid level
real*8  :: tolcg = 1.d-11      ! tolerance energy for Poisson equation
real*8  :: weigrd= 0.4d0       ! preconditioner
integer :: nd2v  = 6           ! order of numerical differentiation
integer :: msrhmx = 14         ! mesh for serial calculation
integer :: msrhmy = 14         !
integer :: msrhmz = 14         !

integer :: ifmd = 0            ! 0:non, 1:CG, 2:NVE-MD, 3:NVT-MD, 4:NPT-MD,
                               ! 10:MSST (multiscale shock technique)
real*8  :: dtmd = 50.d0        ! time step in [a.u.]  Caution !! [Rydberg units]
integer :: nstop= 10           ! total step
integer :: nstep_ini= 0        ! initial step number only for lstart == .false.
real*8  :: treq = 300.d0       ! temperature in [K]
logical :: liscale = .false.   ! .true. = check temperature
integer :: iscnum = 25         ! number of temperature check
integer :: iscstp = 20         ! skip step
logical :: lmomzero = .false.  ! .true. = make momentum zero
integer :: ichest = 3          ! # of previous steps for c.d. estimation
integer :: ihest  = 3          ! # of previous steps for w.f. estimation
logical :: laspc_chg = .true.  ! .true. = the ASPC predictor / .false. = usual extrapotation
logical :: laspc_wv  = .true.  ! .true. = the ASPC predictor / .false. = subspace alignment
                               ! These are for backward compatibility
logical :: laspc = .false.     ! .true. = the ASPC method / .false. = CGMD
logical :: lxlbomd = .false.   ! .true. = time reversible XL-BOMD
integer :: kxlbomd = 0         ! # of previous values in XL-BOMD, if 0, kxlbomd = ihest
real*8  :: aspc_corr =  1.d0   ! ASPC corrector (0: default)
real*8  :: aslh_aspc =  0.d0   ! Charge mixing in the ASPC method (0.0 -> set aslh_aspc = aslh)
real*8  :: bslh_aspc =  0.d0   ! Charge mixing in the ASPC method (0.0 -> set bslh_aspc = bslh)
integer :: iscfmx_aspc = 1     ! No. of corrector step
real*8  :: tolres_aspc = 1.0d-06 ! tolerance for average residual in the ASPC method
integer :: nskp_aspc = 0       ! correction skip step (0: no correction)

integer :: ioptmze             ! for structural optimization (default value is given in remd_in.f90)
integer :: ioptmze_cell        ! for supercell  optimization (default value is given in remd_in.f90)

integer :: nshockv(1:3) = (/1,0,0/) ! direction with respect to super cell
real*8  :: trsmatrix(3,3)           ! coordinate transformation natrix
real*8  :: hcell_rec(3,3)           ! store original supercell

logical :: lsave = .true.      ! .true. = save data
logical :: lsreal8 = .true.    ! .true. = in real*8 data

logical :: lmulken = .false.   ! .true. = Mulliken analysis
integer :: nskpmulk= 5         ! skip step
logical :: ldecmpovp = .false. ! .true. = decompose overlap population
logical :: lspdatom = .false.  ! .true. = output weights for each band associated
                               !          with each atom
logical :: lpmchgat = .false.  ! .true. = output m-decomposed atomic charge
logical :: lmdecmpdos = .false.! .true. = output m-decomposed PDOS

logical :: lsphexp = .false.   ! .true. = spherical harmonics expansion
integer :: nskpsphexp = 5      ! skip step
real*8  :: radsphexp = 0.d0    ! default radius [a.u.],
                               ! if radsphexp = 0, slater_r are used

logical :: leda    = .false.   ! .true. = Energy density analysis
integer :: nskpeda = 5         ! skip step
real*8  :: radgeda = 10.d0     ! radius for grid EDA [a.u.]
!      real*8  :: rdelaunay = 20.d0   ! radius for Delaunay tetrahedralization [a.u.]

logical :: lwannier= .false.   ! .true. = maximully localized Wannier func.
integer :: iwbnd1=0, iwbnd2=0  ! band index to be transformed by Unitary matrix
integer :: iterwan = 200       ! maximum No. of iterations
real*8  :: tolwan  = 1.d-06    ! tolerance
integer :: iwstt1=1, iwstt2=0  ! orbit index
integer :: natwan  = 0         ! No. of focused atoms
integer, allocatable, dimension(:) :: natwno, nwanaa  ! atom No., No. of WF
integer :: nskpwann= 5         ! skip step
logical :: loutuni = .false.   ! .true. = output unitary matrix

logical :: lconduct = .false.  ! .true. = conductivity calculation
integer :: nskpconduct = 5     ! skip step
real*8  :: wgconduct = 0.1d0   ! width of Gaussian used alternative to Delta func.
real*8  :: tempconduct = 300d0 ! Temperature
logical :: ldcconduct = .false.! .true. = Spatial distribution of DC current density
real*8  :: efconduct(3) = 0d0  ! Constant electric-field vector
integer :: icband1  = 0, icband2  = 0     ! range of band indices for total current
integer :: ichband1 = 0, ichband2 = 0     ! range of band indices for hole current
integer :: iceband1 = 0, iceband2 = 0     ! range of band indices for electron current
logical :: lacconduct = .false.! .true. = Frequency dependence of conductivity
real*8  :: freqacmx  = 1d0     ! maximum frequency
integer :: immtband1 = 0, immtband2 = 0   ! range of band indices for mometum matrix
!- icband1, icband2, ichband1, ichband2, iceband1, iceband2, immtband1, immtband2 = 0
!- gives an automatic selection from energy (wgconduct, freqacmx) and temperature (tempconduct)

logical :: lstress = .false.   ! .true. = stress calculation
integer :: nskip_stress = 5    ! skip step

logical :: lintchg = .false.   ! .true. = atomic charge
integer :: nskip_intchg = 5    ! skip step

logical :: ldpchg = .false.    ! .true. = dump charge density
integer :: nskip_dpchg = 5     ! skip step
real*8  :: dc_xmin = 1.d0, dc_xmax = 0.d0 ! scaled x coordinates of output area
real*8  :: dc_ymin = 1.d0, dc_ymax = 0.d0 ! scaled y coordinates of output area
real*8  :: dc_zmin = 1.d0, dc_zmax = 0.d0 ! scaled z coordinates of output area
logical :: ldpchg_paw = .false.  ! .true. = dump soft-charge density ( only for PAW )
logical :: lcompress_dpchg = .true. ! .true. = compressed format
integer :: ndigit_dpchg = 3         ! significant figures

logical :: ldpwav = .false.    ! .true. = dump wave functions
integer :: nskip_dpwav = 5     ! skip step
integer :: ibstt1=0, ibstt2=0  ! band index
real*8  :: wv_xmin = 1.d0, wv_xmax = 0.d0 ! scaled x coordinates of output area
real*8  :: wv_ymin = 1.d0, wv_ymax = 0.d0 ! scaled y coordinates of output area
real*8  :: wv_zmin = 1.d0, wv_zmax = 0.d0 ! scaled z coordinates of output area
logical :: lcompress_dpwav = .true. ! .true. = compressed format
integer :: ndigit_dpwav = 3         ! significant figures

logical :: ldppot = .false.    ! .true. = dump local potential
integer :: nskip_dppot = 5     ! skip step
integer :: nav_dppot = 1       ! 1:xy, 2:yz, 3:zx
real*8  :: pt_xmin = 1.d0, pt_xmax = 0.d0 ! scaled x coordinates of output area
real*8  :: pt_ymin = 1.d0, pt_ymax = 0.d0 ! scaled y coordinates of output area
real*8  :: pt_zmin = 1.d0, pt_zmax = 0.d0 ! scaled z coordinates of output area
logical :: lcompress_dppot = .true. ! .true. = compressed format
integer :: ndigit_dppot = 3         ! significant figures

logical :: lhoteh = .false.    ! .true. = hot electron/hole relaxation rate
integer :: nskip_hoteh = 5     ! skip step
integer :: ibhoth = 1          ! # of states below Ef
integer :: ibhote = 1          ! # of states above Ef

real*8, dimension(3,3) :: h_MD  = 0.d0    ! MD cell vectors
real*8, dimension(3,3) :: hcell = 0.d0    ! supercell vectors
real*8, dimension(3,3) :: h_bkup = 0.d0   ! supercell vectors in the previous step
real*8, dimension(3,3) :: hci             ! (hci)^transpose = (inverse of hcell)
real*8, dimension(3,3) :: rdelg           ! supercell vectors / n_i
logical :: lorthrhmbc = .true.            ! .true. = orthorhombic supercell
real*8, dimension(3) :: rdel              ! the grid spacing,
                                          ! only valid for orthorhombic cells
real*8 :: rdelv                           ! volume elements for dense grids
real*8 :: volume                          ! volume of supercell

logical, dimension(3) :: lvacuum= .false. ! set vacuum regions
real*8,  dimension(3) :: vacuum = 0.d0
real*8  :: alpha_ldouble_grid_recip = 0.d0 ! the parameter in the error function
logical :: lshdep = .false.

real*8  :: ecut     = 0.d0     ! cutoff energy for wavefunctions
real*8  :: ecutdens = 0.d0     ! cutoff energy for charge density
real*8  :: ecutsoft = 0.d0     ! cutoff energy for soft part of density
real*8  :: ecutorth = 0.d0     ! cutoff energy for orthogonalization in CG iteration
real*8  :: ecutc    = 0.d0     ! cutoff energy for charge density mixing
! cutoff energies for long-range Coulomb potential
real*8  :: ecutlong = 0.d0        ! for non-periodic direction
real*8  :: ecutlong_small = 0.d0  ! for     periodic direction
integer :: iosp = 3            ! order of cardinal B spline in the double-grid method
                               !  iosp must be an odd integer

logical :: lvshape = .false.   ! .true. = variable-shape method
real*8  :: pwscale = 1.2d0     ! scale factor for # of PW's in the (NPT) MD

integer :: noband = 0          ! No. of occupied bands
integer :: nband  = 0          ! total No. of electronic bands
integer :: lfermi = 3          ! 1:nonmetallic, 2:Fermi, 3:Gaussian,
                               ! lfermi > 3:Methfessel & Paxton, order = lfermi - 3
real*8  :: tfermi = 2000.d0    ! electronic temp.(K), if metallic

logical :: lspin  = .false.    ! .true. = take into account spin polarization
logical :: lnoncollinear = .false. ! .true. = noncollinear magnetism
                                   ! .false.=    collinear magnetism
integer :: ncscale = 1         ! scale factor for # of PW's
logical :: lfixud = .true.     ! .true. = fix diff. between N_up & N_down
real*8  :: diffud = 0.d0       ! diff. between N_up & N_down
logical :: lwfrand = .false.   ! .false. = same as up-spin,
                               ! .true.  = by random No.
logical :: lduplicate = .false. ! .true. = duplicate the up-spin state
logical :: latmref = .false.    ! .true. = the direction of the magnetization of each atom is restricted.
logical :: lclearamag = .false. ! .true. = the magnetization of each atom is cleared.
real*8  :: refmx = 0.d0, refmy = 0.d0, refmz = 0.d0 
                                   ! reference direction for noncollinear magnetism
logical :: lclearmagne_in_x = .false.
logical :: lclearmagne_in_y = .false.
logical :: lclearmagne_in_z = .false.
! If lclearmagne_in_x, y, z == true., clear magnetization
! in the x, y, and z direction, respectively, to zero.
logical :: lfixmagne_in_x = .false.
logical :: lfixmagne_in_y = .false.
logical :: lfixmagne_in_z = .false.
real(8) ::  fixmagne_in_x = 0d0
real(8) ::  fixmagne_in_y = 0d0
real(8) ::  fixmagne_in_z = 0d0
! If lfixmagne_in_x, y, z == true., fix the sum of magnetization
! in the x, y, and z direction, respectively.

!-----parameters for electric field
logical :: lefield = .false.      ! .true. = apply uniform electric field
logical :: lefield_start = .true. ! .true. = restart
real*8  ::  efield(3) = 0.d0      ! electric field vector
logical :: lsawtooth = .false.    ! .true.  for isolated/periodic systems
                                  ! .false. for periodic insulator systems
logical :: lsawtooth_shape = .false.  ! .true.  for periodic sawtooth potential
                                      ! .false. for non-periodic sawtooth potential
logical :: lsawtooth_xyz(3) = .false. !
!# lsawtooth = .true. is only for rectangular supercell
!#                        and for the efield is parallel to the x or y or z direction.
!#   periodic sawtooth potential = x E_0     for   0 <= x < L/2
!#                               = (L-x) E_0 for L/2 <= x < L
! if lsawtooth_xyz(n) == .true., efield // x or y or z for n = 1, 2, 3, respectively
logical :: lefield_islts       ! = lefield .and. .not.lsawtooth
logical :: loutpolarization = .false. ! .true. = output polarization through SCF iterations
logical :: lconstraintD = .false. 
!# if lconstraintD = .true., the electric displacement D is variable (fixed),
!#                           and D is given by efield(1:3).
!# otherwise the electric field E is variable (fixed).

integer :: ntype  = 0          ! the number of atomic species
integer :: natom  = 0          ! the number of atoms
integer :: mxl    = 0          ! the maximum angular momentum
integer :: nylmmx              ! = mxl*(mxl+2) + 1
logical :: lplcom = .false.    ! .true. = place the c.m. on the center of cell


integer :: node_c = 1          ! No. of nodes for Cholesky decomposition
integer :: node_r = 1          ! No. of nodes for eigenvalue problem


real*8,  dimension(3) :: rmax
real*8,  dimension(3) :: cofmas   ! the center of supercell

integer, dimension(3) :: nd1v     ! the number of FD  meshes
integer, dimension(3) :: nd1vks   ! the number of FFT meshes

integer :: nel                    ! the number of electrons
integer :: ncion = 0              ! the charge number of ion

integer :: nhit1 = 1              ! order of preconditioning to solve Poisson eq.

!------for double grid method for solving Poisson equation in real space
logical :: ldouble_grid            ! .true. = use double grid method
real*8,  dimension(3) ::  disp_org ! displacement of double cells
integer, dimension(3) :: ndisp_org ! displacement (unit of mesh) of double cells
logical :: lsphere = .false.       ! .true. = use spherical region in cluster calculation

!------for double grid method for solving Poisson equation in reciprocal space
logical :: ldouble_grid_recip  = .false. ! .true. = use double grid method
integer :: idouble_grid_method = 0       ! 1 = surface, 2 = wire, 3 = cluster geometries

real*8  :: dblock = 1.d+07     ! size of work array in [bytes] used in untryt and sbalin
integer :: iblock              ! for work area used in untryt and sbalin
                               ! If you have large memory, increase iblock to get faster.

integer :: mx1loc = 2**13      ! dimension for local-pseudopotential tables
integer :: mx1    = 2**13      ! dimension for other tables

logical :: lwell = .false.     ! well potential outside the atomic sphere
real*8  :: wellheight = 1.d0   ! height of well potential in Ryd.

!---variables for k points
integer :: nkpnt = 1                ! the number of k points
integer :: nknod1 = 1, nknod2 = 1, nknod = 1    ! k-point index in each node
real*8,  allocatable, dimension(:,:) :: bzk
real*8,  allocatable, dimension(:)   :: wbzk
logical :: lgamma = .true.
logical, allocatable, dimension(:)   :: lgammak
integer :: npkx(3) = 1              ! the number of division for reciprocal vectors
integer, allocatable, dimension(:,:) :: kpsymop
logical :: ltetra = .false.         ! for BZ-zone integration
                     ! Note: Modified tetrahedron method (PRB49,16223(1994)) is used.
                     !       Fermi energy is initially estimated from lfermi.
logical :: lksgamma = .false.       !  .true. = Gamma centered sampling
                                    ! .false. = Monkhorst-Pack sampling
integer :: mxkpdupnm = 1

!---restrictin for MD cell
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


integer :: natom_alloc = 0     ! the number of atoms
integer :: natnod_alloc  = 0   ! the number of atoms per node

save


end module




module param_inispin
!-----------------------------------------------------------------------
! type declaration and initialization of variables for spin initialization
!-----------------------------------------------------------------------
implicit none

integer :: inispin = 1         ! initial spin density
                               !  1: uniformly polarized (default)
                               !  2: ferromagnetic alignment
                               !  3: antiferromagnetic random alignment
                               !  4: specify spin polarization
integer :: nspinat = 0         ! # of atomic types to be polarized for inispin = 2 & 3
                               ! # of atoms to be polarized for inispin = 4
integer, allocatable, dimension(:) :: iatspin  ! type for inispin = 2 & 3
                                               ! atom # for inispin = 4
real*8,  allocatable, dimension(:) :: atmmagne ! magnetic moment per atom
real*8,  allocatable, dimension(:) :: umx, umy, umz ! mignetization direction

save

end module




module param_atom
!-----------------------------------------------------------------------
! type declaration and initialization of variables for atoms
!-----------------------------------------------------------------------
implicit none


real*8,  allocatable, dimension(:) :: zatom    ! atomic number

logical :: lkbpp, lvand, lvfull, llocl         ! pseudopotential type
logical, allocatable, dimension(:) :: lkbppi   ! normconserving pp.
logical, allocatable, dimension(:) :: lvandi   ! ultrasoft pp.
logical, allocatable, dimension(:) :: lvflag   ! full ultrasoft pp.
logical, allocatable, dimension(:) :: llocli   ! local pp
logical :: lkbpp_r, lvand_r, lnlpp_r, llclpp_r, lpcc_r ! .true. = calculate in real space
logical :: lkbpp_g, lvand_g, lnlpp_g, llclpp_g, lpcc_g ! .true. = calculate in reciprocal space
logical, allocatable, dimension(:) :: lking    ! smoothing for nonlocal pp.
real*8,  allocatable, dimension(:) :: rking, gkgmax, gkgexct
logical, allocatable, dimension(:) :: llking   ! smoothing for local pp.
real*8,  allocatable, dimension(:) :: rlking, glkgmax, glkgexct

real*8 :: rctflc = 9.5d0   ! cutoff length for local pp.

logical :: lpcc                                ! partial core correction
logical, allocatable, dimension(:) :: lpcci    ! partial core correction for each atom
real*8,  allocatable, dimension(:) :: rpcc     ! radius for PCC
logical, allocatable, dimension(:) :: lpking   ! smoothing for PCC
real*8,  allocatable, dimension(:) :: rpking, gpkgmax, gpkgexct

integer, allocatable, dimension(:) :: nhk      ! No. of atoms
integer, allocatable, dimension(:) :: nhk1, nhk2
integer, allocatable, dimension(:) :: icscale  ! 1:scaled, 2:real coordinates
logical, allocatable, dimension(:) :: vrandom  ! .true.=random velocities

!---     ratm = real   coordinates for cluster calculations
!---            scaled coordinates for bulk    calculations
real*8,  allocatable, dimension(:,:) :: ratm   ! atomic coordinates
real*8,  allocatable, dimension(:,:) :: realatm ! atomic real coordinates 
                                                ! transferred from MD nodes
real*8,  allocatable, dimension(:,:) :: vatm    ! velocities transferred from MD nodes
real*8, dimension(3)                 :: disp_md ! displacement between QM & MD atoms
real*8,  allocatable, dimension(:,:) :: wrka   ! work array
real*8,  allocatable, dimension(:) :: watom    ! mass number
real*8,  allocatable, dimension(:) :: covrad   ! covalent radius


integer, parameter :: mxlx = 3

integer, parameter :: mulaox = 8
integer, allocatable, dimension(:) :: nmulao   ! No. of atomic orbitals
integer, allocatable, dimension(:,:) :: lmulao ! l for atomic orbital
real*8,  allocatable, dimension(:,:,:) :: rdecmp ! decomposition radius
real*8,  allocatable, dimension(:,:) :: rxmulk ! maximum radius
                                               ! for Mulliken analysis
logical, allocatable, dimension(:,:) :: lmulbsis ! how to generate atomic orbital

real*8,  allocatable, dimension(:) :: rsphexp  ! radius for Ylm expansion

real*8,  allocatable, dimension(:) :: radeda   ! radius for EDA

real*8,  allocatable, dimension(:) :: rintchg  ! radius for atomc charge

real*8,  allocatable, dimension(:) :: zv       ! valence
integer, allocatable, dimension(:) :: lmax     ! maximum angular momnetum (l)
integer, allocatable, dimension(:) :: lclno    ! l for local potential
logical, allocatable, dimension(:,:) :: lchk   ! for each l
integer, allocatable, dimension(:) :: lsomax   ! > 0 = relativistic data set


!-----for ultrasoft pseudopotentials
real*8,  allocatable, dimension(:,:) :: vdrcut ! cutoff length
real*8 :: rctmax = 0.d0                        ! max. of vdrcut
integer, allocatable, dimension(:,:) :: nrefe  ! the number of reference E.
integer :: mxref = 0                           ! max. of nrefe

!-----for the PAW method only
real*8  :: frac_rcomp = 1.3d0                  ! default value for frac_rcomp_it
real*8,  allocatable, dimension(:) :: frac_rcomp_it
                                               ! compensation-charge cutoff = max(r_c)/frac_rcomp_it

!-----for the mean-field Hubbard model (DFT+U)
logical, allocatable, dimension(:) :: lplusU_at  ! on/off DFT+U for each atom
integer, allocatable, dimension(:) :: plusU_l    ! angular momentum l, usually = 2
real*8,  allocatable, dimension(:) :: plusU_U    ! Hubbard parameter U
real*8,  allocatable, dimension(:) :: plusU_J    ! the screened exchange energy J

!-----for an empirical correction to non-local pp. (DFT+C)
logical, allocatable, dimension(:) :: lplusC_at  ! on/off DFT+U for each atom
real*8,  allocatable, dimension(:,:) :: plusC_r    ! cutoff radius
real*8,  allocatable, dimension(:,:) :: plusC_e    ! energy shift

save


end module




subroutine get_pwscale( pwscale_ )
!-----------------------------------------------------------------------
!  get nstep
!-----------------------------------------------------------------------
use param
implicit none
real*8  :: pwscale_

pwscale_ = pwscale

return
end subroutine




subroutine get_lnoncollinear( lnoncollinear_ )
!-----------------------------------------------------------------------
!  get nstep
!-----------------------------------------------------------------------
use param
implicit none
logical :: lnoncollinear_

lnoncollinear_ = lnoncollinear

return
end subroutine




subroutine param_atom_alloc( nfile, myid, nodes,  &
& alloc_mem, ntype, natom, natom_alloc )
!-----------------------------------------------------------------------
!     allocate memory for atoms
!-----------------------------------------------------------------------
use param_atom
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: ntype, natom, natom_alloc

!-----declare local variables
integer :: status
real*8  :: the_mem
integer :: it, i, j


!------allocate memory
allocate( zatom(ntype), lkbppi(ntype), lvandi(ntype),  &
& lvflag(ntype), llocli(ntype),  &
& lking(ntype),  rking(ntype),  gkgmax(ntype),  gkgexct(ntype),  &
& llking(ntype), rlking(ntype), glkgmax(ntype), glkgexct(ntype),  &
& lpcci(ntype), rpcc(ntype), lpking(ntype),  &
& rpking(ntype), gpkgmax(ntype), gpkgexct(ntype),  &
& nhk(ntype), nhk1(ntype), nhk2(ntype), icscale(ntype),  &
& vrandom(ntype), ratm(3,natom), realatm(3,natom_alloc), vatm(3,natom),  &
& wrka(3,natom), watom(ntype), covrad(ntype),  &
& nmulao(ntype), lmulao(mulaox,ntype), rdecmp(2,mulaox,ntype),  &
& rxmulk(mulaox,ntype), lmulbsis(mulaox,ntype), &
& rsphexp(ntype), radeda(ntype),  &
& rintchg(ntype),  &
& zv(ntype), lmax(ntype), lclno(ntype), lsomax(ntype), frac_rcomp_it(ntype), &
& lplusU_at(ntype), plusU_l(ntype), plusU_U(ntype), plusU_J(ntype), &
& lplusC_at(ntype), plusC_r(0:mxlx,ntype), plusC_e(0:mxlx,ntype), &
& stat=status )

the_mem =  &
&   8.d0 * ntype  &
& + 1.d0 * ntype*4  &
& + 1.d0 * ntype + 8.d0 * ntype*3  &
& + 1.d0 * ntype + 8.d0 * ntype*3  &
& + 1.d0 * ntype + 8.d0 * ntype  &
& + 1.d0 * ntype + 8.d0 * ntype*3  &
& + 4.d0 * ntype*4 + 1.d0 * ntype  &
& + 8.d0 * ( 3*natom*4 + ntype*2 )  &
& + 4.d0 * ( mulaox + 1 )*ntype  &
& + 8.d0 * ( 3*mulaox*ntype + ntype*5 )  &
& + 4.d0 * ntype*2  &
& + 1.d0 * ( ntype + mulaox*ntype )  &
& + 1.d0 * ntype + 4.d0 * ntype + 8.d0 * ntype*2

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'param_atom_alloc', .true. )


!-----set initial values
do it = 1, ntype
   zatom(it)  = 0.d0     ! atomic number
   lkbppi(it) = .true.   ! normconserving pp.
   lvandi(it) = .false.  ! ultrasoft pp.
   lvflag(it) = .false.  ! full ultrasoft pp.
   llocli(it) = .false.  ! local pp
   lking(it)  = .false.  ! smoothing for nonlocal pp.
   llking(it) = .false.  ! smoothing for local pp.
   lpcci(it)  = .false.  ! partial core correction
   rpcc(it)   = 0.d0     ! radius for PCC
   lpking(it) = .false.  ! smoothing for PCC
   nhk(it)    = 0          ! No. of atoms
   nhk1(it)   = 0
   nhk2(it)   = 0
   icscale(it) = 1       ! 1:scaled, 2:real coordinates
   vrandom(it) = .true.  ! .true.=random velocities
   watom(it)  = 0.d0     ! mass number
   covrad(it) = 0.d0     ! covalent radius
   rsphexp(it) = 0.d0    ! radius for Ylm expansion
   radeda(it) = 0.d0     ! radius for EDA
   rintchg(it) = 0.d0    ! radius to calculate atomic charge
   zv(it)     = 0.d0     ! valence
   lmax(it)   = 0        ! maximum angular momnetum (l)
   lclno(it)  = 0        ! l for local potential
   lsomax(it) = 0
   frac_rcomp_it(it) = frac_rcomp
   lplusU_at(it) = .false.
   plusU_l(it) = 2
   plusU_U(it) = 0.d0
   plusU_J(it) = 0.d0
   lplusC_at(it) = .false.
   plusC_r(0:mxlx,it) = 0.d0
   plusC_e(0:mxlx,it) = 0.d0
end do

do it = 1, ntype
   nmulao(it) = 4        ! No. of atomic orbitals
   do j = 1, nmulao(it)
      lmulao(j,it) = (j-1)*10 ! l for atomic orbitals
      rxmulk(j,it) = 10.d0    ! maximum radius
      rdecmp(1,j,it) = 0.d0 ! decomposition radius
      rdecmp(2,j,it) = 0.d0 ! decomposition radius
      lmulbsis(j,it) = .true.
   end do
end do

do j = 1, 3
do i = 1, natom
   ratm(j,i) = 0.d0      ! atomic coordinates
   vatm(j,i) = 0.d0      ! atomic velocities
   wrka(j,i) = 0.d0      ! work array
end do
end do


!-----allocate memory for work arrays in ppkb.f
call psvariables_alloc( nfile, myid, nodes,  &
& alloc_mem, ntype )


return
end subroutine




subroutine mxl_alloc( nfile, myid, nodes )
!-----------------------------------------------------------------------
!     allocate memory related to mxl
!-----------------------------------------------------------------------
use param
use param_atom
implicit none
integer :: nfile(*), myid, nodes

!-----declare local variables
integer :: it, l
integer :: status
real*8  :: the_mem


!------allocate memory
allocate( lchk(0:mxl,ntype), nrefe(0:mxl,ntype), stat=status )
if( lvand .and. status == 0 )  &
& allocate( vdrcut(0:mxl,ntype), stat=status )

the_mem = 1.d0 * (mxl+1)*ntype + 4.d0 * (mxl+1)*ntype
if( lvand )  &
& the_mem = the_mem + 8.d0 * (mxl+1)*ntype

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'mxl_alloc', .true. )


!-----allocate memory for work arrays in kbpp.f
call psvariables_mxl_alloc( nfile, myid, nodes )


if( lvand ) then
    do it = 1, ntype
    do l = 0, lmax(it)
       vdrcut(l,it) = 0.d0
       nrefe(l,it)  = 0
    end do
    end do
end if


return
end subroutine




subroutine set_lvsend(  nfile, myid, nodes, lvsend )
!-----------------------------------------------------------------------
!    set lvsend
!-----------------------------------------------------------------------
use param
implicit none
integer :: nfile(*), myid, nodes
logical :: lvsend

lvsend = ltddft_fssh .and. lfssh_vscale

return
end subroutine




module pwlda_variables
!-----------------------------------------------------------------------
! type declaration of shared variables in pwlda.f
!-----------------------------------------------------------------------
implicit none

integer :: myx, myy, myz
integer, dimension(6) :: nn
integer, dimension(3) :: myparity

logical :: lclust = .false.  ! .true. = cluster calculations
real*8  :: gamma             ! Ewald parameter
real*8  :: rccc2
real*8  :: dgalpha           ! convergence parameter in the double-grid method

!----- variables only for bulk calculations
logical :: lhfull    ! .false. = minimum-image calculation in Ewald method

real*8,  dimension(3,3) :: rba

integer :: nplw5, nplw, nplw2, nplw3, nplwc, nplw7, nplwcs
integer :: nplw5ex, nplwex, nplw2ex, nplw3ex, nplwcex, nplw7ex
integer :: kfft0d
integer :: kfft1,  kfft2,  kfft3,  kfft0
integer :: kfft1b, kfft2b, kfft3b, kfft0b
integer :: kmax1d, kmax2d, kmax3d
integer :: kmax1, kmax2, kmax3
integer :: kmax1cs, kmax2cs, kmax3cs
real*8  :: rvol
integer :: nbnod
integer :: natnod
integer :: npnod
integer :: npnod7
integer :: nultg

integer :: nstepCG = 0
integer :: nstepMD = 0
integer :: nstep = 0
integer :: jgcycl = 0

integer :: keypwv
integer :: keypcd
integer :: keypmg

real*8  :: cshdep

real*8  :: sume      ! QM energy

real*8, parameter  :: pi =  3.14159265358979323846d0

character(50):: fname_ion   ! data-file name for the next run
character(50):: fname_eig
character(50):: fname_eigk
character(50):: fname_tddft
character(50):: fname_cds
character(50):: fname_hrt
character(50):: fname_pcds
character(50):: fname_peig
character(50):: fname_peigk

character(50):: dname_wann  ! base-file name to dump Wannier func.
character(50):: dname_eig   ! base-file name to dump wave functions
character(50):: dname_cds   ! base-file name to dump wave functions
character(50):: dname_pot   ! base-file name to dump local potential
character(50):: dname_potav ! base-file name to dump local potential

real*8 :: entrpy = 0.d0
real*8, dimension(2) :: feneud = 0.d0
real*8 :: fermie = 0.d0

save


end module




module pwlda_pp
!-----------------------------------------------------------------------
! type declaration of variables for pseudopotentials
!-----------------------------------------------------------------------
implicit none


!-----variables for spherical harmonics
real*8, allocatable, dimension(:,:) :: ylmr, ylmi

!-----variables for error function
real*8, allocatable, dimension(:) :: tberf, tberfa
real*8 :: drcut
real*8, allocatable, dimension(:) :: tbeff, tbeffa
real*8 :: drcutf

!-----variables for local pseudopotential
real*8, allocatable, dimension(:,:) :: tablc, tablca
real*8, allocatable, dimension(:)   :: dltlc, rmxlc
real*8, allocatable, dimension(:,:) :: tbflc, tbflca
real*8, allocatable, dimension(:)   :: dltflc, rmxflc

integer :: nvlcl
real*8, allocatable, dimension(:) :: vlocli, xitgrd, rr

save


end module




module pwlda_atom
!-----------------------------------------------------------------------
! type declaration of variables for atoms and pseudopotentials
!-----------------------------------------------------------------------
implicit none


integer, allocatable, dimension(:) :: iatoit     ! reference to atomic species
integer, allocatable, dimension(:) :: iatmpt_nod ! atom index in each node
integer :: nion_nod                              ! the number of atoms in each node
real*8,  allocatable, dimension(:)  :: bufatm, bufatmr  ! buffer array
real*8, allocatable, dimension(:,:) :: pratm     ! 


!------- variables for energy, forces & stresses
real*8 :: ecolmb                    ! direct Coulomb energy
real*8, allocatable, dimension(:,:) :: floc     ! local potential force
real*8, allocatable, dimension(:,:) :: fnlc     ! nonlocal potential force
real*8, allocatable, dimension(:,:) :: fclm     ! direct Coulomb force
real*8, allocatable, dimension(:,:) :: fpcc     ! force from PCC
real*8, allocatable, dimension(:,:) :: fhfs     ! short-ranged HF force
real*8, allocatable, dimension(:,:) :: fefd     ! electric field force
real*8, allocatable, dimension(:,:) :: frc      ! total force
real*8, allocatable, dimension(:,:,:) :: prevr
real*8, dimension(3,3) :: strtot = 0.d0    ! total stress
real*8, dimension(3,3) :: strloc    ! local potential stress
real*8, dimension(3,3) :: strnlc    ! nonlocal potential stress
real*8, dimension(3,3) :: strclm    ! direct Coulomb stress
real*8, dimension(3,3) :: strgga = 0.d0   ! stress by density gradient in GGA
real*8, dimension(3,3) :: strpcc    ! partial core correction stress
real*8, dimension(3,3) :: strhfs    ! short-ranged HF stress
real*8, dimension(3,3) :: strefd    ! electric field stress

integer, allocatable, dimension(:) :: nhk1r, nhk2r       ! work array in ivelct
integer, allocatable, dimension(:) :: nhk1_nod, nhk2_nod ! for paralellization
integer, allocatable, dimension(:) :: nhk1_nat, nhk2_nat ! related to natnod
integer, allocatable, dimension(:) :: ioa, ioag, ioas, ioasg


!--- for EDA
real*8, allocatable, dimension(:) :: atom_ecolmb ! atomic direct Coulomb energy

save


end module




module pwlda_proc
!-----------------------------------------------------------------------
! type declaration of variables depending on the number of processors
!-----------------------------------------------------------------------
implicit none

integer,allocatable,dimension(:) :: nplcnt, npldsp ! for G decomp.
integer,allocatable,dimension(:) :: ncgcnt, ncgdsp ! for unifying wavefunctions
integer,allocatable,dimension(:) :: nbncnt, nbndsp ! for band decomp.
integer,allocatable,dimension(:) :: natcnt, natdsp ! for atom decomp.
integer,allocatable,dimension(:) :: lbncnt, lbndsp ! work area
integer,allocatable,dimension(:) :: mbncnt, mbndsp ! work area
integer,allocatable,dimension(:) :: idstnd, jdstnd ! for all-to-all commun.
integer,allocatable,dimension(:) :: nplcnt7, npldsp7 ! for G decomp.
integer,allocatable,dimension(:) :: ncgcnt7, ncgdsp7 ! for unifying wavefunctions
save

end module




module pwlda_pw
!-----------------------------------------------------------------------
! type declaration of variables for electronic bands and plane waves
!-----------------------------------------------------------------------
implicit none

real*8,  allocatable, dimension(:) :: eig     ! eigenvalues
real*8,  allocatable, dimension(:) :: occ     ! occupancies
real*8,  allocatable, dimension(:,:) :: eigtmp ! temporal eigenvalues
integer :: nocc                               ! the number of occupied states
real*8,  allocatable, dimension(:) :: dmtrxr, bijr, bm1r, bx0r
real*8,  allocatable, dimension(:) :: dmtrxi, biji, bm1i, bx0i
integer, allocatable, dimension(:) :: norder  ! order of states
integer :: nbnod1, nbnod2                     ! band index for each node
integer :: nbxxxx
integer :: natnod1, natnod2                     ! band index for each node
integer :: npnod1, npnod2
integer :: nspnod
integer :: npnnbnx
integer :: npnod71, npnod72
integer :: nspnod7

!-----work variables
integer, allocatable, dimension(:) :: nriter
real*8,  allocatable, dimension(:) :: prod, prodr, prods, betk, bzan
real*8,  allocatable, dimension(:) :: w1r, w1i, ew1
real*8,  allocatable, dimension(:,:) :: tmpr


!--- for the plane wave method -----------------------------------------
integer :: nfftk

real*8,  allocatable, dimension(:) :: recnrm
real*8,  allocatable, dimension(:) :: gx, gy, gz
integer, allocatable, dimension(:) :: nga, ngb, ngc
real*8,  allocatable, dimension(:) :: recnrmex ! for the double-grid method

real*8,  allocatable, dimension(:) :: dkgnrm
integer, allocatable, dimension(:) :: ijkgd
integer, allocatable, dimension(:) :: ijkg

!-----for noncollinear magnetism
!real*8,  allocatable, dimension(:) :: dkrecnrm
!real*8,  allocatable, dimension(:) :: dkgx, dkgy, dkgz
!integer, allocatable, dimension(:) :: ndkga, ndkgb, ndkgc
!integer, allocatable, dimension(:) :: ncijkg

!-----for ultrasoft pseudopotential
real*8,  allocatable, dimension(:) :: gboxx, gboxy, gboxz
integer, allocatable, dimension(:) :: ngboxa, ngboxb, ngboxc
real*8,  allocatable, dimension(:) :: abgvc3
integer, allocatable, dimension(:) :: ijkgb
integer :: kmax1b, kmax2b, kmax3b


#if CCREAL4
real*8,  allocatable, dimension(:) :: cgjr
real*4,  allocatable, dimension(:) :: cgjr4
real*4,  allocatable, dimension(:) :: gdcr4
real*4,  allocatable, dimension(:) :: bufcr4
#else
real*8,  allocatable, dimension(:) :: cgjr
#endif
real*8,  allocatable, dimension(:) :: rhcr, gdcr
real*8,  allocatable, dimension(:) :: bufcr
real*8,  allocatable, dimension(:) :: hcsr
real*8,  allocatable, dimension(:) :: sck
real*8,  allocatable, dimension(:) :: pgdcr

real*8,  allocatable, dimension(:) :: prcd, gnk, hnk
real*8,  allocatable, dimension(:) :: ekib
integer, allocatable, dimension(:) :: iod,  iods
integer, allocatable, dimension(:) :: iodg, iodsg
real*8,  allocatable, dimension(:)  :: fft3x, fft3y
complex*16, allocatable, dimension(:) :: fftwork
integer :: ntotfd
integer, allocatable, dimension(:) :: mftnod, mftdsp
integer, allocatable, dimension(:) :: mfd2ft
integer, allocatable, dimension(:) :: mshglb
real*8,  allocatable, dimension(:) :: glocal
real*8,  allocatable, dimension(:) :: apk
real*8,  allocatable, dimension(:) :: apki
real*8,  allocatable, dimension(:) :: eigr
real*8,  allocatable, dimension(:) :: eigi
real*8,  allocatable, dimension(:) :: scwr
real*8,  allocatable, dimension(:) :: scwi
integer, allocatable, dimension(:) :: mftwrk

!-----for FFT mesh on coarse grid
integer :: ntotfd_c
real*8,  dimension(3)   :: rdel_c
real*8,  dimension(3,3) :: rdelg_c
real*8  :: rdelv_c
integer, dimension(3)   :: nd1vks_c
integer, allocatable, dimension(:) :: mfd2ft_c
integer, allocatable, dimension(:) :: mshglb_c
integer, allocatable, dimension(:) :: mshgnx, mshgny, mshgnz

!--- for charge density
real*8,  allocatable, dimension(:)   :: rhgr, rinr, thrhgr
real*8,  allocatable, dimension(:)   :: rhgrsoft

real*8,  allocatable, dimension(:)   :: pdbuf

integer :: nprvmx
real*8,  allocatable, dimension(:) :: prveig

integer :: ncprvmx
real*8,  allocatable, dimension(:) :: prvrho

!-----sine, cosine functions for atoms
real*8, allocatable, dimension(:) :: ycos, ysin
real*8, allocatable, dimension(:) :: DCQ1, DCQ2, DCQ3
real*8, allocatable, dimension(:) :: DSQ1, DSQ2, DSQ3


!--- for spin-polarized calculations -----------------------------------
integer :: nspnmx            ! = 2 for spin-polarized calculations
integer :: nspnmx2           ! = 2 for spin-polarized calculations
                             ! = 4 for noncollinear magnetism
logical :: lcgjsv = .true.

real*8,  allocatable, dimension(:) :: gdcrsv
real*8,  allocatable, dimension(:) :: rhcrsv
real*8,  allocatable, dimension(:) :: hcsrsv
real*8,  allocatable, dimension(:) :: pgdcrsv

real*8,  allocatable, dimension(:) :: wegud
real*8,  allocatable, dimension(:) :: egvlud
integer, allocatable, dimension(:) :: nsordr
real*8,  allocatable, dimension(:) :: hlocal
real*8,  allocatable, dimension(:) :: svtrxr
real*8,  allocatable, dimension(:) :: svtrxi

!--- for linear-response TDDFT
real*8,  allocatable, dimension(:,:) :: rhgrph

!--- for hybrid functionlas
real*8,  allocatable, dimension(:) :: cgjr_m, cgjr_buf, hfrhcr_m
integer, allocatable, dimension(:) :: iod_m, iod_buf
real*8,  allocatable, dimension(:) :: hfchg,  hfv,  hfwav
real*8,  allocatable, dimension(:) :: hfchgi, hfvi, hfwavi
real*8,  allocatable, dimension(:) :: hfrhgr
real*8,  allocatable, dimension(:,:) :: hfrhcr
!real*8,  allocatable, dimension(:) :: thehfrhcr
real*8,  allocatable, dimension(:) :: hfv_c, hfvi_c
real*8,  allocatable, dimension(:) :: hfrecnrm, hfstrrecnrm
real*8,  allocatable, dimension(:) :: hfeiki, hfeikj, hfeikj_tmp
real*8,  allocatable, dimension(:) :: hfcos, hfsin
integer :: khffftd
integer :: nplw5ex1

integer :: nknodx

save


end module




module pwlda_grid
!-----------------------------------------------------------------------
! type declaration of variables for multigrid
!-----------------------------------------------------------------------
implicit none


integer :: nnd2v

integer :: mshmx, mshmy, mshmz
integer :: lnmx
integer :: ndata
integer :: msnodx, msnody, msnodz
integer :: noddatx
integer :: muldatx

integer :: msndx1, msndy1, msndz1
integer :: msndx2, msndy2, msndz2
integer :: nodext
integer :: nodyzx, nodzxx, nodxyx
integer :: mulyzx, mulzxx, mulxyx
integer :: mulext
integer :: nodexe
integer :: nodexf
integer :: mulexf
integer :: mshndv


integer ::  nmm
integer, allocatable, dimension(:) :: nmmcnt, nmmdsp
integer, allocatable, dimension(:) :: nmmrc1, nmmrc2
integer, allocatable, dimension(:) :: meshx1, meshx
integer, allocatable, dimension(:) :: meshy1, meshy
integer, allocatable, dimension(:) :: meshz1, meshz

integer, allocatable, dimension(:) :: mshnod
integer, allocatable, dimension(:) :: mshnew
integer, allocatable, dimension(:) :: mshxyz
integer, allocatable, dimension(:) :: mulpit, mulx1, mulx2,  &
&                               muly1, muly2, mulz1, mulz2
integer, allocatable, dimension(:) :: mshnx, mshny, mshnz
integer, allocatable, dimension(:) :: mulpms
integer, allocatable, dimension(:) :: mshx1, mshy1, mshz1,  &
&                                     mshx,  mshy,  mshz
integer, allocatable, dimension(:) :: mulpem

integer :: mulnpx
integer :: mulnpy
integer :: mulnpz
integer, allocatable, dimension(:) :: ismx
integer, allocatable, dimension(:) :: ismy
integer, allocatable, dimension(:) :: ismz
integer, allocatable, dimension(:) :: ismshx
integer, allocatable, dimension(:) :: ismshy
integer, allocatable, dimension(:) :: ismshz
integer, allocatable, dimension(:) :: mulpsx, mulpsy, mulpsz,  &
&                                     mulyzs, mulzxs, mulxys
integer, allocatable, dimension(:) :: irmx
integer, allocatable, dimension(:) :: irmy
integer, allocatable, dimension(:) :: irmz
integer, allocatable, dimension(:) :: irmshx
integer, allocatable, dimension(:) :: irmshy
integer, allocatable, dimension(:) :: irmshz

integer :: nodyzh
integer :: nodzxh
integer :: nodxyh
integer, allocatable, dimension(:) :: jsmshx
integer, allocatable, dimension(:) :: jsmshy
integer, allocatable, dimension(:) :: jsmshz
integer, allocatable, dimension(:) :: mulprx, mulpry, mulprz,  &
&                                     mulyzr, mulzxr, mulxyr
integer, allocatable, dimension(:) :: ismfpx
integer, allocatable, dimension(:) :: ismfpy
integer, allocatable, dimension(:) :: ismfpz
integer, allocatable, dimension(:) :: irmfpx
integer, allocatable, dimension(:) :: irmfpy
integer, allocatable, dimension(:) :: irmfpz
integer, allocatable, dimension(:) :: idbuf, idbufr

integer :: msndx4
integer :: msndy4
integer :: msndz4
real*8,  allocatable, dimension(:)   :: fmshx
integer, allocatable, dimension(:,:) :: ncomct
integer, allocatable, dimension(:)   :: mulnd2

integer :: ndv2x3
real*8, allocatable,  dimension(:,:) :: amcdv2
real*8, allocatable,  dimension(:,:) :: cdv1, cdv2
#ifdef VECTOR
real*8, allocatable, dimension(:) :: xxx
#endif

integer, allocatable, dimension(:) :: mulgmx
integer, allocatable, dimension(:) :: mulgmy
integer, allocatable, dimension(:) :: mulgmz
!-----for serial calculations
integer :: muls
integer, dimension(3) :: nd1vs
integer, allocatable, dimension(:) :: kulnd2

!-----for second grid
integer :: lmuls
integer, dimension(3) :: lnd1vs
integer, allocatable, dimension(:) :: nulnd2
integer, allocatable, dimension(:) :: lulnd2


!------- variables for potential
real*8, allocatable, dimension(:) :: rho       ! electron density in real space
real*8, allocatable, dimension(:) :: vhar      ! hartree potential
real*8, allocatable, dimension(:) :: vext      ! local pp potential
real*8, allocatable, dimension(:) :: vexc      ! exchange-correlation potential
real*8, allocatable, dimension(:) :: eeexc     ! exchange-correlation energy
real*8, allocatable, dimension(:) :: hdiag     ! = vhar + vext + vexc
real*8, allocatable, dimension(:) :: rhocore   ! core charge
!---for non-self consistent (NSC) forces
real*8, allocatable, dimension(:) :: vhar_out  ! hartree potential
real*8, allocatable, dimension(:) :: delrho    ! rho_out - rho_in
real*8, allocatable, dimension(:) :: xexc      ! 
real*8, allocatable, dimension(:) :: ddelrho   ! 
real*8, allocatable, dimension(:) :: vwell     ! well potential
real*8, allocatable, dimension(:) :: vefield   ! uniform electric field potential


!------- variables for conjugate gradient method
real*8, allocatable, dimension(:) :: x
real*8, allocatable, dimension(:) :: rk, vk
real*8, allocatable, dimension(:) :: dbuf, dbufr
real*8, allocatable, dimension(:) :: xx
real*8, allocatable, dimension(:) :: tmpk
real*8, allocatable, dimension(:) :: tmpl, tmpm, tmpn

integer, allocatable, dimension(:) :: nrstrct
real*8,  allocatable, dimension(:) :: prevres


!--- for shape denpendent terms due to dipole moment of the system
real*8, dimension(3) :: dpion, dprho
real*8, allocatable, dimension(:) :: vlshdp, vhshdp

!--- for spin-polarized calculations
real*8, allocatable, dimension(:) :: rhoud
real*8, allocatable, dimension(:) :: vlocud

!--- for linear-response TDDFT
real*8,  allocatable, dimension(:,:) :: rhoph
real*8,  allocatable, dimension(:) :: rhophe
real*8,  allocatable, dimension(:) :: vexcph      ! exchange-correlation potential
real*8,  allocatable, dimension(:) :: eeexcph     ! exchange-correlation energy
real*8,  allocatable, dimension(:) :: rhoudph
real*8,  allocatable, dimension(:) :: vlocudph
real*8,  allocatable, dimension(:) :: hdiagph

!-----for noncollinear magnetism
!real*8,  allocatable, dimension(:) :: rhomx, rhomy, rhomz
!real*8,  allocatable, dimension(:) :: prvrhom

save


end module




module ncmagne_variables
!-----------------------------------------------------------------------
! type declaration of variables for noncollinear magnetism
!-----------------------------------------------------------------------
implicit none

!-----for electronic bands and plane waves
real*8,  allocatable, dimension(:) :: dkrecnrm
real*8,  allocatable, dimension(:) :: dkgx, dkgy, dkgz
integer, allocatable, dimension(:) :: ndkga, ndkgb, ndkgc
integer, allocatable, dimension(:) :: ncijkg
real*8,  allocatable, dimension(:) :: dkdkgnrm

!-----for magnetic moment
real*8,  allocatable, dimension(:) :: rhomx, rhomy, rhomz
!real*8,  allocatable, dimension(:) :: prvrhom

real*8,  allocatable, dimension(:) :: rhomxini, rhomyini, rhomzini
real*8,  allocatable, dimension(:) :: eig2r, eig2i

real*8,  allocatable, dimension(:,:,:) :: wvratio, wvmagxy

real*8,  allocatable, dimension(:,:) :: totmag  ! total magnetic moment of each atom
real*8,  allocatable, dimension(:) :: fugomag   ! up or down of each atom's spin

logical, allocatable, dimension(:) :: leachatmref, leachclearamag

save

end module




subroutine get_nstep( nstep_ )
!-----------------------------------------------------------------------
!  get nstep
!-----------------------------------------------------------------------
use pwlda_variables
implicit none
integer :: nstep_

nstep_ = nstep

return
end subroutine




subroutine set_rciprl( nfile, myid, nodes,  &
& nplw5, nplw, nplw2, nplwc,  &
& gs1, gs2, gs3, gss, ktmp1, ktmp2, ktmp3, ix, ifftk )
!-----------------------------------------------------------------------
!     allocate memory for reciprocal vectors
!-----------------------------------------------------------------------
use param
use pwlda_pw
use ncmagne_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: nplw5, nplw, nplw2, nplwc
integer :: ifftk
integer, dimension(ifftk) :: ktmp1, ktmp2, ktmp3
real*8,  dimension(ifftk) :: gs1, gs2, gs3, gss
integer, dimension(ifftk) :: ix

!-----declare local variables
integer :: status
real*8  :: the_mem
integer :: i, j, ig, ig0, nplw0, igr, igi
integer :: ngenh = 0
integer :: lplh  = 0
integer :: lpl2h = 0
integer :: lplch = 0
save ngenh, lplh, lpl2h, lplch


if( nplw5 > ngenh .or.  &
&   nplw  > lplh  .or.  &
&   nplw2 > lpl2h .or.  &
&   nplwc > lplch     ) then

    !-----if already allocated, deallocate arrays
    if( allocated(nga) ) then

        the_mem =  &
& + 4.d0 * ( size(nga) + size(ngb) + size(ngc) )  &
& + 8.d0 * ( size(gx) + size(gy) + size(gz) + size(recnrm) + size(dkgnrm) )  &
& + 4.d0 * ( size(ijkgd) + size(ijkg) )  &
& + 8.d0 * ( size(recnrmex) )

        deallocate( nga, ngb, ngc, gx, gy, gz, recnrm,  &
& dkgnrm, ijkgd, ijkg, recnrmex,  &
& stat=status )

        if( status == 0 .and. allocated(hfrecnrm) ) then
            the_mem = the_mem  &
&                   + 8.d0 * ( size(hfrecnrm) + size(hfstrrecnrm) )
            deallocate( hfrecnrm, hfstrrecnrm, stat=status )
        end if

        if( status == 0 .and. allocated(dkgx) ) then
            the_mem = the_mem  &
& + 8.d0 * ( size(dkgx) + size(dkgy) + size(dkgz) + size(dkrecnrm) + size(dkdkgnrm) )  &
& + 4.d0 * ( size(ndkga) + size(ndkgb) + size(ndkgc) + size(ncijkg) )
            deallocate( dkgx, dkgy, dkgz, dkrecnrm, ndkga, ndkgb, ndkgc,  &
& ncijkg, dkdkgnrm, stat=status )
        end if

        !------error trap
        call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'set_rciprl', .true. )
    end if

    ngenh = nplw5 * pwscale
    lplh  = nplw  * pwscale
    lpl2h = nplw2 * pwscale
    lplch = nplwc * pwscale

    !------allocate memory
    allocate( nga(0:ngenh), ngb(0:ngenh), ngc(0:ngenh),  &
& gx(0:ngenh), gy(0:ngenh), gz(0:ngenh), recnrm(0:ngenh),  &
& dkgnrm((lplch+1)*2),  &
& ijkgd(2*ngenh+1), ijkg(2*lpl2h+1), recnrmex(0:ngenh),  &
& stat=status )

    the_mem =  &
& + 4.d0 * ( size(nga) + size(ngb) + size(ngc) )  &
& + 8.d0 * ( size(gx) + size(gy) + size(gz) + size(recnrm) + size(dkgnrm) )  &
& + 4.d0 * ( size(ijkgd) + size(ijkg) )  &
& + 8.d0 * ( size(recnrmex) )

    if( status == 0 .and. jhybrid /= 0 ) then
        if( lgamma ) then
            allocate( hfrecnrm(0:ngenh), hfstrrecnrm(0:ngenh), stat=status )
        else
            allocate( hfrecnrm(ngenh*2+1), hfstrrecnrm(ngenh*2+1), stat=status )
        end if
        the_mem = the_mem  &
&               + 8.d0 * ( size(hfrecnrm) + size(hfstrrecnrm) )
    end if

    if( status == 0 ) then
        if( lnoncollinear ) then
            allocate( dkgx(0:lplh), dkgy(0:lplh), dkgz(0:lplh), dkrecnrm(0:lplh),  &
& ndkga(0:lplh), ndkgb(0:lplh), ndkgc(0:lplh), ncijkg(0:lplh),  &
& dkdkgnrm((lplh+1)*2), stat=status )
        else
            allocate( dkgx(1), dkgy(1), dkgz(1), dkrecnrm(1),  &
& ndkga(1), ndkgb(1), ndkgc(1), ncijkg(1), dkdkgnrm(1), stat=status )
        end if
        the_mem = the_mem  &
& + 8.d0 * ( size(dkgx) + size(dkgy) + size(dkgz) + size(dkrecnrm) + size(dkdkgnrm) )  &
& + 4.d0 * ( size(ndkga) + size(ndkgb) + size(ndkgc) + size(ncijkg) )
    end if

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'set_rciprl', .true. )

end if


   nga(0) = 0
   ngb(0) = 0
   ngc(0) = 0
    gx(0) = 0.d0
    gy(0) = 0.d0
    gz(0) = 0.d0
recnrm(0) = 0.d0
do i = 1, nplw5
   j = ix(i)
   nga(i) = ktmp1(j)
   ngb(i) = ktmp2(j)
   ngc(i) = ktmp3(j)
   gx(i)  = gs1(j)
   gy(i)  = gs2(j)
   gz(i)  = gs3(j)
   recnrm(i) = gss(j)
end do


do i = 0, nplwc
   igr = 2*i + 1
   igi = igr + 1
   dkgnrm(igr) = recnrm(i)
   dkgnrm(igi) = recnrm(i)
end do


if( .not.lnoncollinear ) then
    nplw0 = nplw
!    do ig = 0, nplw
!       dkgx(ig) = gx(ig)
!       dkgy(ig) = gy(ig)
!       dkgz(ig) = gz(ig)
!   dkrecnrm(ig) = recnrm(ig)
!    end do
else
    !-----noncollinear magnetism
    nplw0 = nplw/2

    ig = 0
     dkgx(ig) = gx(ig)
     dkgy(ig) = gy(ig)
     dkgz(ig) = gz(ig)
    ndkga(ig) = nga(ig)
    ndkgb(ig) = ngb(ig)
    ndkgc(ig) = ngc(ig)
dkrecnrm(ig) = recnrm(ig)
    do ig = 1, nplw, 2
       ig0 = (ig+1)/2
       dkgx(ig) = gx(ig0)
       dkgy(ig) = gy(ig0)
       dkgz(ig) = gz(ig0)
      ndkga(ig) = nga(ig0)
      ndkgb(ig) = ngb(ig0)
      ndkgc(ig) = ngc(ig0)
   dkrecnrm(ig) = recnrm(ig0)
    end do
    do ig = 2, nplw, 2
       ig0 = (ig+1)/2
       dkgx(ig) = -gx(ig0)
       dkgy(ig) = -gy(ig0)
       dkgz(ig) = -gz(ig0)
      ndkga(ig) = -nga(ig0)
      ndkgb(ig) = -ngb(ig0)
      ndkgc(ig) = -ngc(ig0)
   dkrecnrm(ig) = recnrm(ig0)
    end do

    do i = 0, nplw
       igr = 2*i + 1
       igi = igr + 1
       dkdkgnrm(igr) = dkrecnrm(i)
       dkdkgnrm(igi) = dkrecnrm(i)
    end do
end if


!------memory allocation for k-point sampling
!call set_rciprlk( nfile, myid, nodes,  &
!& nplw0, pwscale, jhybrid, mxkpdupnm )


!------memory allocation for TDDFT
!call set_rciprl_tddft( nfile, myid, nodes, nplw, pwscale )


return
end subroutine




subroutine nplw3_alloc( nfile, myid, nodes, nplw3 )
!-----------------------------------------------------------------------
!     allocate memory for variables related to nplw3
!-----------------------------------------------------------------------
use param
use pwlda_pw
implicit none
integer :: nfile(*), myid, nodes
integer :: nplw3

!-----declare local variables
integer :: status
real*8  :: the_mem
integer :: lpl3h = -1
save lpl3h


if( nplw3 > lpl3h ) then

    !-----if already allocated, deallocate arrays
    if( allocated(gboxx) ) then
        deallocate( gboxx, gboxy, gboxz,  &
& ngboxa, ngboxb, ngboxc, abgvc3, ijkgb,  &
& stat=status )

        the_mem =  &
& + 8.d0 * ( lpl3h + 1 )*3  &
& + 4.d0 * ( lpl3h + 1 )*3  &
& + 8.d0 * ( lpl3h + 1 )  &
& + 4.d0 * ( 2*lpl3h + 1 )

        !------error trap
        call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'nplw3_alloc', .true. )
    end if

    !------allocate memory
    lpl3h = nplw3 * pwscale
    allocate( gboxx(0:lpl3h), gboxy(0:lpl3h), gboxz(0:lpl3h),  &
& ngboxa(0:lpl3h), ngboxb(0:lpl3h), ngboxc(0:lpl3h),  &
& abgvc3(0:lpl3h), ijkgb(2*lpl3h+1),  &
& stat=status )

    the_mem =  &
& + 8.d0 * ( lpl3h + 1 )*3  &
& + 4.d0 * ( lpl3h + 1 )*3  &
& + 8.d0 * ( lpl3h + 1 )  &
& + 4.d0 * ( 2*lpl3h + 1 )

    !------error trap
    call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'nplw3_alloc', .true. )

end if


return
end subroutine




subroutine get_nplw3( nfile, myid, nodes,  &
& nplw3, nplw5, crsbox, ecutdens )
!-----------------------------------------------------------------------
!     get nplw3
!-----------------------------------------------------------------------
use pwlda_pw
implicit none
integer :: nfile(*), myid, nodes
integer :: nplw3, nplw5
real*8,  dimension(3,3) :: crsbox
real*8  :: ecutdens

!-----declare local variables
integer :: ig
real*8  :: g1, g2, g3, ggs


nplw3 = 0
do ig = 1, nplw5
   G1 = crsbox(1,1)*dble(nga(ig)) + crsbox(1,2)*dble(ngb(ig))  &
&     + crsbox(1,3)*dble(ngc(ig))
   G2 = crsbox(2,1)*dble(nga(ig)) + crsbox(2,2)*dble(ngb(ig))  &
&     + crsbox(2,3)*dble(ngc(ig))
   G3 = crsbox(3,1)*dble(nga(ig)) + crsbox(3,2)*dble(ngb(ig))  &
&     + crsbox(3,3)*dble(ngc(ig))
   ggs = G1*G1 + G2*G2 + G3*G3
   if( ggs <= ecutdens  ) then
       nplw3 = nplw3 + 1
   end if
end do


return
end subroutine




subroutine pwlda_pp_alloc( nfile, myid, nodes,  &
& alloc_mem, ntype, mxl, mx1, mx1loc, lpcc )
!-----------------------------------------------------------------------
!     allocate memory for variables for atoms and pseudopotentials
!-----------------------------------------------------------------------
use pwlda_pp
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: ntype
integer :: mxl
integer :: mx1, mx1loc
logical :: lpcc

!-----declare local variables
integer :: status
real*8  :: the_mem


nvlcl = max( mx1, mx1loc )
!------allocate memory
allocate( ylmr(-mxl:mxl,0:mxl), ylmi(-mxl:mxl,0:mxl),  &
& tberf(0:mx1loc), tberfa(0:mx1loc),  &
& tbeff(0:mx1loc), tbeffa(0:mx1loc),  &
& tablc(0:mx1loc,ntype), tablca(0:mx1loc,ntype), dltlc(ntype),  &
& rmxlc(ntype), tbflc(0:mx1loc,ntype), tbflca(0:mx1loc,ntype),  &
& dltflc(ntype), rmxflc(ntype),  &
& vlocli(0:nvlcl), xitgrd(0:nvlcl), rr(0:nvlcl),  &
& stat=status )

the_mem =  &
&  8.d0 * ( (2*mxl+1)*(mxl+1)*2 + (mx1loc+1)*4  &
&         + (mx1loc+1)*ntype*2 + ntype*2  &
&         + (mx1loc+1)*ntype*2 + ntype*2 + (nvlcl+1)*3 )

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'pwlda_pp_alloc', .true. )


if( lpcc ) call pcc_variables_alloc( nfile, myid, nodes,  &
& alloc_mem, ntype, mx1 )


return
end subroutine




subroutine pwlda_atom_alloc( nfile, myid, nodes )
!-----------------------------------------------------------------------
!     allocate memory for variables for atoms and pseudopotentials
!-----------------------------------------------------------------------
use param
use param_atom
use pwlda_variables
use pwlda_proc
use pwlda_pw
use pwlda_atom
implicit none
integer :: nfile(*), myid, nodes

!-----declare local variables
integer :: status
real*8  :: the_mem
integer :: natmnod
integer :: i, j


!natom = nhk2(ntype)
natmnod = ( natom + nodes - 1 )/nodes
!nylmmx = mxl*(mxl+2) + 1

!------allocate memory
allocate( iatoit(natom),  &
& bufatm(3*natom), bufatmr(3*natom), pratm(3,natom),  &
& floc(3,natom), fnlc(3,natom), fclm(3,natom), fpcc(3,natom),  &
& frc(3,natom), prevr(3,natom,3),  &
& nhk1r(ntype), nhk2r(ntype), nhk1_nod(ntype), nhk2_nod(ntype),  &
& nhk1_nat(ntype), nhk2_nat(ntype),  &
& ioa(natmnod), ioag(natom), ioas(natmnod), ioasg(natom),  &
& atom_ecolmb(natom),  &
& stat=status )

the_mem =  &
&  4.d0 * ( natom )  &
& + 8.d0 * ( 3*natom*3 )  &
& + 8.d0 * ( 3*natom*5 + 3*natom*3 )  &
& + 4.d0 * ( ntype*6 + natmnod*4 ) &
& + 8.d0 * ( natom )

if( status == 0 .and. jhybrid /= 0 ) then
    allocate( fhfs(3,natom), stat=status )

    the_mem = the_mem + 8.d0 * size(fhfs)
end if

if( status == 0 .and. lefield ) then !lefield_islts ) then
    allocate( fefd(3,natom), stat=status )

    the_mem = the_mem + 8.d0 * size(fefd)
end if

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'pwlda_atom_alloc', .true. )


if( lkbpp_r ) then
    !-----allocate memory for variables for nonlocal pp
    call nlkbppr_variables_alloc( nfile, myid, nodes,  &
& alloc_mem, natom, mxl, nylmmx, nbnod, nspnmx )
end if


if( lkbpp_g .or. lkbpp_r ) then
    !-----allocate memory for variables for nonlocal pp
    call nlkbpp_variables_alloc( nfile, myid, nodes,  &
& alloc_mem, natom, mxl, nylmmx, nbnod, nspnmx )

    !-----first-order SO effects by purturbation calculation
!    call nlkbppSOp_variables_alloc( nfile, myid, nodes,  &
!& alloc_mem, natom, nband, nbnod )
end if


!-----clear variables
do i = 1, 3
do j = 1, natom
   prevr(1,j,i) = 0.d0
   prevr(2,j,i) = 0.d0
   prevr(3,j,i) = 0.d0
end do
end do


return
end subroutine




subroutine iatmpt_nod_alloc( nfile, myid, nodes,  &
& alloc_mem, ntype, nhk1, nhk2 )
!-----------------------------------------------------------------------
!     allocate memory for variables for atoms and pseudopotentials
!-----------------------------------------------------------------------
use pwlda_atom
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: ntype
integer, dimension(ntype) :: nhk1, nhk2

!-----declare local variables
integer :: status
real*8  :: the_mem
integer :: i, it
integer :: nion_nodx, nion_nod_alloc = -1
save nion_nod_alloc


nion_nodx = nion_nod
call gimax(nion_nodx)

if( nion_nodx > nion_nod_alloc ) then
    if( allocated(iatmpt_nod) ) then
        the_mem = 4.d0 * ( nion_nod_alloc )
        deallocate( iatmpt_nod, stat=status )
        !------error trap
        call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'iatmpt_nod_alloc', .true. )
    end if

    nion_nod_alloc = max( nion_nodx, 1 )
    !------allocate memory
    allocate( iatmpt_nod(nion_nod_alloc), stat=status )

    the_mem = 4.d0 * ( nion_nod_alloc )

    !------error trap
    call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'iatmpt_nod_alloc', .true. )
end if


do i = 1, nion_nod
   iatmpt_nod(i) = iatoit(i)
end do

do it = 1, ntype
do i = nhk1(it), nhk2(it)
   iatoit(i) = it
end do
end do


return
end subroutine




subroutine pwlda_proc_alloc( nfile, myid, nodes, alloc_mem )
!-----------------------------------------------------------------------
!     allocate memory for processor-dependent variables
!-----------------------------------------------------------------------
use pwlda_proc
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem

!-----declare local variables
integer :: status
real*8  :: the_mem
integer :: nproc


nproc = nodes
!------allocate memory
allocate( nplcnt(nproc), npldsp(nproc), ncgcnt(nproc),  &
& ncgdsp(nproc), nbncnt(nproc), nbndsp(nproc), lbncnt(nproc),  &
& lbndsp(nproc), mbncnt(nproc), mbndsp(nproc), idstnd(nproc),  &
& natcnt(nproc), natdsp(nproc), jdstnd(nproc),  &
& nplcnt7(nproc), npldsp7(nproc), ncgcnt7(nproc), ncgdsp7(nproc),  &
& stat=status )

the_mem = 4.d0 * nproc * 16

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'pwlda_proc_alloc', .true. )


return
end subroutine




subroutine pwlda_iod_alloc( nfile, myid, nodes,  &
& alloc_mem, nband )
!-----------------------------------------------------------------------
!     allocate memory for iod, iodg, iods, iodsg
!-----------------------------------------------------------------------
use pwlda_pw
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: nband

!-----declare local variables
integer :: nbnodx
integer :: status
real*8  :: the_mem


nbnodx = ( nband + nodes - 1 )/nodes

!------allocate memory
allocate( iod(nbnodx), iods(nbnodx), iodg(nband), iodsg(nband),  &
& stat=status )

the_mem = 4.d0 * ( nbnodx + nband ) * 2

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'pwlda_iod_alloc', .true. )


return
end subroutine




subroutine set_nprvmx( ichest, ihest, lxlbomd, kxlbomd )
!-----------------------------------------------------------------------
!     allocate memory for the plane wave method
!-----------------------------------------------------------------------
use pwlda_pw
implicit none
integer :: ichest, ihest
logical :: lxlbomd
integer :: kxlbomd

nprvmx  = max( ihest,  1 )
ncprvmx = max( ichest, 1 )

if( lxlbomd ) then
    if( kxlbomd == 0 ) kxlbomd = ihest
    nprvmx  = max( nprvmx,  kxlbomd + 1 )
    ncprvmx = max( ncprvmx, kxlbomd + 1 )

!    call set_nprvmx_in_nlpaw( kxlbomd )
end if

return
end subroutine




subroutine pwlda_alloc( nfile, myid, nodes )
!-----------------------------------------------------------------------
!     allocate memory for the plane wave method
!-----------------------------------------------------------------------
use param
use param_atom
use pwlda_variables
use pwlda_pw
use ncmagne_variables
implicit none
integer :: nfile(*), myid, nodes

!-----declare local variables
integer :: mshmx, mshmy, mshmz
integer :: msnodx, msnody, msnodz
integer :: noddatx
integer :: lpnbndx, nband2, lpnbndx2
integer :: status, ierror
real*8  :: the_mem
integer :: i
integer :: npnodx
integer :: nbnodx
integer :: nfftkd, nfftke  ! array size for FFT variables
integer :: lplnfftk   = 0
integer :: lplnfftkd  = 0
integer :: lplnfftke  = 0
integer :: lplnoddatx = 0
integer :: ngenhex    = 0    ! nplw5ex
integer :: lplhex0    = 0    ! nplwex
integer :: lplchex    = 0    ! nplwcex
integer :: lpl7hex0   = 0    ! nplw7ex
integer :: lplnpnod0  = 0    ! npnod
integer :: lplnpnodx  = 0    ! npnodx
integer :: lplhex, lpl7hex, lplnpnod
logical :: licall = .true.
logical :: lgamma_
integer :: nkpnt_ !, nknodx
real*8  :: pkscale
!logical :: ltddft, ltddft_fssh, lnoncollinear
integer :: nnspin !, nspnmx2
save lplnfftk, lplnfftkd, lplnfftke, lplnoddatx, ngenhex, lplhex,  &
& lpl7hex, lplchex, lplnpnod, lplnpnodx, licall, lpnbndx, &
 lplhex0, lpl7hex0, lplnpnod0


mshmx = nd1vks(1)
mshmy = nd1vks(2)
mshmz = nd1vks(3)
nfftk = mshmx*mshmy*mshmz

nfftkd = nfftk
nfftke = (mshmx+2)*mshmy*mshmz
#ifdef ASLSX
nfftke = (mshmx+4)*(mshmy+1)*(mshmz+1)
nfftke = max( nfftke + (mshmx+1)/2 + mshmy*mshmz,  &
&             (mshmx+1)*(mshmy+1)*(mshmz+1) + mshmx + mshmy + mshmz + 1 )
nfftkd = (mshmx+1)*(mshmy+1)*(mshmz+1)
#endif
!-----check array size for FFT variables
!if( ldouble_grid_recip )  &
! call check_fftsize( nfile, myid, nodes, nfftkd, nfftke )

msnodx = ( mshmx + npx - 1 )/npx
msnody = ( mshmy + npy - 1 )/npy
msnodz = ( mshmz + npz - 1 )/npz
noddatx = msnodx*msnody*msnodz


#if PCHOLESKY
nband2 = nbnod*nband
#else
nband2 = nband*nband
#endif

npnodx = npnod
nbnodx = nbnod
call gimax(npnodx)
call gimax(nbnodx)


!---get k-point sampling information
call get_kpsi( lgamma_, nkpnt_, pkscale )
call get_nknodx( nknodx )

!---get ltddft
!call get_ltddft( ltddft, ltddft_fssh )

!lnoncollinear = ncscale == 2

if( lspin ) then
    nnspin = 2
  else
    nnspin = 1
end if

if( licall ) then

    !------allocate memory
    allocate( eig(nband*nkpnt), eigtmp(nband*nkpnt,nnspin), occ(nband*nkpnt),  &
& dmtrxr(nband2*nknodx), bijr(nband2),  &
& bm1r(nband2), bx0r(nband2), norder(nband*nkpnt), nriter(nband),  &
& prod(nband), prodr(nband), prods(nband), betk(nband),  &
& bzan(nband), w1r(nband), w1i(nband), ew1(nband),  &
& tmpr(iblock,nband), ekib(nbnod*nnspin*nkpnt),  &
& mftnod(nodes), mftdsp(nodes),  &
& stat=status )

    if( status == 0 ) then
        if( lspin ) then
            allocate(  &
& wegud(nband*2*nkpnt), egvlud(nband*2*nkpnt),  &
& nsordr(nband*2*nkpnt),  &
& svtrxr(nband2*nknodx), stat=status )
            svtrxr = 0.d0
          else
            allocate(  &
& wegud(1), egvlud(1), nsordr(1), svtrxr(1),  &
& stat=status )
        end if
    end if

    if( status == 0 ) then
        if( .not.lgamma .or. ltddft .or. lnoncollinear ) then
            allocate( dmtrxi(nband2*nknodx), biji(nband2),  &
& bm1i(nband2), bx0i(nband2),  &
& stat=status )
          else
            allocate( dmtrxi(1), biji(1), bm1i(1), bx0i(1),  &
& stat=status )
        end if
    end if

    if( status == 0 ) then
        if( .not.lgamma .and. lspin ) then
            allocate( svtrxi(nband2*nknodx),  &
& stat=status )
            svtrxi = 0.d0
          else
            allocate( svtrxi(1),  &
& stat=status )
        end if
    end if

    the_mem =  &
& + 8.d0 * ( nband*(2+nnspin) * nkpnt + nband2*(nknodx+3) )  &
& + 4.d0 * ( nband*( nkpnt + 1 ) )  &
& + 8.d0 * ( nband*8 + iblock*nband + nbnod*nkpnt*nnspin )  &
& + 4.d0 * ( nodes*2 )
    if( lspin ) then
        the_mem = the_mem  &
& + 8.d0 * ( nband*2 * 2 * nkpnt )  &
& + 4.d0 * ( nband*2 * nkpnt )  &
& + 8.d0 * ( nband2 * nknodx )
    end if
    if( .not.lgamma .or. ltddft ) then
        the_mem = the_mem  &
& + 8.d0 * ( nband2 * (nknodx+3) )
    end if
    if( .not.lgamma .and. lspin ) then
        the_mem = the_mem  &
& + 8.d0 * ( nband2 * nknodx )
    end if


    if( status == 0 ) then
        if( jhybrid /= 0 ) then
            allocate( iod_m(nbnodx), iod_buf(nbnodx),  &
& stat=status )
            the_mem = the_mem  &
& + 4.d0 * ( size(iod_m) + size(iod_buf) )
        end if
    end if

    if( status == 0 ) then
        if( lnoncollinear ) then
            allocate( wvratio(nband,2,nkpnt), wvmagxy(nband,2,nkpnt),  &
& stat=status )
        else
            allocate( wvratio(1,1,1), wvmagxy(1,1,1), stat=status )
        end if
        the_mem = the_mem + 8.d0 * ( size(wvratio) + size(wvmagxy) )
    end if

    !------error trap
    call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'pwlda_alloc', .true. )

    wvratio(:,:,:) = 0.d0
    wvmagxy(:,:,:) = 0.d0

    !-----first-order SO effects by purturbation calculation
!    call SOp_variables_alloc( nfile, myid, nodes )

end if


if( nfftk   > lplnfftk   .or.  &
&   nfftkd  > lplnfftkd  .or.  &
&   nfftke  > lplnfftke  .or.  &
&   noddatx > lplnoddatx .or.  &
&   nplw5ex > ngenhex    .or.  &
&   nplwex  > lplhex0    .or.  &
&   nplwcex > lplchex    .or.  &
&   nplw7ex > lpl7hex0   .or.  &
&   npnodx  > lplnpnod0  .or.  &
&   npnodx  > lplnpnodx      ) then

    !-----if already allocated, deallocate arrays
    if( allocated(gdcr) ) then

        !-----save wavefunctions in temporal variables or a scratch file
        call savewv_in_scratch( nfile, myid, nodes, ierror )
        !------error trap
        status = abs(ierror)
        call gimax(status)
        if( status /= 0 ) call fstop( nfile, myid, nodes,  &
& 'open scratch file error in pwlda_alloc' )

        the_mem =  &
& + 8.d0 * ( size(gdcr) + size(rhgr) + size(prveig) + size(prvrho) )  &
& + 8.d0 * size(cgjr)  &
#if CCREAL4
& + 4.d0 * ( size(cgjr4) + size(gdcr4) + size(bufcr4) )  &
#endif
& + 8.d0 * ( size(rhcr) + size(bufcr) + size(prcd) + size(gnk)  &
&          + size(hnk) + size(fft3x) + size(fft3y) )  &
& + 16.d0* ( size(fftwork) )  &
& + 4.d0 * ( size(mfd2ft) + size(mshglb) + size(mftwrk) )  &
& + 8.d0 * ( size(glocal) + size(apk) + size(eigr) + size(rinr)  &
&          + size(thrhgr) + size(pdbuf) )

        !-----deallocate arrays
        deallocate( gdcr, rhgr, prveig, prvrho,  &
#if CCREAL4
& cgjr, cgjr4, gdcr4, bufcr4,  &
#else
& cgjr,  &
#endif
& rhcr, bufcr, prcd, gnk, hnk, fft3x, fft3y, fftwork,  &
& mfd2ft, mshglb, glocal, apk, eigr, mftwrk, rinr, thrhgr, pdbuf,  &
& stat=status )

        if( allocated(gdcrsv) .and. status == 0 ) then
            the_mem = the_mem + 8.d0 * ( size(gdcrsv) + size(rhcrsv) + size(hlocal) )
            deallocate( gdcrsv, rhcrsv, hlocal, stat=status )
        end if

        if( allocated(hcsr) .and. status == 0 ) then
            the_mem = the_mem + 8.d0 * ( size(hcsr) + size(sck) )
            deallocate( hcsr, sck, stat=status )
        end if

        if( allocated(hcsrsv) .and. status == 0 ) then
            deallocate( hcsrsv, stat=status )
            the_mem = the_mem + 8.d0 * ( lpnbndx )
        end if

        if( allocated(rhgrsoft) .and. status == 0 ) then
            deallocate( rhgrsoft, stat=status )
            the_mem = the_mem + 8.d0 * ( ngenhex*nspnmx )
        end if

        if( allocated(eigi) .and. status == 0 ) then
            deallocate( eigi, apki, scwr, scwi, stat=status )
            the_mem = the_mem + 8.d0 * (lplnfftk+1)*4
        end if

        if( allocated(pgdcr) .and. status == 0 ) then
            the_mem = the_mem + 8.d0 * ( size(pgdcr) )
            deallocate( pgdcr, stat=status )
        end if

        if( allocated(pgdcrsv) .and. status == 0 ) then
            the_mem = the_mem + 8.d0 * ( size(pgdcrsv) )
            deallocate( pgdcrsv, stat=status )
        end if

        if( allocated(rhgrph) .and. status == 0 ) then
            the_mem = the_mem + 8.d0 * ( size(rhgrph) )
            deallocate( rhgrph, stat=status )
        end if

        if( allocated(hfrhcr) .and. status == 0 ) then
            the_mem = the_mem  &
& + 8.d0 * ( size(hfrhcr) + size(hfv) )
            deallocate( hfrhcr, hfv, stat=status )
        end if

        if( allocated(cgjr_buf) .and. status == 0 ) then
            the_mem = the_mem  &
& + 8.d0 * ( size(cgjr_m) + size(cgjr_buf) + size(hfchg) + size(hfwav)  &
&          + size(hfrhgr) + size(hfrhcr_m) )
!&          + size(hfrhgr) + size(hfrhcr) + size(thehfrhcr) + size(hfrhcr_m) )
            deallocate( cgjr_m, cgjr_buf, hfchg, hfwav, hfrhgr,  &
& hfrhcr_m, stat=status )
!& hfrhcr, thehfrhcr, hfrhcr_m, stat=status )
        end if

        if( allocated(hfchgi) .and. status == 0 ) then
            the_mem = the_mem  &
& + 8.d0 * ( size(hfchgi) + size(hfvi) + size(hfwavi) + size(hfeiki) + size(hfeikj)  &
& + size(hfeikj_tmp) + size(hfcos) + size(hfsin) )
            deallocate( hfchgi, hfvi, hfwavi, hfeiki, hfeikj, hfeikj_tmp,  &
& hfcos, hfsin, stat=status )
        end if

        if( allocated(eig2r) .and. status == 0 ) then
            the_mem = the_mem + 8.d0 * ( size(eig2r) + size(eig2i) )
            deallocate( eig2r, eig2i, stat=status )
        end if

        !------error trap
        call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'pwlda_alloc(2)', .true. )

    end if


    !-----allocate arrays
    lplnfftk   = nfftk  * pwscale
    lplnfftkd  = nfftkd * pwscale
    lplnfftke  = nfftke * pwscale
    lplnoddatx = noddatx* pwscale
    ngenhex    = nplw5ex* pwscale
    lplhex0    = nplwex * pwscale
    lplchex    = nplwcex* pwscale
    lpl7hex0   = nplw7ex* pwscale
    lplnpnod0  = npnodx * pwscale
    lplnpnodx  = npnodx * pwscale

    lplhex     = lplhex0   * pkscale * dble(ncscale)
    lpl7hex    = lpl7hex0  * pkscale * dble(ncscale)
    lplnpnod   = lplnpnod0 * pkscale !* dble(ncscale)  --- npnodx is already doubled.

    lpnbndx = max( 2*lplnpnod*nband, lplhex*nbnod )
    npnnbnx = lplnpnodx*nbnodx*2*pkscale

    nspnmx2 = nspnmx*ncscale*ncscale
    allocate( gdcr(lpnbndx), rhgr(ngenhex*nspnmx2),  &
& prveig(lplnpnod*nband*2*nprvmx*nspnmx),  &
& prvrho(ngenhex*ncprvmx*nspnmx2),  &
#if CCREAL4
& cgjr(lpnbndx), cgjr4(lpl7hex*nband),  &
& gdcr4(lpnbndx), bufcr4(npnnbnx),  &
#else
& cgjr(lplhex*nband),  &
#endif
& rhcr(lpnbndx), bufcr(npnnbnx),  &
& prcd(lplhex), gnk(lplhex), hnk(lplhex),  &
& fft3x(lplnfftkd), fft3y(lplnfftkd), fftwork(lplnfftke),  &
& mfd2ft(lplnfftk),  &
& mshglb(lplnfftk),  &
& apk(lplnfftk+1), eigr(lplnfftk+1), mftwrk(3*lplnoddatx),  &
& rinr(ngenhex*nspnmx2), thrhgr(ngenhex),  &
& pdbuf(lplnpnod*2),  &
& stat=status )

    the_mem =  &
& + 8.d0 * ( size(gdcr) + size(rhgr) + size(prveig) + size(prvrho) )  &
& + 8.d0 * size(cgjr)  &
#if CCREAL4
& + 4.d0 * ( size(cgjr4) + size(gdcr4) + size(bufcr4) )  &
#endif
& + 8.d0 * ( size(rhcr) + size(bufcr) + size(prcd) + size(gnk)  &
&          + size(hnk) + size(fft3x) + size(fft3y) )  &
& + 16.d0* ( size(fftwork) )  &
& + 4.d0 * ( size(mfd2ft) + size(mshglb) + size(mftwrk) )  &
& + 8.d0 * ( size(apk) + size(eigr) + size(rinr)  &
&          + size(thrhgr) + size(pdbuf) )

    if( status == 0 ) then
        if( .not.lnoncollinear ) then
            allocate( glocal(lplnfftk), stat=status )
        else
            !-----noncollinear magnetism
            allocate( glocal(lplnfftk*4), stat=status )
        end if
        the_mem = the_mem + 8.d0 * size(glocal)
    end if

    if( status == 0 ) then
        if( lspin ) then
            allocate( gdcrsv(lpnbndx),  &
& rhcrsv(lpnbndx),  &
& hlocal(lplnfftk), stat=status )
          else
            allocate( gdcrsv(1), rhcrsv(1), hlocal(1),  &
& stat=status )
        end if
        the_mem = the_mem  &
& + 8.d0 * ( size(gdcrsv) + size(rhcrsv) + size(hlocal) )
    end if

    if( status == 0 ) then
        if( lvand ) then
            allocate( hcsr(lpnbndx), sck(lplhex),  &
& stat=status )
            if( lspin ) then
                allocate( hcsrsv(lpnbndx), stat=status )
              else
                allocate( hcsrsv(1), stat=status )
            end if
        else
            if( jhybrid /= 0 .or. lefield_islts .or. lnoncollinear ) then
                allocate( hcsr(lplhex), sck(lplhex), hcsrsv(1), stat=status )
            else
                allocate( hcsr(1), sck(1), hcsrsv(1), stat=status )
            end if
        end if
        the_mem = the_mem  &
& + 8.d0 * ( size(hcsr) + size(sck) + size(hcsrsv) )
    end if

!    if( status == 0 ) then
!        if( lpaw ) then
!            allocate( rhgrsoft(ngenhex*nspnmx), stat=status )
!            the_mem = the_mem + 8.d0 * ( ngenhex*nspnmx )
!        end if
!    end if

    if( status == 0 ) then
        if( .not.lgamma .or. lnoncollinear ) then
            allocate( eigi(lplnfftk+1), apki(lplnfftk+1),  &
& scwr(lplnfftk+1), scwi(lplnfftk+1),  &
& stat=status )
          else
            allocate( eigi(1), apki(1), scwr(1), scwi(1),  &
& stat=status )
        end if
        the_mem = the_mem  &
& + 8.d0 * ( size(eigi) + size(apki) + size(scwr) + size(scwi) )
    end if

    if( status == 0 ) then
!        if( laspc ) then
        if( laspc .and. .false. ) then !---not used
            allocate( pgdcr(lpnbndx), stat=status )
            if( lspin ) then
                allocate( pgdcrsv(lpnbndx), stat=status )
              else
                allocate( pgdcrsv(1), stat=status )
            end if
        else
            allocate( pgdcr(1), pgdcrsv(1), stat=status )
        end if
        the_mem = the_mem  &
& + 8.d0 * ( size(pgdcr) + size(pgdcrsv) )
    end if

    if( status == 0 ) then
        if( lrtddft ) then
            allocate( rhgrph(ngenhex,2), stat=status )
            the_mem = the_mem + 8.d0 * ( size(rhgrph) )
        end if
    end if

    if( status == 0 ) then
        if( jhybrid /= 0 .or. lefield_islts ) then
            allocate( hfrhcr(lpnbndx,1:nspnmx), hfv(lplnfftk), stat=status )
            the_mem = the_mem  &
& + 8.d0 * ( size(hfrhcr) + size(hfv) )
        else
            if( .not.allocated(hfrhcr) ) then
                allocate( hfrhcr(1,1:nspnmx), hfv(1), stat=status )
            end if
        end if
    end if

    if( status == 0 ) then
        if( jhybrid /= 0 ) then
            lpnbndx2 = max( lpnbndx, lplhex*nbnodx )
            allocate( cgjr_m(lpnbndx2*mxkpdupnm), cgjr_buf(lpnbndx2*mxkpdupnm),  &
& hfchg(lplnfftk), hfwav(lplnfftk+1),  &
& hfrhcr_m(lpnbndx2), stat=status )
!& hfrhcr(lpnbndx,1:nspnmx), thehfrhcr(lpnbndx), hfrhcr_m(lpnbndx2), stat=status )
            if( lgamma ) then
                allocate( hfrhgr(ngenhex), stat=status )
            else
                allocate( hfrhgr(ngenhex*2+1), stat=status )
            end if
            the_mem = the_mem  &
& + 8.d0 * ( size(cgjr_m) + size(cgjr_buf) + size(hfchg) + size(hfwav)  &
&          + size(hfrhgr) + size(hfrhcr_m) )
!&          + size(hfrhgr) + size(hfrhcr) + size(thehfrhcr) + size(hfrhcr_m) )
            if( .not.lgamma .and. status == 0 ) then
                khffftd = max(mshmx, mshmy, mshmz)
                allocate( hfchgi(lplnfftk), hfvi(lplnfftk), hfwavi(lplnfftk+1),  &
& hfeiki(lplnfftk*2), hfeikj(lplnfftk*2*mxkpdupnm), hfeikj_tmp(lplnfftk*2*mxkpdupnm),  &
& hfcos(khffftd*2*3), hfsin(khffftd*2*3),  &
& stat=status )
                the_mem = the_mem + 8.d0 * ( size(hfchgi) + size(hfvi) + size(hfwavi)  &
& + size(hfeiki) + size(hfeikj) + size(hfeikj_tmp) + size(hfcos) + size(hfsin) )
            end if
        else
            if( .not.allocated(cgjr_m) ) then
                allocate( cgjr_m(1), hfvi(1), stat=status )
            end if
        end if
    end if

    if( status == 0 ) then
        if( lnoncollinear ) then
            allocate( eig2r(lplnfftk+1), eig2i(lplnfftk+1), stat=status )
            the_mem = the_mem + 8.d0*( size(eig2r) + size(eig2i) )
        end if
    end if


    !------error trap
    call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'pwlda_alloc(2)', .true. )

    gdcr(1:lpnbndx) = 0.d0
    cgjr(1:lpnbndx) = 0.d0
    prveig(1:lplnpnod*nband*2*nprvmx*nspnmx) = 0.d0
    prvrho(1:ngenhex*ncprvmx*nspnmx2) = 0.d0
!    if( laspc ) then
    if( laspc .and. .false. ) then !---not used
        pgdcr(1:lpnbndx) = 0.d0
    end if
    if( lspin ) then
        gdcrsv(1:lpnbndx)  = 0.d0
        hlocal(1:lplnfftk) = 0.d0
    end if
    if( jhybrid /= 0 ) then
        cgjr_m(1:lpnbndx2*mxkpdupnm) = 0.d0
    end if
    if( jhybrid /= 0 .or. lefield_islts ) then
        hfrhcr(1:lpnbndx,1:nspnmx) = 0.d0
    end if

    !-----read wavefunctions in temporal variables or a scratch file
    call readwv_in_scratch( nfile, myid, nodes, ierror )

    !------error trap
    status = abs(ierror)
    call gimax(status)
    if( status /= 0 ) call fstop( nfile, myid, nodes,  &
& 'read error from scratch file in pwlda_alloc' )

end if


if( lspin ) then
    do i = 1, lpnbndx
       rhcrsv(i) = 0.d0
    end do
    if( lvand ) then
        do i = 1, lpnbndx
           hcsrsv(i) = 0.d0
        end do
    end if
!    if( laspc ) then
    if( laspc .and. .false. ) then !---not used
        do i = 1, lpnbndx
           pgdcrsv(i) = 0.d0
        end do
    end if
end if

bzan = 1.d0

if( licall ) then
    !-----allocate memory for work arrays in fermi.f
    if( lspin .or. lrela .and. lrela_sop .or. lsoc .and. .not.lsoc_full ) then
        call fermi_alloc( nfile, myid, nodes, alloc_mem,  &
& nband*nkpnt*2 )
      else
        call fermi_alloc( nfile, myid, nodes, alloc_mem,  &
& nband*nkpnt )
    end if
end if


#ifdef LIBFFTW3
!---memory allocation for FFTW3
if( jhybrid == 0 ) then
    call fftw3_alloc( nfile, myid, nodes,  &
& alloc_mem, nfftkd, pwscale, lgamma, lnoncollinear, kfft0 )
else
    call fftw3_alloc( nfile, myid, nodes,  &
& alloc_mem, nfftkd, pwscale, lgamma, lnoncollinear, nfftkd )
end if
#endif


!---memory allocation for charge density mixing
call mxchgg_alloc( nfile, myid, nodes, &
& alloc_mem, lspin, nplwcex, imxchg, itratn, imxspin, nmxspin,  &
& cmetric, spinmetric, wkerker, pwscale, lnoncollinear )


!---memory allocation for k-point sampling
!call pwlda_k_alloc( nfile, myid, nodes,  &
!& alloc_mem, lspin, lvand,  &
!& nband, nbnod, nplwex, npnodx, nprvmx, nspnmx,  &
!& pwscale, jhybrid, ncscale )


!---memory allocation for TDDFT
!call pwlda_tddft_alloc( nfile, myid, nodes,  &
!& alloc_mem, nband, nbnod, nplwex, npnodx, nspnmx, pwscale )

!---memory allocation for TDDFT-FSSH
call tddft_fssh_pw_alloc( nfile, myid, nodes, &
& alloc_mem, nspnmx, nplw5ex, pwscale, lpnbndx )


licall = .false.

return
end subroutine




subroutine cfftmsh_alloc( nfile, myid, nodes,  &
& alloc_mem, lnlpp_r, kfft1, kfft2, kfft3, kfft0, pwscale,  &
& jhybrid, mxkpdupnm, nbnod, lefield_islts )
!-----------------------------------------------------------------------
!     allocate memory for the plane wave method
!-----------------------------------------------------------------------
use pwlda_pw
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
logical :: lnlpp_r
integer :: kfft1, kfft2, kfft3, kfft0
real*8  :: pwscale
integer :: jhybrid, mxkpdupnm, nbnod
logical :: lefield_islts

!-----declare local variables
integer :: status
real*8  :: the_mem
integer :: i, m1, ix, iy, iz
logical :: lgamma
integer :: lplntotfd_c = 0    ! ntotfd_c
integer :: lplkfft0    = 0    ! kfft0
save lplntotfd_c, lplkfft0


call get_lgamma( lgamma )

if( ntotfd_c > lplntotfd_c  .or.  &
&   kfft0    > lplkfft0         ) then

    !-----if already allocated, deallocate arrays
    if( allocated(mfd2ft_c) ) then
        !------deallocate memory
        the_mem = 4.d0 * size(mfd2ft_c)

        deallocate( mfd2ft_c, stat=status )

        if( allocated(mshglb_c) .and. status == 0 ) then
            the_mem = the_mem + 4.d0 * ( size(mshglb_c) + size(mshgnx)*3 )
            deallocate( mshglb_c, mshgnx, mshgny, mshgnz,  &
& stat=status )
        end if

        if( allocated(hfv_c) .and. status == 0 ) then
            the_mem = the_mem + 8.d0 * size(hfv_c)
            deallocate( hfv_c, stat=status )
        end if

        if( allocated(hfvi_c) .and. status == 0 ) then
            the_mem = the_mem + 8.d0 * size(hfvi_c)
            deallocate( hfvi_c, stat=status )
        end if

        !------error trap
        call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'cfftmsh_alloc', .true. )
    end if

    lplntotfd_c = ntotfd_c * pwscale
    lplkfft0    = kfft0    * pwscale

    !------allocate memory
    allocate( mfd2ft_c(lplntotfd_c),  &
& stat=status )

    the_mem = 4.d0 * size(mfd2ft_c)


    if( status == 0 ) then
        if( lnlpp_r ) then
            allocate( mshglb_c(lplkfft0),  &
& mshgnx(lplntotfd_c), mshgny(lplntotfd_c), mshgnz(lplntotfd_c),  &
& stat=status )
          else
            allocate( mshglb_c(1), mshgnx(1), mshgny(1), mshgnz(1),  &
& stat=status )
        end if

        the_mem = the_mem + 4.d0 * ( size(mshglb_c) + size(mshgnx)*3 )
    end if


    if( status == 0 ) then
        if( jhybrid /= 0 .or. lefield_islts ) then
            if( .not.lgamma ) then
                allocate( hfv_c(lplkfft0*nbnod), stat=status )
!                allocate( hfv_c(lplkfft0*mxkpdupnm*nbnod), stat=status )
            else
                allocate( hfv_c(lplkfft0*nbnod), stat=status )
            end if
        else
            allocate( hfv_c(1), stat=status )
        end if

        the_mem = the_mem + 8.d0 * size(hfv_c)

        if( .not.lgamma ) then
            if( jhybrid /= 0 ) then
                allocate( hfvi_c(lplkfft0*nbnod), stat=status )
!                allocate( hfvi_c(lplkfft0*mxkpdupnm*nbnod), stat=status )
            else
                allocate( hfvi_c(1), stat=status )
            end if
            the_mem = the_mem + 8.d0 * size(hfvi_c)
        end if
    end if

    !------error trap
    call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'cfftmsh_alloc', .true. )

end if


!-----set mfd2ft_c
do i = 1, kfft0
   m1 = mshglb(i)
   if( m1 >= 1 ) mfd2ft_c(m1) = i
end do


!-----set mesh to calculate nonlocal pp on real-space grid
if( lnlpp_r ) then
    do i = 1, kfft0
       mshglb_c(i) = mshglb(i)
    end do

    i = 0
    do iz = 1, kfft3
    do iy = 1, kfft2
    do ix = 1, kfft1
       i = i + 1
       m1 = mshglb_c(i)
       if( m1 >= 1 ) then
!                 mfd2ft_c(m1) = ix + kfft1*( (iy-1) + kfft2*(iz-1) )
           mshgnx(m1) = ix
           mshgny(m1) = iy
           mshgnz(m1) = iz
       end if
    end do
    end do
    end do

    !--- for k-point sampling
!    call eikr_alloc( nfile, myid, nodes,  &
!& alloc_mem, kfft1, kfft2, kfft3, kfft0, pwscale, jhybrid )

end if

hfv_c(:) = 0.d0
if( .not.lgamma ) hfvi_c(:) = 0.d0


return
end subroutine




subroutine pwlda_grid_alloc( nfile, myid, nodes,  &
& alloc_mem, multg, nultg )
!-----------------------------------------------------------------------
!     allocate memory for multigrid
!-----------------------------------------------------------------------
use pwlda_grid
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: multg, nultg

!-----declare local variables
integer :: status
real*8  :: the_mem


!------allocate memory
allocate( mulnd2(multg), kulnd2(multg),  &
&         nulnd2(nultg), lulnd2(nultg),  &
& stat=status )

the_mem = 4.d0 * ( multg + nultg )*2

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'pwlda_grid_alloc', .true. )


return
end subroutine




subroutine pwlda_mgrid_alloc( nfile, myid, nodes, npx, npy, npz,  &
& alloc_mem, nd1vks, nd2v, multg, nultg, lspin, pwscale, lrtddft,  &
& ltddft_fssh, lnoncollinear, ncprvmx, lwell, lefield, lsawtooth )
!-----------------------------------------------------------------------
!     allocate memory for multigrid
!-----------------------------------------------------------------------
use pwlda_grid
use ncmagne_variables
implicit none
integer :: nfile(*), myid, nodes
integer :: npx, npy, npz
real*8  :: alloc_mem
integer, dimension(3) :: nd1vks
integer :: nd2v
integer :: multg
integer :: nultg
logical :: lspin
real*8  :: pwscale
logical :: lrtddft, ltddft_fssh
logical :: lnoncollinear
integer :: ncprvmx
logical :: lwell, lefield, lsawtooth

!-----declare local variables
real*8  :: the_mem
integer :: status, ierror
integer :: incxyz, incrx, incry, incrz
integer :: tmshmx, tmshmy, tmshmz, tsnodx, tsnody, tsnodz
integer :: tsndx1, tsndy1, tsndz1, tsndx2, tsndy2, tsndz2
integer :: mul
integer :: i, msnd43
integer :: nspnmx
integer :: lplmulexf  = 0    ! mulexf
integer :: lplmuldatx = 0    ! muldatx
integer :: lplmulyzx  = 0    ! mulyzx
integer :: lplmulzxx  = 0    ! mulzxx
integer :: lplmulxyx  = 0    ! mulxyx
integer :: lplnodyzh  = 0    ! nodyzh
integer :: lplnodzxh  = 0    ! nodzxh
integer :: lplnodxyh  = 0    ! nodxyh
integer :: lplmshndv  = 0    ! mshndv
integer :: lplmsnd43  = 0    ! msnd43
integer :: lplnoddatx = 0    ! noddatx
integer :: lplmulext  = 0    ! mulext
integer :: lplnodexe  = 0    ! nodexe
logical :: licall = .true.
save lplmulexf, lplmuldatx,lplmulyzx, lplmulzxx, lplmulxyx,  &
&    lplnodyzh, lplnodzxh, lplnodxyh, lplmshndv, lplmsnd43,  &
&    lplnoddatx,lplmulext, lplnodexe
save licall


if( lspin ) then
    nspnmx = 2
  else
    nspnmx = 1
end if

if( licall ) then
    mulnpx = 2*multg*(npx+1)
    mulnpy = 2*multg*(npy+1)
    mulnpz = 2*multg*(npz+1)
    ndv2x3 = nd2v*(nd2v+2)*3

    !------allocate memory
    allocate( nmmcnt(nodes), nmmdsp(nodes), nmmrc1(0:nodes-1),  &
& nmmrc2(0:nodes-1), meshx1(0:npx-1), meshx(0:npx-1),  &
& meshy1(0:npy-1), meshy(0:npy-1), meshz1(0:npz-1), meshz(0:npz-1),  &
& mshnod(multg), mshnew(multg), mulpit(multg),  &
& mulx1(multg), mulx2(multg), muly1(multg), muly2(multg),  &
& mulz1(multg), mulz2(multg),  &
& mulpms(multg), mshx1(multg), mshy1(multg),  &
& mshz1(multg), mshx(multg),  mshy(multg),  mshz(multg),  &
& mulpem(multg),  &
& mulpsx(multg), mulpsy(multg), mulpsz(multg),  &
& mulyzs(multg), mulzxs(multg), mulxys(multg),  &
& mulprx(multg), mulpry(multg), mulprz(multg),  &
& mulyzr(multg), mulzxr(multg), mulxyr(multg),  &
& ismx(mulnpx), ismy(mulnpy), ismz(mulnpz),  &
& irmx(mulnpx), irmy(mulnpy), irmz(mulnpz),  &
& ismfpx(mulnpx), ismfpy(mulnpy), ismfpz(mulnpz),  &
& irmfpx(mulnpx), irmfpy(mulnpy), irmfpz(mulnpz),  &
& ncomct(12,multg), amcdv2(ndv2x3,multg),  &
& cdv1(ndv2x3,nultg), cdv2(ndv2x3,nultg),  &
& mulgmx(0:npx-1), mulgmy(0:npy-1), mulgmz(0:npz-1),  &
& nrstrct(multg), prevres(multg),  &
& stat=status )

    the_mem =  &
& 4.d0 * ( nodes*4 + (npx + npy + npz)*2 + multg*29  &
&       + ( mulnpx + mulnpy + mulnpz )*4 + 12*multg )  &
& + 8.d0 * ( ndv2x3*multg*1 + ndv2x3*nultg*2 )  &
& + 4.d0 * ( npx + npy + npz )  &
& + 8.d0 * multg + 4.d0 * multg

    !------error trap
    call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'pwlda_mgrid_alloc', .true. )

end if


mshmx = nd1vks(1)
mshmy = nd1vks(2)
mshmz = nd1vks(3)
lnmx  = mshmx*mshmy*mshmz
ndata = lnmx
msnodx = ( mshmx + npx - 1 )/npx
msnody = ( mshmy + npy - 1 )/npy
msnodz = ( mshmz + npz - 1 )/npz
noddatx = msnodx*msnody*msnodz

msndx1 = - nd2v + 1
msndy1 = - nd2v + 1
msndz1 = - nd2v + 1
msndx2 = msnodx + nd2v
msndy2 = msnody + nd2v
msndz2 = msnodz + nd2v
incrx = msndx2 - msndx1 + 1
incry = msndy2 - msndy1 + 1
incrz = msndz2 - msndz1 + 1
incxyz = incrx*incry*incrz

nodyzx = msnody*msnodz*nd2v*2
nodzxx = (msnodz*msnodx*nd2v + msnodz*2)*2
nodxyx = ( msnodx*msnody*nd2v+(msnodx+msnody)*2+4 )*2

!-----define array dimensions for multigrid : muldatx
tmshmx = mshmx
tmshmy = mshmy
tmshmz = mshmz
muldatx = noddatx
mulexf  = incxyz
mulyzx  = nodyzx
mulzxx  = nodzxx
mulxyx  = nodxyx
do mul = 2, multg
   tmshmx = ( tmshmx + 1 )/2
   tmshmy = ( tmshmy + 1 )/2
   tmshmz = ( tmshmz + 1 )/2
   tsnodx = ( tmshmx + npx - 1 )/npx
   tsnody = ( tmshmy + npy - 1 )/npy
   tsnodz = ( tmshmz + npz - 1 )/npz
   muldatx = muldatx + tsnodx*tsnody*tsnodz

   tsndx1 = - mulnd2(mul) + 1
   tsndy1 = - mulnd2(mul) + 1
   tsndz1 = - mulnd2(mul) + 1
   tsndx2 = tsnodx + mulnd2(mul)
   tsndy2 = tsnody + mulnd2(mul)
   tsndz2 = tsnodz + mulnd2(mul)
   incrx = tsndx2 - tsndx1 + 1
   incry = tsndy2 - tsndy1 + 1
   incrz = tsndz2 - tsndz1 + 1
   incxyz = incrx*incry*incrz
   mulexf = mulexf + incxyz

   mulyzx = mulyzx + tsnody*tsnodz*mulnd2(mul)*2
   mulzxx = mulzxx + ( tsnodz*tsnodx*mulnd2(mul) + tsnodz*2 )*2
   mulxyx = mulxyx + ( tsnodx*tsnody*mulnd2(mul)  &
&                    + ( tsnodx + tsnody )*2 + 4 )*2
end do

nodext = noddatx  &
&      + ( msnody*msnodz + msnodz*msnodx + msnodx*msnody )*nd2v*2
mulext = muldatx + mulyzx + mulzxx + mulxyx
nodexe = nodext + ( msnodx + msnody + msnodz )*4 + 8
nodexf = nodext + ( ( msnodx + msnody + msnodz )*4  &
&                                          + 8*nd2v )*nd2v*nd2v

mshndv = max( msnody*msnodz*nd2v,  &
&             msnodz*msnodx*nd2v + 2*msnodz,  &
&             msnodx*msnody*nd2v + 2*(msnodx+msnody) + 8 ) * 3

nodyzh = nodyzx/2
nodzxh = nodzxx/2
nodxyh = nodxyx/2
msndx4 = msnodx+4
msndy4 = msnody+4
msndz4 = msnodz+4
msnd43 = msndx4*msndy4*msndz4


if( mulexf  > lplmulexf  .or.  &
&   muldatx > lplmuldatx .or.  &
&   mulyzx  > lplmulyzx  .or.  &
&   mulzxx  > lplmulzxx  .or.  &
&   mulxyx  > lplmulxyx  .or.  &
&   nodyzh  > lplnodyzh  .or.  &
&   nodzxh  > lplnodzxh  .or.  &
&   nodxyh  > lplnodxyh  .or.  &
&   mshndv  > lplmshndv  .or.  &
&   msnd43  > lplmsnd43  .or.  &
&   noddatx > lplnoddatx .or.  &
&   mulext  > lplmulext  .or.  &
&   nodexe  > lplnodexe      ) then

    !-----if already allocated, deallocate arrays
    if( allocated(mshxyz) ) then

        the_mem =  &
& 4.d0 * ( lplmulexf + lplmuldatx*3  &
&       + lplmulyzx + lplmulzxx + lplmulxyx  &
&       + lplmulyzx + lplmulzxx + lplmulxyx  &
&       + 3*(lplnodyzh+lplnodzxh+lplnodxyh)*2  &
&       + lplmshndv*2 + lplmsnd43 )  &
#ifdef VECTOR
& + 8.d0 * lplmulexf  &
#endif
& + 8.d0 * ( lplnoddatx*9 + lplmuldatx*7 + lplmulext  &
&         + lplnodexe*2 + 1 + size(vhar_out) + size(delrho) + size(xexc)  &
&         + size(ddelrho) )

        !------deallocate memory
        deallocate( mshxyz, mshnx, mshny, mshnz,  &
& ismshx, ismshy, ismshz, irmshx, irmshy, irmshz,  &
& jsmshx, jsmshy, jsmshz, idbuf, idbufr, fmshx,  &
#ifdef VECTOR
& xxx,  &
#endif
& rho, vhar, vext, vexc, eeexc, hdiag, rhocore, x, rk, vk,  &
& dbuf, dbufr,  &
& xx, tmpk, tmpl, tmpm, tmpn, vlshdp, vhshdp, vhar_out, delrho, xexc,  &
& ddelrho,  &
& stat=status )

        if( allocated(rhoud) .and. status == 0 ) then
            the_mem = the_mem + 8.d0 * ( size(rhoud) + size(vlocud) )
            deallocate( rhoud, vlocud, stat=status )
        end if

        if( allocated(rhoph) .and. status == 0 ) then
            the_mem = the_mem  &
& + 8.d0 * ( size(rhoph) + size(rhophe) + size(vexcph) + size(eeexcph) )
            deallocate( rhoph, rhophe, vexcph, eeexcph, stat=status )
        end if

        if( allocated(rhoudph) .and. status == 0 ) then
            the_mem = the_mem  &
& + 8.d0 * ( size(rhoudph) + size(vlocudph) + size(hdiagph) )
            deallocate( rhoudph, vlocudph, hdiagph, stat=status )
        end if

        if( allocated(rhomx) .and. status == 0 ) then
!            !-----save wavefunctions in temporal variables or a scratch file
!            call savemagne_in_scratch( nfile, myid, nodes, ierror )
!            !------error trap
!            status = abs(ierror)
!            call gimax(status)
!            if( status /= 0 ) call fstop( nfile, myid, nodes,  &
!& 'open scratch file error in pwlda_mgrid_alloc' )

            the_mem = the_mem  &
& + 8.d0 * ( size(rhomx) + size(rhomy) + size(rhomz)  & !+ size(prvrhom)  &
&          + size(rhomxini) + size(rhomyini) + size(rhomzini) )
            deallocate( rhomx, rhomy, rhomz,  & !prvrhom,  &
& rhomxini, rhomyini, rhomzini, stat=status )
        end if

        if( allocated(vwell) .and. status == 0 ) then
            the_mem = the_mem + 8.d0 * ( size(vwell) )
            deallocate( vwell, stat=status )
        end if

        if( allocated(vefield) .and. status == 0 ) then
            the_mem = the_mem + 8.d0 * ( size(vefield) )
            deallocate( vefield, stat=status )
        end if

        !------error trap
        call check_dealloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'pwlda_mgrid_alloc(2)', .true. )
    end if

    lplmulexf  = mulexf  * pwscale
    lplmuldatx = muldatx * pwscale
    lplmulyzx  = mulyzx  * pwscale*pwscale
    lplmulzxx  = mulzxx  * pwscale*pwscale
    lplmulxyx  = mulxyx  * pwscale*pwscale
    lplnodyzh  = nodyzh  * pwscale*pwscale
    lplnodzxh  = nodzxh  * pwscale*pwscale
    lplnodxyh  = nodxyh  * pwscale*pwscale
    lplmshndv  = mshndv  * pwscale*pwscale
    lplmsnd43  = msnd43  * pwscale
    lplnoddatx = noddatx * pwscale
    lplmulext  = mulext  * pwscale
    lplnodexe  = nodexe  * pwscale

    !------allocate memory
    allocate( mshxyz(lplmulexf),  &
& mshnx(lplmuldatx), mshny(lplmuldatx), mshnz(lplmuldatx),  &
& ismshx(lplmulyzx), ismshy(lplmulzxx), ismshz(lplmulxyx),  &
& irmshx(lplmulyzx), irmshy(lplmulzxx), irmshz(lplmulxyx),  &
& jsmshx(3*lplnodyzh*2), jsmshy(3*lplnodzxh*2),  &
& jsmshz(3*lplnodxyh*2),  &
& idbuf(lplmshndv), idbufr(lplmshndv), fmshx(lplmsnd43),  &
#ifdef VECTOR
& xxx(lplmulexf),  &
#endif
& rho(lplnoddatx), vhar(lplmuldatx), vext(lplnoddatx),  &
& vexc(lplnoddatx), eeexc(lplnoddatx),  &
& hdiag(lplmuldatx), rhocore(lplnoddatx),  &
& x(lplmulext), rk(lplmuldatx), vk(lplmuldatx), dbuf(lplnoddatx),  &
& dbufr(lplnoddatx), xx(0:lplnodexe), tmpk(lplnodexe),  &
& tmpl(lplmuldatx), tmpm(lplmuldatx),  &
& tmpn(lplmuldatx), vlshdp(lplnoddatx), vhshdp(lplnoddatx),  &
& stat=status )

    the_mem =  &
& 4.d0 * ( lplmulexf + lplmuldatx*3  &
&       + lplmulyzx + lplmulzxx + lplmulxyx  &
&       + lplmulyzx + lplmulzxx + lplmulxyx  &
&       + 3*(lplnodyzh+lplnodzxh+lplnodxyh)*2  &
&       + lplmshndv*2 + lplmsnd43 )  &
#ifdef VECTOR
& + 8.d0 * lplmulexf  &
#endif
& + 8.d0 * ( lplnoddatx*9 + lplmuldatx*7 + lplmulext  &
& + lplnodexe*2 + 1 )

    if( status == 0 ) then
        if( lspin .or. lnoncollinear ) then
            allocate( rhoud(lplnoddatx), stat=status )
            if( lspin ) then
                allocate( vlocud(lplnoddatx*2), stat=status )
            else ! if( lnoncollinear ) then
                allocate( vlocud(lplnoddatx*4), stat=status )
            end if
          else
            allocate( rhoud(1), vlocud(1), stat=status )
        end if
        the_mem = the_mem + 8.d0 * ( size(rhoud) + size(vlocud) )
    end if

    if( status == 0 ) then
        !-----noncollinear magnetism
        if( lnoncollinear ) then
            allocate( rhomx(lplnoddatx), rhomy(lplnoddatx), rhomz(lplnoddatx),  &
!& prvrhom(lplnoddatx*ncprvmx*3),  &
& rhomxini(lplnoddatx), rhomyini(lplnoddatx), rhomzini(lplnoddatx),  &
& stat=status )
            the_mem = the_mem  &
& + 8.d0 * ( size(rhomx) + size(rhomy) + size(rhomz)  & !+ size(prvrhom)  &
&          + size(rhomxini) + size(rhomyini) + size(rhomzini) )
          else
            allocate( rhomx(1), rhomy(1), rhomz(1),  & !prvrhom(1),  &
& rhomxini(1), rhomyini(1), rhomzini(1),  &
&stat=status )
        end if
    end if

    if( status == 0 ) then
        if( lrtddft ) then
            allocate( rhoph(lplnoddatx,2), rhophe(lplnoddatx),  &
& vexcph(lplnoddatx), eeexcph(lplnoddatx), stat=status )
            the_mem = the_mem  &
& + 8.d0 * ( size(rhoph) + size(rhophe) + size(vexcph) + size(eeexcph) )
            if( lspin ) then
                allocate( rhoudph(lplnoddatx), vlocudph(lplnoddatx*2),  &
& hdiagph(lplmuldatx), stat=status )
                the_mem = the_mem  &
& + 8.d0 * ( size(rhoudph) + size(vlocudph) + size(hdiagph) )
            end if
        end if
    end if

    if( status == 0 ) then
        if( ltddft_fssh ) then
            allocate( vhar_out(lplmuldatx), delrho(lplnoddatx),  &
& xexc(lplnoddatx), ddelrho(lplnoddatx*3), stat=status )
            the_mem = the_mem  &
& + 8.d0 * ( size(vhar_out) + size(delrho) + size(xexc) + size(ddelrho) )
        else
            allocate( vhar_out(1), delrho(1), xexc(1), ddelrho(1), stat=status )
        end if
    end if

    if( status == 0 ) then
        if( lwell ) then
            allocate( vwell(lplnoddatx), stat=status )
            the_mem = the_mem + 8.d0 * ( size(vwell) )
        end if
    end if

    if( status == 0 ) then
        if( lefield .and. lsawtooth ) then
            allocate( vefield(lplnoddatx), stat=status )
            the_mem = the_mem + 8.d0 * ( size(vefield) )
        end if
    end if

    !------error trap
    call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'pwlda_mgrid_alloc(2)', .true. )


    vhar(1:lplmuldatx) = 0.d0
    tmpk(1:lplnodexe) = 0.d0
    if( lnoncollinear ) then
        !-----noncollinear magnetism
        rhomx(1:lplnoddatx) = 0.d0
        rhomy(1:lplnoddatx) = 0.d0
        rhomz(1:lplnoddatx) = 0.d0
!        prvrhom(1:lplnoddatx*ncprvmx*3) = 0.d0

!        !-----read wavefunctions in temporal variables or a scratch file
!        call readmagne_in_scratch( nfile, myid, nodes, ierror )
!
!        !------error trap
!        status = abs(ierror)
!        call gimax(status)
!        if( status /= 0 ) call fstop( nfile, myid, nodes,  &
!& 'read error from scratch file in pwlda_mgrid_alloc' )
    end if


    !---allocate memory for EDA
!    call eda_mesh_alloc( nfile, myid, nodes,  &
!& alloc_mem, lplnoddatx )

    !---allocate memory for vdW
!    call vdw_mesh_alloc( nfile, myid, nodes, lplnoddatx )

    !---allocate memory for TDDFT-FSSH with LR-TDDFT
!    if( lrtddft .and. ltddft_fssh ) &
!& call lrtddft_mgrid_alloc( nfile, myid, nodes, lplnoddatx, nspnmx )

    !-----noncollinear magnetism
!    if( lnoncollinear ) call ncmxchg_alloc2( nfile, myid, nodes, &
!& alloc_mem, lplnoddatx )

end if


licall = .false.

return
end subroutine
