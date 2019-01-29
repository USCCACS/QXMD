Input Parameters: IN.PARAM 
=========================
.. Note:: The name of the main input file (e.g. 'IN.PARAM') must be specified in 'control/filename'.

The main input file is divided into sections corresponding to different
controls. Sections start with ``*SECTION_NAME`` and end with ``*end``. To
disable any section, simply prepend section title with a ``#``.

.. Example:: To disable section SECTION_NAME_2

::

   *SECTION_NAME_1
   ...
   *end

   #*SECTION_NAME_2
   ...
   *end

Each section is futher divided into subsections, given by a title in
parentheses (e.g. (QM-Nodes)), where input control parameters are
defined. To disable any subsection, simply prepend subsection title with
a '#'.

Example: To disable subsection_2 in section \*SECTION_NAME

::

   *SECTION_NAME
   (subsection_1)

   #(subsection_2)

   (subsection_3)
   *end

Note: Any unnecessary sections or subsections may be removed from
input.file, depending on your simulation needs (recommended only for
advanced users).


5.1 Parallel Section: \*parallel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


5.1.1 (QM-Nodes)
^^^^^^^^^^^^^^^^

**npx, npy, npz** = [1|even integer] [1|even integer] [1|even integer]

Default: **npx** = 1 **npy** = 1 **npz** = 1

QM-Nodes is a space delimited set of three integers **npx**, **npy**,
**npz**. **QXMD** uses a hybrid spatial and band decompisition for
parallel computing, and this set of three integers multiplied together
defines the number of MPI ranks. For more detail on parallelization.


5.1.2 (k-points)
^^^^^^^^^^^^^^^^

**npk** = [integer]

Default: **npk** = 1

Defines the number of MPI ranks to parallelize k-point sampling. See
section 5.9 for k-point sampling.


5.1.3 (linear-response TDDFT)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**nplr** = [integer]

Default: **nplr** = 1

Defines the number of MPI ranks for parallelization of matrix
computations for linear response TDDFT calculations.


5.1.4 (MD-nodes)
^^^^^^^^^^^^^^^^

**md_npx, md_npy, md_npz** = [1|even integer] [1|even integer] [1|even
integer]

MD-Nodes is a space delimited set of three integers: **md_npx, md_npy,
md_npz**. These variables are used only for spatial decomposition in
classical MD and divide-and-conqure QMD simulations. For linear scaling
DFT please see reference `[1-7] <#shimojo_2014>`__

Reference:
^^^^^^^^^^

1. "A divide-conquer-recombine algorithmic paradigm for large
spatiotemporal quantum molecular dynamics simulations," F. Shimojo, S.
Hattori, R. K. Kalia, M. Kunaseth, W. Mou, A. Nakano, K. Nomura, S.
Ohmura, P. Rajak, K. Shimamura & P. Vashishta J. Chem. Phys. 140, 18A529
(2014) 2. "Linear-scaling density-functional-theory calculations of
electronic structure based on real-space grids: design, analysis, and
scalability test of parallel algorithms," F. Shimojo, R. K. Kalia, A.
Nakano & P. Vashishta Comput. Phys. Commun. 140, 303 (2001) 3. "Embedded
divide-and-conquer algorithm on hierarchical real-space grids: parallel
molecular dynamics simulation based on linear-scaling density functional
theory," F. Shimojo, R. K. Kalia, A. Nakano & P. Vashishta, Comput.
Phys. Commun. 167, 151 (2005) 4. "Divide-and-conquer density functional
theory on hierarchical real-space grids: parallel implementation and
applications," F. Shimojo, R. K. Kalia, A. Nakano & P. Vashishta, Phys.
Rev. 77 085103 (2008) 5. "Scalable atomistic simulation algorithms for
materials research," A. Nakano, R. K. Kalia, P. Vashishta, T. J.
Campbell, S. Ogata, F. Shimojo & S. Saini Proc. Supercomputing, SC01
(ACM/IEEE, 2001) 6. "Metascalable quantum molecular dynamics simulations
of hydrogen-on-demand," K. Nomura, R. K. Kalia, A. Nakano, P. Vashishta,
K. Shimamura, F. Shimojo, M. Kunaseth, P. C. Messina& N. A. Romero Proc.
Supercomputing, SC14 (IEEE/ACM, 2014) 7. "Large nonadiabatic quantum
molecular dynamics simulations on parallel computer," F. Shimojo, S.
Ohmura, W. Mou, R. K. Kalia, A. Nakano & P. Vashishta Comput. Phys.
Commun. 184, 1 (2013)


5.2 Start(on/off) Section
~~~~~~~~~~~~~~~~~~~~~~~~~


5.2.1 (On/Off)
^^^^^^^^^^^^^^

**lstart** = [boolean]

Default: **lstart** = .FALSE.

.FALSE. : Start self-consistent field (SCF) iteration using a random
initial wavefunction .TRUE. : Continue SCF iteration from a previous
run's output wavefunction. Note: QM_*??? output files must be present in
the **data/** directory from the previous calculation.

Determines whether you would like to restart a successfully completed
simulation. Note that if you are restarting a NAQMD simulation, you will
also want to set **ltddft_start** to .TRUE. in the \*TDDFT-MD section to
properly restart a successfully completed simulation.


5.3 TDDFT-MD
~~~~~~~~~~~~


5.3.1 (On/Off)
^^^^^^^^^^^^^^

**ltddft** = [boolean]

Default: **ltddft** = .FALSE.

.FALSE. : Execute Adiabtic QMD based on density functional theory (DFT)
.TRUE. : Execute non-adiabtic QMD (NAQMD) based on time-dependent
density functional theory (TDDFT) `[1] <#Gross_1990>`__.

Determines whether to run QMD simulation under adiabatic or
non-adiabatic methods. Adiabatic methods simulate thermodynamic system
in ground state equilibrium, while non-adiabatic methods simulate
electronic excitations.


5.3.2 (FSSH)
^^^^^^^^^^^^

**ltddft_fssh** = [boolean]

Default: **ltddft_fssh** = .TRUE.

.TRUE. : Perform NAQMD based on Fewest Switches Surface Hopping (FSSH)
method `[2] <#Tully_1990>`__\  .FALSE. : Perform NAQMD based on based on
Ehrenfest dynamics (not yet implemented).

Determines the implementation method for electron state dynamics in
NAQMD. Fewest Switches surface hopping method is proposed by J. Tully
`[2] <#Tully_1990>`__ molecular dynamics simulation of the processes
including electronic transition. Next few flags ask for the
specification of FSSH method


5.3.3 (FSSH-switch)
^^^^^^^^^^^^^^^^^^^

**lfssh_switch** = [boolean]

Default: **lfssh_switch** = .TRUE.

.TRUE. : Allow electrons to move between excited states. .FALSE. : Keep
electronic occupations fixed.

Determines whether electronic occupations can change throughout the
NAQMD simulation.


5.3.4 (FSSH-ground-state-SCF)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lfssh_gsscf** = [boolean]

Default: **lfssh_gsscf** = .FALSE.

.TRUE. : SCF iterations performed based on the ground state .FALSE. :
SCF with the excited state

This parameter should always be set to **true** to obtain convergence.


5.3.5 (FSSH-charge mixing)
^^^^^^^^^^^^^^^^^^^^^^^^^^

**imxchg_fssh** = 0 \| 1 \| 2 \| 3 **aslh_fssh, bslh_fssh** = [real]
[real]

Default: **imxchg_fssh** = 1


=========== ===================
imxchg_fssh Charge Mixing Method                              Notes for aslh_fssh, bslh_fssh
=========== ===================
**0**       No Mixing                                         N/A
**1**       Pulay [`3- <#Pulay_1980>`__\ `4 <#Payne_1992>`__] 0.9, 0.6 (recommended values)
**2**       Anderson                                         
**3**       Simple                                            bslh_fssh not used
=========== ===================

**imxchg_fssh** is used to specify the method for charge mixing during
SCF iterations, while **aslh_fssh** and **bslh_fssh** are used as tuning
parameters. Note that this section is only considered when
**lfssh_gsscf** is set to .TRUE.


5.3.6 (FSSH-random-initialize)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lfssh_random** = [boolean] **rseed_fssh** = [real]

Default: **lfssh_random** = .FALSE.

**lfssh_random** = .TRUE. : Automatically seed random number generator.
**lfssh_random** = .FALSE. : Specify the seed for the random number
generator with value given by **rseed_fssh**.

Determines how to specify the seed for the random number generator used
by the FSSH method.


5.3.7 (Boltzmann factor for upward transition)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lfssh_boltzmn** = [boolean]

Default: **lfssh_boltzmn** = .TRUE.

.TRUE. : Multiply electronic transition probability by the Boltzmann
factor. .FALSE. : Leave electronic transtion probability unaltered.

Determines whether or not to multiply the electronic transition
probability by the Boltzmann factor when the electronic excitation
energy increases due to the transition. This is used in order to
approximately satisfy the detailed balance condition.


5.3.8 (velocity scaling)
^^^^^^^^^^^^^^^^^^^^^^^^

**lfssh_vscale** = [boolean] **tminimum** = [real]

Default: **lfssh_vscale** = .FALSE.

.TRUE. : Rescale atomic velocities. .FALSE. : Do not rescae atomic
velocities.

Determines whether or not to rescale atomic velocities upon electronic
excitation. **tminimum** gives the minimum temperature in [K] and is
used to constrain velocity scaling.


5.3.9 (time step)
^^^^^^^^^^^^^^^^^

**dttddft** = [real]

Default: **dttddft** = 0.02d0

Gives the time step in [a.u.] for numerically integrating the TDDFT
equations.


5.3.10 (parallel calculation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lfssh_parallel** = [boolean]

Default: **lfssh_parallel** = .TRUE.

.TRUE. : Solves time-dependent K-S equations in parallel .FALSE. :
Solves time-dependent K-S equations serially.

Determines whether or not to perform TDDFT calculations in parallel.


5.3.11 (restart)
^^^^^^^^^^^^^^^^

**ltddft_start** = [boolean]

Default: **ltddft_start** = .FALSE.

.FALSE. : Initialize electronic occupations as specified in the
(occupations) subsection .TRUE. : Initialize electronic occupations with
their values from the last step of a previously completed simulation.
All files beginning with QM_\* and MD_\* must be present in the
**data/** directory.

Determines how to initialize electronic excitations for an NAQMD
simulation. Note that if you are continuing a calculation by setting
**lstart** to .TRUE. in the \*start(on/off), you usually want
**ltddft_start** to be set to .TRUE. as well to properly restart a
successfully completed simulation.


5.3.12 (initial exciton)
^^^^^^^^^^^^^^^^^^^^^^^^

**nexciton** = [integer] **iband_hole, iband_electron, ldegenerate** =
[integer] [integer] [boolean]

Default: **nexciton** = 0

**nexciton** = total number of excitons **iband_hole** = band index of
the hole. **iband_electron** = band index of the electron.
**ldegenerate** = .TRUE. : triplet **ldegenerate** = .FALSE. : singlet

This set of variables are used to define an exciton in linear-reponse
TDDFT. Thus, these variables are only read when **lrtddft** is set to
.TRUE. in the \*linear-response TDDFT section.

**nexciton** gives the total number of excitons to be initially created
and should be followed by that many lines of space delimted values for
**iband_hole, iband_electron, ldegenerate**, where each set specifies
one exciton. **nexciton** should be only be set to 0 or 1, as higher
numbers of excitons is not gaurenteed to work.

Example:

::

   (initial exciton)
       1
       10  11  .FALSE.


5.3.13 (ground state force)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lfssh_gsfrc** = [boolean]

Default: **lfssh_gsfrc** = .FALSE.

.TRUE. : ground state forces are used in FSSH .FALSE. : the excited
state forces are used.

Note that this variable is only read if **lfssh_gsscf** is set to .TRUE.
in the (FSSH-ground-state-SCF) subsection.


5.3.14 (NSC force)
^^^^^^^^^^^^^^^^^^

**ltddft_nscforce** = [boolean]

Default: **ltddft_nscforce** = .TRUE.

.TRUE. : On .FALSE.: Off

Determines whether or not to calculate excited state forces using a
non-self-consistent (NSC) method.


5.3.15 (occupations)
^^^^^^^^^^^^^^^^^^^^

**nocc_change** = [integer] **num_band, occ_new** = [integer] [real]
[real]

**nocc_change** = total number of bands with electronic occupation
number to be changed. **num_band** = band index for occupation change.
**occ_new** = new occupation number for the corresponding band for up
(and down) electrons.

**nocc_change** gives the total number of bands that will undergo a
change in electronic occupation number and should be followed by that
many lines of space delimted values for **num_band, occ_new**, where
each set specifies the band index undergoing a change in its electron
occupation values. If not using spin polarization, only the first number
in the **occ_new** varibale will be read. If using spin-polarization,
you must set **lspin** to .TRUE. in the \*spin polarization section,
then, you may specify the number of up electrons to move with the first
number in the **occ_new** varibale set and the number of down electrons
to move with the second number in the **occ_new** varibale set.

Example for no spin-polarization: This example shows how to excited both
electrons in band index 10 and both electrons in band index 11, to band
indices 12 and 13 (four total electrons changing occupation). This
assumes that in the ground state band indicies 10 and 11 are fully
occupied and band indices 12 and 13 are empty.

::

   (occupations)
   4
   10   0.0   0.0
   11   0.0   0.0
   12   2.0   0.0
   13   2.0   0.0

Example for spin-polarization: This example shows how to excite the up
electron in band index 10 and the down electron in band index 11 to band
indices 12 and 13, respectively (two total electrons changing
occupation). This assumes that in the ground state band indicies 10 and
11 are fully occupied and band indices 12 and 13 are empty.

::

   (occupations)
   2
   10   0.0   1.0
   11   1.0   0.0
   12   1.0   0.0
   13   0.0   1.0


5.3.16 (broadening)
^^^^^^^^^^^^^^^^^^^

**tdbroad** = [real]

Default: **tdbroad** = 0.0

Determines the width of Gaussian broadening of the Fermi surface in [K].
Note: **tdbroad** = 0.0 denotes no broadening.


5.3.17 (DISH)
^^^^^^^^^^^^^

**lfssh_dish** = [boolean] **ndishpair** = [integer] **ndishi, ndishj,
decoherence_rate** = [integer] [integer] [real]

Default: **lfssh_dish** = .FALSE. **ndishpair** = 0

**lfssh_dish** = .TRUE. : Enables Decoherence-Induced Surface Hopping
(DISH) **lfssh_dish** = .FALSE. : Disables DISH **ndishpair** = the
number of state pairs **ndishi, ndishj, decoherence_rate** = the two
band indices between which to define the decoherence rate in [a.u.] for
DISH

Example:

::

   (DISH)
    .true.
      6
      23   24  5.063109E-03
      23   25  5.147713E-03
      23   26  4.596093E-03
      24   25  7.877069E-03
      24   26  7.426337E-03
      25   26  2.768402E-03

The decoherence rate for each pair of states is given by
`[5] <#Jaeger_2012>`__:

::

        rate = sqrt(alpha),

where alpha is a parameter of gaussian

::

        gaus(alpha,t) = exp(-alpha*t*t)

which is fitted to the dephasing function

::

       dij(t) = exp(-gij(t))

with

::

       gij(t) = int^t_0 intg(t') dt'

and

::

       int g(t) = int_t^0 Cij(t)

Cij(t) is an autocorrelation function of the energy gap between two
states

::

       Cij(t) = <(Eij(t)-Eij_ave)(Eij(0)-Eij_ave)>
              = <(Eij(t)*Eij(0)> - Eij_ave*Eij_ave

::

     where Eij(t) = Ei(t) - Ej(t)


Reference:
^^^^^^^^^^

1. "Time-dependent density-functional theory", Gross, E. K. U., and W.
Kohn. Adv. Quantum Chem. 21, 255-291, (1990) 2."Molecular dynamics with
electronic transitions." Tully, John C. J. Chem. Phys. 93.2, 1061-1071
(1990) 3."Convergence acceleration of iterative sequences. the case of
scf iteration." Pulay, P Chem. Phys. Lett. 72.3, 393-398 (1980)
4."Iterative minimization techniques for ab initio total-energy
calculations: molecular dynamics and conjugate gradients" Payne, M. C.,
Teter, M. P., Allan, D. C. ,Arias, T. A., and Joannapoulos, J. D. Rev.
Mod. Phys. 64, 1045 (1992) 5. "Decoherence-induced surface hopping."
Jaeger, H. M., Fischer, S., & Prezhdo, O. V. J. Chem. Phys. 137, 22A545
(2012)


5.4 linear-response TDDFT
~~~~~~~~~~~~~~~~~~~~~~~~~


5.4.1 (On/Off)
^^^^^^^^^^^^^^

**lrtddft** = [boolean]

Default: **lrtddft** = .FALSE.

.TRUE. : Execute NAQMD based on linear response time-dependent density
functional theory (LR-TDDFT) `[1] <#Casida_1995>`__ .FALSE. : Do not use
linear response theory for NAQMD.

Determines whether or not to use linear reponse theory, which involves
calculation the Casida coupling matrix, for NAQMD simulations. For more
detail, please see paper by Casida et al. `[1] <#Casida_1995>`__


5.4.2 (whether to specify states)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lrspecific** = [boolean]

Default: **lrspecific** = .FALSE.

.TRUE. : Specify the occupied and unoccupied electronic states to be
considered in the (specific states) subsection. .FALSE. : Specify the
occupied and unoccupied electronic states to be considered in the
(energy difference) subsection.

Determines how the occupied and unoccupied states to be considered in
the linear response TDDFT calculations. When **lrspecific** is set to
.TRUE. the occupied and unoccupied electronic states should be specified
in the (specific states) subsection by their band index numbers. When
**lrspecific** is set to .FALSE. the occupied and unoccupied electronic
states should be specified in the (energy difference) subsection by
giving the maximum energy difference between them.


5.4.3 (specific states)
^^^^^^^^^^^^^^^^^^^^^^^

**nlrstates** = [boolean]

**ihole, iparticle, ispin, ikpts** = [integer] [integer] [integer]
[integer] OR **ihband, ipband** = [integer] [integer]

Default: **nlrstates** = 0 **ihband** = 0 **ipband** = 0

**nlrstates** = the number of states to specify **ihole** = band index
of hole (occupied states) **iparticle** = band index of particle
(unoccupied states) **ispin** = spin index (1 or 2) **ikpts** = k-point
index **ihband** = band index of hole (occupied states) **ipband** =
band index of particle (unoccupied states)

If **nlrstates** > 0, then states are specified with **ihole, iparticle,
ispin, ikpts**. If **nlrstates** = 0, then states are specified with
**ihband, ipband**, where states between **ihband** and **ipband** will
be considered. Note that this subsection is only read if **lrspecific**
is set to .TRUE. in the (whether to specify states) subsection.


5.4.4 (unit of energy)
^^^^^^^^^^^^^^^^^^^^^^

(ry) or (hr) or (ev)

Default: (ry)

The units of energy to be used for the variable **enediff** given below.


5.4.5 (energy difference)
^^^^^^^^^^^^^^^^^^^^^^^^^

**enediff** = [real]

Default: **enediff** = 0.3

**enediff** = the maximum energy difference in the units given above
between occupied and unoccupied states considered in linear response
TDDFT. Note that this subsection is only read if **lrspecific** is set
to .FALSE. in the (whether to specify states) subsection.


5.4.6 (finite-difference mesh)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**fdmeshxc** = [real]

Default: **fdmeshxc** = 1.d-06

**fdmeshxc** is the finite difference mesh size for the 2nd derivative
of Exc


5.4.7 (threshold)
^^^^^^^^^^^^^^^^^

**threshrpa** = [real]

Default: **threshrpa** = 0.0

**threshrpa** is the threshold value for RPA1


5.4.8 (transition-dipole threshold)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**threshdipole** = [real]

Default: **threshdipole** = 0.0

**threshdipole** threshold value for transition dipole.


5.4.9 (energy-difference threshold)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**threshexcite** = [real]

Default: **threshexcite** = 1.d-04

**threshexcite** is the threshold value for the difference between the
energies of hole and particle


5.4.10 (long-range exchange scheme)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**llcexchange** = [boolean]

Default: **llcexchange** = .FALSE.

.TRUE. : Use long-range correction for exchange. .FALSE. : Do not use.

Determines whether or not to use the long-range correction for exchange
functional `[2] <#Tawada_2004>`__. Note that if **lvacuum[1|2|3]** is
set to .TRUE. **alpha_ldouble_grid_recip** needs to be specified to
determine the ratio of the long-range part.


5.4.11 (parameters to divide 1/r)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**alpha_llcexchange** = [real] **large_cell_llcexchange** = [real]

Default: **alpha_llcexchange** = 0.17 **large_cell_llcexchange** = 2.0

Note 1: The cell-vector length of the CUBIC large cell is given by
large_cell_llcexchange \* max(L1,L2,L3). Note 2: (iosp) & (ecutlong)
have to be set in \*double-grid method only if llcexchange &
.not.ldouble_grid_recip


5.4.12 (hybrid mixing parameter)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**hybrid_mixing** = [real]

Default: **hybrid_mixing** = 0.0

**hybrid_mixing** is the mixing paramter for the Hybrid method. Note
that this is only read when **llcexchange** is set to .FALSE.


5.4.13 (Foerster: width of Gaussian)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**foersigma** = [real]

Default: **foersigma** = 0.1d0

**foersigma** gives the Foerster transfer rate


5.4.14 (diagonal approximation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**ldiagonal** = [boolean]

Default:

**ldiagonal** = .TRUE.

.TRUE. : The diagonal approximation is made, where matrix elements are
calculated only for **ihole_1 == ihole_2** and **iparticle_1 ==
iparticle_2** .FALSE. : Diagonal approximation is not used.

Determines whether or not to use the diagonal approximation. Note: If
**ltddft** is set to .TRUE. .and. **ltddft_fssh** is also set to .TRUE.,
set **ldiagonal** equal to .TRUE.


5.4.15 (scissors correction)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lscissors_lrtddft** = [boolean] **ecorr_triplet** = [real]
**ecorr_singlet** = [real]

Default: **lscissors_lrtddft** = .FALSE.

**lscissors_lrtddft** = .TRUE. : Implement scissors correction.
**lscissors_lrtddft** = .FALSE. : Do no implement correction.

Note: The values for **ecorr_triplet** and **ecorr_singlet** depend on
'unit of energy' in the \*linear-response TDDFT section. These
corrections will be added to the excitation energies.


5.4.16 (singlet-fission-rate calculation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lfission_rate** = [boolean]

Default: **lfission_rate** = .FALSE.

.TRUE. : Calculate the singlet fission rate based on time-dependent
perturbation theory. .FALSE. : Do not calculate rate.


5.4.17 (spontaneous emission rate)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lemission_rate** = [boolean] **refractive** = [real]

Default: **lemission_rate** = .TRUE. **refractive** = 1.5

**lemission_rate** = .TRUE. : Compute spontaneous emission rate based on
transition dipole approximation. **lemission_rate** = .FALSE. : Do not
calculate rate.

Determines whether or not to calculate spontaneous emission rate.
**refractive** is the refractive index of the system.


Reference:
^^^^^^^^^^

1. "Time-Dependent Density Functional Response Theory for Molecules" M.
E. Casida, Recent Advances in Density Functional Methods (Part I),
edited by D. P. Chong (World Scientific, Singapore, pp. 155-192 (1995))
2. "A long-range-corrected time-dependent density functional theory."
Tawada, Yoshihiro, et al. J. Chem. Phys. 120, 8425-8433 (2004)


5.5 PAW
~~~~~~~


5.5.1 (On/Off)
^^^^^^^^^^^^^^

**lpaw** = [boolean]

Default: **lpaw** = .FALSE.

.TRUE. : Use the PAW method `[1] <#Blochl_1994>`__\  .FALSE.: Use other
pseudopotential method.


5.5.2 (non-spherical symmetry)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lpaw_sym** = [boolean]

Default: **lpaw_sym** = .FALSE.

.TRUE. : Full symmetry .FALSE. : Spherical symmetry


5.5.3 (onsite charge mixing)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**paw_mix** = [real]

Default: **paw_mix** = 0.5

Gives the onsite charge mixing parameter in the PAW method.


Reference:
^^^^^^^^^^

1. "Projector augmented-wave method" P. E. Blochl Phys Rev B 50, 17953
(1994)


5.6 Cluster
~~~~~~~~~~~


5.6.1 (On/Off)
^^^^^^^^^^^^^^

**lclust** = [boolean]

Default: **lclust** = .FALSE.

.TRUE. : Perform cluster calculation. .FALSE. : Perform bulk
calculation.

Determines whether to perform a cluster or bulk calculation.


5.7 Approximation for Exchange
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


5.7.1 (approximation)
^^^^^^^^^^^^^^^^^^^^^

**jgga** = 1|2|3|4|5|6

Default: **jgga** = 2

======== ==================================
**jgga** Excahnge functional
======== ==================================
**1**    LD\ `[1] <#Perdew_1981>`__
**2**    GGA(PBE)\ `[2] <#Perdew_1996>`__
**3**    GGA(RPBE)\ `[3] <#Hammer_1999>`__
**4**    GGA(revPBE)\ `[4] <#Zhang_1998>`__
**5**    vdW-DF\ `[5] <#Dion_2004>`__
**6**    vdW-DF2 `[6] <#Lee_2010>`__
======== ==================================

**jgga** is used to specify the exchange correction functional.


5.7.2 (DFT-D)
^^^^^^^^^^^^^

**ldftd** = [boolean]

Default: **ldftd** = .FALSE.

.TRUE. : Employ an empirical vdW correction .FASLE. : Do not use an
empirical vdW correction

Determines whether to use an empirical correction for the van der Waals
interaction proposed by S. Grimme `[7] <#grimme_2006>`__


5.7.3 (kernel factorization) ------ jgga=5 or 6 (vdW-DF or vdW-DF2)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lvdw_factorize** = [boolean]

Default: **lvdw_factorize** = .FALSE.

.TRUE. : Take into account. .FASLE. : Do not

only used for for vdW-DF or vdW-DF2 (jgga=5 or 6 )


5.7.4 (q-mesh range) ------ jgga=5 or 6 (vdW-DF or vdW-DF2)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lvdw_factorize** = [boolean] **q_cut** = [real] **dlambda** = [real]
**nqmesh** = [integer]

Default: **q_cut** = 5.0 **dlambda** = 1.2 **nqmesh** = 20

These parameters are used to saturate the original function q_0(n,|Delta
n|) by redefining q_0^{sat} = h[q_0,q_c].

h(x,x_c) = x_c[1 - exp( - sum_{m=1}^{m_c} (x/x_c)^m/m)]

**q_cut** is corresponding to q_c and **dlambda** and **nqmesh** are a
width and the number of a logarithmic mesh.


5.7.5 (kernel cutoff parameters) ------ jgga=5 or 6 (vdW-DF or vdW-DF2)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**d_soft** = [real] **phi_at_0** = [real]

Default: **d_soft** = 1.0 **phi_at_0** = 0.5

Because the kernel phi(d_1,d_2) has a logarithmic divergence when d_1
and d_2 tend to zero, the modified kernel is used for small d_1 and d_2:

phi(d_1,d_2) = phi_0 + phi_2 d^2 + phi_4 d^4 (d < d_s)

where d = sqrt(d_1^2 + d_2^2) and phi_2 and phi_4 are fitting
parameters. d_s and phi_0 are corresponding to **d_soft** and
**phi_at_0**


5.7.6 (cutoff length) ------ lvdw_factorize = .FALSE.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**rvdwmax** = [real]

Default: **rvdwmax** = 10.d0

Cutoff length in [a.u]


5.7.7 (integration skip) ------ lvdw_factorize = .FALSE.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**imod1** = [even integer] **imod2** = [integer]

Default: **imod1** = 2 **imod2** = 3

Integration of r' on every **imod1** meshes Integration of r on every
**imod2** meshes


5.7.8 (order of cardinal B spline) ------ lvdw_factorize = .FALSE.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**iosp** = [odd integer]

Default: **iosp** = 7

order of cardinal B spline in the double-grid method
`[8] <#Martyna_1999>`__


5.7.9 (DFT+U)
^^^^^^^^^^^^^

**lplusU** = [boolen]

Default: **lplusU** = .FALSE.

.TRUE. : Employ DFT+U correction of the mean-field Hubbard model .FASLE.
: Do not use DFT+U correction Determines whether to use DFT+U correction
of the mean-field Hubbard model `[9] <Dudarev_1998>`__. The value of U
and J parameter is set to in the \*atoms section. Note that this
subsection is only considered when **lpaw** is set to .TRUE.


5.7.10 (onsite charge mixing)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**plusU_mix** = [real]

Default: **plusU_mix** = 0.5

**plusU_mix** is used as turning parameter for onsite charge mixing.
Note that this subsection is only considered when both **lpaw** and
**lplusU** are set to .TRUE.


5.7.11 (DFT+C)
^^^^^^^^^^^^^^

**lplusC** = [boolen]

Default: **lplusC** = .FALSE.

.TRUE. : Employ an empirical correction to non-local pseudopotential
`[10] <#Li_2005>`__ .FASLE. : Do not


5.7.12 (hybrid functional)
^^^^^^^^^^^^^^^^^^^^^^^^^^

**jhybrid** = [0|1]

Default: **jhybrid** = 0

=========== ============================
**jhybrid** Hybrid functional
=========== ============================
**0**       Do not use
**1**       HSE `[11-13] <#Heyd_2003>`__
=========== ============================

**jhybrid** is used to specify the range-separated hybrid exchange
functionals. Note that this subsection is only considered when **jgga**
is set to 2(GGA(PBE)).


5.7.13 (hybrid mixing parameter)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**hmixing** = [real]

Default: **hmixing** = 0.25

The HSE exchange correlation functional is constructed as a linear
combination of the Hartree-Fock exact exchange functional and PBE
functional:

Exc^HSE = a Ex^HF + (1-a) Ex^PBE-SR + Ex^PBE-LR + Ec^PBE

where **a** is defined by **hmixing**.


5.7.14 (range-separation parameter)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**hrange** = [real]

Default: **hrange** = 0.11

The HSE exchange correlation functional spits the Coulomb operator into
short-range and long-range components by an error function and
**hrange**, omega is an adjustable parameter:

1/r = erfc(omega r)/r + erf(omega r)/r


5.7.15 (grid-reduction factor)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**mgridred** = [integer]

Default: **mgridred** = 1

**mgridred** gives the grid-reduction factor. Note: The number of
k-point grids should be the multiple of **mgridred**.


Reference:
^^^^^^^^^^

1. "Self-interaction correction to density-functional approximations for
many-electron systems" J. P. Perdew and Alex Zunger Phys. Rev. B 23 ,
5048(1981) 2. "Generalized Gradient Approximation Made Simple" J. P.
Perdew, K. Burke and M. Ernzerhof Phys. Rev. Lett. 78, 3865 (1996) 3.
"Improved adsorption energetics within density-functional theory using
revised Perdew-Burke-Ernzerhof functionals" B. Hammer, L. B. Hansen, and
J. K. Nørskov Phys. Rev. B 59, 7413 (1999) 4. "Comment on “Generalized
Gradient Approximation Made Simple” " Yingkai Zhang and Weitao Yang
Phys. Rev. Lett. 80, 890 (1998) 5. "Van der Waals Density Functional for
General Geometries " M. Dion, H. Rydberg, E. Schröder, D. C. Langreth,
and B. I. Lundqvist Phys. Rev. Lett. 92, 246401 (2004) 6.
"Higher-accuracy van der Waals density functional" Kyuho Lee, Éamonn D.
Murray, Lingzhu Kong, Bengt I. Lundqvist, and David C. Langreth Phys.
Rev. B 82, 081101 (2010) 7. "Semiempirical GGA‐type density functional
constructed with a long‐range dispersion correction" Stefan Grimme J.
Comput. Chem. 27, 1787 (2006) 8. "A reciprocal space based method for
treating long range interactions in ab initio and force-field-based
calculations in clusters" J. Chem. Phys. 110, 2810 (1990) 9.
"Electron-energy-loss spectra and the structural stability of nickel
oxide:  An LSDA+U study" Phys. Rev. B 57, 1505 (1998) 10.
"Band-structure-corrected local density approximation study of
semiconductor quantum dots and wires" Phys. Rev. B 72, 125325 (2005) 11.
"Hybrid functionals based on a screened Coulomb potential." J. Heyd and
G. E. ScuseriaJ. Chem. Phys. 118, 8207-8215 (2003) 12. "Efficient hybrid
density functional calculations in solids: Assessment of the
Heyd–Scuseria–Ernzerhof screened Coulomb hybrid functional " J. Heyd and
G. E. ScuseriaJ. Chem. Phys. 121, 1187 (2004) 13. "Influence of the
exchange screening parameter on the performance of screened hybrid
functionals " A.V. Krukau, O. A. Vydrov, A. F. Izmaylov, and G. E.
ScuseriaJ. Chem. Phys. 125, 224106 (2006)


5.8 scissors correction at donor/acceptor interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


5.8.1 (On/Off)
^^^^^^^^^^^^^^

**lscissors** = [boolean]

Default: **lscissors** = .FALSE.

.TRUE. : Take into account .FASLE. : Do not


5.8.2 (unit of energy)
^^^^^^^^^^^^^^^^^^^^^^

(ry) or (hr) or (ev)


5.8.3 (corrections)
^^^^^^^^^^^^^^^^^^^

**ecorr_acceptor** = [real] **ecorr_donor** = [real]

Default **ecorr_acceptor** = 0.0 **ecorr_donor** = 0.0

Corrections to be subtracted from the eigenenergies.


5.8.4 (unit of length)
^^^^^^^^^^^^^^^^^^^^^^

(bohr) or (ang)


5.8.5 (boundary)
^^^^^^^^^^^^^^^^

**zboundary** = [real]

Default: **zboundary** = 0.0 Defines the boundary between acceptor and
donor regions as follows: z-coordinate < **zboundary**: acceptor region
z-coordinate > **zboundary**: donor region


5.9 k-points
~~~~~~~~~~~~


5.9.1 (number of k points)
^^^^^^^^^^^^^^^^^^^^^^^^^^

**nkpnt** = [integer]

Default: (gamma point) **nkpnt** = 1

**nkpnt** is the number of k points. If **nkpnt** is greater than zero,
you set k vectors in (k vectors). If **nkpnt** is zero, you set the
number of division for the reciprocal vectors in (division number).


5.9.2 (k vectors)
^^^^^^^^^^^^^^^^^

**(bzk,wbzk)** = [real] [real] [real] [real]

Default: (gamma point) **bzk(1:3)** = (/0.0 0.0 0.0/) **wbzk** = 1.0

**bzk(1:3)** are k vectors and **wbzk** is the weight of each k-vector.
Note that this subsection is only considered when **nkpnt** is greater
than zero.Also, Please set K vectors between -0.5 and 0.5.


5.9.3 (division number)
^^^^^^^^^^^^^^^^^^^^^^^

**npkx(1:3)** = [integer] [integer] [integer]

Default: **npkx(1:3)** = (/ 1 1 1 /)

**npkx(1:3)** is the number of division for reciprocal vectors and
Monkhorst-Pack method is used. `[1] <#Monkhorst_1976>`__


5.9.4 (tetrahedron method)
^^^^^^^^^^^^^^^^^^^^^^^^^^

**ltetra** = [boolean]

Default: **ltetra** = .FALSE.

.TRUE. : Use the tetrahedron method .FALSE. : Do not use

Determines whether to use the modified tetrahedron method
`[2] <#Blochl_1994>`__ for the Brillouin Zone integration.


Reference:
^^^^^^^^^^

1. "Special points for Brillouin-zone integrations" H. J. Monkhorst and
J. D. Pack, O. Jepsen, and O. K. Andersen. Phys. Rev. B 13, 5188 (1976)
2. "Improved tetrahedron method for Brillouin-zone integrations." P. E.
Blöchl, O. Jepsen, and O. K. Andersen. Phys. Rev. B 49, 16223 (1994)


5.10 electric field
~~~~~~~~~~~~~~~~~~~


5.10.1 (On/Off)
^^^^^^^^^^^^^^^

**lefield** = [boolean]

Default: **lefield** = .FALSE.

.TRUE. : Apply uniform electric field. .FALSE. : No electric field.


5.10.2 (electric field)
^^^^^^^^^^^^^^^^^^^^^^^

**efield(1:3)** = [real] [real] [real]

Determines the uniform electric field in atomic units.
`[1] <#Umari_2001>`__


Reference:
^^^^^^^^^^

19. "Ab initio molecular dynamics in a finite homogeneous electric
field." P. Umari, and A. Pasquarello. Phys. Rev. Lett. 89,157602 (2002)


5.11 spin polarization
~~~~~~~~~~~~~~~~~~~~~~


5.11.1 (On/Off)
^^^^^^^^^^^^^^^

**lspin** = [boolean]

Default: **lspin** = .FALSE.

.TRUE. : Apply spin polarization .FALSE. : Do not use

Determines whether to take into account spin polarization.


5.11.2 (noncollinear magnetism)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lnoncollinear** = [boolean]

Default: **lnoncollinear** = .FALSE.

.TRUE. : Noncollinear magnetism .FALSE. : Collinear magnetism

Note that the DFT calculation by noncollinear magnetism is still under
development (SCF iterations might not converge).


5.11.3 (fix spin polarization)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lfixud** = [boolean]

Default: **lfixud** = .TRUE.

.TRUE. : Fix the value of spin polarization .FALSE. : Do not fix

**lfixud** set to .TRUE. will be fixed the value of spin polarization.
Note that this subsection is only considered when **lnoncollinear** is
.FALSE.


5.11.4 (spin polarization)
^^^^^^^^^^^^^^^^^^^^^^^^^^

**diffud** = [real]

Default: **diffud** = 0.0

**diffud** is the value of spin polarization to fix total magnetic
moment. Note that this subsection is only considered when **lfixud** is
.TRUE.


5.11.5 (initial wave functions)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lwfrand** = [boolean]

Default: **lfixud** = .FALSE.

.FALSE. : Down spin is given as same as up spin .TRUE. : Down spin is
given by random number

Determines how to give the down spin.


5.11.6 (initial spin density)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**inispin** = [0|1|2|3|4]

Default: **inispin** = 1

=========== ==================================
**inispin** Initial spin density
=========== ==================================
**0**       from LSDA (lnoncollinear = .TRUE.)
**1**       uniformly polarized
**2**       ferromagnetic alignment
**3**       antiferromagnetic random alignment
**4**       specify spin polarization par atom
=========== ==================================

**inispin** is used to specify the method for initial spin density. The
detail parameter is set in another subsection.


5.11.7 (uniform spin density) ------ **inispin** = 1 & **lnoncollinear** = .TRUE.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**atmmagne** = [real] **umx**, **umy**, **umz** = [real] [real] [real]

If **inispin** is one, spin polarized uniformly. **atmmagne** is total
magnetic moment and (**umx**, **umy**, **umz**) are a magnetization
direction.


5.11.8 (ferromagnetic spin density) ------ **inispin** = 2 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**nspinat** = [integer] **iatspin** = [integer] **atmmagne** = [real]
(**umx**, **umy**, **umz**) = [real] [real] [real]

Default: **nspinat** = 1

All **itspin**-type atoms are polarized with the magnetic moment of
**atmmagne**. (**umx**, **umy**, **umz**) are available only when
**lnoncollinear** is .TRUE. **nspinat** is the number of atom types to
be polarized.

Example:

::

   (ferromagnetic spin density)
      2                          : (nspinat)
      1    2.d0  0.d0 0.d0 1.d0  : (iatspin, atmmagne, umx, umy, umz)
      2   -4.d0  0.d0 0.d0 1.d0  : (iatspin, atmmagne, umx, umy, umz)


5.11.9 (antiferromagnetic spin density) ------ **inispin** = 3 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**nspinat** = [integer] **iatspin** = [integer] **atmmagne** = [real]
(**umx**, **umy**, **umz**) = [real] [real] [real]

Default: **nspinat** = 1

Half of **itspin**-type atoms are polarized with the magnetic moment of
**atmmagne**, and the rest are polarized with the magnetic moment of
**-atmmagne**. (**umx**, **umy**, **umz**) are available only when
**lnoncollinear** is .TRUE. **nspinat** is the number of atom types to
be polarized.

Example:

::

   (antiferromagnetic spin density)
      2                          : (nspinat)
      1    2.d0  0.d0 0.d0 1.d0  : (iatspin, atmmagne, umx, umy, umz)
      2   -4.d0  0.d0 0.d0 1.d0  : (iatspin, atmmagne, umx, umy, umz)


5.11.10 (atomic spin density) ------ **inispin** = 4
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**nspinat** = [integer] **iatspin** = [integer] **atmmagne** = [real]
(**umx**, **umy**, **umz**) = [real] [real] [real]

Default: **nspinat** = 1

You can specify the magnetic moment of **atmmagne** par atom.
**iatspin** is the index of atom number. (**umx**, **umy**, **umz**) are
available only when **lnoncollinear** is .TRUE. **nspinat** is the
number of atom types to be polarized.

Example:

::

   (atomic spin density)
      4                          : (nspinat)
     10    2.d0  0.d0 0.d0 1.d0  : (iatspin, atmmagne, umx, umy, umz)
     11    2.d0  0.d0 0.d0 1.d0  : (iatspin, atmmagne, umx, umy, umz)
     23    2.d0  0.d0 0.d0 1.d0  : (iatspin, atmmagne, umx, umy, umz)
     24    2.d0  0.d0 0.d0 1.d0  : (iatspin, atmmagne, umx, umy, umz)


5.11.11 (reference direction) ------ **lnoncollinear** = .TRUE.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**refmx**, **refmy**, **refmz** = [real] [real] [real]

Default: **refmx**, **refmy**, **refmz** = (/ 0.0 0.0 1.0/)

When **lnoncollinear** is .TRUE., **refmx**, **refmy**, **refmz** shows
the reference direction.


5.11.12 (duplicate the up-spin state) ------ **lnoncollinear** = .FALSE.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lduplicate** = [boolean]

Default: **lduplicate** = .FALSE.

.TRUE. : Duplicate wavefunction of the up spin state .FALSE. : Do not
duplicate


5.12 SCF iterations
~~~~~~~~~~~~~~~~~~~


5.12.1 (global iterations)
^^^^^^^^^^^^^^^^^^^^^^^^^^

**iscfmx** = [integer]

Default: **iscfmx** = 100

Gives the maximum number of global iterations of SCF to be performed.


5.12.2 (trial global iterations)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**itrial** = [integer]

Default: **itrial** = 0

Gives the number of trial global iterations


5.12.3 (trial charge mixing)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**imxchg_trial** = [0|1] **aslh_trial, bslh_trial** = [real] [real]

Default: **imxchg_trial** = 0 **aslh_trial** = 0.1 **bslh_trial** = 0.64

**imxchg_trial** = 0 : No trial charge mixing **imxchg_trial** = 1 :
Trial charge mixing

Determines whether to perform trial charge mixing, and if so, the tuning
parameters, **aslh_trial, bslh_trial**, to be used.


5.12.4 (trial spin-density mixing)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**imxspin_trial** = [0|1] **amxspin_trial, bmxspin_trial** = [real]
[real]

Default: **imxspin_trial** = 0 **amxspin_trial** = 0.1 **bmxspin_trial**
= 0.64


5.12.5 (tolerances)
^^^^^^^^^^^^^^^^^^^

**tolpot** = [real] **tolres** = [real]

Default: **tolpot** = 5.0d-09 **tolres** = 5.0d-09

**tolpot** gives the tolerance for change in total energy. Once the
change in total energy between subsequent SCF iterations is smaller than
the tolerance, the SCF iterations are considered to have converged.
**tolres** similarly gives the tolerance for average residual.


5.12.6 (HC products)
^^^^^^^^^^^^^^^^^^^^

**lhcunt** = [boolean]

Default: **lhcunt** = .TRUE.

.TRUE. : HC products by Unitary tr. .FALSE. : ???

???


5.12.7 (charge mixing)
^^^^^^^^^^^^^^^^^^^^^^

**imxchg** = [0|1|2|3|4|5|6] **aslh, bslh** = [real] [real]

Default: **imxchg** = 1

====== ===========================
imxchg Charge Mixing Method
====== ===========================
**0**  No Mixing
**1**  Pulay
**2**  Anderson
**3**  Simple
**4**  Srivastava
**5**  Johnson
**6**  Johnson w/ variable weights
====== ===========================

Determines which charge mixing method will be used, along with tuning
parameters for that method.


5.12.8 (number of mixing)
^^^^^^^^^^^^^^^^^^^^^^^^^

**itratn** = [integer]

Default: **itratn** = 10

Determines how many charge densities from previous SCF iterations to use
for charge mixing. Not available for imxchg = [0|3].


5.12.9 (metric in charge mixing)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**cmetric** = [real]

Default: **cmetric** = 20.0

Gives the metric for charge mixing. Used as the weighting factor w(g) in
the calculation of residual is given by:

w(g) = (cmetric-1)*G_min^2*\ G_max^2/(G_max^2 - cmetric*G_min^2)

Not available for imxchg = [0|3].


5.12.10 (magnetic moment mixing)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lmixreal** = [boolean]

Default: **lmixreal** = .TRUE.

.TRUE. :real-space mixing .FALSE. : reciprocal-spcace mixing

This is only read when **lnoncollinear** in the (noncollinear magnetism)
subsection of the spin polarization section is set to .TRUE.


5.12.11 (spin-density mixing)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**imxspin** = [0|1|2|3|4|5|6] **amxspin, bmxspin** = [real] [real]

Default: **imxspin** = 0 **amxspin** = 0.8 **bmxspin** = 0.64

======= ===========================
imxspin Charge Mixing Method
======= ===========================
**0**   No Mixing
**1**   Pulay
**2**   Anderson
**3**   Simple
**4**   Srivastava
**5**   Johnson
**6**   Johnson w/ variable weights
======= ===========================

This is only calculated when spin polarization is enabled in the \*spin
polarization section.


5.12.12 (# of spin-density mixing)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**nmxspin** = [integer]

Default: **nmxspin** = 10

Determines how many spin densities from previous SCF iterations to use
for spin density mixing.


5.12.13 (spin metric in spin-density mixing)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**spinmetric** = [real]

Default: **spinmetric** = 1.0

The metric in spin-density mixing. This variable is only read when
**lnoncollinear** == .FALSE. or **lmixreal** == .FALSE.


5.12.14 (Kerker-mixing weight in spin-density mixing)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**wkerker** = [real]

Default: **wkerker** = 1.0

Kerker-mixing weight in spin-density mixing. This variable is only read
when **lnoncollinear** == .FALSE. or **lmixreal** == .FALSE.

.. warning:: Sigma_new(G) = sigma_in(G) + K(G)\ ``*(sigma_out(G) -sigma_in(G)), where sigma is the spin density, and K(G) = amxspin*\ (1 - wkerker + wkerker*G~2/(G^2 + bmxspin)). Note that sigma(G=0) is not updated when wkerker = 1.d0, i.e., the difference between the numbers of up- and down-spin is fixed.`` 

Not available for imxchg = [0|3].


5.12.15 (output eigenvalues)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**louteig** = [boolean]

Default: **louteig** = .FALSE.

.TRUE. : Output eigenvalues to file .FALSE. : Do not output eigenvalues

Whether or not to output eigenvalues through SCF iterations.


5.12.16 (output energy parts)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**loutenergy** = [boolean]

Default: **loutenergy** = .FALSE.

.TRUE. : Output energy parts to file .FALSE. : Do not output energy
parts

Whether or not to output energy parts through SCF iterations


5.12.17 (output residuals)
^^^^^^^^^^^^^^^^^^^^^^^^^^

**loutzansa** = [boolean]

Default: **loutzansa** = .FALSE.

.TRUE. : Output residuals (standard deviation) to file .FALSE. : Do not
output residuals.

Whether or not to output residuals through SCF iterations


5.12.18 (convergence stabilizer)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**nstabi** = [integer] **xstabi** = [real]

Default: **nstabi** = 10 **xstabi** = 10.0

If zansa2 < zansa2/xstabi at nstabi iterations before, set lreset =
.FALSE.


5.13 well potenial
~~~~~~~~~~~~~~~~~~


5.13.1 (On/Off)
^^^^^^^^^^^^^^^

**lwell** = [boolean]

Default: **lwell** = .FALSE.

Determines whether to use a well potential outide atomic sphere.


5.13.2 (height of potential)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**wellheight** = [real]

Default: **wellheight** = 1.0

Gives the height of the well potenital in Rydbergs.


5.13.3 (radius around each atom)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**nwellat** = [integer] **iatwell, radwell** = [integer] [real]

Default: **nwellat** = 0

**nwellat** gives the total number of atom types for which to make a
well potential around, followed by that many lines of **iatwell**,
giving the atom type index and **radwell**, giving the radius around
each atom in bohrs. If not specified (atom type 2 in this case), the
covalent radius is used.

Example: This defines the radius around atom types 1,3, a total of 2
atom types, to be 1.9 bohrs. Atom type 2, for which no radius is given,
will have its covalent radius used.

::

   (radius around each atom)
      2
      1    1.9d0
      3    1.9d0


5.14 Kohn-Sham equation
~~~~~~~~~~~~~~~~~~~~~~~


5.14.1 (On/Off)
^^^^^^^^^^^^^^^

**ihldam** = [1|2]

Default: **ihldam** = 1

1 : Conjugate-Gradient (CG) Method 2 : Residual minimization scheme,
direct inversion in the iterative subspace (RMM-DIIS) Method

Determines the method used in the iterative minimization of energy as a
functional of wavefunctions.


5.14.2 (tolerance)
^^^^^^^^^^^^^^^^^^

**toleig** = [real]

Default: **toleig** = 1.d-10

Gives the tolerance for the Kohn-Sham equations.


5.14.3 (threshold for w.f. direction)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**threinn** = [real]

Default: **threinn** = 0.0

Gives the threshold for direction of new wavefunction


5.14.4 (iteration)
^^^^^^^^^^^^^^^^^^

**itermx** = [integer]

Default: **itermx** = 4

Gives the maximum number of iterations.


5.14.5 (empty-band iteration)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**iteremx** = [integer]

Default: **iteremx** = 4

Gives the maximum number of iterations for empty bands.


5.14.6 (trial global iterations)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**kstrial** = [integer]

Default: **kstrial** = 0

Gives the number of trial global iterations with itermx = iteremx = 1


5.14.7 (CG method)
^^^^^^^^^^^^^^^^^^

**methodcg** = [1|2]

Default: **methodcg** = 2

1:line minimization 2:BKL [ref]

Gives the method used of CG This variable is only read if **ihldam** ==
1


5.15 Poisson equation
~~~~~~~~~~~~~~~~~~~~~

Note this section is only used for cluster calculations.


5.15.1 (multigrid level)
^^^^^^^^^^^^^^^^^^^^^^^^

**multg** = [integer]

Default: **multg** = 2

??


5.15.2 (tolerance)
^^^^^^^^^^^^^^^^^^

**tolcg** = [real]

Default: **tolcg** = 1.d-11

Gives the tolerance energy in a.u. for one electron


5.15.3 (preconditioner)
^^^^^^^^^^^^^^^^^^^^^^^

**weigrd** = [real]

Default: **weigrd** = 0.4

???


5.15.4 (differentiation)
^^^^^^^^^^^^^^^^^^^^^^^^

**nd2v** = [integer]

Default: **nd2v** = 6

Gives the order of numerical differentiation.


5.15.5 (mesh for serial calculation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**msrhmx, msrhmy, msrhmz** = [integer] [integer] [integer]

Default: **msrhmx, msrhmy, msrhmz** = 14 14 14

Gives the mesh size for serial calculation.


5.16 molecular dynamics
~~~~~~~~~~~~~~~~~~~~~~~


5.16.1 (On/Off)
^^^^^^^^^^^^^^^

**ifmd** = [0|1|2|3|4|5|10]

Default: **ifmd** = 0

====== =================
ifmd   Type of Dynamics
====== =================
**0**  Single Step
**1**  Optimization
**2**  NVE
**3**  NVT
**4**  NPT
**5**  NVT for each atom
**10** MSST [ref]
====== =================

Determines the type of QMD simulation to run. Add short description of
each type.


5.16.2 (time step)
^^^^^^^^^^^^^^^^^^

**dtmd, nstop** = [real] [integer]

Default: **dtmd** = 50.0 **nstop** = 10

**dtmd** gives the time step for QMD simulation in [a.u.], while
**nstop** gives the total number of time steps to simulate. Thus, the
total simulation time will equal (**nstop** \* **dtmd**)


5.16.3 (initial step number)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**nstep_ini** = [integer]

Default: **nstep_ini** = 0

Gives the initial step number. This varibale is ignored for **lstart**
== .TRUE.


5.16.4 (temperature)
^^^^^^^^^^^^^^^^^^^^

**treq** = [real]

Default: **treq** = 300.0

Gives the initial temperature in [K] for NVE, NVT, and NPT QMD
simulations.


5.16.5 (check temperature)
^^^^^^^^^^^^^^^^^^^^^^^^^^

**liscale** = [boolean] **iscnum** = [integer] **iscstp** = [integer]

Default: **liscale** = .FALSE. **iscnum** = 25 **iscstp** = 20

**liscale** = .FALSE. : Do not scale temperature. **liscale** = .TRUE. :
Scale temperature a total of **iscstp** times, with **iscnum** steps in
between each scaling.

Determines whether and how to scale temperature to keep it near the
initial given temperature.


5.16.6 (make total momentum zero)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lmomzero** = [boolean]

Deafault: **lmomzero** = .FALSE.

.TRUE. : Make the total momentum zero. .FALSE. : Do not


5.16.7 (optimization)
^^^^^^^^^^^^^^^^^^^^^

**ioptmze** = [-1|0|1|2|3|10]

Default: **ioptmze** = 2

======= =====================================
ioptmze Type of Structural Optimization
======= =====================================
**-1**  Do not optimize atomic coords
**0**   Conjugate gradient
**1**   Projected velocity Verlet
**2**   Quasi-Newton method with BFGS formula
**3**   RFO saddle points search
**10**  Harmonic-mode analysis
======= =====================================

Determines method for structural optimization of atomic coordinates.
This varibale is onlt read when **ifmd** == 1.


5.16.8 (cell optimization)
^^^^^^^^^^^^^^^^^^^^^^^^^^

**ioptmze_cell** = [-1|0|1|2]

Default: **ioptmze_cell** = -1

============ =====================================
ioptmze_cell Type of Cell Optimization
============ =====================================
**-1**       Do not optimize supercell
**0**        Conjugate gradient
**1**        Not yet implemented
**2**        Quasi-Newton method with BFGS formula
============ =====================================

Determines method for optimization of (super)cell. This varibale is only
read when **ifmd** == 1, and is not not read when **ioptmze** == 1 or
10.


5.16.9 (cell CG time step)
^^^^^^^^^^^^^^^^^^^^^^^^^^

**dtcellcg** = [real]

Default: **dtcellcg** = 0.1

only for Conjugate gradient method (**ifmd** == 1 & **ioptmze_cell** ==
0).


5.16.10 (stabilizer for quasi-Newton)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**gammamin** = [real]

Default: **gammamin** = 0.1

only for quasi-Newton method (**ifmd** == 1 & **ioptmze** == 2).


5.16.11 (clear Hessian)
^^^^^^^^^^^^^^^^^^^^^^^

**ibfgsclear** = [0|1]

Default:

**ibfgsclear** = 0

1: Clear Hessian after every ibfgsclear step 0: Hessian is not cleared


5.16.12 (clear cell Hessian)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**iclearcellh** = [integer]

Default: **iclearcellh** = 0

Clear Hessian every **iclearcellh** step. If **iclearcellh** == 0, the
Hessian is not cleared.


5.16.13 (hybrid optimization)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lhybridopt** = [boolean] **nstep_hybrid** = [integer]
**nstep_hybrid_cell** = [integer]

Default: **lhybridopt** = .TRUE. **nstep_hybrid** = 10
**nstep_hybrid_cell** = 10

**lhybridopt** = .TRUE. : Perform structural optimization first, then
cell optimization. **lhybridopt** = .FALSE. : Perform cell optimization
first, then structural optimization.

**nstep_hybrid** is the time step for structural optimization.
**nstep_hybrid_cell** is the time step for cell optimization. These
variables are only read for optimization calculations(\ **ifmd** == 1 &
**ioptmze** >= 0 & **ioptmze_cell** >= 0 )


5.16.14 (displacement)
^^^^^^^^^^^^^^^^^^^^^^

**hmadisp** = [real]

Default: **hmadisp** = 0.01

Atomic displacement in [a.u.] for Harmonic analysis by finite
difference. This is only used for Harmonic-mode analysis (**ifmd** == 1
& **ioptmze** == 10).


5.16.15 (order of differentiation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**nhmaord** = [integer]

Default: **nhmaord** = 1

Order of differentiation only defined for Harmonic-mode analysis
(**ifmd** == 1 & **ioptmze** == 10).


5.16.16 (thermostat parameters)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**nnos** = [integer] **nresn** = [integer] **nyosh** = [integer]
**tomega** = [real]

Default: **nnos** = 1 **nresn** = 1 **nyosh** = 1 **tomega** = 5500.0

[add references] Only for NVT & NPT-MD (**ifmd** == 3, 4) Notes for
(tomega): The frequency for the thermostat will be given by omega [a.u.]
= 2*pi / tomega. The tomega will be determined from the period of VAF.
If you are not aware of it, try a value, 100 \* dtmd in [a.u.], as an
initial guess.


5.16.17 (atom type with thermostat)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**nathermo** = [integer] **ithermo** = [integer]

Only used for NVT for each atom (**ifmd** == 5)

**nathermo** gives the total number of atom types to be given a
thermostat, followed by that many lines of **ithermo**, which specifies
the atom type.

Example:

::

      2  : (nathermo) # of atom types with thermostat
      1  : (ithermo) atom type
      4  : (ithermo) atom type


5.16.18 (atom thermostat parameters)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**nnos_atom** = [integer] **nresn_atom** = [integer] **nyosh_atom** =
[integer] **tomega_atom** = [real]

Default: **nnos_atom** = 1 **nresn_atom** = 1 **nyosh_atom** = 1
**tomega_atom** = 5500.0

**nnos_atom** gives the number of thermostats, **nresn_atom** gives the
MTS steps for heat bath, **nyosh_atom** gives the Yoshida-Suzuki
decomposition step, and **tomega_atom** a period of fluctuation in
[a.u.]. This is only used for NVT for each atom (**ifmd** == 5).


5.16.19 (pressure)
^^^^^^^^^^^^^^^^^^

**hpext** = [real]

Default: **hpext** = 0.0

Defines the pressure in [GPa] for NPT-MD & MSST (**ifmd** == 4, 10) &
**ioptmze_cell** >= 0.


5.16.20 (barostat parameters)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**tbomega** = [real] **blkmod** = [real]

Default: **tbomega** = 5500.0 **blkmod** = 250.0

Barostat parameters for NPT (**ifmd** == 4). **tbomega** gives the time
scale for barostat in [a.u.], and **blkmod** gives the bulk modulus in
[GPa].


5.16.21 (restriction for MD cell)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**irstrct** = [0|1|2|3|4|5|6|10|11|12|13|14|15|16]

Default: **irstrct** = 0

======== ============================================================================
irstrct   Type of Restriction on cell
======== ============================================================================
**0**      No restriction
**1**     **Cubic**: a=b=c; alpha=beta=gamma=90
**2**     **orthorhombic**: a /= b/= c; alpha=beta=gamma=90
**3**     **Tetragonal**: a = b /= c; alpha=beta=gamma=90
**4**     **Monoclinic**: a /= b /= c; alpha=beta=90; gamma /= 90
**5**     **Hexagonal**: a = b /= c; alpha=beta=90; gamma=120
**6**     **Trigonal**: a=b=c; alpha = beta = gamma /= 90
**10**   symmetry opearations (enabled automatically when lsymop=.true. & irstrct==0)
======== ============================================================================
======== ============================================================================

.. note:: Add 10 to the above values to enforce the lattice shape with symmetry opearation. For example, set **irstrct** == 12 for an orthorhombic cell with symmetry opearations.

Set restrictions for geometry of MD cell in NPT-MD simulations.


5.16.22 (sub-restriction for MD cell)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**irstrct_sub** = [integer]

Default: **irstrct_sub** = 1

======= =========== ===========================================
irstrct irstrct_sub Type of Sub-restriction on cell
======= =========== ===========================================
**3**   **1**       a = b /= c, alpha = beta = gamma = 90
\       **2**       b = c /= a, alpha = beta = gamma = 90
\       **other**   c = a /= b, alpha = beta = gamma = 90
**4**   **1**       a /= b /= c, alpha = beta = 90; gamma /= 90
\       **2**       a /= b /= c, beta = gamma = 90; alpha /= 90
\       **other**   a /= b /= c, gamma = alpha = 90; beta /= 90
**5**   **1**       a = b /= c, alpha = beta = 90; gamma = 120
\       **2**       b = c /= a, beta = gamma = 90; alpha = 120
\       **other**   c = a /= b, gamma = alpha = 90; beta = 120
======= =========== ===========================================

Defines sub-options for restriction for MD cell in NPT-MD simulations.


5.16.23 (MD cell edge restriction)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lcell_rstrct(1)** = [boolean] **lcell_rstrct(2)** = [boolean]
**lcell_rstrct(3)** = [boolean]

.TRUE. = fix the length of L[1,2,3] .FALSE. = do not fix length

Determines whether or not to fix the lengths of MD cell edges.

.. _51624-shock-wave-velocity:

5.16.24 (shock wave velocity)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**shockspeed** = [real] **nshockv(1:3)** = [0|1, 0|1, 0|1]

Default: **shockspeed** = 2000.0 **nshockv(1:3)** = (/1,0,0/)

These variables are only read for MSST (**ifmd** == 10). They define the
shockspeed and the shock direction. For example:

**nshockv(1:3)** = 1 0 0: L1 direction **nshockv(1:3)** = 0 1 0: L2
direction **nshockv(1:3)** = 0 0 1: L3 direction **nshockv(1:3)** = 1 1
0: L1 + L2 direction

.. _51625-clear-barostat-velocity:

5.16.25 (clear barostat velocity)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lmsstscale** = [boolean] **msstscnum** = [integer] **msstscstp** =
[integer]

Default: **lmsstscale** = .FALSE. **msstscnum** = 20 **msstscstp** = 500

These variables are only read for MSST (**ifmd** == 10).

.. _51626-initial-barostat-velocity:

5.16.26 (initial barostat velocity)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**vxmsst** = [real]

Default: **vxmsst** = 0.0

Defines the initial barostat velocity. Only read for MSST (**ifmd** ==
10).

.. _51627-output-data:

5.16.27 (output data)
^^^^^^^^^^^^^^^^^^^^^

**ioskip** = [integer] **locoor, ioskipcoor** = [boolean, integer]

Defines what data to output for MD-nodes.

.. _51628-charge-estimation:

5.16.28 (charge estimation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**ichest** = [integer]

Default: **ichest** = 3

Number of previous steps

.. _51629-wavefunction-estimation:

5.16.29 (wavefunction estimation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**ihest** = [integer]

Number of previous steps

.. _51630-aspc-charge-estimation:

5.16.30 (ASPC charge estimation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**laspc_chg** = [boolean]

.TRUE. = ASPC predictor .FALSE. = usual extrapolation

.. _51631-aspc-wavefunction-estimation:

5.16.31 (ASPC wavefunction estimation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**laspc_wv** = [boolean]

.TRUE. = ASPC predictor .FALSE. = subspace alignment

.. _51632-aspc-acceleration:

5.16.32 (ASPC acceleration)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**laspc** = [boolean]

.TRUE. = do ASPC

.. _51633-extended-lagrangian-scheme:

5.16.33 (Extended Lagrangian scheme)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lxlbomd** = [boolean]

.TRUE. : time reversible XL-BOMD

.. _51634-previous-values:

5.16.34 (previous values)
^^^^^^^^^^^^^^^^^^^^^^^^^

**kxlbomd** = [integer]

Number of previous values. If 0, kxlbomd = ihest.

.. _51635-aspc-parameter-for-corrector:

5.16.35 (ASPC parameter for corrector)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**aspc_corr** = [real]

Default: **aspc_corr** = 0.0

.. _51636-aspc-mixing-charge:

5.16.36 (ASPC mixing charge)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**aslh_aspc, bslh_aspc** = [real, real]

.. _51637-aspc-iterations:

5.16.37 (ASPC iterations)
^^^^^^^^^^^^^^^^^^^^^^^^^

**iscfmx_aspc** = [integer]

Gives the number of iterations.

.. _51638-aspc-tolerance:

5.16.38 (ASPC tolerance)
^^^^^^^^^^^^^^^^^^^^^^^^

**tolres_aspc** = [real]

Gives the tolerance for the residual.

.. _51639-aspc-correction-skip-step:

5.16.39 (ASPC correction skip step)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**nskp_aspc** = [integer]

.. _51640-planewave-expansion:

5.16.40 (planewave expansion)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**pwscale** = [real]

only for NPT-MD (ifmd == 4 )

.. _51641-statistical-calculation:

5.16.41 (statistical calculation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lstat** = [boolean]

.. _51642-tolerance:

5.16.42 (tolerance)
^^^^^^^^^^^^^^^^^^^

**tol_energy** = [real] **tol_force** = [real]

tolerance for CG optimization (ifmd == 1 ). **tol_energy** given in
energy/atom in [a.u.]. **tol_force** gives maximum force in [a.u.].

.. _517-save-data:

5.17 save data
~~~~~~~~~~~~~~

.. _5171-onoff:

5.17.1 (On/Off)
^^^^^^^^^^^^^^^

**lsave** = [boolean]

Default: **lsave** = .TRUE.

Whether or not to save data...which data?

.. _5172-data-type:

5.17.2 (data type)
^^^^^^^^^^^^^^^^^^

**lsreal8** = [boolean]

Default: **lsreal8** = .TRUE.

.TRUE. : Saves data as real\ *8 .FALSE. : Saves data as integer*\ 2

Determines the data type for saved data

.. _518-soft-walls:

5.18 soft walls
~~~~~~~~~~~~~~~

.. _5181-walls:

5.18.1 (walls)
^^^^^^^^^^^^^^

**nwalls** = [integer] **wallp** = [real, real, real] **wallv** = [real,
real, real] **wallf** = [real] **nwallp** = [integer]

Default: **nwalls** = 0

**nwalls** gives the number of soft walls to be inserted into the
system, followed by that many sets of **wallp**, **wallv**, **wallf**,
**nwallp**, which define each soft wall. **wallp** gives any point on
the wall in [a.u.]. **wallv** gives the normal vector to the wall,
**wallf** gives the strength of repulsive potential, and **nwallp**
gives the power of the repulsive potential.

Example: This defines two soft walls

::

   (walls)
      2
    10.0, 0.0, 0.0
     1.0, 0.0, 0.0
     1.d-04
     12
    0.0, 10.0, 0.0
    0.0,  1.0, 0.0
     1.d-04
     12

.. _519-gravitational-field:

5.19 gravitational field
~~~~~~~~~~~~~~~~~~~~~~~~

.. _5191-onoff:

5.19.1 (On/Off)
^^^^^^^^^^^^^^^

**lgravi** = [boolean]

Default: **lgravi** = .FALSE.

Determines whether to apply a gravitational field.

.. _5192-magnitude:

5.19.2 (magnitude)
^^^^^^^^^^^^^^^^^^

**gravmag** = [real]

Default: **gravmag** = 0.0

Determines the stength of the gravitational field in units of g = 9.8
[m/s^2]

.. _5193-direction:

5.19.3 (direction)
^^^^^^^^^^^^^^^^^^

**gravdir** = [real, real, real]

Default: **gravdir** = [0.0, 0.0, 0.0]

Determines the direction of the gravitational field.

.. _520-mulliken-analysis:

5.20 Mulliken analysis
~~~~~~~~~~~~~~~~~~~~~~

.. _52021-onoff:

5.20.21 (On/Off)
^^^^^^^^^^^^^^^^

**lmulken** = [boolean]

.. _52022-skip-step:

5.20.22 (skip step)
^^^^^^^^^^^^^^^^^^^

**nskpmulk** = [integer]

.. _52023-decompose-overlap-population:

5.20.23 (decompose overlap population)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**ldecmpovp** = [boolean]

.. _52024-weights-associated-with-each-atom:

5.20.24 (weights associated with each atom)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lspdatom** = [boolean]

.. _52025-bond-charge:

5.20.25 (bond charge)
^^^^^^^^^^^^^^^^^^^^^

**nbondc** = [integer] **ibondc, ilbondc, jbondc, jlbondc** = [integer
integer integer integer]

.. _521-spherical-harmonics-expansion:

5.21 Spherical Harmonics expansion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _5211-onoff:

5.21.1 (On/Off)
^^^^^^^^^^^^^^^

**lsphexp** = [boolean]

Default: **lsphexp** = .FALSE.

Determines whether or not to use spherical harmonics expansion. Note
this is only read for molecular dynamics simulations.

.. _5212-skip-step:

5.21.2 (skip step)
^^^^^^^^^^^^^^^^^^

**nskpsphexp** = [integer]

Default: **nskpsphexp** = 5

.. _5213-default-radius:

5.21.3 (default radius)
^^^^^^^^^^^^^^^^^^^^^^^

**radsphexp** = [real]

Default: **radsphexp** = 0.0

Gives the default radius in a.u.. If **radsphexp** = 0, slater_r are
used.

.. _522-eda:

5.22 EDA
~~~~~~~~

.. _5221-onoff:

5.22.1 (On/Off)
^^^^^^^^^^^^^^^

**leda** = [boolean]

Default: **leda** = .FALSE.

Determines whether or not to perform energy density analysis (EDA). Note
this can only be used in molecular dynamics.

.. _5222-skip-step:

5.22.2 (skip step)
^^^^^^^^^^^^^^^^^^

**nskpeda** = [intger]

Default: **nskpeda** = 5

Will peform EDA every **nskpeda** steps.

.. _5223-radius-for-grid-eda:

5.22.3 (radius for grid EDA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**radgeda** = [real]

Default: **radgeda** = 10.0

Gives the radius in a.u. for EDA calculation.

.. _523-wannier-function:

5.23 Wannier function
~~~~~~~~~~~~~~~~~~~~~

.. _5231-onoff:

5.23.1 (On/Off)
^^^^^^^^^^^^^^^

**lwannier** = [boolean]

Default: **lwannier** = .FALSE.

Determines whether or not to compute maximully localized Wannier
functions.

.. _5232-band-index:

5.23.2 (band index)
^^^^^^^^^^^^^^^^^^^

**iwbnd1, iwbnd2** = [integer] [integer]

Default: **iwbnd1, iwbnd2** = 0 0

Gives the range of band indices to be transformed by Unitary matrix. If
**iwbnd1, iwbnd2** = 0 0, then all occupied bands will be transformed.

.. _5233-iteration:

5.23.3 (iteration)
^^^^^^^^^^^^^^^^^^

**iterwan** = [integer]

Default: **iterwan** = 200

Gives the maximum number of iterations.

.. _5234-tolerance:

5.23.4 (tolerance)
^^^^^^^^^^^^^^^^^^

**tolwan** = [real]

Default: **tolwan** = 1.d-06

Gives the tolerance???

.. _5235-dump-orbitals:

5.23.5 (dump orbitals)
^^^^^^^^^^^^^^^^^^^^^^

**iwstt1, iwstt2** = [integer] [integer]

Default: **iwstt1, iwstt2** = 1 0

Gives the orbital indices to dump data for. If **iwstt1, iwstt2** = 0 0,
then all orbitals are dumped. If **iwstt1, iwstt2** = 1 0, then no
orbitals are dumped.

.. _5236-orbitals-around-atoms:

5.23.6 (orbitals around atoms)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**natwan** = [integer] **natwno, nwanaa** = [integer] [integer]

Default: **natwan** = 0

**natwan** gives the number of focused atoms, followed by that many
lines of **natwno, nwanaa** pairs, where **natwno** is the atom number
and **nwanaa** is the wavefunction number.

.. _5237-skip-step:

5.23.7 (skip step)
^^^^^^^^^^^^^^^^^^

**nskpwann** = [integer]

Default: **nskpwann** = 5

Gives the number of steps to skip.

.. _5238-output-unitary-matrix:

5.23.8 (output unitary matrix)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**loutuni** = [boolean]

Default: **loutuni** = .FALSE.

Determines whether or not to output the unitary matrix.

.. _524-conductivity:

5.24 Conductivity
~~~~~~~~~~~~~~~~~

.. _5241-onoff:

5.24.1 (On/Off)
^^^^^^^^^^^^^^^

**lconduct** = [boolean]

Default: **lconduct** = .FALSE.

.TRUE. : Execute the conductivity calculation .FASLE. : Do not execute
the conductivity calculation

Determines whether to run the conductivity calculation. At least one of
**ldcconduct** or **lacconduct** is set to .TRUE.

.. _5242-skip-step:

5.24.2 (skip step)
^^^^^^^^^^^^^^^^^^

**nskpconduct** = [integer]

Default: **nskpconduct** = 5

**nskpconduct** is skip step for the conductivity calculation.

.. _5243-unit-of-energy:

5.24.3 (unit of energy)
^^^^^^^^^^^^^^^^^^^^^^^

(ry) or (hr) or (ev)

Default: (ry)

.. _5244-width-of-energy:

5.24.4 (width of energy)
^^^^^^^^^^^^^^^^^^^^^^^^

**wgconduct** = [real]

Default: **wgconduct** = 0.1

The gaussian broadening is used, alternative to the Delta function.

.. _5245-temperature:

5.24.5 (temperature)
^^^^^^^^^^^^^^^^^^^^

**tempconduct** = [real]

Default: **tempconduct** = 300.0

Temperature in [K]

.. _5246-dc-current-density:

5.24.6 (DC current density)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lconduct** = [boolean]

Default: **ldcconduct** = .FALSE.

.TRUE. : .FASLE. :

???

.. _5247-electric-field:

5.24.7 (electric field)
^^^^^^^^^^^^^^^^^^^^^^^

**efconduct** = (/ [real] [real] [real] /)

Default: **efconduct(1:3)** = (/ 0.0 0.0 0.0/)

**efconduct(1:3)** is the constant electric-field vector. The unit of
**efconduct(1:3)** is atomic unit.

.. _5248-current-states:

5.24.8 (current states)
^^^^^^^^^^^^^^^^^^^^^^^

(**icband1**, **icband2**) = [integer] [integer]

Default: (**icband1**, **icband2**) = (/ 0, 0 /)

Range of band indices for total current. If both **icband1** and
**icband2** are set to zero, the range is determined automatically.

.. _5249-hole-states:

5.24.9 (hole states)
^^^^^^^^^^^^^^^^^^^^

(**ichband1**, **ichband**) = [integer] [integer]

Default: (**ichband1**, **ichband2**) = (/ 0, 0 /)

Range of band indices for hole current. If both **ichband1** and
**ichband2** are set to zero, the range is determined automatically.

.. _52410-electron-states:

5.24.10 (electron states)
^^^^^^^^^^^^^^^^^^^^^^^^^

(**iceband1**, **iceband**) = [integer] [integer]

Default: (**iceband1**, **iceband2**) = (/ 0, 0 /)

Tange of band indices for electron current. If both **iceband1** and
**iceband2** are set to zero, the range is determined automatically.
Note:

icband1 <= ichband1 <= ichband2 <= iceband1 <= iceband2 <= icband2

.. _52411-output-momentum-matrix:

5.24.11 (output momentum matrix)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lacconduct** = [boolean]

Default: **lacconduct** = .FALSE.

.TRUE. : .FASLE. :

??? to caluculate optical conductivity

.. _52412-unit-of-energy:

5.24.12 (unit of energy)
^^^^^^^^^^^^^^^^^^^^^^^^

(ry) or (hr) or (ev)

Default: (ry)

.. _52413-frequency-range:

5.24.13 (frequency range)
^^^^^^^^^^^^^^^^^^^^^^^^^

**freqacmx** = [real]

Default: **freqacmx** = 1.0

??? occupied: e_HOMO - **freqacmx** ~ e_HOMO unoccupied: e_LUMO ~ e_LUMO
+ **freqacmx**

.. _52414-matrix-states:

5.24.14 (matrix states)
^^^^^^^^^^^^^^^^^^^^^^^

(**iceband1**, **iceband**) = [integer] [integer]

Default: (**immtband1**, **immtband2**) = (/ 0, 0 /)

??? range of band indices for mometum matrix Both of **immtband1** and
**immtband2** are zero, the range is determined automatically.

.. _525-stress-calculation--------only-for-bulk-calculations:

5.25 stress calculation ------ only for bulk calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _5251-onoff:

5.25.1 (On/Off)
^^^^^^^^^^^^^^^

**lstress** = [boolean]

Default: **lstress** = .TRUE.

.TRUE. : Estimate the stress .FASLE. : Do not estimate the stress

Determines whether to estimate the stress.

.. _5252-skip-step:

5.25.2 (skip step)
^^^^^^^^^^^^^^^^^^

**nskip_stress** = [integer]

Default: **nskip_stress** = 5

**nskip_stress** is skip step to estimate the stress.

.. _526-dump-charge-density:

5.26 dump charge density
~~~~~~~~~~~~~~~~~~~~~~~~

.. _5261-onoff:

5.26.1 (On/Off)
^^^^^^^^^^^^^^^

**ldpchg** = [boolean]

Default: **ldpchg** = .FALSE.

.TRUE. : Output charge density .FASLE. : Do not output

Determines whether to output the charge density.

.. _5262-skip-step:

5.26.2 (skip step)
^^^^^^^^^^^^^^^^^^

**nskip_dpchg** = [integer]

Default: **nskip_dpchg** = 5

**nskip_dpchg** is skip step to output the charge density.

.. _5263-output-area:

5.26.3 (output area)
^^^^^^^^^^^^^^^^^^^^

**x_min** **x_max** = [real] [real] **y_min** **y_max** = [real] [real]
**z_min** **z_max** = [real] [real]

Default: **x_min** **x_max** = 1.0 0.0 **y_min** **y_max** = 1.0 0.0
**z_min** **z_max** = 1.0 0.0

These parameter specify the area to output the charge density. They are
reduced coordinates. If **[x|y|z]_min** > **[x|y|z]_max**, output are is
the whole space.

.. _527-dump-wavefunctions:

5.27 dump wavefunctions
~~~~~~~~~~~~~~~~~~~~~~~

.. _5271-onoff:

5.27.1 (On/Off)
^^^^^^^^^^^^^^^

**ldpwav** = [boolean]

Default: **ldpwav** = .FALSE.

.TRUE. : Output wavefunction .FASLE. : Do not output

Determines whether to output the wavefunction.

.. _5272-bands:

5.27.2 (bands)
^^^^^^^^^^^^^^

(**ibstt1**, **ibstt2**) = [integer] [integer]

Default: (**ibstt1**, **ibstt2**) = (/ 0, 0 /)

The wavefunction will be output from **ibstt1** to **ibstt2**.
**ibstt1** and **ibstt2** are the index of bands. If both **ibstt1** and
**ibstt2** are zero, all bands will be output.

.. _5273-skip-step:

5.27.3 (skip step)
^^^^^^^^^^^^^^^^^^

**nskip_dpwav** = [integer]

Default: **nskip_dpwav** = 5

**nskip_dpwav** is skip step to output the wavefunction.

.. _5274-output-area:

5.27.4 (output area)
^^^^^^^^^^^^^^^^^^^^

**x_min** **x_max** = [real] [real] **y_min** **y_max** = [real] [real]
**z_min** **z_max** = [real] [real]

Default: **x_min** **x_max** = 1.0 0.0 **y_min** **y_max** = 1.0 0.0
**z_min** **z_max** = 1.0 0.0

These parameter specify the area to output the wavefunction. They are
reduced coordinates. If **[x|y|z]_min** > **[x|y|z]_max**, output are is
the whole space.

.. _528-dump-local-potential:

5.28 dump local potential
~~~~~~~~~~~~~~~~~~~~~~~~~

.. _5281-onoff:

5.28.1 (On/Off)
^^^^^^^^^^^^^^^

**ldppot** = [boolean]

Default: **ldppot** = .FALSE.

.TRUE. : Output local potential .FASLE. : Do not output

Determines whether to output the local potential.

.. _5282-average-plane:

5.28.2 (average plane)
^^^^^^^^^^^^^^^^^^^^^^

**nav_dppot** = [1|2|3]

Default: **nav_dppot** = 1

============= =====
**nav_dppot** Plane
============= =====
**1**         xy
**2**         yz
**3**         zx
============= =====

The plane-averaged local potential will be output.

.. _5283-skip-step:

5.28.3 (skip step)
^^^^^^^^^^^^^^^^^^

**nskip_dppot** = [integer]

Default: **nskip_dppot** = 5

**nskip_dppot** is skip step to output the local potential.

.. _5284-output-area:

5.28.4 (output area)
^^^^^^^^^^^^^^^^^^^^

**x_min** **x_max** = [real] [real] **y_min** **y_max** = [real] [real]
**z_min** **z_max** = [real] [real]

Default: **x_min** **x_max** = 1.0 0.0 **y_min** **y_max** = 1.0 0.0
**z_min** **z_max** = 1.0 0.0

These parameter specify the area to output the local potential. They are
reduced coordinates. If **[x|y|z]_min** > **[x|y|z]_max**, output are is
the whole space.

.. _529-supercell:

5.29 supercell
~~~~~~~~~~~~~~

.. _5291-unit-of-length:

5.29.1 (unit of length)
^^^^^^^^^^^^^^^^^^^^^^^

(bohr) or (ang)

Default: (bohr)

Determines the unit of supercell (bohr or angstrom). And the supercell
is specified by either cell vector or 'lengths and angles'. You need
comment out the subsection you do not specify. Note that for cluster
calculations, the supercell has to be ORTHORHOMBIC.

.. _5292-cell-vector:

5.29.2 (cell vector)
^^^^^^^^^^^^^^^^^^^^

L1 = [real] [real] [real] L2 = [real] [real] [real] L3 = [real] [real]
[real]

Example:

::

    10.68   0.0    0.0        : super cell vector L1
     0.0   10.68   0.0        : super cell vector L2
     0.0    0.0   10.68       : super cell vector L3

.. _5293-lengths--angles:

5.29.3 (lengths & angles)
^^^^^^^^^^^^^^^^^^^^^^^^^

Example:

::

    10.68,  10.68,  10.68  :  lengths of cell vectors
    90.00,  90.00,  90.00  :  angles between cell vec. in [deg.]

.. _530-vacuum:

5.30 vacuum
~~~~~~~~~~~

This section is mainly used for the double grid method.

.. _5301-unit-of-length:

5.30.1 (unit of length)
^^^^^^^^^^^^^^^^^^^^^^^

(bohr) or (ang)

Defaubr> (bohr)

Determines the unit of vacuum (bohr or angstrom).

.. _5302-onoff:

5.30.2 (On/Off)
^^^^^^^^^^^^^^^

**lvacuum(1)**, **vacuum(1)** [boolean] [real] **lvacuum(2)**,
**vacuum(2)** [boolean] [real] **lvacuum(3)**, **vacuum(3)** [boolean]
[real]

Default: **lvacuum(1:3)** = .FALSE. **vacuum(1:3)** = 0.0

If **lvacuum([1|2|3])** is .TRUE., the vaccum is **vacuum([1|2|3])**\ <
[x|y|z] < [L1|L2|L3].

.. _5303-parameter-in-error-function:

5.30.3 (parameter in error function)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**alpha_ldouble_grid_recip** = [real]

Default: **alpha_ldouble_grid_recip** = 0.0

If **alpha_ldouble_grid_recip** is zero, **alpha_ldouble_grid_recip** is
set automatically. **alpha_ldouble_grid_recip** should be about 5.0/(the
length of cell).

.. _531-spherical-region:

5.31 spherical region
~~~~~~~~~~~~~~~~~~~~~

.. _5311-onoff:

5.31.1 (On/Off)
^^^^^^^^^^^^^^^

**lsphere** = [boolean]

Default: **lsphere** = .FALSE.

.TRUE. : Use spherical region in cluster calculations .FASLE. : Do not
use

.. _532-planewaves:

5.32 planewaves
~~~~~~~~~~~~~~~

.. _5321-unit-of-cutoff-energy:

5.32.1 (unit of cutoff energy)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

(ry) or (hr) or (ev)

Default: (ry)

Determines the unit of cutoff energy.

.. _5322-for-wavefuctions:

5.32.2 (for wavefuctions)
^^^^^^^^^^^^^^^^^^^^^^^^^

**ecut** = [real]

**ecut** is the cutoff energy to expand wavefunctions.

.. _5323-for-electron-density:

5.32.3 (for electron density)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**ecutdens** = [real]

Default: **ecutdens** = 0.0

**ecutdens** is the cutoff energy to expand electron density.
**ecutdens** must be greater than **ecutsoft**. if .not.lvand,
**ecutdens** must be equal to **ecutsoft**.

.. _5324-for-soft-part-of-density:

5.32.4 (for soft part of density)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**ecutsoft** = [real]

Default: **ecutorth** = 0.0

**ecutsoft** is the cutoff energy to expand soft part of electron
density. It must be greater than ecut, and smaller than ecut*4.

.. _5325-for-orthogonalization:

5.32.5 (for orthogonalization)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**ecutorth** = [real]

Default: **ecutorth** = 0.0

**ecutorth** is the cutoff energy for orthogonalization in CG iteration.

.. _5326-for-charge-density-mixing:

5.32.6 (for charge density mixing)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**ecutc** = [real]

Default: **ecutc** = 0.0

**ecutc** is the cutoff energy for charge density mixing.

.. _533-double-grid-method:

5.33 double-grid method
~~~~~~~~~~~~~~~~~~~~~~~

[add ref]

.. _5331-order-of-cardinal-b-spline:

5.33.1 (order of cardinal B spline)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**iosp** = [integer]

Default: *iosp* = 3. This number must be an odd integer

.. _5332-unit-of-cutoff-energy:

5.33.2 (unit of cutoff energy)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

(ry) or (hr) or (ev)

Default: (ry)

.. _5333-for-non-periodic-direction:

5.33.3 (for non-periodic direction)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**ecutlong** = [real]

Default: **ecutlong** = 0.0

Cutoff energies for the long-range Coulomb interaction

.. _5333-for-periodic-direction:

5.33.3 (for periodic direction)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**ecutlong_small** = [real]

Default: **ecutlong_small** = 0.0

Cutoff energies for the long-range Coulomb interaction in the periodic
direction (for wire and surface geometries)

.. _534-electronic-bands:

5.34 electronic bands
~~~~~~~~~~~~~~~~~~~~~

.. _5341-occupied-bands:

5.34.1 (occupied bands)
^^^^^^^^^^^^^^^^^^^^^^^

**noband** = [integer]

Default: **noband** = 0

The number of occupied bands should at least be sufficient to
accommodate all the valence electrons in the simulation.

.. _5342-empty-bands:

5.34.2 (empty bands)
^^^^^^^^^^^^^^^^^^^^

**neband** = [integer]

Default: **neband** = 0

A good rule of thumb is to set neband equal to 10% of occupied bands.

.. _5343-broadening:

5.34.3 (broadening)
^^^^^^^^^^^^^^^^^^^

**lfermi** = [integer] flag **tfermi** = [real]

Default: **lfermi** = 3 **tfermi** = 2000.0

lfermi = 1 for Non-metallic systems, lfermi = 2 for Fermi smearing,
lfermi = 3 for Gaussian smearing and lfermi > 3 for Methfessel-Paxton
smearing scheme of order lfermi-3

.. _5344-charge-number-of-ion:

5.34.4 (charge number of ion)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**ncion** = [integer]

Default: **ncion** = 0. For calculation of electronically non-neutral
systems. The number of electrons is given by ``nel-ncion``.
Charged-periodic calculations are perfomed with the uniform background.

.. _535-symmetry-operations:

5.35 Symmetry operations
~~~~~~~~~~~~~~~~~~~~~~~~

.. _5351-onoff:

5.35.1 (On/Off)
^^^^^^^^^^^^^^^

**lsymop** = [boolean]

Default: **lsymop** = .FALSE.

The flag **lsymop** defines if symmetry operations are checked.

.. _5352-unit-cell:

5.35.2 (unit cell)
^^^^^^^^^^^^^^^^^^

**nsymm** = [integer][integer][integer]

Default: **nsymm** = 1 1 1

**nsymm** defines the number of unit cells in the x, y and z directions.

The transformation matrix is specified through either (transformation
matrix) or (inverse transformation matrix) as follows.
``trsfrm(j,i) <-- trsfrm(j,i)/nsymm(i)`` or
``trsfri(i,j) <-- trsfri(i,j)*nsymm(i)``

.. _5353-transformation-matrix:

5.35.3 (transformation matrix)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**trsfrm(1,i)** = [real][real][real] **trsfrm(2,i)** =
[real][real][real] **trsfrm(3,i)** = [real][real][real]

Default: **trsfrm** = No defaults

.. _5353-inverse-transformation-matrix:

5.35.3 (inverse transformation matrix)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**trsfri(1,i)** = [real][real][real] **trsfri(2,i)** =
[real][real][real] **trsfri(3,i)** = [real][real][real]

Default: **trsfri** = No defaults

.. _5354-origin-of-symmetry-op:

5.35.4 (origin of symmetry op.)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**symm** = [real][real][real]

Default: **symm** = 0.0 0.0 0.0

.. _5355-symmetry-operations:

5.35.5 (symmetry operations)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**nmsyop** = [integer] ``???``

Default: **nmsyop** = 1

.. _5356-denominator-of-trans-vector:

5.35.6 (denominator of trans. vector)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**denomt** = [real]

Default: **denomt** = No default

.. _5357-inversion-operations:

5.35.7 (inversion operations)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**linvsn** = [boolean]

Default: **linvsn** = .FALSE.

.. _5358-radius-to-identify-atoms:

5.35.8 (radius to identify atoms)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**radsym** = [real]

Default: **radsym** = 1.0

Radius to identify symmetrically equivalent atoms

.. _5359-for-charge-density:

5.35.9 (for charge density)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lsymcd** = [boolean]

Default: **lsymcd** = .FALSE.

.. _53510-for-spin-density:

5.35.10 (for spin density)
^^^^^^^^^^^^^^^^^^^^^^^^^^

**lsymspin** = [boolean]

Default: **lsymspin** = No default

.. _536-atoms:

5.36 Atoms
~~~~~~~~~~

.. _5361-species:

5.36.1 (species)
^^^^^^^^^^^^^^^^

**ntype** = [integer]

Default: **ntype** = 0

Number of atom types. The subsequent sections must be duplicated
``ntype`` times.

.. _5362-atomic-number:

5.36.2 (atomic number)
^^^^^^^^^^^^^^^^^^^^^^

**zatom** = [integer]

Default: **zatom** = No default

Atomic number of the current atom type

.. _5363-pseudopotential:

5.36.3 (pseudopotential)
^^^^^^^^^^^^^^^^^^^^^^^^

kbpp or uspp or vand or local

Default: No default

Type of pseudopotential defined for the current atom type.

.. _5364-nonlocal-potential:

5.36.4 (nonlocal potential)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lking** = [boolean] **rking** = [real] **gkgmax** = [real] **gkgexct**
= [real]

Default: **lking** = No default **rking** = No default **gkgmax** = No
default **gkgexct** = No default

These flags specify if nonlocal potntial is computed in real space or
reciprocal space. We follow the formalism by [ R. D. King-smith et al.
]. lking=.false. allows calculation to proceed in reciporocal space. In
such case, next three parameter are read but not use. If lking=.true.,
caculation will proceed in real space. In such case, rking, gkgmax and
gkgexct are used. rking corresponds to real space cutoff radius. gkgmax
and gkgexct correspond to accuracy of non-local potentials. For more
detail, please follow paper by R. D. King-smith et al. Accuracy of real
space must be tested if input parameters are not provided.

Recommendation: **lking** = false for small system and true for large
system **rking** = 1.5- 2 **gkgmax** = 1.25 **gkgexct** = 0.98

.. _5365-local-potential:

5.36.5 (local potential)
^^^^^^^^^^^^^^^^^^^^^^^^

**llking** = [boolean] **rlking** = [real] **glkgmax** = [real]
**glkgexct** = [real]

Default: **llking** = No default **rlking** = No default **glkgmax** =
No default **glkgexct** = No default

These flags specify if local potntial is computed in real space or
reciprocal space. It is recommended to use reciprocal space reciprocal
space for the calculation. In reciprocal space rlking, glkgmax and
glkgexct are not used but values are read.

Recommendation: **llking** = false

.. _5366-partial-core-correction:

5.36.6 (partial core correction)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lpcc** = [boolean] **rpcc** = [real] **lpking** = [boolean] **rpking**
= [real] **gpkgmax** = [real] **gpkgexct** = [real]

Default: **lpcc** = No default **rpcc** = No default **lpking** = No
default **rpking** = No default **gpkgmax** = No default **gpkgexct** =
No default

Psuedopotential formalism treat core and valence electron seperately.
Partial core correction includes the effect of core charge as a
perturbation which increase the transferrability of the code. It is
recommended to use partial core correction. Please see paper by Louie et
al.

.. _5367-compensation-charge-cutoff:

5.36.7 (compensation charge cutoff)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**frac_rcomp_it** = [real]

Default: **frac_rcomp_it** = 1.3

Compensation charge must be used in PAW method to match all electron
wave function. A recommended value is between 1.2 - 1.3. For more detail
please see paper by Kresse et al.

.. _5368-the-number-of-atoms:

5.36.8 (the number of atoms)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**nhk** = [integer]

Default: **nhk** = 0

Number of atoms of the current type. If ``nhk`` is zero, then atomic
configuration is read from the positions file

.. _5369-unit-of-length:

5.36.9 (unit of length)
^^^^^^^^^^^^^^^^^^^^^^^

(bohr) or (ang)

Default: No default

.. _53610-position-file:

5.36.10 (position file)
^^^^^^^^^^^^^^^^^^^^^^^

**fname** = [String] **scaledflag**\ \* = [integer] **keyword**\ \* =
[integer]

Default: **fname** = No default **scaledflag**\ \* = No default
**keyword**\ \* = No default

Filename to read atomic positions from. Filenames are given with
reference to the current working directory. **scaledflag** determines if
the aomtic positions in the coordinate file are in real or scaled units.
Atoms of this type are denoted by the integer **keyword** in the
corrdinates file.

The structure of 'CONFIGURATION_FILE' is as follows.

::

      ntotal
      n_1      x_1      y_1      z_1
      n_2      x_2      y_2      z_2
      ...      ...      ...      ...
      n_ntotal x_ntotal y_ntotal z_ntotal

The coordinates for (n_i == keyword) are selected.

.. _53611-positions:

5.36.11 (positions)
^^^^^^^^^^^^^^^^^^^

**nunitcells** = [integer][integer][integer] **scaledflag**\ \* =
[integer] ``nhk`` lines of **positions** = [real][real][real]

Default: **nunitcells** = No defaults **scaledflag**\ \* = No defaults

**nunitcells** gives the number of unit cells along the x, y and z
directions. The flag **scaledflag** defines if the atomic positions are
defined in real or scaled coordinates (**scaledflag** = 1 is scaled
coordinates, **scaledflag** = 2 is real coordinates).

If '1 (scaled)' is selected, the scaled coordinates ( 0.0 < x,y,z < 1.0
) have to be given above. If '2 (real)' is selected, the real
coordinates ( in a.u. ) have to be given like below.

::

    2                            : 1:scaled, 2:real coordinates
     0.00d0  0.00d0  0.00d0      : If '2 (real)' is selected,
    10.34d0 10.34d0 10.34d0      :  the real coordinates ( in a.u. )
    10.34d0  0.00d0 10.34d0      :  have to be given.
   0.00d0 10.34d0 10.34d0

.. _53612-fix-positions:

5.36.12 (fix positions)
^^^^^^^^^^^^^^^^^^^^^^^

**lfixion** = [boolean]

Default: **lfixion** = .FALSE.

If .TRUE., atoms of current type are held fixed by setting atomic
velocities to zero.

.. _53613-fix-displacement:

5.36.13 (fix displacement)
^^^^^^^^^^^^^^^^^^^^^^^^^^

**ldispion** = [boolean] **dvector** = [real] [real] [real]

Default: **ldispion** = .FALSE. **dvector** = No default

If .TRUE., positions of atoms of current type are updated by adding
**dvector** to the current position. ``???``

.. _53614-fictitious-mass:

5.36.14 (fictitious mass)
^^^^^^^^^^^^^^^^^^^^^^^^^

**ficmass** = [real]

Default: **ficmass** = 1.0

If set, the mass of atoms of this type are scaled by **ficmass**

.. _53615-velocity-file:

5.36.15 (velocity file)
^^^^^^^^^^^^^^^^^^^^^^^

**fname** = [String] **keyword**\ \* = [integer]

Default: **fname** = No default **keyword**\ \* = No default

**fname** is the filename to read atomic velocities from. Filename is
given with reference to the current working directory. Atoms of this
type are denoted by the integer **keyword** in the velocities file.

The structure of 'VELOCITY_FILE' is as follows.

::

      ntotal
      vmax                              : maximum velocity in (a.u.)
      n_1      vx_1      vy_1      vz_1 : an integer & scaled velocity
      n_2      vx_2      vy_2      vz_2
      ...      ...       ...       ...
      n_ntotal vx_ntotal vy_ntotal vz_ntotal

The velocities for (n_i == keyword) are selected.

.. _53616-velocities:

5.36.16 (velocities)
^^^^^^^^^^^^^^^^^^^^

**nunitcells** = [integer][integer][integer] **scalingfactor**\ \* =
[real] ``nhk`` lines of scaled **velocities** = [real][real][real] .
Ensure that -1.0 < **velocities** < 1.0

Default: **nunitcells** = No defaults **scalingfactor**\ \* = No
defaults

**nunitcells** gives the number of unit cells along the x, y and z
directions. Atomic velocities are initialized by multiplying the vector
**velocities** by the **scalingfactor**.

.. _53617-spherical-harmonics-expansion:

5.36.17 (Spherical Harmonics expansion)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**rsphexp** = [real]

Default: **rsphexp** = No default

Radius for *Ylm* expansion

.. _53618-eda:

5.36.18 (EDA)
^^^^^^^^^^^^^

**radeda** = [real]

Default: **radeda** = No default

Radius for EDA

.. _53619-atomic-charge:

5.36.19 (atomic charge)
^^^^^^^^^^^^^^^^^^^^^^^

**rintchg** = [real]

Default: **rintchg** = No default

Radius of integration to calculate atomic charge

.. _53620-dftu:

5.36.20 (DFT+U)
^^^^^^^^^^^^^^^

**lplusU_at** = [boolean] **plusU_l** = [integer] **plusU_U** = [real]
**plusU_J** = [real]

Default: **lplusU_at** = .FALSE. **plusU_l** = 2 **plusU_U** = 0.0
**plusU_J** = 0.0

The flag **lplusU_at** determines if the DFT+U method is used.
**plusU_l** defines the angular momentum quantum number of the orbital
of interest. **plusU_U** and **plusU_J** define the Hubbard parameter U
and the screened exchange energy in electron-volts.

.. _53621-dftc:

5.36.21 (DFT+C)
^^^^^^^^^^^^^^^

**lplusC_at** = [boolean] **plusC_r_e_s** = [real] [real]
**plusC_r_e_p** = [real] [real] **plusC_r_e_d** = [real] [real]
**plusC_r_e_f** = [real] [real]

Default: **lplusC_at** = .FALSE. **plusC_r_e_s** = 0.0 0.0
**plusC_r_e_p** = 0.0 0.0 **plusC_r_e_d** = 0.0 0.0 **plusC_r_e_f** =
0.0 0.0

The flag **lplusC_at** determines if the DFT+C method is used.
**plusC_r_e** defines the ``plusC_r`` and ``plusC_e`` vales for the s,
p, d and f subshells.

The form of correction

``modified V_nl(r) = V_nl(r) + plusC_e * sin(a*r)/(a*r), a = pi/plusC_r``

where ``plusC_r`` is cutoff radius in [a.u.] and ``plusC_e`` is energy
shift in [Ryd.]

.. _537-constraint-conditions:

5.37 Constraint conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _5371-unit-of-length:

5.37.1 (unit of length)
^^^^^^^^^^^^^^^^^^^^^^^

(bohr) or (ang)

Default: (ang)

.. _5371-bond-length-constraint:

5.37.1 (bond length constraint)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**ncbonds** = [integer]

Followed by **ncbonds** lines of constraints.

**ncbatm1** = [integer] **ncbatm2** = [integer] **cblength** = [real]

For each line, the bond length between **ncbatm1** and **ncbatm2** is
held fixed at **cblength** Angstrom.

Default: **ncbonds** = 0 **ncbatm1** = No default **ncbatm2** = No
default **cblength** = No default

.. _538-virtual-molecular-dynamics:

5.38 Virtual Molecular Dynamics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _539-cholesky-decomposition:

5.39 Cholesky decomposition
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _5391-the-number-of-nodes:

5.39.1 (the number of nodes)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**node_c** = [integer]

Default: **node_c** = 4

.. _540-eigenvalue-problem:

5.40 Eigenvalue problem
~~~~~~~~~~~~~~~~~~~~~~~

.. _5401-the-number-of-nodes:

5.40.1 (the number of nodes)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**node_r** = [integer]

Default: **node_r** = 4

.. _541-work-array:

5.41 Work Array
~~~~~~~~~~~~~~~

.. _5411-size-of-array:

5.41.1 (size of array)
^^^^^^^^^^^^^^^^^^^^^^

**dblock** = [real]

Default: **dblock** = 1.d+07

.. _542-table-dimension:

5.42 Table Dimension
~~~~~~~~~~~~~~~~~~~~

.. _5421-local-potential:

5.42.1 (local potential)
^^^^^^^^^^^^^^^^^^^^^^^^

**mx1loc** = [integer]

Default: **mx1loc** = 2**13

.. _5422-others:

5.42.2 (others)
^^^^^^^^^^^^^^^

**mx1** = [integer]

Default: **mx1** = 2**13

.. raw:: html

   <br>

--------------
