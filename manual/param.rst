Input Parameters: IN.PARAM 
==========================
.. Note:: The name of the main input file (e.g. 'IN.PARAM') must be specified in 'control/filename'.

The main input file is divided into sections corresponding to different
controls. Sections start with ``*SECTION_NAME`` and end with ``*end``. To
disable any section, simply prepend section title with a ``#``.

To disable section SECTION_NAME_2

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

.. Note:: Any unnecessary sections or subsections may be removed from input.file, depending on your simulation needs (recommended only for advanced users).


5.1 Parallel Section: \*parallel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


5.1.1 (QM-Nodes)
^^^^^^^^^^^^^^^^

**npx, npy, npz** = [1|even integer] [1|even integer] [1|even integer]

::

   Default **npx** = 1 **npy** = 1 **npz** = 1

QM-Nodes is a space delimited set of three integers **npx**, **npy**,
**npz**. **QXMD** uses a hybrid spatial and band decompisition for
parallel computing, and this set of three integers multiplied together
defines the number of MPI ranks. 


5.1.2 (k-points)
^^^^^^^^^^^^^^^^

**npk** = [integer]

::
 
   Default **npk** = 1

Defines the number of MPI ranks to parallelize k-point sampling. See
section 5.8 for k-point sampling.



5.1.3 (MD-nodes)
^^^^^^^^^^^^^^^^

**md_npx, md_npy, md_npz** = [1|even integer] [1|even integer] [1|even
integer]

::

   Default: **md_npx =1 , md_npy =1, md_npz** = 1

MD-Nodes is a space delimited set of three integers: **md_npx, md_npy,
md_npz**. 

These variables are used only for spatial decomposition inclassical MD and divide-and-conqure QMD simulations. For linear scaling DFT please see reference. :cite:`RN19` :cite:`RN20` :cite:`RN22` :cite:`RN23` :cite:`RN24` :cite:`RN25`


5.2 Start(on/off) Section
~~~~~~~~~~~~~~~~~~~~~~~~~


5.2.1 (On/Off)
^^^^^^^^^^^^^^

**lstart** = [boolean]

::

   Default **lstart** =  ``.FALSE.`` 


``.FALSE.``  : Start self-consistent field (SCF) iteration using a random initial wavefunction 

``.TRUE.`` : Continue SCF iteration from a previous run's output wavefunction. 

.. Note:: QM_* output files must be present in the **data/** directory from the previous calculation.

Determines whether you would like to restart a successfully completed
simulation. Note that if you are restarting a NAQMD simulation, you will
also want to set **ltddft_start** to .TRUE. in the \*TDDFT-MD section to
properly restart a successfully completed simulation.


5.3 TDDFT-MD
~~~~~~~~~~~~


5.3.1 (On/Off)
^^^^^^^^^^^^^^

**ltddft** = [boolean]

::

   Default **ltddft** =  ``.FALSE.`` 


``.FALSE.``  : Execute Adiabtic QMD based on density functional theory (DFT)

``.TRUE.`` : Execute non-adiabtic QMD (NAQMD) based on time-dependent density functional theory (TDDFT). :cite:`RN26`

Determines whether to run QMD simulation under adiabatic or
non-adiabatic methods. Adiabatic methods simulate thermodynamic system
in ground state equilibrium, while non-adiabatic methods simulate
electronic excitations.


5.3.2 (FSSH)
^^^^^^^^^^^^

**ltddft_fssh** = [boolean]

::
	
   Default **ltddft_fssh** = .TRUE.


``.TRUE.`` : Perform NAQMD based on Fewest Switches Surface Hopping (FSSH) method  

``.FALSE.`` : :cite:`RN27`. Perform NAQMD based on based on

.. note:: Ehrenfest dynamics (not yet implemented).

Determines the implementation method for electron state dynamics in
NAQMD. Fewest Switches surface hopping method is proposed by J. Tully :cite:`RN27` molecular dynamics simulation of the processes
including electronic transition. Next few flags ask for the specification of FSSH method


5.3.3 (FSSH-switch)
^^^^^^^^^^^^^^^^^^^

**lfssh_switch** = [boolean]

::

   Default **lfssh_switch** = .TRUE.

``.TRUE.`` : Allow electrons to move between excited states.  

``.FALSE.``: Keep electronic occupations fixed.

Determines whether electronic occupations can change throughout the
NAQMD simulation.


5.3.4 (FSSH-ground-state-SCF)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lfssh_gsscf** = [boolean]

::

   Default **lfssh_gsscf** =  ``.FALSE.`` 

``.TRUE.`` : SCF iterations performed based on the ground state  

``.FALSE.``:SCF with the excited state

.. note:: This parameter should always be set to  ``.TRUE.``  to obtain convergence.

5.3.5 (FSSH-charge mixing)
^^^^^^^^^^^^^^^^^^^^^^^^^^

**imxchg_fssh** = 0 \| 1 \| 2 \| 3 **aslh_fssh, bslh_fssh** = [real]
[real]

::

   Default **imxchg_fssh** = 1


=========== ===============================================
imxchg_fssh Charge Mixing Method                           
=========== ===============================================
**0**       No Mixing                                      
**1**       Pulay :cite:`RN28` 0.9, 0.6(Recommended values)
**2**       Anderson                                       
**3**       Simple                                         
=========== ===============================================

**imxchg_fssh** is used to specify the method for charge mixing during
SCF iterations, while **aslh_fssh** and **bslh_fssh** are used as tuning
parameters. Note that this section is only considered when
**lfssh_gsscf** is set to .TRUE.


5.3.6 (FSSH-random-initialize)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lfssh_random** = [boolean] **rseed_fssh** = [real]

::

   Default **lfssh_random** = .FALSE.

**lfssh_random** = ``.TRUE.`` : Automatically seed random number generator.
**lfssh_random** = ``.FALSE.`` : Specify the seed for the random number
generator with value given by **rseed_fssh**.

Determines how to specify the seed for the random number generator used
by the FSSH method.


5.3.7 (Boltzmann factor for upward transition)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lfssh_boltzmn** = [boolean]

::

   Default **lfssh_boltzmn** = .TRUE.

``.TRUE.`` : Multiply electronic transition probability by the Boltzmann
factor. 

``.FALSE.`` : Leave electronic transtion probability unaltered.

Determines whether or not to multiply the electronic transition
probability by the Boltzmann factor when the electronic excitation
energy increases due to the transition. This is used in order to
approximately satisfy the detailed balance condition.


5.3.8 (velocity scaling)
^^^^^^^^^^^^^^^^^^^^^^^^

**lfssh_vscale** = [boolean] **tminimum** = [real]

::

   Default **lfssh_vscale** = .FALSE.

``.TRUE.`` : Rescale atomic velocities. 

``.FALSE.`` : Do not rescae atomicvelocities.

Determines whether or not to rescale atomic velocities upon electronic
excitation. **tminimum** gives the minimum temperature in [K] and is
used to constrain velocity scaling.


5.3.9 (time step)
^^^^^^^^^^^^^^^^^

**dttddft** = [real]

::
 
   Default **dttddft** = 0.02d0

Gives the time step in [a.u.] for numerically integrating the TDDFT
equations.


5.3.10 (parallel calculation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lfssh_parallel** = [boolean]

::
 
   Default **lfssh_parallel** = .TRUE.

.TRUE. : Solves time-dependent K-S equations in parallel .FALSE. :
Solves time-dependent K-S equations serially.

Determines whether or not to perform TDDFT calculations in parallel.


5.3.11 (restart)
^^^^^^^^^^^^^^^^

**ltddft_start** = [boolean]

::
 
   Default **ltddft_start** = .FALSE.

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

::
 
   Default **nexciton** = 0

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

::
 
   Default **lfssh_gsfrc** = .FALSE.

.TRUE. : ground state forces are used in FSSH .FALSE. : the excited
state forces are used.

Note that this variable is only read if **lfssh_gsscf** is set to .TRUE.
in the (FSSH-ground-state-SCF) subsection.


5.3.14 (NSC force)
^^^^^^^^^^^^^^^^^^

**ltddft_nscforce** = [boolean]

::
 
   Default **ltddft_nscforce** = .TRUE.

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

::
 
   Default **tdbroad** = 0.0

Determines the width of Gaussian broadening of the Fermi surface in [K].
Note: **tdbroad** = 0.0 denotes no broadening.


5.3.17 (DISH)
^^^^^^^^^^^^^

**lfssh_dish** = [boolean] **ndishpair** = [integer] **ndishi, ndishj,
decoherence_rate** = [integer] [integer] [real]

::
 
   Default **lfssh_dish** = .FALSE. **ndishpair** = 0

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
      23   25  5.117713E-03
      23   26  4.5.5.83E-03
      24   25  7.877069E-03
      24   26  7.426337E-03
      25   26  2.768402E-03

The decoherence rate for each pair of states is given by Jaeger :cite:`RN30`:

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


5.4 Approximation for Exchange
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


5.4.1 (approximation)
^^^^^^^^^^^^^^^^^^^^^

**jgga** = 1|2|3|4|5.5

::
 
   Default **jgga** = 2

======== ==================================
**jgga** Excahnge functional
======== ==================================
**1**    LDA :cite:`RN31`
**2**    GGA(PBE)\ :cite:`RN42` 
======== ==================================

**jgga** is used to specify the exchange correction functional.


5.4.2 (DFT-D)
^^^^^^^^^^^^^

**ldftd** = [boolean]

::
 
   Default **ldftd** = .FALSE.

.TRUE. : Employ an empirical vdW correction .FASLE. : Do not use an
empirical vdW correction

Determines whether to use an empirical correction for the van der Waals
interaction proposed by S. Grimme :cite:`RN37` :cite:`RN38`



5.5 SCF iterations
~~~~~~~~~~~~~~~~~~~


5.5.1 (global iterations)
^^^^^^^^^^^^^^^^^^^^^^^^^^

**iscfmx** = [integer]

::
 
   Default **iscfmx** = 100

Gives the maximum number of global iterations of SCF to be performed.


5.5.2 (tolerances)
^^^^^^^^^^^^^^^^^^^

**tolpot** = [real] **tolres** = [real]

::
 
   Default **tolpot** = 5.0d-09 **tolres** = 5.0d-09

**tolpot** gives the tolerance for change in total energy. Once the
change in total energy between subsequent SCF iterations is smaller than
the tolerance, the SCF iterations are considered to have converged.
**tolres** similarly gives the tolerance for average residual.


5.5.3 (charge mixing)
^^^^^^^^^^^^^^^^^^^^^^

**imxchg** = [0|1|2|3|4|5.5] **aslh, bslh** = [real] [real]

::
 
   Default **imxchg** = 1

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


5.5.4 (number of mixing)
^^^^^^^^^^^^^^^^^^^^^^^^^

**itratn** = [integer]

::
 
   Default **itratn** = 10

Determines how many charge densities from previous SCF iterations to use
for charge mixing. Not available for imxchg = [0|3].



5.6 Kohn-Sham equation
~~~~~~~~~~~~~~~~~~~~~~~


5.6.1 (On/Off)
^^^^^^^^^^^^^^^

**ihldam** = [1|2]

::
 
   Default **ihldam** = 1

1 : Conjugate-Gradient (CG) Method 2 : Residual minimization scheme,
direct inversion in the iterative subspace (RMM-DIIS) Method

Determines the method used in the iterative minimization of energy as a
functional of wavefunctions.


5.6.2 (tolerance)
^^^^^^^^^^^^^^^^^^

**toleig** = [real]

::
 
   Default **toleig** = 1.d-10

Gives the tolerance for the Kohn-Sham equations.


5.6.3 (threshold for w.f. direction)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**threinn** = [real]

::
 
   Default **threinn** = 0.0

Gives the threshold for direction of new wavefunction


5.6.4 (iteration)
^^^^^^^^^^^^^^^^^^

**itermx** = [integer]

::
 
   Default **itermx** = 4

Gives the maximum number of iterations.


5.6.5 (empty-band iteration)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**iteremx** = [integer]

::
 
   Default **iteremx** = 4

Gives the maximum number of iterations for empty bands.


5.6.6 (CG method)
^^^^^^^^^^^^^^^^^^

**methodcg** = [1|2]

::
 
   Default **methodcg** = 2

1:line minimization 2:BKL [ref]

Gives the method used of CG This variable is only read if **ihldam** ==
1


5.7 molecular dynamics
~~~~~~~~~~~~~~~~~~~~~~~


5.7.1 (On/Off)
^^^^^^^^^^^^^^^

**ifmd** = [0|1|2|3|4|5.7]

::
 
   Default **ifmd** = 0

====== =================
ifmd   Type of Dynamics
====== =================
**0**  Single Step
**1**  Optimization
**2**  NVE
**3**  NVT :cite:`RN45` 
**4**  NPT :cite:`RN45` 
**5**  NVT for each atom
**10** MSST :cite:`RN47` :cite:`RN48`
====== =================

Determines the type of QMD simulation to run. Add short description of
each type.


5.7.2 (time step)
^^^^^^^^^^^^^^^^^^

**dtmd, nstop** = [real] [integer]

::
 
   Default **dtmd** = 50.0 **nstop** = 10

**dtmd** gives the time step for QMD simulation in [a.u.], while
**nstop** gives the total number of time steps to simulate. Thus, the
total simulation time will equal (**nstop** \* **dtmd**)


5.7.3 (initial step number)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**nstep_ini** = [integer]

::
 
   Default **nstep_ini** = 0

Gives the initial step number. This varibale is ignored for **lstart**
== .TRUE.


5.7.4 (temperature)
^^^^^^^^^^^^^^^^^^^^

**treq** = [real]

::
 
   Default **treq** = 300.0

Gives the initial temperature in [K] for NVE, NVT, and NPT QMD
simulations.


5.7.5 (check temperature)
^^^^^^^^^^^^^^^^^^^^^^^^^^

**liscale** = [boolean] **iscnum** = [integer] **iscstp** = [integer]

::
 
   Default **liscale** = .FALSE. **iscnum** = 25 **iscstp** = 20

**liscale** = .FALSE. : Do not scale temperature. **liscale** = .TRUE. :
Scale temperature a total of **iscstp** times, with **iscnum** steps in
between each scaling.

Determines whether and how to scale temperature to keep it near the
initial given temperature.


5.7.6 (make total momentum zero)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lmomzero** = [boolean]

Deafault: **lmomzero** = .FALSE.

.TRUE. : Make the total momentum zero. .FALSE. : Do not


5.7.7 (optimization)
^^^^^^^^^^^^^^^^^^^^^

**ioptmze** = [-1|0|1|2|3|10]

::
 
   Default **ioptmze** = 2

======= =====================================
ioptmze Type of Structural Optimization
======= =====================================
**-1**  Do not optimize atomic coords
**0**   Conjugate gradient
**1**   Projected velocity Verlet
**2**   Quasi-Newton method with BFGS formula
======= =====================================

Determines method for structural optimization of atomic coordinates.
This varibale is onlt read when **ifmd** == 1.


5.7.8 (cell optimization)
^^^^^^^^^^^^^^^^^^^^^^^^^^

**ioptmze_cell** = [-1|0|1|2]

::
 
   Default **ioptmze_cell** = -1

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


5.7.9 (cell CG time step)
^^^^^^^^^^^^^^^^^^^^^^^^^^

**dtcellcg** = [real]

::
 
   Default **dtcellcg** = 0.1

only for Conjugate gradient method (**ifmd** == 1 & **ioptmze_cell** ==
0).


5.7.10 (stabilizer for quasi-Newton)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**gammamin** = [real]

::
 
   Default **gammamin** = 0.1

only for quasi-Newton method (**ifmd** == 1 & **ioptmze** == 2).


5.7.11 (clear Hessian)
^^^^^^^^^^^^^^^^^^^^^^^

**ibfgsclear** = [0|1]

::
 
   Default

**ibfgsclear** = 0

1: Clear Hessian after every ibfgsclear step 0: Hessian is not cleared


5.7.12 (clear cell Hessian)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**iclearcellh** = [integer]

::
 
   Default **iclearcellh** = 0

Clear Hessian every **iclearcellh** step. If **iclearcellh** == 0, the
Hessian is not cleared.


5.7.13 (hybrid optimization)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lhybridopt** = [boolean] **nstep_hybrid** = [integer]
**nstep_hybrid_cell** = [integer]

::
 
   Default **lhybridopt** = .TRUE. **nstep_hybrid** = 10
**nstep_hybrid_cell** = 10

**lhybridopt** = .TRUE. : Perform structural optimization first, then
cell optimization. **lhybridopt** = .FALSE. : Perform cell optimization
first, then structural optimization.

**nstep_hybrid** is the time step for structural optimization.
**nstep_hybrid_cell** is the time step for cell optimization. These
variables are only read for optimization calculations(\ **ifmd** == 1 &
**ioptmze** >= 0 & **ioptmze_cell** >= 0 )


5.7.14 (pressure)
^^^^^^^^^^^^^^^^^^

**hpext** = [real]

::
 
   Default **hpext** = 0.0

Defines the pressure in [GPa] for NPT-MD & MSST (**ifmd** == 4, 10) &
**ioptmze_cell** >= 0.


5.7.15 (barostat parameters)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**tbomega** = [real] **blkmod** = [real]

::
 
   Default **tbomega** = 5500.0 **blkmod** = 250.0

Barostat parameters for NPT (**ifmd** == 4). **tbomega** gives the time
scale for barostat in [a.u.], and **blkmod** gives the bulk modulus in
[GPa].



5.7.16 (shock wave velocity)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**shockspeed** = [real] **nshockv(1:3)** = [0|1, 0|1, 0|1]

::
 
   Default **shockspeed** = 2000.0 **nshockv(1:3)** = (/1,0,0/)

These variables are only read for MSST (**ifmd** == 10). They define the
shockspeed and the shock direction. For example:

**nshockv(1:3)** = 1 0 0: L1 direction **nshockv(1:3)** = 0 1 0: L2
direction **nshockv(1:3)** = 0 0 1: L3 direction **nshockv(1:3)** = 1 1
0: L1 + L2 direction

.. _51025-clear-barostat-velocity:

5.7.17 (clear barostat velocity)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lmsstscale** = [boolean] **msstscnum** = [integer] **msstscstp** =
[integer]

::
 
   Default **lmsstscale** = .FALSE. **msstscnum** = 20 **msstscstp** = 500

These variables are only read for MSST (**ifmd** == 10).

.. _51026-initial-barostat-velocity:

5.7.18 (initial barostat velocity)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**vxmsst** = [real]

::
 
   Default **vxmsst** = 0.0

Defines the initial barostat velocity. Only read for MSST (**ifmd** ==
10).

.. _51027-output-data:

5.7.19 (output data)
^^^^^^^^^^^^^^^^^^^^^

**ioskip** = [integer] **locoor, ioskipcoor** = [boolean, integer]

Defines what data to output for MD-nodes.

.. _51028-charge-estimation:

5.7.20 (charge estimation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**ichest** = [integer]

::
 
   Default **ichest** = 3

Number of previous steps

.. _51029-wavefunction-estimation:

5.7.21 (wavefunction estimation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**ihest** = [integer]

Number of previous steps

.. _51030-aspc-charge-estimation:


5.7.22 (tolerance)
^^^^^^^^^^^^^^^^^^^

**tol_energy** = [real] **tol_force** = [real]

tolerance for CG optimization (ifmd == 1 ). **tol_energy** given in
energy/atom in [a.u.]. **tol_force** gives maximum force in [a.u.].

.. _511-save-data:

5.8 save data
~~~~~~~~~~~~~~

.. _5.8-onoff:

5.8.1 (On/Off)
^^^^^^^^^^^^^^^

**lsave** = [boolean]

::
 
   Default **lsave** = .TRUE.

Whether or not to save data...which data?

.. _5.9-data-type:

5.8.2 (data type)
^^^^^^^^^^^^^^^^^^

**lsreal8** = [boolean]

::
 
   Default **lsreal8** = .TRUE.

.TRUE. : Saves data as real\ *8 .FALSE. : Saves data as integer*\ 2

Determines the data type for saved data



5.9 stress calculation ------ only for bulk calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _5121-onoff:

5.9.1 (On/Off)
^^^^^^^^^^^^^^^

**lstress** = [boolean]

::
 
   Default **lstress** = .TRUE.

.TRUE. : Estimate the stress .FASLE. : Do not estimate the stress

Determines whether to estimate the stress.

.. _5122-skip-step:

5.9.2 (skip step)
^^^^^^^^^^^^^^^^^^

**nskip_stress** = [integer]

::
 
   Default **nskip_stress** = 5

**nskip_stress** is skip step to estimate the stress.

.. _513-dump-charge-density:

5.10 dump charge density
~~~~~~~~~~~~~~~~~~~~~~~~

.. _5131-onoff:

5.10.1 (On/Off)
^^^^^^^^^^^^^^^

**ldpchg** = [boolean]

::
 
   Default **ldpchg** = .FALSE.

.TRUE. : Output charge density .FASLE. : Do not output

Determines whether to output the charge density.

.. _5.102-skip-step:

5.10.2 (skip step)
^^^^^^^^^^^^^^^^^^

**nskip_dpchg** = [integer]

::
 
   Default **nskip_dpchg** = 5

**nskip_dpchg** is skip step to output the charge density.

.. _5133-output-area:

5.10.3 (output area)
^^^^^^^^^^^^^^^^^^^^

**x_min** **x_max** = [real] [real] **y_min** **y_max** = [real] [real]
**z_min** **z_max** = [real] [real]

::
 
   Default **x_min** **x_max** = 1.0 0.0 **y_min** **y_max** = 1.0 0.0
**z_min** **z_max** = 1.0 0.0

These parameter specify the area to output the charge density. They are
reduced coordinates. If **[x|y|z]_min** > **[x|y|z]_max**, output are is
the whole space.

.. _514-dump-wavefunctions:

5.11 dump wavefunctions
~~~~~~~~~~~~~~~~~~~~~~~

.. _5141-onoff:

5.11.1 (On/Off)
^^^^^^^^^^^^^^^

**ldpwav** = [boolean]

::
 
   Default **ldpwav** = .FALSE.

.TRUE. : Output wavefunction .FASLE. : Do not output

Determines whether to output the wavefunction.

.. _5142-bands:

5.11.2 (bands)
^^^^^^^^^^^^^^

(**ibstt1**, **ibstt2**) = [integer] [integer]

::
 
   Default (**ibstt1**, **ibstt2**) = (/ 0, 0 /)

The wavefunction will be output from **ibstt1** to **ibstt2**.
**ibstt1** and **ibstt2** are the index of bands. If both **ibstt1** and
**ibstt2** are zero, all bands will be output.

.. _5143-skip-step:

5.11.3 (skip step)
^^^^^^^^^^^^^^^^^^

**nskip_dpwav** = [integer]

::
 
   Default **nskip_dpwav** = 5

**nskip_dpwav** is skip step to output the wavefunction.

.. _5144-output-area:

5.11.4 (output area)
^^^^^^^^^^^^^^^^^^^^

**x_min** **x_max** = [real] [real] **y_min** **y_max** = [real] [real]
**z_min** **z_max** = [real] [real]

::
 
   Default **x_min** **x_max** = 1.0 0.0 **y_min** **y_max** = 1.0 0.0
**z_min** **z_max** = 1.0 0.0

These parameter specify the area to output the wavefunction. They are
reduced coordinates. If **[x|y|z]_min** > **[x|y|z]_max**, output are is
the whole space.

.. _515-dump-local-potential:

5.12 dump local potential
~~~~~~~~~~~~~~~~~~~~~~~~~

.. _5151-onoff:

5.12.1 (On/Off)
^^^^^^^^^^^^^^^

**ldppot** = [boolean]

::
 
   Default **ldppot** = .FALSE.

.TRUE. : Output local potential .FASLE. : Do not output

Determines whether to output the local potential.

.. _5152-average-plane:

5.12.2 (average plane)
^^^^^^^^^^^^^^^^^^^^^^

**nav_dppot** = [1|2|3]

::
 
   Default **nav_dppot** = 1

============= =====
**nav_dppot** Plane
============= =====
**1**         xy
**2**         yz
**3**         zx
============= =====

The plane-averaged local potential will be output.

.. _5153-skip-step:

5.12.3 (skip step)
^^^^^^^^^^^^^^^^^^

**nskip_dppot** = [integer]

::
 
   Default **nskip_dppot** = 5

**nskip_dppot** is skip step to output the local potential.

.. _5154-output-area:

5.12.4 (output area)
^^^^^^^^^^^^^^^^^^^^

**x_min** **x_max** = [real] [real] **y_min** **y_max** = [real] [real]
**z_min** **z_max** = [real] [real]

::
 
   Default **x_min** **x_max** = 1.0 0.0 **y_min** **y_max** = 1.0 0.0
**z_min** **z_max** = 1.0 0.0

These parameter specify the area to output the local potential. They are
reduced coordinates. If **[x|y|z]_min** > **[x|y|z]_max**, output are is
the whole space.

.. _5.5-supercell:

5.13 supercell
~~~~~~~~~~~~~~

.. _5.51-unit-of-length:

5.13.1 (unit of length)
^^^^^^^^^^^^^^^^^^^^^^^

(bohr) or (ang)

::
 
   Default (bohr)

Determines the unit of supercell (bohr or angstrom). And the supercell
is specified by either cell vector or 'lengths and angles'. You need
comment out the subsection you do not specify. Note that for cluster
calculations, the supercell has to be ORTHORHOMBIC.

.. _5.52-cell-vector:

5.13.2 (cell vector)
^^^^^^^^^^^^^^^^^^^^

L1 = [real] [real] [real] L2 = [real] [real] [real] L3 = [real] [real]
[real]

Example:

::

    10.68   0.0    0.0        : super cell vector L1
     0.0   10.68   0.0        : super cell vector L2
     0.0    0.0   10.68       : super cell vector L3

.. _5.53-lengths--angles:

5.13.3 (lengths & angles)
^^^^^^^^^^^^^^^^^^^^^^^^^

Example:

::

    10.68,  10.68,  10.68  :  lengths of cell vectors
    90.00,  90.00,  90.00  :  angles between cell vec. in [deg.]


5.14 planewaves
~~~~~~~~~~~~~~~

.. _5.61-unit-of-cutoff-energy:

5.14.1 (unit of cutoff energy)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

(ry) or (hr) or (ev)

::
 
   Default (ry)

Determines the unit of cutoff energy.

.. _5.62-for-wavefuctions:

5.14.2 (for wavefuctions)
^^^^^^^^^^^^^^^^^^^^^^^^^

**ecut** = [real]

**ecut** is the cutoff energy to expand wavefunctions.

.. _5.63-for-electron-density:

5.14.3 (for electron density)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**ecutdens** = [real]

::
 
   Default **ecutdens** = 0.0

**ecutdens** is the cutoff energy to expand electron density.
**ecutdens** must be greater than **ecutsoft**. if .not.lvand,
**ecutdens** must be equal to **ecutsoft**.

.. _5.64-for-soft-part-of-density:

5.14.4 (for soft part of density)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**ecutsoft** = [real]

::
 
   Default **ecutorth** = 0.0

**ecutsoft** is the cutoff energy to expand soft part of electron
density. It must be greater than ecut, and smaller than ecut*4.

.. _5.65-for-orthogonalization:

5.14.5 (for orthogonalization)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**ecutorth** = [real]

::
 
   Default **ecutorth** = 0.0

**ecutorth** is the cutoff energy for orthogonalization in CG iteration.

.. _5.66-for-charge-density-mixing:

5.14.5 (for charge density mixing)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**ecutc** = [real]

::
 
   Default **ecutc** = 0.0

**ecutc** is the cutoff energy for charge density mixing.




5.15 electronic bands
~~~~~~~~~~~~~~~~~~~~~

.. _5.71-occupied-bands:

5.15.1 (occupied bands)
^^^^^^^^^^^^^^^^^^^^^^^

**noband** = [integer]

::
 
   Default **noband** = 0

The number of occupied bands should at least be sufficient to
accommodate all the valence electrons in the simulation.

.. _5.72-empty-bands:

5.15.2 (empty bands)
^^^^^^^^^^^^^^^^^^^^

**neband** = [integer]

::
 
   Default **neband** = 0

A good rule of thumb is to set neband equal to 10% of occupied bands.

.. _5.73-broadening:

5.15.3 (broadening)
^^^^^^^^^^^^^^^^^^^

**lfermi** = [integer] flag **tfermi** = [real]

::
 
   Default **lfermi** = 3 **tfermi** = 2000.0

lfermi = 1 for Non-metallic systems, lfermi = 2 for Fermi smearing,
lfermi = 3 for Gaussian smearing and lfermi > 3 for Methfessel-Paxton
smearing scheme of order lfermi-3

.. _5.74-charge-number-of-ion:

5.15.4 (charge number of ion)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**ncion** = [integer]

::
 
   Default **ncion** = 0. For calculation of electronically non-neutral
systems. The number of electrons is given by ``nel-ncion``.
Charged-periodic calculations are perfomed with the uniform background.


.. _5.8-atoms:

5.16 Atoms
~~~~~~~~~~

.. _5.81-species:

5.16.1 (species)
^^^^^^^^^^^^^^^^

**ntype** = [integer]

::
 
   Default **ntype** = 0

Number of atom types. The subsequent sections must be duplicated
``ntype`` times.

.. _5.82-atomic-number:

5.16.2 (atomic number)
^^^^^^^^^^^^^^^^^^^^^^

**zatom** = [integer]

::
 
   Default **zatom** = No default

Atomic number of the current atom type

.. _5.83-pseudopotential:


5.16.3 (partial core correction)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lpcc** = [boolean] **rpcc** = [real] **lpking** = [boolean] **rpking**
= [real] **gpkgmax** = [real] **gpkgexct** = [real]

::
 
   Default **lpcc** = No default **rpcc** = No default **lpking** = No
default **rpking** = No default **gpkgmax** = No default **gpkgexct** =
No default

Psuedopotential formalism treat core and valence electron seperately.
Partial core correction includes the effect of core charge as a
perturbation which increase the transferrability of the code. It is
recommended to use partial core correction. Please see paper by Louie et
al.

.. _5.87-compensation-charge-cutoff:

.. _5.88-the-number-of-atoms:

5.16.4 (the number of atoms)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**nhk** = [integer]

::
 
   Default **nhk** = 0

Number of atoms of the current type. If ``nhk`` is zero, then atomic
configuration is read from the positions file

.. _5.89-unit-of-length:

5.16.5 (unit of length)
^^^^^^^^^^^^^^^^^^^^^^^

(bohr) or (ang)

::
 
   Default No default

.. _5.810-position-file:

5.16.6 (position file)
^^^^^^^^^^^^^^^^^^^^^^^

**fname** = [String] **scaledflag**\ \* = [integer] **keyword**\ \* =
[integer]

::
 
   Default **fname** = No default **scaledflag**\ \* = No default
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

.. _5.811-positions:

5.16.7 (positions)
^^^^^^^^^^^^^^^^^^^

**nunitcells** = [integer][integer][integer] **scaledflag**\ \* =
[integer] ``nhk`` lines of **positions** = [real][real][real]

::
 
   Default **nunitcells** = No defaults **scaledflag**\ \* = No defaults

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

.. _5.812-fix-positions:


5.16.8 (velocity file)
^^^^^^^^^^^^^^^^^^^^^^^

**fname** = [String] **keyword**\ \* = [integer]

::
 
   Default **fname** = No default **keyword**\ \* = No default

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

.. _5.816-velocities:

5.16.9 (velocities)
^^^^^^^^^^^^^^^^^^^^

**nunitcells** = [integer][integer][integer] **scalingfactor**\ \* =
[real] ``nhk`` lines of scaled **velocities** = [real][real][real] .
Ensure that -1.0 < **velocities** < 1.0

::
 
   Default **nunitcells** = No defaults **scalingfactor**\ \* = No
defaults

**nunitcells** gives the number of unit cells along the x, y and z
directions. Atomic velocities are initialized by multiplying the vector
**velocities** by the **scalingfactor**.

.. _5.817-spherical-harmonics-expansion:

.. bibliography:: param.bib
   :style: unsrt



.. raw:: html

   <br>


