Input/Output files 
==================


4.1 Input Files
~~~~~~~~~~~~~~~

4.1.1 filename
^^^^^^^^^^^^^^

**filename** is a simple, single line text file containing the path to
the main input file, relative to the QXMD executable, in between single
quotes. Below shows the entire contents of a sample **filename**

4.1.2 CONFIG
^^^^^^^^^^^^^^^

An initial configuration file must be included for any new QXMD run. It
is best practice to name this file **CONFIG**. Below shows a sample
configuration file for a water molecule.

::

   3
   1 0.579029 0.305583 0.264540
   1 0.420969 0.305577 0.264541
   2 0.500002 0.244625 0.264685

The first line should contain the total number of atoms in the
simulation system. In this case, the water molecule has 3 atoms: 2
hydrogen and 1 oxygen. Each of the following lines (one line for each
atom) contains a keyword representing the species the atom, as well as
the atom's **x**, **y**, and **z** coordinates. These may be given in
real space of fractional coordinates (you define which kind are provided
in the IN.PARAM file). In the example above, the second and third lines
give the spatial coordinates for the hydrogen atoms, using '1' as the
keyword for hydrogen, while the last line give the spatial coordinates
for the oxygen atom, using '2' as the keyword for oxygen. Integer
numbers should be used for keywords, though it is arbitrary which
numbers are chosen for each species, as long as what you choose is used
consistently with atomic data entered in the main input file, IN.PARAM.

4.1.3 IN.PARAM
^^^^^^^^^^^^^^

A main input file is required for all QXMD runs, which contains all the
parameters settings. It is best practice to name the main input file
**IN.PARAM**. There are on the order of a few hundred input parameters
that may be tailored for various kinds of QXMD simulations, though most
parameters have default settings, many of which need not be changed for
most routine QXMD simulations. Details of these parameters may be found
in section 5.

4.1.4 VELOC
^^^^^^^^^^^^^^

A file defining the intial velocities of the atoms may optionally be
provided. It is best practice to name this file **VELOC**. If this is
not provided, random inital velocities will be assigned to each atom, in
accordance with the initial system temperature. Below shows a sample
initial velocity file for a water molecule.

::

   3
     2.7271510E-03
    1  0.058157  0.043210  0.003117
    2 -1.000000 -0.688863 -0.057829
    2  0.076852  0.002967  0.008360

The first line should contain the total number of atoms in the
simulation system. In this case, the water molecule has 3 atoms: 2
hydrogen and 1 oxygen. The next line is a scaling factor by which all of
the following veloctiy components should be multiplied. Each of the
following lines (one line for each atom) contains a keyword representing
the species the atom, as well as the atom's **x**, **y**, and **z**
velocity components in atomic units. The integer keywords for the atomic
species should be consistent with those provided in the IN.CONFIG and
IN.PARAM files.


4.2 Output Files
~~~~~~~~~~~~~~~~

All output data files are dumped to the **data** directory during QXMD
runs. Output files with names beginning with a lower case letter are
human-readable text files detailing simulation data, while output files
beginning with upper case letters are binary files to be used for
restarting a QXMD run. These files must be present in the **data**
directory when attempting to restart a job. After a simulation, all
output files that you wish to save should be moved out of the **data**
directory and store elsewhere as a new QXMD run will write over all
files present in **data**.

4.2.1 md_box.d
^^^^^^^^^^^^^^

This file gives the simulation box size for a classical MD simulation. A
sample of output file md_box.d is shown below.

::

   # MD cell vectors
   #        L_1           L_2           L_3          angle(2-3) angle(3-1) angle(1-2)
         0 3.7280514E+01 3.7280514E+01 3.0613561E+01  90.000000  90.000000 120.000000

L_1, L_2, and L_3 are the lengths of three cell vectors in Bohr units,
while angles between these vectors, in degrees, follow. The '0' at the
beginning of the last line represents the MD step number. In this
example, the simulation box size is held fixed, and thus this is the
only line that appears in the output. If these numbers were updated at
any point, a new line would appear under the last line with the step
number in which the box size was altered, followed by the new box
dimensions.

Note that there is also an output file named **qm_box.d**, which
provides the supercell size for a quantum MD simulation.


4.2.2 md_cel.d
^^^^^^^^^^^^^^

Similar to **md_box.d**, **md_cel.d** holds the cell vector data for the
simulation cell in which molecular dynamics take place, however instead
of giving the lengths and angles of the vectors, this file gives the
**x**, **y**, and **z** components of the three cell vectors. A sample
of output file md_cel.d is shown below.

::

   # MD cell vectors
   #        L_(1:3,1:3)
         0  3.7280514E+01  0.0000000E+00  0.0000000E+00 -1.8640257E+01  3.2285872E+01  0.0000000E+00  1.8745400E-15  3.2467985E-15  3.0613561E+01

After the two comment lines, the file lists the step number, in this
case '0', and then **x**, **y**, and **z** components of the first cell
vector, followed by those of the second and third cell vectors. In this
case, there were no changes in the cell vectors during the simulation,
however, if changes do occur, additional lines will be written beginning
with the step number in which the cell vectors were changed, followed by
their new values.

4.2.3 md_eng.d
^^^^^^^^^^^^^^

This file contains the Hamiltonian, potential, and kinetic energies
[Hartree] of the total system, as well as the temperature [K] at every
time step. A sample of the output is shown below.

::

   # step    H [hartree]     P.E. [hartree]  K.E. [hartree]     T[K]
         0 -1.4792707116E+04 -1.47928610E+04  1.53902180E-01   300.0000
         1 -1.4792706554E+04 -1.47928611E+04  1.54538843E-01   301.2410
         2 -1.4792705828E+04 -1.47928609E+04  1.55091556E-01   302.3184
         3 -1.4792704962E+04 -1.47928605E+04  1.55563674E-01   303.2387
         4 -1.4792703950E+04 -1.47928599E+04  1.55970999E-01   304.0327
         5 -1.4792703306E+04 -1.47928596E+04  1.56301244E-01   304.6765
         6 -1.4792702504E+04 -1.47928591E+04  1.56548702E-01   305.1588
         7 -1.4792701327E+04 -1.47928580E+04  1.56710424E-01   305.4741
         8 -1.4792700224E+04 -1.47928570E+04  1.56777438E-01   305.6047
         9 -1.4792699401E+04 -1.47928562E+04  1.56763969E-01   305.5785
        10 -1.4792698443E+04 -1.47928551E+04  1.56675161E-01   305.4053


4.2.4 md_log
^^^^^^^^^^^^

This file provides of log of the simulation, including values of input
parameters set and lengths of computation time for force calculations
for each time step.

4.2.5 MD_mts0
^^^^^^^^^^^^^

This is a binary file that is required for restarting a QXMD simulation.

4.2.5 MD_mts0
^^^^^^^^^^^^^

This is a binary file that is required for restarting a QXMD simulation.


4.2.6 md_spc.d
^^^^^^^^^^^^^^

This file lists the species keyword for each atom at each time step. A
sample of the output is shown below.

::

   #  Atomic species               
         2    42  34
         0     12
    1 1 1 1 2 2 2 2 2 2 2 2
         1     12
    1 1 1 1 2 2 2 2 2 2 2 2
         2     12
    1 1 1 1 2 2 2 2 2 2 2 2
         3     12
    1 1 1 1 2 2 2 2 2 2 2 2
         4     12
    1 1 1 1 2 2 2 2 2 2 2 2
         5     12
    1 1 1 1 2 2 2 2 2 2 2 2

In the second line, '2' represents the total number of species in the
system. The following numbers,'42' and '34', are the atomic numbers of
the two atomic species in the system ordered by the keywords
representing those species. In this example, molybdenum has atomic
number 42, and is represented by keyword '1', while selenium has atomic
number 34, and is represented by keyword '2'. The next line lists the
step number, followed by the total number of atoms in the system, in
this case '12'.


4.2.7 md_str.d
^^^^^^^^^^^^^^

This file outputs components of the stress tensor [GPa] of the system,
computed classically, at given time steps intervals (the interval is
defined in IN.PARAM). This file will only be written if dumping stress
data is set to **true** in IN.PARAM. Below is sample output.

::

   # Stress in [GPa]               
   #        Pxx             Pyy             Pzz            Pyz             Pzx             Pxy     
   0  1.33846254E+01  1.32248154E+01  5.38055547E+00  2.74846400E-02  4.09607470E-02 -9.31727041E-01
   5  1.37099707E+01  1.34067898E+01  6.81653796E+00  1.21177457E-02  2.32549066E-01 -1.03833531E+00
   10  1.52252127E+01  1.57211516E+01  9.64813968E+00  3.34044020E-02  3.45617372E-01 -8.50283630E-01
   15  1.68630437E+01  1.59302869E+01  1.19623442E+01  1.09971448E-01  4.41730570E-01  7.93058521E-01
   20  1.78613247E+01  1.62869327E+01  1.36732341E+01  1.16372513E-01  4.25461913E-01  1.08703422E+00

After the first two comment lines, the following lines list the step
number followed by the components of the stress tensor for that time
step. Here, stress data is dumped every 5 time steps. Note that
contributions from the kinetic energy of the ions are included in these
values, while in qm_str.d they are not.


4.2.8 md_str_diag.d
^^^^^^^^^^^^^^^^^^^

This file outputs the diagonalized stress tensor, along with the
respective eigenvectors at given time step intervals. Below is sample
output.

::

   # Diagonalized stress in [GPa]  
   #        Pxx             Pyy             Pzz            Eigenvectors             
         0  1.42398826E+01  1.23699022E+01  5.38021140E+00  7.3672E-01 -6.7620E-01  1.3084E-03  6.7618E-01  7.3671E-01  6.8597E-03 -5.6024E-03 -4.1690E-03  9.9998E-01
         5  1.46113490E+01  1.25135917E+01  6.80835783E+00  7.5737E-01 -6.5263E-01  2.1581E-02  6.5207E-01  7.5764E-01  2.8228E-02 -3.4773E-02 -7.3072E-03  9.9937E-01
        10  1.46051510E+01  1.63637815E+01  9.62557143E+00  7.9488E-01  6.0384E-01  5.9491E-02 -6.0340E-01  7.9698E-01 -2.7090E-02 -6.3770E-02 -1.4364E-02  9.9786E-01
        15  1.73523452E+01  1.54808731E+01  1.19224564E+01  8.6781E-01  4.9024E-01  8.1122E-02 -4.8915E-01  8.7153E-01 -3.4170E-02 -8.7452E-02 -1.0028E-02  9.9612E-01
        20  1.84552850E+01  1.57357781E+01  1.36304284E+01  8.8838E-01  4.5019E-01  8.9996E-02 -4.4823E-01  8.9293E-01 -4.2081E-02 -9.9304E-02 -2.9550E-03  9.9505E-01

After the first two comment lines, the following lines list the step
number followed by the components of the diagonalized stress tensor and
the **x**, **y**, and **z** components of the three eigenvectors for
that time step. Here, diagonalized stress data is dumped every 5 time
steps.


4.2.9 MD_Tocontinue
^^^^^^^^^^^^^^^^^^^

This is a binary file that is required in the **data** directory to
restart a QXMD simulation.


4.2.10 md_vel.d
^^^^^^^^^^^^^^^

This file contains scaled atomic velocities for each atom at each time
step. Below is sample output.

::

   #  Atomic scaled velocities     
         0     12
    1.7316159E-06
   -6.94686 -2.92359 -1.82411  5.50162 -1.28242  2.31888  0.35217 -0.34956  1.04083
   -2.29616  2.29092  0.53871 -0.03100 -5.38400 -1.78311  2.52055  1.50530  1.29961
    1.94346  1.81671 -2.48381  5.25972  2.70791  1.54716  5.36706  8.68984  1.05675
   -5.95609 -9.99900 -0.80332 -3.99628  1.60613 -0.32696 -0.98935  1.80875 -1.02671

The second line gives the step number, in this case '0', and the total
number of atoms in the system, in this case, '12'. The next line gives
the number by which to scale the components of velocity given in the
following lines (i.e. multiply the following numbers to get the absolute
velocity values). The following lines give the **x**, **y**, and **z**
components of velocity for the first atom, followed by those for the
second atom, ending with those for the last atom.


4.2.11 qm_box.d
^^^^^^^^^^^^^^^

This file holds the cell vector data for the supercell for a quantum
molecular dynamics simulation. A sample of output file is shown below.

::

   #  supercell (FFT cell) vectors (lengths & angles)
   #        L_1           L_2           L_3          angle(2-3) angle(3-1) angle(1-2)                   
         0 1.2426838E+01 1.2426838E+01 3.1114338E+01  90.000000  90.000000 120.000000

L_1, L_2, and L_3 are the lengths of three supercell vectors in Bohr
units, while angles between these vectors, in degrees, follow. The '0'
at the beginning of the last line represents the step number. In this
example, the simulation box size is held fixed, and thus this is the
only line that appears in the output. If the supercell size were
changed, a new line would appear under the last line with the step
number in which the size was altered, followed by the new supercell
dimensions.

Note that there is also an output file named **md_box.d**, which holds
the box size for a classical MD simulation. For purely quantum MD
simulations, **md_box.d** and **qm_box.d** will hold the same
information (except in cases where the double-grid method is used).
However, it is possible to perform a hybrid classical-quantum MD
simulation in which case the system is treated classically, except for a
defined area inside which is treted quantum mechanically. In this case,
**md_box.d** gives the box size of the entire system, while **qm_box.d**
defines the supercell area that is to be treated with QM.


4.2.12 QM_cds
^^^^^^^^^^^^^

This is a binary file that is required in the **data** directory to
restart a QXMD simulation, containing charge density data for the system
at the last step of the previous run.


4.2.13 qm_cel.d
^^^^^^^^^^^^^^^

Similar to **qm_box.d**, **qm_cel.d** holds the cell vector data for the
supercell, however instead of giving the lengths and angles of the
vectors, this file gives the **x**, **y**, and **z** components of the
three supercell vectors. Below is sample output.

::

   #  supercell (FFT cell) vectors
   #        L_(1:3,1:3)
         0  1.2426838E+01  0.0000000E+00  0.0000000E+00 -6.2134191E+00  1.0761957E+01  0.0000000E+00  1.9052037E-15  3.2999097E-15  3.1114338E+01

After the two comment lines, the file lists the step number, in this
case '0', and then **x**, **y**, and **z** components of the first
supercell vector, followed by those of the second and third supercell
vectors. In this case, there were no changes in the supercell vectors
during the simulation, however, if changes do occur, additional lines
will be written beginning with the step number in which the supercell
vectors were changed, followed by their new values.


4.2.14 QM_cell
^^^^^^^^^^^^^^

This is a binary file that is required in the **data** directory to
restart a QXMD simulation, containing simulation cell data for the
system at the last step of the previous run.


4.2.15 QM_eig
^^^^^^^^^^^^^

This is a binary file that is required in the **data** directory to
restart a QXMD simulation, containing energy eigenvalue data for the
system at the last step of the previous run.


4.2.16 qm_eig.d
^^^^^^^^^^^^^^^

This file lists the energy eigenvalues, along with their electronic
occupation number for each time step. Below is sample output.

::

   #  Eigenvalues                   
         0      7     10
        1 -1.79706E+00 2.000
        2 -8.95274E-01 2.000
        3 -6.37670E-01 2.000
        4 -4.86501E-01 2.000
        5 -9.48240E-02 0.000
        6  1.26780E-01 0.000
        7  1.34093E-01 0.000
        8  1.69431E-01 0.000
        9  1.94460E-01 0.000
       10  2.79168E-01 0.000
         1     11     10
        1 -1.79215E+00 2.000
        2 -8.86033E-01 2.000
        3 -6.40897E-01 2.000
        4 -4.85737E-01 2.000
        5 -9.72281E-02 0.000
        6  1.20947E-01 0.000
        7  1.34115E-01 0.000
        8  1.68074E-01 0.000
        9  1.93969E-01 0.000
       10  2.75769E-01 0.000

The first number in the second line gives the step number, in this case
'0'. The second number gives the cumulative number of SCF iterations
performed, in this case '7' iterations were performed in step 0.
Finally, the third number, '10' represents the total number of energy
eigenvalues and occupation numbers to follow. This number corresponds to
the number of energy bands for the system, defined in IN.PARAM. Thus the
next ten lines list the band index number, the energy in eV, and the
number of electrons which occupy that energy band. In the second step, 4
SCF iterations were performed since the number '11' follows the step
number '1' (7 SCF iterations were performed in the first step, followed
by 4 SCF iterations, for a total of 11 SCF iterations after the first
two time steps.) In this example, the eight electrons of the system
occupy the four lowest energy bands in the first time two steps.


4.2.17 qm_eng.d
^^^^^^^^^^^^^^^

This file provides various components of the system's energy [Rydberg]
at each time step. Below is sample output.

::

   #  Total potential energy and energy parts in [Ryd.] units
   #              Total(HF)         Total(KS)         Kinetic      External      Hartree       Exchange      Correlation   ------        Entropy       Onsite E.     --------      ------        Ewald E.      DFT-D
        0     7   -4.4003610122E+01 -4.4003640357E+01  1.370805E+01 -6.080189E+01  2.800936E+01 -7.786772E+00 -5.786949E-01  0.000000E+00  0.000000E+00 -1.533928E+01  0.000000E+00  0.000000E+00 -1.214255E+00 -1.610012E-04
        1     11   -4.4002961457E+01 -4.4003371148E+01  1.367946E+01 -6.058249E+01  2.791215E+01 -7.771542E+00 -5.778955E-01  0.000000E+00  0.000000E+00 -1.532657E+01  0.000000E+00  0.000000E+00 -1.336327E+00 -1.631839E-04
        2    15   -4.4001366424E+01 -4.4002003283E+01  1.365563E+01 -6.038887E+01  2.782594E+01 -7.758172E+00 -5.771782E-01  0.000000E+00  0.000000E+00 -1.531663E+01  0.000000E+00  0.000000E+00 -1.442551E+00 -1.647241E-04
        3    18   -4.3999479244E+01 -4.3999841762E+01  1.363668E+01 -6.023643E+01  2.775745E+01 -7.747680E+00 -5.766015E-01  0.000000E+00  0.000000E+00 -1.530791E+01  0.000000E+00  0.000000E+00 -1.525182E+00 -1.650145E-04
        4    22   -4.3997932717E+01 -4.3997243098E+01  1.362538E+01 -6.013514E+01  2.771149E+01 -7.740722E+00 -5.762123E-01  0.000000E+00  0.000000E+00 -1.530223E+01  0.000000E+00  0.000000E+00 -1.579643E+00 -1.636864E-04
        5    25   -4.3997129513E+01 -4.3997300982E+01  1.362022E+01 -6.009034E+01  2.769035E+01 -7.737485E+00 -5.760705E-01  0.000000E+00  0.000000E+00 -1.530042E+01  0.000000E+00  0.000000E+00 -1.603399E+00 -1.606360E-04

After the first two comment lines, each line gives the step number, the
cumulative number of SCF iterations completed up to that step number,
followed by various components of system energy as labeled by the column
titles.


4.2.18 qm_fer.d
^^^^^^^^^^^^^^^

This file gives the Fermi energy of the system at each time step. Below
is sample output.

::

   #  Fermi energy
         0      7 -2.90493E-01
         1     11 -2.90662E-01
         2     15 -2.91483E-01
         3     18 -2.92223E-01
         4     22 -2.92859E-01
         5     25 -2.93001E-01

The first column gives the step number, the second column gives the
cumulative number of SCF iterations up to that step number, while the
third column gives the Fermi energy in eV.


4.2.19 qm_frc.d
^^^^^^^^^^^^^^^

This file gives the three components of force on each atom at each time
step. Below is sample output for a monolayer of MoSe2 with 12 total
atoms.

::

   #  Atomic forces in [a.u.]       
         0      2      4      8
    3.7133826E-02
    0.14634 1.41724-0.01118-0.15415 0.36847-0.00294 0.87782-0.64926 0.00080
   -0.89534-0.92856 0.00524 1.28410 0.95522 8.33389-1.27916 1.17417 9.99839
   -0.17476-0.37722 6.37526 0.18425-1.85265 8.15680 1.29026 0.94518-8.33006
   -1.28567 1.15640-9.99900-0.17251-0.36263-6.37149 0.17883-1.84636-8.15571
         1      2      4      8
    3.5751399E-02
   -1.40927-0.53578 0.00521 1.21420-0.62553 0.03759 0.23010-0.18164-0.00182
   -0.27261 1.29953 0.02477-0.13760-1.72196 9.40196 0.26416 0.31938 6.59092
   -1.83555 0.61826 8.18927 1.83320 0.80236 9.98594-0.10811-1.77374-9.39287
    0.19630 0.32076-6.62970-1.84943 0.69302-8.21226 1.87463 0.78534-9.99900

After the comment line, the next line gives the step number, the total
number of atomic species in the system ('2' for molybdenum and
selenium), the number of atoms of the species corresponding to keyword
'1' (in this case, keyword '1' was used for molybdenum and there are '4'
molybdenum atoms), followed by the number of atoms of the species
corresponding to keyword '2' (in this case, keyword '2' was used for
selenium and there are '8' selenium atoms). The number in the third line
is a scaling factor for the following force vector components (multiply
all of the following force components by this number to get the true
values for force). The following lines give the **x**, **y**, and **z**
components of force on the first atom, followed by those for the second
atom, and so on.


4.2.20 qm_fsshprob_***to***-u.d
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This set of files are only written during a Non-Adiabatic QMD simulation
(i.e. TD-DFT set to .true.). These files give the probabilities for
electrons to hop from one band to another, where these band indices will
be given in the title of the file in place of the '***'. A sample file
name from this set could be qm_fsshprob_29to32-u.d, which will give the
probability for an electron to transition from band index 29 to band
index 32 at each time step that has a non-zero probability. The '-u' at
the end of the file name indicates that the data is for spin-up
electrons in the case that spin polarization is used in the QXMD
simulation. In this case, spin-down electron data will be stored in
files named 'qm_fsshprob_***to***-d.d' If spin polarization is not used,
all electron transition probabilities will be stored in the spin-up data
files by default. Below is a sample of this output from this file.

::

   # step  probability  accumulation
       14  8.19528E-07  8.19528E-07
       16  4.54580E-06  4.54580E-06
       17  1.01866E-05  1.47324E-05
       18  5.77833E-06  2.05108E-05

As the column titles suggest, the first column gives the step number at
which there exists a finite probability, given by the second column, of
the electron hopping from, in this case band index 29 to 32.
Probabilities are only written for step numbers for which there is a
finite (non-zero) probability. Thus, in the first 13 steps of the
simulation there was zero probability for the electron to hop from band
29 to 32, as well as time step 15. However, at time steps 14, 16, 17,
and 18, there were non-zero probabilities for this electronic
transition. The accumulation column refers to ???


4.2.21 QM_hrt
^^^^^^^^^^^^^

This is a binary file that is required in the **data** directory to
restart a QXMD simulation.


4.2.22 QM_ion
^^^^^^^^^^^^^

This is a binary file that is required in the **data** directory to
restart a QXMD simulation, containing atomic position data for the
system at the last step of the previous run.


4.2.23 qm_ion.d
^^^^^^^^^^^^^^^

This files gives the three components of the positions for all the atoms
at each time step. It can be used for visualization of the atomic
trajectories throughout the simulation. Below is sample output for
monolayer MoSe2 with a total of 12 atoms.

::

   #  Atomic scaled coordinates     
         0      2      4      8
    1.0000000E-01
    1.66624 3.33424 2.49974 6.66751 3.33318 2.49997 1.66758 8.33529 2.50021
    6.66486 8.33089 2.50009 3.33548 1.66849 1.46790 8.33121 1.66281 1.46807
    3.33216 6.66685 1.46770 8.33431 6.66561 1.46808 3.33503 1.66953 3.53207
    8.33384 1.66623 3.53173 3.33184 6.66746 3.53228 8.33294 6.66581 3.53215
         1      2      4      8
    1.0000000E-01
    1.66036 3.33190 2.49817 6.67225 3.33212 2.50197 1.66794 8.33491 2.50111
    6.66273 8.33275 2.50055 3.33569 1.66400 1.46680 8.33330 1.66428 1.46971
    3.33378 6.66836 1.46589 8.33873 6.66766 1.46984 3.33989 1.67716 3.53255
    8.32863 1.65779 3.53052 3.32835 6.66879 3.53167 8.33197 6.66709 3.53084

After the comment line, the next line gives the step number, the total
number of atomic species in the system ('2' for molybdenum and
selenium), the number of atoms of the species corresponding to keyword
'1' (in this case, keyword '1' was used for molybdenum and there are '4'
molybdenum atoms), followed by the number of atoms of the species
corresponding to keyword '2' (in this case, keyword '2' was used for
selenium and there are '8' selenium atoms). The number in the third line
is a scaling factor for the following spatial coordinates of each atom
(multiply all of the following spatial coordinates by this number to get
the true atomic positions). The following lines give the **x**, **y**,
and **z** coordinates of the first atom, followed by those for the
second atom, and so on.


4.2.24 qm_log
^^^^^^^^^^^^^

This file provides of log of simulation details, including values of
input parameters set, time statistics and energies computed for each SCF
iteration for every simulation time step, as well as total computation
time. This file is the best place to look to determine why a QXMD
simulation fails/crashes, and general debugging.


4.2.25 qm_mul.d
^^^^^^^^^^^^^^^

This file provides Mulliken analysis data. Note that this file will only
be written if Mulliken analysis is set to .true. in IN.PARAM. Sample
output for a water molecule is shown below.

::

   #  Mulliken analysis : s, p, & d population for each atom      
         0      2      2      1     -1
     2   0 1
     1   0
       1  1.8974  5.5777  7.4751
       2  0.2624  0.2624
       3  0.2625  0.2625

After the first comment line, the second line gives the step number, the
total number of atomic species in the system ('2' for hydrogen and
oxygen), the number of atoms of the species corresponding to keyword '1'
(in this case, keyword '1' was used for hydrogen and there are '2' H
atoms), followed by the number of atoms of the species corresponding to
keyword '2' (in this case, keyword '2' was used for oxygen and there is
'1' O atom). The '-1' at the end of the line may be ignored as is it
only included for backward compatibility with older versions of the
code. The next set of two lines gives the number of different angular
momenta (orbital types) each species has. Since oxygen was assigned
keyword '1', the first of these lines (line 3) starts with '2'
indicating that oxygen has 2 types of orbitals which are identified as
'0' (meaning 's-type') and '1' (meaning 'p-type'). The second of these
line corresponds to hydrogen, which only has one orbtial, '0' or
's-type'. Finally, the partial charges in each orbital type, followed by
the total partial charge around each atom is given in the last three
lines. For example, in line 5 of the sample output, '1' refers to atom
1, in this case an oxygen atom, the next number is the partial charge in
the s-type orbital, the next number is the partial charge in the p-type
orbitals, while the last number is the total partial charge around the
oxygen atom. The final two lines give the partial charge in the s-type
orbital in atoms '2' and '3', the two H atoms, along with the total
partial charge around the H atoms which is the same as the s-type
orbital partial charge, since there is only one orbital.


4.2.26 qm_ovp.d
^^^^^^^^^^^^^^^

This file provides the overlap charge densities between atoms at the
given time step intervals. Note that this file is only written if
Mulliken analysis is set to **true**. Below is sample output for a
monolayer of MoSe2 with 12 total atoms.

::

   #  Mulliken analysis : overlap population between atoms        
         0      2      4      8     -1
     3   0 1 2
     3   0 1 2
       9
   2 -0.0137 -0.0003  0.0171 -0.0002 -0.0983 -0.0017  0.0156 -0.0018 -0.2547
   3 -0.0057 -0.0003  0.0163 -0.0002 -0.2605  0.0009  0.0143  0.0009 -0.0488
   4 -0.0457 -0.0003  0.0207 -0.0003 -0.0089  0.0004  0.0232  0.0005 -0.0259
   5  0.0041  0.0555 -0.0021 -0.0127  0.0859 -0.0138 -0.0024  0.2150  0.0596
   6  0.0067  0.0514 -0.0028 -0.0130  0.1780 -0.0220 -0.0024  0.2105  0.0451
   7  0.0045  0.0520 -0.0027 -0.0126  0.1058 -0.0041 -0.0024  0.2361  0.0196
   9  0.0051  0.0562 -0.0024 -0.0125  0.0870 -0.0137 -0.0019  0.2140  0.0593
   10  0.0050  0.0513 -0.0024 -0.0126  0.1794 -0.0219 -0.0025  0.2126  0.0452
   11  0.0068  0.0523 -0.0026 -0.0128  0.1068 -0.0042 -0.0020  0.2331  0.0193

After the comment line, the next line gives the step number, the total
number of atomic species in the system ('2' for molybdenum and
selenium), the number of atoms of the species corresponding to keyword
'1' (in this case, keyword '1' was used for molybdenum and there are '4'
molybdenum atoms), followed by the number of atoms of the species
corresponding to keyword '2' (in this case, keyword '2' was used for
selenium and there are '8' selenium atoms). The '-1' at the end of the
line may be ignored as is it only included for backward compatibility
with older versions of the code. The next set of two lines give the
total number of angular momenta (orbital types) for each species, in
this case both Mo and Se have s-type ('0'), p-type ('1'), and d-type
('2') orbitals. The '9' in the next line signifies that the following
lines will give partial charge overlaps between atoms number 9 and the
other atoms. Lines 6-14 begin the the atom number, and then give the
overlap between the all combinations of orbital types of atom 9 and all
orbital types of the line's respective atom number. Line 6, for example,
gives the partial charge overlaps between atom 9 and atom 2. The 9
charge overlap numbers give overlaps between: s-type(atom 9)/s-type(atom
2), s-type(atom 9)/p-type(atom 2), s-type(atom 9)/d-type(atom 2),
p-type(atom 9)/s-type(atom 2), p-type(atom 9)/p-type(atom 2),
p-type(atom 9)/d-type(atom 2), d-type(atom 9)/s-type(atom 2),
d-type(atom 9)/p-type(atom 2), d-type(atom 9)/d-type(atom 2).


4.2.27 QM_pcds
^^^^^^^^^^^^^^

This is a binary file that is required in the **data** directory to
restart a QXMD simulation.


4.2.28 qm_pds.d
^^^^^^^^^^^^^^^

This file gives the different angular momemta contributions to each
energy band. Sample output for a water molecule is given below.

::

   #  Mulliken analysis : s, p, & d contribution to each band     
         0      2     10     -1
     2   0 1
     1   0
       1  0.8239  0.1159  0.0602  1.0000
       2 -0.0000  0.7977  0.2023  1.0000
       3  0.1248  0.8752  0.0000  1.0000
       4  0.0000  1.0000  0.0000  1.0000
       5  0.0130 -0.0058  0.4393  0.4465
       6 -0.0000 -0.0071  0.2595  0.2524
       7 -0.0103  0.0017  0.1011  0.0925
       8  0.0011 -0.0010  0.0401  0.0402
       9  0.0099 -0.0033  0.0676  0.0742
      10  0.0000  0.0306  0.0655  0.0961

After the comment line, the next line gives the step number, the total
number of atomic species in the system ('2' for hydrogen and oxygen),
and the total number energy bands, in this case '10'. The '-1' refers to
\***. The next two lines give the number of diffrent angular momenta
(orbital types) for each atomic species. In this case, the species
corresponding to keyword '1', oxygen, has s-type '0' and p-type '1'
orbitals, while the species corresponding to keyword '2', hydrogen, has
only s-type orbitals. The next 10 lines give the different angular
momemta contributions to each of the energy bands. Looking at the 4th
line, for example, the first number, '1', represents the energy band
index, the next three numbers give the contributions from oxygen s-type
and p-type orbitals and hydrogen s-type orbitals. The final number is a
sum of these three numbers. In this case, the sum is '1' since the first
energy band is comprised entirely of oxygen s-type and p-type orbitals
and hydrogen s-type orbitals. However, higher energy bands may have
some, even majority contribution from higher angular momenta, in which
case the last number in the line will not be '1'.


4.2.28 QM_peig
^^^^^^^^^^^^^^

This is a binary file that is required in the **data** directory to
restart a QXMD simulation.


4.2.29 qm_str.d
^^^^^^^^^^^^^^^

This file outputs components of the stress tensor [GPa] of the system,
computed quantum mechanically, at given time steps intervals (the
interval is defined in IN.PARAM). This file will only be written if
dumping stress data is set to **true** in IN.PARAM. Below is sample
output.

::

   #  Stress in [GPa]               
   #        Pxx             Pyy             Pzz            Pyz             Pzx             Pxy          
        0  1.33128802E+01  1.31392351E+01  5.29606443E+00 -3.16079663E-03  5.31442065E-04 -9.35079164E-01
        5  1.36386034E+01  1.33425150E+01  6.22908157E+00 -7.05568167E-03  1.18151950E-01 -1.04635711E+00
       10  1.51764215E+01  1.56807679E+01  8.35762389E+00  2.34764126E-02  2.40917207E-01 -8.65032053E-01
       15  1.68410502E+01  1.59107220E+01  1.10723229E+01  1.08462786E-01  4.08898975E-01  7.83659509E-01
       20  1.78493502E+01  1.62776114E+01  1.32497640E+01  1.19138247E-01  4.21580509E-01  1.08288678E+00

After the first two comment lines, the following lines list the step
number followed by the components of the stress tensor for that time
step. Here, stress data is dumped every 5 time steps. Note that qm_str.d
varies slightly from md_str.d due since md_str.d includes contributions
from the kinetic energy of the ions, while qm_str.d does not.


4.2.30 QM_tddftfssh
^^^^^^^^^^^^^^^^^^^

This is a binary file that is required in the **data** directory to
restart a Non-Adiabtic QXMD simulation, containing electronic transition
probability data for the system at the last step of the previous run.


4.2.31 qm_td_eig.d
^^^^^^^^^^^^^^^^^^

This file gives the enery eigenvalues and the electron occupation
numbers for each energy band for each time step in a Non-Adiabtic QMD
simulation. As such, this file is only written when TD-DFT is set to
**true**. Sample output for a water molecule with two electrons excited
from the highest occupied molecular orbital (HOMO) to the lowest
unoccupied molecular orbital (LUMO) is shown below.

::

   #  Eigenvalues of GS & occupations of Excited States 
         0     19     10
        1 -2.12398E+00 2.000
        2 -1.11576E+00 2.000
        3 -6.58981E-01 2.000
        4 -5.43674E-01 0.000
        5 -5.55956E-02 2.000
        6  1.62120E-01 0.000
        7  1.63390E-01 0.000
        8  1.87276E-01 0.000
        9  2.13296E-01 0.000
       10  3.57972E-01 0.000

After the comment line, the second line gives the step number, the
cumulative number of SCF iterations performed up to that step number,
and the number of energy bands in the system. The following lines give
the energy band index number, the energy eigenvalue for that band, and
the electronic occupation number for that band. In this case, spin
polarization was not used, however, if it is turned on, there will be
one column for spin-up electrons and one column for spin-down electrons.
As can be seen, the lowest three energy bands are fully occupied, the
fourth band had two electrons removed and placed in the fifth band,
simulation an electronic excitation in a NAQMD simulation.


4.2.32 qm_zan.d
^^^^^^^^^^^^^^^

This file gives energy differences per electron and residuals for each
time step. Note that the word "zansa" means residual in Japanese.

::

   # Difference of tot E/el. E(HF)-E(KS) & maximum & average residuals
   #              difene      difene2     zansa1      zansa2      bfzansa1    bfzansa2
        0    19  1.4700E-07  3.0082E-05  2.0454E-05  2.1813E-08  7.7861E-06  1.7108E-08
        1    32  8.4920E-07  5.7235E-05  5.8779E-05  3.3700E-08  5.3446E-05  1.5236E-08
        2    42  9.8762E-07  2.2899E-04  2.5317E-05  4.9784E-08  3.4752E-06  1.0213E-08
        3    52  8.5924E-07  1.1135E-04  3.0278E-05  2.1025E-08  2.3697E-05  1.2596E-08

After the two comment lines, each line gives the time step number, the
total number of SCF iterations completed up to that time step, the
difference in energy between the current and previous iteration upon
convergence, the difference between the Harris-Foulkes and Kohn-Sham
energies, the maximum (zansa1) and average (zansa2) residuals before the
the Kohn-Sham equations are appoximately solved, and the maximum
(bfzansa1) and average (bfzansa2) residuals after the the Kohn-Sham
equations are appoximately solved.

