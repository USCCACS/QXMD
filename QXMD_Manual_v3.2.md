# QXMD : *ab initio* Molecular Dynamics Simulation Package


## <a name="ht">Table of Contents </a>

- [QXMD : *ab initio* Molecular Dynamics Simulation Package](#qxmd--ab-initio-molecular-dynamics-simulation-package)
  - [<a name="ht">Table of Contents </a>](#a-name%22ht%22table-of-contents-a)
  - [<a name="h0">Introduction </a>](#a-name%22h0%22introduction-a)
  - [<a name="h1"> 1. Prerequisites </a>](#a-name%22h1%22-1-prerequisites-a)
  - [<a name="h2"> 2. Installation </a>](#a-name%22h2%22-2-installation-a)
    - [2.1 Download QXMD](#21-download-qxmd)
    - [2.2 Set Working Directory](#22-set-working-directory)
    - [2.3 Configure Makefiles](#23-configure-makefiles)
    - [2.4 Precompiler variables](#24-precompiler-variables)
    - [2.5 Adding New Machine configuration in Makefile](#25-adding-new-machine-configuration-in-makefile)
    - [2.6 Build QXMD](#26-build-qxmd)
  - [<a name="h3"> 3. Program Structure </a>](#a-name%22h3%22-3-program-structure-a)
    - [3.1 Control Directory](#31-control-directory)
    - [3.2 Data directory](#32-data-directory)
    - [3.3 QXMD executable](#33-qxmd-executable)
  - [<a name="h4"> 4. Files used by **QXMD** </a>](#a-name%22h4%22-4-files-used-by-qxmd-a)
    - [4.1 Input Files](#41-input-files)
    - [4.2 Output Files](#42-output-files)
  - [<a name="h5"> 5. Main input file: IN.PARAM </a>](#a-name%22h5%22-5-main-input-file-inparam-a)
    - [5.1 Parallel Section: *parallel](#51-parallel-section-parallel)
    - [5.2 Start(on/off) Section](#52-startonoff-section)
    - [5.3 TDDFT-MD](#53-tddft-md)
    - [5.4 linear-response TDDFT](#54-linear-response-tddft)
    - [5.5 PAW](#55-paw)
    - [5.6 Cluster](#56-cluster)
    - [5.7 Approximation for Exchange](#57-approximation-for-exchange)
    - [5.8 scissors correction at donor/acceptor interface](#58-scissors-correction-at-donoracceptor-interface)
    - [5.9 k-points](#59-k-points)
    - [5.10 electric field](#510-electric-field)
    - [5.11 spin polarization](#511-spin-polarization)
    - [5.12 SCF iterations](#512-scf-iterations)
    - [5.13 well potenial](#513-well-potenial)
    - [5.14 Kohn-Sham equation](#514-kohn-sham-equation)
    - [5.15 Poisson equation](#515-poisson-equation)
    - [5.16 molecular dynamics](#516-molecular-dynamics)
    - [5.17 save data](#517-save-data)
    - [5.18 soft walls](#518-soft-walls)
    - [5.19 gravitational field](#519-gravitational-field)
    - [5.20 Mulliken analysis](#520-mulliken-analysis)
    - [5.21 Spherical Harmonics expansion](#521-spherical-harmonics-expansion)
    - [5.22 EDA](#522-eda)
    - [5.23 Wannier function](#523-wannier-function)
    - [5.24 Conductivity](#524-conductivity)
    - [5.25 stress calculation](#525-stress_compute)
    - [5.26 dump charge density](#526-dump-charge-density)
    - [5.27 dump wavefunctions](#527-dump-wavefunctions)
    - [5.28 dump local potential](#528-dump-local-potential)
    - [5.29 supercell](#529-supercell)
    - [5.30 vacuum](#530-vacuum)
    - [5.31 spherical region](#531-spherical-region)
    - [5.32 planewaves](#532-planewaves)
    - [5.33 double-grid method](#533-double-grid-method)
    - [5.34 electronic bands](#534-electronic-bands)
    - [5.35 Symmetry operations](#535-symmetry-operations)
    - [5.36 Atoms](#536-atoms)
    - [5.37 Constraint conditions](#537-constraint-conditions)
    - [5.38 Virtual Molecular Dynamics](#538-virtual-molecular-dynamics)
    - [5.39 Cholesky decomposition](#539-cholesky-decomposition)
    - [5.40 Eigenvalue problem](#540-eigenvalue-problem)
    - [5.41 Work Array](#541-work-array)
    - [5.42 Table Dimension](#542-table-dimension)
  - [<a name="h6"> 6. QXMD Simulation Quickstart </a> ](#a-name%22h7%22-6-qxmd-simulation-quickstart-a)
    - [6.1 Minimal IN.PARAM file](#61-minimal-minimal-inparam-file)
    - [6.2 Adding parameters to IN.PARAM](#62-adding-parameters-to-inparam)
  - [<a name="h7"> 7. Utilities </a>](#a-name%22h7%22-7-utilities-a)
    - [7.1 Creating PDB file from output](#71-creating-pdb-file-from-output)
    - [7.2 Picking last configuration from simulation](#72-picking-last-configuration-from-simulation)
    - [7.3 Create gaussian cube file to visualize wave function](#73-create-gaussian-cube-file-to-visualize-wave-function)
    - [7.4 Kohn Sham eigenvalues](#74-kohn-sham-eigenvalues)
  - [References](#references)

------------------------------------------
 
## <a name="h0">Introduction </a>  

QXMD is a Quantum Molecular Dynamics (QMD) simulation software with various eXtensions.
QMD follows the trajectories of all atoms while computing interatomic forces quantum mechanically
in the framework of density functional theory (DFT)
[<a href="https://journals.aps.org/pr/abstract/10.1103/PhysRev.136.B864">P. Hohenberg & W. Kohn, <i>Phys. Rev.</i> <b>136</b>, B864 (1964)</a>;
 <a href="https://journals.aps.org/pr/abstract/10.1103/PhysRev.140.A1133">W. Kohn & L. J. Sham, <i>Phys. Rev.</i> <b>140</b>, A1133 (1965)</a>].
The QXMD software has been developed by Fuyuki Shimojo since 1994 [[1]](#ref1). Since 1999, various extensions have been developed in collaboration with
Rajiv Kalia, Aiichiro Nakano and Priya Vashishta [[2]](#ref2).
<br><br>
The basic QXMD code is based on a plane-wave basis to represent electronic wave functions and pseudopotential (PP) methods to describe electron-ion interaction.
Supported PPs include norm-conserving PP
[<a href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.43.1993">N. Troullier & J. L. Martins, <i>Phys. Rev. B</i> <b>41</b>, 1993 (1991)</a>]
and ultrasoft PP
[<a href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.41.7892">D. Vanderbilt, <i>Phys. Rev. B</i> <b>41</b>, 7892 (1991)</a>],
as well as an all-electron projector augmented-wave (PAW) method
[<a href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.50.17953">P. E. Blochl, <i>Phys Rev B</i> <b>50</b>, 17953 (1994)</a>].
Electron-electron interaction beyond the mean-field Hartree approximation is included using various exchange-correlation functionals,
with and without spin polarization:
generalized gradient approximation (GGA)
[<a href="https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.77.3865">J. P. Perdew, K. Burke & M. Ernzerhof, <i>Phys. Rev. Lett.</i> <b>77</b>, 3865 (1996)</a>],
DFT+U method for transition metals
[<a href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.52.R5467">A. I. Liechtenstein, V. I. Anisimov & J. Zaanen, <i>Phys. Rev. B</i> <b>52</b>, R5467 (1995)</a>],
van der Waals (vDW) functional for molecular crystals and layered materials
[<a href="http://onlinelibrary.wiley.com/doi/10.1002/jcc.20078/abstract">S. Grimme, <i>J. Comput. Chem.</i> <b>25</b>, 1463 (2004)</a>],
nonlocal correlation functional
[<a href="https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.92.246401">M. Dion et al., <i>Phys. Rev. Lett.</i> <b>92</b>, 246401 (2004)</a>],
and range-separated exact-exchange functional
[<a href="http://aip.scitation.org/doi/10.1063/1.1564060">J. Heyd, G. E. Scuseria & M. Ernzerhof, <i>J. Chem. Phys.</i> <b>118</b>, 8207 (2003)</a>].
<br><br>
Various unique capabilities included in the QXMD code (some of which are described in [[2]](#ref2)) include:
<ul>

1. Linear-scaling DFT algorithms [[3]](#ref3),[[4]](#ref4),[[5]](#ref5)
2. Scalable algorithms on massively parallel computers [[6]](#ref6),[[7]](#ref7)
3. Nonadiabatic quantum molecular dynamics (NAQMD) to describe excitation dynamics [[8]](#ref8)
4. Omni-directional multiscale shock technique (OD-MSST) to study shock response of materials [[9]](#ref9)
</ul>
<br>


------------------------------------------

Go to **[Table of Contents](#ht)**, **[Introduction](#h0)**,  **[1.Prerequisites](#h1)**,  **[2.Installation](#h2)**,  **[3.Program Structure](#h3)**, **[4.I/O Files](#h4)** ,**[5.Main input file](#h5)**,  **[6.Utility files](#h6)**

------------------------------------------

## <a name="h1"> 1. Prerequisites </a> 

**QXMD** parallel version requires Fortran compiler, FFT and MPI libraries. There are no other library requirement for the **QXMD**

------------------------------------------

Go to **[Table of Contents](#ht)**, **[Introduction](#h0)**,  **[1.Prerequisites](#h1)**,  **[2.Installation](#h2)**,  **[3.Program Structure](#h3)**, **[4.I/O Files](#h4)** ,**[5.Main input file](#h5)**,  **[6.Utility files](#h6)**

------------------------------------------

## <a name="h2"> 2. Installation </a>

### 2.1 Download QXMD
To get started,  clone this repository to your computer.
```
~$ git clone https://github.com/USCCACS/QXMD_Course.git
```

### 2.2 Set Working Directory

First, change the working directory to **qxmd/**
```
~$ cd qxmd
```
You will see following files and directories when you use the 'ls' command:

```
qxmd $ ls
Include/    Makefile    Qxmd/    Sources/     util/
```

**Include/** contains necessary libraries for compilation, **Makefile** will be used to build the program, **Qxmd** contains default variables for many input parameters, **Sources/** contains all **QXMD** source codes,and **util/** contains helpful codes for post-processing output data from **QXMD**.

### 2.3 Configure Makefiles
You have to choose makefile based on the configuration of your machine. For most of the facility supercomputers makefile is already configured. Complete list of the preconfigured makefile can be obtain by typing **make help**.

````
$make help 

make dec        : DEC Alpha Degital Fortran
make gnu        : GNU Fortran (only for classical MD)
make ffc        : PC LINUX (Fujitsu Fortran&C)
make ffc-v64    : PC LINUX (Fujitsu Fortran&C) on VT-64
make ifc        : PC LINUX (Intel Fortran)
make ifc-x86    : PC LINUX (Intel Fortran) on x86-64
make ifc-v64    : PC LINUX (Intel Fortran) on VT-64
make hpc        : PC LINUX (Intel Fortran) on x86-64 at USC-HPC
make hp         : HP FORTRAN 90
make sp         : IBM SP
make origin     : SGI ORIGIN 2000
make origin26   : SGI ORIGIN 2600
make t3e        : CRAY T3E
make sr2201     : HITACHI SR2201
make sr8000     : HITACHI SR8000
make sr11000    : HITACHI SR11000 at ISSP
make ha8000     : HITACHI HA8000-tc/HT210
make vpp        : Fujitsu VPP5000
make sx7        : NEC SX-7
make altix      : SGI Altix3700/1280 at ISSP
make p5         : IBM eServer model p5
make primequest : FUJITSU PRIMEQUEST 580
make primergy   : FUJITSU PRIMERGY RX200S3
make cx400      : FUJITSU PRIMERGY CX400
make primehpc   : FUJITSU PRIMEHPC FX10
make sekirei    : Intel Fortran on SGI ICE XA/UV (systemB) at ISSP
````
If you are not using aforementioned machine, you might need to modify the **Makefile** according to your computing environment.

### 2.4 Precompiler variables
Precompiler flags are defined in CPPDEFS as:

````
CPPDEFS= -DFLAG1 -DFLAG2
````
Following flags are defined in **QXMD**
````
-LIBFFTW = link FFTW library 
-LIBFFTW3 = link FFTW3 library 
-DPOINTER64 = used for 64-bit machines 
-DSSL2VP = used on Fujitsu machine 
-DVECTOR = special implementation. It may result in faster computation on some machine. 
````
There are several flags are defined for future development. Currently, they are not used. Several flags are specific to compiler. Please check compiler manual for further details.

### 2.5 Adding New Machine configuration in Makefile
**Makefile** defines which compiler you will use to build the **QXMD** executable.There are several predefined compiler setting which must be set. First, open Makefile available in QXMD directory and define the machine name. For example, theta machine at Argonne National Lab

````
theta: $(WDIR) $(DDIR) $(SOURCEFILE)
        sed "s/^#THETA#//" $(SOURCEFILE) > $(WORKFILE)
````
This process will ensure that make will understand the target. Next, we need to define the macro such as **LINKER**, **LDFLAGS** and **FFTLIB** in the makefile included in Source directory. For example, theta machine definition are as follow 

````
#THETA#     LINKER        = ftn
#THETA#     FFLAGS        = -c -fast -warn none                                                                     
#THETA#     LDFLAGS       = -O3 -ipo -no-prec-div 
#THETA#     FFLAGS        = -c -O3 -ipo -no-prec-div  -warn none                                              
#THETA#     CPPDEFS       = -traditional $(CPPDEF) -DLIBFFTW3 -DIFC -DVECTOR                                         
#THETA#     MKLPATH       = $(MKLROOT)/intel/lib64
#THETA#     FFTLIB        = -L$(MKLPATH) -Wl,--start-group -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -Wl,--end-group

````
Definition of each macro is as follow:

```
LINKER= MPI Fortran compiler 
FFLAGS= Compiler flags
LDFLAGS= linking libraries
CPPDEFS= Preprocessing MACROS
MKLPATH= Variable defined for MKL libraries
FFTLIB= FFT library 
```
It is up to the user to define a new macro for better readability. Once you have defined all the macros you are ready to build code for new architecture 

### 2.6 Build QXMD
After tailoring the Makefile for your computing environment, type the command below to build the **QXMD** executable.

```
qxmd $ make $(MACHINE_NAME)
qxmd $ make qxmd
```
To create a **parallel** version of QXMD, please use the following steps 

```
qxmd $ make $(MACHINE_NAME)
qxmd $ make qxmd_mpi
```

Check to see if you the **QXMD** executable has been created, and if so, you are ready to start a simulation.

```
qxmd $ ls
Include/    Makefile    Qxmd/    qxmd_mpi Sources/     util/
```

Parallel build is activated by default. The program will use all available core on the machine to build the program efficiently. 


------------------------------------------

Go to **[Table of Contents](#ht)**, **[Introduction](#h0)**,  **[1.Prerequisites](#h1)**,  **[2.Installation](#h2)**,  **[3.Program Structure](#h3)**, **[4.I/O Files](#h4)** ,**[5.Main input file](#h5)**,  **[6.Utility files](#h6)**

------------------------------------------
## <a name="h3"> 3. Program Structure </a>
Below is a tree diagram of the general structure for the directories and files needed to run a **QXMD** simulation.
```
.
├── control
│   ├── filename
│   ├── IN.CONFIG
│   ├── IN.PARAM
│   ├── IN.VELOC
│   └── PAW
│       ├── H.PBE
│       └── O.PBE
├── data
└── qxmd_mpi
```

An overview is given below, while a more detailed explanation can be found in section 4.1.
### 3.1 Control Directory
The **control** directory is where QXMD will search for all input data.  The main input file, in this case named **IN.PARAM**, contains the settings for all input parameters for running QXMD.  It is best practice to name this file 'IN.PARAM', though it can, in principle, be named anything.  Thus, it is required to have a file named 'filename' which is a simple, one-line text file holding the path (relative to the QXMD executable) to the main input file, in this case 'control/IN.PARAM'.  Most QXMD runs will also require a file detailing the initial configuration of the simulation system (this may be omitted if you are restarting a QXMD run from a previous run).  It is best practice to name the input configuration file 'IN.CONFIG', as shown above in the tree diagram.  Optionally, a file defining the inital velocities of all the atoms may be defined. It is best practice to name the input velocity file 'IN.VELOC', as shown above in the tree diagram.  If not intial velocities are provided, random initial velocities will be used corresponding to the intial temperature of the system.  Finally, a directory holding pseudopotential information for all atomic species present in your simulation system is required.  The **PAW** directory, shown above, holds pseudopotential information for hydrogen in **H.PBE** and oxygen in **O.PBE**, which would both be required for simulating, say, a water molecule.  

### 3.2 Data directory
A directory named **data** is required for any QXMD run.  All output data will be dumped to this directory.  If there are files already present in this directory at the start of a new QXMD run, they will be overwritten.  Thus, it is recommended to move files in **data** from a previous run to a new directory before beginning a new QXMD run.  If you are restarting a QXMD run, there are some required files that must be present in **data**, which will be detailed in section 4.2.

### 3.3 QXMD executable
The QXMD executable, in this case **qxmd_mpi** should be at the same directory level as the **control** and **data** directories.  

------------------------------------------

Go to **[Table of Contents](#ht)**, **[Introduction](#h0)**,  **[1.Prerequisites](#h1)**,  **[2.Installation](#h2)**,  **[3.Program Structure](#h3)**, **[4.I/O Files](#h4)** ,**[5.Main input file](#h5)**,  **[6.Utility files](#h6)**

------------------------------------------

## <a name="h4"> 4. Files used by **QXMD** </a>
### 4.1 Input Files
#### 4.1.1 filename
**filename** is a simple, single line text file containing the path to the main input file, relative to the QXMD executable, in between single quotes. Below shows the entire contents of a sample **filename**

```
'control/IN.PARAM'
```

#### 4.1.2 IN.CONFIG
An initial configuration file must be included for any new QXMD run.  It is best practice to name this file **IN.CONFIG**.  Below shows a sample configuration file for a water molecule.

```
3
1 0.579029 0.305583 0.264540
1 0.420969 0.305577 0.264541
2 0.500002 0.244625 0.264685
```

The first line should contain the total number of atoms in the simulation system.  In this case, the water molecule has 3 atoms: 2 hydrogen and 1 oxygen.  Each of the following lines (one line for each atom) contains a keyword representing the species the atom, as well as the atom's **x**, **y**, and **z** coordinates.  These may be given in real space of fractional coordinates (you define which kind are provided in the IN.PARAM file).  In the example above, the second and third lines give the spatial coordinates for the hydrogen atoms, using '1' as the keyword for hydrogen,  while the last line give the spatial coordinates for the oxygen atom, using '2' as the keyword for oxygen.  Integer numbers should be used for keywords, though it is arbitrary which numbers are chosen for each species, as long as what you choose is used consistently with atomic data entered in the main input file, IN.PARAM.  

#### 4.1.3 IN.PARAM
A main input file is required for all QXMD runs, which contains all the parameters settings. It is best practice to name the main input file **IN.PARAM**.   There are on the order of a few hundred input parameters that may be tailored for various kinds of QXMD simulations, though  most parameters have default settings, many of which need not be changed for most routine QXMD simulations.  Details of these parameters may be found in section 5.

#### 4.1.4 IN.VELOC
A file defining the intial velocities of the atoms may optionally be provided.  It is best practice to name this file **IN.VELOC**.  If this is not provided, random inital velocities will be assigned to each atom, in accordance with the initial system temperature.  Below shows a sample initial velocity file for a water molecule.

```
3
  2.7271510E-03
 1  0.058157  0.043210  0.003117
 2 -1.000000 -0.688863 -0.057829
 2  0.076852  0.002967  0.008360
```

The first line should contain the total number of atoms in the simulation system.  In this case, the water molecule has 3 atoms: 2 hydrogen and 1 oxygen.  The next line is a scaling factor by which all of the following veloctiy components should be multiplied.  Each of the following lines (one line for each atom) contains a keyword representing the species the atom, as well as the atom's **x**, **y**, and **z** velocity components in atomic units. The integer keywords for the atomic species should be consistent with those provided in the IN.CONFIG and IN.PARAM files.

### 4.2 Output Files
All output data files are dumped to the **data** directory during QXMD runs. Output files with names beginning with a lower case letter are human-readable text files detailing simulation data, while output files beginning with upper case letters are binary files to be used for restarting a QXMD run.  These files must be present in the **data** directory when attempting to restart a job.  After a simulation, all output files that you wish to save should be moved out of the **data** directory and store elsewhere as a new QXMD run will write over all files present in **data**. 

#### 4.2.1 md_box.d
This file gives the simulation box size for a classical MD simulation.  A sample of output file md_box.d is shown below.

```
# MD cell vectors
#        L_1           L_2           L_3          angle(2-3) angle(3-1) angle(1-2)
      0 3.7280514E+01 3.7280514E+01 3.0613561E+01  90.000000  90.000000 120.000000
```
L\_1, L\_2, and L_3 are the lengths of three cell vectors in Bohr units, while angles between these vectors, in degrees, follow.  The '0' at the beginning of the last line represents the MD step number.  In this example, the simulation box size is held fixed, and thus this is the only line that appears in the output.  If these numbers were updated at any point, a new line would appear under the last line with the step number in which the box size was altered, followed by the new box dimensions.

Note that there is also an output file named **qm\_box.d**, which provides the supercell size for a quantum MD simulation.  

#### 4.2.2 md_cel.d
Similar to **md\_box.d**, **md\_cel.d** holds the cell vector data for the simulation cell in which molecular dynamics take place, however instead of giving the lengths and angles of the vectors, this file gives the **x**, **y**, and **z** components of the three cell vectors. A sample of output file md_cel.d is shown below.

```
# MD cell vectors
#        L_(1:3,1:3)
      0  3.7280514E+01  0.0000000E+00  0.0000000E+00 -1.8640257E+01  3.2285872E+01  0.0000000E+00  1.8745400E-15  3.2467985E-15  3.0613561E+01
```
After the two comment lines, the file lists the step number, in this case '0', and then **x**, **y**, and **z** components of the first cell vector, followed by those of the second and third cell vectors.  In this case, there were no changes in the cell vectors during the simulation, however, if changes do occur, additional lines will be written beginning with the step number in which the cell vectors were changed, followed by their new values.

#### 4.2.3 md_eng.d
This file contains the Hamiltonian, potential, and kinetic energies [Hartree] of the total system, as well as the temperature [K] at every time step.  A sample of the output is shown below.

```
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
```

#### 4.2.4 md_log
This file provides of log of the simulation, including values of input parameters set and lengths of computation time for force calculations for each time step.

#### 4.2.5 MD_mts0
This is a binary file that is required for restarting a QXMD simulation.

#### 4.2.6 md_spc.d
This file lists the species keyword for each atom at each time step.  A sample of the output is shown below.

```
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
```
In the second line, '2' represents the total number of species in the system. The following numbers,'42' and '34', are the atomic numbers of the two atomic species in the system ordered by the keywords representing those species.  In this example, molybdenum has atomic number 42, and is represented by keyword '1',  while selenium has atomic number 34, and is represented by keyword '2'. The next line lists the step number, followed by the total number of atoms in the system, in this case '12'.

#### 4.2.7 md_str.d
This file outputs components of the stress tensor [GPa] of the system, computed classically, at given time steps intervals (the interval is defined in IN.PARAM).  This file will only be written if dumping stress data is set to **true** in IN.PARAM. Below is sample output.
```
# Stress in [GPa]               
#        Pxx             Pyy             Pzz            Pyz             Pzx             Pxy     
0  1.33846254E+01  1.32248154E+01  5.38055547E+00  2.74846400E-02  4.09607470E-02 -9.31727041E-01
5  1.37099707E+01  1.34067898E+01  6.81653796E+00  1.21177457E-02  2.32549066E-01 -1.03833531E+00
10  1.52252127E+01  1.57211516E+01  9.64813968E+00  3.34044020E-02  3.45617372E-01 -8.50283630E-01
15  1.68630437E+01  1.59302869E+01  1.19623442E+01  1.09971448E-01  4.41730570E-01  7.93058521E-01
20  1.78613247E+01  1.62869327E+01  1.36732341E+01  1.16372513E-01  4.25461913E-01  1.08703422E+00
```
After the first two comment lines, the following lines list the step number followed by the components of the stress tensor for that time step.  Here, stress data is dumped every 5 time steps.  Note that contributions from the kinetic energy of the ions are included in these values, while in qm\_str.d they are not.

#### 4.2.8 md\_str_diag.d
This file outputs the diagonalized stress tensor, along with the respective eigenvectors at given time step intervals.  Below is sample output.

```
# Diagonalized stress in [GPa]  
#        Pxx             Pyy             Pzz            Eigenvectors             
      0  1.42398826E+01  1.23699022E+01  5.38021140E+00  7.3672E-01 -6.7620E-01  1.3084E-03  6.7618E-01  7.3671E-01  6.8597E-03 -5.6024E-03 -4.1690E-03  9.9998E-01
      5  1.46113490E+01  1.25135917E+01  6.80835783E+00  7.5737E-01 -6.5263E-01  2.1581E-02  6.5207E-01  7.5764E-01  2.8228E-02 -3.4773E-02 -7.3072E-03  9.9937E-01
     10  1.46051510E+01  1.63637815E+01  9.62557143E+00  7.9488E-01  6.0384E-01  5.9491E-02 -6.0340E-01  7.9698E-01 -2.7090E-02 -6.3770E-02 -1.4364E-02  9.9786E-01
     15  1.73523452E+01  1.54808731E+01  1.19224564E+01  8.6781E-01  4.9024E-01  8.1122E-02 -4.8915E-01  8.7153E-01 -3.4170E-02 -8.7452E-02 -1.0028E-02  9.9612E-01
     20  1.84552850E+01  1.57357781E+01  1.36304284E+01  8.8838E-01  4.5019E-01  8.9996E-02 -4.4823E-01  8.9293E-01 -4.2081E-02 -9.9304E-02 -2.9550E-03  9.9505E-01
```
After the first two comment lines, the following lines list the step number followed by the components of the diagonalized stress tensor and the **x**, **y**, and **z** components of the three eigenvectors for that time step.  Here, diagonalized stress data is dumped every 5 time steps.

#### 4.2.9 MD_Tocontinue
This is a binary file that is required in the **data** directory to restart a QXMD simulation.

#### 4.2.10 md_vel.d
This file contains scaled atomic velocities for each atom at each time step.  Below is sample output.
```
#  Atomic scaled velocities     
      0     12
 1.7316159E-06
-6.94686 -2.92359 -1.82411  5.50162 -1.28242  2.31888  0.35217 -0.34956  1.04083
-2.29616  2.29092  0.53871 -0.03100 -5.38400 -1.78311  2.52055  1.50530  1.29961
 1.94346  1.81671 -2.48381  5.25972  2.70791  1.54716  5.36706  8.68984  1.05675
-5.95609 -9.99900 -0.80332 -3.99628  1.60613 -0.32696 -0.98935  1.80875 -1.02671
```
The second line gives the step number, in this case '0', and the total number of atoms in the system, in this case, '12'.  The next line gives the number by which to scale the components of velocity given in the following lines (i.e. multiply the following numbers to get the absolute velocity values).  The following lines give the **x**, **y**, and **z** components of velocity for the first atom, followed by those for the second atom, ending with those for the last atom.

#### 4.2.11 qm_box.d
This file holds the cell vector data for the supercell for a quantum molecular dynamics simulation.  A sample of output file is shown below.

```
#  supercell (FFT cell) vectors (lengths & angles)
#        L_1           L_2           L_3          angle(2-3) angle(3-1) angle(1-2)                   
      0 1.2426838E+01 1.2426838E+01 3.1114338E+01  90.000000  90.000000 120.000000
```
L\_1, L\_2, and L_3 are the lengths of three supercell vectors in Bohr units, while angles between these vectors, in degrees, follow.  The '0' at the beginning of the last line represents the step number.  In this example, the simulation box size is held fixed, and thus this is the only line that appears in the output.  If the supercell size were changed, a new line would appear under the last line with the step number in which the size was altered, followed by the new supercell dimensions.

Note that there is also an output file named **md\_box.d**, which holds the box size for a classical MD simulation.  For purely quantum MD simulations, **md\_box.d** and **qm\_box.d** will hold the same information (except in cases where the double-grid method is used).  However, it is possible to perform a hybrid classical-quantum MD simulation in which case the system is treated classically, except for a defined area inside which is treted quantum mechanically.  In this case, **md\_box.d** gives the box size of the entire system, while **qm\_box.d** defines the supercell area that is to be treated with QM.

#### 4.2.12 QM_cds
This is a binary file that is required in the **data** directory to restart a QXMD simulation, containing charge density data for the system at the last step of the previous run.

#### 4.2.13 qm_cel.d
Similar to **qm\_box.d**, **qm\_cel.d** holds the cell vector data for the supercell, however instead of giving the lengths and angles of the vectors, this file gives the **x**, **y**, and **z** components of the three supercell vectors. Below is sample output.

```
#  supercell (FFT cell) vectors
#        L_(1:3,1:3)
      0  1.2426838E+01  0.0000000E+00  0.0000000E+00 -6.2134191E+00  1.0761957E+01  0.0000000E+00  1.9052037E-15  3.2999097E-15  3.1114338E+01
```
After the two comment lines, the file lists the step number, in this case '0', and then **x**, **y**, and **z** components of the first supercell vector, followed by those of the second and third supercell vectors.  In this case, there were no changes in the supercell vectors during the simulation, however, if changes do occur, additional lines will be written beginning with the step number in which the supercell vectors were changed, followed by their new values.

#### 4.2.14 QM_cell
This is a binary file that is required in the **data** directory to restart a QXMD simulation, containing simulation cell data for the system at the last step of the previous run.

#### 4.2.15 QM_eig
This is a binary file that is required in the **data** directory to restart a QXMD simulation, containing energy eigenvalue data for the system at the last step of the previous run.

#### 4.2.16 qm_eig.d
This file lists the energy eigenvalues, along with their electronic occupation number for each time step. Below is sample output.
```
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
```
The first number in the second line gives the step number, in this case '0'.  The second number gives the cumulative number of SCF iterations performed, in this case '7' iterations were performed in step 0. Finally, the third number, '10' represents the total number of energy eigenvalues and occupation numbers to follow.  This number corresponds to the number of energy bands for the system, defined in IN.PARAM.  Thus the next ten lines list the band index number, the energy in eV, and the number of electrons which occupy that energy band.  In the second step, 4 SCF iterations were performed since the number '11' follows the step number '1' (7 SCF iterations were performed in the first step, followed by 4 SCF iterations, for a total of 11 SCF iterations after the first two time steps.) In this example, the eight electrons of the system occupy the four lowest energy bands in the first time two steps.

#### 4.2.17 qm_eng.d
This file provides various components of the system's energy [Rydberg] at each time step.  Below is sample output.
```
#  Total potential energy and energy parts in [Ryd.] units
#              Total(HF)         Total(KS)         Kinetic      External      Hartree       Exchange      Correlation   ------        Entropy       Onsite E.     --------      ------        Ewald E.      DFT-D
     0     7   -4.4003610122E+01 -4.4003640357E+01  1.370805E+01 -6.080189E+01  2.800936E+01 -7.786772E+00 -5.786949E-01  0.000000E+00  0.000000E+00 -1.533928E+01  0.000000E+00  0.000000E+00 -1.214255E+00 -1.610012E-04
     1     11   -4.4002961457E+01 -4.4003371148E+01  1.367946E+01 -6.058249E+01  2.791215E+01 -7.771542E+00 -5.778955E-01  0.000000E+00  0.000000E+00 -1.532657E+01  0.000000E+00  0.000000E+00 -1.336327E+00 -1.631839E-04
     2    15   -4.4001366424E+01 -4.4002003283E+01  1.365563E+01 -6.038887E+01  2.782594E+01 -7.758172E+00 -5.771782E-01  0.000000E+00  0.000000E+00 -1.531663E+01  0.000000E+00  0.000000E+00 -1.442551E+00 -1.647241E-04
     3    18   -4.3999479244E+01 -4.3999841762E+01  1.363668E+01 -6.023643E+01  2.775745E+01 -7.747680E+00 -5.766015E-01  0.000000E+00  0.000000E+00 -1.530791E+01  0.000000E+00  0.000000E+00 -1.525182E+00 -1.650145E-04
     4    22   -4.3997932717E+01 -4.3997243098E+01  1.362538E+01 -6.013514E+01  2.771149E+01 -7.740722E+00 -5.762123E-01  0.000000E+00  0.000000E+00 -1.530223E+01  0.000000E+00  0.000000E+00 -1.579643E+00 -1.636864E-04
     5    25   -4.3997129513E+01 -4.3997300982E+01  1.362022E+01 -6.009034E+01  2.769035E+01 -7.737485E+00 -5.760705E-01  0.000000E+00  0.000000E+00 -1.530042E+01  0.000000E+00  0.000000E+00 -1.603399E+00 -1.606360E-04
```
After the first two comment lines, each line gives the step number, the cumulative number of SCF iterations completed up to that step number, followed by various components of system energy as labeled by the column titles.

#### 4.2.18 qm_fer.d
This file gives the Fermi energy of the system at each time step.  Below is sample output.
```
#  Fermi energy
      0      7 -2.90493E-01
      1     11 -2.90662E-01
      2     15 -2.91483E-01
      3     18 -2.92223E-01
      4     22 -2.92859E-01
      5     25 -2.93001E-01
```
The first column gives the step number, the second column gives the cumulative number of SCF iterations up to that step number, while the third column gives the Fermi energy in eV.

#### 4.2.19 qm_frc.d
This file gives the three components of force on each atom at each time step.  Below is sample output for a monolayer of MoSe2 with 12 total atoms.
```
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
```
After the comment line, the next line gives the step number, the total number of atomic species in the system ('2' for molybdenum and selenium), the number of atoms of the species corresponding to keyword '1' (in this case, keyword '1' was used for molybdenum and there are '4' molybdenum atoms), followed by the number of atoms of the species corresponding to keyword '2' (in this case, keyword '2' was used for selenium and there are '8' selenium atoms).  The number in the third line is a scaling factor for the following force vector components (multiply all of the following force components by this number to get the true values for force).  The following lines give the **x**, **y**, and **z** components of force on the first atom, followed by those for the second atom, and so on.  

#### 4.2.20 qm\_fsshprob\_\*\*\*to\*\*\*-u.d
This set of files are only written during a Non-Adiabatic QMD simulation (i.e. TD-DFT set to .true.).  These files give the probabilities for electrons to hop from one band to another, where these band indices will be given in the title of the file in place of the '***'. A sample file name from this set could be qm\_fsshprob\_29to32-u.d, which will give the probability for an electron to transition from band index 29 to band index 32 at each time step that has a non-zero probability.  The '-u' at the end of the file name indicates that the data is for spin-up electrons in the case that spin polarization is used in the QXMD simulation.  In this case, spin-down electron data will be stored in files named 'qm\_fsshprob\_\*\*\*to\*\*\*-d.d' If spin polarization is not used, all electron transition probabilities will be stored in the spin-up data files by default.  Below is a sample of this output from this file.
```
# step  probability  accumulation
    14  8.19528E-07  8.19528E-07
    16  4.54580E-06  4.54580E-06
    17  1.01866E-05  1.47324E-05
    18  5.77833E-06  2.05108E-05
```
As the column titles suggest, the first column gives the step number at which there exists a finite probability, given by the second column, of the electron hopping from, in this case band index 29 to 32. Probabilities are only written for step numbers for which there is a finite (non-zero) probability.  Thus, in the first 13 steps of the simulation there was zero probability for the electron to hop from band 29 to 32, as well as time step 15.  However, at time steps 14, 16, 17, and 18, there were non-zero probabilities for this electronic transition.  The accumulation column refers to ???


#### 4.2.21 QM_hrt
This is a binary file that is required in the **data** directory to restart a QXMD simulation.

#### 4.2.22 QM_ion
This is a binary file that is required in the **data** directory to restart a QXMD simulation, containing atomic position data for the system at the last step of the previous run.

#### 4.2.23 qm_ion.d
This files gives the three components of the positions for all the atoms at each time step.  It can be used for visualization of the atomic trajectories throughout the simulation.  Below is sample output for monolayer MoSe2 with a total of 12 atoms.
```
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
 ```
After the comment line, the next line gives the step number, the total number of atomic species in the system ('2' for molybdenum and selenium), the number of atoms of the species corresponding to keyword '1' (in this case, keyword '1' was used for molybdenum and there are '4' molybdenum atoms), followed by the number of atoms of the species corresponding to keyword '2' (in this case, keyword '2' was used for selenium and there are '8' selenium atoms).  The number in the third line is a scaling factor for the following spatial coordinates of each atom (multiply all of the following spatial coordinates by this number to get the true atomic positions).  The following lines give the **x**, **y**, and **z** coordinates of the first atom, followed by those for the second atom, and so on. 

 #### 4.2.24 qm_log
 This file provides of log of simulation details, including values of input parameters set, time statistics and energies computed for each SCF iteration for every simulation time step, as well as total computation time.  This file is the best place to look to determine why a QXMD simulation fails/crashes, and general debugging.

#### 4.2.25 qm_mul.d
This file provides Mulliken analysis data.  Note that this file will only be written if Mulliken analysis is set to .true. in IN.PARAM.  Sample output for a water molecule is shown below.
```
#  Mulliken analysis : s, p, & d population for each atom      
      0      2      2      1     -1
  2   0 1
  1   0
    1  1.8974  5.5777  7.4751
    2  0.2624  0.2624
    3  0.2625  0.2625

```
After the first comment line, the second line gives the step number, the total number of atomic species in the system ('2' for hydrogen and oxygen), the number of atoms of the species corresponding to keyword '1' (in this case, keyword '1' was used for hydrogen and there are '2' H atoms), followed by the number of atoms of the species corresponding to keyword '2' (in this case, keyword '2' was used for oxygen and there is '1' O atom).  The '-1' at the end of the line may be ignored as is it only included for backward compatibility with older versions of the code. <br>
The next set of two lines gives the number of different angular momenta (orbital types) each species has.  Since oxygen was assigned keyword '1', the first of these lines (line 3) starts with '2' indicating that oxygen has 2 types of orbitals which are identified as '0' (meaning 's-type') and '1' (meaning 'p-type').  The second of these line corresponds to hydrogen, which only has one orbtial, '0' or 's-type'. <br>
Finally, the partial charges in each orbital type, followed by the total partial charge around each atom is given in the last three lines.  For example, in line 5 of the sample output, '1' refers to atom 1, in this case an oxygen atom, the next number is the partial charge in the s-type orbital, the next number is the partial charge in the p-type orbitals, while the last number is the total partial charge around the oxygen atom.  The final two lines give the partial charge in the s-type orbital in atoms '2' and '3', the two H atoms, along with the total partial charge around the H atoms which is the same as the s-type orbital partial charge, since there is only one orbital.


#### 4.2.26 qm_ovp.d
This file provides the overlap charge densities between atoms at the given time step intervals.  Note that this file is only written if Mulliken analysis is set to **true**. Below is sample output for a monolayer of MoSe2 with 12 total atoms.
```
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
```
After the comment line, the next line gives the step number, the total number of atomic species in the system ('2' for molybdenum and selenium), the number of atoms of the species corresponding to keyword '1' (in this case, keyword '1' was used for molybdenum and there are '4' molybdenum atoms), followed by the number of atoms of the species corresponding to keyword '2' (in this case, keyword '2' was used for selenium and there are '8' selenium atoms). The '-1' at the end of the line may be ignored as is it only included for backward compatibility with older versions of the code. <br>
The next set of two lines give the total number of angular momenta (orbital types) for each species, in this case both Mo and Se have s-type ('0'), p-type ('1'), and d-type ('2') orbitals.
The '9' in the next line signifies that the following lines will give partial charge overlaps between atoms number 9 and the other atoms.  Lines 6-14 begin the the atom number, and then give the overlap between the all combinations of orbital types of atom 9 and all orbital types of the line's respective atom number. <br>
Line 6, for example, gives the partial charge overlaps between atom 9 and atom 2.  The 9 charge overlap numbers give overlaps between: s-type(atom 9)/s-type(atom 2), s-type(atom 9)/p-type(atom 2), s-type(atom 9)/d-type(atom 2), p-type(atom 9)/s-type(atom 2), p-type(atom 9)/p-type(atom 2), p-type(atom 9)/d-type(atom 2), d-type(atom 9)/s-type(atom 2), d-type(atom 9)/p-type(atom 2), d-type(atom 9)/d-type(atom 2).


#### 4.2.27 QM_pcds
This is a binary file that is required in the **data** directory to restart a QXMD simulation.

#### 4.2.28 qm_pds.d
This file gives the different angular momemta contributions to each energy band.  Sample output for a water molecule is given below.
```
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
```
After the comment line, the next line gives the step number, the total number of atomic species in the system ('2' for hydrogen and oxygen), and the total number energy bands, in this case '10'.  The '-1' refers to ***.
The next two lines give the number of diffrent angular momenta (orbital types) for each atomic species.  In this case, the species corresponding to keyword '1', oxygen, has s-type '0' and p-type '1' orbitals, while the species corresponding to keyword '2', hydrogen, has only s-type orbitals.  The next 10 lines give the different angular momemta contributions to each of the energy bands.  Looking at the 4th line, for example, the first number, '1', represents the energy band index, the next three numbers give the contributions from oxygen s-type and p-type orbitals and hydrogen s-type orbitals.  The final number is a sum of these three numbers.  In this case, the sum is '1' since the first energy band is comprised entirely of oxygen s-type and p-type orbitals and hydrogen s-type orbitals.  However, higher energy bands may have some, even majority contribution from higher angular momenta, in which case the last number in the line will not be '1'.

#### 4.2.28 QM_peig
This is a binary file that is required in the **data** directory to restart a QXMD simulation.

#### 4.2.29 qm_str.d
This file outputs components of the stress tensor [GPa] of the system, computed quantum mechanically, at given time steps intervals (the interval is defined in IN.PARAM).  This file will only be written if dumping stress data is set to **true** in IN.PARAM. Below is sample output.
```
#  Stress in [GPa]               
#        Pxx             Pyy             Pzz            Pyz             Pzx             Pxy          
     0  1.33128802E+01  1.31392351E+01  5.29606443E+00 -3.16079663E-03  5.31442065E-04 -9.35079164E-01
     5  1.36386034E+01  1.33425150E+01  6.22908157E+00 -7.05568167E-03  1.18151950E-01 -1.04635711E+00
    10  1.51764215E+01  1.56807679E+01  8.35762389E+00  2.34764126E-02  2.40917207E-01 -8.65032053E-01
    15  1.68410502E+01  1.59107220E+01  1.10723229E+01  1.08462786E-01  4.08898975E-01  7.83659509E-01
    20  1.78493502E+01  1.62776114E+01  1.32497640E+01  1.19138247E-01  4.21580509E-01  1.08288678E+00
```
After the first two comment lines, the following lines list the step number followed by the components of the stress tensor for that time step.  Here, stress data is dumped every 5 time steps.  Note that qm\_str.d varies slightly from md\_str.d due since md\_str.d includes contributions from the kinetic energy of the ions, while qm\_str.d does not.

#### 4.2.30 QM_tddftfssh
This is a binary file that is required in the **data** directory to restart a Non-Adiabtic QXMD simulation, containing electronic transition probability data for the system at the last step of the previous run.  

#### 4.2.31 qm\_td\_eig.d
This file gives the enery eigenvalues and the electron occupation numbers for each energy band for each time step in a Non-Adiabtic QMD simulation.  As such, this file is only written when TD-DFT is set to **true**. Sample output for a water molecule with two electrons excited from the highest occupied molecular orbital (HOMO) to the lowest unoccupied molecular orbital (LUMO) is shown below.
```
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
```
After the comment line, the second line gives the step number, the cumulative number of SCF iterations performed up to that step number, and the number of energy bands in the system.
The following lines give the energy band index number, the energy eigenvalue for that band, and the electronic occupation number for that band.  In this case, spin polarization was not used, however, if it is turned on, there will be one column for spin-up electrons and one column for spin-down electrons.  As can be seen, the lowest three energy bands are fully occupied, the fourth band had two electrons removed and placed in the fifth band, simulation an electronic excitation in a NAQMD simulation.

#### 4.2.32 qm_zan.d
This file gives energy differences per electron and residuals for each time step. Note that the word "zansa" means residual in Japanese.
```
# Difference of tot E/el. E(HF)-E(KS) & maximum & average residuals
#              difene      difene2     zansa1      zansa2      bfzansa1    bfzansa2
     0    19  1.4700E-07  3.0082E-05  2.0454E-05  2.1813E-08  7.7861E-06  1.7108E-08
     1    32  8.4920E-07  5.7235E-05  5.8779E-05  3.3700E-08  5.3446E-05  1.5236E-08
     2    42  9.8762E-07  2.2899E-04  2.5317E-05  4.9784E-08  3.4752E-06  1.0213E-08
     3    52  8.5924E-07  1.1135E-04  3.0278E-05  2.1025E-08  2.3697E-05  1.2596E-08
```
After the two comment lines, each line gives the time step number, the total number of SCF iterations completed up to that time step, the difference in energy between the current and previous iteration upon convergence, the difference between the Harris-Foulkes and Kohn-Sham energies, the maximum (zansa1) and average (zansa2) residuals before the the Kohn-Sham equations are appoximately solved, and the maximum (bfzansa1) and average (bfzansa2) residuals after the the Kohn-Sham equations are appoximately solved.


------------------------------------------

Go to **[Table of Contents](#ht)**, **[Introduction](#h0)**,  **[1.Prerequisites](#h1)**,  **[2.Installation](#h2)**,  **[3.Program Structure](#h3)**, **[4.I/O Files](#h4)** ,**[5.Main input file](#h5)**,  **[6.Utility files](#h6)**

------------------------------------------


## <a name="h5"> 5. Main input file: IN.PARAM </a>
Note: The name of the main input file (e.g. 'IN.PARAM') must be specified in 'control/filename'.

The main input file is divided into sections corresponding to different controls. Sections start with '\*SECTION_NAME' and end with '*end'.  To disable any section, simply prepend section title with a '#'. <br>

Example: To disable section SECTION_NAME_2

```
*SECTION_NAME_1
...
*end

#*SECTION_NAME_2
...
*end
```

 Each section is futher divided into subsections, given by a title in parentheses (e.g. (QM-Nodes)), where input control parameters are defined.  To disable any subsection, simply prepend subsection title with a '#'.

 Example: To disable subsection\_2 in section \*SECTION_NAME

```
*SECTION_NAME
(subsection_1)
#(subsection_2)
(subsection_3)
*end
```


Note: Any unnecessary sections or subsections may be removed from input.file, depending on your simulation needs (recommended only for advanced users).

### 5.1 Parallel Section: *parallel
#### 5.1.1 (QM-Nodes)
**npx, npy, npz** = [1|even integer] [1|even integer] [1|even integer]

Default:<br />
**npx** = 1<br>
**npy** = 1<br>
**npz** = 1

QM-Nodes is a space delimited set of three integers **npx**, **npy**, **npz**.  **QXMD** uses a hybrid spatial and band decompisition for parallel computing, and this set of three integers multiplied together defines the number of MPI ranks. For more detail on parallelization, please see [[1]](#shimojo_2014).

#### 5.1.2 (k-points)
**npk** = [integer]

Default:<br />
**npk** = 1

Defines the number of MPI ranks to parallelize k-point sampling.  See section 5.9 for k-point sampling.

#### 5.1.3 (linear-response TDDFT)
**nplr**  = [integer]

Default:<br />
**nplr**  = 1

Defines the number of MPI ranks for parallelization of matrix computations for linear response TDDFT calculations.

#### 5.1.4 (MD-nodes)
**md\_npx, md\_npy, md\_npz** = [1|even integer] [1|even integer] [1|even integer]

MD-Nodes is a space delimited set of three integers: **md_npx, md_npy, md_npz**.  These variables are used only for spatial decomposition in classical MD and divide-and-conqure QMD simulations. For linear scaling DFT please see reference [[1-7]](#shimojo_2014)

#### Reference:
<a name="shimojo_2014">1. "A divide-conquer-recombine algorithmic paradigm for large spatiotemporal quantum molecular dynamics simulations,"
   F. Shimojo, S. Hattori, R. K. Kalia, M. Kunaseth, W. Mou, A. Nakano, K. Nomura, S. Ohmura, P. Rajak, K. Shimamura & P. Vashishta
   <a href="http://aip.scitation.org/doi/abs/10.1063/1.4869342?journalCode=jcp"><i>J. Chem. Phys.</i> <b>140</b>, 18A529 (2014)</a><br></a>
<a name="Shimojo_2001">2. "Linear-scaling density-functional-theory calculations of electronic structure based on real-space grids: design, analysis, and scalability test of parallel algorithms,"
   F. Shimojo, R. K. Kalia, A. Nakano & P. Vashishta
   <a href="https://www.sciencedirect.com/science/article/pii/S0010465501002478"><i>Comput. Phys. Commun.</i> <b>140</b>, 303 (2001)</a><br></a>
<a name="Shimojo_2005">3. "Embedded divide-and-conquer algorithm on hierarchical real-space grids: parallel molecular dynamics simulation based on linear-scaling density functional theory,"
   F. Shimojo, R. K. Kalia, A. Nakano & P. Vashishta,
   <a href="http://www.sciencedirect.com/science/article/pii/S0010465505000688"><i>Comput. Phys. Commun.</i> <b>167</b>, 151 (2005)</a><br></a>
<a name="Shimojo_2008">4. "Divide-and-conquer density functional theory on hierarchical real-space grids: parallel implementation and applications,"
   F. Shimojo, R. K. Kalia, A. Nakano & P. Vashishta,
   <a href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.77.085103"><i>Phys. Rev.</i> <b>77</b>  085103 (2008)</a><br></a>
<a name="Nakano_2001">5. "Scalable atomistic simulation algorithms for materials research,"
   A. Nakano, R. K. Kalia, P. Vashishta, T. J. Campbell, S. Ogata, F. Shimojo & S. Saini
   <a href="http://ieeexplore.ieee.org/document/1592786/"><i>Proc. Supercomputing, SC01</i> (ACM/IEEE, 2001)</a><br></a>
<a name="Nomura_2014">6. "Metascalable quantum molecular dynamics simulations of hydrogen-on-demand,"
   K. Nomura, R. K. Kalia, A. Nakano, P. Vashishta, K. Shimamura, F. Shimojo, M. Kunaseth, P. C. Messina& N. A. Romero
   <a href="http://ieeexplore.ieee.org/document/7013041/"><i>Proc. Supercomputing, SC14</i> (IEEE/ACM, 2014)</a><br></a>
<a name="Shimojo_2013">7. "Large nonadiabatic quantum molecular dynamics simulations on parallel computer,"
   F. Shimojo, S. Ohmura, W. Mou, R. K. Kalia, A. Nakano & P. Vashishta
   <a href="https://www.sciencedirect.com/science/article/pii/S0010465512002548"><i>Comput. Phys. Commun.</i> <b>184</b>, 1 (2013)</a><br></a>


### 5.2 Start(on/off) Section
#### 5.2.1 (On/Off)
**lstart** = [boolean]

Default: <br />
**lstart** = .FALSE.

.FALSE. : Start self-consistent field (SCF) iteration using a random initial wavefunction <br/>
.TRUE. : Continue SCF iteration from a previous run's output wavefunction.  Note: QM_*??? output files must be present in the **data/** directory from the previous calculation.

Determines whether you would like to restart a successfully completed simulation.  Note that if you are restarting a NAQMD simulation, you will also want to set **ltddft_start** to .TRUE. in the *TDDFT-MD section to properly restart a successfully completed simulation.

### 5.3 TDDFT-MD
#### 5.3.1 (On/Off)
**ltddft** = [boolean]

Default: <br />
**ltddft** = .FALSE.

.FALSE. : Execute Adiabtic QMD based on density functional theory (DFT) <br/>
.TRUE. : Execute non-adiabtic QMD (NAQMD) based on time-dependent density functional theory (TDDFT) [[1]](#Gross_1990).

Determines whether to run QMD simulation under adiabatic or non-adiabatic methods.  Adiabatic methods simulate thermodynamic system in ground state equilibrium, while non-adiabatic methods simulate electronic excitations.

#### 5.3.2 (FSSH)
**ltddft\_fssh**  = [boolean]

Default: <br />
**ltddft\_fssh** = .TRUE.

.TRUE. : Perform NAQMD based on Fewest Switches Surface Hopping (FSSH) method [[2]](#Tully_1990)<br/>
.FALSE. : Perform NAQMD based on based on Ehrenfest dynamics (not yet implemented).

Determines the implementation method for electron state dynamics in NAQMD. Fewest Switches surface hopping method is proposed by J. Tully [[2]](#Tully_1990) molecular dynamics simulation of the processes including electronic transition. Next few flags ask for the specification of FSSH method

#### 5.3.3 (FSSH-switch)
**lfssh\_switch**  = [boolean]

Default: <br />
**lfssh\_switch** = .TRUE.

.TRUE. : Allow electrons to move between excited states. <br/>
.FALSE. : Keep electronic occupations fixed.

Determines whether electronic occupations can change throughout the NAQMD simulation.

#### 5.3.4 (FSSH-ground-state-SCF)
**lfssh\_gsscf** = [boolean]

Default: <br>
**lfssh\_gsscf** = .FALSE.

.TRUE. : SCF iterations performed based on the ground  state <br>
.FALSE. : SCF with the excited state

This parameter should always be set to **true** to obtain convergence. 


#### 5.3.5 (FSSH-charge mixing)
**imxchg\_fssh** = 0 | 1 | 2 | 3 <br/>
**aslh\_fssh, bslh\_fssh** = [real] [real]

Default: <br/>
**imxchg\_fssh** = 1


| imxchg\_fssh | Charge Mixing Method | Notes for aslh\_fssh, bslh\_fssh |
| ------------ | :------------------: | :------------------------------: |
| **0**        | No Mixing            |  N/A                             |
| **1**        | Pulay [[3-](#Pulay_1980)[4](#Payne_1992)]               | 0.9, 0.6 (recommended values)    |
| **2**        | Anderson             |                                  |
| **3**        | Simple               |  bslh\_fssh not used             |

**imxchg\_fssh** is used to specify the method for charge mixing during SCF iterations, while **aslh\_fssh** and **bslh\_fssh** are used as tuning parameters.
Note that this section is only considered when **lfssh\_gsscf** is set to .TRUE.

#### 5.3.6 (FSSH-random-initialize)
**lfssh\_random** = [boolean] <br/>
**rseed\_fssh** = [real]

Default: <br/>
**lfssh\_random** = .FALSE.

**lfssh\_random** = .TRUE. : Automatically seed random number generator. <br/>
**lfssh\_random** = .FALSE. : Specify the seed for the random number generator with value given by **rseed\_fssh**.

Determines how to specify the seed for the random number generator used by the FSSH method.

#### 5.3.7 (Boltzmann factor for upward transition)
**lfssh\_boltzmn** = [boolean]

Default: <br/>
**lfssh\_boltzmn** = .TRUE.

.TRUE. : Multiply electronic transition probability by the Boltzmann factor. <br/>
.FALSE. : Leave electronic transtion probability unaltered.

Determines whether or not to multiply the electronic transition probability by the Boltzmann factor when the electronic excitation energy increases due to the transition. This is used in order to approximately satisfy the detailed balance condition.

#### 5.3.8 (velocity scaling)
**lfssh\_vscale** = [boolean] <br/>
**tminimum** = [real]

Default: <br/>
**lfssh\_vscale** = .FALSE.

.TRUE. : Rescale atomic velocities. <br>
.FALSE. : Do not rescae atomic velocities.

Determines whether or not to rescale atomic velocities upon electronic excitation.  **tminimum** gives the minimum temperature in [K] and is used to constrain velocity scaling.

#### 5.3.9 (time step)
**dttddft** = [real]

Default: <br/>
**dttddft** = 0.02d0

Gives the time step in [a.u.] for numerically integrating the TDDFT equations.

#### 5.3.10 (parallel calculation)
**lfssh_parallel** = [boolean]

Default: <br/>
**lfssh_parallel** = .TRUE.

.TRUE. : Solves time-dependent K-S equations in parallel <br>
.FALSE. : Solves time-dependent K-S equations serially.

Determines whether or not to perform TDDFT calculations in parallel.

#### 5.3.11 (restart)
**ltddft_start** = [boolean]

Default: <br/>
**ltddft_start** = .FALSE.

.FALSE. : Initialize electronic occupations as specified in the (occupations) subsection <br>
.TRUE. : Initialize electronic occupations with their values from the last step of a previously completed simulation.  All files beginning with QM\_* and MD\_* must be present in the **data/** directory.

Determines how to initialize electronic excitations for an NAQMD simulation.  Note that if you are continuing a calculation by setting **lstart** to .TRUE. in the *start(on/off), you usually want **ltddft_start** to be set to .TRUE. as well to properly restart a successfully completed simulation.

#### 5.3.12 (initial exciton)
**nexciton** = [integer] <br/>
**iband\_hole, iband\_electron, ldegenerate** = [integer] [integer] [boolean]

Default: <br/>
**nexciton** = 0

**nexciton** = total number of excitons <br/>
**iband\_hole** = band index of the hole. <br/>
**iband\_electron** = band index of the electron. <br/>
**ldegenerate** = .TRUE. : triplet <br/>
**ldegenerate** = .FALSE. : singlet

This set of variables are used to define an exciton in linear-reponse TDDFT.  Thus, these variables are only read when **lrtddft** is set to .TRUE. in the *linear-response TDDFT section.<br/>

**nexciton** gives the total number of excitons to be initially created and should be followed by that many lines of space delimted values for **iband\_hole, iband\_electron, ldegenerate**, where each set specifies one exciton. **nexciton** should be only be set to 0 or 1, as higher numbers of excitons is not gaurenteed to work.

Example: <br/>
```
(initial exciton)
    1
    10  11  .FALSE.
```

#### 5.3.13 (ground state force)
**lfssh_gsfrc** = [boolean]

Default: <br/>
**lfssh_gsfrc** = .FALSE.

.TRUE. : ground state forces are used in FSSH <br/>
.FALSE. : the excited state forces are used. <br/>

Note that this variable is only read if **lfssh_gsscf** is set to .TRUE. in the (FSSH-ground-state-SCF) subsection.

#### 5.3.14 (NSC force)
**ltddft_nscforce** = [boolean]

Default: <br/>
**ltddft_nscforce** = .TRUE.

.TRUE. : On <br>
.FALSE.: Off

Determines whether or not to calculate excited state forces using a non-self-consistent (NSC) method.

#### 5.3.15 (occupations)
**nocc_change** = [integer] <br/>
**num\_band, occ\_new** = [integer] [real] [real]

**nocc_change** = total number of bands with electronic occupation number to be changed. <br/>
**num_band** = band index for occupation change. <br/>
**occ_new** = new occupation number for the corresponding band for up (and down) electrons. <br/>

**nocc_change** gives the total number of bands that will undergo a change in electronic occupation number and should be followed by that many lines of space delimted values for **num\_band, occ\_new**, where each set specifies the band index undergoing a change in its electron occupation values.  If not using spin polarization, only the first number in the **occ\_new** varibale will be read.  If using spin-polarization, you must set **lspin** to .TRUE. in the *spin polarization section, then, you may specify the number of up electrons to move with the first number in the **occ\_new** varibale set and the number of down electrons to move with the second number in the **occ\_new** varibale set.

Example for no spin-polarization:
This example shows how to excited both electrons in band index 10 and both electrons in band index 11, to band indices 12 and 13 (four total electrons changing occupation).  This assumes that in the ground state band indicies 10 and 11 are fully occupied and band indices 12 and 13 are empty.
```
(occupations)
4
10   0.0   0.0
11   0.0   0.0
12   2.0   0.0
13   2.0   0.0
```

Example for spin-polarization:
This example shows how to excite the up electron in band index 10 and the down electron in band index 11 to band indices 12 and 13, respectively (two total electrons changing occupation).  This assumes that in the ground state band indicies 10 and 11 are fully occupied and band indices 12 and 13 are empty.
```
(occupations)
2
10   0.0   1.0
11   1.0   0.0
12   1.0   0.0
13   0.0   1.0
```

#### 5.3.16 (broadening)
**tdbroad** = [real]

Default:
**tdbroad** = 0.0

Determines the width of Gaussian broadening of the Fermi surface in [K].
Note: **tdbroad** = 0.0 denotes no broadening.

#### 5.3.17 (DISH)
**lfssh\_dish** = [boolean] <br/>
**ndishpair** = [integer] <br/>
**ndishi, ndishj, decoherence\_rate** = [integer] [integer] [real]

Default: <br/>
**lfssh\_dish** = .FALSE.
**ndishpair** = 0

**lfssh\_dish** = .TRUE. : Enables Decoherence-Induced Surface Hopping (DISH) <br/>
**lfssh\_dish** = .FALSE. : Disables DISH <br/>
**ndishpair** = the number of state pairs <br/>
**ndishi, ndishj, decoherence\_rate** = the two band indices between which to define the decoherence rate in [a.u.] for DISH

Example:
```
(DISH)
 .true.
   6
   23   24  5.063109E-03
   23   25  5.147713E-03
   23   26  4.596093E-03
   24   25  7.877069E-03
   24   26  7.426337E-03
   25   26  2.768402E-03
```

The decoherence rate for each pair of states is given by [[5]](#Jaeger_2012):
```
     rate = sqrt(alpha),
```
  where alpha is a parameter of gaussian
```
     gaus(alpha,t) = exp(-alpha*t*t)
```
  which is fitted to the dephasing function
```
    dij(t) = exp(-gij(t))
```
  with
```
    gij(t) = int^t_0 intg(t') dt'
```
  and
```
    int g(t) = int_t^0 Cij(t)
```
  Cij(t) is an autocorrelation function of the energy gap between two states
```
    Cij(t) = <(Eij(t)-Eij_ave)(Eij(0)-Eij_ave)>
           = <(Eij(t)*Eij(0)> - Eij_ave*Eij_ave
```
      where Eij(t) = Ei(t) - Ej(t)

#### Reference:
<a name="Gross_1990">1. "Time-dependent density-functional theory", Gross, E. K. U., and W. Kohn. <a href="https://www.sciencedirect.com/science/article/pii/S0065327608606000"> <i>Adv. Quantum Chem. </i> <b>21</b>, 255-291, (1990)</a><br></a>
<a name="Tully_1990">2."Molecular dynamics with electronic transitions." Tully, John C. <a href="http://aip.scitation.org/doi/abs/10.1063/1.459170"><i>J. Chem. Phys. </i> <b>93.2</b>, 1061-1071 (1990)</a><br></a>
<a name="Pulay_1980">3."Convergence acceleration of iterative sequences. the case of scf iteration." Pulay, P <a href="https://doi.org/10.1016/0009-2614(80)80396-4"><i> Chem. Phys. Lett. </i> <b>72.3</b>, 393-398 (1980)</a><br></a>
 <a name="Payne_1992">4."Iterative minimization techniques for ab initio total-energy calculations: molecular dynamics and conjugate gradients" Payne, M. C., Teter, M. P., Allan, D. C. ,Arias, T. A., and Joannapoulos, J. D. <a href="https://doi.org/10.1103/RevModPhys.64.1045"><i> Rev. Mod. Phys. </i> <b>64</b>, 1045
 (1992)</a><br></a>
<a name="Jaeger_2012">5. "Decoherence-induced surface hopping." Jaeger, H. M., Fischer, S., & Prezhdo, O. V. <a href="https://aip.scitation.org/doi/abs/10.1063/1.4757100"><i>J. Chem. Phys. </i> <b>137</b>, 22A545 (2012)</a><br></a>


### 5.4 linear-response TDDFT
#### 5.4.1 (On/Off)
**lrtddft** = [boolean]

Default: <br />
**lrtddft** = .FALSE.

.TRUE. : Execute NAQMD based on linear response time-dependent density functional theory (LR-TDDFT) [[1]](#Casida_1995) <br/>
.FALSE. : Do not use linear response theory for NAQMD.

Determines whether or not to use linear reponse theory, which involves calculation the Casida coupling matrix, for NAQMD simulations. For more detail, please see paper by Casida et al. [[1]](#Casida_1995)

#### 5.4.2 (whether to specify states)
**lrspecific** = [boolean]

Default: <br/>
**lrspecific** = .FALSE.

.TRUE. : Specify the occupied and unoccupied electronic states to be considered in the (specific states) subsection. <br/>
.FALSE. : Specify the occupied and unoccupied electronic states to be considered in the (energy difference) subsection.

Determines how the occupied and unoccupied states to be considered in the linear response TDDFT calculations.  When **lrspecific** is set to .TRUE. the occupied and unoccupied electronic states should be specified in the (specific states) subsection by their band index numbers.  When **lrspecific** is set to .FALSE. the occupied and unoccupied electronic states should be specified in the (energy difference) subsection by giving the maximum energy difference between them.

#### 5.4.3 (specific states)
**nlrstates** = [boolean] <br/>

**ihole, iparticle, ispin, ikpts** = [integer] [integer] [integer] [integer]<br/>
OR <br/>
**ihband, ipband** = [integer] [integer]

Default: <br>
**nlrstates** = 0 <br>
**ihband** = 0 <br>
**ipband** = 0


**nlrstates** = the number of states to specify <br/>
**ihole** = band index of hole (occupied states) <br/>
**iparticle** = band index of particle (unoccupied states) <br/>
**ispin** = spin index (1 or 2) <br/>
**ikpts** = k-point index <br/>
**ihband** = band index of hole (occupied states) <br/>
**ipband** = band index of particle (unoccupied states) <br/>

If **nlrstates** > 0, then states are specified with **ihole, iparticle, ispin, ikpts**.  If **nlrstates** = 0, then states are specified with **ihband, ipband**, where states between **ihband** and **ipband** will be considered. Note that this subsection is only read if **lrspecific** is set to .TRUE. in the (whether to specify states) subsection.

#### 5.4.4 (unit of energy)
(ry) or (hr) or (ev)

Default: <br/>
(ry)

The units of energy to be used for the variable **enediff** given below.

#### 5.4.5 (energy difference)
**enediff** = [real]

Default: <br>
**enediff** = 0.3

**enediff** = the maximum energy difference in the units given above between occupied and unoccupied states considered in linear response TDDFT.  Note that this subsection is only read if **lrspecific** is set to .FALSE. in the (whether to specify states) subsection.

#### 5.4.6 (finite-difference mesh)
**fdmeshxc** = [real] <br/>

Default: <br>
**fdmeshxc** = 1.d-06

**fdmeshxc** is the finite difference mesh size for the 2nd derivative of Exc 

#### 5.4.7 (threshold)
**threshrpa** = [real] <br>

Default: <br>
**threshrpa** = 0.0

**threshrpa** is the threshold value for RPA1

#### 5.4.8 (transition-dipole threshold)
**threshdipole** = [real] <br/>

Default: <br>
**threshdipole** = 0.0

**threshdipole** threshold value for transition dipole.

#### 5.4.9 (energy-difference threshold)
**threshexcite** = [real] <br>

Default: <br>
**threshexcite** = 1.d-04

**threshexcite** is the threshold value for the difference between the energies of hole and particle

#### 5.4.10 (long-range exchange scheme)
**llcexchange** = [boolean]

Default: <br/>
**llcexchange** = .FALSE.

.TRUE. : Use long-range correction for exchange. <br>
.FALSE. : Do not use.

Determines whether or not to use the long-range correction for exchange functional [[2]](#Tawada_2004).  Note that if **lvacuum[1|2|3]** is set to .TRUE. **alpha\_ldouble\_grid_recip** needs to be specified to determine the ratio of the long-range part.


#### 5.4.11 (parameters to divide 1/r)
**alpha_llcexchange** = [real] <br>
**large\_cell\_llcexchange** = [real]

Default: <br>
**alpha_llcexchange** = 0.17 <br>
**large\_cell\_llcexchange** = 2.0

Note 1: The cell-vector length of the CUBIC large cell is given by large_cell_llcexchange * max(L1,L2,L3). <br/>
Note 2: (iosp) & (ecutlong) have to be set in *double-grid method
only if llcexchange & .not.ldouble_grid_recip

#### 5.4.12 (hybrid mixing parameter)
**hybrid_mixing** = [real]

Default: <br>
**hybrid_mixing** = 0.0

**hybrid_mixing** is the mixing paramter for the Hybrid method. Note that this is only read when **llcexchange** is set to .FALSE. 

#### 5.4.13 (Foerster: width of Gaussian)
**foersigma** = [real] <br/>

Default: <br>
**foersigma** = 0.1d0

**foersigma** gives the Foerster transfer rate

#### 5.4.14 (diagonal approximation)
**ldiagonal** = [boolean]

Default: <br>

**ldiagonal** = .TRUE.

.TRUE. :  The diagonal approximation is made, where matrix elements are calculated only for **ihole\_1 == ihole\_2** and **iparticle\_1 == iparticle\_2** <br/>
.FALSE. : Diagonal approximation is not used.

Determines whether or not to use the diagonal approximation. <br/>
Note: If **ltddft** is set to .TRUE. .and. **ltddft_fssh** is also set to .TRUE., set **ldiagonal** equal to .TRUE. <br/>

#### 5.4.15 (scissors correction)
**lscissors\_lrtddft** = [boolean] <br>
**ecorr\_triplet** = [real] <br>
**ecorr\_singlet** = [real]

Default: <br>
**lscissors\_lrtddft** = .FALSE.

**lscissors\_lrtddft** = .TRUE. : Implement scissors correction.
**lscissors\_lrtddft** = .FALSE. : Do no implement correction.


Note: The values for **ecorr\_triplet** and **ecorr\_singlet** depend on 'unit of energy' in the *linear-response TDDFT section.  These corrections will be added to the excitation energies.

#### 5.4.16 (singlet-fission-rate calculation)
**lfission_rate** = [boolean]

Default:
**lfission_rate** = .FALSE.

.TRUE. : Calculate the singlet fission rate based on time-dependent perturbation theory. <br>
.FALSE. : Do not calculate rate.

#### 5.4.17 (spontaneous emission rate)
**lemission_rate** = [boolean] <br>
**refractive** = [real]

Default: <br>
**lemission_rate** = .TRUE. <br>
**refractive** = 1.5

**lemission_rate** = .TRUE. : Compute spontaneous emission rate based on transition dipole approximation. <br>
**lemission_rate** = .FALSE. : Do not calculate rate.

Determines whether or not to calculate spontaneous emission rate.  **refractive**  is the refractive index of the system.

#### Reference: 
<a name="Casida_1995">1. "Time-Dependent Density Functional Response Theory for Molecules" M. E. Casida, <a href="http://homepage.univie.ac.at/mario.barbatti/papers/method/tddft_casida_1995.pdf"> Recent Advances in Density Functional Methods (Part I), edited by D. P. Chong (World Scientific, Singapore, pp. 155-192 (1995))</a><br></a>
<a name="Tawada_2004">2. "A long-range-corrected time-dependent density functional theory." Tawada, Yoshihiro, et al.  <a href="https://aip.scitation.org/doi/abs/10.1063/1.1688752"><i> J. Chem. Phys. </i> <b>120</b>, 8425-8433 (2004) </a><br></a>


### 5.5 PAW
#### 5.5.1 (On/Off)
**lpaw** = [boolean]

Default: <br>
**lpaw** = .FALSE.

.TRUE. : Use the PAW method [[1]](#Blochl_1994)<br>
.FALSE.: Use other pseudopotential method.

#### 5.5.2 (non-spherical symmetry)
**lpaw_sym** = [boolean]

Default: <br>
**lpaw_sym** = .FALSE.

.TRUE. : Full symmetry <br>
.FALSE. : Spherical symmetry

#### 5.5.3 (onsite charge mixing)
**paw_mix** = [real]

Default: <br>
**paw_mix** = 0.5

Gives the onsite charge mixing parameter in the PAW method.

#### Reference:
<a name="Blochl_1994">1. "Projector augmented-wave method" P. E. Blochl <a href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.50.17953"><i>Phys Rev B</i> <b>50</b>, 17953 (1994)</a><br></a>

### 5.6 Cluster
#### 5.6.1 (On/Off)
**lclust** = [boolean]

Default: <br>
**lclust** = .FALSE.

.TRUE. : Perform cluster calculation. <br>
.FALSE. : Perform bulk calculation.

Determines whether to perform a cluster or bulk calculation.

### 5.7 Approximation for Exchange
#### 5.7.1 (approximation)
**jgga** = 1|2|3|4|5|6

Default: <br>
**jgga** = 2

| **jgga**  | Excahnge functional |
| --------- | :-----------------: |
| **1**     | LD[[1]](#Perdew_1981)              |
| **2**     | GGA(PBE)[[2]](#Perdew_1996)        |
| **3**     | GGA(RPBE)[[3]](#Hammer_1999)       |
| **4**     | GGA(revPBE)[[4]](#Zhang_1998)     |
| **5**     | vdW-DF[[5]](#Dion_2004)          |
| **6**     | vdW-DF2 [[6]](#Lee_2010)        |


**jgga** is used to specify the exchange correction functional.

#### 5.7.2 (DFT-D)
**ldftd** = [boolean]

Default: <br>
**ldftd** = .FALSE.

.TRUE.  : Employ an empirical vdW correction <br>
.FASLE. : Do not use an empirical vdW correction

Determines whether to use an empirical correction for the van der Waals interaction proposed by S. Grimme [[7]](#grimme_2006)

#### 5.7.3 (kernel factorization)  ------  jgga=5 or 6 (vdW-DF or vdW-DF2)
**lvdw_factorize** = [boolean]

Default: <br>
**lvdw_factorize** = .FALSE.

.TRUE.  : Take into account. <br>
.FASLE. : Do not

only used for for vdW-DF or vdW-DF2 (jgga=5 or 6 )

#### 5.7.4 (q-mesh range)  ------  jgga=5 or 6 (vdW-DF or vdW-DF2)
**lvdw_factorize** = [boolean] <br>
**q_cut** = [real] <br>
**dlambda** = [real] <br>
**nqmesh** = [integer]

Default: <br>
**q_cut** = 5.0 <br>
**dlambda** = 1.2 <br>
**nqmesh** = 20 

These parameters are used to saturate the original function q_0(n,|Delta n|) by redefining q_0^{sat} = h[q_0,q_c].

h(x,x_c) = x_c[1 - exp( - sum_{m=1}^{m_c} (x/x_c)^m/m)]

**q_cut** is corresponding to q_c and **dlambda** and **nqmesh** are a width and the number of a logarithmic mesh.


#### 5.7.5 (kernel cutoff parameters)  ------  jgga=5 or 6 (vdW-DF or vdW-DF2)
**d_soft** = [real] <br>
**phi_at_0** = [real]

Default: <br>
**d_soft** = 1.0 <br>
**phi_at_0** = 0.5

Because the kernel phi(d\_1,d\_2) has a logarithmic divergence when d\_1 and d\_2 tend to zero, the modified kernel is used for small d\_1 and d\_2:

phi(d\_1,d\_2) = phi\_0 + phi\_2 d^2 + phi\_4 d^4    (d < d_s)

where d = sqrt(d\_1^2 + d\_2^2) and phi\_2 and phi\_4 are fitting parameters. d\_s and phi\_0 are corresponding to **d\_soft** and **phi\_at\_0**


#### 5.7.6 (cutoff length)  ------  lvdw_factorize = .FALSE.
**rvdwmax** = [real]

Default: <br>
**rvdwmax** = 10.d0

Cutoff length in [a.u]

#### 5.7.7 (integration skip)  ------  lvdw_factorize = .FALSE.
**imod1** = [even integer] <br>
**imod2** = [integer]

Default: <br>
**imod1** = 2 <br>
**imod2** = 3

Integration of r' on every **imod1** meshes <br>
Integration of r  on every **imod2** meshes


#### 5.7.8 (order of cardinal B spline)  ------  lvdw_factorize = .FALSE.
**iosp** = [odd integer]

Default: <br>
**iosp** = 7


order of cardinal B spline in the double-grid method [[8]](#Martyna_1999)


#### 5.7.9 (DFT+U)
**lplusU** = [boolen]

Default: <br>
**lplusU** = .FALSE.

.TRUE.  : Employ DFT+U correction  of the mean-field Hubbard model <br>
.FASLE. : Do not use DFT+U correction <br>
Determines whether to use DFT+U correction  of the mean-field Hubbard model [[9]](Dudarev_1998). The value of U and J parameter is set to in the \*atoms section. Note that this subsection is only considered when **lpaw** is set to .TRUE.


#### 5.7.10 (onsite charge mixing)
**plusU_mix** = [real]

Default: <br>
**plusU_mix** = 0.5

**plusU_mix** is used as turning parameter for onsite charge mixing. Note that this subsection is only considered when both **lpaw** and **lplusU** are set to .TRUE.


#### 5.7.11 (DFT+C)
**lplusC** = [boolen]

Default: <br>
**lplusC** = .FALSE.

.TRUE.  : Employ an empirical correction to non-local pseudopotential [[10]](#Li_2005) <br>
.FASLE. : Do not


#### 5.7.12 (hybrid functional)
**jhybrid** = [0|1]

Default: <br>
**jhybrid** = 0

| **jhybrid**  | Hybrid functional   |
| ------------ | :-----------------: |
| **0**        | Do not use          |
| **1**        | HSE  [[11-13]](#Heyd_2003)          |

**jhybrid** is used to specify the range-separated hybrid exchange functionals. Note that this subsection is only considered when **jgga** is set to 2(GGA(PBE)).


#### 5.7.13 (hybrid mixing parameter)
**hmixing** = [real]

Default: <br>
**hmixing** = 0.25

The HSE exchange correlation functional is constructed as a linear combination of the Hartree-Fock exact exchange functional and PBE functional:

Exc^HSE = a Ex^HF + (1-a) Ex^PBE-SR + Ex^PBE-LR + Ec^PBE

where **a** is defined by **hmixing**.


#### 5.7.14 (range-separation parameter)
**hrange** = [real]

Default: <br>
**hrange** = 0.11

The HSE exchange correlation functional spits the Coulomb operator into short-range and long-range components by an error function and **hrange**, omega is an adjustable parameter:

1/r = erfc(omega r)/r + erf(omega r)/r


#### 5.7.15 (grid-reduction factor)
**mgridred** = [integer]

Default: <br>
**mgridred** = 1

**mgridred** gives the grid-reduction factor.
Note: The number of k-point grids should be the multiple of **mgridred**.

#### Reference:
<a name="Perdew_1981">1. "Self-interaction correction to density-functional approximations for many-electron systems" J. P. Perdew and Alex Zunger <a href="https://doi.org/10.1103/PhysRevB.23.5048"> <i>Phys. Rev. B</i> <b> 23 </b>, 5048(1981)</a><br></a>
<a name="Perdew_1996">2. "Generalized Gradient Approximation Made Simple" J. P. Perdew, K. Burke and M. Ernzerhof <a href="https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.77.3865"> <i>Phys. Rev. Lett.</i> <b>78</b>, 3865 (1996)</a><br></a>
<a name="Hammer_1999">3. "Improved adsorption energetics within density-functional theory using revised Perdew-Burke-Ernzerhof functionals" B. Hammer, L. B. Hansen, and J. K. Nørskov <a href="https://doi.org/10.1103/PhysRevB.59.7413"> <i>Phys. Rev. B</i> <b>59</b>, 7413 (1999)</a><br></a>
<a name="Zhang_1998">4. "Comment on “Generalized Gradient Approximation Made Simple” " Yingkai Zhang and Weitao Yang <a href="https://doi.org/10.1103/PhysRevLett.80.890"> <i>Phys. Rev. Lett.</i> <b>80</b>, 890 (1998)</a><br></a>
<a name="Dion_2004">5. "Van der Waals Density Functional for General Geometries " M. Dion, H. Rydberg, E. Schröder, D. C. Langreth, and B. I. Lundqvist <a href="https://doi.org/10.1103/PhysRevLett.92.246401"> <i>Phys. Rev. Lett.</i> <b>92</b>, 246401 (2004)</a><br></a>
<a name="Lee_2010">6. "Higher-accuracy van der Waals density functional" Kyuho Lee, Éamonn D. Murray, Lingzhu Kong, Bengt I. Lundqvist, and David C. Langreth <a href="https://doi.org/10.1103/PhysRevB.82.081101"> <i>Phys. Rev. B</i> <b>82</b>, 081101 (2010)</a><br></a>
<a name="Grimme_2006">7. "Semiempirical GGA‐type density functional constructed with a long‐range dispersion correction" Stefan Grimme <a href="https://doi.org/10.1002/jcc.20495"> <i>J. Comput. Chem.</i> <b>27</b>, 1787 (2006)</a><br></a>
<a name="Martyna_1999">8. "A reciprocal space based method for treating long range interactions in ab initio and force-field-based calculations in clusters" <a href="https://doi.org/10.1063/1.477923"> <i>J. Chem. Phys.</i> <b>110</b>, 2810 (1990)</a><br></a>
<a name="Dudarev_1998">9. "Electron-energy-loss spectra and the structural stability of nickel oxide:  An LSDA+U study" <a href="https://doi.org/10.1103/PhysRevB.57.1505"> <i>Phys. Rev. B</i> <b>57</b>, 1505 (1998)</a><br></a>
<a name="Dudarev_1998">10. "Band-structure-corrected local density approximation study of semiconductor quantum dots and wires" <a href="https://doi.org/10.1103/PhysRevB.72.125325"> <i>Phys. Rev. B</i> <b>72</b>, 125325 (2005)</a><br></a>
<a name="Heyd_2003">11. "Hybrid functionals based on a screened Coulomb potential." J. Heyd and G. E. Scuseria<a href="https://aip.scitation.org/doi/abs/10.1063/1.1564060"><i>J. Chem. Phys.</i> <b>118</b>, 8207-8215 (2003)</a></br></a>
<a name="Heyd_2004_1">12. "Efficient hybrid density functional calculations in solids: Assessment of the Heyd–Scuseria–Ernzerhof screened Coulomb hybrid functional " J. Heyd and G. E. Scuseria<a href="https://doi.org/10.1063/1.1760074"><i>J. Chem. Phys.</i> <b>121</b>, 1187 (2004)</a></br></a>
<a name="Heyd_2004_2">13. "Influence of the exchange screening parameter on the performance of screened hybrid functionals
" A.V. Krukau, O. A. Vydrov, A. F. Izmaylov, and G. E. Scuseria<a href="https://doi.org/10.1063/1.2404663"><i>J. Chem. Phys.</i> <b>125</b>, 224106 (2006)</a></br></a>



### 5.8 scissors correction at donor/acceptor interface
#### 5.8.1 (On/Off)
**lscissors** = [boolean]

Default: <br>
**lscissors** = .FALSE.

.TRUE.  : Take into account <br>
.FASLE. : Do not

#### 5.8.2 (unit of energy)
(ry) or (hr) or (ev)


#### 5.8.3 (corrections)
**ecorr_acceptor** = [real] <br>
**ecorr_donor** = [real]

Default <br>
**ecorr_acceptor** = 0.0 <br>
**ecorr_donor** = 0.0

Corrections to be subtracted from the eigenenergies.

#### 5.8.4 (unit of length)
(bohr) or (ang)

#### 5.8.5 (boundary)
**zboundary** = [real]

Default: <br>
**zboundary** = 0.0
Defines the boundary between acceptor and donor regions as follows: <br>
z-coordinate < **zboundary**: acceptor region <br>
z-coordinate > **zboundary**: donor region

### 5.9 k-points
#### 5.9.1 (number of k points)
**nkpnt** = [integer]

Default: (gamma point)<br>
**nkpnt** = 1

**nkpnt** is the number of k points. If **nkpnt** is greater than zero, you set k vectors in (k vectors). If **nkpnt** is zero, you set the number of division for the reciprocal vectors in (division number).


#### 5.9.2 (k vectors)
**(bzk,wbzk)** = [real] [real] [real] [real]

Default: (gamma point)<br>
**bzk(1:3)** = (/0.0 0.0 0.0/) <br>
**wbzk** = 1.0

**bzk(1:3)** are k vectors and **wbzk** is the weight of each k-vector. Note that this subsection is only considered when **nkpnt** is greater than zero.Also, Please set K vectors between -0.5 and 0.5. 


#### 5.9.3 (division number)
**npkx(1:3)** = [integer] [integer] [integer]

Default: <br>
**npkx(1:3)** = (/ 1 1 1 /)

**npkx(1:3)** is the number of division for reciprocal vectors and Monkhorst-Pack method is used. [[1]](#Monkhorst_1976) 


#### 5.9.4 (tetrahedron method)
**ltetra** = [boolean]

Default: <br>
**ltetra** = .FALSE.

.TRUE.  : Use the tetrahedron method <br>
.FALSE. : Do not use

Determines whether to use the modified tetrahedron method [[2]](#Blochl_1994) for the Brillouin Zone integration.

#### Reference: 
<a name="Monkhorst_1976">1. "Special points for Brillouin-zone integrations" H. J. Monkhorst and J. D. Pack, O. Jepsen, and O. K. Andersen. <a href="https://doi.org/10.1103/PhysRevB.13.5188"><b>Phys. Rev. B</b> <i>13</i>, 5188 (1976)</a></br></a>
<a name="Blochl_1994">2. "Improved tetrahedron method for Brillouin-zone integrations."  P. E. Blöchl, O. Jepsen, and O. K. Andersen. <a href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.49.16223"><b>Phys. Rev. B</b> <i>49</i>, 16223 (1994)</a></br></a>

### 5.10 electric field
#### 5.10.1 (On/Off)
**lefield** = [boolean]

Default: <br>
**lefield** = .FALSE.

.TRUE. : Apply uniform electric field. <br>
.FALSE. : No electric field.

#### 5.10.2 (electric field)
**efield(1:3)** = [real] [real] [real]

Determines the uniform electric field in atomic units. [[1]](#Umari_2001)

#### Reference:
<a name="ref19">19. "Ab initio molecular dynamics in a finite homogeneous electric field."   P. Umari, and A. Pasquarello. <a href="https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.89.157602"><b> Phys. Rev. Lett.</b> <i>89</i>,157602 (2002)</a></br></a>

### 5.11 spin polarization
#### 5.11.1 (On/Off)
**lspin** = [boolean]

Default: <br>
**lspin** = .FALSE.

.TRUE.  : Apply spin polarization
.FALSE. : Do not use

Determines whether to take into account spin polarization.

#### 5.11.2 (noncollinear magnetism)
**lnoncollinear** = [boolean]

Default: <br>
**lnoncollinear** = .FALSE.

.TRUE.  : Noncollinear magnetism
.FALSE. :    Collinear magnetism

Note that the DFT calculation by noncollinear magnetism is still under development (SCF iterations might not converge).


#### 5.11.3 (fix spin polarization)
**lfixud** = [boolean]

Default: <br>
**lfixud** = .TRUE.

.TRUE.  : Fix the value of spin polarization <br>
.FALSE. : Do not fix

**lfixud** set to .TRUE. will be fixed the value of spin polarization. Note that this subsection is only considered when **lnoncollinear** is .FALSE.


#### 5.11.4 (spin polarization)
**diffud** = [real]

Default: <br>
**diffud** = 0.0

**diffud** is the value of spin polarization to fix total magnetic moment. Note that this subsection is only considered when **lfixud** is .TRUE.


#### 5.11.5 (initial wave functions)
**lwfrand** = [boolean]

Default: <br>
**lfixud** = .FALSE.

.FALSE. : Down spin is given as same as up spin <br>
.TRUE.  : Down spin is given by random number

Determines how to give the down spin.

#### 5.11.6 (initial spin density)
**inispin** = [0|1|2|3|4]

Default: <br>
**inispin** = 1

| **inispin**  | Initial spin density               |
| ------------ | :--------------------------------: |
| **0**        | from LSDA (lnoncollinear = .TRUE.) |
| **1**        | uniformly polarized                |
| **2**        | ferromagnetic alignment            |
| **3**        | antiferromagnetic random alignment |
| **4**        | specify spin polarization par atom |

**inispin** is used to specify the method for initial spin density. The detail parameter is set in another subsection.


#### 5.11.7 (uniform spin density)  ------  **inispin** = 1 & **lnoncollinear** = .TRUE.

**atmmagne** = [real] <br>
**umx**, **umy**, **umz** = [real] [real] [real]

If **inispin** is one, spin polarized uniformly. **atmmagne** is total magnetic moment and (**umx**, **umy**, **umz**) are a magnetization direction.


#### 5.11.8 (ferromagnetic spin density)  ------  **inispin** = 2 <br>
**nspinat** = [integer] <br>
**iatspin** = [integer] <br>
**atmmagne** = [real] <br>
(**umx**, **umy**, **umz**) = [real] [real] [real]

Default: <br>
**nspinat** = 1 <br>

All **itspin**-type atoms are polarized with the magnetic moment of **atmmagne**. (**umx**, **umy**, **umz**) are available only when **lnoncollinear** is .TRUE. **nspinat** is the number of atom types to be polarized.

Example:
```
(ferromagnetic spin density)
   2                          : (nspinat)
   1    2.d0  0.d0 0.d0 1.d0  : (iatspin, atmmagne, umx, umy, umz)
   2   -4.d0  0.d0 0.d0 1.d0  : (iatspin, atmmagne, umx, umy, umz)
```


#### 5.11.9 (antiferromagnetic spin density)  ------  **inispin** = 3 <br>
**nspinat** = [integer] <br>
**iatspin** = [integer] <br>
**atmmagne** = [real] <br>
(**umx**, **umy**, **umz**) = [real] [real] [real]

Default: <br>
**nspinat** = 1 <br>

Half of **itspin**-type atoms are polarized with the magnetic moment of **atmmagne**, and the rest are polarized with the magnetic moment of **-atmmagne**. (**umx**, **umy**, **umz**) are available only when **lnoncollinear** is .TRUE. **nspinat** is the number of atom types to be polarized.

Example:
```
(antiferromagnetic spin density)
   2                          : (nspinat)
   1    2.d0  0.d0 0.d0 1.d0  : (iatspin, atmmagne, umx, umy, umz)
   2   -4.d0  0.d0 0.d0 1.d0  : (iatspin, atmmagne, umx, umy, umz)
```


#### 5.11.10 (atomic spin density)  ------  **inispin** = 4
**nspinat** = [integer] <br>
**iatspin** = [integer] <br>
**atmmagne** = [real] <br>
(**umx**, **umy**, **umz**) = [real] [real] [real]

Default: <br>
**nspinat** = 1 <br>

You can specify the magnetic moment of **atmmagne** par atom. **iatspin** is the index of atom number. (**umx**, **umy**, **umz**) are available only when **lnoncollinear** is .TRUE. **nspinat** is the number of atom types to be polarized.

Example:
```
(atomic spin density)
   4                          : (nspinat)
  10    2.d0  0.d0 0.d0 1.d0  : (iatspin, atmmagne, umx, umy, umz)
  11    2.d0  0.d0 0.d0 1.d0  : (iatspin, atmmagne, umx, umy, umz)
  23    2.d0  0.d0 0.d0 1.d0  : (iatspin, atmmagne, umx, umy, umz)
  24    2.d0  0.d0 0.d0 1.d0  : (iatspin, atmmagne, umx, umy, umz)
```


#### 5.11.11 (reference direction)  ------  **lnoncollinear** = .TRUE.
**refmx**, **refmy**, **refmz** = [real] [real] [real]

Default: <br>
**refmx**, **refmy**, **refmz** = (/ 0.0 0.0 1.0/)

When **lnoncollinear** is .TRUE., **refmx**, **refmy**, **refmz** shows the reference direction.


#### 5.11.12 (duplicate the up-spin state)  ------  **lnoncollinear** = .FALSE.
**lduplicate** = [boolean]

Default: <br>
**lduplicate** = .FALSE.

.TRUE.  : Duplicate wavefunction of the up spin state <br>
.FALSE. : Do not duplicate

### 5.12 SCF iterations
#### 5.12.1 (global iterations)
**iscfmx** = [integer]

Default: <br>
**iscfmx** = 100

Gives the maximum number of global iterations of SCF to be performed.

#### 5.12.2 (trial global iterations)
**itrial** = [integer]

Default: <br>
**itrial** = 0

Gives the number of trial global iterations

#### 5.12.3 (trial charge mixing)
**imxchg_trial** = [0|1] <br>
**aslh\_trial, bslh\_trial** = [real] [real]

Default: <br>
**imxchg_trial** = 0 <br>
**aslh\_trial** = 0.1 <br>
**bslh\_trial** = 0.64

**imxchg_trial** = 0 : No trial charge mixing <br>
**imxchg_trial** = 1 : Trial charge mixing

Determines whether to perform trial charge mixing, and if so, the tuning parameters, **aslh\_trial, bslh\_trial**, to be used.

#### 5.12.4 (trial spin-density mixing)
**imxspin_trial** = [0|1] <br>
**amxspin\_trial, bmxspin\_trial** = [real] [real]

Default:
**imxspin_trial** = 0 <br>
**amxspin\_trial** = 0.1 <br>
**bmxspin\_trial** = 0.64


#### 5.12.5 (tolerances)
**tolpot** = [real] <br>
**tolres** = [real]

Default: <br>
**tolpot** = 5.0d-09 <br>
**tolres** = 5.0d-09

**tolpot** gives the tolerance for change in total energy.  Once the change in total energy between subsequent SCF iterations is smaller than the tolerance, the SCF iterations are considered to have converged.  **tolres** similarly gives the tolerance for average residual. 

#### 5.12.6 (HC products)
**lhcunt** = [boolean]

Default: <br>
**lhcunt** = .TRUE.

.TRUE. : HC products by Unitary tr. <br>
.FALSE. : ???

???

#### 5.12.7 (charge mixing)
**imxchg** = [0|1|2|3|4|5|6]
**aslh, bslh** = [real] [real]

Default: <br>
**imxchg** = 1


| imxchg | Charge Mixing Method | 
| ------ | :------------------: | 
| **0**  | No Mixing            |  
| **1**  | Pulay                |
| **2**  | Anderson             |
| **3**  | Simple               |
| **4**  | Srivastava           |
| **5**  | Johnson              |
| **6**  | Johnson w/ variable weights |

Determines which charge mixing method will be used, along with tuning parameters for that method. 

#### 5.12.8 (number of mixing)
**itratn** = [integer]

Default: <br>
**itratn** = 10

Determines how many charge densities from previous SCF iterations to use for charge mixing. Not available for imxchg = [0|3].

#### 5.12.9 (metric in charge mixing)
**cmetric** = [real]

Default: <br>
**cmetric** = 20.0

Gives the metric for charge mixing. Used as the weighting factor w(g) in the calculation of residual is given by:

 w(g) = (cmetric-1)*G_min^2*G_max^2/(G_max^2 - cmetric*G_min^2)

Not available for imxchg = [0|3].

#### 5.12.10 (magnetic moment mixing)
**lmixreal** = [boolean]

Default: <br>
**lmixreal** = .TRUE.

.TRUE. :real-space mixing <br>
.FALSE. : reciprocal-spcace mixing <br>

This is only read when **lnoncollinear** in the (noncollinear magnetism) subsection of the spin polarization section is set to .TRUE.

#### 5.12.11 (spin-density mixing)
**imxspin** = [0|1|2|3|4|5|6] <br>
**amxspin, bmxspin** = [real] [real]

Default: <br>
**imxspin** = 0 <br>
**amxspin** = 0.8 <br>
**bmxspin** = 0.64

| imxspin | Charge Mixing Method |
| ------  | :------------------: | 
| **0**   | No Mixing            |
| **1**   | Pulay                |
| **2**   | Anderson             |
| **3**   | Simple               |
| **4**   | Srivastava           |
| **5**   | Johnson              |
| **6**   | Johnson w/ variable weights | 

This is only calculated when spin polarization is enabled in the *spin polarization section.

#### 5.12.12 (# of spin-density mixing)
**nmxspin** = [integer]

Default: <br>
**nmxspin** = 10

Determines how many spin densities from previous SCF iterations to use for spin density mixing.

#### 5.12.13 (spin metric in spin-density mixing)
**spinmetric** = [real]

Default: <br>
**spinmetric** = 1.0

The metric in spin-density mixing.  This variable is only read when **lnoncollinear** == .FALSE. or **lmixreal** == .FALSE.

#### 5.12.14 (Kerker-mixing weight in spin-density mixing)
**wkerker** = [real]

Default: <br>
**wkerker** = 1.0

Kerker-mixing weight in spin-density mixing.  This variable is only read when **lnoncollinear** == .FALSE. or **lmixreal** == .FALSE.

  Caution!!  sigma_new(G) = sigma_in(G) + K(G)*(sigma_out(G) - sigma_in(G)),
             where sigma is the spin density, and
             K(G) = amxspin*(1 - wkerker + wkerker*G~2/(G^2 + bmxspin)).
  Note that sigma(G=0) is not updated when wkerker = 1.d0,
  i.e., the difference between the numbers of up- and down-spin is fixed.

Not available for imxchg = [0|3].

#### 5.12.15 (output eigenvalues)
**louteig** = [boolean]

Default: <br>
**louteig** = .FALSE.

.TRUE. : Output eigenvalues to file <br>
.FALSE. : Do not output eigenvalues

Whether or not to output eigenvalues through SCF iterations.

#### 5.12.16 (output energy parts)
**loutenergy** = [boolean]

Default: <br>
**loutenergy** = .FALSE.

.TRUE. : Output energy parts to file <br>
.FALSE. : Do not output energy parts

Whether or not to output energy parts through SCF iterations

#### 5.12.17 (output residuals)
**loutzansa** = [boolean]

Default: <br>
**loutzansa** = .FALSE.

.TRUE. : Output residuals (standard deviation) to file <br>
.FALSE. : Do not output residuals.

Whether or not to output residuals through SCF iterations

#### 5.12.18 (convergence stabilizer)
**nstabi** = [integer] <br>
**xstabi** = [real]

Default: <br>
**nstabi** = 10 <br>
**xstabi** = 10.0

If zansa2 < zansa2/xstabi at nstabi iterations before, set lreset = .FALSE.

### 5.13 well potenial
#### 5.13.1 (On/Off)
**lwell** = [boolean]

Default: <br>
**lwell** = .FALSE.

Determines whether to use a well potential outide atomic sphere.

#### 5.13.2 (height of potential)
**wellheight** = [real]

Default: <br>
**wellheight** = 1.0

Gives the height of the well potenital in Rydbergs.

#### 5.13.3 (radius around each atom)
**nwellat** = [integer] <br>
**iatwell, radwell** = [integer] [real]

Default: <br>
**nwellat** = 0

**nwellat** gives the total number of atom types for which to make a well potential around, followed by that many lines of **iatwell**, giving the atom type index and **radwell**, giving the radius around each atom in bohrs.  If not specified (atom type  2 in this case), the covalent radius is used.

Example: This defines the radius around atom types 1,3, a total of 2 atom types, to be 1.9 bohrs.  Atom type 2, for which no radius is given, will have its covalent radius used.
```
(radius around each atom)
   2
   1    1.9d0
   3    1.9d0
```

### 5.14 Kohn-Sham equation
#### 5.14.1 (On/Off)
**ihldam** = [1|2]

Default: <br>
**ihldam** = 1

1 : Conjugate-Gradient (CG) Method
2 : Residual minimization scheme, direct inversion in the iterative subspace (RMM-DIIS) Method

Determines the method used in the iterative minimization of energy as a functional of wavefunctions.

#### 5.14.2 (tolerance)
**toleig** = [real]

Default: <br>
**toleig** = 1.d-10

Gives the tolerance for the Kohn-Sham equations.

#### 5.14.3 (threshold for w.f. direction)
**threinn** = [real]

Default: <br>
**threinn** = 0.0

Gives the threshold for direction of new wavefunction

#### 5.14.4 (iteration)
**itermx** = [integer]

Default: <br>
**itermx** = 4

Gives the maximum number of iterations.


#### 5.14.5 (empty-band iteration)
**iteremx** = [integer]

Default: <br>
**iteremx** = 4

Gives the maximum number of iterations for empty bands.

#### 5.14.6 (trial global iterations)
**kstrial** = [integer]

Default: <br>
**kstrial** = 0

Gives the number of trial global iterations with itermx = iteremx = 1

#### 5.14.7 (CG method)
**methodcg** = [1|2]

Default: <br>
**methodcg** = 2

1:line minimization <br>
2:BKL [ref]

Gives the method used of CG
This variable is only read if **ihldam** == 1

### 5.15 Poisson equation
Note this section is only used for cluster calculations.
#### 5.15.1 (multigrid level)
**multg** = [integer]

Default: <br>
**multg** = 2

??

#### 5.15.2 (tolerance)
**tolcg** = [real]

Default: <br>
**tolcg** = 1.d-11

Gives the tolerance energy in a.u. for one electron

#### 5.15.3 (preconditioner)
**weigrd** = [real]

Default: <br>
**weigrd** = 0.4

???

#### 5.15.4 (differentiation)
**nd2v** = [integer]

Default: <br>
**nd2v** = 6

Gives the order of numerical differentiation.

#### 5.15.5 (mesh for serial calculation)
**msrhmx, msrhmy, msrhmz** = [integer] [integer] [integer]

Default: <br>
**msrhmx, msrhmy, msrhmz** = 14 14 14

Gives the mesh size for serial calculation.


### 5.16 molecular dynamics
#### 5.16.1 (On/Off)
**ifmd** = [0|1|2|3|4|5|10]

Default: <br>
**ifmd** = 0

| ifmd   | Type of Dynamics  |
| ---    | :-------------:   |
| **0**  | Single Step       |
| **1**  | Optimization      |
| **2**  | NVE               |
| **3**  | NVT               |
| **4**  | NPT               |
| **5**  | NVT for each atom |
| **10** | MSST [ref]        |

Determines the type of QMD simulation to run.  Add short description of each type.

#### 5.16.2 (time step)
**dtmd, nstop** = [real] [integer]

Default: <br>
**dtmd** = 50.0 <br>
**nstop** = 10

**dtmd** gives the time step for QMD simulation in [a.u.], while **nstop** gives the total number of time steps to simulate.  Thus, the total simulation time will equal (**nstop** * **dtmd**)


#### 5.16.3 (initial step number)
**nstep_ini** = [integer]

Default: <br>
**nstep_ini** = 0

Gives the initial step number.  This varibale is ignored for **lstart** == .TRUE.

#### 5.16.4 (temperature)
**treq** = [real]

Default: <br>
**treq** = 300.0

Gives the initial temperature in [K] for NVE, NVT, and NPT QMD simulations.

#### 5.16.5 (check temperature)
**liscale** = [boolean] <br>
**iscnum** = [integer] <br>
**iscstp** = [integer]

Default: <br>
**liscale** = .FALSE. <br>
**iscnum** = 25 <br>
**iscstp** = 20 <br>

**liscale** = .FALSE. : Do not scale temperature. <br>
**liscale** = .TRUE.  : Scale temperature a total of **iscstp** times, with **iscnum** steps in between each scaling.

Determines whether and how to scale temperature to keep it near the initial given temperature.

#### 5.16.6 (make total momentum zero)
**lmomzero** = [boolean]

Deafault: <br>
**lmomzero** = .FALSE.

.TRUE. : Make the total momentum zero.
.FALSE. : Do not

#### 5.16.7 (optimization)
**ioptmze** = [-1|0|1|2|3|10]

Default: <br>
**ioptmze** = 2

| ioptmze | Type of Structural Optimization       |
| ---     | :-------------------------------:     |
| **-1**  | Do not optimize atomic coords         |
| **0**   | Conjugate gradient                    |
| **1**   | Projected velocity Verlet             |
| **2**   | Quasi-Newton method with BFGS formula |
| **3**   | RFO saddle points search              |
| **10**  | Harmonic-mode analysis                |

Determines method for structural optimization of atomic coordinates.  This varibale is onlt read when **ifmd** == 1.

#### 5.16.8 (cell optimization)
**ioptmze_cell** = [-1|0|1|2]

Default: <br>
**ioptmze_cell** = -1

| ioptmze_cell | Type of Cell Optimization             |
| ------------ | :---------------------------------:   |
| **-1**       | Do not optimize supercell             |
| **0**        | Conjugate gradient                    |
| **1**        | Not yet implemented                   |
| **2**        | Quasi-Newton method with BFGS formula |

Determines method for optimization of (super)cell.  This varibale is only read when **ifmd** == 1, and is not not read when **ioptmze** == 1 or 10.

#### 5.16.9 (cell CG time step)
**dtcellcg** = [real]

Default: <br>
**dtcellcg** = 0.1

only for Conjugate gradient method (**ifmd** == 1 & **ioptmze_cell** == 0).

#### 5.16.10 (stabilizer for quasi-Newton)
**gammamin** = [real]

Default: <br>
**gammamin** = 0.1

only for quasi-Newton method (**ifmd** == 1 & **ioptmze** == 2).

#### 5.16.11 (clear Hessian)
**ibfgsclear** = [0|1]

Default: <br>

**ibfgsclear** = 0

1: Clear Hessian after every ibfgsclear step
0: Hessian is not cleared

#### 5.16.12 (clear cell Hessian)
**iclearcellh** = [integer]

Default: <br>
**iclearcellh** = 0

Clear Hessian every **iclearcellh** step.  If **iclearcellh** == 0, the Hessian is not cleared.

#### 5.16.13 (hybrid optimization)
**lhybridopt** = [boolean] <br>
**nstep_hybrid** = [integer] <br>
**nstep\_hybrid\_cell** = [integer]

Default: <br>
**lhybridopt** = .TRUE. <br>
**nstep_hybrid** = 10 <br>
**nstep\_hybrid\_cell** = 10

**lhybridopt** = .TRUE. : Perform structural optimization first, then cell optimization. <br>
**lhybridopt** = .FALSE. : Perform cell optimization first, then structural optimization.

**nstep_hybrid** is the time step for structural optimization. <br>
**nstep\_hybrid\_cell** is the time step for cell optimization. <br>
These variables are only read for optimization calculations(**ifmd** == 1 & **ioptmze** >= 0 & **ioptmze_cell** >= 0 )

#### 5.16.14 (displacement)
**hmadisp** = [real]

Default: <br>
**hmadisp** = 0.01

Atomic displacement in [a.u.] for Harmonic analysis by finite difference.  This is only used for Harmonic-mode analysis (**ifmd** == 1 & **ioptmze** == 10).

#### 5.16.15 (order of differentiation)
**nhmaord** = [integer]

Default: <br>
**nhmaord** = 1

Order of differentiation only defined for Harmonic-mode analysis (**ifmd** == 1 & **ioptmze** == 10).

#### 5.16.16 (thermostat parameters)
**nnos** = [integer] <br>
**nresn** = [integer] <br>
**nyosh** = [integer] <br>
**tomega** = [real]

Default: <br>
**nnos** = 1
**nresn** = 1
**nyosh** = 1
**tomega** = 5500.0

[add references]  <br>
Only for NVT & NPT-MD (**ifmd** == 3, 4) <br>
Notes for (tomega): <br>
The frequency for the thermostat will be given by omega [a.u.] = 2*pi / tomega.  The tomega will be determined from the period of VAF.  If you are not aware of it, try a value,  100 * dtmd in [a.u.], as an initial guess.

#### 5.16.17 (atom type with thermostat)
**nathermo** = [integer]
**ithermo** = [integer]

Only used for NVT for each atom (**ifmd** == 5)

**nathermo** gives the total number of atom types to be given a thermostat, followed by that many lines of **ithermo**, which specifies the atom type.

Example:
```
   2  : (nathermo) # of atom types with thermostat
   1  : (ithermo) atom type
   4  : (ithermo) atom type
```

#### 5.16.18 (atom thermostat parameters)
**nnos_atom** = [integer] <br>
**nresn_atom** = [integer] <br>
**nyosh_atom** = [integer] <br>
**tomega_atom** = [real]

Default: <br>
**nnos_atom** = 1 <br>
**nresn_atom** = 1 <br>
**nyosh_atom** = 1 <br>
**tomega_atom** = 5500.0

**nnos\_atom** gives the number of thermostats, **nresn\_atom** gives the MTS steps for heat bath, **nyosh\_atom** gives the Yoshida-Suzuki decomposition step, and **tomega_atom** a period of fluctuation in [a.u.].  This is only used for NVT for each atom (**ifmd** == 5).

#### 5.16.19 (pressure)
**hpext** = [real]

Default: <br>
**hpext** = 0.0

Defines the pressure in [GPa] for NPT-MD & MSST (**ifmd** == 4, 10) & **ioptmze_cell** >= 0.

#### 5.16.20 (barostat parameters)
**tbomega** = [real] <br>
**blkmod** = [real]

Default: <br>
**tbomega** = 5500.0 <br>
**blkmod** = 250.0

Barostat parameters for NPT (**ifmd** == 4).  **tbomega** gives the time scale for barostat in [a.u.], and **blkmod** gives the bulk modulus in [GPa].

#### 5.16.21 (restriction for MD cell)
**irstrct** = [0|1|2|3|4|5|6|10|11|12|13|14|15|16]

Default: <br>
**irstrct** = 0

| irstrct | Type of Restriction on cell                         |
| ------- | :---------------------------------:                 |
| **0**   | No restriction                                      |
| **1**   | **Cubic**: a=b=c; alpha=beta=gamma=90                   |
| **2**   | **orthorhombic**: a /= b/= c; alpha=beta=gamma=90       |
| **3**   | **Tetragonal**: a = b /= c; alpha=beta=gamma=90         |
| **4**   | **Monoclinic**: a /= b /= c; alpha=beta=90; gamma /= 90 |
| **5**   | **Hexagonal**: a = b /= c; alpha=beta=90; gamma=120     |
| **6**   | **Trigonal**: a=b=c; alpha = beta = gamma /= 90         |
| **10**  | Set by symmetry opearations (enabled automatically when lsymop=.true. & irstrct==0)|
Add 10 to the above values to enforce the lattice shape with symmetry opearation.  For example, set **irstrct** == 12 for an orthorhombic cell with symmetry opearations.

Set restrictions for geometry of MD cell in NPT-MD simulations.

#### 5.16.22 (sub-restriction for MD cell)
**irstrct_sub** = [integer]

Default: <br>
**irstrct_sub** = 1

| irstrct | irstrct_sub | Type of Sub-restriction on cell         |
| ------- | :---------: | :----------------------------:          |
| **3**   | **1**       |  a  = b /= c, alpha = beta = gamma = 90 |
|         | **2**       |  b  = c /= a, alpha = beta = gamma = 90 |
|         | **other**   |  c  = a /= b, alpha = beta = gamma = 90 |
| **4**   | **1**       |  a  /= b /= c, alpha = beta = 90; gamma /= 90 |
|         | **2**       |  a  /= b /= c, beta = gamma = 90; alpha /= 90 |
|         | **other**   |  a  /= b /= c, gamma = alpha = 90; beta /= 90 |
| **5**   | **1**       |  a  = b /= c, alpha = beta = 90; gamma = 120 |
|         | **2**       |  b  = c /= a, beta = gamma = 90; alpha = 120 |
|         | **other**   |  c  = a /= b, gamma = alpha = 90; beta = 120 |

Defines sub-options for restriction for MD cell in NPT-MD simulations.

#### 5.16.23 (MD cell edge restriction)
**lcell_rstrct(1)** = [boolean] <br>
**lcell_rstrct(2)** = [boolean] <br>
**lcell_rstrct(3)** = [boolean] 

.TRUE. = fix the length of L[1,2,3] <br>
.FALSE. = do not fix length

Determines whether or not to fix the lengths of MD cell edges.

#### 5.16.24 (shock wave velocity)
**shockspeed** = [real]
**nshockv(1:3)** = [0|1, 0|1, 0|1]

Default: <br>
**shockspeed** = 2000.0 <br>
**nshockv(1:3)** = (/1,0,0/)

These variables are only read for MSST (**ifmd** == 10).  They define the shockspeed and the shock direction.  For example:

**nshockv(1:3)** = 1 0 0:  L1 direction <br>
**nshockv(1:3)** = 0 1 0:  L2 direction <br>
**nshockv(1:3)** = 0 0 1:  L3 direction <br>
**nshockv(1:3)** = 1 1 0:  L1 + L2 direction

#### 5.16.25 (clear barostat velocity)
**lmsstscale** = [boolean] <br>
**msstscnum** = [integer] <br>
**msstscstp** = [integer]

Default: <br>
**lmsstscale** = .FALSE. <br>
**msstscnum** = 20 <br>
**msstscstp** = 500

These variables are only read for MSST (**ifmd** == 10).  

#### 5.16.26 (initial barostat velocity)
**vxmsst** = [real]

Default: <br>
**vxmsst** = 0.0

Defines the initial barostat velocity.  Only read for MSST (**ifmd** == 10). 

#### 5.16.27 (output data)
**ioskip** = [integer]
**locoor, ioskipcoor** = [boolean, integer]

Defines what data to output for MD-nodes.

#### 5.16.28 (charge estimation)
**ichest** = [integer]

Default: <br>
**ichest** = 3

Number of previous steps

#### 5.16.29 (wavefunction estimation)
**ihest** = [integer]

Number of previous steps

#### 5.16.30 (ASPC charge estimation)
**laspc_chg** = [boolean]

.TRUE. = ASPC predictor <br>
.FALSE. = usual extrapolation

#### 5.16.31 (ASPC wavefunction estimation)
**laspc_wv** = [boolean]

.TRUE. = ASPC predictor <br>
.FALSE. = subspace alignment

#### 5.16.32 (ASPC acceleration)
**laspc** = [boolean]

.TRUE. = do ASPC 

#### 5.16.33 (Extended Lagrangian scheme)
**lxlbomd** = [boolean]

.TRUE. : time reversible XL-BOMD

#### 5.16.34 (previous values)
**kxlbomd** = [integer]

Number of previous values.  If 0, kxlbomd = ihest.

#### 5.16.35 (ASPC parameter for corrector)
**aspc_corr** = [real]

Default: <br>
**aspc_corr** = 0.0

#### 5.16.36 (ASPC mixing charge)  
**aslh_aspc, bslh_aspc** = [real, real]

#### 5.16.37 (ASPC iterations)
**iscfmx_aspc** = [integer]

Gives the number of iterations.

#### 5.16.38 (ASPC tolerance)
**tolres_aspc** = [real]

Gives the tolerance for the residual.

#### 5.16.39 (ASPC correction skip step)
**nskp_aspc** = [integer]

#### 5.16.40 (planewave expansion)
**pwscale** = [real]

only for NPT-MD (ifmd == 4 )

#### 5.16.41 (statistical calculation)
**lstat** = [boolean]

#### 5.16.42 (tolerance)
**tol_energy** = [real]
**tol_force** = [real]

tolerance for CG optimization (ifmd == 1 ).  **tol_energy** given in energy/atom in [a.u.].  **tol\_force** gives maximum force  in [a.u.].


### 5.17 save data
#### 5.17.1 (On/Off)
**lsave** = [boolean]

Default: <br>
**lsave** = .TRUE.

Whether or not to save data...which data?

#### 5.17.2 (data type)
**lsreal8** = [boolean]

Default: <br>
**lsreal8** = .TRUE.

.TRUE. : Saves data as real*8 <br>
.FALSE. : Saves data as integer*2

Determines the data type for saved data


### 5.18 soft walls
#### 5.18.1 (walls)
**nwalls** = [integer] <br>
**wallp** = [real, real, real] <br>
**wallv** = [real, real, real] <br>
**wallf** = [real] <br>
**nwallp** = [integer]

Default: <br>
**nwalls** = 0

**nwalls** gives the number of soft walls to be inserted into the system, followed by that many sets of **wallp**, **wallv**, **wallf**, **nwallp**, which define each soft wall.  **wallp** gives any point on the wall in [a.u.].  **wallv** gives the normal vector to the wall, **wallf** gives the strength of repulsive potential, and **nwallp** gives the power of the repulsive potential.

Example: This defines two soft walls
```
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
```

### 5.19 gravitational field
#### 5.19.1 (On/Off)
**lgravi** = [boolean]

Default: <br>
**lgravi** = .FALSE.

Determines whether to apply a gravitational field.

#### 5.19.2 (magnitude)
**gravmag** = [real]

Default:
**gravmag** = 0.0

Determines the stength of the gravitational field in units of g = 9.8 [m/s^2]

#### 5.19.3 (direction)
**gravdir** = [real, real, real]

Default: <br>
**gravdir** = [0.0, 0.0, 0.0]

Determines the direction of the gravitational field.

### 5.20 Mulliken analysis
#### 5.20.21 (On/Off)
**lmulken** = [boolean]

#### 5.20.22 (skip step)
**nskpmulk** = [integer]

#### 5.20.23 (decompose overlap population)
**ldecmpovp** = [boolean]

#### 5.20.24 (weights associated with each atom)
**lspdatom** = [boolean]

#### 5.20.25 (bond charge)
**nbondc** = [integer]
**ibondc, ilbondc, jbondc, jlbondc** = [integer integer integer integer]


### 5.21 Spherical Harmonics expansion
#### 5.21.1 (On/Off)
**lsphexp** = [boolean]

Default: <br>
**lsphexp** = .FALSE.

Determines whether or not to use spherical harmonics expansion.  Note this is only read for molecular dynamics simulations.

#### 5.21.2 (skip step)
**nskpsphexp** = [integer]

Default: <br>
**nskpsphexp** = 5

#### 5.21.3 (default radius)
**radsphexp** = [real]

Default: <br>
**radsphexp** = 0.0

Gives the default radius in a.u.. If **radsphexp** = 0, slater_r are used.


### 5.22 EDA
#### 5.22.1 (On/Off)  
**leda** = [boolean]

Default: <br>
**leda** = .FALSE.

Determines whether or not to perform energy density analysis (EDA).  Note this can only be used in molecular dynamics.

#### 5.22.2 (skip step)
**nskpeda** = [intger]

Default: <br>
**nskpeda** = 5

Will peform EDA every **nskpeda** steps.

#### 5.22.3 (radius for grid EDA)
**radgeda** = [real]

Default: <br>
**radgeda** = 10.0

Gives the radius in a.u. for EDA calculation.

### 5.23 Wannier function
#### 5.23.1 (On/Off)
**lwannier** = [boolean]

Default: <br>
**lwannier** = .FALSE.

Determines whether or not to compute maximully localized Wannier functions.

#### 5.23.2 (band index)
**iwbnd1, iwbnd2** = [integer] [integer]

Default: <br>
**iwbnd1, iwbnd2** = 0 0

Gives the range of band indices to be transformed by Unitary matrix.  If **iwbnd1, iwbnd2** = 0 0, then all occupied bands will be transformed.

#### 5.23.3 (iteration)
**iterwan** = [integer]

Default: <br>
**iterwan** = 200

Gives the maximum number of iterations.

#### 5.23.4 (tolerance)
**tolwan** = [real]

Default: <br>
**tolwan** = 1.d-06

Gives the tolerance???

#### 5.23.5 (dump orbitals)
**iwstt1, iwstt2** = [integer] [integer]

Default: 
**iwstt1, iwstt2** = 1 0

Gives the orbital indices to dump data for.  If **iwstt1, iwstt2** = 0 0, then all orbitals are dumped. If **iwstt1, iwstt2** = 1 0, then no orbitals are dumped.

#### 5.23.6 (orbitals around atoms)
**natwan** = [integer]
**natwno, nwanaa** = [integer] [integer]

Default: <br>
**natwan** = 0

**natwan** gives the number of focused atoms, followed by that many lines of **natwno, nwanaa** pairs, where **natwno** is the atom number and **nwanaa** is the wavefunction number.

#### 5.23.7 (skip step)
**nskpwann** = [integer]

Default: <br>
**nskpwann** = 5

Gives the number of steps to skip.

#### 5.23.8 (output unitary matrix)
**loutuni** = [boolean]

Default: <br>
**loutuni** = .FALSE.

Determines whether or not to output the unitary matrix.


### 5.24 Conductivity
#### 5.24.1 (On/Off)
**lconduct** = [boolean]

Default: <br>
**lconduct** = .FALSE.

.TRUE.  : Execute the conductivity calculation <br>
.FASLE. : Do not execute the conductivity calculation

Determines whether to run the conductivity calculation. At least one of **ldcconduct** or **lacconduct** is set to .TRUE.


#### 5.24.2 (skip step)
**nskpconduct** = [integer]

Default: <br>
**nskpconduct** = 5

**nskpconduct** is skip step for the conductivity calculation.


#### 5.24.3 (unit of energy)
(ry) or (hr) or (ev)

Default: <br>
(ry)


#### 5.24.4 (width of energy)
**wgconduct** = [real]

Default: <br>
**wgconduct** = 0.1

The gaussian broadening is used, alternative to the Delta function.


#### 5.24.5 (temperature)
**tempconduct** = [real]

Default: <br>
**tempconduct** = 300.0

Temperature in [K]


#### 5.24.6 (DC current density)
**lconduct** = [boolean]

Default: <br>
**ldcconduct** = .FALSE.

.TRUE.  : <br>
.FASLE. :

???


#### 5.24.7 (electric field)
**efconduct** = (/ [real] [real] [real] /)

Default: <br>
**efconduct(1:3)** = (/ 0.0 0.0 0.0/)

**efconduct(1:3)** is the constant electric-field vector. The unit of **efconduct(1:3)** is atomic unit.


#### 5.24.8 (current states)
(**icband1**, **icband2**) = [integer] [integer]

Default: <br>
(**icband1**, **icband2**) = (/ 0, 0 /)

Range of band indices for total current.  If both **icband1** and **icband2** are set to zero, the range is determined automatically.


#### 5.24.9 (hole states)
(**ichband1**, **ichband**) = [integer] [integer]

Default: <br>
(**ichband1**, **ichband2**) = (/ 0, 0 /)

Range of band indices for hole current.  If both **ichband1** and **ichband2** are set to zero, the range is determined automatically.


#### 5.24.10 (electron states)
(**iceband1**, **iceband**) = [integer] [integer]

Default: <br>
(**iceband1**, **iceband2**) = (/ 0, 0 /)

Tange of band indices for electron current. If both **iceband1** and **iceband2** are set to zero, the range is determined automatically.
Note:

icband1  <= ichband1 <= ichband2 <= iceband1 <= iceband2 <= icband2


#### 5.24.11 (output momentum matrix)
**lacconduct** = [boolean]

Default: <br>
**lacconduct** = .FALSE.

.TRUE.  : <br>
.FASLE. :

???
to caluculate optical conductivity


#### 5.24.12 (unit of energy)
(ry) or (hr) or (ev)

Default: <br>
(ry)


#### 5.24.13 (frequency range)
**freqacmx** = [real]

Default: <br>
**freqacmx** = 1.0

???
occupied:   e_HOMO - **freqacmx** ~ e_HOMO
unoccupied: e_LUMO ~ e_LUMO + **freqacmx**


#### 5.24.14 (matrix states)
(**iceband1**, **iceband**) = [integer] [integer]

Default: <br>
(**immtband1**, **immtband2**) = (/ 0, 0 /)

???
range of band indices for mometum matrix
Both of **immtband1** and **immtband2** are zero, the range is determined automatically.


### 5.25 stress calculation  ------  only for bulk calculations
#### 5.25.1 (On/Off)
**lstress** = [boolean]

Default: <br>
**lstress** = .TRUE.

.TRUE.  : Estimate the stress <br>
.FASLE. : Do not estimate the stress

Determines whether to estimate the stress.
#### 5.25.2 (skip step)
**nskip_stress** = [integer]

Default: <br>
**nskip_stress** = 5

**nskip_stress** is skip step to estimate the stress.


### 5.26 dump charge density
#### 5.26.1 (On/Off)
**ldpchg** = [boolean]

Default: <br>
**ldpchg** = .FALSE.

.TRUE.  : Output charge density <br>
.FASLE. : Do not output

Determines whether to output the charge density.


#### 5.26.2 (skip step)
**nskip_dpchg** = [integer]

Default: <br>
**nskip_dpchg** = 5

**nskip_dpchg** is skip step to output the charge density.


#### 5.26.3 (output area)
**x\_min** **x\_max** = [real] [real] <br>
**y\_min** **y\_max** = [real] [real] <br>
**z\_min** **z\_max** = [real] [real] <br>

Default: <br>
**x\_min** **x\_max** = 1.0 0.0 <br>
**y\_min** **y\_max** = 1.0 0.0 <br>
**z\_min** **z\_max** = 1.0 0.0

These parameter specify the area to output the charge density. They are reduced coordinates. If **[x|y|z]\_min** > **[x|y|z]\_max**, output are is the whole space.

### 5.27 dump wavefunctions
#### 5.27.1 (On/Off)
**ldpwav** = [boolean]

Default: <br>
**ldpwav** = .FALSE.

.TRUE.  : Output wavefunction <br>
.FASLE. : Do not output

Determines whether to output the wavefunction.


#### 5.27.2 (bands)
(**ibstt1**, **ibstt2**) = [integer] [integer]

Default: <br>
(**ibstt1**, **ibstt2**) = (/ 0, 0 /)

The wavefunction will be output from **ibstt1** to **ibstt2**. **ibstt1** and **ibstt2** are the index of bands. If both **ibstt1** and **ibstt2** are zero, all bands will be output.

#### 5.27.3 (skip step)
**nskip_dpwav** = [integer]

Default: <br>
**nskip_dpwav** = 5

**nskip_dpwav** is skip step to output the wavefunction.


#### 5.27.4 (output area)
**x_min** **x_max** = [real] [real] <br>
**y_min** **y_max** = [real] [real] <br>
**z_min** **z_max** = [real] [real] <br>

Default: <br>
**x_min** **x_max** = 1.0 0.0 <br>
**y_min** **y_max** = 1.0 0.0 <br>
**z_min** **z_max** = 1.0 0.0 

These parameter specify the area to output the wavefunction. They are reduced coordinates. If **[x|y|z]\_min** > **[x|y|z]\_max**, output are is the whole space.

### 5.28 dump local potential
#### 5.28.1 (On/Off)
**ldppot** = [boolean]

Default: <br>
**ldppot** = .FALSE.

.TRUE.  : Output local potential <br>
.FASLE. : Do not output

Determines whether to output the local potential.


#### 5.28.2 (average plane)
**nav_dppot** = [1|2|3]

Default: <br>
**nav_dppot** = 1

| **nav_dppot**  | Plane |
| -------------- | :---: |
| **1**          | xy    |
| **2**          | yz    |
| **3**          | zx    |

The plane-averaged local potential will be output.

#### 5.28.3 (skip step)
**nskip_dppot** = [integer]

Default: <br>
**nskip_dppot** = 5

**nskip_dppot** is skip step to output the local potential.


#### 5.28.4 (output area)
**x_min** **x_max** = [real] [real] <br>
**y_min** **y_max** = [real] [real] <br>
**z_min** **z_max** = [real] [real] <br>

Default: <br> 
**x_min** **x_max** = 1.0 0.0 <br>
**y_min** **y_max** = 1.0 0.0 <br>
**z_min** **z_max** = 1.0 0.0 <br>

These parameter specify the area to output the local potential. They are reduced coordinates. If **[x|y|z]\_min** > **[x|y|z]\_max**, output are is the whole space.

### 5.29 supercell
#### 5.29.1 (unit of length)
(bohr) or (ang)

Default: <br>
(bohr)

Determines the unit of supercell (bohr or angstrom). And the supercell is specified by either cell vector or 'lengths and angles'. You need comment out the subsection you do not specify. Note that for cluster calculations, the supercell has to be ORTHORHOMBIC.


#### 5.29.2 (cell vector)
L1 = [real] [real] [real] <br>
L2 = [real] [real] [real] <br>
L3 = [real] [real] [real]

Example:
```
 10.68   0.0    0.0        : super cell vector L1
  0.0   10.68   0.0        : super cell vector L2
  0.0    0.0   10.68       : super cell vector L3
```


#### 5.29.3 (lengths & angles)

Example:
```
 10.68,  10.68,  10.68  :  lengths of cell vectors
 90.00,  90.00,  90.00  :  angles between cell vec. in [deg.]
```


### 5.30 vacuum
This section is mainly used for the double grid method.
#### 5.30.1 (unit of length)
(bohr) or (ang)

Defaubr>
(bohr)

Determines the unit of vacuum (bohr or angstrom).


#### 5.30.2 (On/Off)
**lvacuum(1)**, **vacuum(1)** [boolean] [real] <br>
**lvacuum(2)**, **vacuum(2)** [boolean] [real] <br>
**lvacuum(3)**, **vacuum(3)** [boolean] [real] <br>

Default: <br>
**lvacuum(1:3)** = .FALSE. <br>
**vacuum(1:3)** = 0.0

If **lvacuum([1|2|3])** is .TRUE., the vaccum is **vacuum([1|2|3])**< [x|y|z] < [L1|L2|L3].


#### 5.30.3 (parameter in error function)
**alpha_ldouble_grid_recip** = [real]

Default: <br>
**alpha_ldouble_grid_recip** = 0.0

If **alpha_ldouble_grid_recip** is zero, **alpha_ldouble_grid_recip** is set automatically. **alpha_ldouble_grid_recip** should be about 5.0/(the length of cell).


### 5.31 spherical region
#### 5.31.1 (On/Off)
**lsphere** = [boolean]

Default: <br>
**lsphere** = .FALSE.

.TRUE.  : Use spherical region in cluster calculations <br>
.FASLE. : Do not use


### 5.32 planewaves
#### 5.32.1 (unit of cutoff energy)
(ry) or (hr) or (ev)

Default: <br>
(ry)

Determines the unit of cutoff energy.


#### 5.32.2 (for wavefuctions)
**ecut** = [real]

**ecut** is the cutoff energy to expand wavefunctions.


#### 5.32.3 (for electron density)
**ecutdens** = [real]

Default: <br>
**ecutdens** = 0.0  <br>

**ecutdens** is the cutoff energy to expand electron density. **ecutdens** must be greater than **ecutsoft**.  if .not.lvand, **ecutdens** must be equal to **ecutsoft**.


#### 5.32.4 (for soft part of density)
**ecutsoft** = [real]

Default: <br>
**ecutorth** = 0.0  <br>

**ecutsoft** is the cutoff energy to expand soft part of electron density.  It must be greater than ecut, and smaller than ecut*4.


#### 5.32.5 (for orthogonalization)
**ecutorth** = [real]

Default: <br>
**ecutorth** = 0.0 

**ecutorth** is the cutoff energy for orthogonalization in CG iteration.


#### 5.32.6 (for charge density mixing)
**ecutc** = [real]

Default: <br>
**ecutc** = 0.0

**ecutc** is the cutoff energy for charge density mixing. 



### 5.33 double-grid method
[add ref]
#### 5.33.1 (order of cardinal B spline)
**iosp** = [integer]

Default: <br />
*iosp* = 3. This number must be an odd integer

#### 5.33.2 (unit of cutoff energy)
(ry) or (hr) or (ev)

Default: <br/>
(ry)

#### 5.33.3 (for non-periodic direction)
**ecutlong** = [real]

Default: <br/>
**ecutlong** = 0.0

Cutoff energies for the long-range Coulomb interaction

#### 5.33.3 (for periodic direction)
**ecutlong_small** = [real]

Default:<br />
**ecutlong_small** = 0.0

Cutoff energies for the long-range Coulomb interaction in the periodic direction (for wire and surface geometries)

### 5.34 electronic bands
#### 5.34.1 (occupied bands)
**noband** = [integer]

Default:<br />
**noband** = 0

The number of occupied bands should at least be sufficient to accommodate all the valence electrons in the simulation.

#### 5.34.2 (empty bands)
**neband** = [integer]

Default:<br />
**neband** = 0

A good rule of thumb is to set neband equal to 10% of occupied bands.

#### 5.34.3 (broadening)
**lfermi** = [integer] flag <br>
**tfermi** = [real]

Default:<br />
**lfermi** = 3<br />
**tfermi** = 2000.0

lfermi = 1 for Non-metallic systems, lfermi = 2 for Fermi smearing, lfermi = 3 for Gaussian smearing and lfermi > 3 for Methfessel-Paxton smearing scheme of order lfermi-3

#### 5.34.4 (charge number of ion)
**ncion** = [integer]

Default:<br />
**ncion** = 0. For calculation of electronically non-neutral systems. The number of electrons is given by `nel-ncion`. Charged-periodic calculations are perfomed with the uniform background.




### 5.35 Symmetry operations
#### 5.35.1 (On/Off)
**lsymop** = [boolean]

Default:<br />
**lsymop** = .FALSE.


The flag **lsymop** defines if symmetry operations are checked.


#### 5.35.2 (unit cell)
**nsymm** = [integer][integer][integer]<br />

Default:<br />
**nsymm** = 1 1 1

**nsymm** defines the number of unit cells in the x, y and z directions.


The transformation matrix is specified through either (transformation matrix) or (inverse transformation matrix) as follows.<br />
`trsfrm(j,i) <-- trsfrm(j,i)/nsymm(i)` or `trsfri(i,j) <-- trsfri(i,j)*nsymm(i)`


#### 5.35.3 (transformation matrix)
**trsfrm(1,i)** = [real][real][real]<br />
**trsfrm(2,i)** = [real][real][real]<br />
**trsfrm(3,i)** = [real][real][real]<br />

Default:<br />
**trsfrm** = No defaults


#### 5.35.3 (inverse transformation matrix)
**trsfri(1,i)** = [real][real][real]<br />
**trsfri(2,i)** = [real][real][real]<br />
**trsfri(3,i)** = [real][real][real]<br />

Default:<br />
**trsfri** = No defaults


#### 5.35.4 (origin of symmetry op.)
**symm** = [real][real][real]

Default: <br />
**symm** = 0.0 0.0 0.0


#### 5.35.5 (symmetry operations)
**nmsyop** = [integer]
`???`

Default: <br />
**nmsyop** = 1


#### 5.35.6 (denominator of trans. vector)
**denomt** = [real]

Default: <br />
**denomt** = No default


#### 5.35.7 (inversion operations)
**linvsn** = [boolean]

Default: <br />
**linvsn** = .FALSE.




#### 5.35.8 (radius to identify atoms)
**radsym** = [real]

Default: <br />
**radsym** = 1.0

Radius to identify symmetrically equivalent atoms


#### 5.35.9 (for charge density)
**lsymcd** = [boolean]

Default: <br />
**lsymcd** = .FALSE.


#### 5.35.10 (for spin density)
**lsymspin** = [boolean]

Default: <br />
**lsymspin** = No default








### 5.36 Atoms
#### 5.36.1 (species)
**ntype** = [integer]

Default:<br />
**ntype** = 0

Number of atom types. The subsequent sections must be duplicated `ntype` times.

#### 5.36.2 (atomic number)
**zatom** = [integer]

Default:<br />
**zatom** = No default

Atomic number of the current atom type

#### 5.36.3 (pseudopotential)
kbpp or uspp or vand or local

Default: <br/>
No default

Type of pseudopotential defined for the current atom type.

#### 5.36.4 (nonlocal potential)
**lking** = [boolean]<br />
**rking** = [real]<br />
**gkgmax** = [real]<br />
**gkgexct** = [real]


Default:<br />
**lking** = No default<br />
**rking** = No default<br />
**gkgmax** = No default<br />
**gkgexct**  = No default

These flags specify if nonlocal potntial is computed in real space or reciprocal space. We follow the formalism by [<a href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.44.13063"> R. D. King-smith et al. </a>]. lking=.false. allows calculation to proceed in reciporocal space. In such case, next three parameter are read but not use. If lking=.true., caculation will proceed in real space. In such case, rking, gkgmax and gkgexct are used. rking corresponds to real space cutoff radius. gkgmax and gkgexct correspond to accuracy of non-local potentials. For more detail, please follow paper by R. D. King-smith et al. Accuracy of real space must be tested if input parameters are not provided. 

Recommendation:<br />
**lking** = false for small system and true for large system <br />
**rking** = 1.5- 2<br />
**gkgmax** = 1.25<br />
**gkgexct**  = 0.98

#### 5.36.5 (local potential)
**llking** = [boolean]<br />
**rlking** = [real]<br />
**glkgmax** = [real]<br />
**glkgexct** = [real]

Default:<br />
**llking** = No default<br />
**rlking** = No default<br />
**glkgmax** = No default<br />
**glkgexct**  = No default

These flags specify if local potntial is computed in real space or reciprocal space. It is recommended to use reciprocal space reciprocal space for the calculation. In reciprocal space rlking, glkgmax and glkgexct are not used but values are read. 

Recommendation:<br />
**llking** = false<br />


#### 5.36.6 (partial core correction)
**lpcc** = [boolean]<br />
**rpcc** = [real]<br />
**lpking** = [boolean]<br />
**rpking** = [real]<br />
**gpkgmax** = [real]<br />
**gpkgexct** = [real]

Default:<br />
**lpcc** = No default<br />
**rpcc** = No default<br />
**lpking** = No default<br />
**rpking** = No default<br />
**gpkgmax** = No default<br />
**gpkgexct** = No default

Psuedopotential formalism treat core and valence electron seperately. Partial core correction includes the effect of core charge as a perturbation which increase the transferrability of the code. It is recommended to use partial core correction. Please see paper by <a href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.26.1738"> Louie et al. </a>

#### 5.36.7 (compensation charge cutoff)
**frac_rcomp_it** = [real]

Default:<br />
**frac_rcomp_it** = 1.3

Compensation charge must be used in PAW method to match all electron wave function. A recommended value is between 1.2 - 1.3. For more detail please see paper by <a href="https://journals.aps.org/prb/pdf/10.1103/PhysRevB.59.1758"> Kresse et al. </a>
#### 5.36.8 (the number of atoms)
**nhk** = [integer]

Default:<br />
**nhk** = 0

Number of atoms of the current type. If `nhk` is zero, then atomic configuration is read from the positions file

#### 5.36.9 (unit of length)
(bohr) or (ang)

Default:<br />
No default


#### 5.36.10 (position file)
**fname** = [String]<br />
**scaledflag*** = [integer]<br />
**keyword*** = [integer]<br />

Default:<br />
**fname** = No default<br />
**scaledflag*** = No default<br />
**keyword*** = No default

Filename to read atomic positions from. Filenames are given with reference to the current working directory. **scaledflag** determines if the aomtic positions in the coordinate file are in real or scaled units. Atoms of this type are denoted by the integer **keyword** in the corrdinates file.

The structure of 'CONFIGURATION_FILE' is as follows.
```
   ntotal
   n_1      x_1      y_1      z_1
   n_2      x_2      y_2      z_2
   ...      ...      ...      ...
   n_ntotal x_ntotal y_ntotal z_ntotal

```

The coordinates for (n_i == keyword) are selected.

#### 5.36.11 (positions)
**nunitcells** = [integer][integer][integer]<br />
**scaledflag*** = [integer]<br />
`nhk` lines of **positions** = [real][real][real]

Default: <br />
**nunitcells** = No defaults <br />
**scaledflag*** = No defaults

**nunitcells** gives the number of unit cells along the x, y and z directions. The flag **scaledflag** defines if the atomic positions are defined in real or scaled coordinates (**scaledflag** = 1 is scaled coordinates, **scaledflag** = 2 is real coordinates).


If '1 (scaled)' is selected, the scaled coordinates ( 0.0 < x,y,z < 1.0 ) have to be given above. <br />
If '2 (real)' is selected, the real coordinates ( in a.u. ) have to be given like below.

```
 2                            : 1:scaled, 2:real coordinates
  0.00d0  0.00d0  0.00d0      : If '2 (real)' is selected,
 10.34d0 10.34d0 10.34d0      :  the real coordinates ( in a.u. )
 10.34d0  0.00d0 10.34d0      :  have to be given.
0.00d0 10.34d0 10.34d0
```

#### 5.36.12 (fix positions)
**lfixion** = [boolean]

Default: <br />
**lfixion** = .FALSE.

If .TRUE., atoms of current type are held fixed by setting atomic velocities to zero.


#### 5.36.13 (fix displacement)
**ldispion** = [boolean]<br />
**dvector** = [real] [real] [real]

Default: <br />
**ldispion** = .FALSE. <br />
**dvector** = No default

If .TRUE., positions of atoms of current type are updated by adding **dvector** to the current position.
`???`

#### 5.36.14 (fictitious mass)
**ficmass** = [real]

Default: <br />
**ficmass** = 1.0

If set, the mass of atoms of this type are scaled by **ficmass**


#### 5.36.15 (velocity file)
**fname** = [String]<br />
**keyword*** = [integer]<br />

Default:<br />
**fname** = No default<br />
**keyword*** = No default

**fname** is the filename to read atomic velocities from. Filename is given with reference to the current working directory. Atoms of this type are denoted by the integer **keyword** in the velocities file.

The structure of 'VELOCITY_FILE' is as follows.

```
   ntotal
   vmax                              : maximum velocity in (a.u.)
   n_1      vx_1      vy_1      vz_1 : an integer & scaled velocity
   n_2      vx_2      vy_2      vz_2
   ...      ...       ...       ...
   n_ntotal vx_ntotal vy_ntotal vz_ntotal
```

The velocities for (n_i == keyword) are selected.



#### 5.36.16 (velocities)
**nunitcells** = [integer][integer][integer]<br />
**scalingfactor*** = [real]<br />
`nhk` lines of scaled **velocities** = [real][real][real] . Ensure that -1.0 < **velocities** < 1.0

Default: <br />
**nunitcells** = No defaults <br />
**scalingfactor*** = No defaults

**nunitcells** gives the number of unit cells along the x, y and z directions. Atomic velocities are initialized by multiplying the vector **velocities** by the **scalingfactor**.


#### 5.36.17 (Spherical Harmonics expansion)
**rsphexp** = [real]

Default: <br />
**rsphexp** = No default

Radius for *Ylm* expansion


#### 5.36.18 (EDA)
**radeda** = [real]

Default: <br />
**radeda** = No default

Radius for EDA


#### 5.36.19 (atomic charge)
**rintchg** = [real]

Default: <br />
**rintchg** = No default

Radius of integration to calculate atomic charge


#### 5.36.20 (DFT+U)
**lplusU_at** = [boolean]<br />
**plusU_l** = [integer]<br />
**plusU_U** = [real]<br />
**plusU_J** = [real]

Default: <br />
**lplusU_at** = .FALSE.<br />
**plusU_l** = 2 <br />
**plusU_U** = 0.0<br />
**plusU_J** = 0.0

The flag **lplusU_at** determines if the DFT+U method is used. **plusU_l** defines the angular momentum quantum number of the orbital of interest. **plusU_U** and **plusU_J** define the Hubbard parameter U and the screened exchange energy in electron-volts.


#### 5.36.21 (DFT+C)
**lplusC_at** = [boolean]<br />
**plusC_r_e_s** = [real] [real]<br />
**plusC_r_e_p** = [real] [real]<br />
**plusC_r_e_d** = [real] [real]<br />
**plusC_r_e_f** = [real] [real]

Default: <br />
**lplusC_at** = .FALSE.<br />
**plusC_r_e_s** = 0.0 0.0<br />
**plusC_r_e_p** = 0.0 0.0<br />
**plusC_r_e_d** = 0.0 0.0<br />
**plusC_r_e_f** = 0.0 0.0

The flag **lplusC_at** determines if the DFT+C method is used. **plusC_r_e** defines the `plusC_r` and `plusC_e` vales for the s, p, d and f subshells.

The form of correction

`modified V_nl(r) = V_nl(r) + plusC_e * sin(a*r)/(a*r), a = pi/plusC_r`

where `plusC_r` is cutoff radius in [a.u.] and `plusC_e` is energy shift in [Ryd.]

### 5.37 Constraint conditions

#### 5.37.1 (unit of length)
(bohr) or (ang)

Default: <br/>
(ang)

#### 5.37.1 (bond length constraint)
**ncbonds** = [integer]<br />

Followed by **ncbonds** lines of constraints.

**ncbatm1** = [integer]<br />
**ncbatm2** = [integer]<br />
**cblength** = [real]

For each line, the bond length between **ncbatm1** and **ncbatm2** is held fixed at **cblength** Angstrom.

Default:<br />
**ncbonds** = 0<br />
**ncbatm1** = No default<br />
**ncbatm2** = No default<br />
**cblength** = No default


### 5.38 Virtual Molecular Dynamics


### 5.39 Cholesky decomposition

#### 5.39.1 (the number of nodes)
**node_c** = [integer]

Default:<br />
**node_c** = 4



### 5.40 Eigenvalue problem

#### 5.40.1 (the number of nodes)
**node_r** = [integer]

Default:<br />
**node_r** = 4


### 5.41 Work Array
#### 5.41.1 (size of array)
**dblock** = [real]

Default:<br />
**dblock** = 1.d+07

### 5.42 Table Dimension
#### 5.42.1 (local potential)
**mx1loc** = [integer]

Default:<br />
**mx1loc** = 2**13

#### 5.42.2 (others)
**mx1** = [integer]

Default:<br />
**mx1** = 2**13

<br>

------------------------------------------

Go to **[Table of Contents](#ht)**, **[Introduction](#h0)**,  **[1.Prerequisites](#h1)**,  **[2.Installation](#h2)**,  **[3.Program Structure](#h3)**, **[4.I/O Files](#h4)** ,**[5.Main input file](#h5)**,  **[6.Utility files](#h6)**

------------------------------------------
## <a name="h6"> 6. QXMD Simulation Quick-Start
The **Examples** directory holds examples of various QXMD simulations for a variety of systems.  Examining the input files in these examples may be the best way to begin to get an intuition for how to run your first QXMD simulation.  While all the input parameters are given in section 5, many of these parameters have default settings that need not be changed, and many of the input sections can be left out of the input file for many standard (simple) QXMD simulations. Thus, what follows can be viewed as a quick-start guide for getting your first QXMD simulation running as soon as possible.  To run more complex simulations, the reader should refer to input parameter settings in Section 5 and sample input files given in the **Examples** directory.

### 6.1 Minimal IN.PARAM file
Given below is an sample input file for the optimization of a water molecule which includes all the most important input parameters for running standard QXMD simulations.  Most users will find this input file sufficient for running optimization and quantum molecular dynamics jobs.  Of course, input parameters given in section 5 can be added to this input file (an example of how to add a new section of variables will be shown in the next section on running a non-adiabatic QMD simulation) to tailor a QXMD job to the user's specifications. 
``` 
---(Input-parameter file for QXMD)---
                              :
*parallel                     :
(QM-nodes)                    :
  1  1  1                     : (npx, npy, npz)
(MD-nodes)                    :
  1  1  1                     : (md_npx, md_npy, md_npz)
*end                          :
                              :
                              :
*start(on/off)                :
(how of it)                   :
 .false.                       : (lstart) .true. = restart
*end                          :
                              :
*PAW                          :
(how of it)                   :
 .true.                       : (lpaw)  .true.  = PAW method
                              :         .false. = pseudopotential method
(non-spherical symmetry)      :
 .false.                      : (lpaw_sym) .true. = full      symmetry
                              :           .false. = spherical symmetry
                              :
(onsite charge mixing)        :
    0.5d0                     : (paw_mix)
*end                          :
                              :
                              :
                              :
*approximation for Exc        :
(approximation)               :
       2                      : (jgga) = 1:LDA, 2:GGA(PBE), 3:GGA(RPBE),
                              :          4:GGA(revPBE), 5:vdW-DF, 6:vdW-DF2
                              :
(DFT-D)                       : an empirical correction for the vdW interaction
 .true.                       : (ldftd) .true. = on, .false. = off
                              :
*end                          :
                              :
                              :
                              :
*SCF iterations               :
(global iterations)           :
     100                      : (iscfmx) maximum No. of global iterations
(tolerances)                  :
      3.0d-08                 : (tolpot) tolerance for total energy
      5.0d-08                 : (tolres) tolerance for average residual
*end                          :
                              :
                              :
*molecular dynamics           :
(how of it)                   :
   1                          : (ifmd)
                              :   0:non, 1:optimization, 2:NVE, 3:NVT, 4:NPT
(time step)                   :
   50.0d0   1000               : (dtmd, nstop) time step, total step
(temperature)                 : only for real dynamics (NVE-, NVT-, NPT-MD )
   300.d0                     : (treq) temperature in [K]
(check temperature)           :
  .true.                     : (liscale) .true. = Do it !
   25                         : (iscnum)  number of temperature check
   20                         : (iscstp)  skip step
                              :
(optimization)                : only for structural optimization (ifmd == 1 )
    2                         : (ioptmze)
                              :  -1: do not optimize atomic coordinates
                              :   0: Conjugate gradient
                              :   1: Projected velocity Verlet
                              :   2: Quasi-Newton method with BFGS formula
                              :
(stabilizer for quasi-Newton) : only for quasi-Newton method (ifmd==1 & ioptmze==2)
   0.1d0                      : (gammamin)
                              :
(clear Hessian)               : only for quasi-Newton method (ifmd==1 & ioptmze==2)
    0                         : (ibfgsclear) clear Hessian every ibfgsclear step
                              :           if ibfgsclear == 0, Hessian is not cleared.
                              :
(atomic stress & energy)      : only for MD nodes
 .true.                       : (latomic) .true. = output atomic stress & energy
                              : Note that nskip_atomic = ioskip
                              : When lstress = false, atomic stress is not output.
                              :
(tolerance)                   : tolerance for CG optimization (ifmd == 1 )
  1.d-07                      : (tol_energy) energy/atom in [a.u.]
  5.d-04                      : (tol_force ) max. force  in [a.u.]
                              :
(output data)                 : only for MD nodes
    1                         : (ioskip)  skip step
 .true.    1                  : (locoor, ioskipcoor) scaled coordinates
 .false.    1                  : (lovelo, ioskipvelo) scaled velocities
 .false.    1                  : (loforc, ioskipforc) scaled forces
*end                          :
------------------------------:---------------------------------------
*supercell                    :
(unit of length)              :
(ang)                         : (bohr) or (ang)
                              :
(lengths & angles)            :
7.00d0, 7.00d0, 5.0d0         :  lengths of cell vectors
90.000,  90.000,  90.000  :  angles between cell vec. in [deg.]
*end                          :
                              :
*planewaves                   :
(unit of cutoff energy)       :
(ry)                          : (ry) or (hr) or (ev)
(for wavefunctions)           :
  30.0                        : (ecut)
(for electron density)        :
 250.0                        : (ecutdens)
(for soft part of density)    :
 70.0                         : (ecutsoft)
*end                          :
                              :
*electronic bands             :
(occupied bands)              :
      8                     : (noband)  No. of occupied bands
(empty bands)                 :
      2                     : (neband)  No. of empty bands
                              :           total No.= noband + neband
(broadening)                  :
       3    500.d0            : (lfermi) = 1:nonmetallic, 2:Fermi, 3:Gaussian,
                              :   lfermi(>3):Methfessel & Paxton, order=lfermi-3
                              : (tfermi) = electronic temp.(K), if metallic
*end                          :
                              :
------------------------------:---------------------------------------
*atoms                        :
(species)                     :
  2                           : (ntype) No. of atomic species
==============================:=======================================
(atomic number)               :
  8.0                         : (zatom)
(pseudopotential)             :
uspp                          : kbpp .or. uspp .or. vand
(nonlocal potential)          :
 .true. 1.5d0 1.25d0 0.8d0    : (lking) .true. = on, (rking, gkgmax, gkgexct)
                              : smoothing parameters
(local potential)             :
 .false. 1.5d0 1.15d0 0.8d0   : (llking) .true. = on, (rlking, glkgmax, glkgexct)
                              : smoothing parameters
(partial core correction)     :
 .true.   1.4d0               : (lpcc) .true. = on, (r_cut) in [a.u.]
 .true. 1.1d0 1.15d0 0.8d0    : (lpking) .true. = on, (rpking, gpkgmax, gpkgexct)
                              : smoothing parameters
                              :
(unit of length)              : only for positions
(ang)                         : (bohr) or (ang)
                              :
(position file)               : 
'control/IN.CONFIG'           :
 2                            : 1:scaled, 2:real coordinates
 1                            : (keyword)
(end)                         :
==============================:=======================================
(atomic number)               :
  1.0                         : (zatom)
(pseudopotential)             :
uspp                          : kbpp .or. uspp .or. vand
(nonlocal potential)          :
 .true. 1.5d0 1.25d0 0.8d0    : (lking) .true. = on, (rking, gkgmax, gkgexct)
                              : smoothing parameters
(local potential)             :
 .false. 1.5d0 1.15d0 0.8d0   : (llking) .true. = on, (rlking, glkgmax, glkgexct)
                              : smoothing parameters
(partial core correction)     :
 .false.   1.4d0              : (lpcc) .true. = on, (r_cut) in [a.u.]
 .true. 1.1d0 1.15d0 0.8d0    : (lpking) .true. = on, (rpking, gpkgmax, gpkgexct)
                              : smoothing parameters
                              :
(unit of length)              : only for positions
(ang)                         : (bohr) or (ang)
                              :
(position file)               : 
'control/IN.CONFIG'           :
 2                            : 1:scaled, 2:real coordinates
 2                            : (keyword)
(end)                         :
*end                          : end of setting *atoms
```

In fact, many of these parameters need not be defined for an optimization run, but are included so that only two things need to be changed to make this input file run a quantum molecular dynamics (QMD) simulation.  First and foremost, the '(how of it)' parameter (ifmd) in the \*molecular dynamics section must be changed.  Above, this parameter is set to 1, which will execute an optimization run, but setting this parameter to 3 (as shown below) will execute a QMD simulation in the NVT ensemble (this can also be set to '2' for the NVE ensemble or '4' for the NPT ensemble).
```
*molecular dynamics           :
(how of it)                   :
   3                          : (ifmd)
                              :   0:non, 1:optimization, 2:NVE, 3:NVT, 4:NPT
```

The other thing that must be changed is the input configuration file.  Usually, you will want to start a QMD simulation from the optimized geometry.  The optimized atomic positions are obtained after an optimization run using the 'pick_config.f90' utility file (see section 7.2).  One can either give the optimized geometry configuration file a new name and update this in the (position file) parameter in the *atoms section, or replace IN.CONFIG with the new geometry. Once the '(how of it)' parameter in the \*molecular dynamics section and the input configuration file is adjusted to the optimized geometry, the input file is complete for execution of a standard QMD simulation.

### 6.2 Adding parameters to IN.PARAM
Here, we give an example of how to add input parameters to the basic IN.PARAM file given above in section 6.1.  As described in the beginning of Section 5, the main input file is divided into sections corresponding to different groups of input parameters. Sections start with '\*SECTION_NAME' and end with '*end'. If you would like to add a parameter from a section that already exists in the input file, you main simply add in the parameter in any position within the section (i.e. after the *section name and before the corresponding *end). In order to add a new input parameter from a section that does not already exist in the input file, you must add that section (i.e. add the *SECTION_NAME and *end lines with your new parameter in between) along with the parameter.  You need not add all paramters in the newly added section.  Any parameters in the section that are not included will simply be set to their default values.  In this section, we show an example for altering the above IN.PARAM file in Section 6.1 to execute a non-adiabtic QMD (NAQMD) simulation.  In this case, we must add the *TDDFT-MD section (see section 5.3), as time-dependent DFT is required to run a NAQMD simulation.  NAQMD simulations involve exciting electrons from lower lying energy bands to higher lying ones, which is defined in the *TDDFT-MD section shown below.

``` 
---(Input-parameter file for QXMD)---
                              :
*parallel                     :
(QM-nodes)                    :
  1  1  1                     : (npx, npy, npz)
(MD-nodes)                    :
  1  1  1                     : (md_npx, md_npy, md_npz)
*end                          :
                              :
                              :
*start(on/off)                :
(how of it)                   :
 .false.                       : (lstart) .true. = restart
*end                          :
                              :
*PAW                          :
(how of it)                   :
 .true.                       : (lpaw)  .true.  = PAW method
                              :         .false. = pseudopotential method
(non-spherical symmetry)      :
 .false.                      : (lpaw_sym) .true. = full      symmetry
                              :           .false. = spherical symmetry
                              :
(onsite charge mixing)        :
    0.5d0                     : (paw_mix)
*end                          :
                              :
                              :
                              :
*approximation for Exc        :
(approximation)               :
       2                      : (jgga) = 1:LDA, 2:GGA(PBE), 3:GGA(RPBE),
                              :          4:GGA(revPBE), 5:vdW-DF, 6:vdW-DF2
                              :
(DFT-D)                       : an empirical correction for the vdW interaction
 .true.                       : (ldftd) .true. = on, .false. = off
                              :
*end                          :
                              :
                              :
                              :
*SCF iterations               :
(global iterations)           :
     100                      : (iscfmx) maximum No. of global iterations
(tolerances)                  :
      3.0d-08                 : (tolpot) tolerance for total energy
      5.0d-08                 : (tolres) tolerance for average residual
*end                          :
                              :
*TDDFT-MD                     :
(how of it)                   :
 .true.                       : (ltddft) .true. = execute MD based on TDDFT
                              :
(FSSH)                        :
 .true.                       : (ltddft_fssh) .true. = FSSH, .false. = Ehrenfest
                              :
(FSSH-switch)                 :
 .true.                       : (lfssh_switch) .true.  = switching aveilable
                              :                .false. = cccupations are fixed
(FSSH-ground-state-SCF)       :
 .true.                       : (lfssh_gsscf) .true.  = SCF with the ground  state
                              :               .false. = SCF with the excited state
(FSSH-mixing charge)          : only for lfssh_gsscf = .true.
   0.8d0  0.13d0              : (aslh_fssh, bslh_fssh)
			      :
(time step)                   :
   0.04d0                     : (dttddft) time step in [Hartree a.u.] in TDDFT-FSSH
                              :
(restart)                     :
 .false.                      : (ltddft_start) .true. = restart
                              :
(occupations)                 : for lrtddft = .false.
    2                         : (nocc_change) # of occupations to be changed
   4   0.0  0.0              : (numband, occ_new) band index, occupations(up&down)
   5   2.0  0.0              : (numband, occ_new) band index, occupations(up&down)

*end                          :
                              :
*molecular dynamics           :
(how of it)                   :
   3                          : (ifmd)
                              :   0:non, 1:optimization, 2:NVE, 3:NVT, 4:NPT
(time step)                   :
   50.0d0   1000               : (dtmd, nstop) time step, total step
(temperature)                 : only for real dynamics (NVE-, NVT-, NPT-MD )
   300.d0                     : (treq) temperature in [K]
(check temperature)           :
  .false.                     : (liscale) .true. = Do it !
   25                         : (iscnum)  number of temperature check
   20                         : (iscstp)  skip step
                              :
(optimization)                : only for structural optimization (ifmd == 1 )
    2                         : (ioptmze)
                              :  -1: do not optimize atomic coordinates
                              :   0: Conjugate gradient
                              :   1: Projected velocity Verlet
                              :   2: Quasi-Newton method with BFGS formula
                              :
(stabilizer for quasi-Newton) : only for quasi-Newton method (ifmd==1 & ioptmze==2)
   0.1d0                      : (gammamin)
                              :
(clear Hessian)               : only for quasi-Newton method (ifmd==1 & ioptmze==2)
    0                         : (ibfgsclear) clear Hessian every ibfgsclear step
                              :           if ibfgsclear == 0, Hessian is not cleared.
                              :
(atomic stress & energy)      : only for MD nodes
 .true.                       : (latomic) .true. = output atomic stress & energy
                              : Note that nskip_atomic = ioskip
                              : When lstress = false, atomic stress is not output.
                              :
(tolerance)                   : tolerance for CG optimization (ifmd == 1 )
  1.d-07                      : (tol_energy) energy/atom in [a.u.]
  5.d-04                      : (tol_force ) max. force  in [a.u.]
                              :
(output data)                 : only for MD nodes
    1                         : (ioskip)  skip step
 .true.    1                  : (locoor, ioskipcoor) scaled coordinates
 .false.    1                  : (lovelo, ioskipvelo) scaled velocities
 .false.    1                  : (loforc, ioskipforc) scaled forces
*end                          :
------------------------------:---------------------------------------
*supercell                    :
(unit of length)              :
(ang)                         : (bohr) or (ang)
                              :
(lengths & angles)            :
7.00d0, 7.00d0, 5.0d0         :  lengths of cell vectors
90.000,  90.000,  90.000  :  angles between cell vec. in [deg.]
*end                          :
                              :
*planewaves                   :
(unit of cutoff energy)       :
(ry)                          : (ry) or (hr) or (ev)
(for wavefunctions)           :
  30.0                        : (ecut)
(for electron density)        :
 250.0                        : (ecutdens)
(for soft part of density)    :
 70.0                         : (ecutsoft)
*end                          :
                              :
*electronic bands             :
(occupied bands)              :
      8                     : (noband)  No. of occupied bands
(empty bands)                 :
      2                     : (neband)  No. of empty bands
                              :           total No.= noband + neband
(broadening)                  :
       3    500.d0            : (lfermi) = 1:nonmetallic, 2:Fermi, 3:Gaussian,
                              :   lfermi(>3):Methfessel & Paxton, order=lfermi-3
                              : (tfermi) = electronic temp.(K), if metallic
*end                          :
                              :
------------------------------:---------------------------------------
*atoms                        :
(species)                     :
  2                           : (ntype) No. of atomic species
==============================:=======================================
(atomic number)               :
  8.0                         : (zatom)
(pseudopotential)             :
uspp                          : kbpp .or. uspp .or. vand
(nonlocal potential)          :
 .true. 1.5d0 1.25d0 0.8d0    : (lking) .true. = on, (rking, gkgmax, gkgexct)
                              : smoothing parameters
(local potential)             :
 .false. 1.5d0 1.15d0 0.8d0   : (llking) .true. = on, (rlking, glkgmax, glkgexct)
                              : smoothing parameters
(partial core correction)     :
 .true.   1.4d0               : (lpcc) .true. = on, (r_cut) in [a.u.]
 .true. 1.1d0 1.15d0 0.8d0    : (lpking) .true. = on, (rpking, gpkgmax, gpkgexct)
                              : smoothing parameters
                              :
(unit of length)              : only for positions
(ang)                         : (bohr) or (ang)
                              :
(position file)               : 
'control/IN.CONFIG'           :
 2                            : 1:scaled, 2:real coordinates
 1                            : (keyword)
(end)                         :
==============================:=======================================
(atomic number)               :
  1.0                         : (zatom)
(pseudopotential)             :
uspp                          : kbpp .or. uspp .or. vand
(nonlocal potential)          :
 .true. 1.5d0 1.25d0 0.8d0    : (lking) .true. = on, (rking, gkgmax, gkgexct)
                              : smoothing parameters
(local potential)             :
 .false. 1.5d0 1.15d0 0.8d0   : (llking) .true. = on, (rlking, glkgmax, glkgexct)
                              : smoothing parameters
(partial core correction)     :
 .false.   1.4d0              : (lpcc) .true. = on, (r_cut) in [a.u.]
 .true. 1.1d0 1.15d0 0.8d0    : (lpking) .true. = on, (rpking, gpkgmax, gpkgexct)
                              : smoothing parameters
                              :
(unit of length)              : only for positions
(ang)                         : (bohr) or (ang)
                              :
(position file)               : 
'control/IN.CONFIG'           :
 2                            : 1:scaled, 2:real coordinates
 2                            : (keyword)
(end)                         :
*end                          : end of setting *atoms
```

Note that the *TDDFT_MD section can be inserted in between any two sections.  Also note that (check temperature) should be set to .FALSE. for NAQMD runs.  Otherwise, everything else in the IN.PARAM file for a regular (adiabatic) QMD run may remain the same for an NAQMD run, albeit with with *TDDFT-MD section added and set to .TRUE.


## <a name="h7"> 7. Utilities </a>                                                                    
Utilities program are  provided in the util directory in the download. Util directory consists of multiple small fortran and C/C++ program. These small program are used to convert output file generated by QXMD to a formatted structure. This section gives an overveiw of utility program available with QXMD. It should be noted that this list is not comprehensive and any user can create there own utility program to extract data from QXMD output. output file stucture can be referred to section 4.2. 


All the utility program can be compiled using following step.

````
util $ ifort $(util_program).f -o $(util_program)
````
All the utility programs take argument **-h** for help. -h flag prints the argument taken by utility programs. 
````
util $ ./$(util_program) -h 
````

Each program takes path of output data directory as a command line argument. To run the program please eter following steps

````
util $ ./$(util_program) -d $(path of output data directory)
````
Utility files may also take more than one argument depending on their calculation. Please check below about the command line section for each analysis. 

NOTE: All the utility program are written for INTEL Fortran compiler. It is possible the your program will not compile with gfortran. 


### 7.1 Creating PDB file from output 
Program **toPDBcell.f** is used to create the PDB file from QXMD output. This program requires following files in the output directory. 

#### Input file  
````
qm_ion.d 
qm_box.d
md_spc.d 
````
#### Output file 
````
config.pdb 
````
The program will read complete trajectory after a simulation and create a PDB trajectory file. PDB file format can be read by multiple visualization software including VMD, OVITO. 

#### Executing program  
The program takes one argument.
```` 
-d $(Path of the output data directory)
````
Path of the output data directory is mandatory.

To run the program, please use following step
````
Example: 
./toPDBcell -d data
````
It will create a file name config.pdb. It will also print timestep on the terminal 

````
 open : 
 data/qm_ion.d                                                                  
  
 open : 
 data/qm_box.d                                                                  
  
 open : 
 data/md_spc.d                                                                  
  
           0
           1
           2
           3
           4
           5
           6
           7
           ..
           ..
           ..
````
### 7.2 Picking last configuration from simulation 

Program **pick_config.f90** is used to obtain the last configuration from a trajectory. This program is useful to restart the job in case you have not saved restart files. The program requires following files.
#### Input file 
````
qm_ion.d 
qm_box.d
md_spc.d 
qm_cel.d
md_vel.d (optional)
````
#### Output file 
````
IN.VELOC (Velocity data)
IN.CONFIG (Coordinate of each atom)
Lattice_vector.dat (Supercell size)
````

**md_vel.d** file is optional, if **md_vel.d** file exist. Otherwise, program will igonre the option. Program will create three file IN.VELOC, IN. CONFIG and Lattice_vector.dat. IN.VELOC contains velocity data for each atoms. IN.CONFIG contains coordiates of all the atoms in the system and Lattice_vecotr.dat has supercell data.

#### Executing program  
The program takes two argument. 
````
-d $(Path of the output data directory)
-n $(TIMESTEP)
````
Path of the output data directory is mandatory. TIMESTEP selection is optional. If you do not specify -n TIMESTEP option, the program will create output from last snanshot of the trajectory. If you would like to pick certain step from the complete trajectory, specify the timestep. 

To run the program, please use following step.

````
Example:
To pick 5th step from 100 step trajectory
  ./pick_config -d ../examples/01_Water/optimization/data -n 5
To pick last step from 100 step trajecotry 
  ./pick_config -d ../examples/01_Water/optimization/data
````

If the program is executed correctely. It will show output on terminal as following.
````
 open  :
 ../examples/01_Water/optimization/data//md_spc.d                               
                      
 open  :
 ../examples/01_Water/optimization/data//qm_ion.d                               
                      
 open  :
 ../examples/01_Water/optimization/data//qm_box.d                               
                      
 open  :
 ../examples/01_Water/optimization/data//qm_cel.d                               
                      
 open  :
 ../examples/01_Water/optimization/data//md_vel.d                               
                      
Pick up the configuration & velocity files at    100 step
BOX at      0 step [A]
  7.0000,   7.0000,   7.0000  :  angles between cell vec. in [deg.]
 90.0000,  90.0000,  90.0000  :  lengths of cell vectors
CELL at      0 step [A]
  7.00000  0.00000  0.00000   : super cell vector L1
  0.00000  7.00000  0.00000   : super cell vector L2
  0.00000  0.00000  7.00000   : super cell vector L3
````

First 14 lines shows all the file indentified and opened by the program. Next line 
````
Pick up the configuration & velocity files at    100 step
```` 
should correspond to your selected step. 

### 7.3 Create gaussian cube file to visualize wave function
This program is used to create gaussian cube file to visualize wave function. It is to be noted that input files are this program are not outputed as default. You must set **dump wavefunctions** true in obtain the input files **qm_eigv.d**. Wavefunction dumping is extremely I/O intesive and each files are very large. You should dump minimal amount of wavefunction files. 

#### Input file 
````
qm_ion.d 
qm_box.d
md_spc.d 
qm_eigv.d.***
````
#### Output file 
````
state.****.cube 
````
#### Executing program 
The program can take upto 5 argument.
````
-d $(Path_of_the_output_data_directory) 
-n $(frequency_of_snapshot)
-ib $(create cube file from state number)
-eb $(create cube file till state number)
-w  $(selector between wave function or square of wave function )
````

Path of the output data directory is mandatory. **-n** option allows you to choose frequency of cube file created. If you have printed several cube files for each iteration, **-ib** and  **-eb** allows you to choose create dump file in one run. **-w** is optional. Adding -w flag will print wave function &#936;. If you would like to print &#936;<sup>2</sup>, please do not include **-w** flag.  

To run the program, use the 
````
Example:
Create cube files for wave function every 5th step for band 3, 4 and 5 
  ./gcube -d ../examples/01_Water/adiabatic_qmd/data -n 5 -ib 3 -eb 5 -w 
Create cube files for square of wave function every 10th step for band 2, 3, 4 and 5  
  ./gcube -d ../examples/01_Water/adiabatic_qmd/data -n 10 -ib 2 -eb 5  
````

If the program is executed correctely. It will show output on terminal as following.

```
 Open:
 ../examples/01_Water/optimization/data/qm_ion.d                                
                                                                                
                                           
           2
           2
 Open: 
 ./state.02.000110.cube                                                         
                                                                                
                                           
           3
           3
 Open: 
 ./state.03.000110.cube                                                         
                                                                                
                                           
           4
           4
 Open: 
 ./state.04.000110.cube                                                         
                                                                                
                                           
           5
           5
 Open: 
 ./state.05.000110.cube                                                         
                                                                                
                                           
           2
           2
 Open: 
 ./state.02.000120.cube
 ....
 ....
 ....
 ....
```

### 7.4 Kohn Sham eigenvalues  
This program is used to provide Kohn Sham eigenvales in table format. 

#### Input file 
````
qm_eig.d 
qm_fer.d
````
#### Output file 
````
EIG.dat
````
#### Executing program 
The program can take upto 1 argument. Providing the data path is mandotory 
````
-d $(Path_of_the_output_data_directory) 
````

Path of the output data directory is mandatory which is passed as a argument by -d flag. 

To run the program, use the 
````
Example:
Create cube files for wave function every 5th step for band 3, 4 and 5 
  ./eig -d ../examples/01_Water/adiabatic_qmd/data  
````

If the program is executed correctely. It will show output on terminal as following.

```
 open : 
 ../examples/01_Water/optimization/data/qm_eig.d                                
  
 open : 
 ../examples/01_Water/optimization/data/qm_fer.d                                
  
           0
           1
           2
           3
           4
           5
           6
           7
 ....
 ....
 ....
 ....
```

The program will create 'EIG.dat' file. EIG.dat file contains two column. First column represent to TIMESTEP and second column corresponds to eigenvalues for timestep. For example, in case of single water molecule in box example in example directory, we have used 10 total bands. Thus, first row represent to 10 eigenvalues. Timesteps will remain same for first column of first 10 row. Afterthat, it will change to step 1 and next 10 rows will corresponds to 10 eigenvalues for timestep 1. 

````
Format of EIG.dat:
     0 -0.2318E+01
     0 -0.4377E+00
     0 -0.3899E+00
     0 -0.3797E+00
     0  0.3797E+00
     0  0.5791E+00
     0  0.6430E+00
     0  0.6445E+00
     0  0.6462E+00
     0  0.6562E+00
     1 -0.2320E+01
     1 -0.4662E+00
     1 -0.3909E+00
     1 -0.3763E+00
     1  0.3763E+00
 ....
 ....
 ....
 ....
````




## References
<b>References</b><br>
<a name="ref1">1. "First-principles molecular-dynamics simulation of expanded liquid rubidium,"
   F. Shimojo, Y. Zempo, K. Hoshino & M. Watabe,
   <a href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.52.9320"><i> Phys. Rev. B</i>, <b>52</b>, 9320 (1995)</a><br></a>
<a name="ref2">2. "A divide-conquer-recombine algorithmic paradigm for large spatiotemporal quantum molecular dynamics simulations,"
   F. Shimojo, S. Hattori, R. K. Kalia, M. Kunaseth, W. Mou, A. Nakano, K. Nomura, S. Ohmura, P. Rajak, K. Shimamura & P. Vashishta
   <a href="http://aip.scitation.org/doi/abs/10.1063/1.4869342?journalCode=jcp"><i>J. Chem. Phys.</i> <b>140</b>, 18A529 (2014)</a><br></a>
<a name="ref3">3. "Linear-scaling density-functional-theory calculations of electronic structure based on real-space grids: design, analysis, and scalability test of parallel algorithms,"
   F. Shimojo, R. K. Kalia, A. Nakano & P. Vashishta
   <a href="https://www.sciencedirect.com/science/article/pii/S0010465501002478"><i>Comput. Phys. Commun.</i> <b>140</b>, 303 (2001)</a><br></a>
<a name="ref4">4. "Embedded divide-and-conquer algorithm on hierarchical real-space grids: parallel molecular dynamics simulation based on linear-scaling density functional theory,"
   F. Shimojo, R. K. Kalia, A. Nakano & P. Vashishta,
   <a href="http://www.sciencedirect.com/science/article/pii/S0010465505000688"><i>Comput. Phys. Commun.</i> <b>167</b>, 151 (2005)</a><br></a>
<a name="ref5">5. "Divide-and-conquer density functional theory on hierarchical real-space grids: parallel implementation and applications,"
   F. Shimojo, R. K. Kalia, A. Nakano & P. Vashishta,
   <a href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.77.085103"><i>Phys. Rev.</i> <b>B</b> 77, 085103 (2008)</a><br></a>
<a name="ref6">6. "Scalable atomistic simulation algorithms for materials research,"
   A. Nakano, R. K. Kalia, P. Vashishta, T. J. Campbell, S. Ogata, F. Shimojo & S. Saini
   <a href="http://ieeexplore.ieee.org/document/1592786/"><i>Proc. Supercomputing, SC01</i> (ACM/IEEE, 2001)</a><br></a>
<a name="ref7">7. "Metascalable quantum molecular dynamics simulations of hydrogen-on-demand,"
   K. Nomura, R. K. Kalia, A. Nakano, P. Vashishta, K. Shimamura, F. Shimojo, M. Kunaseth, P. C. Messina& N. A. Romero
   <a href="http://ieeexplore.ieee.org/document/7013041/"><i>Proc. Supercomputing, SC14</i> (IEEE/ACM, 2014)</a><br></a>
<a name="ref8">8. "Large nonadiabatic quantum molecular dynamics simulations on parallel computer,"
   F. Shimojo, S. Ohmura, W. Mou, R. K. Kalia, A. Nakano & P. Vashishta
   <a href="https://www.sciencedirect.com/science/article/pii/S0010465512002548"><i>Comput. Phys. Commun.</i> <b>184</b>, 1 (2013)</a><br></a>
<a name="ref9">9."Crystalline anisotropy of shock-induced phenomena: omni-directional multiscale shock technique,"
   K. Shimamura, M. Misawa, S. Ohmura, F. Shimojo, R. K. Kalia, A. Nakano & P. Vashishta
   <a href="http://aip.scitation.org/doi/abs/10.1063/1.4942191?journalCode=apl"><i>Appl. Phys. Lett.</i> <b>108</b>, 071901 (2016)</a><br></a>
<a name="ref10">10."Molecular dynamics with electronic transitions." Tully, John C. <a href="http://aip.scitation.org/doi/abs/10.1063/1.459170"><i>J. Chem. Phys. </i> <b>93.2</b>, 1061-1071 (1990)</a><br></a>
<a name="ref11">11. "Time-dependent density-functional theory", Gross, E. K. U., and W. Kohn. <a href="https://www.sciencedirect.com/science/article/pii/S0065327608606000"> <i>Adv. Quantum Chem. </i> <b>21</b>, 255-291, (1990)</a><br></a>
<a name="ref12">12. "Decoherence-induced surface hopping." Jaeger, H. M., Fischer, S., & Prezhdo, O. V. <a href="https://aip.scitation.org/doi/abs/10.1063/1.4757100"><i>J. Chem. Phys. </i> <b>137</b>, 22A545 (2012)</a><br></a>
<a name="ref13">13. "Time-Dependent Density Functional Response Theory for Molecules" M. E. Casida, <a href="http://homepage.univie.ac.at/mario.barbatti/papers/method/tddft_casida_1995.pdf"> Recent Advances in Density Functional Methods (Part I), edited by D. P. Chong (World Scientific, Singapore, pp. 155-192 (1995))</a><br></a>
<a name="ref14">14. "A long-range-corrected time-dependent density functional theory." Tawada, Yoshihiro, et al.  <a href="https://aip.scitation.org/doi/abs/10.1063/1.1688752"><i> J. Chem. Phys. </i> <b>120.18</b>, 8425-8433 (2004) </a><br></a>
<a name="ref15">15. "Projector augmented-wave method" P. E. Blochl <a href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.50.17953"><i>Phys Rev B</i> <b>50</b>, 17953 (1994)</a><br></a>
<a name="ref16">16. "Generalized Gradient Approximation Made Simple" J. P. Perdew, K. Burke & M. Ernzerhof <a href="https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.77.3865"> <i>Phys. Rev. Lett.</i> <b>77</b>, 3865 (1996)</a><br></a>
<a name="ref17">17. "Hybrid functionals based on a screened Coulomb potential." <a href="https://aip.scitation.org/doi/abs/10.1063/1.1564060"><i>J. Chem. Phys.</i> <b>118</b>, 8207-8215 (2003)</a></br></a>
<a name="ref18">18. "Improved tetrahedron method for Brillouin-zone integrations." Blöchl, Peter E., Ove Jepsen, and Ole Krogh Andersen. <a href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.49.16223"><b>Phys. Rev. B</b> <i>49</i>, 16223 (1994)</a></br></a>
<a name="ref19">19. "Ab initio molecular dynamics in a finite homogeneous electric field."  Umari, P., and Alfredo Pasquarello. <a href="https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.89.157602"><b> Phys. Rev. Lett.</b> <i>89</i>,157602 (2002)</a></br></a>
1.  R. D. King-Smith, M. C. Payne, and J. S. Lin "Real-space implementation of nonlocal pseudopotentials for first-principles total-energy calculations" Physical Review B 44, (1991) :13063
