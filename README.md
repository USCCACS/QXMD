# QXMD : ab initio Molecular Dynamics Simulation Package
**QXMD** is a scalable, parallel program for Quantum Molecular Dynamics simulations with various eXtensions. Its simulation engine is based on (time-dependent) density functional theory using pseudopotentials and plane-wave basis, while extensions include nonadiabatic electron-nuclei dynamics and multiscale shock technique. **QXMD** serves as a community-development platform for new methods and algorithms, a research platform on high-end parallel supercomputers, and an educational platform for hands-on training.
 
## 0. Prerequisites
The parallel version of **QXMD** requires a Fortran compiler, as well as FFT and MPI libraries. There are no other library requirements for the QXMD.

## 1. Getting Started
To get started, clone this repository to your computer.
```
~$ git clone https://github.com/USCCACS/QXMD.git
```

## 2. How to build QXMD

### 2.1 Working Directory
First, change working directory to QXMD/
```
~$ cd QXMD
```

You will see following files and directories when you use the 'ls' command:
```
QXMD $ ls
Include/   Examples/   Lib/    Makefile    QXMD_Manual_v3.2/    Sources/     util/
```

Include/ contains necessary libraries for compilation, Examples contains sample input files for various example simulations, Lib contains pseudopotential files, Makefile will be used to build the program, QXMD_Manual_v3.2 is the manual for this software packgae, Sources/ contains all QXMD source codes,and util/ contains helpful codes for post-processing output data from QXMD.

### 2.2 Configure Makefile
You have to choose makefile based on the configuration of your machine. For most of the facility supercomputers makefile is already configured. Complete list of the preconfigured makefile can be obtain by typing make help.
```
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
```

If you are not using aforementioned machine, you might need to modify the Makefile according to your computing environment.

### 2.4 Precompiler variables
Precompiler flags are defined in CPPDEFS as:
```
CPPDEFS= -DFLAG1 -DFLAG2
```

Following flags are defined in QXMD
```
-LIBFFTW = link FFTW library 
-LIBFFTW3 = link FFTW3 library 
-DPOINTER64 = used for 64-bit machines 
-DSSL2VP = used on Fujitsu machine 
-DVECTOR = special implementation. It may result in faster computation on some machine. 
```

There are several flags are defined for future development. Currently, they are not used. Several flags are specific to compiler. Please check compiler manual for further details.

### 2.5 Adding New Machine configuration in Makefile
Makefile defines which compiler you will use to build the QXMD executable.There are several predefined compiler setting which must be set. First, open Makefile available in QXMD directory and define the machine name. For example, theta machine at Argonne National Lab
```
theta: $(WDIR) $(DDIR) $(SOURCEFILE)
        sed "s/^#THETA#//" $(SOURCEFILE) > $(WORKFILE)
```

This process will ensure that make will understand the target. Next, we need to define the macro such as LINKER, LDFLAGS and FFTLIB in the makefile included in Source directory. For example, theta machine definition are as follow
```
#THETA#     LINKER        = ftn
#THETA#     FFLAGS        = -c -fast -warn none                                                                     
#THETA#     LDFLAGS       = -O3 -ipo -no-prec-div 
#THETA#     FFLAGS        = -c -O3 -ipo -no-prec-div  -warn none                                              
#THETA#     CPPDEFS       = -traditional $(CPPDEF) -DLIBFFTW3 -DIFC -DVECTOR                                         
#THETA#     MKLPATH       = $(MKLROOT)/intel/lib64
#THETA#     FFTLIB        = -L$(MKLPATH) -Wl,--start-group -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -Wl,--end-group
```

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
After tailoring the Makefile for your computing environment, type the command below to build the QXMD executable.
```
QXMD $ make $(MACHINE_NAME)
QXMD $ make qxmd
```

To create a parallel version of QXMD, please use the following steps
```
QXMD $ make $(MACHINE_NAME)
QXMD $ make qxmd_mpi
```

Check to see if you the QXMD executable has been created, and if so, you are ready to start a simulation.

```
qxmd $ ls
Include/    Makefile    Qxmd/    qxmd_mpi Sources/     util/
```

Parallel build is activated by default. The program will use all available core on the machine to build the program efficiently.

## 3. How to run
In order to run QXMD, there are a few mandatory directories and files that must be present and in the correct hierarchy (Fig. 1). The working directory, from which a **QXMD** job is launched must contain a directory called “data/”, where all output data will be dumped, as well as a directory called “control/”, which must contain the following:
```
•	CONFIG: a configuration file detailing the atomic coordinates of the system to be simulated
•	filename: a simple text file containing the path to the main input file
•	IN.PARAM: an input parameter file with various settings for the dynamics simulation 
•	VELOC: an optional initial velocity files containing the initial three component velocities for each atom in the system
•	NCPP/ or PAW/ or USPP/: a directory containing pseudopotential files for each atomic species in the system
```

There are many example input files for various types of simulations in the Examples/ directory, including optimization of water, adiabatic QMD of water in the canonical ensemble, non-adiabatic QMD of monolayer MoSe2 in the microcanonical ensemble, and a MSST simulation of Si.  The NAQMD and MSST examples are explained in more detail in the sections below. 

To learn more about QXMD, please refer to [QXMD manual](https://github.com/USCCACS/QXMD/blob/master/QXMD_Manual_v3.2.md).

## 4. License
This project is licensed under the GPU 3.0 license - see the [LICENSE](https://github.com/USCCACS/QXMD/blob/master/LICENSE) file for details

## 5. Publications

1) F. Shimojo, Y. Zempo, K. Hoshino, and M. Watabe, "First-principles molecular-dynamics simulation of expanded liquid rubidium," Physical Review B, vol. 52, pp. 9320-9329, Oct 1 1995. <br>

2) F. Shimojo, S. Hattori, R. K. Kalia, M. Kunaseth, W. Mou, A. Nakano, et al., "A divide-conquer-recombine algorithmic paradigm for large spatiotemporal quantum molecular dynamics simulations," Journal of Chemical Physics, vol. 140, 18A529, May 14, 2014.


