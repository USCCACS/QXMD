Program Structure
=================

Below is a tree diagram of the general structure for the directories and
files needed to run a **QXMD** simulation.

::

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

An overview is given below, while a more detailed explanation can be
found in section 4.1.

3.1 Control Directory
~~~~~~~~~~~~~~~~~~~~~

The **control** directory is where QXMD will search for all input data.
The main input file, in this case named **IN.PARAM**, contains the
settings for all input parameters for running QXMD. It is best practice
to name this file 'IN.PARAM', though it can, in principle, be named
anything. Thus, it is required to have a file named 'filename' which is
a simple, one-line text file holding the path (relative to the QXMD
executable) to the main input file, in this case 'control/IN.PARAM'.
Most QXMD runs will also require a file detailing the initial
configuration of the simulation system (this may be omitted if you are
restarting a QXMD run from a previous run). It is best practice to name
the input configuration file 'IN.CONFIG', as shown above in the tree
diagram. Optionally, a file defining the inital velocities of all the
atoms may be defined. It is best practice to name the input velocity
file 'IN.VELOC', as shown above in the tree diagram. If not intial
velocities are provided, random initial velocities will be used
corresponding to the intial temperature of the system. Finally, a
directory holding pseudopotential information for all atomic species
present in your simulation system is required. The **PAW** directory,
shown above, holds pseudopotential information for hydrogen in **H.PBE**
and oxygen in **O.PBE**, which would both be required for simulating,
say, a water molecule.


3.2 Data directory
~~~~~~~~~~~~~~~~~~

A directory named **data** is required for any QXMD run. All output data
will be dumped to this directory. If there are files already present in
this directory at the start of a new QXMD run, they will be overwritten.
Thus, it is recommended to move files in **data** from a previous run to
a new directory before beginning a new QXMD run. If you are restarting a
QXMD run, there are some required files that must be present in
**data**, which will be detailed in section 4.2.


3.3 QXMD executable
~~~~~~~~~~~~~~~~~~~

The QXMD executable, in this case **qxmd_mpi** should be at the same
directory level as the **control** and **data** directories.
