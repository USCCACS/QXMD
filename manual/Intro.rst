Introduction
============

QXMD is a Quantum Molecular Dynamics (QMD) simulation software with
various eXtensions. QMD follows the trajectories of all atoms while
computing interatomic forces quantum mechanically in the framework of
density functional theory (DFT). :cite:`RN1` The
QXMD software has been developed by Fuyuki Shimojo since 1994. :cite:`RN2`
Since 1999, various extensions have been developed in collaboration with
Rajiv Kalia, Aiichiro Nakano and Priya Vashishta. :cite:`RN3` The basic QXMD
code is based on a plane-wave basis to represent electronic wave
functions and pseudopotential (PP) methods to describe electron-ion
interaction. Supported PPs include norm-conserving PP :cite:`RN4` and ultrasoft PP :cite:`RN5`
. Electron-electron interaction beyond the mean-field Hartree approximation is included using various exchange-correlation
functionals, with and without spin polarization: generalized gradient approximation (GGA) :cite:`RN6`, DFT+U method for transition metals :cite:`RN7`, van der Waals (vDW) functional for molecular crystals and
layered materials :cite:`RN8`,
nonlocal correlation functional :cite:`RN9`, and range-separated exact-exchange functional :cite:`RN10` :cite:`RN99`. Various
unique capabilities included in the QXMD code (some of which are
described in :cite:`RN12`) include:

.. raw:: html

   <ul>

1. Linear-scaling DFT algorithms :cite:`RN12` :cite:`RN13` :cite:`RN14` 
   
2. Scalable algorithms on massively parallel computers :cite:`RN15` :cite:`RN16`
   
3. Nonadiabatic quantum molecular dynamics (NAQMD) to describe excitation dynamics :cite:`RN3` 

4. Omni-directional multiscale shock technique (OD-MSST) to study shock response of materials :cite:`RN17` :cite:`RN18`


.. bibliography:: intro.bib
   :style: unsrt
