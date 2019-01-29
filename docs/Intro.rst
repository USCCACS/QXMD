Introduction
============

QXMD is a Quantum Molecular Dynamics (QMD) simulation software with
various eXtensions. QMD follows the trajectories of all atoms while
computing interatomic forces quantum mechanically in the framework of
density functional theory (DFT). :cite:`RN114` The
QXMD software has been developed by Fuyuki Shimojo since 1994. :cite:`RN334`
Since 1999, various extensions have been developed in collaboration with
Rajiv Kalia, Aiichiro Nakano and Priya Vashishta. :cite:`RN182` The basic QXMD
code is based on a plane-wave basis to represent electronic wave
functions and pseudopotential (PP) methods to describe electron-ion
interaction. Supported PPs include norm-conserving PP :cite:`RN116` and ultrasoft PP :cite:`RN211` , as well as an all-electron
projector augmented-wave (PAW) method :cite:`RN263`. Electron-electron interaction beyond the mean-field
Hartree approximation is included using various exchange-correlation
functionals, with and without spin polarization: generalized gradient
approximation (GGA) :cite:`RN277`, DFT+U method for transition metals :cite:`RN322`, van der Waals (vDW) functional for molecular crystals and
layered materials :cite:`RN21`,
nonlocal correlation functional :cite:`RN126`, and range-separated exact-exchange functional :cite:`RN324` :cite:`RN99`. Various
unique capabilities included in the QXMD code (some of which are
described in :cite:`RN181`) include:

.. raw:: html

   <ul>

1. Linear-scaling DFT algorithms :cite:`RN181` :cite:`RN184` :cite:`RN128` 
   
2. Scalable algorithms on massively parallel computers :cite:`RN139` :cite:`RN139` :cite:`RN64`
   
3. Nonadiabatic quantum molecular dynamics (NAQMD) to describe excitation dynamics :cite:`RN329` 

4. Omni-directional multiscale shock technique (OD-MSST) to study shock response of materials :cite:`RN250` :cite:`RN251`


.. bibliography:: references.bib
   :style: unsrt
