.. _introduction:

*******************************************
Introduction of pydlpoly's general concepts
*******************************************

.. _motivation:

Motivation
==========

The ``pydlpoly`` project was initiated out of a need. Before we had used the ``Tinker`` code (V 4.2)
which was not parallel and limited to about 3000-5000 atoms. We looked at ``lammps`` and ``DL_Poly 2``.
Because of the simpler structure and the fact that ``DL_Poly 2`` was not using a domain decompostion to be 
data-parallel we chose it as a basis for wrapping the code with ``f2py``. Since all the data like atom positions,
forces, velocities etc. are available on all the nodes the code is limited to about 30,000 atoms, which is still
an order of magnitude more than what was possible for us with ``Tinker``. Our target is not so much to investigate
larger systems but to improve the potential energy functions for the hybrid systems we are interested in (like MOFs).
We are working on polarizable and "semi-reacive" force fields, which means that the computational cost per atom will 
get larger in the future. We initalliy started with V2.20 of ``DL_Poly``, which was later converted and renamed to 
``DL_Poly Classsic``, distributed under BSD-license. ``pydlpoly`` is a fork of ``DL_Poly CLassic`` and is currently 
available to academic groups on the basis of a "collaborational license".

.. _features:

Features and Missing Features
=============================

The idea behind ``pydlpoly`` is to use the F90 (with a lot of F77 flair :-) MM energy engine as a backend and to do all
the operations like geometry optimization, MD startup, analysis on Python level. Few things have been added into the F90
code, whereas most extensions and features are implemented on the scripting level. Some things suffer from efficiency issues on
Python level and will have to be converted to F90 at some point. Please note that ``pydlpoly`` is a growing and evolving 
project (as is this documentation ... it will always be outdated). A very important point is that ``pydlpoly`` has been developed
with our projects in mind. This means it supports all we needed to do for our work and probably a number of things are missing.
Most importantly ``pydlpoly`` -- even though it is extended in severel ways by the Python frontend -- must be seen as a "cripled"
``DL_Poly Classic`` since all the features implemented in ``Dl_Poly Classic`` which were not needed might be broken. This means some
of these features (like hyperdynamics, solvation, metal potentials etc.) might still work in the remaining F90 routines
(some might not) but none of these things are supported by the Python frontend. If you want to use these features you are probably
better off if you use the ``DL_Poyl Classic`` code itself.

Changes/Additions to ``DL_Poly`` (F90-level)
--------------------------------------------

* vdW lookup tables use cubic spline interpolation and forces are calculated from analytic derivative of cubic
  spline to be consistent with energy (original ``DL_Poly 2`` uses interpolation for energy and force giving inconsistent
  forces leading to non energy conserving microcanonic MD).
* vdW and short ranged Coulomb (either truncated with shift damping or the short range part of Ewald/SPME) are switched to
  give continuous forces at the cutoff.
* Gaussian charges have been added for all types of Coulomb interactions (truncated, Ewald and SPME) with a different
  width for each atom.
* A number of potentials have been added like MM3 type bond and stretch terms, four fold torsion or Feynman-Hibbs potentials force
  vdW interactions.
* Switching of non-bonded interactions (vdW and Coulomb) on a per molecule basis is added (also kinetic energy can be switched off 
  on a per molecule basis).
  
Features
--------

* Access to all internal data structures like positions, forces, velocities.
* Efficient L-BFGS optimizer for molecular structures and a robust minimzer for periodic systems was added.
* Energy calculations are performed from a central Python level function, which allows to add arbitrary energy 
  contributions from other Python modules (like QMMM etc)
* The central MD loop is controlled on Python level, which allows complex scripting.
* The input for the underlying ``DL_Poly`` engine (FIELD, CONFIG, CONTROL) are generated completely automatic
  from input files inspired by the ``Tinker`` format. In particular, the force field is defined using atom types in a general
  force filed file definition (key-file) which is parsed to produce the FIELD input
* Setup allows to automatically add e.g. guest molecules to a periodic system or to fill a box with solvent molecules.
* Simple mechanic embedding QMMM (currently with ``TURBOMOLE`` as QM engine).
* Restart and trajectory information is written to a hdf5 restart file allowing complex scripting by grouping the
  data in "stages"
  
Missing Features (``DL_Poly`` features not supported any more)
--------------------------------------------------------------

* Only Velocity Verlet but no Leap frog propagator from originial code
* Many special potentials (metal potentials etc.) not supported by frontend
* Hpyerdynamics and solvation features not supported
* Metadynamics implemented in ``DL_Poly Classic`` is not available. This can be replaced by an experimental 
  coupling to the ``plumed`` metadynamics library (we have written a ``f2py`` wrapped ``pyplumed`` for this)
* probably more I am not even aware of ...

.. _dependencies-intro:

Dependencies
============

For details see the Installation section. This is just an overview over the major libraries one needs for ``pydlpoly``

* Python 2 (not 3!!) .. I guess you need at least 2.5
* Numpy (which includes ``f2py``)
* mpi4py (Cython based MPI wrapper ... we use OpenMPI)
* h5py (wrapper for hdf5 library)
* fftw2 (``DL_Poly`` is interfaced to fftw2 only and not to fftw3) 

Status and Outlook
==================

This project is constantly developing. Currently we are working on QEq based polarization models.
Most importanly the force field paramter fitting, which was built on the ``Tinker`` engine is currently 
replaced by completely new pythonic ff_generator using ``pydlpoly`` as MM engine.

Contributors
=============

Concept, Idea and most of the work by R. Schmid
Further contributions by

* Chris Spickermann (Gaussian charges, inital QEq)
* Mohammad Alagehmandi (initial versions of GCMD)
* Johannes P. Dürholt (ff_gen infrastructure, addtions to QM/MM)
* Niklas Siemer (FH-potentials)
