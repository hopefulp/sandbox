.. role:: bash(code)
   :language: bash
   
.. role:: makefile(code)
    :language: makefile 

.. _installation:

**********************************
Compiling and Installing pydlpoly
**********************************

Getting the code
================

This assumes you are a member of the CMC-group or have at least read access to our main repository,
which lives on Rochus's workstation legolas.aci.rub.de. Get yourself ``mercurial`` installed (or better
use ``TortoiseHg``, which is a very nice GUI to mercurial) and clone the repo like this:

.. code-block:: bash

    cd
    mkdir sandbox
    cd sandbox
    hg clone ssh://legolas.aci.rub.de//home/repo/repos/pydlpoly
    
This should produce the pydlpoly source files into the ``~/sandbox/pydlpoly`` directory.

If you are not a CMC-group member, you probably have tgz file to unpack. For the further discussion we assume you do this 
in the same directory.

Generating the Docu
===================

For building the docu (these files) you need to have the ``sphinx`` documentation system installed. On Ubuntu variants it is usally
the package ``python-sphinx`` which needs to be installed.
Just go into ``doc/`` and ``make html``. To open open ``doc/_build/html/index.html`` in your favorite browser.

Compiling the code
==================

Dependencies
------------

As overviewed in the Introduction (:ref:`dependencies-intro`) the compilation of ``pydlpoly`` depends on a number of Python
packages. First you need a F90 compiler. We have used both ``gfortran`` and ``ifort`` on our machines, but others will work as well.
Then you need the development packages (libs&includes) of Python and Numpy. Numpy contains ``f2py``, which is needed to generate
the F90/Python bindings.
On Ubuntu and derived distributions it is often sufficient to use the package manager to install the Python wrapper packages in order 
to install also the underlying libraries itself.

* ``mpi4py`` : The MPI bindings for Python - in addition the MPI library incuding the development options will be needed. We
  prefer to use OpenMPI.
* ``h5py`` : The hdf5 library bindings for Python - in addition the hdf5 lib incuding development options will be installed
* ``fftw2`` : You need the Fastest Fourier Transform in the West Version 2 (not 3!). On some distributions only fft3 is available as
  a package. In this case you need to compile it manually. You have at least to add ``-fPIC`` because `f2py`` will generate a 
  shared object loaded at runtime. 
* ``blas`` : You need a ``libblas`` to get linked in. This is not an extremly critical lib in terms of performance so any "vanilla"
  BLAS-library can be used, as long as it is compiled ``-fPIC``.
  
.. _makefile:

Makefile
--------

The F90 sources, the ``f2py`` interface definition (_pydlpoly.pyf) and the Makefile live in the ``source/`` directory. 
Depending on the position of the libraries you might have to adapt the Makefile. The following is our standard Makefile with two
sections, one for the typical LINUX worstation (using gfortran) and a special variant on our cluster (hostname cmcc), where
all packages (including Python/Numpy, OpenMPI, hf5 etc.) are manually compiled with ``ifort`` and using MKL. 

.. code-block:: makefile

    FC=mpif90
    LD=mpif90 -o

    ifneq ($(HOSTNAME),cmcc.cmcc.cluster) 
    #================== GNU Fortran, MPI version ==============================
    # gfortran (SUSE version , fftw2 is called libdfftw
    LDFLAGS=-O2 -ffast-math
    FFLAGS=-c -O2 -ffast-math -fPIC
    LDLIBS=-L/usr/local/lib -lfftw -lblas
    F2PY_COMP=gnu95
    else
    ## cmcc version using ifort and fftw and MKL Blas
    LDFLAGS=-O3 $(LDLIBS)
    LDLIBS=-L/usr/local/lib -lfftw
    FFLAGS= -FI -c -O3 -xhost -fPIC
    F2PY_COMP=intelem
    endif
    
    #=====================================================================

It is important to set the ``F2PY_COMP`` variable to the proper value for your F90 compiler. Check the
manual and help of ``f2py`` to get this right.

Compiling and Installing
------------------------

If your Makefile is properly configured a simple ``make`` in the ``source/`` directory will compile all the fortran code,
generate the binary _pydlpoly.so, which will be copied to the ``so/`` directory. If this fails you have most probably a dependency
problem. In this case please check the error message to fix it.
In order to run the code you need to add the two directories ``py/`` and ``so/`` to your ``PYTHONPATH`` environment variable.
Usually this is doen in .bashrc (provided you use ``bash``) as your shell.

.. code-block:: bash

    export PDLPDIR=/home/rochus/sandbox/pydlpoly
    export PYTHONPATH=$PDLPDIR/so:$PDLPDIR/py:$PYTHONPATH

To test if everything works allright just import the ``pydlpoly`` module in your interactive Python session like this::

    Python 2.7.5+ (default, Feb 27 2014, 19:37:08) 
    [GCC 4.8.1] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import pydlpoly
    Intialized Random with System Time
    >>>
    
If there are no error messages you are ready to go!




