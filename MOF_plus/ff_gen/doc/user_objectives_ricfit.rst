.. molsys documentation master file, created by
   sphinx-quickstart on Mon Aug 21 14:29:21 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _ric_fit:

RIC_FIT
#######

Introduction
============

The ric_fit objective implements the standard objective function used for MOF-FF. 
It optimizes the geometry of the system and calculates afterwards the hessian of the
optimized geometry. Then both geometry and hessian are projected into redundant internal
coordinates using the ric module from Raimondas Galvelis. The set of rics is defined
by the connectivity of the system. Then the fitness is calculated as msd between geometry
and hessian of the model and the reference geometry and hessian.

How to use?
===========

It is assumed, that there is a mol object with an attached ric addon class instance.
The use of the ric addon is described in more detail in the molsys documentation. The 
following code shows how an ric addon can be initialized and attached. This has to be
done after the pydlpoly setup. It is recommended to use the **full=False** flag in the **setup_rics** routine. This forces the program
to use only those internal coordinates as rics which are really used in a FF potential.

.. code-block:: python 

    >> m.addon("ric")
    >> m.ric.setup_rics(full = False)

In the next step the ric_fit objective function module needs to be importet and initialized.

.. code-block:: python

   >> import ff_gen.objectives.ric_fit3 as ric_fit
   >> ric = ric_fit.ric_fit()

There are a lot of setup option to tweak the behavion of the ric_fit instance.

.. autoclass:: ff_gen.objectives.ric_fit3.ric_fit

The ric_fit class reads its reference information by the help of the reference class out
of an hdf5 file. This hdf5 file has to comprise the following data structure in order to be
usable by the ric_fit class: 

    * In the root group has to be a group called **hessians**
    * **hessians** can contain several hessian subgroups storing coordinates
      and hessian on different points of the PES. The name is of a hessian subgroup
      is the tag which is a keyword argument for the init of a ric_fit instance. By defaults
      the name **primary** is assumed
    * A hessian subgroup has to include the Datasets:
        + **coord**: (N,3) array of atomic coordinates
        + **hessian**: (3N, 3N) array storing the actual hessian

After the initialization of the ric_fit instance, the actual FFgen is initialzed and the
the **ric** object is commited.
