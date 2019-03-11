.. molsys documentation master file, created by
   sphinx-quickstart on Mon Aug 21 14:29:21 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


General Usage
#############

FFgen makes it easy to parameterize FFs in respect to arbitrary objective functions.
It uses pydlpoly as molecular mechanics machine and relies in data stored in an hdf5
file. Objectives functions are implemented in objective modules. These has to be programmed
in a predesigned fashion. On the basis of these objective function, FFgen mades a general
objective function available, which can be minmized/maximized by arbitrary optimizers.

Here the general usage is illustrated employing the :ref:`ric_fit` objective function together 
with the :ref:`CMA-ES` optimizer using biphenyl as example. This example can be found in
the FFgen distribution in the directory examples/biphenyl.

.. literalinclude:: ../examples/biphenyl/fit.py
    :language: python
    :linenos:
    :lines: 13-50

In the lines 1-26, the necessary modules are imported. Molsys, Pydlpoly and the ric_fit class
are instanciated. Then the actual FFgen is started. The following options are possible for this:

.. autoclass:: ff_gen.ff_gen.ff_gen

In the shown example the **name** is set to "ph-ph". From this follows that a 
reference file called "ph-ph.hdf5" in the hdf5 file format is available. In this reference file
the objective specific reference information is stored together with a required **system** Group.
In this group two datasets has to be available. An integer called **natoms** specifying the 
number of atoms of the reference system and a list called **elements** with elementy symbols 
of all atoms in the reference system.

In addition the already created ric objective object is handed over together with the pd
object. Since we want to use CMA-ES as optimizer we set **minimize** to True. 

Afterwards the optimizer object is created, and the actual optimization is performed. In
the end it is necessary to call the method **finish**. Since this method writes out the
par and ric file with the optimized parameters. In addition objective specific information
is written. All this information could be found in the run directory after the optimization
is finished.

In the example case the optimizer is started at the point in the search landscape defined
by the initial values of the parameters which should be optimized. These initial variables
are an attribute of the **ff_gen** class and are called **initials**. In addition the 
upper and lower bounds for the parameters are available under **pmin** and **pmax**. The 
optimizer could also be started from a random guess. For this purpose one can use
the method **random_guess**. The method **get_bounds** is resposible to introduce 
weak/hard margin information into the bounds.
