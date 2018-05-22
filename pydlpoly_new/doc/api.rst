.. _api:

********************************************
Documentation of pydlpoly's API
********************************************

Starting up pydlpoly
####################

In order to start up pydlpoly you have to import the module `pydlpoly` and
instantiate an object ot the pyldpoly class like this::

    import pydlpoly
    pd = pydlpoly.pydlpoly("test")
    
This generates a still "empty" pydlpoly object `pd` with the default runname "test". The default
runname is used as a default for input file prefixes etc..


.. currentmodule:: pydlpoly
    
.. autoclass:: pydlpoly
    :members:  __init__, setup

.. note::

    By default a standard CONTROL file is generated when calling :class:`pydlpoly.setup` which is generated using a dictionary
    of the `pydlpoly` class. The default entries are the following::
        
        "timestep" : 0.001 
        "cut"      : 12.0  
        "delr"     : 1.0  
        "shift"    : "damping"
    
    In order to change these settings you can directly change these values or add other options.
    Not all keywords from DL_Poly classic are still used in pydlpoly (like for example temp or pressure are set in `MD_init`).
    In order to use a 2 fs timestep and a 10.0 A cutoff you have to change the values *before* calling `setup` like this::
    
        pd.control["cut"] = 10.0
        pd.control["timestep] = 0.002
        
    A second way is to specifiy an explcit control file in setup using the `control` parameter.

.. warning::

    Internally the `pydlpoly` class imports the `_pydlpoly.so` library, which contains the python wrapped
    F90 code. All the variables like xyz coordinates etc. are held in allocatable F90 module level arrays, which exist only
    once (like common blocks). Therefore it is not possible to generate two instances of the :class:`pydlpoly.pydlpoly` class, since they
    would share the same arrays. To prevent this an error is raised if you try to generate a second `pydlpoly` incarnation.
    

Optimization Methods
####################

The following methods should be used for the structure optimization (periodic or nonperiodic) or the optimization of lattice and 
structure.

.. autoclass:: pydlpoly
    :members: MIN_lbfgs, LATMIN_sd
        
.. note::

    The force used in the LATMIN optimization is symmetrized according to the boundary conditions. Thus, if for example `bcond = 1`,
    which means cubic boundary conditions the stress tensor is symmetrized to be diagonal with averaged values on the diagonal.
    Make sure to use proper boundary conditions. In order to allow a cubic system to become tricilinc you need to explcitly set
    `bcond = 3` in :class:`pydlpoly.setup`.
    
Molecular Dynamics Methods
##########################

The MD simulations (starting form the current coordinates) are setup by two methods. The MD parameters are setup up by `MD_init`
and the simulation itself is run with `MD_run`. You can use mutliple `MD_run` calls to just continue after doing something like
storing the current coordinates or whatever you like. Note that the same can be accomplished by passing a function `do_every_step`
to `MD_run`.

The following ensemble types are available:

* `nve` : Microcanonic Ensemble (no Thermostat/Barostat)
* `nvt` : Canonic ensemble with Thermostat, either `hoover` (Noose-Hoover) or `ber` (Berendsen)
* `npt` : NPT with isotropic pressure (for liquids with cubic bcond)
* `nst` : NsigmaT with non-isotropic pressure (make sure to have bcond = 3)

.. autoclass:: pydlpoly
    :members: MD_init, MD_run
    
.. note::

    Also during MD in the NPT or NST ensemble one has to make sure to use the proper boundary condition with `bcond`.
    See the corresponding note for the Optimization methods.
    
Advanced Methods
################

the following methods can be used in various scripts to do more advanced things as regular optimization and MD.

Energy calculation
------------------

.. autoclass:: pydlpoly
    :members: calc_energy, calc_energy_force
    
Getting and setting various things
----------------------------------

The following methods allow to get and set various properties like atom positions and velocities (of the complete
system or molecular subset) or the cell paramters.

.. autoclass:: pydlpoly
    :members: get_natoms, get_tstep, get_elements, get_atomtypes, get_xyz, get_subset_xyz, set_xyz, set_subset_xyz,\
              get_vel, get_subset_vel, set_vel, set_subset_vel, set_atoms_moved, get_cell, set_cell, get_frac, set_frac, \
              get_cell_volume, get_stress, get_temperature, get_bcond, get_charges, set_charges, get_masses, get_subset_masses
              
Temperatures and Thermostating
------------------------------

.. autoclass:: pydlpoly
    :members: set_thermostat_sigma, get_degfree, set_degfree

Additional systems or potential terms
-------------------------------------

An extra term is an additional energy expression for the molecular system which is evaluated after the regular molecular mechanics
energy evaluation. An extra system also contains additional degrees of freedom and the extra system must be a class with a defined 
interface containing methods to propagate the system with a velocity verlet propagator etc. (see special docu on extra systems).
                            
.. autoclass:: pydlpoly
    :members: set_extra_term, add_extra_system
          
Writing structures
------------------

.. autoclass:: pydlpoly
    :members: write_xyz, write_tinker_xyz, write_config

                                      
.. todo::
    
    Need more details here. Complete the docstrings including Paramters and Returns of all methods to be 
    reasonably called from users in the pydlpoly.py file
    
    
    
    


    
