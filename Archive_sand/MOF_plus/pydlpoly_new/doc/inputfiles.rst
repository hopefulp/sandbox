.. _inputfiles:

**********************************
Format of pydlpoly's Inputfiles
**********************************

The native input files ``CONFIG``, ``FIELD`` and ``CONTROL`` read by the core ``DL_Poly`` engine are still produced on the fly by the Python frontend and
can be checked in the temporary rundirectory. However, the main input files are the <key> and <xyz> file, which contain the force field definition and 
the geometry/connectivity. The are inspired by the Tinker format, but especially the <key> file format is not exactly equal to Tinker and not compatible to it.

.. _xyz-files:

The Geometry file (extension .xyz)
==================================

Here are the first lines of the file mof5.xyz describing the periodic MOF called MOF-5::
    
      424    26.0699    26.0699    26.0699    90.0000    90.0000    90.0000
      1 o      0.003061   -0.000459   -0.000450   165      2      3      4      5 
      2 zn     1.139801   -1.137243   -1.137311   166      1     14     18     21 
      3 zn    -1.133735   -1.137131    1.136441   166      1      6     15     23 
      4 zn    -1.134000    1.136288   -1.136999   166      1     11     12     17 
      5 zn     1.140090    1.136351    1.136079   166      1      8      9     20 
      6 o     -0.796133   -0.799964    3.041162   167      3      7 
      7 c      0.003076    0.000224    3.621042   168      6      8    230 
      8 o      0.802262    0.800295    3.040974   167      5      7 
      9 o      0.804535    3.041271    0.797591   167      5     10 
     10 c      0.003073    3.621184   -0.000322   168      9     11    210 
     11 o     -0.798196    3.041106   -0.798310   167      4     10 
     .
     .
     
Appart from the unit cell parameter information in the first line, this file
is exactly identical to the regular Tinker xyz format and can be opened by ``molden`` 
or ``vmd``.

First line
----------

* The first integer defines the number of atoms (and thus following lines in the file)
* The next three floats define the a, b and c axis length of the unit cell in Angstrom.
* The last three floats define the cell angles alpha, beta and gamma in degree.

If the system is non-periodic the cell parmaters must be omitted (Warning: ``molden`` writes
Tinker xyz files with some text after the number of atoms. You need to remove this before reading with
``pydlpoly``).

Further lines
------------- 
* The first integer number is just a running index starting from 1! (not really clear why but we kept it to keep the format readable by ``molden``)
* The second entry is the element
* The next three floats (col 3-5) are the xyz cartesian coordinates in Angstrom
* The next entry (col 6) is the atom type. ``pydlpoly`` handles these as strings and you can have any name (like "Zn_MOF5" or so). But for compatibility with
  ``molden`` or ``vmd`` the atom types must be numbers as required in Tinker. Alternatively, an ``#atomtypes`` flag can be appended at the end of the xyz file 
  to provide a list of strings, numbers and elements and allow ``pydlpoly`` to work with a mixed input from xyz- and :ref:`key-files`. Example::
      #atomtypes
      165 O_Zn_MOF5 o
      166 Zn_MOF5 zn
      167 O_C_Zn_MOF5 o
      168 C_MOF5 c
      
* All further integers are the indices of the connected atoms (defined back and forth .. so 1 binds to 2 and 2 to 1)

.. _key-files:

The Force Field file (extension .key)
=====================================

The force field parameters are defined in the key file, which is derived from the original Tinker
format. More recently we have abandonned the backwards compatibility to Tinker and thus some 
obsolete information was discarded from the key files. We refer to these files as "version 2.0" files.
The script ``convert_key`` can be used to convert the old file format into the new format. Note that the current ``pydlpoly`` can 
**only** read the new format.

The original Tinker key-file format is described in detail on the `Tinker website <http://dasher.wustl.edu/tinkerwiki/index.php/Main_Page>`_.

The pydlpoly specific key file format is described below. The file is organized in lines starting with a keyword. Comment
lines start with a hash tag "#". There are settings keywords which should appear only once (last occurance will overwrite 
previous settings) and potential term keys. The latter define a specific potential type. Potential type keywords
are always followed by one or more potential type specifiers. Originally these potential types were integer numbers (derived
from the "ancient" MM3 atom types). However, the potential types do not need to be integer numbers but can also be 
strings (see also :ref:`xyz-files`) for better readability. The drawback is that in this case the written xyz files are not
readable by ``molden``. This problem can be circumvented by setting the ``moldenr`` flag of the
:class:`pydlpoly.pydlpoly.write_tinker_xyz` method to ``True``, which automatically replaces the strings with numbers.
For historic reasons, the examples shown below still use integer atom types.

Here is an example for the ``MOF-FF`` (see this `paper <http://dx.doi.org/10.1002/pssb.201248460>`_) force field 
for MOF-5::

    version        2.0
    parameters     none

    bondunit       71.94
    angleunit      0.02191418
    strbndunit     2.51118
    opbendunit     0.02191418
    torsionunit    0.5
    vdwtype        exp6_damped
    vdwdampfact    0.25
    radiusrule     arithmetic
    radiustype     r-min
    radiussize     radius
    epsilonrule    geometric
    a-expterm      184000.0
    b-expterm      12.0
    c-expterm      2.25
    bondtype       mixmorse_bde
    strbndtype     mmff
    opbendtype     mmff
    chargetype     gaussian
    
    # params for phenyl group
    atom           5         h
    atom           2         c
    
    vdw            5               1.5000      0.0200
    vdw            2               1.9600      0.0560
    
    bond           2         2               7.0800      1.3940
    bond           2         5               5.4300      1.0940
    
    angle          2         2         2               0.7410    127.0500
    angle          2         2         5               0.5030    120.3500
    
    torsion        2         2         2         2               0.0000      4.3790      0.0000
    torsion        5         2         2         5               0.0000      5.9720      0.0000
    torsion        2         2         2         5               0.0000      6.3160      0.0000
    
    opbend         2         2         5         2               0.0190
    
    strbnd         2         2         2               0.0470      0.0470      0.4990
    strbnd         2         2         5              -0.1750      0.3720      0.6490
        
    charge         5               0.1200      0.7236
    charge         2              -0.1200      1.1630

    # params for Zn4O unit
    atom           165       o
    atom           166       zn
    atom           167       o
    atom           168       c
    
    vdw            165             1.8200      0.0590
    vdw            166             2.2900      0.2760
    vdw            167             1.8200      0.0590
    vdw            168             1.9400      0.0560
    
    bond           165       166             1.4890      1.9870     50.0000
    bond           166       167             1.6650      1.9170     50.0000
    bond           168       167             8.6320      1.2750
    bond           2         168             4.9370      1.4880
    
    angle          166       165       166             0.6980    103.9920
    angle          165       166       167             0.0000    113.5840
    angle          167       166       167             0.0800    123.1030
    angle          166       167       168             0.0940    135.6060
    angle          167       168       167             1.5550    123.0060
    angle          2         168       167             1.0730    116.3680
    angle          2         2         168             0.8060    117.2960
    
    torsion        166       165       166       167             0.0000      0.0000      0.0000
    torsion        165       166       167       168             0.0000      0.0000      0.0000
    torsion        166       167       168       167             0.0000      0.0760      0.0000
    torsion        166       167       168       2               0.0000      3.0120      0.0000
    torsion        167       166       167       168             0.0000      0.0000      0.0000
    torsion        2         2         168       167             0.0000      1.9020      0.0000
    torsion        2         2         2         168             0.0000      0.0000      0.0000
    torsion        5         2         2         168             0.0000      0.0000      0.0000
    
    opbend         168       167       2         167             0.1880
    opbend         2         168       2         2               0.0880
    
    charge         165            -1.7800      1.1176
    charge         166             1.4200      2.0733
    charge         167            -0.7200      1.1176
    charge         168             0.6100      1.1630
    chargemod      2         168             0.1800      1.1630
    
    strbnd         166       165       166             0.1280      0.1280      0.0530
    strbnd         165       166       167             0.1960      0.0400     -0.1610
    
Settings Keywords
-----------------

.. note::
    
    A remark on the units used in pydlpoly: all energies are given in **kcal/mol**. ``Dl_Poly`` allows differnt
    energies to be used but by default a ``kcal`` is written to the ``FIELD`` file. Note that ``DL_Poly`` internally 
    works with a unit of 10 J/mol and internal energies are converted to kcal/mol in ``pydlpoly`` using ``dlp.engunit``.
    However, for historic reasons our key files use the MM3 convention for force constants (as Tinker usually does) using mdyne/A.

* bondunit:    Factor to convert force constants in bond keyword lines to kcal/molA^2 (also contains 0.5!!)
* angleunit:   dito for angles
* strbndunit:  dito for strbnd terms
* opbendunit:  dito for out -of-plane bend terms
* torsionunit: dito for torsions (we use barriers in kcal/mol as in Tinker so factor is 0.5)

.. note::
    
    These factors also depend on the (partly strange) definition of potential terms in Tinker (we have adopted these historically
    for compatibility). Please check the numbers in the ``FIELD`` file in case you are not sure what is really used in the calcualtions.
    
* vdwtype: defines the type of vdw term used (``exp6_damped``: MOF-FF specific modified Buckingham, ``buckingham`` or ``lennard-jones``)
* vdwdampfact: only needed in case of ``exp6_damped``
* radiusrule: combination rule (either ``arithmetic`` or ``geometric``)
* radiustype: either ``r-min`` or ``sigma``
* radiussize: either ``radius`` or ``diameter``
* epsilonrule: combination rule for epsilon (either ``arithmetic`` or ``geometric``)
* a-term/b-term/c-term: specific params for the MM3 type ``buckingham`` potential (also valid for ``exp6_damped``)
* bondtype: can be omitted, if ``mixmorse`` or ``mixmorse_bde`` is used then the morse potential is employed for any bond term with three instead of two paramters (see below)
* strbndtype: for MOF-FF request ``mmff`` type
* opdendtype: for MOF-FF request ``mmff`` type
* chargetype: using ``gaussian`` here switches to Gaussian type charge distributions (changes format of ``charge`` and ``chargemod`` potential keyword format)
* rigid: list of molecules (defined e.g. via molname) to be kept as rigid units (works only for MD!)
* freeze: list of molecules to be kept frozen


Potential Keywords
------------------

This table is giving an overview about potential keywords (they can show up multiple times).
please see the comments at the end of the table and in particular the examples in the next subsection.

=================  =======  ===========  ===========  ===========  ===========  ==================
Keyword            # atoms  Parameter 1  Parameter 2  Parameter 3  Parameter 4  Comments
=================  =======  ===========  ===========  ===========  ===========  ==================
atom               1        element                                             defines a type
vdw                1        radius       epsilon                                use combination rules
vdwpr              2        radius       epsilon                                explicit params for a pair
bond               2        forceconst   refdist      [alph/bde]                mixmorse: use morse pot if third paramter is given.
bond5              2
angle              3        forceconst   refangle                               
angle5             3
anglef             3
anglef-2           3
torsion            4        V1           V2           V3           [V4]         fourth is optional
torsion5           4
opbend             4        forceconst                                          first atom is central
strbnd             3        strbnd1      strbnd2      strstr                    refval from bond and angle
charge             1        q            [sigma]                                sigma for gaussians
chargemod          2        q            [sigma]                                modifies charge
virtbond           2                                                            no params
molname            n        name         atomtypes                              name is a string!                                                    
restrain-distance  2
restrain-angle     3
restrain-torsion   3
equivalent         2                                                            atom types are equal
=================  =======  ===========  ===========  ===========  ===========  ==================

* aph/bde :  for ``mixmorse`` the third param triggers the use of a morse potential with the given alpha. for ``mixmorse_bde`` the third param is the bond dissociation energy in kcal/mol
* torsion : barriers of torsions with n= 1, 2, 3 and optional 4 in kcal/mol
* charge  : if atomtype is a **negative** integer it refers to the atom index (as in Tinker)
* chargemod : modifies charge if atom is bonded to a specfic other atom type
* virtbond: one of the two atoms must be a virtual atom (element "Xx"), no params needed
* molname: can have an arbitrary number of atomtypes and defines the name for such a system

Examples
--------

- **chargemod**: Assume "Car" is the atomtype of an aromatic system which is usually -0.12. However, if "Car" is the alpha-carbon
  of a carboxylate (bonded to a "Ccarbox") you want it to have the same params as all other "Car" but its charge should be 0.0.
  For pointn charges (no sigma you can define this accordingly::
      
      charge     Har             0.12
      charge     Car            -0.12
      chargemod  Car   Ccarbox   0.0
      
- **molname**: You want to have all benzene molecules to be named "bz". The are built from atomtypes "Cbz" and "Hbz"::
    
    molename bz  Hbz  Cbz   
    
- **equivalent**: In your keyfile you have params for Car and Har and you want to "reuse" them with another name, for example to 
  define a molecule name with it or to keep all params and just change a specific thing (the charge or a bond term etc.) you can 
  make two types equivalent. The first atomtype is the new "alias" and the second mus be defined already. As long as no specific 
  potential type for the new "alias" is found it will be replaced by the already existing params. This means an alternative way
  of defining the charges not using chargemod whould be this::
         
     charge     Har             0.12 
     charge     Car            -0.12
     equivalent Caralp   Car
     charge     Caralp          0.0
        
  For the above definition of the benzen molname one could use the following equivalence::
        
     equivalent   Hbz   Har
     equivalent   Cbz   Har
        
        
      

The obsolete CONTROL File
=========================

.. currentmodule:: pydlpoly

A number of paramters in the original ``DL_Poly`` are set in the ``CONTROL`` file. This is not needed, but can be provided.
If ``control=<Filename>`` is used as a named paramter in :class:`pydlpoly.setup` then instead of the default ``CONTROL`` file 
the provide file named <Filename> is used. However, a large number of settings in the ``CONTROL`` file are no longer read and used.
For example the number of MD steps or the temperature or pressure, which is directly set in :class:`pydlpoly.MD_init`.
On the other hand a number of settings like the cutoff radius or the specific way to compute the electrostatics can still be changed via the ``CONTROL`` files.
The most direct way to do this is via the ``control`` directory structure of the class :class:`pydlpoly.pydlpoly`.
Before calling setup (this is important, since during setup the ``CONTROL`` file is generated) just change or add settings to this directory. Here is an example
how to reduce the default cutoff from 12.0 to 10.0 A.::
    
    import pydlpoly
    
    pd = pydlpoly.pydlpoly("test")
    pd.control["cut"] = 10.0
    pd.setup()
    
It is recommended to check the final ``CONTROL`` file. For details about the settings see the ``DL_Poly`` documentation.

.. todo::
    
    We need a list of CONTROL keywords which can be changed.
    
    


Naming of Molecules
===================


Restarting
==========




 
 