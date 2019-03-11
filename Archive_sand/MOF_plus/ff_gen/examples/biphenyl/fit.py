#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#  Fitting script using the ffgenerator and the standard ric_fit objective, which
#  was also used in the old ffgen. This script uses CMA-ES as optimizer.
#
#  https://www.lri.fr/~hansen/cmaes_inmatlab.html#python
#
################################################################################

import pydlpoly
import molsys
import ff_gen.ff_gen as ff_gen
import ff_gen.objectives.ric_fit3 as ric_fit
from molsys.util import ff2pydlpoly
import cma

# specify a name
name = "ph-ph"

#  initialize a molsys object
m = molsys.mol.fromFile("ph-ph.mfpx")

# add ff addon and read in the fpar and ric file
m.addon("ff")
m.ff.read(name, fit= True)

# setup pydlpoly
wm = ff2pydlpoly.wrapper(m)
pd = pydlpoly.pydlpoly(name)
pd.setup(web=wm)

# initialize ric addon, need to be done after pydlpoly setup
m.addon("ric")
m.ric.setup_rics(full = False)

# initialize ric_fit objective
ric = ric_fit.ric_fit(only_diag = False)

# initialize FFgen
ff = ff_gen.ff_gen(name, pd, objective = ric, minimize = True)

# set maxiter and inititialize CMA-ES on the problem
miter = 10
es = cma.CMAEvolutionStrategy(ff.initials, 0.2, {'bounds': ff.get_bounds(), 
'maxiter':miter})

# optimize the objective function and write optimized parameters
es.optimize(ff.calc_objective,verb_disp=1)
ff.finish(es.result()[0])
