import os, sys, copy
import numpy as np
import random

import CNT
import bgf, bgftools
import nutils as nu

version = "150204"

class Graphene(CNT.Nanotube):
    def __init__(self, *args, **kwargs):
        super(Graphene, self).__init__()
        print "successfully loaded an empty graphene object."
