#!/usr/bin/env python

#from jasp import *
import os
from ase import  Atoms
from ase.build import  bulk
from ase.lattice.cubic import BodyCenteredCubic
from ase.build import bcc100
from ase.calculators.vasp import Vasp

import numpy as np

from ase.units import J, m


current_mee='Nb'

atoms=bulk(current_mee)
print(atoms.cell)
