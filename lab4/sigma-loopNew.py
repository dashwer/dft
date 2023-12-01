#!/usr/bin/env python

import os
from ase import Atoms
from ase.build import bulk
from ase.lattice.cubic import BodyCenteredCubic
from ase.build import bcc100
from ase.calculators.vasp import Vasp

import numpy as np

from ase.units import J, m

current_mee = 'Pb'

currentdir = os.getcwd()
workdir = currentdir + '/' + current_mee + '/';

if not os.path.exists(workdir):
	os.makedirs(workdir)

energy_bulk = 0.0

for i in range(0, 11):
	current_work_dir = workdir + "/k_points=5/result-100/%i" % i

	if not os.path.exists(current_work_dir):
		os.makedirs(workdir + "/k_points=5/result-100/%i" % i)

	os.chdir(current_work_dir)

	if i == 0:
		atoms = bulk(current_mee)
		calc = Vasp(xc='PBE', kpts=(5,5,5), encut=500, ibrion=1, nsw=100, lvtot=True, lvhar=True, command='ASE_VASP_COMMAND')
		atoms.set_calculator(calc)
		energy_bulk = atoms.get_total_energy()
		print("process: " + str(i) + "\t" + "Energy Bulk" + str(energy_bulk) + "\n")

	else:
		atoms = bcc100(current_mee, size=(1,1,i), vacuum=10.0)
		atoms.center()
		calc = Vasp(xc='PBE', kpts=(5,5,1), encut=500, ibrion=1, nsw=100, lvtot=True, lvhar=True, command='ASE_VASP_COMMAND')

		atoms.set_calculator(calc)

		current_energy = atoms.get_total_energy()
		sigma = (current_energy - i * energy_bulk) * 0.5

		cell = atoms.get_cell()
		area = np.linalg.norm(np.cross(cell[0], cell[1]))
		sigma_aem = sigma/area/(J/m**2)

		print("process: " + str(i) + "\t" + str(current_energy) + "\t" + str(sigma) + "\t" + str(sigma_aem))

		with open(workdir+"log.txt", "w+") as logfile:
			logfile.write(str(current_energy) + "\t" + str(i) + str(sigma_aem) + "\n")
			logfile.close

	os.chdir(workdir)

