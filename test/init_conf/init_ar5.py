"""Creates initial configuration xyz file for 5 Ar atoms
They are placed equally spaced but with a small perturbation
"""
import numpy as np
import mds


r = 3.82
dr = 0.001
pos = [ [0, 0, 0.0 + dr],
        [0, 0, r],
        [0, 0, r * 2],
        [0, 0, r * 3],
        [0, 0, r * 4] ]
cell = np.eye(3) * r * 5

atoms = mds.system(['Ar','Ar','Ar','Ar','Ar'], pos, cell, pbc=True)


print("N =", len(atoms))
print(atoms.get_cell())
print(atoms.get_masses())
print(atoms.get_eV_mass())

# ase.io.write("init.xyz", atoms, format="extxyz")