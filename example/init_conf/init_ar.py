"""Produces xyz configurations for different size of system."""

import ase.io

# init_xyz = ase.io.read("ar.conf", format='lammps-dump-text')
# # atoms = system(init_xyz.get_chemical_symbols(), init_xyz.get_positions(), init_xyz.get_cell())
# ase.io.write('ar108.xyz', init_xyz, format='extxyz')

x = ase.io.read('ar108.xyz', format='extxyz')

y = x.repeat((2,1,1))
ase.io.write('ar216.xyz', y, format='extxyz')

y = x.repeat((2,2,1))
ase.io.write('ar432.xyz', y, format='extxyz')

y = x.repeat((2,2,2))
ase.io.write('ar864.xyz', y, format='extxyz')