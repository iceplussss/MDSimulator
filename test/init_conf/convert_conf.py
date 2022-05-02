"""Reads a lammps config file (108 Ar) and save its xyz version."""

import ase.io

init_xyz = ase.io.read("ar.conf", format='lammps-dump-text')
# atoms = system(init_xyz.get_chemical_symbols(), init_xyz.get_positions(), init_xyz.get_cell())
ase.io.write('ar108.xyz', init_xyz, format='extxyz')
