import numpy as np
import ase
from atoms import system

def update(atoms, engine, dt):
    mass = atoms.get_eV_mass()           
    mass = mass.reshape((-1,1))

    f = engine(atoms) 
    p = atoms.get_velocities() * mass + dt/2 * f 

    x = atoms.get_positions() + dt/mass * p
    atoms.set_positions(x)

    f = engine(atoms)
    p = p + dt/2 * f
    atoms.set_velocities(p/mass)
    