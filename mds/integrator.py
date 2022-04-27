import numpy as np
import ase
from atoms import system

# force = eV/A, time = fs, length = A
def update(atoms, engine, dt):
    mass = atoms.get_masses()             # unit: u
    mass = mass * (1.66e-27 / 1.6022e-23) # unit: eV fs^2 / A^2
    mass = mass.reshape((5,1))

    f = engine(atoms) 
    p = atoms.get_velocities() * mass + dt/2 * f 

    x = atoms.get_positions() + dt/mass * p
    atoms.set_positions(x)

    f = engine(atoms)
    p = p + dt/2 * f
    atoms.set_velocities(p/mass)
    
    
    print("positions")
    print(atoms.get_positions())
    print("velocities")
    print(atoms.get_velocities())    
    print("force")
    print(f)