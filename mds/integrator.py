"""Functions used to integrate the Newton's equations of motion.
The velocity Verlet algorithm is implemented: Time reversibility 
and preservation of the symplectic form on phase space are kept.
"""

import numpy as np
import ase

from .atoms import system


def update(atoms, engine, dt):
    """Initialises Simulation class.
    Args:
        atoms: a system object describing its current state
        engine:  force engine for the simulation
        dt: timestep for the simulation
    """
    
    mass = atoms.get_eV_mass()           
    mass = mass.reshape((-1,1))

    f = engine(atoms) 
    p = atoms.get_velocities() * mass + dt/2 * f 

    x = atoms.get_positions() + dt/mass * p
    atoms.set_positions(x)

    f = engine(atoms)
    p = p + dt/2 * f
    atoms.set_velocities(p/mass)
    