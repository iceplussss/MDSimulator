"""An extended Class which contains information of the atoms, 
including their symbols, positions, masses, velocities. 
Also some "get functions" for calculation and output.
Unit: force = eV/A, time = fs, length = A
"""

import ase
import ase.io
import numpy as np


class system(ase.Atoms):

    def __init__(self, symbols, positions, cell, pbc=[True,True,True]):
        """Initializes a "system" object. 
        Args:
           symbols: an array of the atoms symbols  
           positions: an array of the atoms positions
           cell: an array giving the size of simulation box
           pbc: an optional array for boundary conditions of simulation box
        """
        super().__init__(symbols=symbols, positions=positions, cell=cell, pbc=pbc)


    def get_eV_mass(self):
        """Calculates the masses in unit [eV fs^2 / A^2].
        Converts masses from [u] to [eV fs^2 / A^2].
        """

        mass = self.get_masses()              
        mass = mass * (1.66e-27 / 1.6022e-29) 
        return mass


    def create_velocities(self, temp):
        """Creates initial velocities for given temperature by Gaussian  
        """

        kb = 8.617333262e-5 # in eV
        m = self.get_eV_mass()  
        m = m.reshape((-1,1)) 
        N = self.get_positions().shape[0] 
        v = np.random.normal(0, 1, N*3).reshape(N,3)
        self.set_velocities( v * (kb * temp / m)**0.5 )


    def get_ke(self):
        """Calculates the total kinetic energy."""

        m = self.get_eV_mass()  
        m = m.reshape((-1,1)) 
        v = self.get_velocities()
        ke = 1/2 * m * v**2
        total_ke = ke.sum()
        return total_ke


    def get_pe(self, pot_engine):
        """Calculates the total potential energy."""

        total_pe = pot_engine(self)
        return total_pe