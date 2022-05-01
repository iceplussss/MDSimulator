import ase
import numpy as np
import ase.io

# Unit: force = eV/A, time = fs, length = A
class system(ase.Atoms):

    def __init__(self, symbols, positions, cell, pbc=[True,True,True]):
        super().__init__(symbols=symbols, positions=positions, cell=cell, pbc=pbc)

    def get_eV_mass(self):
        mass = self.get_masses()              # unit: u
        mass = mass * (1.66e-27 / 1.6022e-23) # unit: eV fs^2 / A^2
        return mass

    def get_ke(self):
        m = self.get_eV_mass()  
        m = m.reshape((-1,1)) 
        v = self.get_velocities()
        ke = 1/2 * m * v**2
        total_ke = ke.sum()
        return total_ke

    def get_pe(self, pot_engine):
        total_pe = pot_engine(self)
        return total_pe


if __name__ == "__main__":
    r = 3.82
    dr = 0.001
    pos = [ [0, 0, 0.0 + dr],
            [0, 0, r],
            [0, 0, r * 2],
            [0, 0, r * 3],
            [0, 0, r * 4]]
    cell = np.eye(3) * r * 5
    atoms = system(['Ar','Ar','Ar','Ar','Ar'], pos, cell, pbc=True)
    # print(len(atoms))
    # print(atoms.get_cell())
    print(atoms.get_masses())
    print(atoms.get_eV_mass())
    # ase.io.write("init.xyz", atoms, format="extxyz")