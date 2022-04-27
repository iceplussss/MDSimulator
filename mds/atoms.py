import ase
import numpy as np
import ase.io

class system(ase.Atoms):
    def __init__(self, symbols, positions, cell, pbc=[True,True,True]):
        super().__init__(symbols=symbols, positions=positions, cell=cell, pbc=pbc)


if __name__ == "__main__":
    r = 3.82
    dr = 0.001
    pos = [ [0, 0, 0.0 + dr],
            [0, 0, r],
            [0, 0, r*2],
            [0, 0, r*3],
            [0, 0, r*4]]
    cell = np.eye(3) * r*5
    atoms = system(['Ar','Ar','Ar','Ar','Ar'], pos, cell, pbc=True)
    # print(len(atoms))
    # print(atoms.get_cell())
    # print(atoms.get_velocities())
    print(atoms.get_masses())
    ase.io.write("init.xyz", atoms, format="extxyz")

    # print(atoms.get_all_distances(mic=True, vector=True))
    # print(atoms.get_distances(0, 1, mic=True, vector=False)[0])
    # fij = atoms.get_distances(0, 1, mic=True, vector=True)[0]
    # fij = atoms.get_distances(1, 0, mic=True, vector=True)[0]
    # print(fij)
    # print(fij + [0,0,1])

    # positions = atoms.get_positions()
    # forces = positions*0
    # print(forces)
    # forces[0] = [1,2,3]
    # print(forces)



