import numpy as np
import ase
from atoms import system

def lj_force(r, epsilon, sigma, cutoff):
    if r <= cutoff and r > 0:
        sigma_r = sigma / r
        res = 24 * epsilon / r**2 * (2 * sigma_r**12 - sigma_r**6)
        return res
    else: 
        return 0


def lj(atoms, epsilon, sigma, cutoff=4):
    '''
    Args:
        atoms[ase.Atoms]
    Returns:
        forces[np.ndarray]
    '''
    def lj_engine(atoms):
        cell_size = np.diag(atoms.get_cell())
        if cutoff > cell_size.min()/2:
            raise NotImplementedError("cutoff error")

        forces = np.zeros_like(atoms.get_positions())
        N = len(atoms)
        for j in range(N):
            fij = [0,0,0]
            for i in range(N):
                rij = atoms.get_distances(i, j, mic=True, vector=True)[0]      # rj-ri
                r = abs(atoms.get_distances(i, j, mic=True, vector=False)[0])  # |rj-ri|
                fij += lj_force(r, epsilon, sigma, cutoff) * rij               # forces i to j
            forces[j] = fij * 1
        return forces

    return lj_engine
