import numpy as np
import ase
from .atoms import system


def lj(atoms, epsilon, sigma, cutoff=4):
    '''
    '''
    def lj_force_engine_slow(atoms):
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

    def lj_force_engine(atoms):
        cell_size = np.diag(atoms.get_cell())
        if cutoff > cell_size.min()/2:
            raise NotImplementedError("cutoff error")

        forces = np.zeros_like(atoms.get_positions())
        N = len(atoms)
        D = atoms.get_all_distances(mic=True, vector=True)
        for i in range(N):
            rij = D[i]   ## (n_pair, 3)
            r = ((rij**2).sum(-1))**0.5   ## (n_pair)
            rfilter = ((r <= cutoff) & (r > 0))  ## apply cutoff
            rij = rij[rfilter]
            r = r[rfilter]
            sigma_r_6 = (sigma / r)**6
            res = 24 * epsilon / r**2 * (2 * sigma_r_6**2 - sigma_r_6) # (n_neighbor)
            forces[i] = (res.reshape(-1,1)  * (-rij)).sum(0)
        return forces

    def lj_pot_engine(atoms):
        cell_size = np.diag(atoms.get_cell())
        if cutoff > cell_size.min()/2:
            raise NotImplementedError("cutoff error")

        pot = 0
        N = len(atoms)
        for j in range(N):
            for i in range(N):
                if i == j:
                    continue
                r = abs(atoms.get_distances(i, j, mic=True, vector=False)[0])  # |rj-ri|
                pot += lj_energy(r, epsilon, sigma, cutoff)          
        pot = pot/2 # double counting
        return pot

    return lj_force_engine, lj_pot_engine

def lj_force(r, epsilon, sigma, cutoff):
    if r <= cutoff and r > 0:
        sigma_r = sigma / r
        res = 24 * epsilon / r**2 * (2 * sigma_r**12 - sigma_r**6)
        return res
    else: 
        return 0

def lj_energy(r, epsilon, sigma, cutoff):
    if r <= cutoff and r > 0:
        sigma_r = sigma / r
        res = 4 * epsilon * (sigma_r**12 - sigma_r**6)
        return res
    else: 
        return 0