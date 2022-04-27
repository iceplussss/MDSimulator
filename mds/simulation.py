import json
from atoms import system
from integrator import update
import ase.io
from force_fields import lj
from matplotlib import pyplot as plt
class simulation:

    def __init__(self, json_path):

        with open(json_path, 'r') as file:
            jdata = json.load(file)

        # Configuration
        init_xyz = ase.io.read(jdata["init_config_path"], format='extxyz')
        self.atoms = system(init_xyz.get_chemical_symbols(), init_xyz.get_positions(), init_xyz.get_cell())

        # Ensemble Type
        if jdata['mode'] == 'nve':
            self.thermostat = False
        elif jdata['mode'] == 'nvt':
            self.thermostat = True
            self.temp = jdata["temp"]
        else:
            raise NotImplementedError('Only NVE and NVT ensembles are supported')

        # Force Fields
        if jdata['ff_style'] == 'lj':
            self.ff = 'lj'
            self.epsilon = jdata['ff_coeff'][0]
            self.sigma = jdata['ff_coeff'][1]
            self.cutoff = jdata['ff_coeff'][2]

        # Simulation Setup
        self.total_steps = jdata["total_steps"]
        self.time_step = jdata["time_step"]
        self.total_time = self.time_step * self.total_steps
        self.step = 0

        # Output Setup
        self.print_freq = jdata["print_freq"]

    def run(self):
        if self.ff == 'lj':
            engine = lj(self.atoms, self.epsilon, self.sigma, self.cutoff)
        else: 
            raise NotImplementedError('Only lj are supported')
        debug=True
        if debug:
            t_list = []
            d_list =[]
        while self.step < self.total_steps:
            if debug:
                t_list.append(self.step * self.time_step)
                d_list.append(self.atoms.get_distance(0,1,mic=True))
            if self.step % self.print_freq == 0:
                # ase.io.write()
                pass
            if self.thermostat:
                raise NotImplementedError
            else:
                update(self.atoms, engine, self.time_step)
            if self.step % 100 == 0:
                self.atoms.wrap()
            

            self.step += 1
        if debug:
            plt.plot(t_list, d_list)
            plt.savefig('traj.png')


if __name__ == "__main__":
    import numpy as np
    np.set_printoptions(precision=5,suppress=True)
    new_sim = simulation("./example.json")
    new_sim.run()
