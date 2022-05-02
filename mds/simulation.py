

import json
import ase.io
from .atoms import system
from .integrator import update
from .force_fields import lj
from matplotlib import pyplot as plt
import logging
from time import time

class dynamics:

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
        if jdata["init_temp"]:
            self.init_temp = jdata["init_temp"]
            self.atoms.create_velocities(self.init_temp)

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
        self.log_freq = jdata["log_freq"]
        self.dump_freq = jdata["dump_freq"]

    def run(self):

        debug=True
        if debug:
            t_list = []
            d_list =[]

        logging.basicConfig(filename='md.log', encoding='utf-8', level=logging.INFO)
        logging.info('============================= SIMULATION =============================')
        logging.info('Step       PotEng          KinEng          TotEng          Volume')
        
        if self.ff == 'lj':
            force_engine, pot_engine = lj(self.atoms, self.epsilon, self.sigma, self.cutoff)
        else: 
            raise NotImplementedError('Only lj are supported')
        
        while self.step < self.total_steps:
            t0 = time()
            if debug:
                t_list.append(self.step * self.time_step)
                d_list.append(self.atoms.get_distance(0,1,mic=True))
            if self.step % self.log_freq == 0:
                ke = self.atoms.get_ke()
                pe = self.atoms.get_pe(pot_engine)
                vol = self.atoms.get_volume()
                logging.info('{:d}    {:e}    {:e}    {:e}    {:e}'.format(
                    self.step, pe, ke, pe+ke, vol ))

            if self.step % self.dump_freq == 0:
                ase.io.write('./traj/traj_{:d}.xyz'.format(self.step),self.atoms, format='extxyz')

            if self.thermostat:
                raise NotImplementedError
            else:
                update(self.atoms, force_engine, self.time_step)
            if self.step % 10 == 0:
                self.atoms.wrap()
                v = self.atoms.get_velocities()
                v_zeromean = v - v.mean(0)
                self.atoms.set_velocities(v_zeromean)
            
            tt = time() - t0
            self.step += 1
        if debug:
            plt.plot(t_list, d_list)
            plt.savefig('traj.png')
            print('used time :{}'.format(tt))


if __name__ == "__main__":
    import numpy as np
    np.set_printoptions(precision=5,suppress=True)
    new_sim = simulation("./example.json")
    new_sim.run()
