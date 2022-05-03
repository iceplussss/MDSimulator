"""A class which reads the setting, runs the simulation and outputs results.
"""

import json
import ase.io
import logging
import os

from .atoms import system
from .integrator import update
from .forcefields import lj


class dynamics:

    def __init__(self, json_path):
        """Reads an json file and initialize a "dynamics" object.
        Args:
            json_path: path to find the json file
        """

        with open(json_path, 'r') as file:
            jdata = json.load(file)

        # Initial Configuration
        init_xyz = ase.io.read(jdata["init_config_path"], format='extxyz')
        self.atoms = system(init_xyz.get_chemical_symbols(), init_xyz.get_positions(), init_xyz.get_cell())

        # Ensemble Info
        if jdata["init_temp"]:
            self.init_temp = jdata["init_temp"]
            self.atoms.create_velocities(self.init_temp)
        if jdata['mode'] == 'nve':
            self.mode = 'nve'
        else:
            raise NotImplementedError('Only NVE ensemble is supported')
        
        # Force Fields
        if jdata['ff_style'] == 'lj':
            self.ff = 'lj'
            self.epsilon = jdata['ff_coeff'][0]
            self.sigma = jdata['ff_coeff'][1]
            self.cutoff = jdata['ff_coeff'][2]
        else:
            raise NotImplementedError('Only LJ potential is supported')

        # Simulation Setup
        self.total_steps = jdata["total_steps"]
        self.time_step = jdata["time_step"]
        self.total_time = self.time_step * self.total_steps
        self.step = 0

        # Output Setup
        self.log_freq = jdata["log_freq"]
        self.dump_freq = jdata["dump_freq"]

    def run(self):
        """Function which runs the simulation.
        Does the whole simulation and performs output (logging and dumping)

        logging: Prints info (step, potential energy, kinetic energy, total energy and volume)
        according to log frequency to a single "md.log" file
        dumping: Prints system info (atomic configuration, velosities and cell size) 
        of current time step according to dump frequency to a "traj" folder 
        """

        logging.basicConfig(filename='md.log', encoding='utf-8', level=logging.INFO)
        logging.info('============================= SIMULATION =============================')
        logging.info('Step       PotEng          KinEng          TotEng          Volume')
        
        if os.path.exists('./traj') is False:
            os.mkdir('traj')

        if self.ff == 'lj':
            force_engine, pot_engine = lj(self.atoms, self.epsilon, self.sigma, self.cutoff)
        
        while self.step < self.total_steps:
            # Logging output info 
            if self.step % self.log_freq == 0:
                ke = self.atoms.get_ke()
                pe = self.atoms.get_pe(pot_engine)
                vol = self.atoms.get_volume()
                logging.info('{:d}    {:e}    {:e}    {:e}    {:e}'.format(
                    self.step, pe, ke, pe+ke, vol ))

            # Dumping trajectories 
            if self.step % self.dump_freq == 0:
                ase.io.write('./traj/traj_{:d}.xyz'.format(self.step),self.atoms, format='extxyz')

            # Integrate one step
            update(self.atoms, force_engine, self.time_step)
            if self.step % 10 == 0:
                self.atoms.wrap()
                v = self.atoms.get_velocities()
                v_zeromean = v - v.mean(0)
                self.atoms.set_velocities(v_zeromean)
            self.step += 1
