# A MD simulator

MD Simulator is written in python and the main dependencies are python packages Numpy and ASE. It is capable of performing molecular dynamics simulation for a collection of atoms with a given simulated ensemble and force field.

## Getting Started 

To run a simulation using MDS, the user should first initialize a *dynamics* object. The argument is the path of the input json file. 
```
my_simulation = mds.dynamics("./input.json")
my_simulation.run()
```

A sample input file is shown below:
```
{
    "mode": "nve",
    "init_temp": 94,

    "time_step": 1,
    "total_steps": 100000,

    "init_config_path": "./ar108.xyz",
    "ff_style": "lj",
    "ff_coeff": [0.01, 3.405, 8], 

    "log_freq": 100,
    "dump_freq": 100
}
```

### Prerequisites 

The main dependencies of MDS are ASE and numpy libraries, which are specified in the installation setup file. 
```
install_requires=[
    'ase>=3.21',
    'numpy>=1.14.5',
]
```

### Installation 

The user could install MDS using pip:
```
$ pip install .
```

### Usage example 

Please check examples in *./example/*. All the simulation results of sample problems are also stored in *backup_results* folder in the corresponding example's directory.
