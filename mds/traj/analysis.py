import ase
import numpy as np
import ase.io
from ase.geometry.analysis import Analysis
import matplotlib.pyplot as plt

num_frames = 100
lst = np.zeros(100)

for i in range(num_frames):
    atoms = ase.io.read("./traj_{}.xyz".format(i*100), format='extxyz')
    x = Analysis(atoms).get_rdf(8, 100)[0]
    lst += x

rdf = lst/num_frames
print(np.shape(rdf))

plt.plot(np.arange(100), rdf)
plt.savefig("new.png")
