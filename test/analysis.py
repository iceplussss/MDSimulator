import ase
import numpy as np
import ase.io
from ase.geometry.analysis import Analysis
import matplotlib.pyplot as plt

num_frames = 350
num_bins = 100
throw = 50
lst = []

for i in range(throw, num_frames):
    atoms = ase.io.read("./traj/traj_{}.xyz".format(i*100), format='extxyz')
    x = Analysis(atoms).get_rdf(8, num_bins)[0]
    lst.append(x)

total_rdf = np.array(lst)
rdf = total_rdf.mean(0)

plt.plot(np.arange(num_bins), rdf)
plt.savefig("new_.png")
