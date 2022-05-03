import ase
import numpy as np
import ase.io
from ase.geometry.analysis import Analysis
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['lines.markersize'] = 2
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
plt.rcParams.update({'font.size':15})

sigma = 3.405
num_bins = 100
num_frames = 1000
throw = 100

lst = []
for i in range(throw, num_frames):
    atoms = ase.io.read("./traj/traj_{}.xyz".format(i*100), format='extxyz')
    x = Analysis(atoms).get_rdf(8, num_bins)[0]
    lst.append(x)
total_rdf = np.array(lst)
rdf = total_rdf.mean(0)

data = np.loadtxt("1964data.csv", dtype=float,delimiter=',')
# print(data)

fig, ax = plt.subplots(figsize=(8,6),dpi=200)

ax.set_xlabel(r'$r/\sigma$')
ax.set_ylabel('g(r)')
ax.plot(np.arange(num_bins)/num_bins * 8 / sigma, rdf, label="NVE")
ax.scatter(data[:,0], data[:,1], label="NVT (A.Rahman 1964)", color='red')
ax.legend()
ax.grid(linestyle = ':')
ax.set_ylim(-0.2, 3)

plt.savefig("RDF.png")