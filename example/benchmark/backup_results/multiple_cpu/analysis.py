import ase
import numpy as np
import ase.io
from ase.geometry.analysis import Analysis
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['lines.markersize'] = 5
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
plt.rcParams.update({'font.size':15})

lst_cpus = [2, 4, 10, 20]
lst_size = [108, 216, 432, 864]
lst_time = []

lst_time.append([5.432769775390625, 33.95189571380615, 145.53084325790405, 304.01843214035034]) 
lst_time.append([5.828061103820801, 34.62122392654419, 144.73593854904175, 315.83761954307556])
lst_time.append([5.158254384994507, 31.35063409805298, 131.58981609344482, 272.59346532821655])
lst_time.append([5.484041452407837, 31.26741647720337, 129.33909106254578, 267.62385296821594])

xp = np.linspace(100, 1000, 100)

fig, ax = plt.subplots(figsize=(8,6),dpi=200)

ax.set_xlabel('Number of atoms in the system')
ax.set_ylabel('Wall time [s]')

for i,cpu in enumerate(lst_cpus):
    z = np.polyfit(lst_size, lst_time[i], 1)
    ax.plot(lst_size, lst_time[i], '-o', label="{} cpus".format(cpu))
    ax.plot(xp, z[0] * xp + z[1], ':')

ax.legend()

plt.savefig("benchmark.png")