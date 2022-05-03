from time import time
import numpy as np
import mds

np.set_printoptions(precision=5, suppress=True)

Nlist = [108, 216, 432, 864]
for N in Nlist:
    new_sim = mds.dynamics("./ar{}.json".format(N))

    t0 = time()
    new_sim.run()
    t1 = time()

    print("Time used for {} Ar:".format(N), t1-t0)