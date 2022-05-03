import numpy as np
import mds

np.set_printoptions(precision=5, suppress=True)
new_sim = mds.dynamics("./ar108.json")
new_sim.run()
