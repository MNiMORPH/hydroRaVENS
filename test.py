# Program below the class
import numpy as np
from numpy.random import poisson
from matplotlib import pyplot as plt
import bucketHydrology as bh
#plt.ion()

rain = poisson(.2, 100) # Convert to an import for real data
dt = 1. # day
        
# Change these parameters with an input file or script, eventually
# Arbitrary units; will be made real in an actual landscape
res_surface = bh.reservoir(t_efold=1., f_to_discharge=0.5, Hmax=10.)
res_deep = bh.reservoir(t_efold=10., f_to_discharge=1., Hmax=np.inf)

strat_column = [res_surface, res_deep]

watershed = bh.buckets(reservoir_list=strat_column, dt=dt)
watershed.run(rain)
watershed.plot()
