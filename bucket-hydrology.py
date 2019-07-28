#! /usr/bin/env python

# Started by A. Wickert
# 25 July 2019

class reservoir(object):
    """
    Generic reservoir. Accepts new water (recharge), and sends it to other
    reservoirs and/or out of the system (discharge) at a rate that is 
    proportional to the amount of water held in the reservoir.
    """

    import numpy as np

    def __init__(self, t, Vmax=np.inf):
        self.Vwater = 0.
        self.Vmax = Vmax
        self.t_efold = t
        self.excess = 0.
        self.Vout = np.nan
    
    def recharge(self, V):
        if self.Vwater+V <= self.Vmax:
            excess = 0.
            self.Vwater += V
        if self.Vwater+V > self.Vmax:
            excess = self.Vwater+V - self.Vmax
            self.Vwater = self.Vmax
        self.excess += excess

    def discharge(self, dt):
        dV = self.Vwater * (1 - np.exp(-dt/self.t_efold))
        self.Vout = dV + self.excess
        self.Vwater -= dV
        self.excess = 0.


# Program below the class

import numpy as np
from numpy.random import poisson
from matplotlib import pyplot as plt
#plt.ion()

rain = poisson(.2, 100) # Convert to an import for real data

dt = 1. # day
        

# Change these parameters with an input file or script, eventually
# Arbitrary units; will be made real in an actual landscape
res_tile = reservoir(0.1, 0.)
res_surface = reservoir(1., 10.)
res_deep = reservoir(10., np.inf)

# Change this parameter with a script as well
f_tile = 0.4 # fraction of land surface covered in drain tile

# Partitioning -- thse should also be external parameters
f_tile_to_discharge = 0.9
f_surface_to_discharge = 0.5


Q = []
# This includes made-up rules about how much of the discharge from each
# bucket goes to the other buckets or to the discharge at the outlet.
# I also use arbitarry units
for ti in range(len(rain)):
    Qi = 0.
    # Tile 
    res_tile.recharge(f_tile * rain[ti])
    res_tile.discharge(dt)
    Qi += f_tile_to_discharge * res_tile.Vout
    # Surface
    res_surface.recharge( (1.-f_tile) * rain[ti] + 
                          (1. - f_tile_to_discharge) * res_tile.Vout )
    res_surface.discharge(dt)
    Qi += f_surface_to_discharge * res_surface.Vout
    # Deep
    res_deep.recharge( (1. - f_surface_to_discharge) * res_surface.Vout )
    res_deep.discharge(dt)
    Qi += res_deep.Vout
    # time series
    Q.append(Qi)
    Qi = 0

plt.figure()
plt.plot(rain, 'g', label='rainfall')
plt.plot(Q, 'b', label='discharge')
plt.legend(fontsize=11)
plt.title('Tile-drained fraction: '+'%.1f' %f_tile)
plt.ylabel('Rainfall or discharge [units arbitrary]', fontsize=14)
plt.xlabel('Time [units arbitrary]', fontsize=14)
#plt.savefig('TDF-'+'%.1f' %f_tile+'.png')
plt.show()

