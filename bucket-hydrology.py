#! /usr/bin/env python

# Started by A. Wickert
# 25 July 2019
# Updated by J. Jones
# Starting 08 Oct 2019

import numpy as np
from matplotlib import pyplot as plt

class reservoir(object):
    """
    Generic reservoir. Accepts new water (recharge), and sends it to other
    reservoirs and/or out of the system (discharge) at a rate that is 
    proportional to the amount of water held in the reservoir.
    """

    import numpy as np

    def __init__(self, t_efold, f_to_discharge=1., Vmax=np.inf):
        """
        t_efold: e-folding time for reservoir depletion
        f_to_discharge: fraction of the water lost during that e-folding time
                        that exfiltrates to river discharge (as opposed to 
                        entering one or more other reservoirs)
        Vmax: Maximum water volume that can be held
        """
        self.Vwater = 0.
        self.Vmax = Vmax
        self.t_efold = t_efold
        self.excess = 0.
        self.Vout = np.nan
        self.f_to_discharge = f_to_discharge
    
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
        self.V_exfiltrated = self.excess + dV * self.f_to_discharge
        self.V_infiltrated = dV * (1 - self.f_to_discharge)
        self.Vwater -= dV
        self.excess = 0.

class buckets(object):
    """
    Incorportates a list of reservoirs into a linear hierarchy that sends water 
    either downwards our out to the surface.
    
    Input: reservoir list, in order from top to bottom
    (surface to deep groundwater)
    """

    def __init__(self, reservoir_list):
        self.reservoirs = reservoir_list
        self.rain = None
        self.Q = [] # discharge
    
    def set_rainfall_time_series(self, rain):
        self.rain = rain
    
    def initialize(self):
        """
        Part of CSDMS BMI
        Initialization handled in __init__
        Nothing more to do
        """
        pass
    
    def update(self, rain_at_timestep):
        """
        Updates water flow for one time step (typically a day)
        
        NOTE FALLACY: recharging before discharging,
        even though during the same day
        consider changing to use half-recharge from each day
        """
        # Top layer is special: interacts with atmosphere
        self.reservoirs[0].recharge(self.rain[ti])
        self.reservoirs[0].discharge()
        Qi += self.reservoirs[0].V_exfiltrated
        for i in range(1, len(self.reservoirs)):
            self.reservoirs[i].recharge(res[i-1].V_infiltrated)
            self.reservoirs[i].discharge()
        return Qi
    
    def run(self, rain=None):
        if rain:
            if self.rain:
                print "Warning: overwriting existing rainfall time series"
            self.set_rainfall_time_series(rain)
        if self.rain is None:
            sys.exit("Please set the rainfall time series")
        for ti in range(len(self.rain)):
            Qi = self.update()
            self.Q.append(Qi)


# Program below the class
import numpy as np
from numpy.random import poisson
from matplotlib import pyplot as plt
#plt.ion()

rain = poisson(.2, 100) # Convert to an import for real data

dt = 1. # day
        

# Change these parameters with an input file or script, eventually
# Arbitrary units; will be made real in an actual landscape
res_surface = reservoir(t_efold=1., f_to_discharge=0.5, Vmax=10.)
res_deep = reservoir(t_efold=10., f_to_discharge=1., Vmax=np.inf)

strat_column = [res_surface, res_deep]

watershed = buckets(strat_column)

Q = []
# This includes made-up rules about how much of the discharge from each
# bucket goes to the other buckets or to the discharge at the outlet.
# I also use arbitarry units
for ti in range(len(rain)):
    Qi = 0.
    # Tile 
    res_tile.recharge(rain[ti])
    res_tile.discharge(dt)
    Qi += res_tile.Vout
    # Surface
    res_surface.recharge( rain[ti] )
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

