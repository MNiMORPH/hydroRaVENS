#! /usr/bin/env python

# Started by A. Wickert
# 25 July 2019
# Updated by J. Jones
# Starting 08 Oct 2019

import numpy as np
from matplotlib import pyplot as plt
import sys

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
    
    reservoir_list: list of subsurface layers in order from top to bottom
                    (surface to deep groundwater)
    dt: time step (typically in days)
    
    """

    def __init__(self, reservoir_list, dt=1):
        self.reservoirs = reservoir_list
        self.dt = dt
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
        even though during the same time step
        consider changing to use half-recharge from each time step
        """
        # Top layer is special: interacts with atmosphere
        self.reservoirs[0].recharge(rain_at_timestep)
        self.reservoirs[0].discharge(self.dt)
        Qi = self.reservoirs[0].V_exfiltrated
        for i in range(1, len(self.reservoirs)):
            self.reservoirs[i].recharge(self.reservoirs[i-1].V_infiltrated)
            self.reservoirs[i].discharge(self.dt)
            Qi += self.reservoirs[i].V_exfiltrated
        return Qi
    
    def run(self, rain=None):
        if rain is not None:
            if self.rain is not None:
                print "Warning: overwriting existing rainfall time series"
            self.set_rainfall_time_series(rain)
        if self.rain is None:
            sys.exit("Please set the rainfall time series")
        for rain_ti in self.rain:
            Qi = self.update(rain_ti)
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
watershed = buckets(reservoir_list=strat_column, dt=dt)

watershed.run(rain)

plt.figure()
plt.plot(watershed.rain, 'g', label='rainfall')
plt.plot(watershed.Q, 'b', label='discharge')
plt.legend(fontsize=11)
plt.ylabel('Rainfall or discharge [units arbitrary]', fontsize=14)
plt.xlabel('Time [units arbitrary]', fontsize=14)
#plt.savefig('TDF-'+'%.1f' %f_tile+'.png')
plt.show()

