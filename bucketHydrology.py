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

    def __init__(self, t_efold, f_to_discharge=1., Hmax=np.inf):
        """
        t_efold: e-folding time for reservoir depletion
        f_to_discharge: fraction of the water lost during that e-folding time
                        that exfiltrates to river discharge (as opposed to 
                        entering one or more other reservoirs)
        Hmax: Maximum water volume that can be held
        """
        self.Hwater = 0.
        self.Hmax = Hmax
        self.t_efold = t_efold
        self.excess = 0.
        self.Hout = np.nan
        self.f_to_discharge = f_to_discharge
    
    def recharge(self, H):
        if self.Hwater+H <= self.Hmax:
            excess = 0.
            self.Hwater += H
        if self.Hwater+H > self.Hmax:
            excess = self.Hwater+H - self.Hmax
            self.Hwater = self.Hmax
        self.excess += excess

    def discharge(self, dt):
        dH = self.Hwater * (1 - np.exp(-dt/self.t_efold))
        self.H_exfiltrated = self.excess + dH * self.f_to_discharge
        self.H_infiltrated = dH * (1 - self.f_to_discharge)
        self.Hwater -= dH
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
        Qi = self.reservoirs[0].H_exfiltrated
        for i in range(1, len(self.reservoirs)):
            self.reservoirs[i].recharge(self.reservoirs[i-1].H_infiltrated)
            self.reservoirs[i].discharge(self.dt)
            Qi += self.reservoirs[i].H_exfiltrated
        return Qi
    
    def run(self, rain=None):
        if rain is not None:
            if self.rain is not None:
                print "Warning: overwriting existing rainfall time series"
            self.set_rainfall_time_series(rain)
        if self.rain is None:
            sys.exit("Please set the rainfall time series")
        self.time = np.arange(len(self.rain)) * self.dt
        for rain_ti in self.rain:
            Qi = self.update(rain_ti)
            self.Q.append(Qi)
        self.rain = np.array(self.rain)
        self.Q = np.array(self.Q)

    def plot(self):
        plt.figure()
        plt.bar(left=self.time, height=self.rain/self.dt, width=1.,
                align='center', label='Rainfall', linewidth=0, alpha=0.5)
        plt.plot(self.time, self.Q/self.dt, 'k',
                label='Unit discharge', linewidth=2)
        plt.legend(fontsize=11)
        plt.ylabel('[mm/day]', fontsize=14)
        plt.xlabel('Time [days]', fontsize=14)
        plt.show()

    def computeNashSutcliffeEfficiency(self, Qdata):
        """
        Compute the NSE of the model outputs vs. a set of supplied data
        """
        _realvalue = np.isfinite(self.Q * Qdata)
        NSE_num = np.sum( (self.Q[_realvalue] - Qdata[_realvalue])**2 )
        NSE_denom = np.sum((Qdata[_realvalue] - np.mean(Qdata[_realvalue]))**2)
        if np.sum(1 - _realvalue):
            print "Calculated with ", np.sum(1 - _realvalue), "no-data points"
        self.NSE = NSE_num / NSE_denom
        
