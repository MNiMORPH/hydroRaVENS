# Started by A. Wickert
# 25 July 2019
# Updated by J. Jones
# Starting 08 Oct 2019

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import matplotlib.dates as mdates
import sys
import warnings
import yaml

class Reservoir(object):
    """
    Generic reservoir. Accepts new water (recharge), and sends it to other
    reservoirs and/or out of the system (discharge) at a rate that is
    proportional to the amount of water held in the reservoir.
    """

    import numpy as np

    def __init__(self, t_efold, f_to_discharge=1., Hmax=np.inf, H0=0.):
        """
        t_efold: e-folding time for reservoir depletion (same units as time
                 steps; typically days)
        f_to_discharge: fraction of the water lost during that e-folding time
                        that exfiltrates to river discharge (as opposed to
                        entering one or more other reservoirs)
        Hmax: Maximum water volume that can be held
        """
        self.Hwater = H0
        self.Hmax = Hmax
        self.t_efold = t_efold
        self.excess = 0.
        self.Hout = np.nan
        self.f_to_discharge = f_to_discharge

        # Check values and note whether they are reasonable
        if t_efold < 0:
            raise ValueError("Negative t_efold nonsensical.")
        if f_to_discharge < 0:
            raise ValueError("Negative f_to_discharge not possible")
        elif f_to_discharge > 1:
            raise ValueError("f_to_discharge: Cannot discharge >100% of water")
        elif f_to_discharge == 0:
            warnings.warn("All water infiltrates when f_to_discharge is 0:"+
                          " you may have created a\n"+
                          "redudnant pass-through water-storage layer")
        if Hmax < 0:
            raise ValueError("Hmax must be >= 0 (and >0 makes more sense)")

    def recharge(self, H):
        """
        Recharge can be positive (precipitation) or negative
        (evapotranspiration)
        """
        excess = 0.
        if self.Hwater < 0:
            # Allowing ET in top layer only.
            self.Hwater = 0
        elif self.Hwater+H <= self.Hmax:
            self.Hwater += H
        elif self.Hwater+H > self.Hmax:
            excess = self.Hwater+H - self.Hmax
            self.Hwater = self.Hmax
        self.excess += excess

    def discharge(self, dt):
        dH = self.Hwater * (1 - np.exp(-dt/self.t_efold))
        self.H_exfiltrated = dH * self.f_to_discharge
        self.H_discharge = self.excess + self.H_exfiltrated
        self.H_infiltrated = dH * (1 - self.f_to_discharge)
        self.Hwater -= dH
        self.excess = 0.

class Snowpack(object):

    def __init__(self, melt_factor=0.5):
        """
        A unique reservoir type that adds and removes water based on
        temperature.

        If included in a list of reservoirs within a watershed model,
        should always be on top.

        The melt factor is given as a positive-degree-day factor [mm/day melt]
        """
        self.Hwater = 0. # SWE
        self.melt_factor = melt_factor

    def set_melt_factor(self, melt_factor):
        """
        Change the melt factor (PDD)
        """
        self.melt_factor = melt_factor

    def set_temperature(self, T):
        self.T = T

    def recharge(self, H):
        """
        If T <= 0, all recharge to snowpack
        Else, recharge magically bypasses snowpack
        """
        if self.T <= 0:
            self.Hwater += H
            self.H_infiltrated = 0.
        else:
            # Incoming precip component; melt sums with this
            self.H_infiltrated = H

    def discharge(self, dt):
        """
        Currently, 50/50 snowmelt split between infiltration and discharge.
        This is arbitrary.
        """
        if self.T > 0:
            dH_melt = np.min((self.Hwater, self.melt_factor * self.T * dt))
            self.H_infiltrated += dH_melt * 1.
        else:
            dH_melt = 0
        self.H_discharge = dH_melt * 0.
        self.Hwater -= dH_melt

class Buckets(object):
    """
    Incorportates a list of reservoirs into a linear hierarchy that sends water
    either downwards our out to the surface.

    reservoir_list: list of subsurface layers in order from top to bottom
                    (surface to deep groundwater)

    """

    def __init__(self, reservoir_list):
        self.snowpack = snowpack() # allow changes to melt factor later
        self.reservoirs = reservoir_list
        self.rain = None
        self.ET = None

        # Check if bottom reservoir discharges all to river:
        # conserve mass
        # But allow through with a warning in case the user wants a
        # deep and non-discharging reservoir (although this could be set up)
        # explicitly too)
        if self.reservoirs[-1].f_to_discharge < 1:
            warnings.warn("f_to_discharge of bottom water-storage layer < 1.\n"+
                          "You are not conserving mass.")

        # Evapotranspiration
        self.Chang_I = 41.47044637
        self.Chang_a_i = 6.75E-7*self.Chang_I**3 \
                         - 7.72E-5*self.Chang_I**2 \
                         + 1.7912E-2*self.Chang_I \
                         + 0.49239

    def set_rainfall_time_series(self, rain):
        self.rain = np.array(rain)

    def set_mean_temperature(self, T):
        """
        Mean temperature each time step for the temperature-index approch to
        basic snowpack modeling
        """
        self.Tmean = np.array(T)

    def export_Hlist(self):
        """
        Export the list of water depths, for reinitialization
        (e.g., to start after a spin-up phase or to restart a paused model run)
        """
        Hlist = []
        for reservoir in self.reservoirs:
            Hlist.append( reservoir.Hwater )
        return Hlist

    def __init__(self):
        # Empty __init__ as another option
        # Perhaps overloading isn't the best way to go here.
        # But I wil hold it here as a bookmark until I change (potentially)
        # the full initialization and instantiation method set
        pass

    def initialize(self, config_file=None):
        if config_file is None:
            warnings.warn("No configuration file provided; all values needed "+
                          "for a model run therefore must be set indpendently.")

        # Import yml configuration file
        with open(config_file, "r") as ymlfile:
            cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)

        # Import dataframe from yml
        self.hydrodata = pd.read_csv(cfg['timeseries']['datafile'])

        # Set variables on reservoirs
        # First, check if all reservoirs have the same length
        for _key in cfg['reservoirs'].keys():
            if len(cfg['reservoirs'][_key]) == \
            len(cfg['initial_conditions']['water_reservoir_effective_depths__m']):
                pass
            else:
                raise ValueError(_key + 'within "reservoirs" contains a\n'+
                                 'different number of entries, implying'+
                                 'a different number of subsurface water\n'+
                                 'reservoirs, than '+
                                 'water_reservoir_effective_depths__m'+
                                 'within "initial_conditions".')

        # If all are the same length, then we will assign a number of reservoirs
        self.n_reservoirs = len(cfg['initial_conditions']['water_reservoir_effective_depths__m'])
        # Using this, we will build a list of reservoir objects
        # and initialize them based on the provided inputs
        self._Reservoir_object_list = [
                Reservoir(
                t_efold = cfg['reservoirs']['e_folding_residence_times__days'][i],
                f_to_discharge = cfg['reservoirs']['exfiltration_fractions'][i],
                Hmax = cfg['reservoirs']['maximum_effective_depths__m'][i],
                H0 = cfg['initial_conditions']['water_reservoir_effective_depths__m'][i]
                )
                for i in range(self.n_reservoirs)]
        # Python note: This is a compact way of a loop and using append()

        # Set variables based on yaml
        """
        self.Hwater = 0.
        self.Hmax = Hmax
        self.t_efold = t_efold
        self.excess = 0.
        self.Hout = np.nan
        self.f_to_discharge = f_to_discharge
        """

    def old_initialize(self, dateTimes, Hlist=None, SWE=None):
        """
        Part of CSDMS BMI
        Can use this to initialize from an old run or a spin-up

        dateTimes: A list/array/etc. of python DateTime objects that go along
                   with the time series of precipitation.
                   IF A FUTURE CAPABILITY FOR NONUNIFORM DATA IS ENACTED
                   The calculation will be exclusive of the outer two data
                   points in this series in order to calculate a time step,
                   *unless* a specific dt is indicated
                   FOR FUTURE: dt: Time step [days]; , dt=None
        Hlist: A list of water depths, in the same order as the reservoirs
               (top to bottom) that can be used to set initial water depths.
               This can be helpful for initial conditions or to return from a
               spin up.
        """
        self.time = np.array(dateTimes)
        self.dt = self.__compute_dt()
        if Hlist is not None:
            i = 0
            for reservoir in self.reservoirs:
                reservoir.Hwater = Hlist[i]
                i += 1
        if SWE is not None:
            self.snowpack.Hwater = SWE

    def __compute_dt(self, scalar_dt=True):
        """
        Calculates the time step from the DateTime series of the input data
        If scalar_dt:
            Tests for and returns a single dt in days
        Else:
            Returns it in days using a centered approach for nonuniformly
            spaced data. THIS IS NOT YET IMPLEMENTED IN THE FULL CODE
        """
        _dt = np.diff(self.time)
        if scalar_dt:
            _dt = float( np.mean(np.diff(self.time)) )
            if (np.diff(self.time).astype(float) == _dt).all():
                _dt /= 86400E9 # convert to days
            else:
                sys.exit("Holes in input data series")
        else:
            dt_timeSeries = np.array( _dt[:-1] + _dt[1:] ).astype(float) / 2.
            dt = dt_timeSeries / 86400E9 # convert to days
        return _dt

    def update(self, rain_at_timestep, ET_at_timestep, T_at_timestep):
        """
        Updates water flow for one time step (typically a day)

        rain_at_timestep: rainfall rate
        ET_at_timestep: evapotranspiration
        T_at_timestep: temperature; default "10" as a number above zero,
                       meaning that there is no snowpack if T is unset

        NOTE FALLACY: recharging before discharging,
        even though during the same time step
        consider changing to use half-recharge from each time step

        FOR LATER: , dt_at_timestep=self.dt
        """
        recharge_at_timestep = rain_at_timestep - ET_at_timestep
        self.snowpack.set_temperature(T_at_timestep)
        self.snowpack.recharge(recharge_at_timestep)
        self.snowpack.discharge(self.dt)
        Qi = self.snowpack.H_discharge
        # First, compute snowpack and direct discharge
        for i in range(0, len(self.reservoirs)):
            # Top layer is special: snowpack
            if i == 0:
                self.reservoirs[i].recharge(self.snowpack.H_infiltrated)
            else:
                self.reservoirs[i].recharge(self.reservoirs[i-1].H_infiltrated)
            self.reservoirs[i].discharge(self.dt)
            Qi += self.reservoirs[i].H_discharge
        # Then passs infiltrated water downwards
        # Using this as a separate step so the water can take only one step
        # per time step -- either down or out. This is quite schematic.
        for i in range(1, len(self.reservoirs)):
            self.reservoirs[i].Hwater += self.reservoirs[i-1].H_infiltrated
        return Qi

    def evapotranspirationChang2019(self, Tmax, Tmin, photoperiod):
        """
        Modified daily Thorntwaite Equation
        """

        Teff = 0.5 * 0.69 * (3*Tmax - Tmin)
        C = photoperiod/360.

        # Simple and inefficient logical implementation
        # Is this ET or PET?
        # I will add/subtract from the top reservoir only.
        self.ET = C*(-415.85 + 32.24*Teff - 0.43*Teff**2) * (Teff >= 26) \
                   + 16.*C * (10.*Teff / self.Chang_I)**self.Chang_a_i \
                     * (Teff > 0) * (Teff < 26)

    def set_evapotranspiration(self, ET):
        """
        User provides ET time series
        """
        self.ET = np.array(ET)

    def run(self, rain=None, ET_flag=False, Tmean_flag=False):
        # SET UP VARIABLES
        self.Q = [] # discharge
        self.SWE = [] # snowpack
        if rain is not None:
            if self.rain is not None:
                print("Warning: overwriting existing rainfall time series.")
            self.set_rainfall_time_series(rain)
        if self.rain is None:
            sys.exit("Please set the rainfall time series")
        if ET_flag:
            ET = self.ET
        else:
            print("Warning: neglecting evapotranspiration.")
            ET = np.zeros(self.rain.shape)
        if Tmean_flag:
            Tmean = self.Tmean
        else:
            print("Warning: neglecting snowpack")
            Tmean = np.ones(self.rain.shape) # >0 for no snowpack
        # RUN
        for ti in range(len(self.rain)):
            Qi = self.update(rain[ti], ET[ti], Tmean[ti])
            self.Q.append(Qi)
            self.SWE.append(self.snowpack.Hwater)
        self.Q = np.array(self.Q)

    def plot(self, Qdata=None):
        """
        Plot rainfall and discharge.
        Optionally pass specific discharge data to plot this as well.
        """
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%y-%m-%d'))
        plt.xlabel('Date', fontsize=14)
        #for tick in ax.get_xticklabels():
        #    tick.set_rotation(45)
        plt.xticks(rotation=45, horizontalalignment='right')
        plt.ylabel('[mm/day]', fontsize=14)
        plt.bar(self.time, height=self.rain/self.dt, width=1., align='center', label='Rainfall', linewidth=0, alpha=0.5)
        #plt.plot(self.time, self.rain/self.dt, 'b-', label='Rainfall', alpha=0.5)
        plt.twinx()
        if Qdata is not None:
            plt.plot(self.time, Qdata/self.dt, 'b',
                    label='Unit discharge data', linewidth=2)
        plt.plot(self.time, self.Q/self.dt, 'k',
                label='Unit discharge model', linewidth=2)
        plt.ylim(0, plt.ylim()[-1])
        plt.legend(fontsize=11)
        plt.ylabel('[mm/day]', fontsize=14)
        plt.tight_layout()
        plt.show()

    def computeNashSutcliffeEfficiency(self, Qdata):
        """
        Compute the NSE of the model outputs vs. a set of supplied data
        """
        _realvalue = np.isfinite(self.Q * Qdata)
        NSE_num = np.sum( (self.Q[_realvalue] - Qdata[_realvalue])**2 )
        NSE_denom = np.sum((Qdata[_realvalue] - np.mean(Qdata[_realvalue]))**2)
        if np.sum(1 - _realvalue):
            print("Calculated with ", np.sum(1 - _realvalue), "no-data points")
        self.NSE = 1 - NSE_num / NSE_denom
