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

        self.excess_ET_error = 0.

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
        Recharge can be positive (precipitation > evapotranspiration) or
        negative (evapotranspiration > precipitation)
        """
        self.H_excess = 0. # Extra water above a maximum cap
        self.H_deficit = 0. # Water that this layer cannot hold and
                            # cannot be passed on to a deeper layer

        # ERROR if water is less than 0 -- may be able to remove
        # this check later
        if self.Hwater < 0:
            raise ValueError("Hwater in reservoir < 0; non-physical")

        # What if more water is lost during "recharge" than exists in reservoir?
        # Create a deficit and bring Hwater to 0
        if self.Hwater+H < 0:
            self.H_deficit += self.Hwater+H # inefficient
            self.Hwater = 0
        # What if more water is added than maximum reservoir capacity?
        # Mark excess (straight to runoff) and bring Hwater to Hmax
        elif self.Hwater+H > self.Hmax:
            self.H_excess += self.Hwater+H - self.Hmax
            self.Hwater = self.Hmax
        # Otherwise, we're in a range in which 0 <= H <= Hmax
        # Yay! Things are easier!
        else:
            self.Hwater += H

    # Split exfiltration and discharge re: self.H_excess?
    def discharge(self, dt):
        dH = self.Hwater * (1 - np.exp(-dt/self.t_efold))
        self.H_exfiltrated = dH * self.f_to_discharge
        self.H_discharge = self.H_excess + self.H_exfiltrated
        self.H_infiltrated = dH * (1 - self.f_to_discharge)
        self.Hwater -= dH

class Snowpack(object):

    def __init__(self, melt_factor=None):
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

        self.H_deficit = 0. # Water deficit with neg ET; just this time step
        # If positive recharge
        if H >= 0:
            if self.T <= 0:
                self.Hwater += H
                self.H_infiltrated = 0.
            else:
                # Incoming precip component; melt sums with this
                # This is then directly passed to the first layer of the
                # set of hydrological reservoirs
                self.H_infiltrated = H
        # If negative recharge, take from snowpack, and then from reservoirs
        # in order from top to bottom. If no water is in the bottom reservoir,
        # this deficit will be handled first, but will not produce discharge

        else:
            # Submimation (effectively) if snow present;
            # Otherwise pass water deficit
            if self.Hwater > -H:
                self.Hwater += H
            else:
                self.H_deficit += H + self.Hwater
                self.Hwater = 0
            self.H_infiltrated = 0.

    def discharge(self, dt):
        """
        Currently, sends all snowmelt to the top layer for discharge.
        This is arbitrary and neglects melt atop frozen soil.
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

    def __init__(self):
        # Evapotranspiration constants
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

    def initialize(self, config_file=None):
        """
        Part of CSDMS BMI
        """
        if config_file is None:
            warnings.warn("No configuration file provided; all values needed "+
                          "for a model run therefore must be set indpendently.")

        # Import yml configuration file
        with open(config_file, "r") as ymlfile:
            self.cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)

        # Import dataframe from yml
        self.hydrodata = pd.read_csv(self.cfg['timeseries']['datafile'],
                            parse_dates=['Date'])

        # Set variables on reservoirs
        # First, check if all reservoirs have the same length
        for _key in self.cfg['reservoirs'].keys():
            if len(self.cfg['reservoirs'][_key]) == \
            len(self.cfg['initial_conditions']['water_reservoir_effective_depths__m']):
                pass
            else:
                raise ValueError(_key + ' within "reservoirs" contains a\n'+
                                 'different number of entries, implying'+
                                 'a different number of subsurface water\n'+
                                 'reservoirs, than '+
                                 'water_reservoir_effective_depths__m'+
                                 ' within "initial_conditions".')

        # If all are the same length, then we will assign a number of reservoirs
        self.n_reservoirs = len(self.cfg['initial_conditions']['water_reservoir_effective_depths__m'])
        # Using this, we will build a list of reservoir objects
        # and initialize them based on the provided inputs
        self.reservoirs = [
            Reservoir(
                t_efold = self.cfg['reservoirs']['e_folding_residence_times__days'][i],
                f_to_discharge = self.cfg['reservoirs']['exfiltration_fractions'][i],
                Hmax = self.cfg['reservoirs']['maximum_effective_depths__m'][i],
                H0 = self.cfg['initial_conditions']['water_reservoir_effective_depths__m'][i]
            )
            for i in range(self.n_reservoirs)]
        # Python note: This is a compact way of a loop and using append()

        # Check if bottom reservoir discharges all to river: conserve mass.
        # But allow through with a warning in case the user wants a
        # deep and non-discharging reservoir (although this could be set up)
        # explicitly too)
        if self.reservoirs[-1].f_to_discharge < 1:
            warnings.warn("f_to_discharge of bottom water-storage layer < 1.\n"+
                          "You are not conserving mass.")

        # Set scalar variables based on yaml
        self.melt_factor = self.cfg['snowmelt']['PDD_melt_factor']
        self.et_method = self.cfg['catchment']['evapotranspiration_method']
        self.water_year_start_month = self.cfg['catchment']['water_year_start_month']
        self.drainage_basin_area__km2 = self.cfg['catchment']['drainage_basin_area__km2']

        # Check if there is a mean temperature column for snowpack.
        # If not, note that no snowpack processes will be included
        if 'Mean Temperature [C]' in self.hydrodata.columns:
            # Instantiate snowpack
            self.snowpack = Snowpack(self.melt_factor) # allow changes to melt factor later
        else:
            warnings.warn('"Mean Temperature [C]" has not been set.'+
                            'No snowpack processes will be simulated.')

        # Initial conditions if resuming from prior run
        self.snowpack.Hwater = self.cfg['initial_conditions']['snowpack__m_SWE']
        # H0 in loop above.

        # Check that dt is 1 day everywhere.
        # Do not work otherwise.
        if (self.hydrodata['Date'].diff()[1:] == pd.Timedelta('1 day') ).all():
            self.dt = 1.
        else:
            raise ValueError("All time steps must be 1 day.")

        # Compute specific discharge from data
        self.hydrodata['Specific Discharge [mm/day]'] = \
                self.hydrodata['Discharge [m^3/s]'] \
                / (self.drainage_basin_area__km2*1E3) * 86400

        # Create columns for model output
        # But not needed here -- don't have for subsurface storage !!!!!
        self.hydrodata['Specific Discharge (modeled) [mm/day]'] = pd.NA
        self.hydrodata['Snowpack (modeled) [mm SWE]'] = pd.NA

        # Start out at first timestep
        # Could modify this to pick up a run in the middle
        # Or start at the beginning of a water year
        # for example
        self._timestep_i = self.hydrodata.index[0]

        # If there is a water-supply deficit (ET > P) that is larger than
        # any existing water reservoirs, note this in this variable.
        # Important: This is the cumulative H_deficit, whereas
        # class Snowpack's H_deficit is for that time step only.
        self.H_deficit = 0.

    def compute_water_year(self):
        """
        Adds a "water year" column to the Pandas DataFrame
        """
        self.hydrodata['Water Year'] = pd.DatetimeIndex(self.hydrodata['Date']).year
        self.hydrodata['Water Year'] += \
            pd.DatetimeIndex(self.hydrodata['Date']).month >= self.water_year_start_month

    def compute_ET_multiplier(self):
        """
        Calculates the total amount of evapotranspiration required to balance
        discharge and precipitation in each water year
        """
        # Originally used "sum", but then used "mean" so the headers would
        # still be sensible
        self.hydromeansWY = self.hydrodata.groupby(self.hydrodata['Water Year']).mean()
        # Not needed, but no real harm in calculating
        self.hydromeansWY['Runoff ratio'] = \
                            self.hydromeansWY['Specific Discharge [mm/day]'] / \
                            self.hydromeansWY['Precipitation [mm/day]']
        _ET_required = - ( self.hydromeansWY['Specific Discharge [mm/day]'] -
                            self.hydromeansWY['Precipitation [mm/day]'] )
        self.hydromeansWY['ET multiplier'] = _ET_required / \
                            self.hydromeansWY['Evapotranspiration [mm/day]']

    def compute_ET(self):
        """
        Computes an evapotranspiration column from the input data.

        Does this in two steps:
        1. Initial ET (provided or from Thorntwaite)
        2. Modifying this over a water year to enforce water balance
        """
        # Nonfunctional update in progress
        if self.et_method == 'datafile':
            _raw_ET = self.hydrodata['Evapotranspiration [mm/day]']
        elif self.et_method == 'ThorntwaiteChang2019':
            _raw_ET = self.evapotranspirationChang2019()
        else:
            raise ValueError('evapotranspiration_method must be "datafile" or '+
                             '"ThorntwaiteChang2019".')
        # There should be a better way to do this fully in an operation
        # rather than adding it to the dataframe + memory
        # But this is pretty straightforward and doesn't use much memory
        self.hydrodata = self.hydrodata.merge(
                                    self.hydromeansWY['ET multiplier'],
                                    on='Water Year' )
        self.hydrodata['ET for model [mm/day]'] = \
                                    _raw_ET * self.hydrodata['ET multiplier']

    def update(self, time_step = None):
        """
        Updates water flow for one time step (typically a day)

        NOTE FALLACY: recharging before discharging,
        even though during the same time step
        consider changing to use half-recharge from each time step

        FOR LATER: , dt_at_timestep=self.dt
        FOR SOONER: WATER-YEAR BALANCE
        """

        # time_step sets the index of the row in the Pandas DataFrame
        # that will be used to update, run, and store the outputs
        if time_step is None:
            time_step = self._timestep_i
            # Advance internal variable if external time step is not selected
            # This should be a different variable and therefore not
            # modify the value of "time_step" by reference
            self._timestep_i += 1

        # If no mean temperature included, no snowpack processes simulated
        if 'Mean Temperature [C]' in self.hydrodata.columns:
            self.snowpack.set_temperature(
                    self.hydrodata['Mean Temperature [C]'][time_step] )
            # Recharge P - ET + running deficit (neg if deficit exists)
            self.snowpack.recharge(
                    self.hydrodata['Precipitation [mm/day]'][time_step] -
                    self.hydrodata['ET for model [mm/day]'][time_step]
                    + self.H_deficit)
            self.snowpack.discharge(self.dt)
            # Specific discharge from snowpack
            qi = self.snowpack.H_discharge
            # An existing deficit will cause H_deficit to be negative
            # H_deficit is updated based on snowpack processes
            self.H_deficit = self.snowpack.H_deficit
        else:
            # Just declare variable at 0 if no snowpack processes
            qi = 0.
        # First, compute snowpack (if being used) and direct discharge
        for i in range(0, len(self.reservoirs)):
            # Top layer is special: snowmelt and/or precip infiltrates
            # immediately (at least on our daily time scales)
            if i == 0:
                if 'Mean Temperature [C]' in self.hydrodata.columns:
                    self.reservoirs[i].recharge(self.snowpack.H_infiltrated
                                                + self.H_deficit)
                else:
                    self.reservoirs[i].recharge(
                        self.hydrodata['Precipitation [mm/day]'][time_step] -
                        self.hydrodata['ET for model [mm/day]'][time_step] +
                        self.H_deficit
                    )
            else:
                # Let water infiltrate to lower layers effectively
                # instantaneously; this isn't quite realistic, but
                # should be a simpler approach for parameter calibration
                # (Plus, this is just the water that did exit that above
                # container, which is already free to discharge, so this
                # seems more self-consistent.)
                # The amount of infiltrated water from above could be
                # negative; this represents ET in excess of what the
                # unsaturated zone ("soil zone"; top reservoir) holds.
                # Deeper loss of water could be due to plants tapping into
                # groundwater, direct lake evaporation, etc. -- or related
                # to this model not being physical or distributed, so just
                # needing to balance mass.
                self.reservoirs[i].recharge(self.reservoirs[i-1].H_infiltrated
                                    + self.reservoirs[i-1].H_deficit)
            self.reservoirs[i].discharge(self.dt)
            qi += self.reservoirs[i].H_discharge
        # At the bottom of the heap, if a deficit still exists, then
        # note this the running water budget
        # (Just add: Is 0 unless )
        # And then we have to deal wtih this in the next round
        # If 0, great. If not 0, passes on
        self.H_deficit = self.reservoirs[i].H_deficit

        # No need to return value anymore; just place it in the data table directly
        # return Qi
        self.hydrodata.at[time_step, 'Specific Discharge (modeled) [mm/day]'] = qi
        self.hydrodata.at[time_step, 'Snowpack (modeled) [mm SWE]'] = self.snowpack.Hwater
        self.hydrodata.at[time_step, 'Subsurface storage (modeled total) [mm]'] = \
                                np.sum([res.Hwater for res in self.reservoirs])

    def evapotranspirationChang2019(self, Tmax = None, Tmin = None,
                                                    photoperiod = None):
        """
        Modified daily Thorntwaite Equation
        """
        if Tmax is None:
            Tmax = self.hydrodata['Maximum Temperature [C]']
        if Tmin is None:
            Tmin = self.hydrodata['Maximum Temperature [C]']
        if photoperiod is None:
            photoperiod = self.hydrodata['Photoperiod [hr]']

        Teff = 0.5 * 0.69 * (3*Tmax - Tmin)
        C = photoperiod/360.

        # Simple and inefficient logical implementation
        # Is this ET or PET?
        # I will add/subtract from the top reservoir only.
        return C*(-415.85 + 32.24*Teff - 0.43*Teff**2) * (Teff >= 26) \
                   + 16.*C * (10.*Teff / self.Chang_I)**self.Chang_a_i \
                     * (Teff > 0) * (Teff < 26)

    def run(self):
        for ti in self.hydrodata.index:
            self.update()

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
        plt.bar(self.hydrodata['Date'], height=self.hydrodata['Precipitation [mm/day]']/self.dt, width=1., align='center', label='Precipitation [mm/day]', linewidth=0, alpha=0.5)
        #plt.plot(self.time, self.rain/self.dt, 'b-', label='Rainfall', alpha=0.5)
        plt.twinx()
        plt.plot(self.hydrodata['Date'], self.hydrodata['Specific Discharge [mm/day]'],
                            'b', label='Specific discharge (data)', linewidth=2)
        plt.plot(self.hydrodata['Date'], self.hydrodata['Specific Discharge (modeled) [mm/day]'],
                'k', label='Specific discharge (model)', linewidth=2)
        plt.ylim(0, plt.ylim()[-1])
        plt.legend(fontsize=11)
        plt.ylabel('[mm/day]', fontsize=14)
        plt.tight_layout()
        plt.show()

    def check_mass_balance(self, time_step=None):
        """
        Add up balance components in mm
        """
        # Additions equals discharge out; set up this way, and can check.
        total_additions = \
            self.hydrodata['Precipitation [mm/day]'][:time_step+1].sum() \
            - self.hydrodata['ET for model [mm/day]'][:time_step+1].sum()
        # Storage reservoirs
        snow_storage = self.hydrodata['Snowpack (modeled) [mm SWE]'][time_step]
        subsurface_storage = self.hydrodata['Subsurface storage (modeled total) [mm]'][time_step]
        # Mass removal
        outlet_discharge = self.hydrodata['Specific Discharge (modeled) [mm/day]'][:time_step+1].sum()
        # Carried-over water deficit
        deficit = self.H_deficit

        # Discrepancy
        excess_mass_in_model = outlet_discharge + subsurface_storage \
                                + snow_storage - total_additions + deficit

        return excess_mass_in_model

    def computeNSE(self, verbose=False):
        """
        Compute the Nash-Sutcliffe Efficiency of measured vs. modeled
        specific discharge
        """

        # Shorthand for fcn
        q_data = self.hydrodata['Specific Discharge [mm/day]']
        q_model = self.hydrodata['Specific Discharge (modeled) [mm/day]']

        # Calculate NSE
        _realvalue = ~q_model.isna() * ~q_data.isna()
        NSE_num = np.sum( (q_model[_realvalue] - q_data[_realvalue])**2 )
        NSE_denom = np.sum( (q_data[_realvalue] -
                                    np.mean(q_data[_realvalue]))**2 )
        if np.sum(1 - _realvalue):
            print("Calculated with ", np.sum(1 - _realvalue), "no-data points")

        self.NSE = 1 - NSE_num / NSE_denom

        if verbose:
            print( "NSE:", self.NSE )
