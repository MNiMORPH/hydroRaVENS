#! /usr/bin/python3

########################################
# Then, methought, the air grew denser #
#         - Edgar Allan Poe            #
#              THE RAVEN               #
########################################

# Started by A. Wickert
# 25 July 2019
# Updated slightly by J. Jones
# 08 Oct 2019
# Significant Update by A. Wickert
# October 2022
# CLI added by A. Wickert
# November 2023

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import matplotlib.dates as mdates
import sys
import warnings
import yaml
import argparse


class Reservoir(object):
    """
    Generic reservoir. Accepts new water (recharge), and sends it to other
    reservoirs and/or out of the system (discharge) at a rate that is
    proportional to the amount of water held in the reservoir.
    """

    def __init__(self, t_efold, f_to_discharge=1., Hmax=np.inf, H0=0.):
        """
        Initialize a linear reservoir.

        Parameters
        ----------
        t_efold : float
            E-folding residence time for reservoir depletion (days, or
            whatever time unit matches the model time steps).
        f_to_discharge : float, optional
            Fraction of water lost each time step that exits as river
            discharge. The remainder (1 - f_to_discharge) infiltrates to
            the next-deeper reservoir. Default 1.0 (all to discharge).
        Hmax : float, optional
            Maximum water depth the reservoir can hold. Default np.inf.
        H0 : float, optional
            Initial water depth at the start of the simulation. Default 0.

        Raises
        ------
        ValueError
            If t_efold <= 0, f_to_discharge < 0 or > 1, or Hmax < 0.
        """
        self.Hwater = H0
        self.Hmax = Hmax
        self.t_efold = t_efold
        self.f_to_discharge = f_to_discharge

        # Initialized here so all instance attributes exist before
        # recharge() and discharge() are first called
        self.H_excess = 0.
        self.H_deficit = 0.
        self.H_exfiltrated = 0.
        self.H_infiltrated = 0.
        self.H_discharge = 0.

        # Check values and note whether they are reasonable
        if t_efold <= 0:
            raise ValueError("t_efold must be > 0.")
        if f_to_discharge < 0:
            raise ValueError("Negative f_to_discharge not possible")
        elif f_to_discharge > 1:
            raise ValueError("f_to_discharge: Cannot discharge >100% of water")
        elif f_to_discharge == 0:
            warnings.warn("All water infiltrates when f_to_discharge is 0:"+
                          " you may have created a\n"+
                          "redundant pass-through water-storage layer")
        if Hmax < 0:
            raise ValueError("Hmax must be >= 0 (and >0 makes more sense)")

    def recharge(self, H):
        """
        Add or remove water from the reservoir.

        Recharge H can be positive (net water input, e.g. P > ET) or
        negative (net deficit, e.g. ET > P). Sets H_excess if the
        reservoir overflows Hmax, or H_deficit if more water is removed
        than the reservoir holds.

        Parameters
        ----------
        H : float
            Depth of water added (positive) or removed (negative).

        Raises
        ------
        ValueError
            If Hwater is already negative before recharge is applied.
        """
        # Extra water above a maximum cap
        self.H_excess = 0.
        # Water that this layer cannot hold and cannot be passed to a deeper layer
        self.H_deficit = 0.

        # ERROR if water is less than 0 -- may be able to remove
        # this check later
        if self.Hwater < 0:
            raise ValueError("Hwater in reservoir < 0; non-physical")

        # What if more water is lost during "recharge" than exists in reservoir?
        # Create a deficit and bring Hwater to 0
        if self.Hwater+H < 0:
            self.H_deficit += self.Hwater+H
            self.Hwater = 0.
        # What if more water is added than maximum reservoir capacity?
        # Mark excess (straight to runoff) and bring Hwater to Hmax
        elif self.Hwater+H > self.Hmax:
            self.H_excess += self.Hwater+H - self.Hmax
            self.Hwater = self.Hmax
        # Otherwise, we're in a range in which 0 <= H <= Hmax
        # Yay! Things are easier!
        else:
            self.Hwater += H

    def discharge(self, dt):
        """
        Discharge water from the reservoir over one time step.

        Computes water lost by exponential decay, partitions it between
        river discharge (H_exfiltrated) and infiltration to the next-deeper
        reservoir (H_infiltrated), and adds overflow from recharge()
        (H_excess) to H_discharge.

        Parameters
        ----------
        dt : float
            Time step duration (same units as t_efold; typically days).
        """
        dH = self.Hwater * (1 - np.exp(-dt/self.t_efold))
        self.H_exfiltrated = dH * self.f_to_discharge
        self.H_discharge = self.H_excess + self.H_exfiltrated
        self.H_infiltrated = dH * (1 - self.f_to_discharge)
        self.Hwater -= dH


class Snowpack(object):
    """
    Snowpack reservoir driven by temperature.

    Accumulates precipitation as snow when mean temperature is at or below
    0 °C. Melts at a positive-degree-day rate when temperature is above 0 °C.
    All melt is routed to the top subsurface reservoir as infiltration; direct
    discharge to the river is not modeled.

    Should precede the subsurface reservoir list in a watershed model.
    The melt factor is a positive-degree-day factor [mm/°C/day].
    """

    def __init__(self, melt_factor=None):
        """
        Initialize an empty snowpack.

        Parameters
        ----------
        melt_factor : float, optional
            Positive-degree-day melt factor (mm SWE °C⁻¹ day⁻¹).
            Can be set or updated later via set_melt_factor().
        """
        self.Hwater = 0.  # SWE
        self.melt_factor = melt_factor
        self.T = 0.
        self.H_infiltrated = 0.
        self.H_discharge = 0.
        self.H_deficit = 0.

    def set_melt_factor(self, melt_factor):
        """
        Set or update the positive-degree-day melt factor.

        Parameters
        ----------
        melt_factor : float
            Melt rate per positive degree-day (mm SWE °C⁻¹ day⁻¹).
        """
        self.melt_factor = melt_factor

    def set_temperature(self, T):
        """
        Set the mean air temperature for the current time step.

        Parameters
        ----------
        T : float
            Mean air temperature (°C).
        """
        self.T = T

    def recharge(self, H):
        """
        Apply net water input or deficit to the snowpack.

        If T <= 0, positive H accumulates as snow (SWE). If T > 0,
        positive H bypasses the snowpack and is passed directly to the
        top subsurface reservoir via H_infiltrated. Negative H (ET > P)
        is removed from the snowpack as sublimation; any remainder that
        exceeds available SWE becomes H_deficit.

        Parameters
        ----------
        H : float
            Net water depth for this time step (mm). Positive = input
            (P - ET > 0); negative = deficit (ET - P > 0).
        """

        self.H_deficit = 0.  # Water deficit with neg ET; just this time step
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
            # Sublimation (effectively) if snow present;
            # Otherwise pass water deficit
            if self.Hwater > -H:
                self.Hwater += H
            else:
                self.H_deficit += H + self.Hwater
                self.Hwater = 0
            self.H_infiltrated = 0.

    def discharge(self, dt):
        """
        Compute snowmelt and route it to the top subsurface reservoir.

        Melt is computed using the positive-degree-day method and added
        to H_infiltrated. All melt infiltrates to the top reservoir;
        direct snowmelt runoff to the river is not modeled. H_discharge
        is always 0 for the snowpack.

        Parameters
        ----------
        dt : float
            Time step duration (days).

        Notes
        -----
        Routing all melt to infiltration is a simplification that
        neglects surface runoff from melt atop frozen or saturated soil.
        """
        if self.T > 0:
            dH_melt = np.min((self.Hwater, self.melt_factor * self.T * dt))
        else:
            dH_melt = 0.
        self.H_infiltrated += dH_melt
        self.H_discharge = 0.
        self.Hwater -= dH_melt


class Buckets(object):
    """
    Incorporates a list of reservoirs into a linear hierarchy that sends water
    either downwards or out to the surface. Reservoirs are ordered from top
    (nearest Earth's surface) to bottom (deepest groundwater); this order
    controls the direction of infiltration between layers.
    """

    def __init__(self, T_monthly_normals=None):
        """
        Initialize the watershed model.

        If using the ThorntwaiteChang2019 ET method, pass
        T_monthly_normals here so that the thermal index I and exponent
        a are computed once from climatological normals and remain fixed
        throughout the simulation.

        Parameters
        ----------
        T_monthly_normals : array-like of length 12, optional
            Long-term mean monthly temperatures (°C) used to compute the
            Thornthwaite thermal index I and exponent a per Chang et al.
            (2019), https://doi.org/10.1002/ird.2309. Required when
            evapotranspiration_method is 'ThorntwaiteChang2019'.
        """
        # Thornthwaite thermal index and exponent, per Chang et al. (2019)
        # https://doi.org/10.1002/ird.2309
        # I is climatologically imposed by the local normal temperature regime
        # and must remain fixed during simulation (not recomputed each timestep).
        if T_monthly_normals is not None:
            self.Chang_I = self._compute_Chang_I(T_monthly_normals)
            self.Chang_a = self._compute_Chang_a(self.Chang_I)

        # Frozen ground index (Molnau & Bissell 1983).  Disabled by default
        # (threshold = inf); set fdd_threshold after initialize() to activate.
        self.fdd_threshold = np.inf  # [°C·day]
        self._fgi          = 0.0    # current frozen ground index [°C·day]

    def _compute_Chang_I(self, T_monthly_normals):
        """
        Compute the Thornthwaite thermal index I from long-term monthly normal
        temperatures, per Chang et al. (2019), Eq. 1.
        https://doi.org/10.1002/ird.2309

        Parameters
        ----------
        T_monthly_normals : array-like, length 12
            Long-term mean monthly temperatures (°C). Negative values are
            treated as 0 per the Thornthwaite convention.

        Returns
        -------
        I : float
            Thermal index (dimensionless).
        """
        Tn = np.maximum(T_monthly_normals, 0)
        return np.sum((0.2 * Tn) ** 1.514)

    def _compute_Chang_a(self, I):
        """
        Compute the Thornthwaite exponent a from thermal index I, per
        Chang et al. (2019), Eq. 1.
        https://doi.org/10.1002/ird.2309

        Parameters
        ----------
        I : float
            Thermal index, as returned by _compute_Chang_I.

        Returns
        -------
        a : float
            Thornthwaite exponent (dimensionless).
        """
        return (6.75e-7 * I**3
                - 7.71e-5 * I**2
                + 1.7912e-2 * I
                + 0.49239)

    def export_Hlist(self):
        """
        Return the current water depths in all subsurface reservoirs.

        Useful for checkpointing reservoir state between a spin-up run
        and the main simulation, or for restarting a paused run.

        Returns
        -------
        list of float
            Water depth in each reservoir, ordered from shallowest
            (index 0) to deepest.
        """
        return [reservoir.Hwater for reservoir in self.reservoirs]

    def initialize(self, config_file=None):
        """
        Set up the model from a YAML configuration file.

        Reads the configuration file, loads the input time series, builds
        the reservoir stack, instantiates snowpack if temperature data are
        present, computes the water-year ET multiplier, and runs any
        requested spin-up cycles. Part of the CSDMS Basic Model Interface.

        Parameters
        ----------
        config_file : str, optional
            Path to the YAML configuration file. If None, all required
            values must be set on the object directly before calling
            update().
        """
        if config_file is None:
            warnings.warn("No configuration file provided; all values needed "+
                          "for a model run therefore must be set independently.")

        # Parse YAML configuration file
        # And assign variables except for optimization bounds and plotting
        if config_file is not None:
            try:
                with open(config_file, "r") as yamlfile:
                    self.cfg = yaml.load(yamlfile, Loader=yaml.FullLoader)
            except FileNotFoundError:
                print("\nConfig file not found:", config_file, "\n")
                sys.exit(2)
            except yaml.YAMLError as e:
                print("\nCould not parse config file:", config_file, "\n", e)
                sys.exit(2)

        # Read input time series from the CSV path specified in the config
        self.hydrodata = pd.read_csv(
            self.cfg['timeseries']['datafile'],
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
        self.n_reservoirs = len(
            self.cfg['initial_conditions']['water_reservoir_effective_depths__m'])
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
        # deep and non-discharging reservoir (although this could be set up
        # explicitly too).
        if self.reservoirs[-1].f_to_discharge < 1:
            warnings.warn("f_to_discharge of bottom water-storage layer < 1.\n"+
                          "You are not conserving mass.")

        # Set scalar variables based on yaml
        self.melt_factor = self.cfg['snowmelt']['PDD_melt_factor']
        self.et_method = self.cfg['catchment']['evapotranspiration_method']
        if self.et_method == 'ThorntwaiteChang2019' and not hasattr(self, 'Chang_I'):
            raise ValueError(
                'ThorntwaiteChang2019 requires long-term monthly temperature normals.\n'
                'Pass T_monthly_normals (array of 12 monthly mean temperatures in °C)\n'
                'to Buckets() before calling initialize().'
            )
        self.water_year_start_month = self.cfg['catchment']['water_year_start_month']
        self.drainage_basin_area__km2 = self.cfg['catchment']['drainage_basin_area__km2']

        # Check if there is a mean temperature column for snowpack.
        # If not, note that no snowpack processes will be included
        self.has_snowpack = 'Mean Temperature [C]' in self.hydrodata.columns
        if self.has_snowpack:
            # Instantiate snowpack
            self.snowpack = Snowpack(self.melt_factor)  # allow changes to melt factor later
        else:
            warnings.warn('"Mean Temperature [C]" has not been set. '+
                          'No snowpack processes will be simulated.')

        # How many times to loop the full time series for the spin-up
        # Maybe I should permit a more sophisticated spin-up at some point!
        self.n_spin_up_cycles = self.cfg['general']['spin_up_cycles']

        # Initial conditions if resuming from prior run
        if self.has_snowpack:
            self.snowpack.Hwater = self.cfg['initial_conditions']['snowpack__m_SWE']
        # H0 in loop above.

        # Check that dt is 1 day everywhere.
        # Do not work otherwise.
        if (self.hydrodata['Date'].diff()[1:] == pd.Timedelta('1 day')).all():
            self.dt = 1.
        else:
            raise ValueError("All time steps must be 1 day.")

        # Compute specific discharge from data
        self.hydrodata['Specific Discharge [mm/day]'] = (
            self.hydrodata['Discharge [m^3/s]']
            / (self.drainage_basin_area__km2*1E3) * 86400)

        # Create columns for model output
        self.hydrodata['Specific Discharge (modeled) [mm/day]'] = pd.NA
        self.hydrodata['Snowpack (modeled) [mm SWE]'] = pd.NA
        self.hydrodata['Subsurface storage (modeled total) [mm]'] = pd.NA

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

        # Compute the water years based on the input month for
        # water-year rollover
        self.compute_water_year()

        # Scale evapotranspiration to enable water balance
        # We use this because ET estimates usually have much more
        # error than discharge, but may in the future want a way to
        # disable it.
        self.compute_ET_multiplier()
        self.compute_ET()

        # Model spin-up, if requested
        for _ in range(self.n_spin_up_cycles):
            self.run()  # Spin-up
            # Later allow user-selected starting and ending time steps
            # For now, just the whole series
            self._timestep_i = self.hydrodata.index[0]  # Restart

    def compute_water_year(self):
        """
        Assign a water-year label to each row in self.hydrodata.

        Adds a 'Water Year' column. A water year begins on
        water_year_start_month and is labelled by the calendar year in
        which it ends. For example, with a start month of October (USGS
        convention), 1 Oct 2020 – 30 Sep 2021 is water year 2021.
        """
        self.hydrodata['Water Year'] = pd.DatetimeIndex(self.hydrodata['Date']).year
        self.hydrodata['Water Year'] += \
            pd.DatetimeIndex(self.hydrodata['Date']).month >= self.water_year_start_month

    def compute_ET_multiplier(self):
        """
        Compute per-water-year ET scaling factors to enforce water balance.

        For each water year, computes the ratio of required ET (P - Q) to
        measured or computed ET, and stores this as 'ET multiplier' in
        self.hydrodata_WY_means. This multiplier is later applied in
        compute_ET() to scale ET so that P - Q - ET ≈ 0 over each water year.
        """
        # Originally used "sum", but then used "mean" so the headers would
        # still be sensible
        self.hydrodata_WY_means = self.hydrodata.groupby(
            self.hydrodata['Water Year']).mean(numeric_only=True)
        # Not needed, but no real harm in calculating
        self.hydrodata_WY_means['Runoff ratio'] = (
            self.hydrodata_WY_means['Specific Discharge [mm/day]']
            / self.hydrodata_WY_means['Precipitation [mm/day]'])
        _ET_required = -(self.hydrodata_WY_means['Specific Discharge [mm/day]'] -
                         self.hydrodata_WY_means['Precipitation [mm/day]'])
        self.hydrodata_WY_means['ET multiplier'] = (
            _ET_required / self.hydrodata_WY_means['Evapotranspiration [mm/day]'])

    def compute_ET(self):
        """
        Build the water-balance-adjusted ET time series used in the model.

        Assembles ET in two steps:

        1. Obtain raw daily ET from the input data file or the
           Thornthwaite–Chang 2019 equation (see evapotranspirationChang2019()).
        2. Scale raw ET by the per-water-year multiplier from
           compute_ET_multiplier() so that P - Q - ET ≈ 0 in each water year.

        The result is stored as 'ET for model [mm/day]' in self.hydrodata.
        """
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
            self.hydrodata_WY_means['ET multiplier'],
            on='Water Year')
        self.hydrodata['ET for model [mm/day]'] = (
            _raw_ET * self.hydrodata['ET multiplier'])

    def update(self, time_step=None):
        """
        Advance the model by one time step.

        Routes precipitation minus ET through the snowpack (if present) and
        then through each subsurface reservoir in order from shallowest to
        deepest. Stores modeled specific discharge, snowpack SWE, and total
        subsurface storage in self.hydrodata for the current time step.
        Part of the CSDMS Basic Model Interface.

        NOTE FALLACY: recharging before discharging,
        even though during the same time step.
        Consider changing to use half-recharge from each time step.

        FOR LATER: , dt_at_timestep=self.dt
        FOR SOONER: WATER-YEAR BALANCE

        Parameters
        ----------
        time_step : int, optional
            Index into self.hydrodata for the time step to update. If None,
            uses and then increments the internal counter self._timestep_i.
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
        if self.has_snowpack:
            self.snowpack.set_temperature(
                    self.hydrodata['Mean Temperature [C]'][time_step])
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
        # Frozen ground index (FGI): accumulates freezing degree-days,
        # decays symmetrically during warming.  When FGI > fdd_threshold the
        # top reservoir's infiltration to deeper layers is blocked (all drainage
        # becomes direct runoff), simulating frozen-soil reduction of
        # infiltration capacity.
        #
        # FGI(t) = max(0, FGI(t-1) - T_mean(t))
        #   T_mean < 0  → FGI rises  (freezing degree-days accumulate)
        #   T_mean > 0  → FGI falls  (thawing erodes the index)
        #
        # References
        # ----------
        # Molnau, M. and Bissell, V. C. (1983). A continuous frozen ground
        #     index for flood forecasting. Proc. 51st Annual Western Snow
        #     Conference, 109–119.
        #     https://westernsnowconference.org/sites/westernsnowconference.org/
        #     PDFs/1983Molnau.pdf
        # Shanley, J. B. and Chalmers, A. (1999). The effect of frozen soil on
        #     snowmelt runoff at Sleepers River, Vermont. Hydrol. Process.,
        #     13(12–13), 1843–1857.
        #     https://doi.org/10.1002/(SICI)1099-1085(199909)13:12/13
        #     <1843::AID-HYP879>3.0.CO;2-G
        # Dunne, T. and Black, R. D. (1971). Runoff processes during snowmelt.
        #     Water Resour. Res., 7(5), 1160–1172.
        #     https://doi.org/10.1029/WR007i005p01160
        _f0_calibrated = self.reservoirs[0].f_to_discharge
        if self.fdd_threshold < np.inf and self.has_snowpack:
            T_mean     = self.hydrodata['Mean Temperature [C]'][time_step]
            self._fgi  = max(0.0, self._fgi - T_mean)
            if self._fgi > self.fdd_threshold:
                # Block infiltration: all top-reservoir drainage to runoff
                self.reservoirs[0].f_to_discharge = 1.0

        # First, compute snowpack (if being used) and direct discharge
        for i in range(len(self.reservoirs)):
            # Top layer is special: snowmelt and/or precip infiltrates
            # immediately (at least on our daily time scales)
            if i == 0:
                if self.has_snowpack:
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
                self.reservoirs[i].recharge(
                    self.reservoirs[i-1].H_infiltrated
                    + self.reservoirs[i-1].H_deficit)
            self.reservoirs[i].discharge(self.dt)
            qi += self.reservoirs[i].H_discharge
        # Restore calibrated f_to_discharge for the top reservoir
        self.reservoirs[0].f_to_discharge = _f0_calibrated

        # Carry any unmet deficit forward to the next time step
        self.H_deficit = self.reservoirs[-1].H_deficit

        self.hydrodata.at[time_step, 'Specific Discharge (modeled) [mm/day]'] = qi
        if self.has_snowpack:
            self.hydrodata.at[time_step, 'Snowpack (modeled) [mm SWE]'] = self.snowpack.Hwater
        self.hydrodata.at[time_step, 'Subsurface storage (modeled total) [mm]'] = (
            np.sum([res.Hwater for res in self.reservoirs]))

    def evapotranspirationChang2019(self, Tmax=None, Tmin=None, photoperiod=None,
                                    k=0.69):
        """
        Modified daily Thornthwaite ET₀ equation.

        Chang et al. (2019), Eq. 1–4. https://doi.org/10.1002/ird.2309

        Parameters
        ----------
        Tmax : array-like
            Daily maximum temperature (°C).
        Tmin : array-like
            Daily minimum temperature (°C).
        photoperiod : array-like
            Photoperiod N (hours), computed from latitude and Julian day
            per Allen et al. (1998), Eqs. 2–4 of Chang et al. (2019).
        k : float
            Calibration coefficient in the T_ef formula. Default 0.69,
            recommended by Pereira & Pruitt (2004) for daily ET₀
            (https://doi.org/10.1016/j.agrformet.2004.01.005).
            Use 0.72 for monthly ET₀ per Camargo et al. (1999).

        Returns
        -------
        ET0 : array-like
            Daily reference evapotranspiration (mm day⁻¹).
        """
        if Tmax is None:
            Tmax = self.hydrodata['Maximum Temperature [C]']
        if Tmin is None:
            Tmin = self.hydrodata['Minimum Temperature [C]']
        if photoperiod is None:
            photoperiod = self.hydrodata['Photoperiod [hr]']

        Tef = 0.5 * k * (3 * Tmax - Tmin)
        C = photoperiod / 360.

        quadratic  = C * (-415.85 + 32.24 * Tef - 0.43 * Tef**2)
        power_law  = 16. * C * (10. * Tef / self.Chang_I) ** self.Chang_a

        ET0 = np.where(Tef > 26,  quadratic,
                       np.where(Tef > 0,   power_law,
                                0.))
        return ET0

    def run(self):
        """
        Advance the model through all time steps in self.hydrodata.

        Calls update() sequentially for every row in self.hydrodata, beginning
        from self._timestep_i. Part of the CSDMS Basic Model Interface.
        """
        for _ in self.hydrodata.index:
            self.update()

    def finalize(self):
        """
        Report model skill and display output plots.

        Calls computeNSE(verbose=True) to print the Nash–Sutcliffe Efficiency
        to stdout, then calls plot() to display a time-series comparison of
        observed and modeled specific discharge. Part of the CSDMS Basic Model
        Interface.
        """
        # Goodness of fit
        # Add options to print and/or save values later
        self.computeNSE(verbose=True)
        # Plot
        # Add flag for plotting (or not) later
        self.plot()

    def plot(self):
        """
        Display a time-series comparison of precipitation and specific discharge.

        Produces a dual-axis figure: precipitation as a bar chart on the left
        axis and both observed and modeled specific discharge as line plots on
        the right axis.
        """
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%y-%m-%d'))
        plt.xlabel('Date', fontsize=14)
        plt.xticks(rotation=45, horizontalalignment='right')
        plt.ylabel('Precipitation [mm/day]', fontsize=14, color='C0')
        plt.bar(self.hydrodata['Date'].values,
                height=self.hydrodata['Precipitation [mm/day]'].values/self.dt,
                width=1., align='center', label='Precipitation [mm/day]',
                linewidth=0, color='C0', alpha=0.5)  # C0 is the default bar-plot color
        plt.twinx()
        plt.plot(self.hydrodata['Date'].values,
                 self.hydrodata['Specific Discharge [mm/day]'].values,
                 'royalblue', label='Data', linewidth=2, alpha=0.8)
        plt.plot(self.hydrodata['Date'].values,
                 self.hydrodata['Specific Discharge (modeled) [mm/day]'].values,
                 'k', label='Model', linewidth=2, alpha=0.8)
        plt.ylim(0, plt.ylim()[-1])
        plt.legend(title='Specific Discharge', fontsize=11,
                   title_fontsize=11, labelcolor='linecolor')
        plt.ylabel('Specific Discharge [mm/day]', fontsize=14, color='0.3')
        plt.tight_layout()
        plt.show()

    def check_mass_balance(self, time_step=None):
        """
        Compute the mass-balance discrepancy at a given time step.

        Compares cumulative inputs (P - ET) from the start of the record
        through time_step with cumulative outputs (discharge) plus current
        storage (snowpack + subsurface reservoirs) and any carried-over
        deficit. Returns the excess mass still in the model; a value near
        zero indicates good mass conservation.

        Parameters
        ----------
        time_step : int, optional
            Row index in self.hydrodata at which to evaluate the balance.
            Defaults to the last row.

        Returns
        -------
        excess_mass_in_model : float
            Excess water remaining in the model budget (mm). Should be ~0
            for a mass-conserving run.
        """
        if time_step is None:
            time_step = self.hydrodata.index[-1]
        # Additions equals discharge out; set up this way, and can check.
        total_additions = \
            self.hydrodata['Precipitation [mm/day]'][:time_step+1].sum() \
            - self.hydrodata['ET for model [mm/day]'][:time_step+1].sum()
        # Storage reservoirs; snowpack is 0 when not simulated
        snow_storage = (self.hydrodata['Snowpack (modeled) [mm SWE]'][time_step]
                        if self.has_snowpack else 0.)
        subsurface_storage = self.hydrodata['Subsurface storage (modeled total) [mm]'][time_step]
        # Mass removal
        outlet_discharge = self.hydrodata[
            'Specific Discharge (modeled) [mm/day]'][:time_step+1].sum()
        # Carried-over water deficit
        deficit = self.H_deficit

        # Discrepancy
        excess_mass_in_model = (outlet_discharge + subsurface_storage
                                + snow_storage - total_additions + deficit)

        return excess_mass_in_model

    def computeNSE(self, return_nse=True, verbose=False):
        """
        Compute the Nash–Sutcliffe Efficiency of the discharge simulation.

        Compares modeled and observed specific discharge for all time steps
        where both values are non-missing. Stores the result as self.NSE.

        Parameters
        ----------
        return_nse : bool, optional
            If True (default), return the NSE value.
        verbose : bool, optional
            If True, print the NSE value to stdout.

        Returns
        -------
        NSE : float or None
            Nash–Sutcliffe Efficiency coefficient. Returns None if
            return_nse is False. A value of 1 indicates perfect agreement;
            values below 0 indicate the model performs worse than the
            observed-mean predictor.
        """

        # Shorthand for fcn
        q_data = self.hydrodata['Specific Discharge [mm/day]']
        q_model = self.hydrodata['Specific Discharge (modeled) [mm/day]']

        # Calculate NSE
        _realvalue = ~q_model.isna() & ~q_data.isna()
        NSE_num = np.sum((q_model[_realvalue] - q_data[_realvalue])**2)
        NSE_denom = np.sum((q_data[_realvalue] -
                            np.mean(q_data[_realvalue]))**2)
        if np.sum(~_realvalue):
            print("Excluded", np.sum(~_realvalue), "no-data points from NSE calculation")

        self.NSE = 1 - NSE_num / NSE_denom

        if verbose:
            print("NSE:", self.NSE)

        if return_nse:
            return self.NSE


def main():
    parser = argparse.ArgumentParser(
        description='Pass the configuration file path to run hydroRaVENS.')
    parser.add_argument('-y', '--configfile', type=str,
                        help='YAML file from which all inputs are read.')

    # Parse args if anything is passed.
    # If nothing is passed, then print help and exit.
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    b = Buckets()
    b.initialize(args.configfile)
    b.run()
    b.finalize()


if __name__ == "__main__":
    main()
