#! /usr/bin/env python

# Started by A. Wickert
# 25 July 2019
# Updated by J. Jones
# 08 Oct 2019

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
import pandas

#plt.ion()

#simulated rainfall series
##rain = poisson(.2, 100) # Convert to an import for real data

#import observed rainfall, currently following very specific format derived from NOAA NCEI (formerly NCDC) daily weather station

colnames = ['date', 'precip', 'snow', 'tmax', 'tmin', 'precip_mm', 'photoperiod', 'tmax_C', 'tmin_C']
rain_test = pandas.read_csv('C:/Users/Jabari Jones/Documents/UMN/Cannon/bucket_hydrology/albert_lea_daily_precip_2000.csv', names = colnames)

#extract rainfall, temperature, and photoperiod
precip_mm = rain_test.precip_mm.tolist()
tmax = rain_test.tmax_C.tolist()
tmin = rain_test.tmin_C.tolist()
photoperiod = rain_test.photoperiod.tolist()

#rain variable for imput
rain = map(float, precip_mm)


#import observed streamflow to compare w/model, currently using USGS daily data w/Q in mm/day

colnames_Q = ['USGS', 'code', 'date', 'Q', 'Q_mm']
Q_test = pandas.read_csv('C:/Users/Jabari Jones/Documents/UMN/Cannon/bucket_hydrology/cannon_welch_daily_2000.csv', names = colnames_Q)
Q_mm = Q_test.Q_mm.tolist()

#Timestep in days
dt = 1.

# Partitioning -- thse should also be external parameters
f_surface_to_discharge = 0.3

###
###I is a climatic parameter based on average monthly temps used to estimate daily PET. Method is detailed in Chang et al. 2019.
##FORECAST OF DAILY REFERENCE EVAPOTRANSPIRATION USING A MODIFIED DAILY THORNTHWAITE EQUATION AND TEMPERATURE FORECASTS.Irrigation and drainage. 68.
I = 41.47044637

#Create empty lists for outputs
Q = []
ET = []
T_eff = []
rain_eff = []

# This includes made-up rules about how much of the discharge from each
# bucket goes to the other buckets or to the discharge at the outlet.

#Model spin up
res_surface = reservoir(ai, 5000.)
res_subsurface = reservoir(bi, np.inf)
f_surface_to_discharge = ci

for ti in range(len(rain)):
    Qi = 0.
    Q_dif = 0.

    #Daily PET using method from Chang et al. 2019

    T_eff_i = 0.5*0.69*(3*tmax[ti]-tmin[ti])
    C = photoperiod[ti]/360
    a_i = 6.75*10**-7*I**3 - 7.72*10**-5*I**2 + 1.7912*10**-2*I + 0.49239

    if T_eff_i >= 26:
        ET_i = C*(-415.85 + 32.24*T_eff_i - 0.43*T_eff_i**2)
    elif 0< T_eff_i< 26:
        ET_i = 16*C*(10*T_eff_i/I)**a_i
    else:
        ET_i = 0.

#subtract ET from rainfall to calculate effective rainfall series, might try simple % reduction for interception
    rain_eff_i = rain[ti]-ET_i
    if rain_eff_i > 0:
        rain_eff_i = rain_eff_i *0.5
    else:
        rain_eff_i = 0.

    # Surface
    res_surface.recharge( rain_eff_i )
    res_surface.discharge(dt)
    Qi += f_surface_to_discharge * res_surface.Vout


    # Deep
    res_deep.recharge( (1. - f_surface_to_discharge) * res_surface.Vout )
    res_deep.discharge(dt)
    Qi += res_deep.Vout

    #Nash Sutcliffe Efficiency to measure model offset/error

    numerator = (Qi-Q_mm[ti])**2
    denominator = (Q_mm[ti]-np.nanmean(Q_mm))**2
    NSE_num.append(numerator)
    NSE_denom.append(denominator)

    # Append values into time series
    Q.append(Qi)

    ET.append(ET_i)
    T_eff.append(T_eff_i)
    rain_eff.append(rain_eff_i)

    T_eff_i = 0
    ET_i = 0


#Pseudocode for this NSE section
#Create empty elements to be filled by outputs

NSE = []
a_all = []
b_all = []
c_all = []
#still need to figure out how to store the elements asi:
#(ai_1, bi_1, ci_1, NSE_1)
#(ai_1, bi_1, ci_2, NSE_2)
counter = 1

for ai in range(1,20, 1):
    for bi in range(100, 1000, 100):
        for ci in np.arange(0.2, 0.7, 0.1):
            print counter
            counter += 1
            # Change these parameters with an input file or script, eventually
            # Units of mm;
            res_surface = reservoir(ai, 5000.)
            res_subsurface = reservoir(bi, np.inf)
            f_surface_to_discharge = ci

            Q = []

            NSE_num = []
            NSE_denom = []

            # Q calculations
            for ti in range(len(rain)):
                Qi = 0.

                #Daily PET using method from Chang et al. 2019

                T_eff_i = 0.5*0.69*(3*tmax[ti]-tmin[ti])
                C = photoperiod[ti]/360
                a_i = 6.75*10**-7*I**3 - 7.72*10**-5*I**2 + 1.7912*10**-2*I + 0.49239

                if T_eff_i >= 26:
                    ET_i = C*(-415.85 + 32.24*T_eff_i - 0.43*T_eff_i**2)
                elif 0< T_eff_i< 26:
                    ET_i = 16*C*(10*T_eff_i/I)**a_i
                else:
                    ET_i = 0.

                #subtract ET from rainfall to calculate effective rainfall series, might try simple % reduction for interception
                rain_eff_i = rain[ti]-ET_i
                if rain_eff_i > 0:
                    rain_eff_i = rain_eff_i *0.5
                else:
                    rain_eff_i = 0.

                # Surface
                res_surface.recharge( rain_eff_i)
                res_surface.discharge(dt)
                Qi += f_surface_to_discharge * res_surface.Vout

                # Deep
                res_subsurface.recharge( (1. - f_surface_to_discharge) * res_surface.Vout )
                res_subsurface.discharge(dt)
                Qi += res_subsurface.Vout

                # Append the time step's discharge to the time series
                Q.append(Qi)

            #Nash Sutcliffe Efficiency to measure model offset/error

                numerator = (Qi-Q_mm[ti])**2
                denominator = (Q_mm[ti]-np.nanmean(Q_mm))**2
                NSE_num.append(numerator)
                NSE_denom.append(denominator)

            # Append values into time series
            NSE.append(1 - (np.nansum(NSE_num)/np.nansum(NSE_denom)))
            a_all.append(ai)
            b_all.append(bi)
            c_all.append(ci)

NSE = np.array(NSE)
NSEmax_i = np.where(NSE == np.max(NSE))[0][0]
print "Max NSE:", np.max(NSE)
print "a =", a_all[NSEmax_i]
print "b =", b_all[NSEmax_i]
print "c =", c_all[NSEmax_i]




'''
                append i, ti, mi, NSE to something (first row of a matrix? Or use i, ti, and mi as
                the indexes for a 3D matrix/array. Actually on second thought, that won't work cause the
                parameters won't necessarily be integers, which I assume Python would want for index #s.
                So instead, we'll go back to that first idea, and have something (it would be a dataframe in R)
                that looks like i = (i_1,ti_1,mi_1,NSE)
                                    (i_1,ti_1,mi_2,NSE):
                just as a heads up, this will be an enormous thing

                find and print the row in which NSE is max/RMSE is min and use that
                '''


           #sub


#Generate p    lot
plt.figure(    )
#plt.plot(rain_eff, 'g', label='rainfall')
#plt.plot(r    ain, 'g', label = 'rain')
plt.plot(Q,     'b', label='discharge')
plt.plot(Q_mm, 'r', label = 'obs Q')
#plt.plot(ET    , 'black', label = 'ET')
plt.legend(    fontsize=11)
plt.ylabel(    'Rainfall or discharge [mm/day]', fontsize=14)
plt.xlabel(    'Time [units arbitrary]', fontsize=14)
#plt.savefi    g('TDF-'+'%.1f' %f_tile+'.png')
plt.show()














