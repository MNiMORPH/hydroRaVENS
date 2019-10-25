# Program below the class
import numpy as np
from numpy.random import poisson
from matplotlib import pyplot as plt
import bucketHydrology as bh
import pandas as pd
#plt.ion()

###############
# IMPORT DATA #
###############

#import observed rainfall, currently following very specific format derived from NOAA NCEI (formerly NCDC) daily weather station
colnames = ['date', 'precip', 'snow', 'tmax', 'tmin', 'precip_mm', 'photoperiod', 'tmax_C', 'tmin_C']
weather = pd.read_csv('albert_lea_daily_precip_2000.csv', names = colnames)

#import observed streamflow to compare w/model, currently using USGS daily data w/Q in mm/day
colnames_Q = ['USGS', 'code', 'date', 'Q', 'Q_mm']
Q_measured = pd.read_csv('cannon_welch_daily_2000.csv', names = colnames_Q)


#################################
# SET RESERVOIRS AND PARAMETERS #
#################################

dt = 1. # day
        
res_surface = bh.reservoir(t_efold=5., f_to_discharge=0.4, Hmax=np.inf)
res_deep = bh.reservoir(t_efold=100., f_to_discharge=1., Hmax=np.inf)

strat_column = [res_surface, res_deep]

#############
# RUN MODEL #
#############

# Initialize
watershed = bh.buckets(reservoir_list=strat_column, dt=dt)
watershed.evapotranspirationChang2019(weather['tmax'], weather['tmin'], weather['photoperiod'])

# Run
watershed.run(rain=weather['precip_mm'], ET=True)

# Finalize
watershed.computeNashSutcliffeEfficiency(Q_measured['Q_mm'])
print "NSE:", watershed.NSE
watershed.plot(Q_measured['Q_mm'])
