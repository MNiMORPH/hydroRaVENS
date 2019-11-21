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
t_P_ET = pd.read_csv('CannonLivneh_Date_P_ET.csv')
t_P_ET['Date'] = t_P_ET['Date'].astype('datetime64[ns]') 
# Just 2000 for comparison with data
years = []
for date in t_P_ET['Date']:
    years.append(date.year)
years = np.array(years)
t_P_ET = t_P_ET[years == 2000]

#import observed streamflow to compare w/model, currently using USGS daily data w/Q in mm/day
colnames_Q = ['USGS', 'code', 'date', 'Q', 'Q_mm']
Q_measured = pd.read_csv('cannon_welch_daily_2000.csv', names = colnames_Q)


#################################
# SET RESERVOIRS AND PARAMETERS #
#################################

res_surface = bh.reservoir(t_efold=15., f_to_discharge=.8, Hmax=np.inf)
res_deep = bh.reservoir(t_efold=100., f_to_discharge=1., Hmax=np.inf)

strat_column = [res_surface, res_deep]

#############
# RUN MODEL #
#############

# Initialize
watershed = bh.buckets(reservoir_list=strat_column)
watershed.initialize(t_P_ET['Date'], Hlist = [9.4706537227959871, 11.24098140709903])
# ET from the model seems too high: Q (data) = 3 * P (data) - ET (reanalysis)
watershed.set_evapotranspiration( np.array(t_P_ET['Evapotranspiration [mm/day]'] * .68) )

# Run
watershed.run(rain=np.array(t_P_ET['Precipitation [mm/day]']), ET_flag=True)

# Finalize
watershed.computeNashSutcliffeEfficiency(Q_measured['Q_mm'])
print("NSE:", watershed.NSE)
Hout = watershed.export_Hlist()
print(Hout)
#watershed.plot()
watershed.plot(Q_measured['Q_mm'])

