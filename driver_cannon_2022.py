#! /usr/bin/python3

# driver_cannon_2022.py

import bucketHydrology
b = bucketHydrology.Buckets()
b.initialize('cannon_cfg_tmp.yml')

# To lump into internal functions later but test here for now
self = b
#self.hydrodata['Mean Temperature [C]'] += 40 # ensure no snowmelt processes
self.compute_water_year()
self.compute_ET_multiplier()
self.compute_ET()
#self.run() #Spin-up
#self._timestep_i = 0. # Restart
self.run()
self.computeNSE(verbose=True)
self.plot()

import numpy as np
for time_step in self.hydrodata.index:
    excess_mass_in_model = self.check_mass_balance(time_step)
    subsurface_storage = self.hydrodata['Subsurface storage (modeled total) [mm]'][time_step]
    print( '%0.3f' %excess_mass_in_model, '%0.3f' %subsurface_storage )
