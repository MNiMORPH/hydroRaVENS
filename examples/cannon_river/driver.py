#! /usr/bin/python3

# driver_cannon_2022.py

import hydroravens

b = hydroravens.Buckets()
b.initialize('cannon_cfg.yml')
b.run()
b.compute_NSE(verbose=True)
b.plot()

#import numpy as np
#for time_step in b.hydrodata.index:
#    excess_mass_in_model = b.check_mass_balance(time_step)
#    subsurface_storage = b.hydrodata['Subsurface storage (modeled total) [mm]'][time_step]
#    print( '%0.3f' %excess_mass_in_model, '%0.3f' %subsurface_storage )
