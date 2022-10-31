#! /usr/bin/python3

# driver_cannon_2022.py

import bucketHydrology
b = bucketHydrology.Buckets()
b.initialize('cannon_cfg.yml')

# To lump into internal functions later but test here for now
self = b
self.compute_water_year()
self.compute_ET_multiplier()
self.compute_ET()
#self.run() #Spin-up
#self._timestep_i = 0. # Restart
self.run()
self.computeNSE(verbose=True)
self.plot()
