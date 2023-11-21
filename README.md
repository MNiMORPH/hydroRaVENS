[![DOI](https://zenodo.org/badge/199317220.svg)](https://zenodo.org/badge/latestdoi/199317220)


# :black_nib: RaVENS
Rain and Variable Evapotranspiration, Nieve, and Streamflow

Package name:
```
hydroravens
```

Simple reservoir-based hydrological model

## Description

[CSDMS model page](https://csdms.colorado.edu/wiki/Model:HydroRaVENS)

Simple "conceptual" hydrological model that may include an arbitrary number of linked linear reservoirs (soil-zone water, groundwater, etc.) as well as snowpack (accumulation from precipitation with T<0; positive-degree-day melt) and evapotranspiration (from external input or Thorntwaite equation).

It also includes a water-balance component to adjust ET (typically the least known input) to ensure that P - Q - ET = 0 over the course of a water year.

Other components plot data and compute the NSE (Nash–Sutcliffe model efficiency coefficient). 
