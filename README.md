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

HydroRaVENS is a simple "conceptual" hydrological model that routes daily precipitation through an optional snowpack stage and then through a cascade of one or more linked linear reservoirs (soil-zone water, groundwater, etc.). Evapotranspiration can be supplied from an external data file or computed from temperature using the Thornthwaite–Chang 2019 equation. A water-balance component scales ET (typically the most uncertain input) so that P − Q − ET ≈ 0 over each water year. Model skill is assessed with the Nash–Sutcliffe Efficiency (NSE).

## Model overview

Water moves through HydroRaVENS in three sequential stages each day.

**Snowpack (optional).** When mean air temperature T ≤ 0 °C, precipitation accumulates as snow (stored as snow-water equivalent, SWE). When T > 0 °C, snowpack melts at a positive-degree-day rate and the meltwater is routed to the top subsurface reservoir. Precipitation on warm days bypasses the snowpack and infiltrates directly.

**Linear reservoir cascade.** One or more subsurface reservoirs are stacked from shallowest (soil zone) to deepest (baseflow groundwater). Each reservoir drains exponentially with a characteristic e-folding residence time τ. At each time step a fraction *f* of drainage exits as river discharge; the remainder infiltrates to the next-deeper reservoir. Overflow above a maximum storage capacity exits immediately as direct runoff.

**Evapotranspiration and water balance.** ET is either read from the input file or estimated with the Thornthwaite–Chang 2019 equation (requiring daily Tmax, Tmin, and photoperiod). In either case, ET is scaled by a per-water-year multiplier so that long-term water balance is enforced.

## Installation

```
pip install hydroravens
```

To install from source for development:

```
git clone https://github.com/MNiMORPH/hydroRaVENS.git
cd hydroRaVENS
pip install -e .
```

## Quick start

**Python API**

```python
import hydroravens

b = hydroravens.Buckets()
b.initialize('config.yml')
b.run()
b.computeNSE(verbose=True)
b.plot()
```

**Command-line interface**

```
hydroravens -y config.yml
```

## Configuration file

The model is fully configured through a YAML file. A minimal two-reservoir example:

```yaml
timeseries:
    datafile: input_data.csv          # path to the CSV input file

initial_conditions:
    water_reservoir_effective_depths__m:
        - 2                           # top (shallowest) reservoir initial depth
        - 400                         # bottom reservoir initial depth
    snowpack__m_SWE: 0                # initial snowpack depth

catchment:
    drainage_basin_area__km2: 3800
    evapotranspiration_method: datafile   # or ThorntwaiteChang2019
    water_year_start_month: 10            # 10 = October (USGS convention)

general:
    scalar_dt: true                   # all time steps must be 1 day
    spin_up_cycles: 1                 # full-record loops before main run

# Reservoirs listed top (shallowest) to bottom (deepest)
reservoirs:
    e_folding_residence_times__days:
        - 16
        - 2000
    exfiltration_fractions:           # fraction of drainage going to river discharge
        - 0.8                         # 0.2 infiltrates to the next layer
        - 1.0                         # bottom layer: all to discharge
    maximum_effective_depths__m:
        - .inf
        - .inf

snowmelt:
    PDD_melt_factor: 1.               # mm SWE per positive °C per day
```

Key configuration options:

| Section | Key | Description |
|---------|-----|-------------|
| `timeseries` | `datafile` | Path to the CSV input file |
| `initial_conditions` | `water_reservoir_effective_depths__m` | Initial water depth in each reservoir (m) |
| `initial_conditions` | `snowpack__m_SWE` | Initial snowpack depth (m SWE) |
| `catchment` | `drainage_basin_area__km2` | Basin area used to convert discharge to specific discharge (mm/day) |
| `catchment` | `evapotranspiration_method` | `datafile` or `ThorntwaiteChang2019` |
| `catchment` | `water_year_start_month` | Month number at which a new water year begins |
| `general` | `spin_up_cycles` | Number of full-record spin-up passes before the main simulation |
| `reservoirs` | `e_folding_residence_times__days` | Exponential residence time for each reservoir (days) |
| `reservoirs` | `exfiltration_fractions` | Fraction of each reservoir's drainage routed to river discharge |
| `reservoirs` | `maximum_effective_depths__m` | Storage cap per reservoir (use `.inf` for no cap) |
| `snowmelt` | `PDD_melt_factor` | Positive-degree-day melt rate (mm SWE °C⁻¹ day⁻¹) |

## Input data format

Input is a daily CSV file. Required and optional columns:

| Column | Units | When required |
|--------|-------|---------------|
| `Date` | — | Always; must be daily with no gaps |
| `Precipitation [mm/day]` | mm/day | Always |
| `Discharge [m^3/s]` | m³/s | Always; used to compute NSE |
| `Mean Temperature [C]` | °C | Snowpack processes |
| `Minimum Temperature [C]` | °C | `ThorntwaiteChang2019` ET |
| `Maximum Temperature [C]` | °C | `ThorntwaiteChang2019` ET |
| `Photoperiod [hr]` | hours | `ThorntwaiteChang2019` ET |
| `Evapotranspiration [mm/day]` | mm/day | `datafile` ET method |

The model raises a `ValueError` if any time step is not exactly 1 day.

## Model description

HydroRaVENS is a lumped, daily-timestep conceptual hydrological model. All fluxes are expressed as depth over the drainage basin (mm/day), and the model conserves mass to within numerical precision over each water year.

### Snowpack

If `Mean Temperature [C]` is present in the input file, snowpack processes are enabled.

- **Accumulation:** When T ≤ 0 °C, P accumulates as SWE.
- **Melt:** When T > 0 °C, melt = min(SWE, PDD_melt_factor × T × dt). All melt is routed to the top reservoir; surface runoff from snowmelt is not modeled.
- **ET deficit:** A negative net water input (ET > P) first sublimates snow; any remaining deficit is passed downward to the reservoirs.

### Linear reservoirs

Each reservoir drains by first-order exponential decay:

    ΔH = H × (1 − exp(−dt / τ))

where H is the current water depth, dt is the time step (days), and τ is the e-folding residence time (days). Of the water drained, a fraction *f* (`exfiltration_fraction`) exits as river discharge, and 1 − *f* infiltrates to the next-deeper reservoir. If inflow exceeds the maximum storage capacity `Hmax`, the excess exits immediately as discharge.

Reservoirs are stacked top to bottom, and infiltration from one layer becomes the recharge input for the layer below. The bottom reservoir must fully discharge to the river (f = 1) to conserve mass; a warning is issued if this is not the case.

### Evapotranspiration

Two methods are supported:

- **`datafile`:** ET is read directly from the input CSV and scaled per water year so that P − Q − ET ≈ 0.
- **`ThorntwaiteChang2019`:** Daily ET₀ is estimated from Tmax, Tmin, and photoperiod using the modified Thornthwaite equation of Chang et al. (2019), https://doi.org/10.1002/ird.2309. Long-term monthly temperature normals must be passed to `Buckets(T_monthly_normals=...)` to compute the thermal index *I* and exponent *a* before `initialize()` is called.

In both cases, ET is multiplied by a per-water-year factor (derived from observed P and Q) so that the annual water balance closes. This corrects for systematic bias in ET estimates without altering the seasonal pattern.

### Model skill

`computeNSE()` computes the Nash–Sutcliffe Efficiency:

    NSE = 1 − Σ(Q_mod − Q_obs)² / Σ(Q_obs − Q̄_obs)²

NSE = 1 is a perfect simulation; NSE < 0 means the model performs worse than the observed mean as a predictor.
