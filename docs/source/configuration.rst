Configuration Reference
========================

HydroRaVENS is configured entirely through a YAML file. This page documents all 
available options.

Configuration File Structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

    timeseries:
        # Data input section
    
    initial_conditions:
        # Starting state of reservoirs and snowpack
    
    catchment:
        # Basin properties
    
    general:
        # Simulation settings
    
    reservoirs:
        # Reservoir cascade configuration
    
    snowmelt:
        # (Optional) Snow processes

The ``timeseries`` section
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Option
     - Type
     - Description
   * - ``datafile``
     - string
     - Path to CSV input file with daily time series

Example:

.. code-block:: yaml

    timeseries:
        datafile: data/streamflow_2010_2020.csv

The ``initial_conditions`` section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 35 15 40
   :header-rows: 1

   * - Option
     - Type
     - Description
   * - ``water_reservoir_effective_depths__m``
     - list of floats
     - Initial water depth (m) in each reservoir, listed top to bottom
   * - ``snowpack__m_SWE``
     - float
     - Initial snowpack depth in m snow-water equivalent (SWE)

Example:

.. code-block:: yaml

    initial_conditions:
        water_reservoir_effective_depths__m:
            - 1.5    # Top (soil) reservoir starts at 1.5 m
            - 200    # Bottom (groundwater) starts at 200 m
        snowpack__m_SWE: 0.5  # 0.5 m SWE initially

The ``catchment`` section
~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 35 15 40
   :header-rows: 1

   * - Option
     - Type
     - Description
   * - ``drainage_basin_area__km2``
     - float
     - Basin area (km²); used for discharge-to-depth conversion
   * - ``evapotranspiration_method``
     - string
     - ``datafile`` or ``ThorntwaiteChang2019``
   * - ``water_year_start_month``
     - int
     - Month (1–12) when water year begins (10 = October for USGS)

Example:

.. code-block:: yaml

    catchment:
        drainage_basin_area__km2: 3800
        evapotranspiration_method: datafile
        water_year_start_month: 10

The ``general`` section
~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Option
     - Type
     - Description
   * - ``scalar_dt``
     - bool
     - Required field; set to ``true``. The model always enforces a 1-day
       time step regardless of this value.
   * - ``spin_up_cycles``
     - int
     - Number of complete passes through data before the main run. Use 0 to
       skip spin-up.

Example:

.. code-block:: yaml

    general:
        scalar_dt: true
        spin_up_cycles: 2  # Run through data twice to initialize

The ``reservoirs`` section
~~~~~~~~~~~~~~~~~~~~~~~~~~

This section defines the cascade of linear reservoirs (1 or more).
All lists must have the same length.

.. list-table::
   :widths: 35 15 40
   :header-rows: 1

   * - Option
     - Type
     - Description
   * - ``e_folding_residence_times__days``
     - list of floats
     - Drainage time constant (days) for each reservoir
   * - ``exfiltration_fractions``
     - list of floats
     - Fraction (0–1) of drainage exiting as discharge
   * - ``maximum_effective_depths__m``
     - list of floats
     - Storage capacity (m) per reservoir; use ``.inf`` for unlimited

Example:

.. code-block:: yaml

    reservoirs:
        e_folding_residence_times__days:
            - 16      # Soil: drains in ~16 days
            - 2000    # Groundwater: drains in ~2000 days
        exfiltration_fractions:
            - 0.8     # 80% to river, 20% to next layer
            - 1.0     # 100% to river (should be 1.0; warning issued if not)
        maximum_effective_depths__m:
            - .inf
            - .inf

The ``snowmelt`` section (optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Only used if mean temperature is in the input CSV.

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Option
     - Type
     - Description
   * - ``PDD_melt_factor``
     - float
     - Melt rate (mm SWE per °C per day)

Example:

.. code-block:: yaml

    snowmelt:
        PDD_melt_factor: 1.0  # 1 mm SWE melts per °C per day

Complete Example
~~~~~~~~~~~~~~~~

.. code-block:: yaml

    timeseries:
        datafile: data/input.csv
    
    initial_conditions:
        water_reservoir_effective_depths__m: [2, 400]
        snowpack__m_SWE: 0
    
    catchment:
        drainage_basin_area__km2: 3800
        evapotranspiration_method: ThorntwaiteChang2019
        water_year_start_month: 10
    
    general:
        scalar_dt: true
        spin_up_cycles: 1
    
    reservoirs:
        e_folding_residence_times__days: [16, 2000]
        exfiltration_fractions: [0.8, 1.0]
        maximum_effective_depths__m: [.inf, .inf]
    
    snowmelt:
        PDD_melt_factor: 1.0

Common Configuration Presets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Alpine watershed with snowpack:**

.. code-block:: yaml

    reservoirs:
        e_folding_residence_times__days: [8, 1500]
        exfiltration_fractions: [0.7, 1.0]
    
    snowmelt:
        PDD_melt_factor: 1.5
    
    catchment:
        evapotranspiration_method: ThorntwaiteChang2019

**Humid lowland catchment:**

.. code-block:: yaml

    reservoirs:
        e_folding_residence_times__days: [25, 3000]
        exfiltration_fractions: [0.9, 1.0]
    
    catchment:
        evapotranspiration_method: datafile

**Arid/semi-arid watershed:**

.. code-block:: yaml

    reservoirs:
        e_folding_residence_times__days: [3, 500]
        exfiltration_fractions: [0.5, 1.0]
    
    snowmelt:
        PDD_melt_factor: 0.5
