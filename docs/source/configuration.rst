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
        # Required section; snowpack processes are only active if
        # Mean Temperature [C] is present in the input CSV

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
   * - ``water_reservoir_effective_depths__mm``
     - list of floats
     - Initial water depth (mm) in each reservoir, listed top to bottom
   * - ``snowpack__mm_SWE``
     - float
     - Initial snowpack depth in mm snow-water equivalent (SWE)

Example:

.. code-block:: yaml

    initial_conditions:
        water_reservoir_effective_depths__mm:
            - 150    # Top (soil) reservoir starts at 150 mm
            - 400    # Bottom (groundwater) starts at 400 mm
        snowpack__mm_SWE: 50  # 50 mm SWE initially

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
   * - ``spin_up_cycles``
     - int
     - Number of complete passes through data before the main run. Use 0 to
       skip spin-up.
   * - ``scale_et``
     - bool
     - Scale ET by a per-water-year multiplier so that P − Q − ET = 0 over
       each water year. Default ``true``. Set to ``false`` only when
       supplying trusted measured ET (e.g. eddy covariance) that should not
       be corrected. Using ``false`` with ``ThorntwaiteChang2019`` will raise
       a warning because Thornthwaite ET carries large systematic biases.

Example:

.. code-block:: yaml

    general:
        spin_up_cycles: 2  # Run through data twice to initialize
        scale_et: true     # default; omit to accept the default

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
   * - ``maximum_effective_depths__mm``
     - list of floats
     - Storage capacity (mm) per reservoir; use ``.inf`` for unlimited

Example:

.. code-block:: yaml

    reservoirs:
        e_folding_residence_times__days:
            - 16      # Soil: drains in ~16 days
            - 2000    # Groundwater: drains in ~2000 days
        exfiltration_fractions:
            - 0.8     # 80% to river, 20% to next layer
            - 1.0     # 100% to river (should be 1.0; warning issued if not)
        maximum_effective_depths__mm:
            - .inf
            - .inf

The ``snowmelt`` section
~~~~~~~~~~~~~~~~~~~~~~~~

This section is **always required** in the configuration file. Snowpack
*processes* are only activated when ``Mean Temperature [C]`` is present in
the input CSV, but the ``PDD_melt_factor`` key must be present regardless.

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
        water_reservoir_effective_depths__mm: [2, 400]
        snowpack__mm_SWE: 0

    catchment:
        drainage_basin_area__km2: 3800
        evapotranspiration_method: ThorntwaiteChang2019  # requires T_monthly_normals in Python
        water_year_start_month: 10

    general:
        spin_up_cycles: 1

    reservoirs:
        e_folding_residence_times__days: [16, 2000]
        exfiltration_fractions: [0.8, 1.0]
        maximum_effective_depths__mm: [.inf, .inf]
    
    snowmelt:
        PDD_melt_factor: 1.0

Illustrative Starting Points
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These are rough starting points only, not calibrated or validated parameter
sets. Actual values depend strongly on the specific catchment.

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
