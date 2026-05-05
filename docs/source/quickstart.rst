Quick Start
===========

This guide will get you running your first HydroRaVENS model.

Prepare Your Data
~~~~~~~~~~~~~~~~~

HydroRaVENS requires daily time series data in CSV format. Required and optional columns:

.. list-table:: Input CSV Columns
   :widths: 35 15 40
   :header-rows: 1

   * - Column Name
     - Units
     - When required
   * - ``Date``
     - YYYY-MM-DD
     - Always; must be continuous daily data with no gaps
   * - ``Precipitation [mm/day]``
     - mm/day
     - Always
   * - ``Discharge [m^3/s]``
     - m³/s
     - Always; used to compute NSE
   * - ``Mean Temperature [C]``
     - °C
     - Snowpack processes
   * - ``Minimum Temperature [C]``
     - °C
     - ``ThorntwaiteChang2019`` ET method
   * - ``Maximum Temperature [C]``
     - °C
     - ``ThorntwaiteChang2019`` ET method
   * - ``Photoperiod [hr]``
     - hours
     - ``ThorntwaiteChang2019`` ET method
   * - ``Evapotranspiration [mm/day]``
     - mm/day
     - ``datafile`` ET method

Example input (first few rows):

.. code-block:: text

    Date,Precipitation [mm/day],Discharge [m^3/s],Mean Temperature [C],Evapotranspiration [mm/day]
    2010-01-01,0.0,15.2,-2.5,0.2
    2010-01-02,2.1,16.8,-1.3,0.2
    2010-01-03,0.5,15.9,0.2,0.3
    2010-01-04,5.8,22.1,3.1,0.5

Create a Configuration File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

HydroRaVENS is configured through a YAML file (``config.yml``):

.. code-block:: yaml

    timeseries:
        datafile: input_data.csv

    initial_conditions:
        water_reservoir_effective_depths__m:
            - 2      # Top (soil) reservoir
            - 400    # Bottom (groundwater) reservoir
        snowpack__m_SWE: 0

    catchment:
        drainage_basin_area__km2: 3800
        evapotranspiration_method: datafile
        water_year_start_month: 10

    general:
        scalar_dt: true
        spin_up_cycles: 1

    reservoirs:
        e_folding_residence_times__days:
            - 16      # Soil: fast response
            - 2000    # Groundwater: slow response
        exfiltration_fractions:
            - 0.8     # 80% to discharge, 20% infiltrates
            - 1.0     # 100% to discharge (bottom layer)
        maximum_effective_depths__m:
            - .inf
            - .inf

    snowmelt:
        PDD_melt_factor: 1.0

Run the Model
~~~~~~~~~~~~~

**Using the Python API:**

.. code-block:: python

    import hydroravens

    model = hydroravens.Buckets()
    model.initialize('config.yml')
    model.run()
    model.computeNSE(verbose=True)
    model.plot()

**Using the command-line interface:**

.. code-block:: bash

    hydroravens -y config.yml

Adjust Parameters
~~~~~~~~~~~~~~~~~

Model performance depends on the reservoir parameters:

**Residence times** (``e_folding_residence_times__days``)
  Larger values = slower response. Typical ranges:

  * Soil zone: 5--50 days (fast response)
  * Groundwater: 100--5000 days (slow, baseflow response)

**Exfiltration fractions** (``exfiltration_fractions``)
  Fraction of each reservoir's drainage going directly to the river.

  * Higher = more direct runoff
  * Lower = more infiltration to deeper layers
  * Bottom layer should be 1.0 (all to discharge)

**Initial depths** (``water_reservoir_effective_depths__m``)
  Starting water content in each reservoir. Spin-up cycles reduce
  sensitivity to these initial values.

**Spin-up cycles** (``spin_up_cycles``)
  Number of passes through the full record before the main run.
  One cycle is usually sufficient.

Next Steps
~~~~~~~~~~

* Read the :doc:`model_description` for the theory behind each component
* Explore :doc:`configuration` for all configuration options
* See the :doc:`api` for full Python class documentation
