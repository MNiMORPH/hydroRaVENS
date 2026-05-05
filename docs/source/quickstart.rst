Quick Start
===========

This 5-minute tutorial will get you running your first HydroRaVENS model.

Prepare Your Data
~~~~~~~~~~~~~~~~~

HydroRaVENS requires daily time series data in CSV format with these columns:

.. list-table:: Required CSV Columns
   :widths: 25 15 40
   :header-rows: 1

   * - Column Name
     - Units
     - Notes
   * - ``Date``
     - YYYY-MM-DD
     - Must be continuous daily data (no gaps)
   * - ``Precipitation [mm/day]``
     - mm/day
     - Always required
   * - ``Discharge [m^3/s]``
     - m³/s
     - Used to compute NSE (model skill)
   * - ``Mean Temperature [C]``
     - °C
     - Required only if snowpack processes needed
   * - ``Minimum Temperature [C]``
     - °C
     - Required for Thornthwaite-Chang ET method
   * - ``Maximum Temperature [C]``
     - °C
     - Required for Thornthwaite-Chang ET method
   * - ``Photoperiod [hr]``
     - hours
     - Required for Thornthwaite-Chang ET method
   * - ``Evapotranspiration [mm/day]``
     - mm/day
     - Alternative to computing ET; if present, used instead

Example input file (``input_data.csv``):

.. code-block:: csv

    Date,Precipitation [mm/day],Discharge [m^3/s],Mean Temperature [C]
    2010-01-01,0.0,15.2,-2.5
    2010-01-02,2.1,16.8,-1.3
    2010-01-03,0.5,15.9,0.2
    2010-01-04,5.8,22.1,3.1

Create a Configuration File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

HydroRaVENS is entirely configured through a YAML file (``config.yml``):

.. code-block:: yaml

    # Minimal configuration example
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
    
    # Create and initialize model
    model = hydroravens.Buckets()
    model.initialize('config.yml')
    
    # Run simulation
    model.run()
    
    # Evaluate model skill
    nse = model.computeNSE(verbose=True)
    print(f"Nash-Sutcliffe Efficiency: {nse:.3f}")
    
    # Generate plots
    model.plot()

**Using the command-line interface:**

.. code-block:: bash

    hydroravens -y config.yml

Expected Output
~~~~~~~~~~~~~~~

.. code-block:: text

    Loading configuration from config.yml...
    Initializing model...
    Running model for 3650 days...
    
    === Model Evaluation ===
    Nash-Sutcliffe Efficiency: 0.823
    Mean Absolute Error: 8.5 mm/day
    Bias: -0.3 mm/day

Adjust Parameters
~~~~~~~~~~~~~~~~~

Model performance depends on reservoir parameters:

**Residence times** (``e_folding_residence_times__days``)
  Larger values = slower response. Typical ranges:
  
  * Soil zone: 5–50 days (fast response)
  * Groundwater: 100–5000 days (slow, baseflow response)

**Exfiltration fractions** (``exfiltration_fractions``)
  What fraction of each reservoir drains to the river?
  
  * Higher = more direct runoff
  * Lower = more infiltration to deeper layers
  * Bottom layer must be 1.0 (all to discharge)

**Initial depths** (``water_reservoir_effective_depths__m``)
  Starting water content in each reservoir.

**Spin-up cycles** (``spin_up_cycles``)
  Number of passes through the data before the main run.
  Increase if initial conditions don't matter; set to 1 for quick runs.

Next Steps
~~~~~~~~~~

* 📖 Read the :doc:`model_description` for theory
* ⚙️ Explore :doc:`configuration` for all options
* 🔬 Check :doc:`examples` for calibration & sensitivity analysis
* 🐍 See the :doc:`api` for Python class documentation
