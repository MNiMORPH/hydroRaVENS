HydroRaVENS Documentation
==========================

.. image:: https://zenodo.org/badge/199317220.svg
   :target: https://zenodo.org/badge/latestdoi/199317220
   :alt: Zenodo DOI

A simple, flexible reservoir-based hydrological model for water balance simulation and streamflow prediction.

**HydroRaVENS** (Rain and Variable Evapotranspiration, Nieve, and Streamflow) is a 
conceptual daily-timestep model that routes precipitation through an optional snowpack 
and a cascade of linear reservoirs. Ideal for ungauged basins, climate impact studies, 
and educational applications.

✨ **Key Features**
~~~~~~~~~~~~~~~~~~~

* 🏔️ **Optional snowpack module** – Positive-degree-day snowmelt calculations
* 💧 **Cascading linear reservoirs** – Stack multiple reservoirs from soil to groundwater
* 🌡️ **Flexible ET modeling** – Read from data or use Thornthwaite-Chang equation
* ⚖️ **Mass-conserving** – Strict water balance closure
* 📊 **Model evaluation** – Built-in Nash-Sutcliffe Efficiency (NSE) computation
* 🐍 **Python + CLI** – Use as library or command-line tool
* ✔️ **Lightweight** – Minimal dependencies (NumPy, Pandas, PyYAML)

Quick Links
~~~~~~~~~~~

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   configuration
   examples
   model_description

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api
   references

Quick Example
~~~~~~~~~~~~~

**Python API:**

.. code-block:: python

    import hydroravens
    
    model = hydroravens.Buckets()
    model.initialize('config.yml')
    model.run()
    print(f"NSE: {model.computeNSE():.3f}")
    model.plot()

**Command-line:**

.. code-block:: bash

    hydroravens -y config.yml

Getting Help
~~~~~~~~~~~~

* 🐛 **Report bugs:** `GitHub Issues <https://github.com/MNiMORPH/hydroRaVENS/issues>`_
* 💬 **Discuss:** `GitHub Discussions <https://github.com/MNiMORPH/hydroRaVENS/discussions>`_
* 📖 **Learn more:** `CSDMS Model Page <https://csdms.colorado.edu/wiki/Model:HydroRaVENS>`_

About
~~~~~

HydroRaVENS is developed and maintained by the `MNiMORPH group <https://github.com/MNiMORPH>`_.
It is published under the GNU General Public License v3.0.

**Citation:**

If you use HydroRaVENS in your research, please cite:

.. code-block:: bibtex

    @software{Wickert2026RaVENS,
      author = {Wickert, Andrew D.},
      doi = {10.5281/zenodo.6787390},
      month = May,
      title = {{RaVENS: Rain and Variable Evapotranspiration, Nieve, and Streamflow}},
      url = {https://github.com/MNiMORPH/RaVENS},
      version = {2.2.0},year = {2026}
    }
