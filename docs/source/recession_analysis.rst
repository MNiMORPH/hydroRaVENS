Recession Analysis
==================

HydroRaVENS includes tools for analysing the nonlinear storage–discharge
relationship directly from observed streamflow records, before running or
calibrating the model.

.. contents:: On this page
   :local:
   :depth: 2

Theory
------

During a baseflow recession (no rainfall or snowmelt), water drains from
storage at a rate set by the storage–discharge relationship.  If discharge
is a power law of storage depth,

.. math::

    Q = c \, H^b,

then substituting :math:`\mathrm{d}H/\mathrm{d}t = -Q` gives

.. math::

    -\frac{\mathrm{d}Q}{\mathrm{d}t} = K \, Q^n,
    \qquad n = \frac{2b - 1}{b},

where :math:`b` is the HydroRaVENS ``recession_exponent`` and :math:`n`
is the slope on a log–log plot of :math:`-\mathrm{d}Q/\mathrm{d}t`
versus :math:`Q` — the Brutsaert & Nieber (1977) recession plot.

Inverting:

.. math::

    b = \frac{1}{2 - n}

The two exponents have the following special values:

.. list-table::
   :widths: 20 20 50
   :header-rows: 1

   * - B–N slope *n*
     - HydroRaVENS *b*
     - Physical interpretation
   * - 1
     - 1
     - Linear reservoir; exponential recession.
   * - 3/2
     - 2
     - Long-time Boussinesq solution for horizontal, unconfined aquifer
       (Brutsaert & Nieber 1977).
   * - → 2
     - → ∞
     - Extreme nonlinearity; storage empties nearly instantaneously.

The conversion :math:`b = 1/(2-n)` is valid only for :math:`n < 2`.
The Boussinesq short-time solution gives :math:`n = 3`, which is outside
this framework; in practice, early fast-response recessions (surface
runoff, tile drains) produce steep tails on the B–N cloud that should not
be used to set the baseflow recession exponent.

Brutsaert–Nieber Analysis
--------------------------

.. autoclass:: hydroravens.BrutsaertNieber
   :members: fit, to_reservoir_exponent, summary, plot
   :member-order: bysource

Workflow
--------

A typical workflow fits the recession curve and uses the result as a
prior for model calibration:

.. code-block:: python

    import pandas as pd
    from hydroravens import BrutsaertNieber

    df = pd.read_csv('streamflow.csv', parse_dates=['Date'])
    Q  = df['Specific Discharge [mm/day]'].values

    bn = BrutsaertNieber(Q, min_recession_days=3).fit()
    bn.summary()
    bn.plot()

    b_prior = bn.to_reservoir_exponent()
    print(f"Suggested recession_exponent prior: {b_prior:.2f}")

The output might look like::

    Brutsaert–Nieber recession analysis
      Recession pairs used : 412
      Fitted slope  n      : 1.5831
      Fitted coeff  a      : 0.0147
      HydroRaVENS   b      : 2.367
      Reference (long-time Boussinesq): n = 1.5, b = 2.0

The lower envelope of the B–N cloud corresponds to long-duration
recessions (groundwater-dominated) and is the most physically meaningful
region.  The upper scatter reflects short, event-driven recessions (tile
drains, surface runoff) that may follow a different exponent.

Interpretation guide
--------------------

* **n ≈ 1.5, b ≈ 2** — groundwater recession consistent with the
  Boussinesq long-time solution.  A reasonable fixed prior for deep
  (karst or fractured-rock) reservoirs.
* **n ≈ 1.3–1.6, b ≈ 1.5–2.5** — typical mixed catchment range.
* **n > 1.7, b > 3** — rapid nonlinear drainage, often reflecting
  near-surface soil processes.  Calibrate rather than fix.
* **n ≥ 2** — unphysical in this framework; refit using only the lower
  envelope or with a longer ``min_recession_days``.

References
----------

Brutsaert, W. and Nieber, J. L. (1977). Regionalized drought flow
hydrographs from a mature glaciated plateau. *Water Resources Research*,
13(3), 637–643. https://doi.org/10.1029/WR013i003p00637

Kirchner, J. W. (2009). Catchments as simple dynamical systems:
Catchment characterization, rainfall-runoff modeling, and doing hydrology
backward. *Water Resources Research*, 45, W02429.
https://doi.org/10.1029/2008WR006912
