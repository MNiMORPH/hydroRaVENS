Model Description
==================

Theory & Mathematical Formulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

HydroRaVENS is a lumped, daily-timestep conceptual hydrological model. Water 
moves through three sequential processes each day: optional snowpack accumulation/melt, 
routing through cascading linear reservoirs, and evapotranspiration.

All fluxes are expressed as depths over the drainage basin (mm/day).
Mass is conserved to within numerical precision.

Daily Water Balance
~~~~~~~~~~~~~~~~~~~

On each day, the model computes:

.. math::

    P + M - E = \Delta S + Q

where:

* :math:`P` = precipitation (mm/day)
* :math:`M` = snowmelt (mm/day)
* :math:`E` = evapotranspiration (mm/day)
* :math:`\Delta S` = change in storage (mm/day)
* :math:`Q` = streamflow (mm/day)

Snowpack Module (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~

If mean air temperature is provided, snowpack processes are enabled.

**Accumulation:**
  When :math:`T \leq 0°C`, all precipitation accumulates as snow (stored as SWE):
  
  .. math::
  
      \text{SWE}_{t+1} = \text{SWE}_t + P_t

**Melt:**
  When :math:`T > 0°C`, melt is computed using the positive-degree-day (PDD) approach:
  
  .. math::
  
      M_t = \min(\text{SWE}_t, \alpha \cdot T_t \cdot \Delta t)
  
  where :math:`\alpha` is the melt factor (mm SWE per °C per day).
  
  All melt is routed directly to the top reservoir.

**Sublimation (ET Deficit):**
  If :math:`E > P + M`, the ET deficit first sublimates snow:
  
  .. math::
  
      \text{Sublimation} = \min(\text{SWE}_t, E - P - M)
  
  Any remaining deficit is passed to the reservoirs.

Linear Reservoir Cascade
~~~~~~~~~~~~~~~~~~~~~~~~

Water drains through a stack of reservoirs (top = shallowest, bottom = deepest).
Each reservoir drains by exponential decay:

.. math::

    Q_i(t) = H_i(t) \cdot (1 - e^{-\Delta t / \tau_i})

where:

* :math:`H_i` = water depth in reservoir :math:`i` (m)
* :math:`\tau_i` = e-folding residence time (days)
* :math:`\Delta t` = time step (1 day)

**Discharge Partitioning:**
  Of the water drained from reservoir :math:`i`, a fraction :math:`f_i` exits as 
  river discharge; the remainder infiltrates to the next layer:
  
  .. math::
  
      Q_{\text{discharge},i} = f_i \cdot Q_i \\
      Q_{\text{infiltrate},i} = (1 - f_i) \cdot Q_i

**Constraint:**
  The bottom reservoir must fully discharge: :math:`f_{\text{bottom}} = 1.0`

**Storage Update:**
  After drainage, infiltration, and ET:
  
  .. math::
  
      H_i(t+1) = \max(0, H_i(t) - Q_i(t) + Q_{\text{recharge},i} - E_i)

Evapotranspiration
~~~~~~~~~~~~~~~~~~

Two methods are supported:

**Method 1: From Data**
  ET is read directly from the input CSV and scaled to close the annual water balance:
  
  .. math::
  
      E_{\text{scaled}} = E_{\text{observed}} \cdot \frac{P_{\text{annual}} - Q_{\text{observed,annual}}}{ET_{\text{observed,annual}}}

**Method 2: Thornthwaite-Chang (2019)**
  Daily reference ET (:math:`ET_0`) is estimated from temperature and photoperiod:
  
  .. math::
  
      ET_0 = \text{f}(T_{\max}, T_{\min}, \text{photoperiod})
  
  Then scaled by water year to match observed P − Q balance.

In both cases, the annual scaling factor is stored and applied to ensure that 
:math:`P - Q - E \approx 0` on each water year.

Model Skill Evaluation
~~~~~~~~~~~~~~~~~~~~~~

The **Nash-Sutcliffe Efficiency (NSE)** quantifies model performance:

.. math::

    \text{NSE} = 1 - \frac{\sum_t (Q_{\text{mod},t} - Q_{\text{obs},t})^2}
                          {\sum_t (Q_{\text{obs},t} - \bar{Q}_{\text{obs}})^2}

* :math:`\text{NSE} = 1`: Perfect simulation
* :math:`\text{NSE} = 0`: Model performs as well as the observed mean
* :math:`\text{NSE} < 0`: Model worse than using the mean as a predictor

**Typical ranges:**
  * :math:`\text{NSE} > 0.75`: Excellent
  * :math:`0.5 < \text{NSE} < 0.75`: Good
  * :math:`0.3 < \text{NSE} < 0.5`: Satisfactory
  * :math:`\text{NSE} < 0.3`: Poor

Model Assumptions & Limitations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Strengths:**

* ✅ Mass-conserving (exact balance)
* ✅ Minimal data requirements (daily P and Q)
* ✅ Fast computation (runs decades in seconds)
* ✅ Physically interpretable parameters
* ✅ Suitable for ungauged basins (transfer parameters)

**Limitations:**

* ⚠️ Lumped model (no spatial variability)
* ⚠️ Linear reservoirs (may not capture threshold behavior)
* ⚠️ Daily timestep only (not suitable for event-scale analysis)
* ⚠️ Simplified groundwater–surface-water interaction
* ⚠️ Snowpack simplified (PDD model; ignores energy balance)
* ⚠️ No representation of lakes or artificial storage
* ⚠️ ET method choice matters; Thornthwaite-Chang is approximate

**Best suited for:**

* Long-term water balance studies (years–decades)
* Climate impact assessments
* Ungauged or poorly-gauged basins
* Parameter transfer to similar basins
* Educational modeling
