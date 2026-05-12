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

Optional Process Modules
~~~~~~~~~~~~~~~~~~~~~~~~~

Several processes can be individually enabled or disabled through the
``modules`` block in the configuration file (see :doc:`configuration`).
All modules that require temperature data are silently inactive when
``Mean Temperature [C]`` is absent from the input CSV, regardless of
their flag setting.

.. list-table::
   :widths: 20 10 60
   :header-rows: 1

   * - Module
     - Default
     - Purpose
   * - ``snowpack``
     - on
     - Accumulation, PDD melt, sublimation, and rain-on-snow.
   * - ``frozen_ground``
     - on
     - Frozen ground index (FGI) that blocks deep infiltration.
       Requires ``snowpack: true``.
   * - ``rain_on_snow``
     - on
     - Sensible-heat contribution of warm rain on snowpack.
       Requires ``snowpack: true``.
   * - ``direct_runoff``
     - off
     - Hortonian-inspired bypass fraction. See below.

Snowpack Module (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~

If mean air temperature is provided, snowpack processes are enabled.

**Accumulation:**
  When :math:`T \leq 0°C`, net water input (precipitation minus ET) accumulates
  as snow (stored as SWE):

  .. math::

      \text{SWE}_{t+1} = \text{SWE}_t + (P_t - E_t)

**Melt:**
  When :math:`T > 0°C`, melt is computed using the positive-degree-day (PDD) approach:

  .. math::

      M_t = \min(\text{SWE}_t,\ \alpha \cdot T_t \cdot \Delta t + M_{\text{ROS},t})

  where :math:`\alpha` is the melt factor (mm SWE °C⁻¹ day⁻¹) and
  :math:`M_{\text{ROS}}` is the rain-on-snow sensible-heat contribution
  (see below). All melt is routed directly to the top reservoir.

**Rain-on-snow (ROS) sensible heat:**
  Rain arriving at temperature :math:`T > 0°C` carries thermal energy that
  can melt additional snow:

  .. math::

      M_{\text{ROS},t} = \frac{c_p}{L_f} \cdot T_t \cdot P_t

  where :math:`c_p / L_f \approx 0.01253\ °C^{-1}` is the ratio of the
  specific heat of water to the latent heat of fusion. This term is
  typically small relative to the PDD term except during warm rain-on-snow
  events (McCabe et al. 2007; Würzer et al. 2016).

**ET deficit:**
  When precipitation minus ET is negative, the deficit first sublimates
  snow:

  .. math::

      \text{Sublimation} = \min\!\left(\text{SWE}_t,\ \max(0,\ E_t - P_t)\right)

  Any deficit exceeding available SWE is passed to the top subsurface
  reservoir and, if still unmet, carried forward to the next time step.

Frozen Ground Module (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When a ``fdd_threshold`` is set, the model tracks a **frozen ground index**
(FGI; Molnau & Bissell 1983) that accumulates freezing degree-days,
decays passively each day, and additionally decays during warm periods:

.. math::

    \text{FGI}_t = \max\!\left(0,\ A \cdot \text{FGI}_{t-1} - T_{\text{eff},t} - D_t\right)

where :math:`A` is the daily decay coefficient (``fgi_decay_coeff``,
default 0.97 following Molnau & Bissell 1983 and all major
implementations), :math:`T_{\text{eff},t} = T_t \cdot e^{-k \cdot
\text{SWE}_t}` is the snow-insulation-adjusted temperature (°C; negative
values increase FGI, positive values reduce it), and :math:`D_t` is an
additional thaw credit described below. The passive decay :math:`(1 - A)
\approx 3\%` per day prevents indefinite accumulation during long cold
spells and sets a finite steady-state
:math:`\text{FGI}^* = |T| / (1 - A)` for sustained temperature
:math:`T < 0`. When :math:`\text{FGI}_t` exceeds ``fdd_threshold``, the
top reservoir's exfiltration fraction is set to 1.0 so that all drainage
becomes direct runoff, simulating frozen-soil blockage of deep
infiltration.

**Coupling snowmelt and frozen-ground thaw via the melt factor:**
  The PDD melt factor :math:`\alpha` has units of mm SWE per °C·day,
  making it a natural conversion factor between the thermal forcing
  (°C·day) and a water-equivalent depth (mm SWE):

  .. math::

      \text{mm SWE} = \alpha \cdot \text{°C·day}
      \qquad \Longleftrightarrow \qquad
      \text{°C·day} = \frac{\text{mm SWE}}{\alpha}

  When total melt energy in a timestep exceeds the available SWE, the
  leftover energy :math:`\Delta E` is expressed in mm SWE. Dividing by
  :math:`\alpha` converts it back to degree-days — the same currency the
  FGI uses — so the residual can be credited toward thawing frozen ground:

  .. math::

      D_t = \frac{\max(0,\ \alpha T_t \Delta t + M_{\text{ROS},t} - \text{SWE}_t)}{\alpha}

  This means that once the snowpack is fully depleted, any remaining
  thermal energy continues to thaw frozen soil rather than being
  discarded. The melt factor thus serves as the bridge between the two
  empirical degree-day representations.

  Note that :math:`\alpha` characterises the snow surface, not the soil.
  Applying the same conversion to the soil implicitly assumes that the
  atmosphere-to-surface thermal coupling is the same in both cases, which
  is a simplification. Within a degree-day framework, however, the
  approach is internally self-consistent.

**Snow insulation:**
  A deep snowpack buffers the soil surface from cold air, reducing the
  effective freezing degree-days that accumulate in the FGI. This is
  represented by an exponential insulation factor:

  .. math::

      T_{\text{eff},t} = T_t \cdot e^{-k \cdot \text{SWE}_t}

  where :math:`k` (mm⁻¹ SWE) is the snow insulation decay constant
  (``snow_insulation_k`` in the ``snowmelt`` config section; default
  0.0). The factor is applied to both freezing and thawing temperature
  forcing; meltwater heat delivery (``excess_dd``) is not scaled because
  meltwater reaches the soil surface directly. The parameterisation
  originates in Molnau & Bissell (1983), and was adopted by LISFLOOD
  (van der Knijff et al. 2010) and GSSHA (Downer & Ogden 2004).

  .. note::

      ``snow_insulation_k`` and ``fdd_threshold`` are correlated: both
      control how much frozen-ground effect the model sees. Calibrating
      them simultaneously from streamflow alone leads to equifinality —
      the optimizer trades a near-zero threshold against moderate
      insulation rather than finding physically meaningful values for
      either. Recommended practice:

      * Fix ``snow_insulation_k`` at a literature or field-derived value
        and calibrate only ``fdd_threshold``, or
      * Leave ``snow_insulation_k = 0`` (default) and treat
        ``fdd_threshold`` as the sole free FGI parameter.

      The insulation term is most useful when independent observations
      (soil temperature, frost depth) are available to constrain
      :math:`k`, or for deep-snowpack alpine catchments where the
      insulation effect is large relative to the threshold uncertainty.

Direct Runoff Bypass (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``direct_runoff: true`` in the ``modules`` block, a fixed fraction
:math:`\gamma` of positive daily recharge is intercepted before it enters
the top reservoir and routed directly to the stream:

.. math::

    Q_{\text{direct},t} = \gamma \cdot \max(0,\ R_t)

    R_{\text{remaining},t} = (1 - \gamma) \cdot R_t

where :math:`R_t` is the recharge available after ET (mm/day) and
:math:`\gamma` is ``direct_runoff_fraction`` in the ``general``
configuration section.

This formulation is *Hortonian-inspired* — conceptually motivated by
infiltration-excess overland flow — but cannot be a rigorous
physical representation at the daily timestep. Without sub-daily
intensity data, it is impossible to determine whether the daily
precipitation total actually exceeds hydraulic conductivity at any
moment. The bypass may nonetheless provide a useful empirical
correction for catchments with significant impervious area, compacted
soils, or urban fractions, and for days dominated by a single intense
storm event (Yilmaz et al. 2008).

The module is off by default (``direct_runoff: false``). Unless the
calibration objective clearly demands it, leaving it off is recommended
to avoid equifinality: the bypass can mimic effects already represented
by the shallow-reservoir timescale, causing parameters to lose their
intended physical interpretation.

Linear Reservoir Cascade
~~~~~~~~~~~~~~~~~~~~~~~~

Water drains through a stack of reservoirs (top = shallowest, bottom = deepest).
Each reservoir first receives its recharge input, then drains by exponential decay:

.. math::

    Q_i(t) = \bigl(H_i(t) + Q_{\text{recharge},i}\bigr) \cdot (1 - e^{-\Delta t / \tau_i})

where:

* :math:`H_i` = water depth in reservoir :math:`i` at the start of the time step (mm)
* :math:`Q_{\text{recharge},i}` = water input to reservoir :math:`i` this time step (mm)
* :math:`\tau_i` = e-folding residence time (days)
* :math:`\Delta t` = time step (1 day)

**Discharge Partitioning:**
  Of the water drained from reservoir :math:`i`, a fraction :math:`f_i` exits as 
  river discharge; the remainder infiltrates to the next layer:
  
  .. math::
  
      Q_{\text{discharge},i} = f_i \cdot Q_i \\
      Q_{\text{infiltrate},i} = (1 - f_i) \cdot Q_i

**Constraint:**
  The bottom reservoir should fully discharge (:math:`f_{\text{bottom}} = 1.0`).
  A warning is issued if not, as this violates mass conservation.

**Storage Update:**
  Recharge is applied first, then exponential drainage:

  .. math::

      H_i(t+1) = \bigl(H_i(t) + Q_{\text{recharge},i}\bigr)\,e^{-\Delta t/\tau_i}

  Overflow above :math:`H_{\max}` exits immediately as direct runoff; any
  deficit is passed to the next-deeper reservoir. ET is not subtracted
  separately at the reservoir level — it is already incorporated into the
  recharge input to the top reservoir.

**Physical interpretation:**
  No reservoir is fixed to a particular process; meaning is set by
  parameter choice. Successive reservoirs naturally span progressively
  longer timescales — interflow (days), soil moisture (months),
  groundwater (years) — but that mapping is the user's choice, analogous
  to the multi-component runoff structure of HBV (Bergström 1976).

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
:math:`P - Q - E = 0` over each water year.

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
* ⚠️ Daily timestep by design — degree-day snowmelt, Thornthwaite ET, and
  linear reservoir drainage are daily-scale parameterisations that lose
  physical meaning at finer resolution
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
