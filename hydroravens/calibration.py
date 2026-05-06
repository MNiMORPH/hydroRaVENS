"""
hydroravens.calibration
~~~~~~~~~~~~~~~~~~~~~~~
Run hydroRaVENS with a given parameter set and return a CalibResult
named tuple containing the goodness-of-fit score, AIC, baseflow index,
flow duration curve, end-of-run reservoir states, and the Buckets object.

Intended use: call run_and_score() from a Dakota driver or any other
optimizer. The YAML config passed to cfg must have spin_up_cycles: 0;
spin-up is managed here so calibrated parameters are in effect before
any simulation runs.

Supported metrics
-----------------
'NSE'    Nash-Sutcliffe Efficiency.  Biased toward high flows because its
         denominator is the variance of observed discharge (dominated by
         peaks).  Use as a baseline or when peaks are the primary concern.

'KGE'    Kling-Gupta Efficiency.  Decomposes fit into correlation (r),
         variability ratio (alpha = std_mod/std_obs), and bias ratio
         (beta = mean_mod/mean_obs), weighting all three equally.
         Better balanced across the full flow range than NSE.

'logKGE' KGE applied to log-transformed flows.  Shifts sensitivity toward
         low flows and base flow; useful when base flow matters as much as
         peaks.  A small epsilon (1 % of mean observed flow) is added
         before taking logs to avoid log(0).

AIC
---
AIC = N * ln(SS_res_log / N) + 2k, where SS_res_log is the sum of squared
residuals on log-transformed flows and k is the number of free parameters
passed to run_and_score().  Log-transforming flows makes the Gaussian
residual assumption more defensible for discharge data.  AIC is intended
for comparing models with different numbers of reservoirs; lower is better.

Baseflow index (BFI)
--------------------
Computed with the Eckhardt (2005) recursive digital filter applied to both
observed and modelled specific discharge within the scoring window.
BFI = baseflow volume / total flow volume.  alpha=0.98 and bfi_max=0.80
are standard values for perennial daily streamflow.

Flow duration curve (FDC)
--------------------------
Discharge at 99 evenly-spaced exceedance probabilities (0.5–99.5 %).
Stored as pd.Series indexed by exceedance probability (%) in both
CalibResult.fdc_obs and CalibResult.fdc_mod.

Chaining decades
----------------
run_and_score() returns final_states, a dict of reservoir water depths and
snowpack SWE at the end of the scored run.  Pass this as initial_states to
the next decade's run_and_score() call with spin_up_cycles=0 so that water
storage is physically continuous across decade boundaries.
"""

from collections import namedtuple

import numpy as np
import pandas as pd

from .hydroravens import Buckets


# ---------------------------------------------------------------------------
# Named tuple for return value
# ---------------------------------------------------------------------------

CalibResult = namedtuple('CalibResult', [
    'score',        # float: KGE / NSE / logKGE (higher is better)
    'aic',          # float: AIC on log-transformed flows (lower is better)
    'bfi_obs',      # float: observed baseflow index
    'bfi_mod',      # float: modelled baseflow index
    'fdc_obs',      # pd.Series: observed flow at exceedance probabilities
    'fdc_mod',      # pd.Series: modelled flow at exceedance probabilities
    'final_states', # dict: {'reservoirs': [...], 'snowpack': float}
    'buckets',      # Buckets object after the final run
])


# ---------------------------------------------------------------------------
# Metric helpers – operate on plain numpy arrays
# ---------------------------------------------------------------------------

def _nse(m, o):
    return float(1.0 - np.sum((m - o) ** 2) / np.sum((o - o.mean()) ** 2))


def _kge(m, o):
    r     = np.corrcoef(m, o)[0, 1]
    alpha = m.std() / o.std()
    beta  = m.mean() / o.mean()
    return float(1.0 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2))


def _log_kge(m, o):
    eps = 0.01 * o.mean()
    return _kge(np.log(m + eps), np.log(o + eps))


def _aic(m, o, k):
    eps        = 0.01 * o.mean()
    ss_res_log = np.sum((np.log(m + eps) - np.log(o + eps)) ** 2)
    n          = len(o)
    return float(n * np.log(ss_res_log / n) + 2 * k)


def _eckhardt_bfi(q, alpha=0.98, bfi_max=0.80):
    """
    Eckhardt (2005) recursive digital filter for baseflow separation.
    Returns baseflow index = baseflow volume / total flow volume.
    alpha   : recession constant (~0.98 for daily perennial streams)
    bfi_max : maximum BFI (0.80 perennial; 0.50 ephemeral)
    """
    b     = np.empty_like(q, dtype=float)
    b[0]  = q[0] * bfi_max
    denom = 1.0 - alpha * bfi_max
    for t in range(1, len(q)):
        b[t] = ((1.0 - bfi_max) * alpha * b[t - 1]
                + (1.0 - alpha) * bfi_max * q[t]) / denom
        if b[t] > q[t]:
            b[t] = q[t]
    return float(b.sum() / q.sum())


_FDC_PROBS = np.arange(0.5, 100.0, 1.0)   # exceedance probabilities (%)


def _fdc(q):
    """pd.Series of discharge at standard exceedance probabilities."""
    flows = np.percentile(q, 100.0 - _FDC_PROBS)
    return pd.Series(flows, index=_FDC_PROBS, name='Specific discharge [mm/day]')


def _nash_cascade(q, N, K, dt=1.0):
    """
    Route a runoff time series through a Nash cascade of N identical
    linear reservoirs, each with storage time constant K [days].

    The cascade impulse response (instantaneous unit hydrograph) is a
    two-parameter gamma distribution:

        h(t) = t^(N-1) * exp(-t/K) / (K^N * Gamma(N))

    with mean travel time N*K [days] and variance N*K^2 [days^2].  For
    N = 1 the response reduces to a single exponential.

    Each reservoir is updated with the exact analytical solution for
    piecewise-constant inflow over a timestep dt:

        S_i(t+dt) = S_i(t) * exp(-dt/K)  +  K * Q_{i-1}(t) * (1 - exp(-dt/K))
        Q_i(t+dt) = S_i(t+dt) / K

    This form is unconditionally stable for any K > 0 and dt > 0, unlike
    the forward-Euler discretisation which requires dt/K < 2.

    Parameters
    ----------
    q : array-like, shape (T,)
        Runoff input time series [mm/day].
    N : int
        Number of linear reservoirs in series (shape parameter).
        N = 2 is typical for medium-sized catchments (area ~ 10^3 km^2).
    K : float
        Storage time constant of each reservoir [days] (scale parameter).
        Mean travel time through the cascade is N * K.
    dt : float, optional
        Timestep [days].  Default 1.0.

    Returns
    -------
    np.ndarray, shape (T,)
        Routed discharge time series [mm/day].

    References
    ----------
    Nash, J. E. (1957). The form of the instantaneous unit hydrograph.
        IAHS Publ. 45, 114–121.
        (Introduced the N-reservoir cascade and its gamma IUH.)

    Dooge, J. C. I. (1959). A general theory of the unit hydrograph.
        J. Geophys. Res., 64(2), 241–256.
        https://doi.org/10.1029/JZ064i002p00241

    Rodriguez-Iturbe, I. and Valdés, J. B. (1979). The geomorphologic
        structure of hydrologic response. Water Resour. Res., 15(6),
        1409–1420.
        https://doi.org/10.1029/WR015i006p01409
        (Shows that the gamma IUH arises naturally from Horton scaling
        laws, justifying its use without explicit flow-path geometry.)
    """
    q     = np.asarray(q, dtype=float)
    N     = int(round(N))
    alpha = np.exp(-dt / K)          # exact decay factor over one timestep
    beta  = K * (1.0 - alpha)        # input gain:  K*(1 - exp(-dt/K))

    S   = np.zeros(N)                # initial storage in each reservoir [mm]
    out = np.empty_like(q)

    for t in range(len(q)):
        inflow = q[t]
        for i in range(N):
            S[i]   = alpha * S[i] + beta * inflow
            inflow = S[i] / K        # outflow from reservoir i → inflow to i+1
        out[t] = inflow              # outflow from the final reservoir

    return out


_METRICS = {'NSE': _nse, 'KGE': _kge, 'logKGE': _log_kge}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_and_score(cfg, t_efold=None, f_to_discharge=None, Hmax=None,
                  melt_factor=None, fdd_threshold=None, initial_states=None,
                  start=None, end=None, spin_up_cycles=3,
                  metric='KGE', routing_N=2, routing_K=None):
    """
    Run hydroRaVENS and return a CalibResult named tuple.

    Parameters
    ----------
    cfg : str
        Path to a YAML configuration file. Should have spin_up_cycles: 0
        so that spin-up is performed here with the calibrated parameters.
    t_efold : list of float, optional
        E-folding residence times [days], one per reservoir (shallowest
        first). Overrides the values in cfg.
    f_to_discharge : list of float, optional
        Exfiltration fractions to stream, one per reservoir except the
        deepest (which is always 1.0). Overrides the values in cfg.
    Hmax : list of float, optional
        Maximum effective water depths [m], one per reservoir. Overrides
        the values in cfg.
    melt_factor : float, optional
        Degree-day snowmelt factor [mm SWE per degC per day]. Overrides
        the value in cfg.
    fdd_threshold : float or None, optional
        Frozen ground index threshold [°C·day].  The frozen ground index
        accumulates freezing degree-days and decays during warming
        (Molnau & Bissell 1983, https://westernsnowconference.org/sites/
        westernsnowconference.org/PDFs/1983Molnau.pdf).  When the index
        exceeds fdd_threshold, infiltration from the top reservoir to
        deeper layers is set to zero for that timestep (all drainage
        becomes direct runoff).  None (default) disables the effect.
    initial_states : dict, optional
        Starting reservoir water depths and snowpack SWE, as returned by
        a previous call's CalibResult.final_states.  Format::

            {'reservoirs': [H_shallow, H_deep, ...],
             'snowpack':    H_snow_SWE}

        When provided, these override the H0 values from cfg.  Use with
        spin_up_cycles=0 when chaining consecutive decades so that water
        storage is physically continuous.
    start : str or datetime-like, optional
        Start of the scoring window (inclusive). Score, AIC, BFI, and FDC
        are all computed within this window. Spin-up still uses the full
        record.
    end : str or datetime-like, optional
        End of the scoring window (inclusive). Same as start.
    spin_up_cycles : int, optional
        Number of times to loop the full record before the scored run.
        Default is 3.  Set to 0 when providing initial_states for chained
        decade runs.
    metric : {'KGE', 'NSE', 'logKGE'}, optional
        Goodness-of-fit metric.  Default is 'KGE'.
    routing_N : int, optional
        Number of identical linear reservoirs in the Nash cascade used
        for channel routing (shape parameter of the gamma IUH).
        Default is 2.  Set routing_K to enable routing; routing_N is
        counted as a free parameter in k only when it is explicitly
        calibrated (the caller must add 1 to the k count via
        run_and_score if routing_N is varied).
    routing_K : float or None, optional
        Storage time constant [days] of each Nash-cascade reservoir
        (scale parameter).  Mean travel time through the cascade is
        routing_N * routing_K.  When None (default), no routing is
        applied and the model output is compared directly to observed
        discharge.  When provided, routing_K is counted as one free
        parameter.

    Returns
    -------
    CalibResult
        Named tuple with fields:
        score        : goodness-of-fit (higher is better)
        aic          : AIC on log flows (lower is better; compare across
                       models with different reservoir counts)
        bfi_obs      : observed baseflow index (Eckhardt filter)
        bfi_mod      : modelled baseflow index (Eckhardt filter)
        fdc_obs      : pd.Series, observed FDC indexed by exceedance %
        fdc_mod      : pd.Series, modelled FDC indexed by exceedance %
        final_states : dict suitable for use as initial_states next decade
        buckets      : Buckets object after the final run

        All scalar fields are np.nan if the scoring window is empty.

    Notes
    -----
    Spin-up runs the full record regardless of start/end so that the
    ET multiplier (computed per water year during initialize()) remains
    valid and the initial storage state reflects long-run climatology.
    """
    if metric not in _METRICS:
        raise ValueError(f"metric must be one of {list(_METRICS)}; got {metric!r}")

    b = Buckets()
    b.initialize(cfg)

    # --- Parameter overrides and free-parameter count ---
    k = 0

    if t_efold is not None:
        for i, val in enumerate(t_efold):
            b.reservoirs[i].t_efold = val
        k += len(t_efold)

    if f_to_discharge is not None:
        for i, val in enumerate(f_to_discharge):
            b.reservoirs[i].f_to_discharge = val
        k += len(f_to_discharge)

    if Hmax is not None:
        for i, val in enumerate(Hmax):
            b.reservoirs[i].Hmax = val
        k += len(Hmax)

    if melt_factor is not None and b.has_snowpack:
        b.snowpack.melt_factor = melt_factor
        k += 1

    if fdd_threshold is not None:
        b.fdd_threshold = fdd_threshold
        k += 1

    if routing_K is not None:
        k += 1

    # --- Optional: override initial storage states (for chained decades) ---
    if initial_states is not None:
        for i, h in enumerate(initial_states['reservoirs']):
            b.reservoirs[i].Hwater = h
        if b.has_snowpack:
            b.snowpack.Hwater = initial_states.get('snowpack', 0.0)

    # --- Spin up on the full record with calibrated parameters ---
    for _ in range(spin_up_cycles):
        b.run()
        b._timestep_i = b.hydrodata.index[0]

    # --- Final scored run ---
    b.run()

    # --- Capture end-of-run states for chaining to next decade ---
    final_states = {
        'reservoirs': [res.Hwater for res in b.reservoirs],
        'snowpack':   b.snowpack.Hwater if b.has_snowpack else 0.0,
    }

    # --- Optional: route total runoff through Nash cascade ---
    # Routing is applied to the full simulation output before the scoring
    # window is applied, so that routing-reservoir state is correct at the
    # window boundaries.  The routed series is written back into the Buckets
    # hydrodata frame so that CalibResult.buckets reflects routed discharge
    # for downstream plotting.
    q_mod = b.hydrodata['Specific Discharge (modeled) [mm/day]']
    if routing_K is not None:
        routed = _nash_cascade(q_mod.to_numpy(), routing_N, routing_K)
        q_mod  = pd.Series(routed, index=q_mod.index, name=q_mod.name)
        b.hydrodata['Specific Discharge (modeled) [mm/day]'] = q_mod

    # --- Mask to scoring window ---
    q_obs = b.hydrodata['Specific Discharge [mm/day]']

    mask = q_mod.notna() & q_obs.notna()
    if start is not None:
        mask &= b.hydrodata['Date'] >= pd.Timestamp(start)
    if end is not None:
        mask &= b.hydrodata['Date'] <= pd.Timestamp(end)

    nan_result = CalibResult(
        score=np.nan, aic=np.nan,
        bfi_obs=np.nan, bfi_mod=np.nan,
        fdc_obs=pd.Series(dtype=float), fdc_mod=pd.Series(dtype=float),
        final_states=final_states, buckets=b,
    )
    if mask.sum() == 0:
        return nan_result

    m = np.asarray(q_mod[mask], dtype=float)
    o = np.asarray(q_obs[mask], dtype=float)

    return CalibResult(
        score       = _METRICS[metric](m, o),
        aic         = _aic(m, o, k),
        bfi_obs     = _eckhardt_bfi(o),
        bfi_mod     = _eckhardt_bfi(m),
        fdc_obs     = _fdc(o),
        fdc_mod     = _fdc(m),
        final_states= final_states,
        buckets     = b,
    )
