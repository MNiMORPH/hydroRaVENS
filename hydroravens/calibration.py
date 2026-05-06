"""
hydroravens.calibration
~~~~~~~~~~~~~~~~~~~~~~~
Run hydroRaVENS with a given parameter set and return a goodness-of-fit
metric, optionally scored over a specific date window.

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
"""

import numpy as np
import pandas as pd

from .hydroravens import Buckets


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
    eps = 0.01 * o.mean()   # 1 % of mean observed flow; keeps log finite
    return _kge(np.log(m + eps), np.log(o + eps))


_METRICS = {
    'NSE':    _nse,
    'KGE':    _kge,
    'logKGE': _log_kge,
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_and_score(cfg, t_efold=None, f_to_discharge=None, Hmax=None,
                  melt_factor=None, start=None, end=None, spin_up_cycles=3,
                  metric='KGE'):
    """
    Run hydroRaVENS and return a goodness-of-fit score.

    Parameters
    ----------
    cfg : str
        Path to a YAML configuration file. Should have spin_up_cycles: 0
        so that spin-up is performed here with the calibrated parameters.
    t_efold : list of float, optional
        E-folding residence times [days], one entry per reservoir in order
        from shallowest to deepest. Overrides the values in cfg.
    f_to_discharge : list of float, optional
        Fraction of each reservoir's exfiltration that goes directly to
        streamflow (vs. percolating to the next reservoir). Provide one
        value per reservoir except the deepest, which must always be 1.0.
        Overrides the values in cfg.
    Hmax : list of float, optional
        Maximum effective water depths [m], one per reservoir. Overrides
        the values in cfg.
    melt_factor : float, optional
        Degree-day snowmelt factor [mm SWE per degC per day]. Overrides
        the value in cfg.
    start : str or datetime-like, optional
        Start of the scoring window (inclusive). The metric is computed
        only over dates >= start. Spin-up still runs the full record.
    end : str or datetime-like, optional
        End of the scoring window (inclusive). The metric is computed only
        over dates <= end. Spin-up still runs the full record.
    spin_up_cycles : int, optional
        Number of times to loop the full record before the scored run.
        Default is 3.
    metric : {'KGE', 'NSE', 'logKGE'}, optional
        Goodness-of-fit metric to return.  Default is 'KGE'.

    Returns
    -------
    score : float
        Goodness-of-fit (higher is better for all metrics).
        Returns np.nan if the scoring window contains no valid data.
    b : Buckets
        The Buckets object after the final run, available for plotting
        or further diagnostics.

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

    # --- Parameter overrides ---
    if t_efold is not None:
        for i, val in enumerate(t_efold):
            b.reservoirs[i].t_efold = val

    if f_to_discharge is not None:
        for i, val in enumerate(f_to_discharge):
            b.reservoirs[i].f_to_discharge = val

    if Hmax is not None:
        for i, val in enumerate(Hmax):
            b.reservoirs[i].Hmax = val

    if melt_factor is not None and b.has_snowpack:
        b.snowpack.melt_factor = melt_factor

    # --- Spin up on the full record with calibrated parameters ---
    for _ in range(spin_up_cycles):
        b.run()
        b._timestep_i = b.hydrodata.index[0]

    # --- Final scored run ---
    b.run()

    # --- Mask to scoring window ---
    q_mod = b.hydrodata['Specific Discharge (modeled) [mm/day]']
    q_obs = b.hydrodata['Specific Discharge [mm/day]']

    mask = q_mod.notna() & q_obs.notna()
    if start is not None:
        mask &= b.hydrodata['Date'] >= pd.Timestamp(start)
    if end is not None:
        mask &= b.hydrodata['Date'] <= pd.Timestamp(end)

    if mask.sum() == 0:
        return np.nan, b

    m = np.asarray(q_mod[mask], dtype=float)
    o = np.asarray(q_obs[mask], dtype=float)

    score = _METRICS[metric](m, o)
    return score, b
