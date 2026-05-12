# Cannon River — forward run

Runs hydroRaVENS for the Cannon River catchment (Minnesota, USA) with a
fixed parameter set and prints goodness-of-fit diagnostics.

**Catchment:** Cannon River near Red Wing, MN (USGS 05355200; 3800 km²)
**Period:** 1992–1995 (daily)

## Requirements

- `hydroravens` installed (`pip install hydroravens`)

## Usage

```bash
python run_forward.py
```

Prints KGE\_logKGE\_logFDC, KGE\_logFDC, AIC, and BFI, then saves
`forward_run.png`.

## Files

| File | Description |
|------|-------------|
| `cannon_cfg.yml` | Model configuration (reservoirs, snowmelt, modules) |
| `CannonTestInput.csv` | Daily forcing and observed discharge |
| `run_forward.py` | Run script |

## Parameter estimation

To calibrate the parameters in `cannon_cfg.yml` see the companion example
in `../cannon_inverse/`.
