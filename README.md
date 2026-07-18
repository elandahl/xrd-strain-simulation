# XRD Strain Simulation

Compute X-ray diffraction rocking curves from depth-dependent strain profiles. Split from the published [thermo-elastic-gaas](https://github.com/elandahl/thermo-elastic-gaas) package for extension to new substrates, reflections, and geometries.

## Current calculator

- **`gaas_004_10kev`** — GaAs (004) reflection, 10 keV, dynamical diffraction (low-memory angle-by-angle evaluation, supports fine angular grids)

## Instrument models

Switch with `--instrument` (see `docs/INSTRUMENTS.md`):

- `notebook` (default) — multi-Gaussian detector model from the original notebook
- `aps_7idc` — Gaussian, configurable FWHM (default 1.8 arcsec, per the Sci. Rep. 2022 experiment)
- `none` — raw dynamical-diffraction curve

## Quick start

```bash
# 1. Generate strain profile (strain-wave-simulation repo)
cd ../strain-wave-simulation
python scripts/run.py --no-show

# 2. Compute XRD rocking curve
cd ../xrd-strain-simulation
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
python scripts/run.py \
  --strain-file ../strain-wave-simulation/results/strain_profile.npz \
  --no-show
```

## Extension points

| Future work | Where to add |
|-------------|--------------|
| Si, Ge, InSb substrates | `src/xrd_strain/crystals/` + `register_crystal()` |
| Other reflections / energies | new crystal calculator modules |
| Grazing-incidence or other geometries | extend calculator interface in `crystals/base.py` |
| Different instrument/detector models | `src/xrd_strain/detector/` + a branch in `pipeline.run_xrd` |

## Input format

Reads `.npz` strain profiles from [strain-wave-simulation](https://github.com/elandahl/strain-wave-simulation). See that repo's `docs/STRAIN_PROFILE.md`.

## Paper reproduction

The full combined pipeline for the published paper remains at [thermo-elastic-gaas](https://github.com/elandahl/thermo-elastic-gaas). Chain this repo with `strain-wave-simulation` to reproduce the XRD step.
