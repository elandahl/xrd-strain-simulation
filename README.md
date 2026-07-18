# XRD Strain Simulation

Compute X-ray diffraction rocking curves from depth-dependent strain profiles. Split from the published [thermo-elastic-gaas](https://github.com/elandahl/thermo-elastic-gaas) package for extension to new substrates, reflections, and geometries.

## Current calculator

- **`gaas_004_10kev`** — GaAs (004) reflection, 10 keV, dynamical diffraction (low-memory angle-by-angle evaluation, supports fine angular grids)

## Instrument models

Switch with `--instrument` (see `docs/INSTRUMENTS.md` for details and the figure guide):

- `empirical` (default) — multi-Gaussian **effective** resolution (instrument + sample + fit) from the original notebook; dominant σ ≈ 22 arcsec. Heavy smoothing; best matches the smoothness of the published Fig. 3 fit curve. (Deprecated alias: `notebook`.)
- `aps_7idc` — single Gaussian at the paper's **stated** APS 7ID-C resolution (default 1.8 arcsec FWHM, configurable). Light blur only.
- `none` — raw dynamical-diffraction curve

Note: the `empirical` model reproduces the paper's *plotted* curve better than the paper's own stated 1.8 arcsec resolution does — the published curve is broader than pure instrument response (sample mosaic, depth averaging, beam divergence, fit smoothing). See `docs/INSTRUMENTS.md`.

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

## Validation

**Internal** physics acceptance (perfect-crystal Darwin curve, strained-layer
Bragg shift, instrument-convolution sanity, frozen-notebook regression):

```bash
python scripts/validate_xrd_physics.py
python -m pytest tests/test_xrd_acceptance.py -q
```

See [`docs/VALIDATION.md`](docs/VALIDATION.md).

**External** benchmarking (formulas → constants → literature-scored checks →
cross-code / data) is planned in
[`docs/EXTERNAL_BENCHMARKS.md`](docs/EXTERNAL_BENCHMARKS.md). Tier 2 constants
provenance for GaAs (004) @ 10 keV:

```bash
python scripts/audit_constants.py
```

See [`docs/CONSTANTS_PROVENANCE.md`](docs/CONSTANTS_PROVENANCE.md).
Tier-3 one-at-a-time sensitivity and Stepanov X0h comparison:
[`docs/CONSTANTS_SENSITIVITY.md`](docs/CONSTANTS_SENSITIVITY.md).

## Paper reproduction

The full combined pipeline for the published paper remains at [thermo-elastic-gaas](https://github.com/elandahl/thermo-elastic-gaas). Chain this repo with `strain-wave-simulation` to reproduce the XRD step.
