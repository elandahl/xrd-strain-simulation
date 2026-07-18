# XRD Strain Simulation

Compute X-ray diffraction rocking curves from depth-dependent strain profiles. Split from the published [thermo-elastic-gaas](https://github.com/elandahl/thermo-elastic-gaas) package for extension to new substrates, reflections, and geometries.

## Current calculator

- **`gaas_004_10kev` (default)** — externally benchmarked GaAs (004),
  10 keV, 300 K, σ-polarized dynamical diffraction. Uses modern lattice and
  scattering factors, species-specific Debye–Waller factors, and distinct
  forward/diffracting structure factors \(F_0,F_h\). Matches Stepanov
  X0h/GID_sl perfect-crystal FWHM and peak reflectivity to ~0.1%.
- **`gaas_004_10kev_legacy`** — notebook-faithful constants and \(F_0=F_h\)
  approximation, retained only for exact archival reproduction.

## Instrument models

Switch with `--instrument` (see `docs/INSTRUMENTS.md` for details and the figure guide):

- `aps_7idc` (default) — single Gaussian at the measured APS 7ID-C resolution
  (default 1.8 arcsec FWHM, configurable). Physical analysis default.
- `empirical` — multi-Gaussian **effective** resolution (instrument + sample +
  fit) from the original notebook; dominant σ ≈ 22 arcsec. Legacy
  paper-figure smoothing option. (Deprecated alias: `notebook`.)
- `none` — raw dynamical-diffraction curve

The `empirical` model reproduces the paper's *plotted* curve better, but it
conflates instrument and sample/fit broadening and is therefore not the
default for new data analysis. See `docs/INSTRUMENTS.md`.

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
Formula mapping and the corrected \(F_0/F_h\) distinction:
[`docs/FORMULA_AUDIT.md`](docs/FORMULA_AUDIT.md).
Tier-4 cross-code benchmark against Stepanov GID_sl on strained-layer
profiles (~1% agreement over 4–5 decades):

```bash
python scripts/benchmark_gid_sl.py
```

See [`docs/GID_SL_BENCHMARK.md`](docs/GID_SL_BENCHMARK.md).
The realistic Figure 3 d'Alembert strain benchmark, including a controlled
comparison using GID_sl's own susceptibilities:

```bash
python scripts/benchmark_fig3_gid_sl.py
```

See
[`docs/FIG3_GID_SL_BENCHMARK.md`](docs/FIG3_GID_SL_BENCHMARK.md).

## Paper reproduction

The full combined pipeline for the published paper remains at [thermo-elastic-gaas](https://github.com/elandahl/thermo-elastic-gaas). Chain this repo with `strain-wave-simulation` to reproduce the XRD step.
