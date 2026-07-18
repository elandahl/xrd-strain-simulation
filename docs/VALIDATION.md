# XRD-calculator validation

Internal, physics-based validation of the GaAs (004) dynamical diffraction
calculator, to be green **before** any external benchmarking (Sci. Rep. Fig. 3
or APS/PLS data). This mirrors the strain repo's `docs/VALIDATION.md`.

Run:

```bash
python scripts/validate_xrd_physics.py           # report + figure, nonzero exit on failure
python -m pytest tests/test_xrd_acceptance.py -q
```

Outputs: `docs/physics_acceptance.json` and `docs/images/xrd_acceptance.png`.
The checks assert on analytic/known quantities, not frozen numerical snapshots
(except the deliberate regression golden), so they document the physics and
catch real regressions.

## 1. Unstrained perfect-crystal rocking curve

An essentially unstrained thick GaAs slab must give the textbook Darwin curve:

- **Peak position** at the kinematic Bragg angle
  `θ_B = asin(λ / 2d_004)` with `d_004 = a/4`. Measured 26.0326° vs 26.0325°
  (+0.18 arcsec dynamical refraction shift). Tolerance: |shift| < 3 arcsec.
- **Peak reflectivity** O(1) for a thick perfect crystal (measured 0.985).
- **Darwin width** a few arcsec (measured ≈ 5.9 arcsec FWHM, consistent with
  GaAs 004 at ~10 keV). Tolerance: 2–15 arcsec.

## 2. Known strain → Bragg shift

A uniformly strained surface layer on the substrate produces a second peak
displaced by the kinematic law `Δθ = -ε·tan(θ_B)` (tensile out-of-plane strain
→ larger d → lower angle):

- ε = 2×10⁻³ gives a layer peak at −196.6 arcsec vs predicted −201.5 arcsec
  (2.4% off, correct sign). Tolerance: 10%.
- the shift **scales linearly with ε**: shift(2e-3)/shift(1e-3) = 2.10
  (expect 2).

This is the core diffraction-physics check: peak separation ↔ strain magnitude.

## 3. Instrument-convolution sanity

Convolution must only *broaden*, never move the peak or invent structure:

- **Broadening** increases FWHM: none 6.0 → aps_7idc(12″) 13.6 → empirical
  18.2 arcsec.
- **No new structure**: the number of local maxima does not increase under
  either kernel.
- **Intensity conservation**: the symmetric normalized `aps_7idc` Gaussian
  conserves integrated intensity (ratio 0.9999). *Note:* the `empirical`
  multi-Gaussian kernel does **not** conserve area (a known normalization
  artifact of the inherited convolution — it rescales the curve; see
  `docs/INSTRUMENTS.md`), so area conservation is only asserted for `aps_7idc`.
- **Peak stability**: `aps_7idc` moves the peak centroid by −0.03 arcsec
  (symmetric kernel). The `empirical` kernel has deliberate offset centers and
  may shift slightly.

![XRD acceptance](images/xrd_acceptance.png)

## 4. Frozen-notebook regression

`tests/test_xrd_acceptance.py::test_frozen_notebook_regression` requires the
current `xrd_slab_gaas` (and the low-memory path) to reproduce a golden rocking
curve **bit-for-bit** (atol 1e-12; observed max abs diff 0.0). The golden data
`tests/data/gaas004_golden.npz` was generated from the archival
thermo-elastic-gaas calculator (tag `paper-v1.0`) on a fixed strain profile,
so the split/refactor is proven to preserve the published numerics without
needing the other repo at test time.

## When to re-run

- After any change to `crystals/gaas_004_dynamical.py`, the detector models, or
  the pipeline.
- Before adding a new substrate/reflection: these checks define the GaAs
  baseline a new calculator should reproduce in structure.
- The checked-in `docs/physics_acceptance.json` is the run-of-record; diff it
  after re-running.

## External next steps

Internal checks do not prove the hard-coded scattering factors and lattice
constant match modern tables. That is Tier 2 of the external plan —
see [`EXTERNAL_BENCHMARKS.md`](EXTERNAL_BENCHMARKS.md) and
[`CONSTANTS_PROVENANCE.md`](CONSTANTS_PROVENANCE.md).
