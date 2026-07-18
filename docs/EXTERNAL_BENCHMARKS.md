# External benchmarking plan

Four-tier plan for validating the XRD (and eventually strain) codes against
sources *outside* this repository. Internal physics acceptance
([VALIDATION.md](VALIDATION.md)) must stay green; these tiers are the next
layer before fitting APS/PLS data.

| Tier | Goal | Status |
|------|------|--------|
| **1. Formulas** | Look up textbook/paper formulas; confirm they match the code and the Sci. Rep. paper | **GaAs (004) audit complete** — corrected \(F_0=F_h\) approximation; see [FORMULA_AUDIT.md](FORMULA_AUDIT.md) |
| **2. Constants** | Provenance table for every hard-coded constant at the correct energy, reflection, material, and temperature | **In progress** — see [CONSTANTS_PROVENANCE.md](CONSTANTS_PROVENANCE.md) |
| **3. Internal checks vs external refs** | Re-score perfect-crystal Darwin curve and strained-layer Bragg shift against literature formulas and Stepanov X0h using *audited* constants | **Complete for GaAs (004), 10 keV, 300 K, σ** — FWHM and peak R match X0h/GID_sl to ~0.1%; see [CONSTANTS_SENSITIVITY.md](CONSTANTS_SENSITIVITY.md) |
| **4. Cross-code / arbitrary strain** | Compare against other dynamical XRD methods (e.g. Stepanov GID_sl / X0h, `xrayutilities`) on synthetic and paper strain fields; eventually APS/PLS data | Planned (later) |

## Suggested order of work

1. Constants provenance table (Tier 2) ← **current**
2. Formula ↔ code ↔ paper notes (Tier 1)
3. Score existing XRD acceptance checks against literature (Tier 3)
4. One online-solver comparison on a uniform strained layer (Tier 4 smoke test)
5. Paper Fig. 3 strain → XRD vs published curve
6. APS / PLS experimental data

## Supporting practices

- **Paper SI audit** — every parameter the paper quotes vs what the code uses
  (partially done for strain in `thermo-elastic-gaas/docs/PAPER_VALIDATION.md`).
- **Sensitivity table** — how much Darwin width / Bragg position / layer shift
  move when \(a\), \(f'\), \(f''\), or \(E\) change by their literature uncertainty.
- Keep instrument convolution **off** (`none`) for cross-code comparisons;
  add `aps_7idc` / `empirical` only when matching published plots or data.

## Scope note

Tier 2 below focuses on the **GaAs (004) 10 keV XRD calculator**. Strain-side
constants (sound speeds, densities, expansion coefficients, G, …) have a
parallel audit path in `strain-wave-simulation`; cross-links will be added as
that table is written.
