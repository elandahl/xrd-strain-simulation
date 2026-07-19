# External benchmarking plan

Four-tier plan for validating the XRD (and eventually strain) codes against
sources *outside* this repository. Internal physics acceptance
([VALIDATION.md](VALIDATION.md)) must stay green; these tiers are the next
layer before fitting APS/PLS data.

| Tier | Goal | Status |
|------|------|--------|
| **1. Formulas** | Look up textbook/paper formulas; confirm they match the code and the Sci. Rep. paper | **GaAs (004) audit complete** — corrected \(F_0=F_h\) approximation; see [FORMULA_AUDIT.md](FORMULA_AUDIT.md) |
| **2. Constants** | Provenance table for every hard-coded constant at the correct energy, reflection, material, and temperature | **Complete for GaAs, Si, Ge, and InSb (004), 10 keV, 300 K** — see the `CONSTANTS_PROVENANCE*.md` files |
| **3. Internal checks vs external refs** | Re-score perfect-crystal Darwin curve and strained-layer Bragg shift against literature formulas and Stepanov X0h using *audited* constants | **Complete for GaAs (004), 10 keV, 300 K, σ** — FWHM and peak R match X0h/GID_sl to ~0.1%; see [CONSTANTS_SENSITIVITY.md](CONSTANTS_SENSITIVITY.md) |
| **4. Cross-code / arbitrary strain** | Compare against other dynamical XRD methods (e.g. Stepanov GID_sl / X0h, `xrayutilities`) on synthetic and paper strain fields; eventually APS/PLS data | **Two independent codes complete for all four substrates.** Ge GID_sl log RMS <= 0.0051; InSb <= 0.0299. xrayutilities matched-χ log RMS <= 0.005 for GaAs/Si/Ge/InSb. See [GID_SL_BENCHMARK.md](GID_SL_BENCHMARK.md), [GID_SL_GE_INSB_BENCHMARK.md](GID_SL_GE_INSB_BENCHMARK.md), [FIG3_GID_SL_BENCHMARK.md](FIG3_GID_SL_BENCHMARK.md), and [XU_BENCHMARK.md](XU_BENCHMARK.md). Combined three-code overview: `docs/images/benchmark_all_codes.png` (`scripts/benchmark_all_codes.py`). |

## Suggested order of work

1. Constants provenance table (Tier 2) — done
2. Formula ↔ code ↔ paper notes (Tier 1) — done
3. Score existing XRD acceptance checks against literature (Tier 3) — done
4. One online-solver comparison on a uniform strained layer (Tier 4 smoke test) — done, plus a two-step profile; see [GID_SL_BENCHMARK.md](GID_SL_BENCHMARK.md)
5. Paper Fig. 3 d'Alembert strain → production XRD vs GID_sl — done; see [FIG3_GID_SL_BENCHMARK.md](FIG3_GID_SL_BENCHMARK.md)
6. **Fig. 3 (Cr/GaAs) closure — done.** Strain morphology, XRD strain-model
   sensitivity, and the instrument-vs-1.8″ finding are summarized in
   [FIG3_CLOSURE.md](FIG3_CLOSURE.md). No new long run needed.
7. **Si (004) calculator + Fig. 2 (Cr/Si) — done (first pass).** See
   [CONSTANTS_PROVENANCE_SI.md](CONSTANTS_PROVENANCE_SI.md) and
   [FIG2_FORWARD.md](FIG2_FORWARD.md). Brillouin sidebands present; fringe
   location ~68″ vs paper ~61″ left for later tuning.
8. Second independent implementation with controlled constants
   (`xrayutilities`, all four substrates) — **done**; see [XU_BENCHMARK.md](XU_BENCHMARK.md)
9. Ge/InSb production calculators + GID_sl/X0h synthetic-layer closure —
   **done**; see [GID_SL_GE_INSB_BENCHMARK.md](GID_SL_GE_INSB_BENCHMARK.md)
10. Material-specific Cr/Ge and Cr/InSb strain profiles, then APS / PLS data
    ← **next**

## Supporting practices

- **Paper SI audit** — every parameter the paper quotes vs what the code uses
  (partially done for strain in `thermo-elastic-gaas/docs/PAPER_VALIDATION.md`).
- **Sensitivity table** — how much Darwin width / Bragg position / layer shift
  move when \(a\), \(f'\), \(f''\), or \(E\) change by their literature uncertainty.
- Keep instrument convolution **off** (`none`) for cross-code comparisons;
  add `aps_7idc` / `empirical` only when matching published plots or data.

## Scope note

The current validated scope is explicitly **GaAs, Si, Ge, and InSb (004),
10 keV, 300 K, sigma, symmetric cut**. Energy, reflection, and crystal-cut
generalization is deferred; see [SCOPE_AND_GENERALIZATION.md](SCOPE_AND_GENERALIZATION.md).

## External tools and references

The two independent dynamical-diffraction codes used for Tier-3 and Tier-4
cross-checks. Both are used with instrument response **off** so that only the
diffraction physics is compared.

### xrayutilities (`simpack.DynamicalModel`)

- Home / docs: https://xrayutilities.sourceforge.io/
- Source: https://github.com/dkriegner/xrayutilities
- Version benchmarked against: **1.7.12**
- Citation: D. Kriegner, E. Wintersberger, and J. Stangl, "xrayutilities: a
  versatile tool for reciprocal space conversion of scattering data recorded
  with linear and area detectors," *J. Appl. Cryst.* **46**, 1162–1170 (2013).
- Independent 2-beam (4-tiepoint) dynamical solver with its own scattering
  database. Note: its default omits the Debye–Waller factor, which produces a
  small flat log-offset versus our audited calculator (most visible for InSb);
  see [XU_BENCHMARK.md](XU_BENCHMARK.md).

### GID_sl / X0h (Sergey Stepanov's X-ray Server, APS)

- GID_sl tool: https://x-server.gmca.aps.anl.gov/GID_sl.html
- Server home: https://x-server.gmca.aps.anl.gov/
- Recursion-matrix method: S. A. Stepanov et al., "Dynamical x-ray diffraction
  of multilayers and superlattices: Recursion matrix extension to grazing
  angles," *Phys. Rev. B* **57**, 4829–4841 (1998).
- Server: S. Stepanov, "X-ray server: an online resource for simulations of
  X-ray diffraction and scattering," *Proc. SPIE* **5536**, 16–26 (2004).
- Cached reference curves used by our tests/figures live under `tests/data/`
  (per-material acceptance data) and `tests/data/gid_sl_combined/` (the shared
  −250…+50 arcsec grid behind `benchmark_all_codes.png`).
