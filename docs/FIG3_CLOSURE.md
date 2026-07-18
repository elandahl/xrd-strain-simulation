# Fig. 3 closure (Sci. Rep. 2022, GaAs)

Summary of what has and has not changed substantially in the modern split-repo
stack, relative to the archival notebook path (`thermo-elastic-gaas`, tag
`paper-v1.0`), for reproducing **Figure 3** of Jo et al., *Sci. Rep.* **12**,
16606 (2022): the Cr/GaAs strain profile and the GaAs (004) rocking curve at
Δt = 1.8 ns.

This is a documentation checkpoint. It draws on results and figures already in
the three repositories; **no new long simulation is required** to reach these
conclusions.

## Validated modern stack

- **Strain:** `ttm_dalembert_cr_gaas` (default), preset `paper_fig3_gaas`
  (`strain-wave-simulation`).
- **XRD:** `gaas_004_10kev` production calculator — distinct \(F_0\)/\(F_h\),
  modern Waasmaier–Kirfel / Henke / Stevenson constants
  (`xrd-strain-simulation`).
- **Instrument kernels:** `none`, `aps_7idc` (1.8″ FWHM, physical APS 7ID-C
  resolution), `empirical` (dominant σ ≈ 22″, effective instrument + sample +
  fit broadening).

## What changed substantially

### 1. Far-field strain: dispersive wake → d'Alembert pulse train

The archival leapfrog solver advances the elastic wave on the thermal-diffusion
time step (acoustic Courant number C ≈ 0.003). Over the ~8.5 μm propagation
path this produces a **numerical-dispersion wake** (~35–55 nm chirped
oscillations) behind the wavefront instead of physics.

The `ttm_dalembert_cr_gaas` model keeps the identical TTM + FD near field but
translates the right-going wave into the far field with an exact,
dispersion-free d'Alembert shift. It recovers the **discrete ~114 nm acoustic
pulse train** (Cr-film round-trip echoes, `2·L_film·v_GaAs/v_Cr`) that the
paper shows. **This is a substantial change and it matches the published Fig. 3
morphology.**

- Figure: `strain-wave-simulation/docs/images/fig3_strain_split_dalembert.png`
  (corrected strain on the paper's split z-axis), and `matrix_strain_far.png`
  (leapfrog vs d'Alembert, near vs far field).

### 2. Production vs legacy XRD constants (absolute scale)

The production `gaas_004_10kev` calculator uses **distinct \(F_0\) and \(F_h\)**
and modern constants, which fixes the **absolute Darwin width and peak
reflectivity** to match the independent Stepanov X0h / GID_sl reference (peak
R and FWHM agree to ~0.1%; see [FIG3_GID_SL_BENCHMARK.md](FIG3_GID_SL_BENCHMARK.md)
and [CONSTANTS_SENSITIVITY.md](CONSTANTS_SENSITIVITY.md)).

This is the correct, well-provenanced absolute normalization. **Its effect on
the heavily broadened Fig. 3 inset is secondary:** once the curve is convolved
with the paper-like instrument kernel, the small absolute-scale/width
corrections are washed out and the inset match is dominated by the instrument
convolution (below).

## What did *not* change substantially

### 3. Rocking curve: strain-model impact is modest at 1.8 ns

Within a fixed instrument, the leapfrog vs d'Alembert strain choice changes the
computed rocking curve by at most ~0.6–0.8 in log₁₀ intensity, and only in the
weak +50…+150″ shoulder. The Bragg peak and the low-angle asymmetry — the
features used for paper validation — are essentially unchanged. The far-field
dispersion artifact never compromised the diffraction validation, because the
X-ray extinction depth weights the near field.

### 4. Instrument kernel dominates the paper-looking match

The **instrument kernel**, not the strain model, controls how much the computed
curve resembles the published inset:

- `aps_7idc` at the paper's stated 1.8″ leaves a sharp Bragg cusp and crisp
  thickness fringes that are **not** present in the published smooth curve.
- `empirical` (σ ≈ 22″) reproduces the **smoothness** of the published inset
  much more closely.

**The empirical model matches the published inset better than the stated 1.8″
does.** The published curve is therefore broader than pure APS 7ID-C
resolution; the extra width is real sample/experiment/fit broadening (mosaic
spread, probe-depth averaging, beam divergence, energy bandwidth, fit
smoothing) that the multi-Gaussian `empirical` model lumps together. See
[INSTRUMENTS.md](INSTRUMENTS.md).

## Existing evidence figures

| Figure | Location | Shows |
|--------|----------|-------|
| `published_fig3_scirep_2022.png` | all three repos' `docs/images/` | The benchmark target (CC BY 4.0). |
| `fig3_strain_split_dalembert.png` | `strain-wave-simulation/docs/images/` | Corrected d'Alembert strain on the paper's split z-axis. |
| `matrix_strain_far.png` | strain + thermo `docs/images/` | Leapfrog vs d'Alembert, near vs far field. |
| `matrix_rocking.png` | xrd + thermo `docs/images/` | 2×2 rocking-curve matrix (panels = instrument, colors = strain model). |
| `fig3_gid_sl_benchmark.png` | `xrd-strain-simulation/docs/images/` | Production vs GID_sl on the Fig. 3 d'Alembert strain (no instrument). |

## Status

**Fig. 3 (Cr/GaAs) closure is complete** at the documentation level with the
validated modern stack. The remaining production-vs-reference differences are
understood (scattering database, not the strain representation or the dynamical
algorithm) and do not require a new multi-hour run.

## Next paper target: Fig. 2 (Cr/Si)

The next paper target is **Figure 2 (Cr/Si substrate)**. That requires a
**Si (004) XRD calculator** (the current calculator is GaAs-only) plus the
Cr/Si strain case. See [EXTERNAL_BENCHMARKS.md](EXTERNAL_BENCHMARKS.md) for the
ordered plan.

![Published Sci. Rep. 2022 Figure 2 (Cr/Si) — next target](images/published_fig2_scirep_2022.png)
