# Instrument response models

A time-resolved rocking curve measurement has **two** instrument responses,
and the forward model must apply both before comparing with data:

| response | physical origin | acts on | implementation |
|---|---|---|---|
| **Angular** | monochromator / analyzer / detector acceptance | the rocking curve, convolution in θ | `XrdConfig(instrument=...)`, this page |
| **Temporal** | x-ray bunch duration (+ timing jitter) | the signal vs pump-probe delay | `xrd_strain.temporal`, see below |

The two are independent, both linear on intensity, and therefore commute —
but both must be applied on *linear* intensity, before any log₁₀.

## Temporal response (delay averaging)

The probe is not a delta function in time. In the APS 24-bunch mode used for
the Sci. Rep. 2022 experiment, the x-ray bunch duration was **~90 ps FWHM**
(`temporal.APS_24BUNCH_FWHM_PS`); timing jitter is negligible in comparison.
During one bunch the acoustic wavefront moves \(v \cdot 90\,\mathrm{ps}\)
(~760 nm in Si), so the measured curve is an incoherent average over delay:

\[ I_{\mathrm{meas}}(\theta; t_0) = \sum_k w_k\, I(\theta; t_k) \]

This cannot be applied to a single strain snapshot. Compute strain profiles
at several delays and average:

```python
from xrd_strain import gaussian_delay_weights, run_xrd_delay_averaged

times, weights = gaussian_delay_weights(0.34e-9, 90e-12, n_samples=9)
# ... run the strain model at each time in `times`, load the profiles ...
result = run_xrd_delay_averaged(profiles, weights, config=config)
```

Consequences for interpretation (demonstrated in
[FIG2_FORWARD.md](FIG2_FORWARD.md), `scripts/fig2_delay_average.py`):

- coherent **depth fringes** (spacing set by the instantaneous strained
  depth) scramble and cancel;
- **Brillouin sidebands** survive, because their angle is set by the acoustic
  pulse-train period, which is delay-independent.

Different sources need different kernels: APS 24-bunch ≈ 90 ps FWHM; the
Pohang Light Source dense-fill data will have its own bunch duration.

## Angular response models

The computed dynamical-diffraction rocking curve can be convolved with an
angular-resolution model before output. Select with `--instrument` on the CLI
or `XrdConfig(instrument=...)`.

| instrument | kernel | meaning |
|---|---|---|
| `aps_7idc` (default) | single Gaussian, `--instrument-fwhm-arcsec` (default 1.8) | The **physical** instrument resolution measured/stated for APS 7ID-C: 0.5 mdeg ≈ 1.8 arcsec. Analysis default. |
| `empirical` | multi-Gaussian, dominant σ ≈ 22 arcsec | **Effective** resolution (instrument + sample + fit) inherited from the original notebook / thermo-elastic-gaas (tag `paper-v1.0`). Heavy smoothing. Legacy paper-figure reproduction option. |
| `none` | — | Raw dynamical-diffraction curve, no convolution. |

> Deprecated alias: `--instrument notebook` still works and maps to
> `empirical` (see `pipeline.resolve_instrument`). It was renamed because it is
> not the physical instrument function — see below.

```bash
python scripts/run.py --strain-file strain.npz \
  --instrument aps_7idc --instrument-fwhm-arcsec 1.8 \
  --angle-min 25.975 --angle-max 26.091 --n-points 581 \
  --save-arrays results/curve.npz --no-show
```

## Why `empirical`, not `notebook` or `instrument`

The original notebook applied a **fitted multi-Gaussian** to match a measured
rocking curve. Its widths (`DetectorConfig`, in degrees of σ) are
approximately 4.8, **21.6**, 5.8, 5.4 arcsec — the dominant component is
~22 arcsec, i.e. **more than 10× wider** than the 1.8 arcsec the paper quotes
as the APS 7ID-C angular resolution.

That width therefore cannot be pure instrument resolution. It is an
*effective* broadening that lumps together several real contributions the
1.8 arcsec number omits:

- crystal mosaic spread / sample imperfection,
- probe-depth averaging over the strained profile,
- beam divergence and energy bandwidth,
- smoothing implicit in the published *fit* line (vs raw data points).

Hence the name `empirical` (empirically-fit effective resolution) rather than
`notebook` (where it came from) or `instrument` (which it is not, exclusively).

## Does the empirical model match the paper better than 1.8 arcsec? Yes.

Comparing against the published Fig. 3 rocking curve:

- `aps_7idc` at the paper's stated 1.8 arcsec is a *light* blur: it leaves a
  sharp Bragg cusp and crisp thickness fringes that are **not** visible in the
  published smooth fit curve.
- `empirical` (σ ≈ 22 arcsec) smooths the peak and fringes and reproduces the
  **smoothness** of the published curve much more closely.

**Conclusion:** the smooth curve plotted in the paper is broader than the
paper's own stated 1.8 arcsec instrument resolution. The extra width comes
from sample/experimental/fit effects (listed above), which `empirical`
captures phenomenologically. So `empirical` is the right **explicit option**
for reproducing the paper's *figure*, while the physical `aps_7idc` response is
the default for new data analysis.

## Reading the validation-matrix rocking-curve figure

Benchmark target — the inset rocking curve of the published Fig. 3 (Jo et al.,
Sci. Rep. **12**, 16606 (2022), CC BY 4.0): log intensity vs θ − θ_B, with the
thermal contribution enhancing the lower-angle side:

![Published Sci. Rep. Fig. 3](images/published_fig3_scirep_2022.png)

Our computed rocking curves:

![rocking curve matrix](images/matrix_rocking.png)

This figure is a 2×2 matrix (2 strain models × 2 instruments) drawn as **two
panels sharing the same two colored lines**:

- **PANELS = instrument.** Left = `empirical` (σ ≈ 22 arcsec, heavy blur);
  right = `aps_7idc` (1.8 arcsec FWHM, light blur).
- **COLORS = strain model.** Blue = leapfrog (`ttm_cr_gaas`); red = d'Alembert
  (`ttm_dalembert_cr_gaas`). **Both** colors appear in **both** panels.

Two things explain why the panels look different, and neither is the strain
model:

1. **Vertical offset (left panel sits ~2 decades higher).** The `empirical`
   convolution ends with `y / len(y)` (a normalization quirk), which shifts the
   absolute log₁₀ scale. Only curve *shape* is meaningful on a log plot — ignore
   the panel-to-panel vertical position.
2. **Peak sharpness.** The wide `empirical` kernel rounds the Bragg peak and
   washes out fringes; the narrow `aps_7idc` kernel leaves a sharper cusp and
   more shoulder structure.

The physics takeaway: within each panel, blue vs red differ by at most
~0.6–0.8 in log₁₀ I and only in the weak +50…+150 arcsec shoulder. The Bragg
peak and low-angle asymmetry — the features used for paper validation — are
unchanged by the strain-model choice.

## Sampling requirement

The convolution acts on the discrete angular grid, so the step must resolve
the kernel: step < σ. For `aps_7idc` (σ ≈ 0.76 arcsec) use ≤ ~0.7 arcsec/pt.
At the historical default (0.1 deg / 100 pts = 3.6 arcsec/pt) a 1.8 arcsec
kernel is effectively an identity, and `apply_gaussian_instrument` degrades
gracefully to almost-no-op.

Fine sampling used to be memory-prohibitive: the original whole-array
`xrd_slab_gaas` allocates O(n_layers × n_angles) 4×4 complex matrices (>1 GB
for ~3700 layers × ~600 angles). The calculator now uses
`xrd_slab_gaas_lowmem`, an angle-by-angle wrapper that is numerically
identical (asserted in `tests/test_smoke.py`) with flat memory use.
