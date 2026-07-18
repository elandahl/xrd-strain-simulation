# Constants provenance — GaAs (004) @ 10 keV

Tier 2 of the [external benchmarking plan](EXTERNAL_BENCHMARKS.md).

Every hard-coded constant in `src/xrd_strain/crystals/gaas_004_dynamical.py`
(and the related instrument defaults) is listed with:

- the **code value**,
- the **paper** value when Jo et al., Sci. Rep. **12**, 16606 (2022) states one,
- an **external literature** reference at the correct energy / reflection /
  material / temperature,
- a **verdict**.

Regenerate the numeric comparisons with:

```bash
python scripts/audit_constants.py
```

## Experimental context (paper ↔ code)

| Quantity | Code | Paper (Methods) | Verdict |
|----------|------|-----------------|---------|
| Photon energy | `E_0 = 10.0` keV | **10 keV** X-rays, APS 7ID-C, diamond (111) DCM | **Match** |
| Reflection | GaAs (004), symmetric Bragg | GaAs (0 0 4) | **Match** |
| Probe temperature (ambient) | implied room temp; strain solver `T0 = 300` K | Ambient / room-temperature experiment | **Consistent** |
| Stated angular resolution | `aps_7idc` default FWHM **1.8 arcsec** (= 0.5 mdeg) | **0.5 millidegrees** | **Match** (physical instrument). The default `empirical` kernel is a separate, broader *effective* resolution — see [INSTRUMENTS.md](INSTRUMENTS.md). |

## Fundamental / conversion constants

| Symbol | Code | Literature | Source | Verdict |
|--------|------|------------|--------|---------|
| Classical electron radius \(r_e\) | `RE = 2.818e-5` Å | \(2.81794032\times10^{-5}\) Å | CODATA 2022 | **OK** (4 significant figures) |
| \(hc\) (for \(\lambda=hc/E\)) | `12.3984` keV·Å | \(12.39841984\ldots\) keV·Å | CODATA / NIST \(hc\) | **OK** |
| \(\lambda\) at 10 keV | `1.23984` Å | \(1.23984198\ldots\) Å | from \(hc/E\) | **OK** |
| m → Å | `1e10` | exact by definition | — | **OK** |

## Lattice geometry

| Symbol | Code | Literature (300 K) | Source | Verdict |
|--------|------|--------------------|--------|---------|
| \(a_\mathrm{GaAs}\) | `A_GAAS = 5.65` Å | \(5.65325\) Å (Vurgaftman et al.); commonly \(5.6533\) Å | Vurgaftman et al., *J. Appl. Phys.* **89**, 5815 (2001); standard RT compilations | **Rounded.** \(\Delta a/a \approx 5.8\times10^{-4}\) → Bragg-angle offset ≈ **58 arcsec** at 10 keV (004). Negligible for *relative* \(\theta-\theta_B\) plots; matters for absolute θ. |
| Reciprocal lattice vector (004) | `H_GAAS = (8π/a)\,[0,0,-1]` | \(\lvert\mathbf{H}_{004}\rvert = 8\pi/a = 4\cdot(2\pi/a)\) | zincblende (004) | **Correct** |
| Surface normal | `[0,0,1]` | symmetric Bragg | — | **OK** |
| Structure-factor phase for (004) | `F_GAAS[0] = F_GAAS[1] = 4(F_Ga+F_As)` | For zincblende, \(h+k+l=4n\) → all atoms in phase; \(F_0 = F_H = 4(f_\mathrm{Ga}+f_\mathrm{As})\) | International Tables / Authier | **Correct** |

Kinematic Bragg angle with the code lattice constant:

\[
\theta_B = \arcsin\!\bigl(\lambda / (2\,d_{004})\bigr),\quad d_{004}=a/4
\]

→ \(\theta_B = 26.0325^\circ\) (matches the internal acceptance suite). With
\(a=5.65325\) Å it would be \(26.0164^\circ\).

## Atomic scattering factors at 10 keV, (004)

The code stores atomic factors as

```text
F_AS = 20.6498 - 1.4384 + 1j * 0.833881
F_GA = 19.3993 - 2.9830 + 1j * 0.670502
```

interpreted as \(f_0(Q) + f'(E) + i\,f''(E)\), then
\(F_{0,H}=4(f_\mathrm{Ga}+f_\mathrm{As})\).

For GaAs (004) with \(a=5.65\) Å:

\[
\sin\theta/\lambda = 1/(2d_{004}) = 0.3540\,\mathrm{Å}^{-1}.
\]

| Piece | Code | External reference | Source | Verdict |
|-------|------|--------------------|--------|---------|
| \(f_0^\mathrm{Ga}(Q)\) | 19.3993 | **19.72** | Waasmaier & Kirfel, *Acta Cryst.* **A51**, 416 (1995); ESRF DABAX `f0_WaasKirf.dat` | **~1.6% low** |
| \(f_0^\mathrm{As}(Q)\) | 20.6498 | **21.07** | same | **~2.0% low** |
| \(f'_\mathrm{Ga}(10\,\mathrm{keV})\) | −2.983 | **−2.974** (\(f_1-Z\), \(f_1=28.03\)) | Henke / CXRO (LBNL) elemental table at 10 000 eV | **Excellent match** |
| \(f''_\mathrm{Ga}(10\,\mathrm{keV})\) | 0.6705 | **0.536** (\(f_2\)) | Henke / CXRO | **~25% high** |
| \(f'_\mathrm{As}(10\,\mathrm{keV})\) | −1.438 | **−1.605** (\(f_1-Z\), \(f_1=31.40\)) | Henke / CXRO | **~10% less negative** |
| \(f''_\mathrm{As}(10\,\mathrm{keV})\) | 0.834 | **0.710** | Henke / CXRO | **~17% high** |

Combined structure factor magnitude \(|F_{004}|\):

| | Real | Imag | \(\lvert F\rvert\) |
|--|------|------|---------------------|
| Code | 142.51 | 6.02 | **142.64** |
| Lit (WK \(f_0\) + Henke \(f',f''\)) | 144.83 | 4.99 | **144.91** |
| Relative difference | | | **~1.6%** |

### Notes on the anomalous parts

- Ga **K-edge ≈ 10.367 keV**. At the paper’s **10.000 keV**, \(f'\) is only
  mildly anomalous (~−3). Near the edge it would plunge (Henke \(f'\sim-11\)
  at 10.367 keV). The code’s Ga \(f'\) is clearly a **10 keV** value, not an
  edge value — consistent with the paper.
- The \(f''\) values in the code are systematically **higher** than Henke.
  Higher \(f''\) → more absorption → slightly broader / lower Darwin peak.
  Source of the notebook values is not yet identified (possibly an older
  Cromer–Liberman / Sasaki table, or a rounded composite). Flagged for
  sensitivity work in Tier 3.
- No Debye–Waller factor is applied; at 300 K for GaAs (004) that would
  reduce \(|F|\) by a few percent. Opposite direction from the \(f_0\)
  underestimate, partially compensating.

## Instrument defaults (not optics constants, but provenance)

| Parameter | Code | Paper | Verdict |
|-----------|------|-------|---------|
| `aps_7idc` FWHM | 1.8 arcsec | 0.5 mdeg ≈ 1.8 arcsec | **Match** |
| `empirical` multi-Gaussian | `DetectorConfig` (dominant σ ≈ 22 arcsec) | not a stated instrument number | **Fit-effective**; see INSTRUMENTS.md |

## Summary of findings

**Solid matches (no action needed for physics):**

- Photon energy 10 keV, GaAs (004), APS 7ID-C geometry
- \(r_e\), \(hc/\lambda\), reciprocal-lattice / structure-factor phase
- Ga \(f'\) vs Henke at 10 keV

**Documented discrepancies (do not change code yet — measure impact in Tier 3):**

1. **Lattice constant** rounded to 5.65 Å → ~58″ absolute Bragg offset vs
   modern \(a(300\,\mathrm{K})\).
2. **Thomson \(f_0\)** ~1.6–2% below Waasmaier–Kirfel at \(\sin\theta/\lambda=0.354\).
3. **Absorption \(f''\)** ~17–25% above Henke for both Ga and As; As \(f'\) ~10% off.
4. Net \(|F_{004}|\) ~1.6% low vs modern tables — likely a small effect on Darwin
   width / peak reflectivity, to be quantified next.

**Not bugs in the sense of wrong energy or wrong reflection** — the
calculator is clearly built for GaAs (004) at 10 keV, matching the paper.
The open question is how much the table-vintage differences move the rocking
curve relative to experimental noise and instrument blur.

## Next steps (Tier 3)

1. Recompute Darwin width and peak reflectivity with Henke+WK constants;
   compare to the current acceptance numbers (~5.9″ FWHM, R≈0.985).
2. Sensitivity: vary \(a\), \(f_0\), \(f''\) one-at-a-time within the
   discrepancies above; record Δ(FWHM), Δ(peak), Δ(layer shift).
3. Only then decide whether to update the hard-coded factors or keep them as
   the archival notebook values with a documented offset.
