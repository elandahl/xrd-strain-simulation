# InSb (004), 10 keV, 300 K constants

Production calculator: `insb_004_10kev`.

## Fixed scope

- zincblende InSb, symmetric (004)
- 10.000 keV (`lambda = 1.239841984 Å`)
- sigma polarization
- 300 K

InSb is the highest-risk material in the present four-substrate set: it has
stronger absorption, a lower Bragg angle, softer acoustics, and a larger,
less precisely known thermal-motion correction.

## Constants

| quantity | In | Sb | source / role |
|---|---:|---:|---|
| lattice constant | 6.47937 Å | — | 300 K reference; `theta_B=22.50131°` |
| `f0(Q004)` | 33.601616 e | 35.023630 e | Waasmaier & Kirfel (1995) |
| `f'` | +0.081548 e | +0.075377 e | Henke/CXRO, interpolated at 10,000 eV |
| `f''` | 3.419065 e | 3.918354 e | Henke/CXRO, interpolated at 10,000 eV |
| effective isotropic B | 1.14730 Å² | 1.14730 Å² | Debye model, 300 K, `theta_D=160 K`, mean In/Sb mass |
| Debye-Waller amplitude | 0.896450 | 0.896450 | `exp(-B [sin(theta)/lambda]^2)` |

For zincblende (004), all four In and four Sb atoms are in phase:

`F0 = 4 [(Z_In + f'_In + i f''_In) + (Z_Sb + f'_Sb + i f''_Sb)]`

`Fh = 4 [(f0_In + f'_In + i f''_In) DW_In + (f0_Sb + f'_Sb + i f''_Sb) DW_Sb]`

## Computed values

- `chi0 = -2.030754e-5 + i 1.487714e-6`
- `chih = -1.250196e-5 + i 1.333661e-6`

## Debye-Waller uncertainty

Reid & Pirie, Acta Cryst. A **39**, 1-13 (1983), show that species-resolved
InSb B factors vary among lattice-dynamical models because their eigenvectors
differ. The current common effective B is physically defensible and explicit,
but should be treated as a sensitivity parameter when fitting high-quality
InSb data. It is not hidden inside the solver.

## External closure

Raw curves (no angular convolution, no temporal averaging):

- Stepanov GID_sl/X0h: log10 RMS = 0.0290 (perfect), 0.0288
  (267 nm uniform layer), 0.0299 (two-step); correlation >= 0.999924.
- xrayutilities: production correlation >= 0.99949; when its own
  susceptibilities are injected into our solver, log10 RMS <= 0.0032.

The larger absolute residual than Ge is a scattering-database / thermal-factor
difference, not a dynamical-solver difference. See
[GID_SL_GE_INSB_BENCHMARK.md](GID_SL_GE_INSB_BENCHMARK.md) and
[XU_BENCHMARK.md](XU_BENCHMARK.md).
