# Ge (004), 10 keV, 300 K constants

Production calculator: `ge_004_10kev`.

## Fixed scope

- diamond Ge, symmetric (004)
- 10.000 keV (`lambda = 1.239841984 Å`)
- sigma polarization
- 300 K

This is intentionally not a general Ge calculator over energy, reflection,
or surface cut.

## Constants

| quantity | value | source / role |
|---|---:|---|
| lattice constant | 5.65785 Å | room-temperature Ge reference value; sets `d004=a/4` and `theta_B=25.99376°` |
| `f0(Q004)` | 20.461447 e | Waasmaier & Kirfel (1995) five-Gaussian form factor |
| `f'` | -2.012203 e | Henke/CXRO table, interpolated at 10,000 eV |
| `f''` | 0.621959 e | Henke/CXRO table, interpolated at 10,000 eV |
| isotropic B | 0.543 Å² | Harada & Sasaki, J. Phys. Soc. Jpn. **38**, 866 (1975), 294 K |
| Debye-Waller amplitude | 0.934400 | `exp(-B [sin(theta)/lambda]^2)` |
| classical electron radius | 2.8179403205e-5 Å | CODATA |

For diamond (004), all eight atoms are in phase:

`F0 = 8 (Z + f' + i f'')`

`Fh = 8 (f0(Q004) + f' + i f'') DW`

The distinction is essential: `F0` controls forward refraction/absorption;
`Fh` controls diffraction.

## Computed values

- `chi0 = -1.826396e-5 + i 3.788015e-7`
- `chih = -1.049933e-5 + i 3.539520e-7`

## External closure

Raw curves (no angular convolution, no temporal averaging):

- Stepanov GID_sl/X0h: log10 RMS = 0.0051 (perfect), 0.0048
  (267 nm uniform layer), 0.0045 (two-step layer); correlation >= 0.999975.
- xrayutilities: production correlation >= 0.99985; when its own
  susceptibilities are injected into our solver, log10 RMS <= 0.0040.

See [GID_SL_GE_INSB_BENCHMARK.md](GID_SL_GE_INSB_BENCHMARK.md) and
[XU_BENCHMARK.md](XU_BENCHMARK.md).
