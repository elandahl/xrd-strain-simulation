# Formula audit — GaAs (004) dynamical diffraction

Tier 1 of the [external benchmarking plan](EXTERNAL_BENCHMARKS.md), focused on
the equations that materially affect new data analysis.

## Susceptibility and structure factors

Standard two-beam dynamical diffraction uses Fourier components

\[
\chi_H=-\frac{r_e\lambda^2}{\pi V_c}F_H,\qquad
F_H=\sum_n f_n(H)\exp(iH\cdot r_n-W_n(H)),
\]

with

\[
f_n(H)=f_{0,n}(H)+f'_n(E)+if''_n(E).
\]

References:

- S. A. Stepanov, [X0h / X-Ray Server](https://x-server.gmca.aps.anl.gov/x0h.html)
- R. R. Lindberg & Yu. V. Shvyd'ko, *Two-beam Bragg diffraction, filtering,
  and delays*, ANL technical report (2022), DOI 10.2172/2007525
- B. W. Batterman & H. Cole, *Rev. Mod. Phys.* **36**, 681 (1964)

### The corrected distinction

\(F_0\) and \(F_h\) are not interchangeable:

- **Forward component \(F_0\):** \(H=0\), so \(f_0(0)=Z\) and
  \(W(0)=0\). It controls mean refraction and absorption.
- **Diffracting component \(F_h\):** \(H=(004)\), so it uses
  \(f_0(Q_h)<Z\) and species-specific Debye–Waller attenuation. It controls
  coupling into the Bragg beam and the Darwin width.

The notebook used `F_GAAS[0] = F_GAAS[1]`, effectively \(F_0=F_h\). This made
\(|\chi_0|\) much too small. The production calculator now uses:

| Quantity | Production | Stepanov X0h |
|---|---:|---:|
| \(\mathrm{Re}\chi_0\) | −1.814×10⁻⁵ | −1.804×10⁻⁵ |
| \(\mathrm{Im}\chi_0\) | +3.805×10⁻⁷ | +3.609×10⁻⁷ |
| \(\mathrm{Re}\chi_h\) | −1.032×10⁻⁵ | −1.036×10⁻⁵ |
| \(\mathrm{Im}\chi_h\) | +3.555×10⁻⁷ | +3.372×10⁻⁷ |

The imaginary differences reflect Henke versus X0h/International-Tables
dispersion databases and are within the documented database spread.

## Zincblende GaAs (004)

GaAs has four Ga and four As atoms in the conventional cubic cell. For (004),
all basis phases are unity:

\[
F_{004}=4(f_\mathrm{Ga}+f_\mathrm{As}).
\]

For the isotropic harmonic model used here, \(F_h=F_{-h}\). This justifies the
current scalar fourth-order solver's single diffracting susceptibility for this
specific reflection. It is **not** a safe generalization to arbitrary
non-centrosymmetric reflections; a generalized future solver must carry
\(\chi_h\) and \(\chi_{-h}\) separately.

## Thermal motion

The production \(F_h\) applies

\[
T_n(H)=\exp[-B_n(\sin\theta/\lambda)^2]
\]

using 300 K values \(B_\mathrm{Ga}=0.622\) and
\(B_\mathrm{As}=0.483\) Å² (Stevenson, *Acta Cryst.* **A50**, 621–632,
1994). \(F_0\) correctly has no Debye–Waller attenuation.

## Geometry and polarization

- Symmetric (004): \(H=(0,0,-8\pi/a)\), surface normal \((0,0,1)\).
- Current solver is scalar **σ polarization**, for which \(P=1\).
- For π polarization, standard two-beam theory uses
  \(P=\cos(2\theta_B)\). This has not been added as an unvalidated switch;
  the current production backend is explicitly 10 keV, 300 K, σ polarized.
- Intensity includes the flux factor \(\gamma_h/\gamma_0\).

## Strain mapping

Each layer uses \(a_z\to a_z+\Delta a_z\), with
\(\epsilon=\Delta a_z/a_z\), shifting the local reciprocal lattice and hence
the Bragg condition. The independent acceptance
\(\Delta\theta=-\epsilon\tan\theta_B\) agrees to 2–3% for finite strained caps
and scales linearly with strain.

## External closure

For a semi-infinite perfect GaAs crystal at 10 keV, σ polarization:

| Metric | Production code | X0h/GID_sl | Difference |
|---|---:|---:|---:|
| Absorbing-curve FWHM | 5.732″ | 5.73686″ | −0.09% |
| Peak reflectivity | 0.97389 | 0.97283 | +0.11% |
| Peak offset from kinematic Bragg | 2.27″ | 2.27″ | sampling-level |

This validates the fourth-order boundary implementation for the present
GaAs(004), 10 keV, σ-polarized scope. It does not yet validate arbitrary
reflection, polarization, energy, or material.

## Remaining formula work before generalization

1. General \(\chi_h,\chi_{-h}\) and basis phases for arbitrary reflections.
2. Explicit σ/π polarization in the solver interface.
3. Energy-dependent scattering-factor provider rather than a fixed audited
   10 keV constants set.
4. Temperature-dependent \(a(T)\) and \(B_n(T)\), with provenance.
5. Cross-code arbitrary-strain comparison against GID_sl or another
   Takagi–Taupin implementation (Tier 4).

