# Current scope and deferred generalization

## Production scope today

The registered production calculators are intentionally explicit:

- `gaas_004_10kev`
- `si_004_10kev`
- `ge_004_10kev`
- `insb_004_10kev`

They all represent 300 K, sigma-polarized, symmetric (004) diffraction at
10 keV from a (001)-cut cubic substrate. The calculator name is the contract;
the implementation must not imply that arbitrary energy, hkl, or cut is
already supported.

## Why not generalize yet

The current backend is safely reusable across these four fixed cases, but
three future dimensions change different physics:

1. **Energy** changes wavelength and energy-dependent anomalous factors
   `f'`/`f''`, especially near absorption edges.
2. **Reflection** changes `Q`, atomic phases, polarization, and potentially
   the relation between `Fh` and `F-h`.
3. **Crystal cut** changes the surface normal, incidence/exit direction
   cosines, and asymmetry factor; `(004)` is symmetric only for the present
   (001) cut.

Prematurely combining these into one flexible configuration would make it
easy to construct unsupported combinations that look valid.

## Future seams

When the data require them:

- replace fixed constants modules with an energy-aware scattering provider;
- separate reflection `hkl` from surface-normal/cut orientation;
- generalize the dynamical backend from its fixed symmetric-(004) reciprocal
  vector and surface normal;
- validate every new energy/reflection/geometry against raw X0h/GID_sl and
  xrayutilities curves before applying angular or temporal instrument models.

The present material calculators and crystal registry remain valid entry
points for that later work; no migration is required now.
