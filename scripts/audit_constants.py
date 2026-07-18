#!/usr/bin/env python3
"""Recompute the XRD constants provenance comparisons.

Prints a machine-readable summary of code vs literature values for GaAs (004)
at 10 keV. See docs/CONSTANTS_PROVENANCE.md for interpretation and sources.

Does not modify any hard-coded constants — audit only.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from xrd_strain.crystals.gaas_004_dynamical import (
    A_GAAS,
    E_0,
    F_AS,
    F_GA,
    F_GAAS,
    LAMB_0,
    RE,
)
from xrd_strain.crystals.gaas_004_10kev_300k import (
    gaas_004_10kev_300k_constants,
)

REPO = Path(__file__).resolve().parents[1]
REPORT = REPO / "docs" / "constants_audit.json"

# --- Literature references (see CONSTANTS_PROVENANCE.md) --------------------
# CODATA 2022 classical electron radius (m → Å).
RE_CODATA_A = 2.8179403205e-15 * 1e10
# hc in keV·Å (CODATA/NIST-derived).
HC_KEV_A = 12.39841984
# GaAs lattice constant at 300 K (Vurgaftman et al. 2001).
A_GAAS_LIT = 5.65325
# Henke / CXRO at 10 000 eV (queried 2026-07-18): f1, f2; f' = f1 - Z.
HENKE = {
    "Ga": {"Z": 31, "f1": 28.0261, "f2": 0.5364},
    "As": {"Z": 33, "f1": 31.3953, "f2": 0.7101},
}
# Waasmaier & Kirfel (1995) five-Gaussian coefficients (ESRF DABAX).
# f0(s) = c + sum a_i exp(-b_i s^2), s = sin(theta)/lambda.
WK = {
    "Ga": {
        "a": [15.758946, 6.841123, 4.121016, 2.714681, 2.395246],
        "c": -0.847395,
        "b": [3.121754, 0.226057, 12.482196, 66.203621, 0.007238],
    },
    "As": {
        "a": [17.025642, 4.503441, 3.715904, 3.937200, 6.790175],
        "c": -2.984117,
        "b": [2.597739, 0.003012, 14.272119, 50.437996, 0.193015],
    },
}


def f0_wk(element: str, s: float) -> float:
    p = WK[element]
    return p["c"] + sum(
        ai * np.exp(-bi * s * s) for ai, bi in zip(p["a"], p["b"])
    )


def decode_code_atomic(f: complex) -> dict:
    """Decode notebook-style f = f0 + f' + i f'' stored as (f0 - |f'|) + i f''.

    The code writes ``20.6498 - 1.4384 + 1j * 0.833881``, i.e. real part is
    already f0+f'. We recover the three pieces by matching the literal
    constants in the source (known a priori).
    """
    # Literal pieces from gaas_004_dynamical.py (must stay in sync with source).
    literals = {
        complex(F_AS): (20.6498, -1.4384, 0.833881),
        complex(F_GA): (19.3993, -2.9830, 0.670502),
    }
    return {
        "f0": literals[complex(f)][0],
        "fp": literals[complex(f)][1],
        "fpp": literals[complex(f)][2],
        "f_real": float(np.real(f)),
        "f_imag": float(np.imag(f)),
    }


def main() -> int:
    production = gaas_004_10kev_300k_constants()
    s = 1.0 / (2.0 * (A_GAAS / 4.0))
    th_b_code = float(np.degrees(np.arcsin(LAMB_0 / (2.0 * (A_GAAS / 4.0)))))
    th_b_lit = float(
        np.degrees(np.arcsin((HC_KEV_A / E_0) / (2.0 * (A_GAAS_LIT / 4.0))))
    )

    ga_code = decode_code_atomic(F_GA)
    as_code = decode_code_atomic(F_AS)

    ga_f0_lit = f0_wk("Ga", s)
    as_f0_lit = f0_wk("As", s)
    ga_fp_lit = HENKE["Ga"]["f1"] - HENKE["Ga"]["Z"]
    as_fp_lit = HENKE["As"]["f1"] - HENKE["As"]["Z"]
    ga_fpp_lit = HENKE["Ga"]["f2"]
    as_fpp_lit = HENKE["As"]["f2"]

    f_ga_lit = (ga_f0_lit + ga_fp_lit) + 1j * ga_fpp_lit
    f_as_lit = (as_f0_lit + as_fp_lit) + 1j * as_fpp_lit
    F_lit = 4.0 * (f_ga_lit + f_as_lit)
    F_code = complex(F_GAAS[0])

    report = {
        "context": {
            "energy_keV": E_0,
            "reflection": "GaAs (004)",
            "sin_theta_over_lambda_Ainv": s,
            "theta_B_code_deg": th_b_code,
            "theta_B_lit_a_deg": th_b_lit,
            "bragg_offset_from_modern_a_arcsec": (th_b_code - th_b_lit) * 3600.0,
        },
        "fundamentals": {
            "re_code_A": RE,
            "re_codata_A": RE_CODATA_A,
            "re_rel_diff": (RE - RE_CODATA_A) / RE_CODATA_A,
            "lambda_code_A": LAMB_0,
            "lambda_codata_A": HC_KEV_A / E_0,
        },
        "lattice": {
            "a_code_A": A_GAAS,
            "a_lit_300K_A": A_GAAS_LIT,
            "a_rel_diff": (A_GAAS - A_GAAS_LIT) / A_GAAS_LIT,
        },
        "atomic_factors": {
            "Ga": {
                "code": ga_code,
                "lit_f0_WaasmaierKirfel": ga_f0_lit,
                "lit_fp_Henke": ga_fp_lit,
                "lit_fpp_Henke": ga_fpp_lit,
                "f0_rel_diff": (ga_code["f0"] - ga_f0_lit) / ga_f0_lit,
                "fp_rel_diff": (ga_code["fp"] - ga_fp_lit) / abs(ga_fp_lit),
                "fpp_rel_diff": (ga_code["fpp"] - ga_fpp_lit) / ga_fpp_lit,
            },
            "As": {
                "code": as_code,
                "lit_f0_WaasmaierKirfel": as_f0_lit,
                "lit_fp_Henke": as_fp_lit,
                "lit_fpp_Henke": as_fpp_lit,
                "f0_rel_diff": (as_code["f0"] - as_f0_lit) / as_f0_lit,
                "fp_rel_diff": (as_code["fp"] - as_fp_lit) / abs(as_fp_lit),
                "fpp_rel_diff": (as_code["fpp"] - as_fpp_lit) / as_fpp_lit,
            },
        },
        "structure_factor_004": {
            "code_real": F_code.real,
            "code_imag": F_code.imag,
            "code_abs": abs(F_code),
            "lit_real": F_lit.real,
            "lit_imag": F_lit.imag,
            "lit_abs": abs(F_lit),
            "abs_rel_diff": (abs(F_code) - abs(F_lit)) / abs(F_lit),
        },
        "production_300K": {
            "f_0_real": production.f_0.real,
            "f_0_imag": production.f_0.imag,
            "f_h_real": production.f_h.real,
            "f_h_imag": production.f_h.imag,
            "chi_0_real": production.chi_0.real,
            "chi_0_imag": production.chi_0.imag,
            "chi_h_real": production.chi_h.real,
            "chi_h_imag": production.chi_h.imag,
            "debye_waller_ga": production.debye_waller_ga,
            "debye_waller_as": production.debye_waller_as,
        },
        "sources": {
            "Henke_CXRO": "https://henke.lbl.gov/ (f1,f2 at 10000 eV; 2026-07-18)",
            "Waasmaier_Kirfel": "Acta Cryst. A51, 416 (1995); ESRF DABAX f0_WaasKirf.dat",
            "CODATA_re": "CODATA 2022 classical electron radius",
            "paper": "Jo et al., Sci. Rep. 12, 16606 (2022), Methods: 10 keV, APS 7ID-C, GaAs (004)",
            "a_GaAs": "Vurgaftman et al., J. Appl. Phys. 89, 5815 (2001)",
        },
    }

    def _jsonable(obj):
        if isinstance(obj, dict):
            return {k: _jsonable(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [_jsonable(v) for v in obj]
        if isinstance(obj, np.generic):
            return obj.item()
        return obj

    payload = _jsonable(report)
    REPORT.write_text(json.dumps(payload, indent=2) + "\n")

    print("GaAs (004) @ 10 keV — constants audit\n")
    print(
        f"  Bragg θ_B code a={A_GAAS} Å: {th_b_code:.4f}°  "
        f"(lit a={A_GAAS_LIT} Å: {th_b_lit:.4f}°, "
        f"offset {(th_b_code - th_b_lit) * 3600:.1f} arcsec)"
    )
    print(
        f"  r_e: code {RE:.6e}  CODATA {RE_CODATA_A:.6e}  "
        f"({100 * report['fundamentals']['re_rel_diff']:+.3f}%)"
    )
    print(
        f"  |F_004|: code {abs(F_code):.2f}  lit {abs(F_lit):.2f}  "
        f"({100 * report['structure_factor_004']['abs_rel_diff']:+.2f}%)"
    )
    print(
        f"  production: F0={production.f_0.real:.2f}+"
        f"{production.f_0.imag:.2f}i, "
        f"Fh={production.f_h.real:.2f}+{production.f_h.imag:.2f}i"
    )
    for el, code, f0l, fpl, fppl in [
        ("Ga", ga_code, ga_f0_lit, ga_fp_lit, ga_fpp_lit),
        ("As", as_code, as_f0_lit, as_fp_lit, as_fpp_lit),
    ]:
        print(
            f"  {el}: f0 {code['f0']:.4f} vs {f0l:.4f} "
            f"({100 * (code['f0'] - f0l) / f0l:+.1f}%); "
            f"f' {code['fp']:.4f} vs {fpl:.4f}; "
            f"f'' {code['fpp']:.4f} vs {fppl:.4f} "
            f"({100 * (code['fpp'] - fppl) / fppl:+.1f}%)"
        )
    print(f"\nReport saved to {REPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
