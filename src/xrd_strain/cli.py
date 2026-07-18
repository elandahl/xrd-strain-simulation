"""CLI entry point."""

import argparse
import json
from pathlib import Path

import numpy as np

from xrd_strain.config import XrdConfig
from xrd_strain.crystals.base import list_crystals
from xrd_strain.io import load_strain_profile
from xrd_strain.pipeline import INSTRUMENTS, run_xrd
from xrd_strain.plotting import plot_xrd


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute XRD rocking curve from a strain profile."
    )
    parser.add_argument(
        "--strain-file",
        type=Path,
        required=True,
        help="Path to strain profile .npz from strain-wave-simulation",
    )
    parser.add_argument(
        "--crystal",
        default="gaas_004_10kev",
        choices=list_crystals(),
        help="Crystal/reflection/energy calculator",
    )
    parser.add_argument(
        "--instrument",
        default="aps_7idc",
        choices=list(INSTRUMENTS) + ["notebook"],
        help="Angular resolution model; 'notebook' is a deprecated alias for "
        "'empirical' (see docs/INSTRUMENTS.md)",
    )
    parser.add_argument(
        "--instrument-fwhm-arcsec",
        type=float,
        default=1.8,
        help="Gaussian FWHM in arcsec for --instrument aps_7idc (default 1.8)",
    )
    parser.add_argument("--angle-min", type=float, default=25.98)
    parser.add_argument("--angle-max", type=float, default=26.08)
    parser.add_argument("--n-points", type=int, default=100)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/xrd_figure.png"),
        help="Path to save rocking curve figure",
    )
    parser.add_argument(
        "--save-arrays",
        type=Path,
        default=None,
        help="Optional path to save angle/intensity arrays as .npz",
    )
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()

    profile = load_strain_profile(args.strain_file)
    config = XrdConfig(
        crystal=args.crystal,
        instrument=args.instrument,
        instrument_fwhm_arcsec=args.instrument_fwhm_arcsec,
        angle_min=args.angle_min,
        angle_max=args.angle_max,
        n_points=args.n_points,
    )
    result = run_xrd(profile, config=config)
    plot_xrd(result, save_path=args.output, show=not args.no_show)

    print(f"Saved XRD figure to {args.output}")
    peak_idx = result.intensity.argmax()
    print(
        f"Peak: {result.angle_deg[peak_idx]:.4f} deg, "
        f"log10 intensity {result.intensity[peak_idx]:.4f} "
        f"(instrument={result.instrument}, strain model={result.strain_model})"
    )

    if args.save_arrays is not None:
        args.save_arrays.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            args.save_arrays,
            angle_deg=result.angle_deg,
            intensity=result.intensity,
            meta_json=json.dumps(
                {
                    "crystal": result.crystal,
                    "instrument": result.instrument,
                    "strain_model": result.strain_model,
                    "substrate_material": result.substrate_material,
                    "instrument_fwhm_arcsec": args.instrument_fwhm_arcsec,
                    "strain_file": str(args.strain_file),
                }
            ),
        )
        print(f"Saved arrays to {args.save_arrays}")


if __name__ == "__main__":
    main()
