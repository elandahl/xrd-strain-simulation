"""CLI entry point."""

import argparse
from pathlib import Path

from xrd_strain.config import XrdConfig
from xrd_strain.crystals.base import list_crystals
from xrd_strain.io import load_strain_profile
from xrd_strain.pipeline import run_xrd
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
        "--output",
        type=Path,
        default=Path("results/xrd_figure.png"),
        help="Path to save rocking curve figure",
    )
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()

    profile = load_strain_profile(args.strain_file)
    config = XrdConfig(crystal=args.crystal)
    result = run_xrd(profile, config=config)
    plot_xrd(result, save_path=args.output, show=not args.no_show)

    print(f"Saved XRD figure to {args.output}")
    peak_idx = result.intensity.argmax()
    print(
        f"Peak: {result.angle_deg[peak_idx]:.4f} deg, "
        f"log10 intensity {result.intensity[peak_idx]:.4f}"
    )


if __name__ == "__main__":
    main()
