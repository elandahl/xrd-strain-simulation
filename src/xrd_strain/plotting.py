"""XRD plotting."""

from pathlib import Path

import matplotlib.pyplot as plt

from xrd_strain.models import XrdResult


def plot_xrd(result: XrdResult, save_path: str | Path | None = None, show: bool = True) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(result.angle_deg, result.intensity, "-k")
    ax.set_xlabel("Incident angle (deg)", size=15)
    ax.set_ylabel("Intensity", size=15)
    ax.tick_params(labelsize=15)
    ax.set_title(f"{result.crystal} / {result.substrate_material} / {result.strain_model}")
    fig.tight_layout()

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150, bbox_inches="tight")

    if show:
        plt.show()

    return fig
