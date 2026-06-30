"""XRD configuration."""

from dataclasses import dataclass


M_TO_ANGSTROM = 1e10


@dataclass
class XrdConfig:
    """Parameters for XRD rocking curve simulation."""

    crystal: str = "gaas_004_10kev"
    angle_min: float = 25.98
    angle_max: float = 26.08
    n_points: int = 100
    strain_eps: float = 1e-6
    log10_intensity: bool = True


@dataclass
class DetectorConfig:
    """Instrument resolution convolution parameters."""

    a1: float = 0.418772
    a2: float = 0.038453
    a3: float = 0.522
    cen1: float = 0.00135
    cen2: float = 0.0070e-2
    cen3: float = -0.7938e-3
    cen4: float = -4.6038e-3
    p1: float = 1.33e-3
    p2: float = 6.0094e-3
    p3: float = 1.60414e-3
    p4: float = 1.50414e-3

    def as_tuple(self) -> tuple:
        return (
            self.a1,
            self.a2,
            self.a3,
            self.cen1,
            self.cen2,
            self.cen3,
            self.cen4,
            self.p1,
            self.p2,
            self.p3,
            self.p4,
        )
