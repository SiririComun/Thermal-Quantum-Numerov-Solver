"""Configuration dataclasses for dimensionless Schrödinger simulations.

The legacy notebook uses reduced (dimensionless) units where:
- $\\hbar = 1$
- $m = 0.5$
- $k_B = 1$
- $L = \\pi$

These values make the infinite-well reference spectrum satisfy
$E_n = n^2$ in units of the ground-state energy scale.
"""

from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True, slots=True)
class PhysicsConfig:
    """Physical constants in dimensionless units from the legacy prototype.

    The values match the normalization used in
    ``legacy/research_prototype.ipynb`` so that the 1D infinite square-well
    reference has energies :math:`E_n = n^2`.
    """

    #: Reduced Planck constant in dimensionless units.  Default: ``1.0``.
    hbar: float = 1.0
    #: Particle mass in dimensionless units.  Default: ``0.5``, which ensures
    #: :math:`\\hbar^2 / (2m) = 1` when ``hbar = 1``.
    mass: float = 0.5
    #: Boltzmann constant in dimensionless thermal units.  Default: ``1.0``,
    #: so that reduced temperature equals :math:`k_B T`.
    kb: float = 1.0
    #: Characteristic well width :math:`L` in dimensionless length units.
    #: Default: :math:`\\pi`.
    well_width: float = math.pi


@dataclass(frozen=True, slots=True)
class NumericalConfig:
    """Numerical discretisation and convergence controls in reduced units.

    The default spatial interval reproduces the legacy grid construction:
    :math:`x \\in [-\\pi/2,\\, 3\\pi/2]` with 600 points.  Convergence tolerance
    is the dimensionless threshold used for thermal-sum stopping criteria.
    """

    #: Total number of grid points.  Dimensionless count.  Default: ``600``.
    n_points: int = 600
    #: Left boundary of the spatial domain in dimensionless length units.
    #: Default: :math:`-\\pi/2`.
    x_min: float = -math.pi / 2.0
    #: Right boundary of the spatial domain in dimensionless length units.
    #: Default: :math:`3\\pi/2`.
    x_max: float = 3.0 * math.pi / 2.0
    #: Maximum absolute change allowed in consecutive thermal-density estimates.
    #: Dimensionless probability-density units.  Default: :math:`10^{-4}`.
    convergence_tolerance: float = 1e-4

    def __post_init__(self) -> None:
        """Validates numerical parameter consistency.

        Raises:
            ValueError: If `n_points < 2`, if `x_max <= x_min`, or if
                `convergence_tolerance <= 0`.
        """
        if self.n_points < 2:
            raise ValueError("n_points must be at least 2.")
        if self.x_max <= self.x_min:
            raise ValueError("x_max must be greater than x_min.")
        if self.convergence_tolerance <= 0.0:
            raise ValueError("convergence_tolerance must be positive.")
