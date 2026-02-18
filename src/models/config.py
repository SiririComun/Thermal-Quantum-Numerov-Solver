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
    legacy/research_prototype.ipynb so that the 1D infinite square-well
    reference has energies $E_n = n^2$.

    Attributes:
        hbar: Reduced Planck constant in dimensionless units. Default is 1.0.
        mass: Particle mass in dimensionless units. Default is 0.5, which
            enforces $\\hbar^2/(2m)=1$ when `hbar=1`.
        kb: Boltzmann constant in dimensionless thermal units. Default is 1.0,
            so reduced temperature is directly represented as $t = k_B T$.
        well_width: Characteristic well width $L$ in dimensionless length units.
            Default is $\\pi$.
    """

    hbar: float = 1.0
    mass: float = 0.5
    kb: float = 1.0
    well_width: float = math.pi


@dataclass(frozen=True, slots=True)
class NumericalConfig:
    """Numerical discretization and convergence controls in reduced units.

    The default spatial interval reproduces the legacy grid construction:
    $x \\in [-\\pi/2,\\, 3\\pi/2]$ with 600 points. Convergence tolerance is the
    dimensionless threshold used for thermal-sum stopping criteria.

    Attributes:
        n_points: Total number of grid points for the spatial discretization.
            Dimensionless count (unitless). Default is 600.
        x_min: Left boundary of the spatial domain in dimensionless length
            units. Default is $-\\pi/2$.
        x_max: Right boundary of the spatial domain in dimensionless length
            units. Default is $3\\pi/2$.
        convergence_tolerance: Maximum allowed absolute change in consecutive
            thermal-density estimates, in dimensionless probability-density
            units. Default is $10^{-4}$.
    """

    n_points: int = 600
    x_min: float = -math.pi / 2.0
    x_max: float = 3.0 * math.pi / 2.0
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
