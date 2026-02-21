"""Abstract base class and concrete implementations for 1D quantum potentials.

This module defines the ``BasePotential`` interface and four concrete
subclasses that faithfully reproduce the potential functions found in
``legacy/research_prototype.ipynb``:

- :class:`InfiniteSquareWell` — Analytical infinite-barrier box.
- :class:`FiniteSquareWell`   — Flat-bottom well with finite barriers.
- :class:`HarmonicPotential`  — Truncated harmonic oscillator (parabolic).
- :class:`VShapedPotential`   — Truncated V-shape (linear) potential.

All subclasses accept a :class:`~src.models.config.PhysicsConfig` object
for their physical constants, enforcing Dependency Injection and eliminating
global-constant coupling.
"""

from __future__ import annotations

import math
from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import NDArray

from src.models.config import PhysicsConfig


# ---------------------------------------------------------------------------
# Abstract Base Class
# ---------------------------------------------------------------------------


class BasePotential(ABC):
    """Abstract interface for all 1D quantum potential definitions.

    Every potential that participates in the Numerov eigenvalue pipeline must
    implement this contract.  The :class:`~src.core.numerov_solver.NumerovSolver`
    depends *only* on this interface, so new geometries (e.g. a double well)
    can be added in a new subclass without ever modifying the solver — in
    strict compliance with the Open-Closed Principle.

    Args:
        config: Frozen dataclass holding the dimensionless physical constants
            (``hbar``, ``mass``, ``kb``, ``well_width``) used throughout the
            simulation.  Defaults to the standard reduced-unit set
            (``hbar=1``, ``mass=0.5``, ``L=pi``) that makes the infinite-well
            spectrum satisfy :math:`E_n = n^2`.
    """

    def __init__(self, config: PhysicsConfig | None = None) -> None:
        """Stores the injected physics configuration.

        Args:
            config: Physics constants dataclass.  When ``None``, the standard
                dimensionless defaults are used automatically.
        """
        self._config: PhysicsConfig = config if config is not None else PhysicsConfig()

    # ------------------------------------------------------------------
    # Public read-only property
    # ------------------------------------------------------------------

    @property
    def config(self) -> PhysicsConfig:
        """PhysicsConfig: The injected physics constants, read-only."""
        return self._config

    # ------------------------------------------------------------------
    # Abstract contract
    # ------------------------------------------------------------------

    @abstractmethod
    def evaluate(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """Evaluates the potential energy :math:`V(x)` at each grid point.

        This is the primary contract method consumed by the Numerov solver.
        All values are returned in dimensionless energy units of :math:`E_1`,
        where :math:`E_1 = \\hbar^2 / (2m)` with the defaults ``hbar=1``,
        ``mass=0.5``.

        Args:
            x: 1D array of spatial coordinates in dimensionless length units
                (e.g. :math:`x/L` or literal position when :math:`L = \\pi`).

        Returns:
            1D float64 array with the same shape as ``x`` containing
            :math:`V(x)` values.  Barrier regions may contain ``np.inf``
            for infinite-wall geometries.
        """
        ...

    @abstractmethod
    def get_domain(self) -> tuple[float, float]:
        """Returns the recommended spatial domain :math:`[x_{min}, x_{max}]`.

        The Numerov solver calls this method to build the spatial grid without
        any hardcoded geometry knowledge.  The domain should comfortably
        contain both the classically allowed region and (for finite barriers)
        the evanescent tunnelling tails.

        Returns:
            A ``(x_min, x_max)`` tuple in dimensionless length units.

        Example:
            A finite well of width :math:`L = \\pi` returns::

                >>> well = FiniteSquareWell(v0=10.0, width=math.pi)
                >>> well.get_domain()
                (-1.5707963267948966, 4.71238898038469)
        """
        ...

    # ------------------------------------------------------------------
    # Concrete helpers
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        """Returns an unambiguous string representation.

        Returns:
            A string in the format ``ClassName(config=<PhysicsConfig>)``.
        """
        return f"{self.__class__.__name__}(config={self._config!r})"


# ---------------------------------------------------------------------------
# Concrete Implementations
# ---------------------------------------------------------------------------


class InfiniteSquareWell(BasePotential):
    """Infinite square-well potential with exact Dirichlet boundary conditions.

    The potential is defined as:

    .. math::

        V(x) =
        \\begin{cases}
            0       & \\text{if } 0 \\le x \\le L \\\\
            +\\infty & \\text{otherwise}
        \\end{cases}

    where :math:`L` = ``config.well_width`` (default :math:`\\pi`).

    This geometry has an exact analytical eigenspectrum given by
    :math:`E_n = n^2` in the dimensionless unit system (``hbar=1``,
    ``mass=0.5``, :math:`L=\\pi`), which is used as the reference baseline
    throughout the project.

    Args:
        config: Physics constants.  Reads ``well_width`` to locate the right
            wall of the box.

    Example:
        >>> well = InfiniteSquareWell()
        >>> import numpy as np
        >>> x = np.array([-0.5, 0.5, np.pi / 2, np.pi, np.pi + 0.1])
        >>> well.evaluate(x)
        array([inf,  0.,  0.,  0., inf])
    """

    def evaluate(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """Returns :math:`V(x)` with :math:`+\\infty` walls at :math:`x=0` and :math:`x=L`.

        Args:
            x: 1D array of spatial coordinates in dimensionless length units.

        Returns:
            Float64 array of the same shape as ``x``.  Values outside
            :math:`[0, L]` are set to ``np.inf``; values inside are ``0.0``.
        """
        L: float = self._config.well_width
        V: NDArray[np.float64] = np.full_like(x, np.inf, dtype=np.float64)
        inside: NDArray[np.bool_] = (x >= 0.0) & (x <= L)
        V[inside] = 0.0
        return V

    def get_domain(self) -> tuple[float, float]:
        """Returns :math:`[0, L]` — the exact physical domain of the box.

        The solver requires Dirichlet conditions at both ends, so no external
        margin is needed.

        Returns:
            The tuple ``(0.0, well_width)``.
        """
        return (0.0, self._config.well_width)


class FiniteSquareWell(BasePotential):
    """Finite square-well potential with flat barriers.

    Faithfully reproduces ``get_potential(x, V0, well_width)`` from
    ``legacy/research_prototype.ipynb``.

    The potential is defined as:

    .. math::

        V(x) =
        \\begin{cases}
            0    & \\text{if } 0 \\le x \\le w \\\\
            V_0  & \\text{otherwise}
        \\end{cases}

    where :math:`w` is the ``width`` parameter and :math:`V_0` is the
    barrier height.

    The returned domain adds a half-width margin on the left and a
    full-width margin on the right — identical to the legacy constants
    ``X_MIN = -L/2``, ``X_MAX = 3L/2`` — to capture the evanescent
    tunnelling tails that extend into the barrier region.

    Args:
        v0: Barrier height in dimensionless energy units of :math:`E_1`.
            The legacy study uses ``v0 ∈ {10, 50, 1000}``.
        width: Well width in dimensionless length units.  Defaults to
            ``config.well_width`` (i.e. :math:`\\pi`) when ``None``.
        config: Physics constants.  Reads ``well_width`` as the fallback
            when ``width`` is not provided.

    Raises:
        ValueError: If ``v0 < 0`` or ``width <= 0``.

    Example:
        >>> import numpy as np
        >>> well = FiniteSquareWell(v0=10.0)
        >>> x = np.array([-0.5, 0.5, np.pi / 2, np.pi, np.pi + 0.5])
        >>> well.evaluate(x)
        array([10. ,  0. ,  0. ,  0. , 10. ])
    """

    def __init__(
        self,
        v0: float,
        width: float | None = None,
        config: PhysicsConfig | None = None,
    ) -> None:
        """Initialises the finite square well.

        Args:
            v0: Barrier height :math:`V_0 \\ge 0`.
            width: Well width.  Defaults to ``config.well_width`` when
                ``None``.
            config: Physics constants dataclass.

        Raises:
            ValueError: If ``v0 < 0`` or the resolved ``width <= 0``.
        """
        super().__init__(config)
        if v0 < 0.0:
            raise ValueError(f"Barrier height v0 must be non-negative, got {v0}.")
        resolved_width: float = width if width is not None else self._config.well_width
        if resolved_width <= 0.0:
            raise ValueError(f"Well width must be positive, got {resolved_width}.")
        self._v0: float = v0
        self._width: float = resolved_width

    @property
    def v0(self) -> float:
        """float: Barrier height in dimensionless energy units."""
        return self._v0

    @property
    def width(self) -> float:
        """float: Well width in dimensionless length units."""
        return self._width

    def evaluate(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """Returns :math:`V(x) = 0` inside :math:`[0, w]`, :math:`V_0` outside.

        Args:
            x: 1D array of spatial coordinates in dimensionless length units.

        Returns:
            Float64 array of the same shape as ``x``.
        """
        V: NDArray[np.float64] = np.full(x.shape, self._v0, dtype=np.float64)
        inside: NDArray[np.bool_] = (x >= 0.0) & (x <= self._width)
        V[inside] = 0.0
        return V

    def get_domain(self) -> tuple[float, float]:
        """Returns the extended domain :math:`[-w/2,\\, 3w/2]`.

        The half-width left margin and full-width right margin replicate the
        legacy grid ``X_MIN = -L/2``, ``X_MAX = 3L/2``, ensuring sufficient
        room for the evanescent wave tails to decay to zero.

        Returns:
            The tuple ``(-width/2, 3*width/2)``.
        """
        return (-self._width / 2.0, 3.0 * self._width / 2.0)

    def __repr__(self) -> str:
        return (
            f"FiniteSquareWell(v0={self._v0}, width={self._width}, "
            f"config={self._config!r})"
        )


class HarmonicPotential(BasePotential):
    """Truncated harmonic oscillator (parabolic) potential.

    Faithfully reproduces ``potential_harmonic_finite(x, V0, omega, center)``
    from ``legacy/research_prototype.ipynb``.

    The potential is defined as:

    .. math::

        V(x) = \\min\\!\\left(\\frac{1}{2} m\\,\\omega^2 (x - x_c)^2,\\; V_{\\rm lim}\\right)

    where :math:`m` = ``config.mass``, :math:`\\omega` is the angular
    frequency, :math:`x_c` is the well centre, and :math:`V_{\\rm lim}` is
    the truncation ceiling.  Truncation ensures finite barrier heights for
    numerical tractability without unphysical divergence.

    Args:
        omega: Angular frequency :math:`\\omega` in dimensionless units.
            Legacy default is ``2.2``.
        v0_limit: Truncation ceiling :math:`V_{\\rm lim}` in dimensionless
            energy units.  Legacy default is ``60.0``.
        center: Position of the parabolic minimum :math:`x_c` in
            dimensionless length units.  Legacy default is :math:`\\pi/2`.
        config: Physics constants.  Reads ``mass`` for the prefactor
            :math:`\\frac{1}{2} m \\omega^2`.

    Raises:
        ValueError: If ``omega <= 0`` or ``v0_limit <= 0``.

    Example:
        >>> import numpy as np, math
        >>> pot = HarmonicPotential(omega=2.2, v0_limit=60.0, center=math.pi/2)
        >>> pot.evaluate(np.array([math.pi / 2]))   # minimum at centre
        array([0.])
    """

    def __init__(
        self,
        omega: float = 2.2,
        v0_limit: float = 60.0,
        center: float = math.pi / 2.0,
        config: PhysicsConfig | None = None,
    ) -> None:
        """Initialises the truncated harmonic potential.

        Args:
            omega: Angular frequency :math:`\\omega > 0`.
            v0_limit: Truncation ceiling :math:`V_{\\rm lim} > 0`.
            center: Parabolic minimum position :math:`x_c`.
            config: Physics constants dataclass.

        Raises:
            ValueError: If ``omega <= 0`` or ``v0_limit <= 0``.
        """
        super().__init__(config)
        if omega <= 0.0:
            raise ValueError(f"omega must be positive, got {omega}.")
        if v0_limit <= 0.0:
            raise ValueError(f"v0_limit must be positive, got {v0_limit}.")
        self._omega: float = omega
        self._v0_limit: float = v0_limit
        self._center: float = center

    @property
    def omega(self) -> float:
        """float: Angular frequency in dimensionless units."""
        return self._omega

    @property
    def v0_limit(self) -> float:
        """float: Truncation ceiling in dimensionless energy units."""
        return self._v0_limit

    @property
    def center(self) -> float:
        """float: Position of the parabolic minimum in dimensionless length units."""
        return self._center

    def evaluate(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """Returns :math:`V(x) = \\min(\\tfrac{1}{2}m\\omega^2(x-x_c)^2,\\, V_{\\rm lim})`.

        Args:
            x: 1D array of spatial coordinates in dimensionless length units.

        Returns:
            Float64 array of the same shape as ``x``.
        """
        m: float = self._config.mass
        V_soft: NDArray[np.float64] = (
            0.5 * m * self._omega**2 * (x - self._center) ** 2
        )
        return np.minimum(V_soft, self._v0_limit)

    def get_domain(self) -> tuple[float, float]:
        """Returns a domain that fully contains both well and barrier regions.

        The truncation radius :math:`r_t = \\sqrt{2 V_{\\rm lim} / (m\\omega^2)}`
        is the distance from the centre at which the parabola reaches
        :math:`V_{\\rm lim}`.  A 50 % safety margin is added to ensure
        the evanescent tail decays to zero within the grid.

        Returns:
            The tuple ``(center - 1.5*r_t, center + 1.5*r_t)``.
        """
        m: float = self._config.mass
        r_trunc: float = math.sqrt(2.0 * self._v0_limit / (m * self._omega**2))
        margin: float = 1.5 * r_trunc
        return (self._center - margin, self._center + margin)

    def __repr__(self) -> str:
        return (
            f"HarmonicPotential(omega={self._omega}, v0_limit={self._v0_limit}, "
            f"center={self._center}, config={self._config!r})"
        )


class VShapedPotential(BasePotential):
    """Truncated V-shape (linear) potential.

    Faithfully reproduces ``potential_vshape_finite(x, V0, slope, center)``
    from ``legacy/research_prototype.ipynb``.

    The potential is defined as:

    .. math::

        V(x) = \\min\\!\\left(s\\,|x - x_c|,\\; V_{\\rm lim}\\right)

    where :math:`s` is the linear slope, :math:`x_c` is the cusp centre,
    and :math:`V_{\\rm lim}` is the truncation ceiling.  The result is a
    tent-shaped well with flat barriers once the slope reaches
    :math:`V_{\\rm lim}`.

    Args:
        slope: Linear slope :math:`s` in dimensionless energy-per-length
            units.  Legacy default is ``28.0``.
        v0_limit: Truncation ceiling :math:`V_{\\rm lim}` in dimensionless
            energy units.  Legacy default is ``60.0``.
        center: Position of the V-cusp :math:`x_c` in dimensionless length
            units.  Legacy default is :math:`\\pi/2`.
        config: Physics constants dataclass.  Not required to evaluate this
            potential (no mass dependency), but accepted for interface
            consistency.

    Raises:
        ValueError: If ``slope <= 0`` or ``v0_limit <= 0``.

    Example:
        >>> import numpy as np, math
        >>> pot = VShapedPotential(slope=28.0, v0_limit=60.0, center=math.pi / 2)
        >>> pot.evaluate(np.array([math.pi / 2]))   # cusp — minimum is zero
        array([0.])
    """

    def __init__(
        self,
        slope: float = 28.0,
        v0_limit: float = 60.0,
        center: float = math.pi / 2.0,
        config: PhysicsConfig | None = None,
    ) -> None:
        """Initialises the truncated V-shaped potential.

        Args:
            slope: Linear slope :math:`s > 0`.
            v0_limit: Truncation ceiling :math:`V_{\\rm lim} > 0`.
            center: Cusp position :math:`x_c`.
            config: Physics constants dataclass.

        Raises:
            ValueError: If ``slope <= 0`` or ``v0_limit <= 0``.
        """
        super().__init__(config)
        if slope <= 0.0:
            raise ValueError(f"slope must be positive, got {slope}.")
        if v0_limit <= 0.0:
            raise ValueError(f"v0_limit must be positive, got {v0_limit}.")
        self._slope: float = slope
        self._v0_limit: float = v0_limit
        self._center: float = center

    @property
    def slope(self) -> float:
        """float: Linear slope in dimensionless energy-per-length units."""
        return self._slope

    @property
    def v0_limit(self) -> float:
        """float: Truncation ceiling in dimensionless energy units."""
        return self._v0_limit

    @property
    def center(self) -> float:
        """float: Cusp position in dimensionless length units."""
        return self._center

    def evaluate(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        """Returns :math:`V(x) = \\min(s\\,|x - x_c|,\\, V_{\\rm lim})`.

        Args:
            x: 1D array of spatial coordinates in dimensionless length units.

        Returns:
            Float64 array of the same shape as ``x``.
        """
        V_soft: NDArray[np.float64] = self._slope * np.abs(x - self._center)
        return np.minimum(V_soft, self._v0_limit)

    def get_domain(self) -> tuple[float, float]:
        """Returns a domain that fully contains the well and flat barriers.

        The truncation radius :math:`r_t = V_{\\rm lim} / s` is the distance
        from the cusp at which the slope reaches :math:`V_{\\rm lim}`.  A
        50 % safety margin is added for the evanescent tails.

        Returns:
            The tuple ``(center - 1.5*r_t, center + 1.5*r_t)``.
        """
        r_trunc: float = self._v0_limit / self._slope
        margin: float = 1.5 * r_trunc
        return (self._center - margin, self._center + margin)

    def __repr__(self) -> str:
        return (
            f"VShapedPotential(slope={self._slope}, v0_limit={self._v0_limit}, "
            f"center={self._center}, config={self._config!r})"
        )
