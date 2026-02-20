"""Thermal statistics engine for 1D quantum systems of identical particles.

This module provides :class:`ThermalEngine`, the service class responsible for
all statistical-mechanics post-processing of a solved :class:`QuantumSystem`.
It refactors and encapsulates two families of functions from
``legacy/research_prototype.ipynb``:

- ``calc_thermal_density`` (single-particle Boltzmann-weighted density)
- ``get_pair_prob`` + ``calc_thermal_2p`` (two-particle symmetrised density)

Design (SOLID)
--------------
- **SRP**: ``ThermalEngine`` does one thing — thermal averaging.  It never
  touches eigenvalue solvers or potential geometry.
- **OCP**: New statistics (e.g. anyons) can be added by extending
  :class:`~src.models.statistics.ParticleType` and the ``_spin_components``
  helper without modifying existing methods.
- **DIP**: All physics input arrives through the :class:`~src.models.states.QuantumSystem`
  interface; no concrete solver class is imported.

Mathematical Summary
--------------------
**Single-particle thermal density**

.. math::

    \\rho_{\\text{th}}(x) = \\frac{1}{Z}
    \\sum_{n=0}^{N-1} e^{-(E_n - E_0)/(k_B T)}\\, |\\psi_n(x)|^2,
    \\qquad Z = \\sum_{n=0}^{N-1} e^{-(E_n - E_0)/(k_B T)}

States are accumulated one at a time; summation stops when the maximum
pointwise change in :math:`\\rho(x)` falls below ``tolerance``.

**Two-particle symmetrised pair density** (for :math:`n \\neq m`)

.. math::

    \\Psi_S(x_1, x_2) = \\frac{1}{\\sqrt{2}}
    \\bigl[\\phi_n(x_1)\\phi_m(x_2) + \\phi_n(x_2)\\phi_m(x_1)\\bigr]
    \\quad (\\text{bosons})

.. math::

    \\Psi_A(x_1, x_2) = \\frac{1}{\\sqrt{2}}
    \\bigl[\\phi_n(x_1)\\phi_m(x_2) - \\phi_n(x_2)\\phi_m(x_1)\\bigr]
    \\quad (\\text{fermions})

The diagonal :math:`\\Psi_A(x, x) = 0` identically (Pauli Exclusion Principle).

**Thermal two-particle pair density**

.. math::

    \\rho_{\\text{pair}}(x_1, x_2) = \\frac{1}{Z}
    \\sum_{n,m} e^{-(E_n + E_m - 2E_0)/(k_B T)}\\,
    \\sum_c f_c\\, |\\Psi_{c,nm}(x_1, x_2)|^2

where :math:`c` indexes symmetry channels and :math:`f_c` are the
spin-degeneracy fractions from :class:`~src.models.statistics.ParticleType`.
"""

from __future__ import annotations

import logging
from math import sqrt
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from src.models.states import QuantumSystem
from src.models.statistics import ParticleType

if TYPE_CHECKING:
    pass


# ---------------------------------------------------------------------------
# Spin-fraction lookup (pure data — no logic)
# ---------------------------------------------------------------------------

# Mapping from ParticleType to [(component_ParticleType, fraction), ...]
# BOLTZMANN is handled as a special branch (no symmetrisation).
# BOSON and FERMION are "pure" — only one component with fraction 1.0.
# SPIN_HALF and SPIN_ONE are spin-averaged mixtures.
#
# Derivations: see module docstring and src/models/statistics.py.
_SPIN_COMPONENTS: dict[ParticleType, list[tuple[ParticleType, float]]] = {
    ParticleType.BOSON: [(ParticleType.BOSON, 1.0)],
    ParticleType.FERMION: [(ParticleType.FERMION, 1.0)],
    # Spin-1/2 fermion pair: singlet (1/4) symmetric + triplet (3/4) antisymmetric
    ParticleType.SPIN_HALF: [(ParticleType.BOSON, 1.0 / 4.0),
                             (ParticleType.FERMION, 3.0 / 4.0)],
    # Spin-1 boson pair: (S=2,S=0) symmetric (6/9) + S=1 antisymmetric (3/9)
    ParticleType.SPIN_ONE: [(ParticleType.BOSON, 2.0 / 3.0),
                            (ParticleType.FERMION, 1.0 / 3.0)],
}

_LOG = logging.getLogger(__name__)


class ThermalEngine:
    """Service class for thermal and multi-particle statistical calculations.

    :class:`ThermalEngine` is a stateless service: all state lives in the
    :class:`~src.models.states.QuantumSystem` objects passed to its methods.
    The engine is safe to instantiate once and reuse across multiple systems
    and temperatures.

    Args:
        log_level: Python logging level (e.g. ``logging.DEBUG``) used for
            internal convergence progress messages.  Defaults to
            ``logging.WARNING`` (silent in normal use).

    Example:
        >>> from src.core import NumerovSolver
        >>> from src.core.statistics import ThermalEngine
        >>> from src.models import InfiniteSquareWell
        >>> from src.models.statistics import ParticleType
        >>>
        >>> system = NumerovSolver().solve(InfiniteSquareWell(), n_states=20)
        >>> engine = ThermalEngine()
        >>> rho    = engine.calculate_thermal_density(system, temperature=5.0)
        >>> rho_2p = engine.calculate_pair_density(system, 5.0, ParticleType.FERMION)
    """

    def __init__(self, log_level: int = logging.WARNING) -> None:
        """Initialises the engine and configures the internal logger.

        Args:
            log_level: Logging verbosity.  Set to ``logging.DEBUG`` to trace
                convergence iteration counts.
        """
        _LOG.setLevel(log_level)
        _LOG.debug("ThermalEngine initialised.")

    # ------------------------------------------------------------------
    # Public methods
    # ------------------------------------------------------------------

    def calculate_thermal_density(
        self,
        system: QuantumSystem,
        temperature: float,
        tolerance: float = 1e-4,
    ) -> NDArray[np.float64]:
        """Computes the single-particle thermal probability density :math:`\\rho_{\\text{th}}(x)`.

        States are added one at a time in ascending energy order.  Summation
        halts as soon as the maximum pointwise change in :math:`\\rho(x)` after
        adding the latest state is smaller than ``tolerance`` — the convergence
        criterion from ``legacy/research_prototype.ipynb``.

        The Boltzmann factor uses :math:`E_0` as the energy reference to
        prevent floating-point underflow for large :math:`E_n / (k_B T)`:

        .. math::

            w_n = e^{-(E_n - E_0)/(k_B T)}

        Args:
            system: A solved :class:`~src.models.states.QuantumSystem`.
            temperature: Absolute temperature :math:`T` in dimensionless units
                (same scale as ``system.potential.config.kb``).
            tolerance: Maximum allowed pointwise change in :math:`\\rho(x)`
                between successive state additions.  Summation stops when
                this threshold is met.  Defaults to ``1e-4``.

        Returns:
            1D float64 array of shape ``(n_points,)`` containing
            :math:`\\rho_{\\text{th}}(x)`, normalised so that its trapezoid
            integral over ``system.x_grid`` is :math:`\\approx 1`.

        Raises:
            ValueError: If ``temperature <= 0``.

        Notes:
            At very low temperatures, only the ground state contributes
            meaningfully.  At high temperatures, the sum may require all
            ``system.n_states`` states before convergence; in that case a
            ``logging.DEBUG`` message is emitted and all states are used.

        Example:
            >>> rho = engine.calculate_thermal_density(system, temperature=10.0)
            >>> float(np.trapz(rho, system.x_grid))  # ≈ 1.0
            1.0
        """
        if temperature <= 0.0:
            raise ValueError(
                f"temperature must be strictly positive, got {temperature}."
            )

        # k_B is stored in the potential's PhysicsConfig (kb=1.0 in default units)
        kb: float = system.potential.config.kb
        kbt: float = kb * temperature       # k_B · T  (dimensionless thermal energy)
        E0: float = system.ground_state_energy

        # ── Accumulate states with convergence check ─────────────────────
        rho_prev: NDArray[np.float64] = np.zeros(len(system.x_grid), dtype=np.float64)
        rho_accum: NDArray[np.float64] = np.zeros(len(system.x_grid), dtype=np.float64)
        Z: float = 0.0
        n_converged: int = system.n_states  # fallback: use all states

        for n in range(system.n_states):
            E_n: float = float(system.energies[n])
            # Boltzmann weight: reference-shifted to avoid underflow
            w_n: float = float(np.exp(-(E_n - E0) / kbt))
            Z += w_n
            # Accumulate weighted probability density for state n
            rho_accum += w_n * system.wavefunctions[:, n] ** 2

            rho_current: NDArray[np.float64] = rho_accum / Z

            # Convergence check: skip the first state (no previous estimate)
            if n > 0:
                max_change: float = float(np.max(np.abs(rho_current - rho_prev)))
                if max_change < tolerance:
                    n_converged = n + 1
                    _LOG.debug(
                        "calculate_thermal_density: converged after %d states "
                        "(max_change=%.2e < tol=%.2e, T=%.4f).",
                        n_converged, max_change, tolerance, temperature,
                    )
                    return rho_current

            rho_prev = rho_current.copy()

        _LOG.debug(
            "calculate_thermal_density: used all %d states without meeting "
            "tolerance %.2e (T=%.4f).  Consider increasing n_states.",
            system.n_states, tolerance, temperature,
        )
        return rho_prev

    def get_pair_density(
        self,
        system: QuantumSystem,
        n1: int,
        n2: int,
        p_type: ParticleType,
    ) -> NDArray[np.float64]:
        """Returns the two-particle spatial density :math:`|\\Psi(x_1, x_2)|^2` for a single eigenstate pair.

        This is the *zero-temperature* building block used by
        :meth:`calculate_pair_density`.  Calling this method directly is
        useful for visualising a specific quantum state pair or verifying
        the Pauli Exclusion Principle on the diagonal.

        **Symmetry Requirement (Quantum Mechanics):**
        For *identical* particles, the total wavefunction must be (anti-)symmetric
        under particle exchange (:math:`x_1 \\leftrightarrow x_2`).  The spatial
        part alone determines the exchange character only when the spin state is
        fixed.  Here, the spatial symmetry is enforced per the ``p_type``
        specification:

        - **BOSON** (:math:`n_1 \\neq n_2`):
          :math:`\\Psi_S = (\\phi_{n_1}(x_1)\\phi_{n_2}(x_2) + \\phi_{n_1}(x_2)\\phi_{n_2}(x_1)) / \\sqrt{2}`
        - **FERMION** (:math:`n_1 \\neq n_2`):
          :math:`\\Psi_A = (\\phi_{n_1}(x_1)\\phi_{n_2}(x_2) - \\phi_{n_1}(x_2)\\phi_{n_2}(x_1)) / \\sqrt{2}`
        - **FERMION** (:math:`n_1 = n_2`):
          :math:`\\Psi_A \\equiv 0` — Pauli Exclusion Principle.
        - **BOSON** (:math:`n_1 = n_2`):
          :math:`\\Psi_S = \\phi_{n_1}(x_1)\\phi_{n_1}(x_2)` (already normalised).

        For :data:`~src.models.statistics.ParticleType.SPIN_HALF` and
        :data:`~src.models.statistics.ParticleType.SPIN_ONE`, the density is
        a weighted sum of the symmetric and antisymmetric components.

        Args:
            system: A solved :class:`~src.models.states.QuantumSystem`.
            n1: Zero-based index of the first eigenstate (must be in
                ``[0, system.n_states)``).
            n2: Zero-based index of the second eigenstate (must be in
                ``[0, system.n_states)``).
            p_type: The :class:`~src.models.statistics.ParticleType` controlling
                the spatial symmetry of the pair wavefunction.

        Returns:
            2D float64 array of shape ``(n_points, n_points)`` where entry
            ``[i, j]`` is :math:`|\\Psi(x_i, x_j)|^2`.

        Raises:
            IndexError: If ``n1`` or ``n2`` is outside ``[0, system.n_states)``.
        """
        if not (0 <= n1 < system.n_states):
            raise IndexError(
                f"n1={n1} is out of range for system with {system.n_states} states."
            )
        if not (0 <= n2 < system.n_states):
            raise IndexError(
                f"n2={n2} is out of range for system with {system.n_states} states."
            )

        psi_n: NDArray[np.float64] = system.wavefunctions[:, n1]
        psi_m: NDArray[np.float64] = system.wavefunctions[:, n2]

        # Outer product: base[i, j] = ψ_n(x_i) · ψ_m(x_j)
        # Shape: (n_points, n_points) — the core of all pair calculations.
        # np.outer is used because it is a vectorized O(N²) BLAS-level operation,
        # far more efficient than any Python loop for constructing the 2D grid.
        base: NDArray[np.float64] = np.outer(psi_n, psi_m)

        return self._compute_pair_density(base, n1, n2, p_type)

    def calculate_pair_density(
        self,
        system: QuantumSystem,
        temperature: float,
        p_type: ParticleType,
    ) -> NDArray[np.float64]:
        """Computes the Boltzmann-weighted two-particle spatial density matrix.

        Sums over all pairs :math:`(n, m)` of eigenstates, weighted by the
        two-particle Boltzmann factor, and applies the symmetrisation required
        by ``p_type``:

        .. math::

            \\rho_{\\text{pair}}(x_1, x_2) = \\frac{1}{Z}
            \\sum_{n=0}^{N-1} \\sum_{m=0}^{N-1}
            e^{-(E_n + E_m - 2 E_0)/(k_B T)}
            \\sum_c f_c\\, |\\Psi_{c,nm}(x_1, x_2)|^2

        The energy reference :math:`2 E_0` (twice the single-particle ground
        state) is subtracted to prevent floating-point underflow.

        Args:
            system: A solved :class:`~src.models.states.QuantumSystem`.
            temperature: Absolute temperature :math:`T` in dimensionless units.
            p_type: The :class:`~src.models.statistics.ParticleType` controlling
                symmetry (and spin-mixture fractions) of the pair wavefunction.

        Returns:
            2D float64 array of shape ``(n_points, n_points)`` where entry
            ``[i, j]`` is the thermal pair probability density at
            :math:`(x_i, x_j)`.

        Raises:
            ValueError: If ``temperature <= 0``.

        Notes:
            **Performance**: The inner loop constructs one ``(N, N)`` outer
            product per state pair using :func:`numpy.outer` (BLAS level-2).
            For ``n_states=20`` and ``N=600``, this is 400 outer products of
            size :math:`600 \\times 600 = 360{,}000` elements each — fully
            vectorised and cache-friendly.

            **Fermionic normalisation**: For :data:`ParticleType.FERMION`, pairs
            with :math:`n_1 = n_2` contribute zero probability density (Pauli
            Exclusion) but their Boltzmann weight still enters :math:`Z`.  This
            mirrors the legacy implementation and preserves relative comparison
            between boson and fermion distributions.

        Example:
            >>> rho_2p = engine.calculate_pair_density(system, 10.0, ParticleType.BOSON)
            >>> rho_2p.shape
            (600, 600)
        """
        if temperature <= 0.0:
            raise ValueError(
                f"temperature must be strictly positive, got {temperature}."
            )

        kb: float = system.potential.config.kb
        kbt: float = kb * temperature
        E0: float = system.ground_state_energy
        E: NDArray[np.float64] = system.energies   # shape (n_states,)
        n_pts: int = len(system.x_grid)
        n_states: int = system.n_states

        P_accum: NDArray[np.float64] = np.zeros((n_pts, n_pts), dtype=np.float64)
        Z: float = 0.0

        for n in range(n_states):
            psi_n: NDArray[np.float64] = system.wavefunctions[:, n]
            for m in range(n_states):
                psi_m: NDArray[np.float64] = system.wavefunctions[:, m]

                # Two-particle Boltzmann weight.
                # Reference: 2·E_0 (total ground-state energy of the pair).
                # This shift is mathematically exact and prevents underflow for
                # high-lying states at low temperatures.
                w_nm: float = float(
                    np.exp(-((E[n] + E[m] - 2.0 * E0) / kbt))
                )
                Z += w_nm

                # Outer product: base[i,j] = ψ_n(x_i)·ψ_m(x_j)
                base: NDArray[np.float64] = np.outer(psi_n, psi_m)

                # Accumulate symmetrised pair density weighted by spin fractions
                P_accum += w_nm * self._compute_pair_density(base, n, m, p_type)

        return P_accum / Z

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _compute_pair_density(
        self,
        base: NDArray[np.float64],
        n1: int,
        n2: int,
        p_type: ParticleType,
    ) -> NDArray[np.float64]:
        """Applies the symmetry requirement to a pre-computed outer product.

        This is the core computational primitive of the symmetrisation algebra.
        It is separated from the public methods so that the outer product (the
        expensive operation) is computed only once per pair and reused across
        the spin-mixture components.

        **Vectorised symmetrisation** (``np.outer`` is a rank-1 update):

        Given ``base[i,j] = ψ_n(x_i) · ψ_m(x_j)``:

        - Symmetric:     ``(base + base.T) / √2``
          → evaluates to zero on the anti-diagonal, enhanced on the diagonal.
        - Antisymmetric: ``(base - base.T) / √2``
          → ``base[i,i] - base.T[i,i] = ψ_n(x_i)ψ_m(x_i) - ψ_m(x_i)ψ_n(x_i) = 0``
          → **exact** zero on the diagonal for all grids, all potentials.

        Args:
            base: 2D float64 array of shape ``(n_points, n_points)`` where
                ``base[i, j] = ψ_n(x_i) · ψ_m(x_j)``.
            n1: Zero-based index of the first eigenstate (used to detect the
                ``n1 == n2`` degenerate case).
            n2: Zero-based index of the second eigenstate.
            p_type: :class:`~src.models.statistics.ParticleType` controlling
                the spatial symmetry.

        Returns:
            2D float64 array of shape ``(n_points, n_points)`` containing
            :math:`|\\Psi(x_1, x_2)|^2` for the given symmetry type.
        """
        n_pts: int = base.shape[0]

        if p_type is ParticleType.BOLTZMANN:
            # Classical / distinguishable particles: no exchange correction.
            # P(x1, x2) = |ψ_n(x1)|² · |ψ_m(x2)|²  (outer product of densities)
            # Using base**2 = ψ_n(x1)²·ψ_m(x2)² is equivalent for real ψ.
            return base ** 2

        # For quantum particles: decompose into symmetry channels weighted by
        # their spin-degeneracy fractions from _SPIN_COMPONENTS.
        components = _SPIN_COMPONENTS[p_type]
        result: NDArray[np.float64] = np.zeros((n_pts, n_pts), dtype=np.float64)

        for component, fraction in components:
            if component is ParticleType.BOSON:
                if n1 == n2:
                    # Same state: Ψ(x1,x2) = ψ_n(x1)·ψ_n(x2) — no √2 factor.
                    # This is already normalised since ∫|ψ_n|² dx = 1.
                    psi_pair: NDArray[np.float64] = base
                else:
                    # ── Bosonic symmetrisation (exchange enhancement) ──────────
                    # Ψ_S = [ψ_n(x1)ψ_m(x2) + ψ_m(x1)ψ_n(x2)] / √2
                    # = (base + base.T) / √2
                    # Diagonal: ψ_n(x)ψ_m(x) + ψ_m(x)ψ_n(x) = 2ψ_n(x)ψ_m(x)
                    # → bosons BUNCH: enhanced probability near x1=x2.
                    psi_pair = (base + base.T) / sqrt(2.0)

            else:  # component is ParticleType.FERMION
                if n1 == n2:
                    # Pauli Exclusion Principle: two identical fermions cannot
                    # occupy the same single-particle state.  The antisymmetric
                    # combination vanishes identically: Ψ_A = 0.
                    result += fraction * np.zeros((n_pts, n_pts), dtype=np.float64)
                    continue
                else:
                    # ── Fermionic antisymmetrisation (exchange hole) ───────────
                    # Ψ_A = [ψ_n(x1)ψ_m(x2) - ψ_m(x1)ψ_n(x2)] / √2
                    # = (base - base.T) / √2
                    #
                    # Diagonal proof (Pauli Exclusion):
                    # (base - base.T)[i,i]
                    #   = ψ_n(x_i)·ψ_m(x_i) - ψ_m(x_i)·ψ_n(x_i)
                    #   = 0  (exactly, for any grid, any potential)
                    # ⟹ |Ψ_A(x,x)|² = 0 for all x — the exchange hole.
                    psi_pair = (base - base.T) / sqrt(2.0)

            # |Ψ|² is the spatial probability density of the pair
            result += fraction * psi_pair ** 2

        return result
