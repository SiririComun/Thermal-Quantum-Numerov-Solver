"""Matrix Numerov eigenvalue solver for the 1D Time-Independent Schrödinger Equation.

This module provides the :class:`NumerovSolver` service class, which refactors
the procedural ``solve_schrodinger_numerov`` function from
``legacy/research_prototype.ipynb`` into a production-grade engine.

Mathematical Background
-----------------------
The TISE in dimensionless units is:

.. math::

    -\\frac{\\hbar^2}{2m} \\frac{d^2\\psi}{dx^2} + V(x)\\psi = E\\psi

which, after substituting the scaled variable :math:`\\varepsilon = 2mE/\\hbar^2`
and :math:`U(x) = (2m/\\hbar^2)V(x)`, becomes:

.. math::

    \\psi'' + (\\varepsilon - U(x))\\psi = 0

**Why O(dx⁴)?  The Numerov Advantage**

A standard central-difference second derivative has :math:`O(dx^2)` local
truncation error:

.. math::

    \\psi''_i = \\frac{\\psi_{i+1} - 2\\psi_i + \\psi_{i-1}}{dx^2} + O(dx^2)

The Numerov method improves this by exploiting the structure of the
differential equation itself.  Starting from the exact Taylor expansion:

.. math::

    \\psi_{i+1} + \\psi_{i-1} = 2\\psi_i + dx^2 \\psi''_i + \\frac{dx^4}{12}
    \\psi''''_i + O(dx^6)

and noting that for the TISE :math:`\\psi'''' = (U-\\varepsilon)^2\\psi + O(dx^2)`,
the fourth-derivative term can be approximated to :math:`O(dx^2)` by a second
finite difference of :math:`(U-\\varepsilon)\\psi`:

.. math::

    \\frac{dx^4}{12}\\psi''''_i \\approx \\frac{dx^2}{12}
    \\left[f_{i+1}\\psi_{i+1} - 2f_i\\psi_i + f_{i-1}\\psi_{i-1}\\right]
    + O(dx^6)

where :math:`f_i = U_i - \\varepsilon`.  Substituting into the Taylor expansion
and rearranging gives the Numerov recurrence, which has a **global truncation
error of** :math:`O(dx^4)` — two orders better than standard FD.

**Matrix Form (Generalized Eigenvalue Problem)**

Applying the Numerov recurrence to all :math:`N-2` interior grid points
(with Dirichlet boundary conditions :math:`\\psi_0 = \\psi_{N-1} = 0`) leads to:

.. math::

    A\\,\\psi = \\varepsilon\\, B\\,\\psi

where:

.. math::

    A = \\underbrace{\\frac{-1}{dx^2} M_{-2}}_{\\text{kinetic}} +
        \\underbrace{\\frac{1}{12} M_{10}\\, \\mathrm{diag}(U)}_{\\text{potential}}
    \\qquad
    B = \\frac{1}{12} M_{10}

:math:`M_{-2}` is the tridiagonal matrix with :math:`-2` on the main diagonal
and :math:`+1` on the first off-diagonals (the standard second-difference
stencil).  :math:`M_{10}` is the tridiagonal matrix with :math:`10` on the
main diagonal and :math:`+1` on the first off-diagonals (the Numerov weight
matrix).  Both are *symmetric*.

After solving for :math:`\\varepsilon`, the physical energies are recovered as:

.. math::

    E_n = \\varepsilon_n \\cdot \\frac{\\hbar^2}{2m} = \\frac{\\varepsilon_n}{f}
    \\qquad f = \\frac{2m}{\\hbar^2}
"""

from __future__ import annotations

import numpy as np
from abc import ABC, abstractmethod
from numpy.typing import NDArray
from scipy.integrate import trapezoid
from scipy.linalg import eigh

from src.models.config import NumericalConfig, PhysicsConfig
from src.models.potentials import BasePotential
from src.models.states import QuantumSystem


# ---------------------------------------------------------------------------
# Solver Interface
# ---------------------------------------------------------------------------


class BaseSolver(ABC):
    """Abstract base class defining the interface for all 1D TISE solvers.

    Any numerical method that solves the Time-Independent Schrödinger Equation
    — be it Matrix Numerov, Finite Element, or a Shooting Method — must
    implement this contract.  High-level "Experiment" code depends *only* on
    this interface, never on a concrete solver implementation.

    This enforces the **Dependency Inversion Principle**: abstractions own the
    policy; concrete solvers own the algorithm.  Swapping from Numerov to FEM
    requires only a new subclass, with zero changes to any consumer.

    Example:
        Define and use a custom solver transparently::

            class ShootingSolver(BaseSolver):
                def solve(
                    self,
                    potential: BasePotential,
                    n_states: int = 20,
                ) -> QuantumSystem:
                    ...  # shooting-method implementation

            experiment: BaseSolver = ShootingSolver()
            system = experiment.solve(InfiniteSquareWell())
    """

    @abstractmethod
    def solve(
        self,
        potential: BasePotential,
        n_states: int = 20,
    ) -> QuantumSystem:
        """Solves the 1D TISE and returns the lowest ``n_states`` eigenpairs.

        Args:
            potential: Any concrete :class:`~src.models.potentials.BasePotential`
                subclass supplying :math:`V(x)` and the spatial domain.
            n_states: Number of lowest-energy eigenvalues and eigenvectors to
                compute.  Must be at least 1.

        Returns:
            A :class:`~src.models.states.QuantumSystem` containing the
            spatial grid, eigenvalues :math:`E_n`, normalised eigenvectors
            :math:`\\psi_n(x)`, and a reference to the originating potential.

        Raises:
            ValueError: If ``n_states < 1``.
        """
        ...


class NumerovSolver(BaseSolver):
    """Service class that solves the 1D TISE via the Matrix Numerov method.

    Encapsulates the full Numerov eigenvalue pipeline:
    (1) auto-build the spatial grid from the potential's own domain,
    (2) assemble the symmetric generalised-eigenvalue matrices A and B,
    (3) solve with ``scipy.linalg.eigh`` (optimised for real-symmetric
        positive-definite pencils),
    (4) recover physical energies and normalise wavefunctions.

    The solver is **closed for modification** and **open for extension** via
    any subclass of :class:`~src.models.potentials.BasePotential` — the
    potential geometry is injected through :meth:`solve`, never hardcoded.

    Args:
        config: Numerical discretisation parameters.  When ``None``,
            defaults are used (600 grid points, convergence tolerance 1e-4).

    Attributes:
        config: The injected :class:`~src.models.config.NumericalConfig`,
            read-only after construction.

    Example:
        >>> from src.models import InfiniteSquareWell
        >>> solver = NumerovSolver()
        >>> energies, wavefunctions = solver.solve(InfiniteSquareWell(), n_states=5)
        >>> # energies ≈ [1, 4, 9, 16, 25] in dimensionless E₁ units
    """

    def __init__(self, config: NumericalConfig | None = None) -> None:
        """Stores the injected numerical configuration.

        Args:
            config: Discretisation constants.  When ``None``, the standard
                600-point grid from :class:`~src.models.config.NumericalConfig`
                is used.
        """
        self._config: NumericalConfig = (
            config if config is not None else NumericalConfig()
        )

    # ------------------------------------------------------------------
    # Public read-only property
    # ------------------------------------------------------------------

    @property
    def config(self) -> NumericalConfig:
        """NumericalConfig: The injected numerical parameters, read-only."""
        return self._config

    # ------------------------------------------------------------------
    # Primary public method
    # ------------------------------------------------------------------

    def solve(
        self,
        potential: BasePotential,
        n_states: int = 20,
    ) -> QuantumSystem:
        """Solves the 1D TISE for the given potential and returns a QuantumSystem.

        The method follows these steps:

        1. **Auto-Grid**: Calls ``potential.get_domain()`` to obtain
           :math:`[x_{min}, x_{max}]` and builds a uniform grid of
           ``config.n_points`` points.  The domain is set by the
           *potential*, not by global constants.
        2. **Matrix assembly**: Constructs the symmetric generalised
           matrices :math:`A` and :math:`B` from the Numerov recurrence
           using only the :math:`N-2` interior grid points (Dirichlet BCs).
        3. **Solve**: Calls ``scipy.linalg.eigh(A, B, subset_by_index=...)``,
           which is optimised for symmetric-positive-definite pencils and
           returns eigenvalues in ascending order without post-sorting.
        4. **Energy recovery**: Converts dimensionless eigenvalues
           :math:`\\varepsilon_n` to physical energies
           :math:`E_n = \\varepsilon_n \\hbar^2 / (2m)`.
        5. **Normalisation**: Calls :meth:`_normalize` on each wavefunction
           so that :math:`\\int |\\psi_n(x)|^2\\, dx = 1`.
        6. **Wrapping**: Bundles the spatial grid, energies, wavefunctions,
           and originating potential into a :class:`~src.models.states.QuantumSystem`
           domain model, enforcing shape consistency and eliminating data
           desynchronisation.

        Args:
            potential: Any concrete :class:`~src.models.potentials.BasePotential`
                subclass.  The solver reads ``potential.config`` for physical
                constants and ``potential.get_domain()`` for the grid domain.
            n_states: Number of lowest eigenvalues/eigenvectors to compute.
                Requesting fewer states is significantly faster for large grids
                because ``eigh`` uses a divide-and-conquer LAPACK driver.
                Defaults to 20.

        Returns:
            A :class:`~src.models.states.QuantumSystem` containing:

            - ``energies`` — 1D float64 array of shape ``(n_states,)``
              containing :math:`E_n` in dimensionless units of :math:`E_1`.
            - ``wavefunctions`` — 2D float64 array of shape
              ``(n_points, n_states)`` where column ``k`` is
              :math:`\\psi_k(x)`, normalised so that
              :math:`\\int |\\psi_k|^2\\, dx = 1`.
            - ``x_grid`` — the uniform spatial coordinate array.
            - ``potential`` — the originating :class:`~src.models.potentials.BasePotential`.

        Raises:
            ValueError: If ``n_states < 1`` or if the grid has fewer than 3
                points (making a 0-sized interior impossible).

        Notes:
            **Symmetrisation of** :math:`A`: The potential term
            :math:`(1/12) M_{10} \\cdot \\mathrm{diag}(U)` is not exactly
            symmetric for spatially varying :math:`U(x)`, because its
            :math:`(i, j)` off-diagonal entry is :math:`M_{10}[i,j] \\cdot U_j`
            while its transpose has :math:`M_{10}[i,j] \\cdot U_i`.  The
            asymmetry at neighbouring nodes is :math:`O(dx)`.  We symmetrise
            :math:`A` explicitly via :math:`A \\leftarrow (A + A^T)/2` before
            passing it to ``eigh``.  By first-order perturbation theory the
            anti-symmetric part has zero quadratic form on real eigenvectors
            (it is anti-Hermitian), so the eigenvalue error from this step is
            at most :math:`O(dx^2)` — well within the :math:`10^{-3}` tolerance
            required for all physics analyses in this project.
        """
        if n_states < 1:
            raise ValueError(f"n_states must be at least 1, got {n_states}.")

        # ── 1. Build spatial grid using the potential's own domain ─────────
        physics: PhysicsConfig = potential.config
        x_min, x_max = potential.get_domain()
        n_pts: int = self._config.n_points

        if n_pts < 3:
            raise ValueError(
                f"n_points must be at least 3 to have interior points, got {n_pts}."
            )

        x: NDArray[np.float64] = np.linspace(x_min, x_max, n_pts)
        dx: float = float(x[1] - x[0])

        # ── 2. Evaluate and scale the potential ─────────────────────────────
        V: NDArray[np.float64] = potential.evaluate(x)

        # Clamp any numerical infinities to a very large finite barrier so that
        # the matrix remains well-conditioned (inf * 0 arithmetic is avoided).
        V_clamped: NDArray[np.float64] = np.where(np.isinf(V), 1.0e20, V)

        # Scaling factor:  f = 2m / ħ²
        # The substitution ε = 2mE/ħ² transforms E (physical) to ε (dimensionless),
        # and U = f·V is the scaled potential entering the recurrence.
        f: float = 2.0 * physics.mass / physics.hbar**2
        U: NDArray[np.float64] = f * V_clamped

        # ── 3. Assemble the Numerov tridiagonal matrices ─────────────────────
        # Interior points only: indices 1 … N-2 (Dirichlet BCs at endpoints).
        n_int: int = n_pts - 2
        U_int: NDArray[np.float64] = U[1:-1]

        # M_{-2}: standard second-difference stencil  ─ diag=-2, off-diag=+1
        # This encodes  psi_{i+1} - 2*psi_i + psi_{i-1},  the O(dx²) approximation
        # of dx²*psi''.  The Numerov correction in B upgrades this to O(dx⁴).
        M_m2: NDArray[np.float64] = np.zeros((n_int, n_int))
        np.fill_diagonal(M_m2, -2.0)
        np.fill_diagonal(M_m2[1:], 1.0)   # lower off-diagonal
        np.fill_diagonal(M_m2[:, 1:], 1.0)  # upper off-diagonal

        # M_{10}: Numerov weight matrix  ─ diag=10, off-diag=+1
        # Encodes the Numerov correction weights (1/12)·[f_{i-1}ψ_{i-1} + 10·f_i·ψ_i + f_{i+1}ψ_{i+1}]
        # that approximate the fourth-derivative term in the Taylor expansion,
        # lifting the scheme to O(dx⁴) global accuracy.
        M_10: NDArray[np.float64] = np.zeros((n_int, n_int))
        np.fill_diagonal(M_10, 10.0)
        np.fill_diagonal(M_10[1:], 1.0)
        np.fill_diagonal(M_10[:, 1:], 1.0)

        # ── 4. Build the generalised eigenvalue matrices A and B ─────────────
        # Kinetic term (symmetric by construction):
        #   T = (-1/dx²) · M_{-2}
        T: NDArray[np.float64] = (-1.0 / dx**2) * M_m2

        # Potential term from Numerov recurrence (column-U weighted; see module docstring):
        #   V_mat = (1/12) · M_{10} · diag(U_int)
        # Off-diagonal element (i, j=i±1) = M_{10}[i,j]·U_j/12 ≠ M_{10}[i,j]·U_i/12.
        # This makes V_mat slightly non-symmetric for spatially varying U(x).
        V_mat: NDArray[np.float64] = (1.0 / 12.0) * (M_10 @ np.diag(U_int))

        # Generalised LHS: A = T + V_mat
        A: NDArray[np.float64] = T + V_mat

        # Explicit symmetrisation — required for eigh (see Notes in docstring).
        # For constant U (e.g. infinite square well), A is already symmetric and
        # this step is a no-op numerically.
        A = (A + A.T) * 0.5

        # Generalised RHS (Numerov mass matrix — symmetric positive definite):
        #   B = (1/12) · M_{10}
        # All eigenvalues of M_{10} are in (8, 12) (Gershgorin), so B is SPD.
        B: NDArray[np.float64] = (1.0 / 12.0) * M_10

        # ── 5. Solve generalised symmetric eigenvalue problem A·ψ = ε·B·ψ ────
        # scipy.linalg.eigh:
        #   - Exploits symmetry and positive-definiteness of B for stability.
        #   - Returns eigenvalues in ascending order — no post-sort needed.
        #   - subset_by_index avoids computing the full spectrum (O(N³) → O(N·k²)).
        n_req: int = min(n_states, n_int)
        epsilons: NDArray[np.float64]
        evecs: NDArray[np.float64]
        epsilons, evecs = eigh(A, B, subset_by_index=[0, n_req - 1])

        # ── 6. Recover physical energies:  E_n = ε_n · ħ² / (2m) = ε_n / f ──
        energies: NDArray[np.float64] = epsilons / f

        # ── 7. Reconstruct full wavefunctions with Dirichlet-zero padding ────
        psi_full: NDArray[np.float64] = np.zeros((n_pts, n_req), dtype=np.float64)
        psi_full[1:-1, :] = evecs

        # ── 8. Normalise each wavefunction using the trapezoidal rule ─────────
        for k in range(n_req):
            psi_full[:, k] = self._normalize(x, psi_full[:, k])

        # ── 9. Wrap into QuantumSystem domain model ───────────────────────────
        # Bundling eliminates data desynchronisation: x, energies, and
        # wavefunctions are created together and validated together.  A consumer
        # can never accidentally pair wavefunctions from one solve call with a
        # grid rebuilt with different parameters.
        return QuantumSystem(
            energies=energies,
            wavefunctions=psi_full,
            x_grid=x,
            potential=potential,
        )

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _normalize(
        self,
        x: NDArray[np.float64],
        psi: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        """Normalises a wavefunction so that :math:`\\int |\\psi(x)|^2 dx = 1`.

        Uses the trapezoidal rule to approximate the integral, consistent with
        the legacy implementation in ``research_prototype.ipynb``.

        Args:
            x: 1D float64 array of spatial grid coordinates.
            psi: 1D float64 array of wavefunction values on the same grid.

        Returns:
            The normalised wavefunction array, same shape as ``psi``.
            If the computed norm is zero or negative (degenerate / null
            state), the original array is returned unchanged to avoid
            division by zero.

        Notes:
            The trapezoidal rule achieves :math:`O(dx^2)` accuracy for the
            normalisation integral.  For smooth wavefunctions on the grids used
            in this project (600+ points), this is numerically exact to machine
            precision relative to the :math:`10^{-4}` convergence tolerance.
        """
        norm: float = float(trapezoid(psi**2, x))
        if norm <= 0.0:
            return psi
        return psi / np.sqrt(norm)

    # ------------------------------------------------------------------
    # Dunder helpers
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        """Returns an unambiguous string representation.

        Returns:
            A string in the format ``NumerovSolver(config=<NumericalConfig>)``.
        """
        return f"NumerovSolver(config={self._config!r})"
