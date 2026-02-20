"""Domain model for the output of a quantum eigenvalue solver.

The :class:`QuantumSystem` class replaces the raw ``tuple[NDArray, NDArray]``
previously returned by :class:`~src.core.solvers.BaseSolver`. It solves the
**Data Desynchronisation** problem: in the legacy prototype, the spatial grid
``x``, the energy array, and the wavefunction matrix were three separate
variables that could silently diverge (e.g. if ``x`` was rebuilt with different
parameters while the wavefunctions were kept from an earlier run).

``QuantumSystem`` bundles all four pieces of a solved quantum problem into a
single, self-validating, immutable object:

- **``x_grid``** — the spatial coordinates on which wavefunctions are defined.
- **``energies``** — the eigenvalues :math:`E_n` for each state.
- **``wavefunctions``** — the normalised :math:`\\psi_n(x)` columns.
- **``potential``** — the :class:`~src.models.potentials.BasePotential` that
  was solved, providing geometry and physical-constant context.

Shape consistency is enforced at construction time, making silent
desynchronisation impossible.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from src.models.potentials import BasePotential


class QuantumSystem:
    """Immutable domain model encapsulating a solved 1D quantum system.

    Bundles the spatial grid, eigenvalues, normalised eigenvectors, and the
    originating potential into a single validated object.  All data is
    accessible through read-only properties; no public attribute permits
    post-construction mutation.

    Args:
        energies: 1D float64 array of shape ``(n_states,)`` containing the
            eigenvalues :math:`E_n` in dimensionless units of :math:`E_1`,
            sorted in ascending order.
        wavefunctions: 2D float64 array of shape ``(n_points, n_states)``
            where column ``k`` is :math:`\\psi_k(x)`, normalised so that
            :math:`\\int |\\psi_k|^2\\, dx = 1`.
        x_grid: 1D float64 array of shape ``(n_points,)`` containing the
            spatial coordinates in dimensionless length units on which every
            column of ``wavefunctions`` is defined.
        potential: The concrete :class:`~src.models.potentials.BasePotential`
            instance that was solved.  Retained for provenance — it carries
            both the geometry and the :class:`~src.models.config.PhysicsConfig`
            used during the solve.

    Raises:
        ValueError: If ``len(energies) != wavefunctions.shape[1]``,
            meaning the number of eigenvalues does not match the number of
            wavefunction columns (state-count desynchronisation).
        ValueError: If ``len(x_grid) != wavefunctions.shape[0]``,
            meaning the spatial grid length does not match the wavefunction
            row count (grid-wavefunction desynchronisation).

    Example:
        >>> from src.core import NumerovSolver
        >>> from src.models import InfiniteSquareWell
        >>> system = NumerovSolver().solve(InfiniteSquareWell(), n_states=5)
        >>> system.n_states
        5
        >>> system.ground_state_energy   # should be ≈ 1.0 in E₁ units
        1.0
        >>> prob = system.get_probability_density(0)  # |ψ₁(x)|²
        >>> prob.shape == system.x_grid.shape
        True
    """

    def __init__(
        self,
        energies: NDArray[np.float64],
        wavefunctions: NDArray[np.float64],
        x_grid: NDArray[np.float64],
        potential: BasePotential,
    ) -> None:
        """Validates shapes and stores all data as private attributes.

        Args:
            energies: Eigenvalue array of shape ``(n_states,)``.
            wavefunctions: Eigenvector matrix of shape ``(n_points, n_states)``.
            x_grid: Spatial coordinate array of shape ``(n_points,)``.
            potential: The :class:`~src.models.potentials.BasePotential` that
                was solved to produce these eigenpairs.

        Raises:
            ValueError: If the energy count does not match the wavefunction
                column count (state-count desynchronisation).
            ValueError: If the grid length does not match the wavefunction
                row count (grid-wavefunction desynchronisation).
        """
        n_states: int = int(energies.shape[0])
        n_wf_states: int = int(wavefunctions.shape[1])
        n_points_wf: int = int(wavefunctions.shape[0])
        n_points_grid: int = int(x_grid.shape[0])

        if n_states != n_wf_states:
            raise ValueError(
                f"State-count desynchronisation: len(energies)={n_states} "
                f"!= wavefunctions.shape[1]={n_wf_states}."
            )
        if n_points_grid != n_points_wf:
            raise ValueError(
                f"Grid-wavefunction desynchronisation: len(x_grid)={n_points_grid} "
                f"!= wavefunctions.shape[0]={n_points_wf}."
            )

        # Store everything as private — all access goes through read-only properties.
        self._energies: NDArray[np.float64] = energies.copy()
        self._wavefunctions: NDArray[np.float64] = wavefunctions.copy()
        self._x_grid: NDArray[np.float64] = x_grid.copy()
        self._potential: BasePotential = potential

    # ------------------------------------------------------------------
    # Read-only properties — raw data access
    # ------------------------------------------------------------------

    @property
    def energies(self) -> NDArray[np.float64]:
        """NDArray[np.float64]: Eigenvalues :math:`E_n`, shape ``(n_states,)``.

        Returned as a read-only view; modifying the returned array does not
        affect the internal state of this object.
        """
        view = self._energies.view()
        view.flags.writeable = False
        return view

    @property
    def wavefunctions(self) -> NDArray[np.float64]:
        """NDArray[np.float64]: Normalised :math:`\\psi_n(x)`, shape ``(n_points, n_states)``.

        Column ``k`` corresponds to the ``k``-th eigenstate.  Returned as a
        read-only view.
        """
        view = self._wavefunctions.view()
        view.flags.writeable = False
        return view

    @property
    def x_grid(self) -> NDArray[np.float64]:
        """NDArray[np.float64]: Spatial coordinates in dimensionless length units, shape ``(n_points,)``.

        The spatial axis shared by all wavefunction columns.  Returned as a
        read-only view.
        """
        view = self._x_grid.view()
        view.flags.writeable = False
        return view

    @property
    def potential(self) -> BasePotential:
        """BasePotential: The originating potential instance, read-only."""
        return self._potential

    # ------------------------------------------------------------------
    # Derived read-only properties — computed on demand
    # ------------------------------------------------------------------

    @property
    def n_states(self) -> int:
        """int: Number of eigenstates contained in this system."""
        return int(self._energies.shape[0])

    @property
    def ground_state_energy(self) -> float:
        """float: Energy of the lowest eigenstate :math:`E_1` in units of :math:`E_1`.

        Equivalent to ``self.energies[0]``, provided as a named property for
        readability in downstream code (e.g. thermal partition functions).
        """
        return float(self._energies[0])

    # ------------------------------------------------------------------
    # Public query methods
    # ------------------------------------------------------------------

    def get_probability_density(self, state_index: int) -> NDArray[np.float64]:
        """Returns the probability density :math:`|\\psi_n(x)|^2` for one eigenstate.

        Args:
            state_index: Zero-based index of the eigenstate.  ``0`` is the
                ground state, ``1`` is the first excited state, and so on.

        Returns:
            1D float64 array of shape ``(n_points,)`` containing
            :math:`|\\psi_n(x)|^2` evaluated on ``self.x_grid``.

        Raises:
            IndexError: If ``state_index`` is outside ``[0, n_states)``.

        Example:
            >>> rho = system.get_probability_density(0)  # ground state density
            >>> from scipy.integrate import trapezoid
            >>> trapezoid(rho, system.x_grid)  # ≈ 1.0 by normalisation
            1.0
        """
        if not (0 <= state_index < self.n_states):
            raise IndexError(
                f"state_index={state_index} is out of range for a system "
                f"with {self.n_states} states (valid: 0 to {self.n_states - 1})."
            )
        return self._wavefunctions[:, state_index] ** 2

    def get_state(self, state_index: int) -> NDArray[np.float64]:
        """Returns the wavefunction :math:`\\psi_n(x)` for one eigenstate.

        Args:
            state_index: Zero-based index of the eigenstate.

        Returns:
            1D float64 array of shape ``(n_points,)`` containing
            :math:`\\psi_n(x)` on ``self.x_grid``.

        Raises:
            IndexError: If ``state_index`` is outside ``[0, n_states)``.
        """
        if not (0 <= state_index < self.n_states):
            raise IndexError(
                f"state_index={state_index} is out of range for a system "
                f"with {self.n_states} states (valid: 0 to {self.n_states - 1})."
            )
        return self._wavefunctions[:, state_index].copy()

    # ------------------------------------------------------------------
    # Dunder helpers
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        """Returns an unambiguous string representation.

        Returns:
            A string summarising the potential type, number of states, and
            ground-state energy.
        """
        return (
            f"QuantumSystem("
            f"potential={self._potential!r}, "
            f"n_states={self.n_states}, "
            f"E_0={self.ground_state_energy:.6f})"
        )
