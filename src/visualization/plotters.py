"""Professional rendering service for 1D quantum-system visualisations.

This module provides :class:`QuantumPlotter`, a stateless service class that
consumes :class:`~src.models.states.QuantumSystem` objects and pre-computed
density arrays and renders publication-quality figures.

Design: Separation of Concerns
-------------------------------
The plotter is a *pure rendering* layer.  It enforces the following contract:

1. **No physics**: The plotter never instantiates a solver, a potential, or a
   thermal engine.  All numerical data arrive as method arguments.
2. **Dependency Injection**: Callers pass in the ``QuantumSystem`` and any
   pre-computed arrays (e.g. ``rho``, ``rho_2p``) produced by
   :class:`~src.core.statistics.ThermalEngine`.  The plotter is blind to how
   those arrays were produced.
3. **Stateless rendering**: Each method returns a ``(Figure, Axes)`` pair and
   optionally saves to disk.  No state is accumulated between calls; the class
   instance only carries style configuration.
4. **Back-end isolation**: ``matplotlib`` is imported only here.  Replacing this
   module with a Plotly or Bokeh implementation requires zero changes to the
   physics core.

Extensibility notes
-------------------
- A ``PlotlyQuantumPlotter`` could implement the same method signatures and be
  swapped in without touching ``NumerovSolver`` or ``ThermalEngine``.
- A headless PDF pipeline can call these methods with ``save_path`` set and
  ``plt.show()`` suppressed (controlled by the ``show`` parameter).
- A Jupyter widget layer can call ``fig.show()`` on the returned ``Figure``
  object after the method returns.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3D projection)
from numpy.typing import NDArray

from src.models.states import QuantumSystem


# ---------------------------------------------------------------------------
# Style constants — kept here so a future alternative back-end can ignore them
# ---------------------------------------------------------------------------

_STYLE: str = "seaborn-v0_8-whitegrid"  # clean academic look, matplotlib ≥ 3.6

# Colour palette: one distinct colour per plotted state (cycles for n > 8)
_STATE_COLOURS: list[str] = [
    "#1f77b4",  # blue       — n=1
    "#d62728",  # red        — n=2
    "#2ca02c",  # green      — n=3
    "#ff7f0e",  # orange     — n=4
    "#9467bd",  # purple     — n=5
    "#8c564b",  # brown      — n=6
    "#e377c2",  # pink       — n=7
    "#7f7f7f",  # grey       — n=8
]

# LaTeX strings used on axes — defined once for consistency
_XLABEL_X: str = r"$x / L$"
_YLABEL_PSI2: str = r"$|\psi_n(x)|^2$"
_YLABEL_RHO: str = r"$\rho_{\mathrm{th}}(x)$"
_XLABEL_X1: str = r"$x_1 / L$"
_YLABEL_X2: str = r"$x_2 / L$"
_ZLABEL_RHO2P: str = r"$\rho(x_1, x_2)$"


class QuantumPlotter:
    """Stateless rendering service for 1D quantum-system visualisations.

    Consumes pre-computed :class:`~src.models.states.QuantumSystem` objects
    and density arrays; never performs any physics calculation.

    All methods return the ``(Figure, Axes)`` pair so callers can make further
    customisations before saving or displaying the figure.

    Args:
        figsize: Default figure size ``(width_inches, height_inches)`` used when
            a method does not override it.  Defaults to ``(9, 5)``.
        dpi: Dots-per-inch for on-screen rendering and for ``save_path`` output.
            Defaults to ``150``.

    Example:
        >>> from src.core import NumerovSolver
        >>> from src.models import InfiniteSquareWell
        >>> from src.visualization import QuantumPlotter
        >>>
        >>> system = NumerovSolver().solve(InfiniteSquareWell(), n_states=5)
        >>> plotter = QuantumPlotter()
        >>> fig, ax = plotter.plot_wavefunctions(system, n_states=3)
        >>> fig.savefig("wavefunctions.png", dpi=150)
    """

    def __init__(
        self,
        figsize: tuple[float, float] = (9.0, 5.0),
        dpi: int = 150,
    ) -> None:
        """Configures global matplotlib style and stores rendering defaults.

        Args:
            figsize: Default ``(width, height)`` in inches for all figures.
            dpi: Resolution used when rendering figures to screen or file.
        """
        self._figsize: tuple[float, float] = figsize
        self._dpi: int = dpi

        # Apply a consistent academic style.  Wrapped in try/except so the
        # plotter degrades gracefully if the style name changes across
        # matplotlib versions.
        try:
            plt.style.use(_STYLE)
        except OSError:
            # Fall back to the built-in 'ggplot' if seaborn-v0_8 is unavailable.
            plt.style.use("ggplot")

        # Global rcParams tweaks — applied once at construction time.
        plt.rcParams.update({
            "font.size": 12,
            "axes.titlesize": 13,
            "axes.labelsize": 12,
            "legend.fontsize": 10,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "figure.dpi": dpi,
            # Use LaTeX-compatible math text without requiring a full TeX install.
            "mathtext.fontset": "dejavusans",
        })

    # ------------------------------------------------------------------
    # Public plotting methods
    # ------------------------------------------------------------------

    def plot_wavefunctions(
        self,
        system: QuantumSystem,
        n_states: int = 3,
        show: bool = False,
        save_path: Optional[str | Path] = None,
    ) -> tuple[Figure, Axes]:
        """Plots the first ``n_states`` normalised probability densities.

        Each eigenstate :math:`|\\psi_n(x)|^2` is plotted on a shared axis
        and *vertically offset* by its energy :math:`E_n`, so the plot reads
        like a textbook energy-level diagram.  The potential profile is drawn
        as a shaded background for geometric context.

        The x-axis is normalised by the well width :math:`L` so the plot is
        unit-agnostic.

        Args:
            system: A solved :class:`~src.models.states.QuantumSystem`
                containing at least ``n_states`` eigenpairs.
            n_states: Number of eigenstates to display, starting from the
                ground state (:math:`n=1`).  Must be at most
                ``system.n_states``.  Defaults to 3.
            show: If ``True``, call ``plt.show()`` before returning.
                Defaults to ``False`` for headless / notebook use.
            save_path: If provided, save the figure to this file path
                before returning.

        Returns:
            Tuple ``(fig, ax)`` for further customisation.

        Raises:
            ValueError: If ``n_states > system.n_states``.

        Example:
            >>> fig, ax = plotter.plot_wavefunctions(system, n_states=5)
            >>> fig.savefig("energy_levels.pdf")
        """
        n_plot: int = min(n_states, system.n_states)
        if n_states > system.n_states:
            raise ValueError(
                f"Requested n_states={n_states} but system only has "
                f"{system.n_states} solved states."
            )

        L: float = system.potential.config.well_width  # characteristic length
        x_norm: NDArray[np.float64] = system.x_grid / L

        fig, ax = plt.subplots(figsize=self._figsize)

        # ── Compute scale FIRST so ylim is derived from actual wavefunction heights ─
        # Scale factor: proportional to the energy gap between consecutive levels
        # so each density curve fills roughly the same vertical slice.
        E0: float = float(system.energies[0])
        E_top_state: float = float(system.energies[n_plot - 1])
        E_range: float = E_top_state - E0
        scale: float = (E_range / n_plot) * 0.8 if n_plot > 1 else 1.0

        # True ceiling of the highest plotted curve: E_n + scale·max(|ψ_n(x)|²).
        # Using the actual wavefunction peak (not a fixed energy multiplier) ensures
        # the top state is never clipped regardless of the scale factor.
        rho_top: NDArray[np.float64] = system.get_probability_density(n_plot - 1)
        wf_ceiling: float = E_top_state + scale * float(np.max(rho_top))
        # Add 20 % breathing room above the highest peak.
        ylim_top: float = wf_ceiling * 1.20

        # ── Background: potential profile ────────────────────────────────
        V: NDArray[np.float64] = system.potential.evaluate(system.x_grid)
        # Clamp to the computed display ceiling so infinite walls don't
        # swamp the energy scale of the plotted states.
        V_display: NDArray[np.float64] = np.clip(V, None, ylim_top)
        ax.fill_between(
            x_norm, V_display, ylim_top,
            color="#dce9f5", alpha=0.5, linewidth=0, label="Potential",
        )
        ax.plot(x_norm, V_display, color="#5a87b5", linewidth=1.2, alpha=0.7)

        # ── Eigenstates: |ψ_n|² offset by E_n ───────────────────────────
        for k in range(n_plot):
            E_k: float = float(system.energies[k])
            colour: str = _STATE_COLOURS[k % len(_STATE_COLOURS)]
            rho_k: NDArray[np.float64] = system.get_probability_density(k)

            # Horizontal baseline at energy level E_k
            ax.axhline(E_k, color=colour, linewidth=0.8, linestyle="--", alpha=0.5)

            # Wavefunction density scaled and offset to sit on its energy level
            ax.plot(
                x_norm,
                E_k + scale * rho_k,
                color=colour,
                linewidth=1.8,
                label=fr"$n={k+1},\ E_{k+1}={E_k:.2f}$",
            )
            # Fill under each density curve for emphasis
            ax.fill_between(
                x_norm,
                E_k,
                E_k + scale * rho_k,
                color=colour,
                alpha=0.18,
            )

        ax.set_xlim(x_norm[0], x_norm[-1])
        ax.set_ylim(bottom=0.0, top=ylim_top)
        ax.set_xlabel(_XLABEL_X)
        ax.set_ylabel(r"$E_n$ (dimensionless)")
        ax.set_title(
            fr"Probability densities $|\psi_n(x)|^2$ — "
            fr"{type(system.potential).__name__}",
        )
        ax.legend(loc="upper right", framealpha=0.85)
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))

        fig.tight_layout()
        self._save_and_show(fig, save_path, show)
        return fig, ax

    def plot_thermal_density(
        self,
        system: QuantumSystem,
        rho: NDArray[np.float64],
        title: str = "Single-particle thermal density",
        show: bool = False,
        save_path: Optional[str | Path] = None,
    ) -> tuple[Figure, Axes]:
        """Plots a single-particle thermal probability density :math:`\\rho_{\\text{th}}(x)`.

        The spatial axis is normalised by the characteristic well width
        :math:`L` from ``system.potential.config.well_width``.

        Args:
            system: The :class:`~src.models.states.QuantumSystem` whose
                ``x_grid`` and ``potential`` provide context (well width
                and potential profile for the background shading).
            rho: 1D float64 array of shape ``(n_points,)`` containing the
                thermal probability density, as returned by
                :meth:`~src.core.statistics.ThermalEngine.calculate_thermal_density`.
            title: Axes title string.  Supports LaTeX math notation.
            show: If ``True``, call ``plt.show()`` before returning.
            save_path: Optional file path to save the figure.

        Returns:
            Tuple ``(fig, ax)``.

        Example:
            >>> rho = engine.calculate_thermal_density(system, temperature=5.0)
            >>> fig, ax = plotter.plot_thermal_density(system, rho, title=r"$T=5.0$")
        """
        L: float = system.potential.config.well_width
        x_norm: NDArray[np.float64] = system.x_grid / L

        fig, ax = plt.subplots(figsize=self._figsize)

        # Shaded potential background (clipped to density scale).
        # Guard: for potentials where all finite values are 0 (e.g.
        # InfiniteSquareWell where V=0 everywhere inside), skip the
        # background shading entirely to avoid a division-by-zero.
        V: NDArray[np.float64] = system.potential.evaluate(system.x_grid)
        rho_max: float = float(np.max(rho)) * 1.35
        finite_vals = V[np.isfinite(V)]
        v_max_finite: float = float(np.max(finite_vals)) if len(finite_vals) > 0 else 0.0
        if v_max_finite > 0.0:
            V_display: NDArray[np.float64] = (
                np.clip(V / v_max_finite, 0.0, 1.0) * rho_max * 0.35
            )
            ax.fill_between(
                x_norm, V_display, rho_max,
                color="#e8f0fb", alpha=0.5, linewidth=0,
            )

        # Thermal density
        ax.plot(x_norm, rho, color="#1f77b4", linewidth=2.2, label=r"$\rho_{\mathrm{th}}(x)$")
        ax.fill_between(x_norm, 0, rho, color="#1f77b4", alpha=0.15)

        ax.set_xlim(x_norm[0], x_norm[-1])
        ax.set_ylim(bottom=0.0, top=rho_max)
        ax.set_xlabel(_XLABEL_X)
        ax.set_ylabel(_YLABEL_RHO)
        ax.set_title(title)
        ax.legend(framealpha=0.85)

        fig.tight_layout()
        self._save_and_show(fig, save_path, show)
        return fig, ax

    def plot_pair_density_heatmap(
        self,
        system: QuantumSystem,
        rho_2p: NDArray[np.float64],
        title: str = r"Two-particle pair density $\rho(x_1, x_2)$",
        show: bool = False,
        save_path: Optional[str | Path] = None,
    ) -> tuple[Figure, Axes]:
        """Renders the two-particle density matrix as a 2D ``pcolormesh``.

        The axes are normalised by :math:`L`.  The diagonal :math:`x_1 = x_2`
        is annotated with a dashed white line so the Pauli exclusion hole (for
        fermions) or the bunching ridge (for bosons) is immediately visible.

        Well boundaries from ``system.potential.get_domain()`` are drawn as
        thin white reference lines on both axes.

        Args:
            system: The solved :class:`~src.models.states.QuantumSystem`
                providing the spatial grid and potential geometry.
            rho_2p: 2D float64 array of shape ``(n_points, n_points)``
                containing :math:`\\rho(x_1, x_2)`, as returned by
                :meth:`~src.core.statistics.ThermalEngine.calculate_pair_density`
                or :meth:`~src.core.statistics.ThermalEngine.get_pair_density`.
            title: Axes title string (LaTeX supported).
            show: If ``True``, call ``plt.show()`` before returning.
            save_path: Optional file path to save the figure.

        Returns:
            Tuple ``(fig, ax)``.
        """
        L: float = system.potential.config.well_width
        x_norm: NDArray[np.float64] = system.x_grid / L

        fig, ax = plt.subplots(figsize=(7.0, 6.0))

        pcm = ax.pcolormesh(
            x_norm, x_norm, rho_2p,
            cmap="viridis",
            shading="auto",
            vmin=0.0,
            vmax=float(np.max(rho_2p)),
        )

        # Colourbar with physical label
        cbar = fig.colorbar(pcm, ax=ax, fraction=0.045, pad=0.04)
        cbar.set_label(_ZLABEL_RHO2P, fontsize=11)
        cbar.ax.tick_params(labelsize=9)

        # Diagonal reference: x1 = x2 (exchange symmetry axis)
        ax.plot(
            [x_norm[0], x_norm[-1]],
            [x_norm[0], x_norm[-1]],
            color="white", linewidth=1.4, linestyle="--", alpha=0.85,
            label=r"$x_1 = x_2$",
        )

        # Well-boundary reference lines from potential domain
        x_lo, x_hi = system.potential.get_domain()
        for boundary in (x_lo / L, x_hi / L):
            ax.axvline(boundary, color="white", linewidth=0.9, linestyle=":", alpha=0.6)
            ax.axhline(boundary, color="white", linewidth=0.9, linestyle=":", alpha=0.6)

        ax.set_aspect("equal")
        ax.set_xlabel(_XLABEL_X1)
        ax.set_ylabel(_YLABEL_X2)
        ax.set_title(title)
        ax.legend(loc="upper left", fontsize=9, framealpha=0.7)

        fig.tight_layout()
        self._save_and_show(fig, save_path, show)
        return fig, ax

    def plot_pair_density_3d(
        self,
        system: QuantumSystem,
        rho_2p: NDArray[np.float64],
        title: str = r"Two-particle pair density $\rho(x_1, x_2)$",
        downsample: int = 4,
        show: bool = False,
        save_path: Optional[str | Path] = None,
    ) -> tuple[Figure, Axes3D]:
        """Renders the two-particle density matrix as a 3D surface plot.

        The full ``(n_points × n_points)`` grid is downsampled by ``downsample``
        before surface rendering to keep the triangle count manageable without
        losing the qualitative structure of the exchange hole / bunching ridge.

        Args:
            system: The solved :class:`~src.models.states.QuantumSystem`.
            rho_2p: 2D float64 array of shape ``(n_points, n_points)``
                containing :math:`\\rho(x_1, x_2)`.
            title: Figure title string (LaTeX supported).
            downsample: Keep every ``downsample``-th grid point along each
                axis when building the surface mesh.  ``1`` keeps all points
                (may be slow for ``n_points=600``).  Defaults to ``4``.
            show: If ``True``, call ``plt.show()`` before returning.
            save_path: Optional file path to save the figure.

        Returns:
            Tuple ``(fig, ax3d)`` where ``ax3d`` is the
            :class:`mpl_toolkits.mplot3d.axes3d.Axes3D` instance.
        """
        L: float = system.potential.config.well_width
        x_norm: NDArray[np.float64] = system.x_grid / L

        # Downsample for performance — surface rendering is O(N²) in triangles.
        x_ds: NDArray[np.float64] = x_norm[::downsample]
        Z_ds: NDArray[np.float64] = rho_2p[::downsample, ::downsample]
        X1, X2 = np.meshgrid(x_ds, x_ds)

        fig = plt.figure(figsize=(8.5, 6.5))
        ax3d: Axes3D = fig.add_subplot(111, projection="3d")

        surf = ax3d.plot_surface(
            X1, X2, Z_ds,
            cmap="viridis",
            linewidth=0,
            antialiased=True,
            alpha=0.92,
        )

        # Colourbar mapped to the surface norm — uses the full rho_2p range.
        norm = Normalize(vmin=0.0, vmax=float(np.max(rho_2p)))
        sm = ScalarMappable(cmap="viridis", norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax3d, shrink=0.55, pad=0.08)
        cbar.set_label(_ZLABEL_RHO2P, fontsize=10)
        cbar.ax.tick_params(labelsize=8)

        ax3d.set_xlabel(_XLABEL_X1, fontsize=10, labelpad=8)
        ax3d.set_ylabel(_YLABEL_X2, fontsize=10, labelpad=8)
        ax3d.set_zlabel(_ZLABEL_RHO2P, fontsize=10, labelpad=8)
        ax3d.set_title(title, pad=12)
        ax3d.view_init(elev=28, azim=40)
        ax3d.tick_params(axis="both", labelsize=8)

        fig.tight_layout()
        self._save_and_show(fig, save_path, show)
        return fig, ax3d  # type: ignore[return-value]

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _save_and_show(
        self,
        fig: Figure,
        save_path: Optional[str | Path],
        show: bool,
    ) -> None:
        """Saves the figure to disk and/or displays it.

        Args:
            fig: The :class:`matplotlib.figure.Figure` to act on.
            save_path: File path for saving (``None`` skips saving).
            show: If ``True``, call ``plt.show()``.
        """
        if save_path is not None:
            path = Path(save_path)
            path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(path, dpi=self._dpi, bbox_inches="tight")
        if show:
            plt.show()
        else:
            plt.close(fig)
