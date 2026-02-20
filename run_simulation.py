"""Master integration entry point for the Thermal-Quantum-Numerov-Solver.

This script is the **front door** of the library — the single file a user
runs to execute the full simulation pipeline end-to-end.  It demonstrates
how the individual LEGO bricks assembled across Phases 1–6 snap together:

    PhysicsConfig / NumericalConfig   (Phase 1 — immutable configuration)
            ↓
    FiniteSquareWell / BasePotential  (Phase 2 — potential hierarchy)
            ↓
    NumerovSolver → QuantumSystem     (Phase 3–4 — Numerov engine + domain model)
            ↓
    ThermalEngine                     (Phase 5 — Boltzmann statistics + symmetrisation)
            ↓
    QuantumPlotter                    (Phase 6 — professional visualisation)

Usage
-----
    python run_simulation.py

All figures are written to ``research_output/`` using :mod:`pathlib` so the
script is safe to run from any working directory.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from src.core import NumerovSolver, ThermalEngine
from src.models import (
    FiniteSquareWell,
    NumericalConfig,
    ParticleType,
    PhysicsConfig,
    QuantumSystem,
)
from src.visualization import QuantumPlotter

# ---------------------------------------------------------------------------
# Output directory — resolved relative to this script's location so it is
# safe to call from any working directory.
# ---------------------------------------------------------------------------
_OUTPUT_DIR: Path = Path(__file__).parent / "research_output"


# ---------------------------------------------------------------------------
# Primary simulation pipeline
# ---------------------------------------------------------------------------


def run_pipeline(
    v0: float = 50.0,
    temperature: float = 5.0,
    n_states: int = 12,
    output_dir: Path = _OUTPUT_DIR,
) -> QuantumSystem:
    """Runs the complete single-potential simulation pipeline.

    Steps:
        1. Build an immutable :class:`~src.models.config.PhysicsConfig` and
           :class:`~src.models.config.NumericalConfig`.
        2. Instantiate a :class:`~src.models.potentials.FiniteSquareWell`
           with barrier height ``v0``.
        3. Solve the TISE with :class:`~src.core.solvers.NumerovSolver`,
           returning a :class:`~src.models.states.QuantumSystem`.
        4. Compute the single-particle thermal density via
           :class:`~src.core.statistics.ThermalEngine`.
        5. Compute the fermion pair density for the lowest two states.
        6. Save three figures (wavefunctions, thermal density, 2D heatmap)
           into ``output_dir``.

    Args:
        v0: Barrier height :math:`V_0` for the finite square well in
            dimensionless energy units.  Defaults to ``50.0``.
        temperature: Dimensionless temperature :math:`T` for Boltzmann
            weighting.  Defaults to ``5.0``.
        n_states: Number of eigenstates to solve for.  Defaults to ``12``.
        output_dir: Directory for saving figures.  Created if absent.

    Returns:
        The solved :class:`~src.models.states.QuantumSystem` for downstream
        use (e.g. in ``run_mass_sweep``).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Step 1: Configuration ────────────────────────────────────────────────
    physics = PhysicsConfig(
        hbar=1.0,
        mass=0.5,
        kb=1.0,
        well_width=math.pi,
    )
    numerics = NumericalConfig(
        n_points=600,
    )

    print(f"  Physics : {physics}")
    print(f"  Numerics: {numerics}")

    # ── Step 2: Potential ────────────────────────────────────────────────────
    potential = FiniteSquareWell(config=physics, v0=v0)
    x_lo, x_hi = potential.get_domain()
    print(f"  Potential: FiniteSquareWell(V0={v0})  domain=[{x_lo:.3f}, {x_hi:.3f}]")

    # ── Step 3: Solve TISE ───────────────────────────────────────────────────
    solver = NumerovSolver(numerics)
    system = solver.solve(potential, n_states=n_states)
    print(f"  Solved {system.n_states} states  |  E_0={system.ground_state_energy:.6f}")
    print(f"  Energies: {np.round(system.energies[:5], 4)}")

    # ── Step 4: Thermal density ──────────────────────────────────────────────
    engine = ThermalEngine()
    rho_th = engine.calculate_thermal_density(system, temperature=temperature)
    print(f"  ∫ρ_th dx ≈ {float(np.trapezoid(rho_th, system.x_grid)):.8f}  (should be 1.0)")

    # ── Step 5: Fermion pair density (n=1, m=2) ──────────────────────────────
    rho_pair = engine.calculate_pair_density(
        system, temperature=temperature, p_type=ParticleType.FERMION
    )
    diag_max = float(np.max(np.abs(np.diag(rho_pair))))
    print(f"  Fermion diagonal max = {diag_max:.2e}  (Pauli exclusion: should be ≈ 0)")

    # ── Step 6: Visualise ────────────────────────────────────────────────────
    plotter = QuantumPlotter(figsize=(9, 5), dpi=150)

    plotter.plot_wavefunctions(
        system,
        n_states=min(5, system.n_states),
        save_path=output_dir / f"fsw_v0{int(v0)}_wavefunctions.png",
    )
    print(f"  ✓ Wavefunctions saved")

    plotter.plot_thermal_density(
        system,
        rho_th,
        title=fr"Finite Square Well ($V_0={v0}$) — Thermal density $T={temperature}$",
        save_path=output_dir / f"fsw_v0{int(v0)}_thermal_density.png",
    )
    print(f"  ✓ Thermal density saved")

    plotter.plot_pair_density_heatmap(
        system,
        rho_pair,
        title=(
            fr"Fermion pair density — FiniteSquareWell($V_0={v0}$), $T={temperature}$"
        ),
        save_path=output_dir / f"fsw_v0{int(v0)}_fermion_heatmap.png",
    )
    print(f"  ✓ Fermion heatmap saved")

    return system


# ---------------------------------------------------------------------------
# Parameter sweep: varying particle mass
# ---------------------------------------------------------------------------


def run_mass_sweep(
    masses: list[float] | None = None,
    v0: float = 50.0,
    n_states: int = 8,
    output_dir: Path = _OUTPUT_DIR,
) -> None:
    """Demonstrates how ground-state energy scales with particle mass.

    Runs :meth:`NumerovSolver.solve` for each mass value, reporting the
    ground-state energy :math:`E_0(m)`.  This illustrates the **Open/Closed**
    principle in action: varying mass requires only a new
    :class:`~src.models.config.PhysicsConfig` — no solver code changes.

    The analytic expectation for the infinite square well is:

    .. math::

        E_1^{(\\infty)} = \\frac{\\hbar^2 \\pi^2}{2 m L^2}

    so :math:`E_0 \\propto 1/m`.

    Args:
        masses: List of dimensionless masses to sweep over.  Defaults to
            ``[0.1, 0.25, 0.5, 1.0, 2.0]``.
        v0: Barrier height for the :class:`~src.models.potentials.FiniteSquareWell`
            used in each sweep step.
        n_states: Number of eigenstates to solve per sweep step.
        output_dir: Output directory (figures are not generated here — the
            sweep is a console report only).
    """
    if masses is None:
        masses = [0.1, 0.25, 0.5, 1.0, 2.0]

    numerics = NumericalConfig(n_points=600)
    solver = NumerovSolver(numerics)

    print(f"\n  {'Mass':>8}  {'E_0 (numeric)':>16}  {'E_0 ∝ 1/m (rel.)':>18}  {'E_0·m (const?)':>16}")
    print(f"  {'-'*8}  {'-'*16}  {'-'*18}  {'-'*16}")

    ref_Em: float | None = None  # E_0 * m at reference mass for relative check
    results: list[tuple[float, float]] = []

    for m in masses:
        # Rebuild config with the new mass — frozen dataclass enforces immutability.
        cfg = PhysicsConfig(mass=m, well_width=math.pi)
        pot = FiniteSquareWell(config=cfg, v0=v0)
        system = solver.solve(pot, n_states=n_states)
        E0 = system.ground_state_energy
        Em = E0 * m  # should be ≈ constant = ħ²π²/(2L²) = 0.5 in these units
        if ref_Em is None:
            ref_Em = Em
        ratio = Em / ref_Em if ref_Em else 1.0
        print(f"  {m:>8.3f}  {E0:>16.6f}  {ratio:>18.6f}  {Em:>16.6f}")
        results.append((m, E0))

    # Qualitative check: E_0 * m should be roughly constant (1/m scaling)
    Em_values = [r[0] * r[1] for r in results]
    spread = (max(Em_values) - min(Em_values)) / np.mean(Em_values)
    trend_ok = spread < 0.15  # finite-well deviates slightly from analytic 1/m
    print(f"\n  E_0·m spread = {spread*100:.1f}%  "
          f"({'≈ constant ✓' if trend_ok else 'large spread (finite-well effect)'})")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")   # headless rendering — no display required

    print("=" * 64)
    print("  Thermal-Quantum-Numerov-Solver — Master Integration")
    print("=" * 64)

    # ── Primary pipeline: FiniteSquareWell V0=50 ────────────────────────────
    print("\n[1/2] Running primary simulation pipeline  (V0=50, T=5.0) …")
    system = run_pipeline(v0=50.0, temperature=5.0)

    # ── Parameter sweep: mass dependence ────────────────────────────────────
    print("\n[2/2] Running mass sweep  (V0=50, masses=[0.1, 0.25, 0.5, 1.0, 2.0]) …")
    run_mass_sweep(masses=[0.1, 0.25, 0.5, 1.0, 2.0], v0=50.0)

    print("\n" + "=" * 64)
    print("  Simulation complete.  Figures written to research_output/")
    print("=" * 64)
