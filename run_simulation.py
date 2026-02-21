"""Thermal-Quantum Numerov Solver — Master Simulation Entry Point.

This script is the **Front Door** of the library.  It demonstrates the complete
simulation pipeline by composing the four layers of the LEGO architecture:

    Config  →  Potential  →  Solver  →  Thermal Engine  →  Plotter

Each layer knows nothing about the layers above it:

- :mod:`src.models.config`        — immutable physical constants
- :mod:`src.models.potentials`    — geometry (what is the well?)
- :mod:`src.core.solvers`         — eigenpairs (what are the states?)
- :mod:`src.core.statistics`      — thermodynamics (how are they populated?)
- :mod:`src.visualization`        — rendering (how does it look?)

Run from the project root::

    python run_simulation.py

All output figures are saved to ``research_output/figures/`` using :mod:`pathlib`
for safe cross-platform path handling.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from src.core.solvers import NumerovSolver
from src.core.statistics import ThermalEngine
from src.models.config import NumericalConfig, PhysicsConfig
from src.models.potentials import FiniteSquareWell
from src.models.statistics import ParticleType
from src.models.states import QuantumSystem
from src.visualization.plotters import QuantumPlotter

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

ROOT: Path = Path(__file__).parent
OUTPUT_DIR: Path = ROOT / "research_output" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Main simulation pipeline
# ---------------------------------------------------------------------------

def run_main_simulation() -> None:
    """Executes the complete single-system simulation pipeline.

    Steps
    -----
    1. Initialise :class:`~src.models.config.PhysicsConfig` and
       :class:`~src.models.config.NumericalConfig` as frozen dataclasses.
    2. Instantiate a :class:`~src.models.potentials.FiniteSquareWell`
       with :math:`V_0 = 50`.
    3. Solve the TISE with :class:`~src.core.solvers.NumerovSolver`.
    4. Compute the Fermion pair density at :math:`T = 5.0` with
       :class:`~src.core.statistics.ThermalEngine`.
    5. Render and save three figures via
       :class:`~src.visualization.plotters.QuantumPlotter`.
    """
    print("=" * 64)
    print("  Thermal-Quantum Numerov Solver — Main Simulation")
    print("=" * 64)

    # ── Step 1: Configuration ────────────────────────────────────────────────
    # Frozen dataclasses: once created, these constants can never be mutated.
    physics = PhysicsConfig(
        hbar=1.0,
        mass=0.5,
        kb=1.0,
        well_width=math.pi,
    )
    numerics = NumericalConfig(
        n_points=600,
        convergence_tolerance=1e-4,
    )
    print(f"\n[1/5] Config       — {physics!r}")
    print(f"                     {numerics!r}")

    # ── Step 2: Potential ────────────────────────────────────────────────────
    # DI: PhysicsConfig is injected; FiniteSquareWell owns no global state.
    potential = FiniteSquareWell(v0=50.0, config=physics)
    print(f"\n[2/5] Potential    — {potential!r}")
    print(f"       Domain:  [{potential.get_domain()[0]:.4f},  {potential.get_domain()[1]:.4f}]")

    # ── Step 3: Solve ────────────────────────────────────────────────────────
    solver = NumerovSolver(numerics)
    system: QuantumSystem = solver.solve(potential, n_states=20)
    print(f"\n[3/5] Solver       — {solver!r}")
    print(f"       States:   {system.n_states}  |  E_0 = {system.ground_state_energy:.6f}")
    print(f"       Energies: {np.round(system.energies[:5], 4)}")

    # ── Step 4: Thermal statistics ───────────────────────────────────────────
    temperature: float = 5.0
    engine = ThermalEngine()
    rho_single = engine.calculate_thermal_density(system, temperature=temperature)
    rho_pair   = engine.calculate_pair_density(system, temperature=temperature,
                                               p_type=ParticleType.FERMION)
    print(f"\n[4/5] Thermal      — T = {temperature}  |  particle type: FERMION")
    print(f"       ∫ρ_single dx = {float(np.trapezoid(rho_single, system.x_grid)):.8f}")

    # ── Step 5: Visualisation ────────────────────────────────────────────────
    plotter = QuantumPlotter(figsize=(9, 5), dpi=150)

    # 5a – Energy-level diagram with probability densities
    out_wf = OUTPUT_DIR / "sim_wavefunctions.png"
    plotter.plot_wavefunctions(system, n_states=5, save_path=out_wf)
    print(f"\n[5/5] Plots saved:")
    print(f"       {out_wf.relative_to(ROOT)}")

    # 5b – Single-particle thermal density
    out_rho = OUTPUT_DIR / "sim_thermal_density.png"
    plotter.plot_thermal_density(
        system, rho_single,
        title=fr"Thermal density — FiniteSquareWell ($V_0=50$, $T={temperature}$)",
        save_path=out_rho,
    )
    print(f"       {out_rho.relative_to(ROOT)}")

    # 5c – Fermion pair density heatmap
    out_hm = OUTPUT_DIR / "sim_pair_heatmap_fermion.png"
    plotter.plot_pair_density_heatmap(
        system, rho_pair,
        title=fr"Fermion pair density — FiniteSquareWell ($V_0=50$, $T={temperature}$)",
        save_path=out_hm,
    )
    print(f"       {out_hm.relative_to(ROOT)}")

    print("\n✓ Main simulation complete.\n")


# ---------------------------------------------------------------------------
# Parameter sweep: mass dependence of ground-state energy
# ---------------------------------------------------------------------------

def run_mass_sweep() -> None:
    """Demonstrates how the ground-state energy :math:`E_0` scales with mass.

    The TISE scaling relation in these dimensionless units is:

    .. math::

        E_0 \\propto \\frac{\\hbar^2}{2m L^2}

    Sweeping over mass values with an otherwise identical configuration
    isolates this dependence cleanly.  The result is printed as a table and
    the wavefunction diagrams for each mass are overlaid in a single saved
    figure.

    Masses tested: :math:`m \\in \\{0.1,\\ 0.5,\\ 1.0\\}` (dimensionless).
    """
    print("=" * 64)
    print("  Mass Sweep — Ground-State Energy vs. Particle Mass")
    print("=" * 64)

    masses: list[float] = [0.1, 0.5, 1.0]
    numerics = NumericalConfig(n_points=600, convergence_tolerance=1e-4)
    plotter  = QuantumPlotter(figsize=(9, 5), dpi=150)

    print(f"\n{'Mass':>8}  {'E_0':>12}  {'E_1':>12}  {'E_2':>12}  {'ratio E_1/E_0':>14}")
    print("─" * 62)

    systems: list[QuantumSystem] = []
    for mass in masses:
        # Frozen config: each iteration creates a new immutable object.
        # Demonstrating the "replace pattern" — we vary only `mass`;
        # all other constants keep their default values.
        physics = PhysicsConfig(
            hbar=1.0,
            mass=mass,
            kb=1.0,
            well_width=math.pi,
        )
        potential = FiniteSquareWell(v0=50.0, config=physics)
        solver    = NumerovSolver(numerics)
        system    = solver.solve(potential, n_states=5)
        systems.append(system)

        E = system.energies
        ratio = E[1] / E[0] if E[0] != 0.0 else float("inf")
        print(f"{mass:>8.1f}  {E[0]:>12.6f}  {E[1]:>12.6f}  {E[2]:>12.6f}  {ratio:>14.4f}")

    # Physical interpretation: E_0 ∝ ħ²/(2mL²), so halving the mass
    # should double E_0.  The table confirms this directly.
    print()
    E0_light = systems[0].ground_state_energy   # m = 0.1
    E0_mid   = systems[1].ground_state_energy   # m = 0.5
    E0_heavy = systems[2].ground_state_energy   # m = 1.0
    ratio_01_05 = E0_light / E0_mid
    ratio_05_10 = E0_mid   / E0_heavy
    print(f"  E_0(m=0.1) / E_0(m=0.5) = {ratio_01_05:.4f}  (ISW theory: {0.5/0.1:.1f}; deviation → finite-barrier tunnelling)")
    print(f"  E_0(m=0.5) / E_0(m=1.0) = {ratio_05_10:.4f}  (ISW theory: {1.0/0.5:.1f}; deviation → finite-barrier tunnelling)")

    # Save overlaid wavefunction plot for all three masses
    # (one per mass, side by side in separate files for clarity)
    for mass, system in zip(masses, systems):
        out = OUTPUT_DIR / f"sweep_wavefunctions_m{mass:.1f}.png"
        plotter.plot_wavefunctions(
            system, n_states=3,
            save_path=out,
        )
    print()
    print("  Wavefunction plots saved to research_output/figures/sweep_wavefunctions_m*.png")
    print("\n✓ Mass sweep complete.\n")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    run_main_simulation()
    run_mass_sweep()
