Thermal-Quantum Numerov Solver
==============================

A production-grade Python library for solving the **1D Time-Independent
Schrödinger Equation** and simulating **identical-particle thermalization**
in the Canonical Ensemble — refactored from a legacy research notebook into
a SOLID-compliant, fully documented scientific library.

----

Key Features
------------

- **Matrix Numerov solver** — :math:`O(dx^4)` global accuracy via a
  generalised eigenvalue problem :math:`A\boldsymbol{\psi}=\varepsilon B\boldsymbol{\psi}`.
- **Four potential geometries** — Infinite/Finite Square Wells, Truncated
  Harmonic, and V-Shaped (linear) potentials.
- **Canonical Ensemble thermodynamics** — Boltzmann-weighted single-particle
  and two-particle densities with configurable temperature.
- **All five particle statistics** — Boltzmann, Boson, Fermion, Spin-1/2,
  and Spin-1 two-particle symmetry rules; Pauli exclusion enforced exactly.
- **Publication-quality plots** — energy-level diagrams, thermal density
  curves, 2-D heatmaps, and 3-D surface renderings via ``QuantumPlotter``.
- **Immutable configuration** — frozen dataclasses guarantee reproducibility
  across every solver and plotter call.

----

Quick Start
-----------

.. code-block:: python

   from src.models import FiniteSquareWell
   from src.core  import NumerovSolver, ThermalEngine
   from src.models.statistics import ParticleType
   from src.visualization import QuantumPlotter

   system  = NumerovSolver().solve(FiniteSquareWell(v0=50.0), n_states=20)
   rho     = ThermalEngine().calculate_thermal_density(system, temperature=5.0)
   rho2p   = ThermalEngine().calculate_pair_density(system, 5.0, ParticleType.FERMION)

   plotter = QuantumPlotter(figsize=(9, 5), dpi=150)
   plotter.plot_wavefunctions(system, n_states=5)
   plotter.plot_pair_density_heatmap(system, rho2p, title="Fermion pair density")

----

API Reference
-------------

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: API Reference

   config_api
   physics_api
   solvers_api
   viz_api

.. grid:: 2

   .. grid-item-card:: ⚙️ Configuration
      :link: config_api
      :link-type: doc

      Immutable frozen dataclasses — ``PhysicsConfig`` and
      ``NumericalConfig``.

   .. grid-item-card:: 🔭 Physics Models
      :link: physics_api
      :link-type: doc

      Potential geometries, ``QuantumSystem`` result container, and
      ``ParticleType`` enum.

.. grid:: 2

   .. grid-item-card:: 🧮 Solvers & Statistics
      :link: solvers_api
      :link-type: doc

      ``NumerovSolver`` (Matrix Numerov, :math:`O(dx^4)`) and
      ``ThermalEngine`` (Canonical Ensemble, pair densities).

   .. grid-item-card:: 📊 Visualization
      :link: viz_api
      :link-type: doc

      ``QuantumPlotter`` — four publication-quality plot methods with
      headless ``save_path`` support.

