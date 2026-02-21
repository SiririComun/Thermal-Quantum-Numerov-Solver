Visualization
=============

``QuantumPlotter`` provides four publication-quality plot methods built on
top of Matplotlib.  Every method returns a ``(Figure, Axes)`` tuple and
accepts an optional ``save_path`` argument — making it suitable for both
interactive notebooks and headless CI pipelines.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Method
     - Output
   * - ``plot_wavefunctions``
     - Energy-level diagram with probability densities overlaid on the
       potential.
   * - ``plot_thermal_density``
     - Single-particle Boltzmann-weighted density :math:`\rho(x)` as a
       filled curve on top of the eigenstate ladder.
   * - ``plot_pair_density_heatmap``
     - Two-particle pair density :math:`\rho^{(2)}(x_1,x_2)` as a 2-D
       colour map.  Fermion diagonal suppression is clearly visible.
   * - ``plot_pair_density_3d``
     - Same data rendered as a 3-D surface plot using ``Axes3D``.

.. rubric:: Usage Example

.. code-block:: python

   from src.models import FiniteSquareWell, PhysicsConfig
   from src.core import NumerovSolver, ThermalEngine
   from src.models.statistics import ParticleType
   from src.visualization import QuantumPlotter
   from pathlib import Path

   system  = NumerovSolver().solve(FiniteSquareWell(v0=50.0), n_states=20)
   engine  = ThermalEngine()
   rho     = engine.calculate_thermal_density(system, temperature=5.0)
   rho2p   = engine.calculate_pair_density(system, 5.0, ParticleType.FERMION)

   plotter = QuantumPlotter(figsize=(9, 5), dpi=150)
   plotter.plot_wavefunctions(system, n_states=5, save_path=Path("wf.png"))
   plotter.plot_pair_density_heatmap(system, rho2p, title="Fermion pair density")

----

.. automodule:: src.visualization.plotters
   :members:
   :undoc-members:
   :show-inheritance:
