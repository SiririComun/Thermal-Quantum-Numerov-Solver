Physics Models
==============

This page documents the **potential geometry** and the **quantum-state
container** — the two model objects that flow through the entire pipeline.

Potential Definitions
---------------------

``BasePotential`` is the **Abstract Base Class** that defines the interface
contract every potential must satisfy.  The ``NumerovSolver`` depends
*only* on this abstract type, never on any concrete subclass — a strict
application of the **Dependency Inversion Principle**.

The four concrete potentials faithfully reproduce the wells defined in
``legacy/research_prototype.ipynb``:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Class
     - Physics
   * - ``InfiniteSquareWell``
     - Hard-wall box.  Analytical spectrum: :math:`E_n = n^2`.
   * - ``FiniteSquareWell``
     - Flat-bottom well of depth :math:`V_0`.  Evanescent decay into barrier.
   * - ``HarmonicPotential``
     - Truncated parabola :math:`V = \tfrac{1}{2}m\omega^2 x^2`.
   * - ``VShapedPotential``
     - Truncated linear (V-shape) :math:`V = \lambda|x - L/2|`.

.. automodule:: src.models.potentials
   :members:
   :undoc-members:
   :show-inheritance:

----

Quantum State Container
-----------------------

``QuantumSystem`` is the **read-only result object** returned by every solver.
It wraps the raw NumPy eigenpairs and the originating potential into a
well-typed, immutable container.  No physics logic lives here — it is pure
data, accessible through properties only.

.. automodule:: src.models.states
   :members:
   :undoc-members:
   :show-inheritance:

----

Particle Statistics Enumeration
--------------------------------

``ParticleType`` selects the quantum-mechanical symmetry requirement applied
to the two-particle spatial wavefunction inside ``ThermalEngine``.

.. automodule:: src.models.statistics
   :members:
   :undoc-members:
   :show-inheritance:
