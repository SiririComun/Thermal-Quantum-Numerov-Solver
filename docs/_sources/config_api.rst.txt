Configuration Layer
===================

Immutable frozen dataclasses that encode the dimensionless unit system
used throughout the library.

.. rubric:: Design Principle

``PhysicsConfig`` and ``NumericalConfig`` are **Frozen Dataclasses**
(``@dataclass(frozen=True, slots=True)``).  Once instantiated they cannot be
mutated, which guarantees that every call to a solver or plotter operates on
a self-consistent, reproducible set of constants.  This is the foundation of
the **Immutability** pillar in the SOLID architecture.

.. rubric:: Dimensionless Units

All default values are chosen so that the 1D infinite square-well reference
spectrum satisfies :math:`E_n = n^2`:

.. math::

   E_n = \frac{n^2 \pi^2 \hbar^2}{2 m L^2}
       \xrightarrow{\hbar=1,\, m=0.5,\, L=\pi} n^2

----

.. automodule:: src.models.config
   :members:
   :undoc-members:
   :show-inheritance:
