Solvers & Statistics Engine
============================

This page documents the **numerical core** of the library: the eigenvalue
solver and the thermal statistics engine.

Matrix Numerov Solver
---------------------

``NumerovSolver`` converts the 1D TISE into a generalised eigenvalue problem
of order :math:`O(dx^4)` using the **Matrix Numerov** method.

.. rubric:: Algorithm Summary

Starting from the three-point Numerov recurrence, the interior grid points
are assembled into two symmetric tridiagonal matrices :math:`A` and :math:`B`:

.. math::

   A\,\boldsymbol{\psi} = \varepsilon\, B\,\boldsymbol{\psi}

The system is solved by ``scipy.linalg.eigh``, which exploits the
real-symmetric positive-definite structure of :math:`B` to return eigenvalues
in ascending order without post-sorting.

.. rubric:: Convergence

Numerov's method achieves :math:`O(dx^4)` global accuracy.  With the default
600-point grid, the infinite-well ground state recovers :math:`E_0 = 1.0000`
to six significant figures.

.. automodule:: src.core.solvers
   :members:
   :undoc-members:
   :show-inheritance:

----

Thermal Statistics Engine
--------------------------

``ThermalEngine`` computes single-particle and two-particle thermal densities
in the **Canonical Ensemble** (Boltzmann distribution).

.. rubric:: Single-Particle Density

.. math::

   \rho(x) = \frac{1}{Z} \sum_{n=0}^{N}
             e^{-E_n / k_B T}\,|\psi_n(x)|^2,
   \qquad Z = \sum_{n=0}^{N} e^{-E_n / k_B T}

The series is truncated when successive terms fall below
``NumericalConfig.convergence_tolerance``.

.. rubric:: Two-Particle Pair Density

For identical particles, the **Symmetry Postulate** requires:

.. math::

   \Psi^{\pm}_{nm}(x_1, x_2)
   = \frac{1}{\sqrt{2}}
     \bigl[\psi_n(x_1)\psi_m(x_2)
     \pm \psi_m(x_1)\psi_n(x_2)\bigr]

The ``−`` branch enforces **Pauli exclusion**: :math:`\rho^{(2)}(x,x) \equiv 0`
for all fermion configurations.

.. automodule:: src.core.statistics
   :members:
   :undoc-members:
   :show-inheritance:
