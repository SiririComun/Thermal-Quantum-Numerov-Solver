"""Particle-statistics enumeration for the two-particle thermal engine.

This module defines :class:`ParticleType`, an enumeration whose members
select the quantum-mechanical symmetry requirement applied to the two-particle
spatial wavefunction inside :class:`~src.core.statistics.ThermalEngine`.

Background
----------
The **Symmetry Postulate** of quantum mechanics (Pauli, 1940) states that the
total wavefunction of a system of *identical* particles must be either
*completely symmetric* (bosons) or *completely antisymmetric* (fermions) under
the exchange of any two particles.  Because the total wavefunction factorises
into a *spatial* part and a *spin* part, different spin-multiplet populations
impose different symmetry requirements on the spatial wavefunction.

The five members of :class:`ParticleType` cover the four cases implemented in
``legacy/research_prototype.ipynb`` plus a classical reference:

.. list-table::
   :header-rows: 1
   :widths: 14 10 30 46

   * - Member
     - Spin
     - Spatial symmetry
     - Physical origin
   * - BOLTZMANN
     - —
     - None (distinguishable)
     - Classical / Maxwell–Boltzmann
   * - BOSON
     - integer
     - Symmetric (:math:`+`)
     - Integer total spin (e.g. :sup:`4` He)
   * - FERMION
     - half-int
     - Antisymmetric (:math:`-`)
     - Half-integer spin, spin-polarised
   * - SPIN_HALF
     - 1/2
     - :math:`\tfrac{1}{4} P_S + \tfrac{3}{4} P_A`
     - 2-electron mixed-spin state (singlet + triplet)
   * - SPIN_ONE
     - 1
     - :math:`\tfrac{2}{3} P_S + \tfrac{1}{3} P_A`
     - 2-boson mixed-spin (spin-1); S=2,0 symmetric, S=1 antisym

**SPIN_HALF derivation** (spin-1/2 fermions, fully mixed spin state):
The two-particle spin Hilbert space has dimension 4.  The singlet state
(:math:`S=0`, antisymmetric spin) pairs with a *symmetric* spatial state so
that the total wavefunction is antisymmetric; the triplet (:math:`S=1`,
3 symmetric spin states) pairs with an *antisymmetric* spatial state.
At infinite temperature, all four spin states are equally likely, giving:

.. math::

    P_{\\text{SPIN-1/2}} = \\frac{1}{4} P_S + \\frac{3}{4} P_A

**SPIN_ONE derivation** (spin-1 bosons, fully mixed spin state):
The spin Hilbert space has dimension 9.  States with total spin
:math:`S=2` (5 states, symmetric spin) and :math:`S=0` (1 state, symmetric
spin) require symmetric spatial parts; :math:`S=1` (3 states, antisymmetric
spin) requires an antisymmetric spatial part.  Hence:

.. math::

    P_{\\text{SPIN-1}} = \\frac{6}{9} P_S + \\frac{3}{9} P_A
                       = \\frac{2}{3} P_S + \\frac{1}{3} P_A
"""

from __future__ import annotations

from enum import Enum, auto


class ParticleType(Enum):
    """Enumeration of two-particle statistics / symmetry types.

    Each member controls how :class:`~src.core.statistics.ThermalEngine`
    symmetrises the two-particle spatial wavefunction.

    Members:
        BOLTZMANN: Classical, distinguishable particles.  No exchange
            symmetry is imposed.  The two-particle density is the direct
            (outer) product of two independent single-particle densities.
            This serves as the classical reference for comparing quantum
            bunching / antibunching effects.
        BOSON: Identical bosons (integer total spin).  The spatial
            wavefunction is fully *symmetric* under exchange:
            :math:`\\Psi_S \\propto \\phi_n(x_1)\\phi_m(x_2)
            + \\phi_n(x_2)\\phi_m(x_1)`.
            Bosons tend to *bunch* (enhanced probability near the diagonal
            :math:`x_1 = x_2`).
        FERMION: Identical fermions in a spin-polarised state (or any
            state with a definite antisymmetric spin part).  The spatial
            wavefunction is fully *antisymmetric* under exchange:
            :math:`\\Psi_A \\propto \\phi_n(x_1)\\phi_m(x_2)
            - \\phi_n(x_2)\\phi_m(x_1)`.
            The Pauli Exclusion Principle is a direct consequence:
            :math:`\\Psi_A(x, x) = 0` exactly, so the diagonal of the
            two-particle density vanishes identically.
        SPIN_HALF: Two spin-1/2 particles in a thermally mixed spin state
            (equal population of all four spin substates).  Results in a
            spatial mixture of :math:`1/4\\, P_S + 3/4\\, P_A`.
        SPIN_ONE: Two spin-1 particles in a thermally mixed spin state
            (equal population of all nine spin substates).  Results in a
            spatial mixture of :math:`2/3\\, P_S + 1/3\\, P_A`.

    Example:
        >>> from src.models.statistics import ParticleType
        >>> ParticleType.SPIN_HALF
        <ParticleType.SPIN_HALF: 4>
        >>> ParticleType.SPIN_HALF.name
        'SPIN_HALF'
    """

    BOLTZMANN = auto()  # classical — no exchange symmetry
    BOSON = auto()      # fully symmetric spatial wavefunction
    FERMION = auto()    # fully antisymmetric spatial wavefunction
    SPIN_HALF = auto()  # 1/4 P_S + 3/4 P_A  (two spin-1/2 fermions)
    SPIN_ONE = auto()   # 2/3 P_S + 1/3 P_A  (two spin-1 bosons)
