"""Professional visualization package for the Numerov quantum library.

Provides :class:`~src.visualization.plotters.QuantumPlotter`, a stateless
rendering service that consumes :class:`~src.models.states.QuantumSystem`
and pre-computed density arrays from :class:`~src.core.statistics.ThermalEngine`.

Separation of Concerns
----------------------
The visualisation layer is deliberately isolated from the physics core:

- It **never** calls a solver or thermal engine internally.
- All data arrives as function arguments (Dependency Injection).
- The plotting backend (``matplotlib``) is a detail of this package only;
  no other package imports it.

This makes the library extensible to alternative front-ends (interactive
Plotly dashboards, headless PDF pipelines, web APIs) by replacing only
this package — the physics core is untouched.
"""

from .plotters import QuantumPlotter

__all__ = ["QuantumPlotter"]
