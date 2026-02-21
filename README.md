# Thermal-Quantum-Numerov-Solver

[![Python Version](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![University](https://img.shields.io/badge/University-UdeA-green)](https://www.udea.edu.co/)
[![Documentation](https://img.shields.io/badge/docs-GitHub_Pages-brightgreen)](https://SiririComun.github.io/Thermal-Quantum-Numerov-Solver/)

**[Leer en Español](README.es.md)** | **[Explore Live Documentation](https://SiririComun.github.io/Thermal-Quantum-Numerov-Solver/)**

A high-performance computational physics engine designed to solve the **1D Time-Independent Schrödinger Equation** and simulate **Identical Particle Thermalization** in a Canonical Ensemble.

---

## 🎯 Project Overview

This project was originally developed for the **Statistical Physics** course at the **Universidad de Antioquia (UdeA)**. It implements a high-order numerical solver to explore how quantum particles behave under different potential constraints and thermal conditions.

### Scientific Core:
- **Numerical Method:** Matrix Numerov (**$O(dx^4)$ global accuracy**) for solving Eigenvalue problems.
- **Statistical Framework:** Canonical Ensemble (Boltzmann distribution) to calculate thermal densities.
- **Identical Particles:** Spatial symmetrization and antisymmetrization for Bosons and Fermions.
- **Spin Statistics:** Implementation of Spin-1/2 (Fermionic) and Spin-1 (Bosonic) mixtures.
- **Potentials:** Infinite/Finite Square Wells, Truncated Harmonic, and V-Shaped (Linear) Wells.

---

## 🛠 The Refactoring Journey (Software Architecture Focus)

The primary goal of this repository is to demonstrate the transition from **Scientific Prototyping** to **Professional Software Architecture**. I have refactored a legacy procedural research notebook into a modular, production-grade Python library.

### Architectural Improvements:
- **Procedural ➔ OOP:** Encapsulating physics entities (Potentials, Solvers, States) into clean class hierarchies.
- **Abstraction & Polymorphism:** Using **Abstract Base Classes (ABCs)** to define "contracts" for both Potentials and Solvers, making the engine algorithm-agnostic.
- **Dependency Inversion (DIP):** Decoupling the high-level simulation logic from low-level constants and specific numerical implementations.
- **Immutability:** Implementing **Frozen Dataclasses** for physical and numerical configurations to ensure scientific reproducibility.
- **Documentation as Code (DaC):** Using **Sphinx with Napoleon** to generate a professional API website that renders LaTeX equations via MathJax.

### 🗺️ Refactoring Roadmap

| Phase | Deliverable | Status |
|---|---|:---:|
| 1 | `PhysicsConfig` & `NumericalConfig` frozen dataclasses | ✅ Complete |
| 2 | `BasePotential` ABC + 4 concrete potentials | ✅ Complete |
| 3 | `NumerovSolver` — Matrix Numerov, O(dx⁴) validated | ✅ Complete |
| 3.5 | `BaseSolver` ABC — Dependency Inversion layer | ✅ Complete |
| 4 | `QuantumSystem` domain model — immutable result container | ✅ Complete |
| 5 | `ParticleType` enum + `ThermalEngine` — Pauli exclusion verified | ✅ Complete |
| 6 | `QuantumPlotter` — 4 publication-quality plot methods | ✅ Complete |
| 6.5 | Quality pass: ylim fix, visual validation suite | ✅ Complete |
| 7 | `run_simulation.py` master pipeline + Sphinx docs (0 warnings) | ✅ Complete |
| 8 | Showcase Jupyter Notebook | 🔜 Next |

### 🔬 Research Showcase

All library results are cross-validated against the original `legacy/research_prototype.ipynb` notebook to ensure numerical equivalence. The `run_simulation.py` entry point reproduces the complete pipeline — from eigenvalue decomposition to thermal pair densities — and saves publication-ready figures to `research_output/figures/`.

| Validation check | Result |
|---|:---:|
| Infinite-well spectrum: $E_n = n^2$ (dimensionless units) | ✅ Verified |
| Thermal density normalization: $\int \rho(x)\,dx = 1$ | ✅ Verified |
| Pauli exclusion: diagonal of Fermion pair density $= 0$ | ✅ Verified |
| Matrix Numerov convergence order $O(dx^4)$ | ✅ Verified |
| Mass scaling $E_0 \propto \hbar^2 / (2mL^2)$ (finite-barrier regime) | ✅ Verified |

### 🤖 AI-Assisted Architecture
This project utilizes a custom **Senior Research Architect** agent (configured in `.github/agents/`). This agent is designed to enforce **SOLID principles**, strict **Type Hinting**, and **Google-style documentation**, ensuring the transition from research to production maintains the highest engineering standards.

---

## 📂 Repository Structure

```text
├── legacy/                  # Original procedural research notebook (The Baseline)
├── research_output/         # Academic Paper (LaTeX), original PDFs, and presentations
│   └── figures/             # ✨ Pipeline-generated publication figures
├── docs/                    # [Deployed] Compiled HTML documentation website
├── docs_site/               # Source files for the Sphinx documentation engine
├── src/                     # Modular OOP Library
│   ├── core/                # Physics Engine (BaseSolver, Numerov, ThermalEngine)
│   ├── models/              # Data Models (Potentials, Configs, States, Statistics)
│   └── visualization/       # Professional plotting and 3D rendering utilities
├── test/
│   └── visual_validation/   # Headless plot validation outputs
├── run_simulation.py        # ✨ Master pipeline entry point
├── .github/agents/          # Custom AI Architect Agent configurations
├── .gitignore               # Multi-language (Python/LaTeX) hygiene rules
├── README.md                # Project documentation (English)
├── README.es.md             # Project documentation (Spanish)
└── requirements.txt         # Reproducibility manifest
```

---

## 🚀 Getting Started

### Prerequisites
- Python 3.10+
- Virtual Environment (recommended)

### Installation & Build
1. Clone the repository and setup environment:
   ```bash
   git clone https://github.com/SiririComun/Thermal-Quantum-Numerov-Solver.git
   cd Thermal-Quantum-Numerov-Solver
   python -m venv .venv
   source .venv/bin/activate  # Or .venv\Scripts\activate on Windows
   pip install -r requirements.txt
   ```
2. Build the documentation locally:
   ```bash
   sphinx-build -b html docs_site/source docs
   ```
3. Run the master simulation pipeline:
   ```bash
   python run_simulation.py
   ```
   Figures are saved to `research_output/figures/`.

---

## 🎓 Credits

- **Pablo Sanchez** - [@SiririComun](https://github.com/SiririComun)
- **Juan Montoya** - [@Juanj27](https://github.com/Juanj27)

Developed as part of the Statistical Physics course at the **Universidad de Antioquia (UdeA)**, Medellín, Colombia. Refactored as a personal project to bridge the gap between Computational Physics and Professional Software Engineering.
