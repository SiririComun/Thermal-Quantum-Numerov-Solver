# Thermal-Quantum-Numerov-Solver

**[Leer en Español](README.es.md)**

[![Python Version](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![University](https://img.shields.io/badge/University-UdeA-green)](https://www.udea.edu.co/)

A high-performance computational physics engine designed to solve the **1D Time-Independent Schrödinger Equation** and simulate **Identical Particle Thermalization** in a Canonical Ensemble.

---

## 🎯 Project Overview

This project was originally developed for the **Statistical Physics** course at the **Universidad de Antioquia (UdeA)**. It implements a high-order numerical solver to explore how quantum particles behave under different potential constraints and thermal conditions.

### Scientific Core:
- **Numerical Method:** Matrix Numerov (O(dx⁴) accuracy) for solving Eigenvalue problems.
- **Statistical Framework:** Canonical Ensemble (Boltzmann distribution) to calculate thermal densities.
- **Identical Particles:** Spatial symmetrization and antisymmetrization for Bosons and Fermions.
- **Spin Statistics:** Implementation of Spin-1/2 (Fermionic) and Spin-1 (Bosonic) mixtures.
- **Potentials:** Infinite/Finite Square Wells, Truncated Harmonic, and V-Shaped (Linear) Wells.

---

## 🛠 The Refactoring Journey (Software Architecture Focus)

The primary goal of this repository is to demonstrate the transition from **Scientific Prototyping** to **Professional Software Architecture**. 

I am currently refactoring the original procedural notebook into a modular, script-based Python library, applying **SOLID principles**, **OOP Design Patterns**, and industry-standard practices.

### Architectural Improvements:
- **Procedural ➔ OOP:** Encapsulating physics entities (Potentials, Solvers, States) into clean class hierarchies.
- **Dependency Injection:** Decoupling the Physics Solver from specific potential implementations.
- **Type Safety:** Full integration of Python Type Hinting for robust scientific computing.
- **Modularity:** Strict separation of concerns between the Physics Engine, Data Models, and Visualization logic.
- **Maintainability:** Transitioning from "linear spaghetti" notebooks to an importable Python package.

---

## 📂 Repository Structure

```text
├── legacy/              # Original procedural research notebook (The Baseline)
├── docs/                # Scientific Paper (LaTeX) and Class Presentation
├── src/                 # [In Progress] New modular OOP Library
│   ├── core/            # Physics solvers and mathematical logic
│   ├── models/          # Data classes and physics entities (Potentials, States)
│   └── visualization/   # Professional plotting and 3D rendering utilities
├── tests/               # Unit tests for physics validation
├── README.md            # Project documentation
└── requirements.txt     # Reproducibility manifest

```
---

## 🚀 Getting Started

### Prerequisites
- Python 3.10+
- Virtual Environment (recommended)

### Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/SiririComun/HotBox-Quantum-Sim.git
   cd HotBox-Quantum-Sim
   ```
2. Create and activate a virtual environment:
   ```bash
   python -m venv .venv
   # Windows:
   source .venv/Scripts/activate
   # macOS/Linux:
   source source .venv/bin/activate
   ```
3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
---

## 🎓 Credits

- **Pablo Sanchez** - [@SiririComun](https://github.com/SiririComun)
- **Juan Montoya** - [@Juanj27](https://github.com/Juanj27)

Developed as part of the Statistical Physics course at the **Universidad de Antioquia (UdeA)**, Medellín, Colombia. Refactored as a personal project to bridge the gap between Computational Physics and Professional Software Engineering.

## 🤖 AI-Assisted Architecture

This project utilizes a custom **Senior Research Architect** agent (configured in `.github/agents/`). This agent is designed to enforce **SOLID principles**, strict **Type Hinting**, and **Google-style documentation**, ensuring the transition from research to production maintains the highest engineering standards.
