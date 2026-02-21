# Thermal-Quantum-Numerov-Solver

[![Versión de Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/downloads/)
[![Licencia: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Universidad](https://img.shields.io/badge/University-UdeA-green)](https://www.udea.edu.co/)
[![Documentación](https://img.shields.io/badge/docs-GitHub_Pages-brightgreen)](https://SiririComun.github.io/Thermal-Quantum-Numerov-Solver/)

**[Read in English](README.md)** | **[Explorar Documentación en Línea](https://SiririComun.github.io/Thermal-Quantum-Numerov-Solver/)**

Motor de física computacional de alto rendimiento diseñado para resolver la **Ecuación de Schrödinger Independiente del Tiempo en 1D** y simular la **Termalización de Partículas Idénticas** en un Ensamble Canónico.

---

## 🎯 Descripción del Proyecto

Este proyecto nació originalmente como un trabajo para el curso de **Física Estadística** en la **Universidad de Antioquia (UdeA)**. Implementa un solucionador numérico de alto orden para explorar el comportamiento de partículas cuánticas bajo diversos potenciales y condiciones térmicas.

### Núcleo Científico:
- **Método Numérico:** Numerov Matricial (**precisión global $O(dx^4)$**) para la resolución de problemas de autovalores.
- **Marco Estadístico:** Ensamble Canónico (distribución de Boltzmann) para el cálculo de densidades térmicas.
- **Partículas Idénticas:** Simetrización y antisimetrización espacial para bosones y fermiones.
- **Estadística de Espín:** Implementación de mezclas para sistemas de Espín-1/2 (fermiónico) y Espín-1 (bosónico).
- **Potenciales:** Pozos cuadrados finitos e infinitos, potencial armónico truncado y pozos lineales (en forma de V).

---

## 🛠 El Proceso de Refactorización (Enfoque en Arquitectura)

El objetivo principal de este repositorio es demostrar la transición de un **Prototipo Científico** a una **Arquitectura de Software Profesional**. He transformado un cuaderno de investigación estructurado de forma procedimental en una librería de Python modular, robusta y lista para producción.

### Mejoras Arquitectónicas:
- **De lo Procedimental a POO:** Encapsulamiento de entidades físicas (potenciales, solucionadores, estados) en jerarquías de clases limpias.
- **Abstracción y Polimorfismo:** Uso de **Clases Base Abstractas (ABCs)** para definir los "contratos" de potenciales y solucionadores, permitiendo que el motor sea agnóstico al algoritmo utilizado.
- **Inversión de Dependencias (DIP):** Desacoplamiento de la lógica de simulación de alto nivel respecto a las constantes físicas y las implementaciones numéricas específicas.
- **Inmutabilidad:** Uso de **Dataclasses inmutables (frozen)** para las configuraciones físicas y numéricas, garantizando la reproducibilidad de los resultados.
- **Documentación como Código (DaC):** Implementación de **Sphinx con Napoleon** para generar un sitio web profesional de la API que renderiza ecuaciones complejas mediante MathJax.

### 🗺️ Hoja de Ruta del Proyecto

| Fase | Entregable | Estado |
|---|---|:---:|
| 1 | Dataclasses inmutables `PhysicsConfig` y `NumericalConfig` | ✅ Completo |
| 2 | ABC `BasePotential` + 4 potenciales concretos | ✅ Completo |
| 3 | `NumerovSolver` — Numerov Matricial, precisión O(dx⁴) validada | ✅ Completo |
| 3.5 | ABC `BaseSolver` — capa de Inversión de Dependencias | ✅ Completo |
| 4 | Modelo de dominio `QuantumSystem` — contenedor de resultados inmutable | ✅ Completo |
| 5 | Enum `ParticleType` + `ThermalEngine` — exclusión de Pauli verificada | ✅ Completo |
| 6 | `QuantumPlotter` — 4 métodos de gráficas de calidad para publicación | ✅ Completo |
| 6.5 | Pase de calidad: corrección de ylim, suite de validación visual | ✅ Completo |
| 7 | Pipeline maestro `run_simulation.py` + documentación Sphinx (0 advertencias) | ✅ Completo |
| 8 | `Showcase.ipynb` — cuaderno ejecutado con salidas incrustadas | ✅ Completo |
| 9 | Pulido Final y Empaquetado — README actualizado, build de documentación | ✅ Completo |

### 🔬 Vitrina de Investigación

Todos los resultados de la librería son validados en comparación con el cuaderno original `legacy/research_prototype.ipynb` para garantizar equivalencia numérica. El punto de entrada `run_simulation.py` reproduce el pipeline completo — desde la descomposición en valores propios hasta las densidades de pares térmicos — y guarda figuras listas para publicación en `research_output/figures/`.

| Verificación de validación | Resultado |
|---|:---:|
| Espectro del pozo infinito: $E_n = n^2$ (unidades adimensionales) | ✅ Verificado |
| Normalización de la densidad térmica: $\int \rho(x)\,dx = 1$ | ✅ Verificado |
| Exclusión de Pauli: diagonal de la densidad de pares fermiónica $= 0$ | ✅ Verificado |
| Orden de convergencia del Numerov Matricial $O(dx^4)$ | ✅ Verificado |
| Escalado de masa $E_0 \propto \hbar^2 / (2mL^2)$ (régimen de barrera finita) | ✅ Verificado |

### 📓 Cuaderno de Demostración

[`Showcase.ipynb`](Showcase.ipynb) (en la raíz del proyecto) es un Jupyter Notebook completamente ejecutado que demuestra el pipeline completo de la librería en **seis secciones autocontenidas**, con todas las figuras incrustadas como salida:

| Sección | Demuestra |
|:--------|:----------|
| §0 — Configuración del Entorno | Una sola celda de importación reemplaza ~20 líneas de constantes globales |
| §1 — Simulación en Tres Líneas | `NumerovSolver().solve(FiniteSquareWell(v0=50))` + diagrama de niveles de energía |
| §2 — Estadística Multi-Partícula | Mapa de calor de densidad de pares fermiónicos — exclusión de Pauli $\rho(x,x)=0$ |
| §3 — Agujero de Intercambio (3D) | Renderizado de superficie con `plot_pair_density_3d` |
| §4 — Barrido de Masa | Barrido paramétrico con `dataclasses.replace()` — 3 configuraciones, cero mutación |
| §5 — Métricas de Boilerplate | Gráfico de barras que evidencia **≈ 98 % de reducción de código** vs. el prototipo legado |

> **Hallazgo clave:** el `research_prototype.ipynb` legado requería ~240 líneas de código
> repetitivo para reproducir lo que esta librería expresa en ~5 llamadas a métodos.

### 🤖 Arquitectura Asistida por IA
Este proyecto emplea un agente personalizado de **Arquitecto de Investigación Senior** (configurado en `.github/agents/`). Este agente supervisa el cumplimiento de los **principios SOLID**, el uso estricto de **Type Hinting** (tipado estático) y la **documentación bajo el estándar de Google**, asegurando que el código mantenga estándares de ingeniería de software de alto nivel.

---

## 📂 Estructura del Repositorio

```text
├── legacy/                  # Código original del proyecto (Línea base)
├── research_output/         # Artículo académico (LaTeX), PDFs originales y presentaciones
│   └── figures/             # ✨ Figuras de publicación generadas por el pipeline
├── docs/                    # [Publicado] Sitio web de documentación (HTML compilado)
├── docs_site/               # Archivos fuente del motor de documentación Sphinx
├── src/                     # Librería modular en POO
│   ├── core/                # Motor de cálculo (BaseSolver, Motor de Numerov, ThermalEngine)
│   ├── models/              # Modelos de datos (Potenciales, Configs, Estados, Estadística)
│   └── visualization/       # Utilidades de graficación y renderizado 3D
├── test/
│   └── visual_validation/   # Salidas de validación visual sin cabecera
├── Showcase.ipynb           # ✨ Cuaderno de demostración completamente ejecutado
├── run_simulation.py        # ✨ Punto de entrada del pipeline maestro
├── .github/agents/          # Configuraciones del Agente Arquitecto (IA)
├── .gitignore               # Reglas de limpieza para Python y LaTeX
├── README.md                # Documentación principal (Inglés)
├── README.es.md             # Documentación principal (Español)
└── requirements.txt         # Manifiesto de dependencias para reproducibilidad
```

---

## � Documentación en Línea

La referencia completa de la API está publicada en **GitHub Pages**, generada desde docstrings estilo Google con Sphinx 9 y el tema PyData Sphinx:

> 🔗 **[https://SiririComun.github.io/Thermal-Quantum-Numerov-Solver/](https://SiririComun.github.io/Thermal-Quantum-Numerov-Solver/)**

| Página | Descripción |
|:-------|:------------|
| [Inicio](https://SiririComun.github.io/Thermal-Quantum-Numerov-Solver/) | Resumen del proyecto y guía rápida |
| [Config API](https://SiririComun.github.io/Thermal-Quantum-Numerov-Solver/config_api.html) | Dataclasses inmutables `PhysicsConfig` y `NumericalConfig` |
| [Physics API](https://SiririComun.github.io/Thermal-Quantum-Numerov-Solver/physics_api.html) | `BasePotential`, potenciales concretos, `ThermalEngine`, `ParticleType` |
| [Solvers API](https://SiririComun.github.io/Thermal-Quantum-Numerov-Solver/solvers_api.html) | `BaseSolver`, `NumerovSolver`, `QuantumSystem` |
| [Visualization API](https://SiririComun.github.io/Thermal-Quantum-Numerov-Solver/viz_api.html) | `QuantumPlotter` — los 4 métodos de graficación |

---

## �🚀 Guía de Inicio

### Requisitos previos
- Python 3.10+
- Entorno virtual (recomendado)

### Instalación y Compilación
1. Clona el repositorio y configura el entorno:
   ```bash
   git clone https://github.com/SiririComun/Thermal-Quantum-Numerov-Solver.git
   cd Thermal-Quantum-Numerov-Solver
   python -m venv .venv
   source .venv/bin/activate  # En Windows usa: .venv\Scripts\activate
   pip install -r requirements.txt
   ```
2. Genera la documentación localmente:
   ```bash
   sphinx-build -b html docs_site/source docs
   ```
3. Ejecuta el pipeline maestro de simulación:
   ```bash
   python run_simulation.py
   ```
   Las figuras se guardan en `research_output/figures/`.

---

## 🎓 Créditos

- **Pablo Sanchez** - [@SiririComun](https://github.com/SiririComun)
- **Juan Montoya** - [@Juanj27](https://github.com/Juanj27)

Proyecto desarrollado para el curso de Física Estadística en la **Universidad de Antioquia (UdeA)**, Medellín, Colombia. Refactorizado para cerrar la brecha entre la Física Computacional y la Ingeniería de Software Profesional.