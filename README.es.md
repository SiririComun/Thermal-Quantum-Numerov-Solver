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

### 🤖 Arquitectura Asistida por IA
Este proyecto emplea un agente personalizado de **Arquitecto de Investigación Senior** (configurado en `.github/agents/`). Este agente supervisa el cumplimiento de los **principios SOLID**, el uso estricto de **Type Hinting** (tipado estático) y la **documentación bajo el estándar de Google**, asegurando que el código mantenga estándares de ingeniería de software de alto nivel.

---

## 📂 Estructura del Repositorio

```text
├── legacy/              # Código original del proyecto (Línea base)
├── research_output/     # Artículo académico (LaTeX), PDFs originales y presentaciones
├── docs/                # [Publicado] Sitio web de documentación (HTML compilado)
├── docs_site/           # Archivos fuente del motor de documentación Sphinx
├── src/                 # Librería modular en POO
│   ├── core/            # Motor de cálculo (BaseSolver, Motor de Numerov)
│   ├── models/          # Modelos de datos (Potenciales, Configs, Estados)
│   └── visualization/   # Utilidades de graficación y renderizado 3D
├── .github/agents/      # Configuraciones del Agente Arquitecto (IA)
├── .gitignore           # Reglas de limpieza para Python y LaTeX
├── README.md            # Documentación principal (Inglés)
├── README.es.md         # Documentación principal (Español)
└── requirements.txt     # Manifiesto de dependencias para reproducibilidad
```

---

## 🚀 Guía de Inicio

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

---

## 🎓 Créditos

- **Pablo Sanchez** - [@SiririComun](https://github.com/SiririComun)
- **Juan Montoya** - [@Juanj27](https://github.com/Juanj27)

Proyecto desarrollado para el curso de Física Estadística en la **Universidad de Antioquia (UdeA)**, Medellín, Colombia. Refactorizado para cerrar la brecha entre la Física Computacional y la Ingeniería de Software Profesional.