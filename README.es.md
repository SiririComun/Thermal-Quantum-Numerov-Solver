# Thermal-Quantum-Numerov-Solver

[![Versión de Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/downloads/)
[![Licencia: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Universidad](https://img.shields.io/badge/Universidad-UdeA-green)](https://www.udea.edu.co/)

Un motor de física computacional de alto rendimiento diseñado para resolver la **Ecuación de Schrödinger Independiente del Tiempo en 1D** y simular la **Termalización de Partículas Idénticas** en un Ensamble Canónico.

---

## 🎯 Descripción del Proyecto

Este proyecto fue desarrollado originalmente para el curso de **Física Estadística** en la **Universidad de Antioquia (UdeA)**. Implementa un solvente numérico de alto orden para explorar el comportamiento de partículas cuánticas bajo diferentes restricciones de potencial y condiciones térmicas.

### Núcleo Científico:
- **Método Numérico:** Numerov Matricial (precisión de orden O(dx⁴)) para resolver problemas de autovalores.
- **Marco Estadístico:** Ensamble Canónico (distribución de Boltzmann) para el cálculo de densidades térmicas.
- **Partículas Idénticas:** Simetrización y antisimetrización espacial para Bosones y Fermiones.
- **Estadística de Espín:** Implementación de mezclas para Espín-1/2 (Fermiónico) y Espín-1 (Bosónico).
- **Potenciales:** Pozos cuadrados finitos/infinitos, armónico truncado y pozos en forma de V (lineales).

---

## 🛠 El Proceso de Refactorización (Enfoque en Arquitectura)

El objetivo principal de este repositorio es demostrar la transición del **Prototipado Científico** a la **Arquitectura de Software Profesional**.

Actualmente, estoy refactorizando el cuaderno (notebook) procedimental original hacia una librería de Python modular basada en scripts, aplicando **principios SOLID**, **Patrones de Diseño de POO** y prácticas estándar de la industria.

### Mejoras Arquitectónicas:
- **Procedimental ➔ POO:** Encapsulamiento de entidades físicas (Potenciales, Solventes, Estados) en jerarquías de clases limpias.
- **Inyección de Dependencias:** Desacoplamiento del solvente físico de las implementaciones específicas de los potenciales.
- **Seguridad de Tipos (Type Safety):** Integración total de *Type Hinting* para un cálculo científico robusto.
- **Modularidad:** Separación estricta de responsabilidades entre el motor físico, los modelos de datos y la lógica de visualización.
- **Mantenibilidad:** Transición de "código espagueti" lineal en cuadernos a un paquete de Python importable y organizado.

---

## 📂 Estructura del Repositorio

```text
├── legacy/              # Cuaderno de investigación original (Línea base)
├── docs/                # Artículo científico (LaTeX) y presentación de clase
├── src/                 # [En progreso] Nueva librería modular en POO
│   ├── core/            # Solventes físicos y lógica matemática
│   ├── models/          # Clases de datos y entidades (Potenciales, Estados)
│   └── visualization/   # Utilidades profesionales de graficación y renderizado 3D
├── tests/               # Pruebas unitarias para validación física
├── README.md            # Documentación en inglés
└── README.es.md         # Documentación en español
```

---

## 🚀 Guía de Inicio

### Requisitos previos
- Python 3.10+
- Entorno virtual (recomendado)

### Instalación
1. Clona el repositorio:
   ```bash
   git clone https://github.com/SiririComun/Thermal-Quantum-Numerov-Solver.git
   cd Thermal-Quantum-Numerov-Solver
   ```
2. Crea y activa un entorno virtual:
   ```bash
   python -m venv .venv
   # Windows:
   .venv\Scripts\activate
   # macOS/Linux:
   source .venv/bin/activate
   ```
3. Instala las dependencias:
   ```bash
   pip install -r requirements.txt
   ```

---

## 🎓 Créditos

- **Pablo Sánchez** - [@SiririComun](https://github.com/SiririComun)
- **Juan Montoya** - [@Juanj27](https://github.com/Juanj27)

Desarrollado como parte del curso de Física Estadística en la **Universidad de Antioquia (UdeA)**, Medellín, Colombia. Refactorizado como un proyecto personal para cerrar la brecha entre la Física Computacional y la Ingeniería de Software Profesional.

