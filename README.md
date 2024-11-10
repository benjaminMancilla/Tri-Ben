# Proyecto de Triangulación Delaunay

Este proyecto implementa una triangulación de Delaunay utilizando C++ y SFML para la visualización gráfica. El objetivo es construir un programa que maneje triangulaciones y permita visualizar los resultados de manera interactiva.

## Requisitos

Antes de ejecutar el proyecto, asegúrate de tener instaladas las siguientes dependencias:

- **C++14** o superior.
- **CMake** 3.10 o superior.
- **SFML** (para gráficos y manejo de ventanas).
  - `libsfml-graphics`
  - `libsfml-window`
  - `libsfml-system`

En sistemas basados en Debian (Ubuntu, etc.), puedes instalar SFML con el siguiente comando:

```bash
sudo apt-get install libsfml-dev
```

En este archivo `README.md`, se incluyen los siguientes pasos adicionales para compilar y ejecutar el proyecto:

1. **Configurar el proyecto con CMake**: `cmake ..`
2. **Compilar el proyecto con `make`**.
3. **Ejecutar el programa**: `./tri-ben.exe`

Este flujo es el típico para proyectos que utilizan CMake y que dependen de bibliotecas externas como SFML.


