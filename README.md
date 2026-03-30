This repository contains Fortran implementations of numerical methods applied to fluid flow problems. 

The project includes:
- Analytical formulations
- Explicit numerical schemes
- A main driver code integrating the approaches

The focus is on understanding the behavior of fluid systems and validating numerical strategies for simplified flow models.

---

## How to run

Compile using a Fortran compiler (e.g., ifort, ifx, gfortran):

```bash
gfortran main.f95 analitica.f95 explicita.f95 -o Burgers
./Burgers
```

---

## Objective

The objective of this project is to compare analytical and numerical solutions for fluid flow problems, providing insight into the accuracy and stability of explicit schemes.
