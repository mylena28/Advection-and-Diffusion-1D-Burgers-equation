# Numerical Solution of Burgers' Equation

This repository contains Fortran implementations of numerical methods applied to fluid flow problems. 

The project includes:
- Analytical formulations
- Explicit numerical schemes
- A main driver code integrating the approaches

The focus is on understanding the behavior of fluid systems and validating numerical strategies for simplified flow models.

---

## 1D Burgers Equation: Advection–Diffusion Problem

The numerical solution of the one-dimensional Burgers' equation, is a nonlinear partial differential equation that combines advection and diffusion effects.

The equation is given by:

$$
\frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x} = \nu \frac{\partial^2 u}{\partial x^2}
$$

where:
- \( u(x,t) \) is the velocity field
- \( \nu \) is the kinematic viscosity

This equation serves as a simplified model for studying nonlinear convection and diffusion phenomena in fluid dynamics.

---

## How to run

Compile using a Fortran compiler (e.g., ifort, ifx, gfortran):

```bash
gfortran main.f95 analitica.f95 explicita.f95 -o Burgers
```

```bash
./Burgers
```

---

## Objective

The objective of this project is to compare analytical and numerical solutions for fluid flow problems, providing insight into the accuracy and stability of explicit schemes.
