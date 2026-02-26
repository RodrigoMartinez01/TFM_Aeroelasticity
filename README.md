# 🚀 Aeroelasticity: Final Thesis Project
This project is my final thesis for the Master's degree in Aeronautical Engineering. The main goal of the project is to study a method to solve a non linear aeroelastic system in the time domain, in order to accouunt for the effect that the potential unsteady drag can have on the behavior of the system. 

## Introduction
The potential unsteady drag, is a non linear aerodynamic force that was originally studied by I.E.Garrick in 1936. Classical methods of flutter analysis do not take into account this force, as it's direct effect on the torsion and bending of the wing is negligible. However, some authors have recently propose new mathods to account for the unsteady drag in the aeroelastic system. 

Some of those include the coupling of the horizontal dynamic equation of the aircraft (i.e. thrust equation) with the aeroelastic system by allowing changes in the freestrem velocity, or by including a new degree of freedom to model the in-plane bending of the wing, usually called surge, so that the unsteady drag can directly influence the aeroelastic response.

The reason why in this repository there are codes for easier cases than the last ones mentioned, is to propose a step by step approach to the solution of the problem, and allow a careful validation of each of the theories and tools implemented throughout the project. 

Moreover, this way the user can easely follow the logical evolution of the project from the begining, to the end of it. 
The matlab scripts are intended to justify every expression used in each code, and to give a more academic perspective of the mathematical developement of the project. 

---

## Brief explanation of each Matlab code

The explanations of the codes is organized following the same time evolution that I have implemented. I strongly recommend to follow this order, as each code is use to understand some tools or phenomena that may be used in the next codes. It is the most logical approach to the problem, from easier and known solutions and methods, to more complex ones.

### 1. Quasi-steady Aerodynamic, K, and PK Flutter codes

These 3 codes are basic implementations of a 2 DoF (degrees of freedom) flutter solver using some of the most classical linear flutter analysis. They will be used lately to validate state-space representations to solve the aeroelastic system evolution in time domain.

### 2. Unsteady Aerodynamic Analysis

The main purpose of this code is to study in detail the different mathematical theories used to model the unsteady aerodynamic forces (lift, moment, and drag). The plots compare, for each aerodynamic load, the time domain representation given by each theory. This different methods are:
1) Quasi-steady theory
2) Theodorsen's and Garrick theory
3) R.T.Jones approximation of Wagner indicial function, and Convolution Integral theory
4) Augmented Aerodynamic States Theory

Note that in this code, there is no aeroelastic system, the user can freely play with the parameters of the problem (reduced frequency, elastic axis position, phase lag between the plunge and the pitch motion, etc) to better understand how each of these parameters influence the aerodynamic loads.

Lastly, the power figure intends to show which configurations (depending on the parameters initially set) produces a worst case scenario for the aeroelastic stability of the system. Positive power means that the aerodynamic loads give energy to the system (unstable behavior), while negative powers substract energy from the system (stable behavior).

### 3. State-Space Flutter solver

This code is done to validate against K or PK flutter solvers, a time domain aeroelastic response of the system, but now represented and solved in a stae-space representation, which allows to check the exactly time evolution of the system to any given initial condition (can be freely set by the user).

### 4. Non Linear Aeroelastic solver with drag

This code is a modification of the previous, but now accounting for the effect of the unsteady drag in the aeroelastic instability, by introducing the equation of thrust to the system. Note that the code solves the equation with ode15i, as now the system is highly non linear, and classical state-space methods cannot lead with it.

### 5. Scripts

These Matlab scripts are intended to support each of the codes with a more academic and self-explanatory content, which can help the reader to better understand the mathematical basis of each methdology and theory used in the coding process. 

---
