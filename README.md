# SYMPLECTIC-SIMULATIONS
This repository contains the codes and results from my bachelor's thesis on variational integrators.

In the main folder you can find the Hamiltonian formulation for simulating: the oscillator, Henon--Heiles, the simple pendulum, the 2-body problem (2-dimensional), the general n-body problem (2 and 3 dimensional), and a restricted case where a particle is influenced by n-1 fixed masses.

In the "methods" folder, codes are given for numerical methods: Euler, symplectic Euler (1,2), Störmer--Verlet (12,21), Heun, Ralston, Midpoint, RK-3, RK-4, RK-5 (all explicit, symplectic methods MUST be used for separable problems).

In the "plotting" folder there are functions for drawing the full trajectories of the n-body problems (for both the 2-dimensional and 3-dimensional case) and the corresponding videos that represent the full trajectories, for plotting the outer solar system case and the restricted case with n-1 fixed masses. Also, the script "plot_config" can be used for adding a better configuration (text and line width, letter style, etc.) for the plots.

In the "examples" folder one can find a functions that returns 3 sets of initial conditions for the n-body problem, namely a chatoic case, a periodic case and the outer solar system case.

All images and videos produced by the codes are to be saved in the folder "results", the user will find some examples of the simulations already there.
