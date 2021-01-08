# e-autocat
Core code for the paper "Fictitious Phase Separation in Li Layered Oxides Driven by Electro-Autocatalysis"

fp_solver is the solver for the Fokker-Planck equation.
Use PS_sim to generate the phase diagram of instability (Fig 5d)

BVSR, BVSRdeta are utility functions that compute the electrochemical reaction rate and its sensitivity with respect to the voltage under Butler-Volmer kinetics with a series resistance.
equilibrate_cavg is a utility function that returns the probability distribution at equilibrium given the average concentration.
fp_jacobian is a utility function that computes the jacobian of the ODE function in fp_solver.
