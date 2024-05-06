"""Calculations for Example 12.7.6 of Reaction Engineering Basics"""

# import libraries
import numpy as np
from score_utils import solve_ivodes
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt

# constants available to all functions
# given
V = 500.0 # cm^3
Vdot_in = 1.0 # cm^3 /s
CA_in = 0.015 # mol /cm^3
CB_in = 0.015 # mol /cm^3
T_in = 50 + 273.15 # K
Cp = 0.35*4.184 # J /g /K
rho = 0.93 # g /cm^3
dH = -20000.0 # J /mol
k0 = 3.24E12 # cm^3 /mol /s
E = 105000.0 # J /mol
# known
R = 8.314 # J /mol /K
# calculated
nA_in = CA_in*Vdot_in
nB_in = CB_in*Vdot_in

# make T available to all functions
T = -1.0

# reactor model
def unknowns(init_guess):
     
	# solve the ATEs
    soln = sp.optimize.root(residuals,init_guess)

    # check that the solution is converged
    if not(soln.success):
        print(f"The initial temperature was NOT found: {soln.message}")

    # return the solution
    return soln.x

# residuals function
def residuals(guess):
    # extract the indiviaual guesses
    nA = guess[0]
    nB = guess[1]
    nY = guess[2]
    nZ = guess[3]

    # calculate the rate
    k = k0*np.exp(-E/R/T)
    CA = nA/Vdot_in
    CB = nB/Vdot_in
    r = k*CA*CB

    # evaluate the residuals
    residual_1 = nA_in - nA - V*r
    residual_2 = nB_in - nB - V*r
    residual_3 = -nY + V*r
    residual_4 = -nZ + V*r

    # return the residuals
    return np.array([residual_1, residual_2, residual_3, residual_4])

# perform the analysis
def perform_the_analysis():
    # allow this function to set T
    global T

	# set an initial guess for the low-conversion steady state
    init_guess = np.array([0.99*nA_in, 0.99*nB_in, 0.01*nA_in, 0.01*nA_in])

    # set a range of T values
    T_range = np.linspace(0., 300.,300) + 273.15

    # allocate storage for Qabs and Qgen
    Qabs = np.zeros(300)
    Qgen = np.zeros(300)

    # calculate Qabs and Qgen for each T
    for iT in range(0,300):
        # set the temperature
        T = T_range[iT]

        # solve the reactor design equations
        solution = unknowns(init_guess)

        # save Qabs and Qgen
        Qabs[iT] = Vdot_in*rho*Cp*(T_range[iT]-T_in)
        Qgen[iT] = -V*k0*np.exp(-E/R/T)*solution[0]/Vdot_in*solution[1]/Vdot_in\
            *dH
        
        # set the guess for the next T
        init_guess = solution

	# create, show and save the graph
    plt.figure(1) 
    plt.plot(T_range - 273.15, Qgen, color = 'r', label = 'heat generated')
    plt.plot(T_range - 273.15, Qabs, color = 'b', label = 'heat absorbed')
    plt.text(52, -20 ,'low')
    plt.text(138, 100 ,'middle')
    plt.text(265, 275 ,'high')
    plt.xlabel("Outlet Temperature (C)")
    plt.ylabel("Heat Generation or Absorption Rate (J s$^{-1}$)")
    plt.legend()
    plt.savefig('reb_12_7_6/python/heat_vs_T.png')
    plt.show()
    
    return

if __name__=="__main__":
    perform_the_analysis()
