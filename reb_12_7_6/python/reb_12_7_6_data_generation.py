"""Data Generation for Example 12.7.6 of Reaction Engineering Basics"""

# import libraries
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt

# constants available to all functions
# given
V = 500 # cm^3
Vdot_in = 1.0 # cm^3 /s
CA_in = 0.015 # mol /cm^3
CB_in = 0.015 # mol /cm^3
Cp = 0.35*4.184 # J /g /K
rho = 0.93 # g /cm^3
dH = -20000 # J /mol
k0 = 3.24E12 # cm^3 /mol /s
E = 105000 # J /mol
# known
R = 8.314 # J /mol /K
# calculated
nA_in = Vdot_in*CA_in
nB_in = Vdot_in*CB_in

# make T available to all functiokns
T = -100.

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
    T_in = guess[4]

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
    residual_5 = -Vdot_in*rho*Cp*(T - T_in) - V*r*dH

    # return the residuals
    return np.array([residual_1, residual_2, residual_3, residual_4
                     , residual_5])

# perform the analysis
def perform_the_analysis():

    # enable this function to set T
    global T
    
    # set a range of outlet temperatures
    T_range = np.linspace(-25.0, 325.0, 351) + 273.15

	# set an initial guess
    init_guess = np.array([0.99*nA_in, 0.99*nB_in, 0.01*nA_in, 0.01*nA_in
            ,T_range[0] - 1.0])

    # allocate storage for the conversions and inlet temperatures
    f_range = np.zeros(351)
    Tin_range = np.zeros(351)

    # calculate the conversion and inlet temperature for each outlet temperature
    for i in range(0,351):

        # set the outlet temperature
        T = T_range[i]

        # solve the design equations
        solution = unknowns(init_guess)

        # save the conversion and inlet temperature
        f_range[i] = 100*(nA_in - solution[0])/nA_in
        Tin_range[i] = solution[4]

        # use the result as the initial guess for the next T
        init_guess = solution

    # generate , show and save graphs
    Tin_range = Tin_range - 273.15
    T_range = T_range -273.15

    # find the extinction and ignition points
    iLast = 0
    need_ignition = True
    for i in range(1,350):
        if need_ignition:
            if Tin_range[i] < Tin_range[iLast]:
                i_ignition = iLast
                need_ignition = False
            iLast = i
        else:
            if Tin_range[i] > Tin_range[iLast]:
                i_extinction = iLast
                break
            iLast = i
    
    # tabulate and save the inlet temperature, conversion, and outlet temperatures
    results_df = pd.DataFrame({'Inlet T':Tin_range, 'Conversion':f_range
                               ,'Outlet T':T_range})
    
    # save the results
    results_df.to_csv('reb_12_7_6/python/steady_state_data.csv',index=False)
	
    return

if __name__=="__main__":
    perform_the_analysis()
