"""Calculations for Example 9.6.1 of Reaction Engineering Basics"""

# import libraries
import math
import numpy as np
import pandas as pd
from score_utils import solve_ivodes
import matplotlib.pyplot as plt

# given and known constants
dH = -101.2E3 # J/mol
k_0 = 5.11e4 * 3600 # L/mol/h
E = 74.8e3 # J/mol
T_0 = 180 + 273.15 # K
V = 1900.0 # L
CA_0 = 2.9 # mol/L
CB_0 = 3.2 # mol/L
Cp = 1.23 * 4.184 # J/g/K
rho = 1.02 * 1000.0 # g/L
t_f = 2.0 # h
Re = 8.314 # J/mol/K

# derivatives function
def derivatives(ind, dep):
	# extract the dependent variables for this integration step
    nA = dep[0]
    nB = dep[1]
    T = dep[4]

    # calculate the rate
    CA = nA/V
    CB = nB/V
    k = k_0*math.exp(-E/Re/T)
    r = k*CA*CB

    # evaluate the derivatives
    dnAdt = -V*r
    dnBdt = -V*r
    dnYdt = V*r
    dnZdt = V*r
    dTdt = -r*dH/rho/Cp

    # return the derivatives
    return [dnAdt, dnBdt, dnYdt, dnZdt, dTdt]

# reactor model
def profiles():
    # set the initial values
    nA_0 = CA_0*V
    nB_0 = CB_0*V
    ind_0 = 0.0
    dep_0 = np.array([nA_0, nB_0, 0.0, 0.0, T_0])

    # define the stopping criterion
    f_var = 0
    f_val = t_f
     
    # solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , derivatives, True)

    # check that a solution was found
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    nA = dep[0,:]
    nB = dep[1,:]
    nY = dep[2,:]
    nZ = dep[3,:]
    T = dep[4,:]

    # return all profiles
    return t, nA, nB, nY, nZ, T

# perform the analysis
def perform_the_analysis():
    # solve the reactor design equations
    t, nA, nB, nY, nZ, T = profiles()

    # calculate the other quantities of interest
    CA = nA/V
    CB = nB/V
    CY = nY/V
    CZ = nZ/V
    k = k_0*np.exp(-E/Re/T)
    r = k*CA*CB
    T_C = T - 273.15
    
    # tabulate the results
    results_df = pd.DataFrame({'t':t , 'nA':nA, 'nB':nB
                   , 'nY':nY, 'nZ':nZ, 'T':T_C, 'r':r})

    # save the results
    results_df.to_csv("reb_9_6_1/python/results/reb_9_6_1_results.csv"
                      , index=False)
    
    # display and save the graphs
    plt.figure() 
    plt.plot(t,CA,label='A')
    plt.plot(t,CB,label='B')
    plt.plot(t,CY,linestyle=':',linewidth=4, label='Y')
    plt.plot(t,CZ,label='Z')
    plt.xlabel("$Time \; (h)$")
    plt.ylabel("$Concentration \; (mol \; L^{-1})$")
    plt.legend()
    plt.savefig('reb_9_6_1/python/results/reb_9_6_1_concentrations.png')
    plt.show()

    plt.figure()
    plt.plot(t,T_C)
    plt.xlabel("$Time \; (h)$")
    plt.ylabel("$Temperature \; (°C)$")
    plt.savefig('reb_9_6_1/python/results/reb_9_6_1_temperature.png')
    plt.show()

    plt.figure(3)
    plt.plot(t,r)
    plt.xlabel("$Time \; (h)$")
    plt.ylabel("$Rate \; (mol \; L^{-1} \; h^{-1})$")
    plt.savefig('reb_9_6_1/python/results/reb_9_6_1_rate.png')
    plt.show()

    return

if __name__=="__main__":
    perform_the_analysis()