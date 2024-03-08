"""Calculations for Example J.7.2 of Reaction Engineering Basics"""

# import libraries
import math
import numpy as np
import scipy as sp
import pandas as pd
from score_utils import solve_ivodes

# given and known constants
T_out = 400. # K
D = 1. # in
k_0 = 7.49E9*61.02 # in^3 /mol /min
E = 15300. # cal /mol
P = 4. # atm
dH = -14500. # cal /mol
Cp_A = 10.9 # cal /mol /K
Cp_Z = 21.8 # cal /mol /K
L = 100. # in
nDot_A_in = 1.5 # mol /min
Re = 1.987 # cal /mol /K
Rw = 0.08206*61.02 # in^3 atm /mol /K

# reactor design equations as derivative expressions
def derivatives(z, dep):
    # extract the dependent variables for the current integration step
    nDot_A = dep[0]
    nDot_Z = dep[1]
    T = dep[2]
    
    # calculate other IVODE variables
    r = other_ivode_variables(T, nDot_A, nDot_Z)

    # evaluate the derivatives
    dnDotAdz = -2*math.pi*D**2/4*r
    dnDotZdz = math.pi*D**2/4*r
    dTdz = -math.pi*D**2/4*r*dH/(nDot_A*Cp_A + nDot_Z*Cp_Z)

    # return the derivatives
    return [dnDotAdz, dnDotZdz, dTdz]

# calculate other IVODE variables
def other_ivode_variables(T, nDot_A, nDot_Z):
    r = k_0*math.exp(-E/Re/T)*(P*nDot_A/Rw/T/(nDot_A + nDot_Z))**2
    return r

# calculate IVODE initial and final values
def initial_and_final_values(T_in_guess):
    # set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nDot_A_in, 0.0, T_in_guess])

    # define the stopping criterion
    f_var = 0
    f_val = L   

    # return the initial and final values     
    return ind_0, dep_0, f_var, f_val

# implicit equation for IVODE initial value as residual
def residual(guess):
    # solve the reactor design equations
    T_in_guess = guess[0]
    z, nDot_A, nDot_Z, T = profiles(T_in_guess)

    # extract Tf
    Tf = T[-1]

    # evaluate the residual
    residual = Tf - T_out

    # return the residual
    return residual

# solve the reactor design equations
def profiles(T_in):
    # get the initial and final values
    ind_0, dep_0, f_var, f_val = initial_and_final_values(T_in)

    # solve the IVODEs
    z, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , derivatives)

    # check that a solution was found
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")
    
    # extract dependent variable profiles
    nDot_A = dep[0,:]
    nDot_Z = dep[1,:]
    T = dep[2,:]

    # return all profiles
    return z, nDot_A, nDot_Z, T

# complete the assignment
def complete_the_assignment():
    # initial guess
    initial_guess = T_out - 100.

    # solve the implict equation for Tin
    soln = sp.optimize.root(residual,initial_guess)

    # check that the solution is converged
    if not(soln.success):
        print(f"The initial temperature was NOT found: {soln.message}")

    # extract the result
    T_in = soln.x[0]

    # solve the reactor design equations
    z, nDot_A, nDot_Z, T = profiles(T_in)

    # tabulate the results
    Tin_results_df = pd.DataFrame(columns=['item','value','units'])
    Tin_results_df.loc[0] = ['T_in' , T_in, 'K']
    profile_results_df = pd.DataFrame({'z':z , 'nDot_A':nDot_A
                                       , 'nDot_Z':nDot_Z, 'T':T})
        
    # display the results
    print(' ')
    print(f'Inlet Temperature: {T_in:.3g} K')
    print(' ')
    print('Molar Flow and Temperature Profiles')
    print(profile_results_df.to_string(index=False))

    # save the results
    Tin_results_file_spec = './reb_J_7_2/python/results/reb_J_7_2_Tin_results.csv'
    profile_results_file_spec = \
            './reb_J_7_2/python/results/reb_J_7_2_profile_results.csv'
    Tin_results_df.to_csv(Tin_results_file_spec, index=False)
    profile_results_df.to_csv(profile_results_file_spec, index=False)

if __name__=="__main__":
    complete_the_assignment()