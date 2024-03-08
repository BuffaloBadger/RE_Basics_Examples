"""Calculations for Example J.7.2 of Reaction Engineering Basics"""

# import libraries
import math
import numpy as np
import scipy as sp
import pandas as pd
from score_utils import solve_ivodes

# given and known constants
T_f = 325. # K
D = 5. # cm
dH = -14000. # cal /mol
Cp = 1.3 # cal /cm^3 /K
nDot_A_in = 1.0 # mol /min
nDot_Z_in = 0.0 # mol /min
T_in = 300. # K
L = 50.0 # cm
k_0 = 4.2E15 # cm^3 /mol /min
E = 18000. # cal /mol
Re = 1.987 # cal /mol /K

# make IVODE constant available to all functions
Vdot = float('nan')

# reactor design equations as derivative expressions
def derivatives(z,dep):
    # extract the dependent variables for the current integration step
    nDot_A = dep[0]
    nDot_Z = dep[1]
    T = dep[2]

    # calculate the rate
    r = other_ivode_variables(T, nDot_A)
    
    # evaluate the derivatives
    dnDotAdz = -math.pi*D**2/4*r
    dnDotZdz = math.pi*D**2/4*r
    dTdz = -math.pi*D**2/4*r*dH/Vdot/Cp

    # return the derivatives
    return [dnDotAdz, dnDotZdz, dTdz]

# calculate other IVODE variables
def other_ivode_variables(T, nDot_A):
    r = k_0*math.exp(-E/Re/T)*nDot_A**2/Vdot**2
    return r

# calculate unknown IVODE constant
def ivode_constant():
    # allow this function to set Vdot
    global Vdot

    # initial guess
    initial_guess = 100.0

    # solve the implicit equation for Vdot
    soln = sp.optimize.root(residual,initial_guess)

    # check that the solution converged
    if not(soln.success):
        print(f"Vdot was NOT found: {soln.message}")

    # set Vdot and return
    Vdot = soln.x[0]
    return Vdot

# implicit equation for the IVODE constant as residual
def residual(guess):
    # allow this function to set Vdot
    global Vdot
    Vdot = guess[0]
  
    # solve the reactor design equations
    z, n_A, n_Z, T = profiles()

    # extract the calculated final T
    T_f_calc = T[-1]

    # evaluate the residual
    residual_1 = T_f - T_f_calc

    # return the residual
    return residual_1

# calculate IVODE initial and final values
def initial_and_final_values():
    ind_0 = 0.0
    dep_0 = np.array([nDot_A_in, nDot_Z_in, T_in])
    f_var = 0
    f_val = L
    return ind_0, dep_0, f_var, f_val

# solve the reactor design equations
def profiles():
    # get the initial values and stopping criterion
    ind_0, dep_0, f_var, f_val = initial_and_final_values()
 
    # solve the IVODEs
    z, dep, success, message = solve_ivodes(ind_0, dep_0
                            , f_var, f_val, derivatives)

    # check that a solution was found
    if not(success):
        print(f"The profiles were NOT found: {message}")

    # extract dependent variable profiles
    nDot_A = dep[0,:]
    nDot_Z = dep[1,:]
    T = dep[2,:]

    # return all profiles
    return z, nDot_A, nDot_Z, T

# calculate other quantities of interest
def other_quantities_of_interest(nDot_A, nDot_Z):
    # extract the final values
    nDot_A_f = nDot_A[-1]
    nDot_Z_f = nDot_Z[-1]

    # return the other quantities of interest
    return nDot_A_f, nDot_Z_f

# complete the assignment
def complete_the_assignment():
    # calculate Vdot
    Vdot = ivode_constant()

    # get the solution of the reactor design equations
    z, nDot_A, nDot_Z, T = profiles()

    # get the quantities of interest
    nDot_A_f, nDot_Z_f = other_quantities_of_interest(nDot_A, nDot_Z)

    # display the results
    print(' ')
    print(f'Volumetric Flow Rate: {Vdot:.3g} cm^3 /min')
    print(f'Final Molar Flow of A: {nDot_A_f:.3g} mol /min')
    print(f'Final Molar Flow of Z: {nDot_Z_f:.3g} mol /min')
    print(' ')

    # save the results
    data = [['$\dot{V}$', f'{Vdot}', 'cm^3^ min^-1^'],
    ['$\dot{n}_{A,f}$', f'{nDot_A_f}', 'mol min^-1^'],
    ['$\dot{n}_{Z,f}$', f'{nDot_Z_f}', 'mol min^-1^']]
    results_df = pd.DataFrame(data, columns=['item','value','units'])
    file_spec = './reb_J_7_2/python/results/reb_J_7_2_results.csv'
    results_df.to_csv(file_spec, index=False)

if __name__=="__main__":
    complete_the_assignment()