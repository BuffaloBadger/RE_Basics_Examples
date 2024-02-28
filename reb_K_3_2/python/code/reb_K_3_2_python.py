# calculations for Reaction Engineering Basics Example K.3.2

# import libraries
import math
import numpy as np
import scipy as sp
import pandas as pd
from score_utils import solve_ivodes

# given and known constants
D = 5. # cm
dH = -14000. # cal /mol
V_flow = 1150. # cm^3 /min
C_p = 1.3 # cal /cm^-3 /K
R_R = 1.3
nDot_A_feed = 1.0 # mol /min
L = 50. # cm
T_feed = 300. # K
k_0 = 4.2E15 # cm^3 /mol /min
E = 18000 # cal /mol
Re = 1.987 # cal /mol /K

# derivatives function
def eval_derivs(cur_ind, cur_dep):
    # extract ind and dep vars for the current integration step
    nDot_A_cur = cur_dep[0]
    nDot_Z_cur = cur_dep[1]
    T_cur = cur_dep[2]

    # calculate the rate
    r = k_0*math.exp(-E/Re/T_cur)*nDot_A_cur*nDot_Z_cur/V_flow**2
    
    # evaluate the derivatives
    dnDotAdz = -math.pi*D**2/4*r
    dnDotZdz = math.pi*D**2/4*r
    dTdz = -math.pi*D**2/4*r*dH/V_flow/C_p

    # return the derivatives
    return [dnDotAdz, dnDotZdz, dTdz]

# residuals function
def eval_residuals(guess):
    # extract the individual guesses
    nAin_guess = guess[0]
    nZin_guess = guess[1]
    Tin_guess = guess[2]

    # define initial values
    ind_0 = 0.0
    dep_0 = np.array([nAin_guess, nZin_guess, Tin_guess])

    # define the stopping criterion
    f_var = 0
    f_val = L

    # solve the IVODEs
    z, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , eval_derivs)

    # check that the solution is valid
    if not(success):
        print(f"A solution of the IVODEs was NOT obtained: {message}")
    
    # extract the outlet values
    nAf = dep[0,-1]
    nZf = dep[1,-1]
    Tf = dep[2,-1]

    # evaluate the residuals
    residual1 = nDot_A_feed + R_R/(1 + R_R)*nAf - nAin_guess
    residual2 = R_R/(1 + R_R)*nZf - nZin_guess
    residual3 = Tin_guess - T_feed - R_R*(Tf - Tin_guess)

    # return the residual
    return [residual1, residual2, residual3]

# initial guess
init_guess = [nDot_A_feed, nDot_A_feed, T_feed + 5]

# solve the ATE
soln = sp.optimize.root(eval_residuals,init_guess)

# Check that the solution is converged
if not(soln.success):
    print(f"A solution of the ATEs was NOT obtained: {soln.message}")

# extract the results
n_dot_A_in = soln.x[0]
n_dot_Z_in = soln.x[1]
T_in = soln.x[2]

# initial values
ind_0 = 0.0
dep_0 = np.array([n_dot_A_in, n_dot_Z_in, T_in])

# define the stopping criterion
f_var = 0
f_val = L

# solve the IVODEs
z, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                    , eval_derivs)

# check that the solution is valid
if not(success):
    print(f"A solution of the IVODEs was NOT obtained: {message}")
    
# extract the outlet values
n_dot_A_out = dep[0,-1]
n_dot_Z_out = dep[1,-1]
T_out = dep[2,-1]

# display the results
print(' ')
print(f'Inlet molar flow of A: {n_dot_A_in:.3g} mol/min')
print(f'Inlet molar flow of Z: {n_dot_Z_in:.3g} mol/min')
print(f'Inlet temperature: {T_in:.3g} K')
print(f'Outlet molar flow of A: {n_dot_A_out:.3g} mol/min')
print(f'Outlet molar flow of Z: {n_dot_Z_out:.3g} mol/min')
print(f'Outlet temperature: {T_out:.3g} K')
