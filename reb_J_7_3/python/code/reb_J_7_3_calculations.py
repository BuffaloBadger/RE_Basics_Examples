""" Calculations for Reaction Engineering Basics Example J.7.3"""

# import libraries
import math
import numpy as np
import scipy as sp
import pandas as pd
from score_utils import solve_ivodes

# set the filespec for the results
Tin_results_file_spec = './reb_J_7_3/python/results/reb_J_7_3_Tin_result.csv'
profile_results_file_spec = \
    './reb_J_7_3/python/results/reb_J_7_3_profile_results.csv'

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

# derivatives function
def eval_derivs(cur_ind, cur_dep):
    # extract ind and dep vars for the current integration step
    nDot_A_cur = cur_dep[0]
    nDot_Z_cur = cur_dep[1]
    T_cur = cur_dep[2]

    # calculate the rate
    r = k_0*math.exp(-E/Re/T_cur)*(P*nDot_A_cur/Rw/T_cur
                                   /(nDot_A_cur + nDot_Z_cur))**2
    
    # evaluate the derivatives
    dnDotAdz = -2*math.pi*D**2/4*r
    dnDotZdz = math.pi*D**2/4*r
    dTdz = -math.pi*D**2/4*r*dH/(nDot_A_cur*Cp_A + nDot_Z_cur*Cp_Z)

    # return the derivatives
    return [dnDotAdz, dnDotZdz, dTdz]

# residuals function
def eval_residuals(guess):
    # define initial values
    ind_0 = 0.0
    dep_0 = np.array([nDot_A_in, 0.0, guess[0]])

    # define the stopping criterion
    f_var = 0
    f_val = L

    # solve the IVODEs
    z, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , eval_derivs)

    # check that the solution is valid
    if not(success):
        print(f"A solution of the IVODEs was NOT obtained: {message}")
    
    # extract Tf
    Tf = dep[2,-1]

    # evaluate the residual
    residual = Tf - T_out

    # return the residual
    return residual

# initial guess
init_guess = T_out - 100.

# solve the ATE
soln = sp.optimize.root(eval_residuals,init_guess)

# Check that the solution is converged
if not(soln.success):
    print(f"A solution of the ATEs was NOT obtained: {soln.message}")

# Extract the result
T_in_soln = soln.x[0]

# define initial values
ind_0 = 0.0
dep_0 = np.array([nDot_A_in, 0.0, T_in_soln])

# define the stopping criterion
f_var = 0
f_val = L

# solve the IVODEs
z, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                    , eval_derivs)

# check that the solution is valid
if not(success):
    print(f"A solution of the IVODEs was NOT obtained: {message}")

# extract the results
n_A = dep[0,:]
n_Z = dep[1,:]
T = dep[2,:]

# tabulate the results
Tin_results_df = pd.DataFrame(columns=['item','value','units'])
Tin_results_df.loc[0] = ['T_in' , T_in_soln, 'K']
profile_results_df = pd.DataFrame({'z':z , 'n_A':n_A, 'n_Z':n_Z, 'T':T})
    
# display the results
print(' ')
print(f'Inlet Temperature: {T_in_soln:.3g} K')
print(' ')
print('Molar Flow and Temperature Profiles')
print(profile_results_df.to_string(index=False))

# save the results
Tin_results_df.to_csv(Tin_results_file_spec, index=False)
profile_results_df.to_csv(profile_results_file_spec, index=False)
