"""Calculations for Reaction Engineering Basics Example J.7.2"""

# import libraries
import math
import numpy as np
import scipy as sp
import pandas as pd
from score_utils import solve_ivodes

# set the filespec for the results
results_file_spec = './reb_J_7_2/python/results/reb_J_7_2_results.csv'

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

# make Vdot available to the derivatives function
Vdot = float('nan')

# derivatives function
def eval_derivs(cur_ind, cur_dep):
    # extract ind and dep vars for the current integration step
    nDot_A_cur = cur_dep[0]
    nDot_Z_cur = cur_dep[1]
    T_cur = cur_dep[2]

    # calculate the rate
    r = k_0*math.exp(-E/Re/T_cur)*nDot_A_cur**2/Vdot**2
    
    # evaluate the derivatives
    dnDotAdz = -math.pi*D**2/4*r
    dnDotZdz = math.pi*D**2/4*r
    dTdz = -math.pi*D**2/4*r*dH/Vdot/Cp

    # return the derivatives
    return [dnDotAdz, dnDotZdz, dTdz]

# residuals function
def eval_residuals(guess):
    # make the guess available to the derivatives function
    global Vdot
    Vdot = guess[0]

    # define initial values
    ind_0 = 0.0
    dep_0 = np.array([nDot_A_in, nDot_Z_in, T_in])

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
    Tf_calc = dep[2,-1]

    # evaluate the residual
    residual = T_f - Tf_calc

    # return the residual
    return residual

# initial guess
init_guess = 100.

# solve the ATE
soln = sp.optimize.root(eval_residuals,init_guess)

# Check that the solution is converged
if not(soln.success):
    print(f"A solution of the ATEs was NOT obtained: {soln.message}")

# define initial values
ind_0 = 0.0
dep_0 = np.array([nDot_A_in, nDot_Z_in, T_in])

# define the stopping criterion
f_var = 0
f_val = L

# solve the IVODEs
z, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                    , eval_derivs)

# check that the solution is valid
if not(success):
    print(f"A solution of the IVODEs was NOT obtained: {message}")

# extract final values
nDot_A_f = dep[0,-1]
nDot_Z_f = dep[1,-1]
    
# display the results
print(' ')
print(f'Volumetric Flow Rate: {soln.x[0]:.3g} cm^3^ min^-1^')
print(f'Final Molar Flow of A: {nDot_A_f:.3g} mol min^-1^')
print(f'Final Molar Flow of Z: {nDot_Z_f:.3g} mol min^-1^')

# Save the results
data = [['$\dot{V}$', f'{Vdot}', 'cm^3^ min^-1^'],
    ['$\dot{n}_{A,f}$', f'{nDot_A_f}', 'mol min^-1^'],
    ['$\dot{n}_{Z,f}$', f'{nDot_Z_f}', 'mol min^-1^']]
results_df = pd.DataFrame(data, columns=['item','value','units'])
results_df.to_csv(results_file_spec, index=False)
