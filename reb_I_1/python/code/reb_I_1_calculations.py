# solve 3 ATEs in Reaction Engineering Basics Example I.6

import math
import scipy as sp
import pandas as pd

# set filepath
filepath_to_results = './reb_I_1/python/results/'

# given and known constants
n_A_in = 500. # mol/h
v_fluid = 500. # L
n_Z_in = 0. # mol/h
cp_vol = 1170. # cal/L/K
v_flow = 250. # L/h
temp_in_K = 423. # K
dh_rxn = 18200. # cal/mol
k0 = 1.14E9 # L/mol/h
e_act = 16200. # cal/mol
gas_const_energy = 1.987 # cal/mol/K

# initial guess the solution
init_guess = [n_A_in/2.0, n_A_in/2.0, temp_in_K - 10.0]

# ATEs to solve, as residuals
def eval_residuals(guess):
    
    # Extract the guess values
    n_A_guess = guess[0]
    n_Z_guess = guess[1]
    temp_guess_K = guess[2]

    # Calculate the concentration of A, eqn 8
    conc_A = n_A_guess/v_flow;

    # Calculate the rate, eqn 7
    r = k0*math.exp(-e_act/gas_const_energy/temp_guess_K)*conc_A**2;
    
    # Evaluate the residuals, eqns 4-6
    residual_1 = n_A_in - n_A_guess - v_fluid*r
    residual_2 = n_Z_in - n_Z_guess + v_fluid*r
    residual_3 = cp_vol*v_flow*(temp_guess_K - temp_in_K) + v_fluid*r*dh_rxn
    return [residual_1, residual_2, residual_3]

# Solve the ATEs
soln = sp.optimize.root(eval_residuals,init_guess)

# Check that the solution is converged
if not(soln.success):
    print(f"A solution was NOT obtained: {soln.message}")

# Extract the results
n_A_out = soln.x[0]
n_Z_out = soln.x[1]
temp_out_K = soln.x[2]

# Display the results
print(' ')
print(f'Flow Rate of A: {n_A_out:.3g} mol/h')
print(f'Flow Rate of Z: {n_Z_out:.3g} mol/h')
print(f'Temperature: {temp_out_K - 273.15:.3g} °C')
print(' ')

# save the results to a .csv file
data = [['Flow Rate of A', f'{n_A_out}', ' mol/h'],
    ['Flow Rate of Z', f'{n_Z_out}', ' mol/h'],
    ['Temperature',f'{temp_out_K - 273.15}', ' °C']]
result = pd.DataFrame(data, columns=['item','value','units'])
result.to_csv(filepath_to_results + 'reb_I_1_results.csv', index=False)
