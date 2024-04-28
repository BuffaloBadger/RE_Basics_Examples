""" Calculations for Reaction Engineering Basics Example I.6"""

# Import libraries
import math
import scipy as sp

# Given and known constants
n_A_in = 500. # mol/h
V_fluid = 500. # L
n_Z_in = 0. # mol/h
Cp_vol = 1170. # cal/L/K
V_flow = 250. # L/h
T_in_K = 423. # K
dH_rxn = 18200. # cal/mol
k_0 = 1.14E9 # L/mol/h
E = 16200. # cal/mol
Re = 1.987 # cal/mol/K

# Residuals function
def eval_residuals(guess):
    
    # Extract the guess values
    n_A_guess = guess[0]
    n_Z_guess = guess[1]
    temp_guess_K = guess[2]

    # Calculate the concentration of A
    C_A = n_A_guess/V_flow

    # Calculate the rate
    r = k_0*math.exp(-E/Re/temp_guess_K)*C_A**2
    
    # Evaluate the residuals
    residual_1 = n_A_in - n_A_guess - V_fluid*r
    residual_2 = n_Z_in - n_Z_guess + V_fluid*r
    residual_3 = Cp_vol*V_flow*(temp_guess_K - T_in_K) + V_fluid*r*dH_rxn

    # Return the residuals
    return [residual_1, residual_2, residual_3]

# Initial guess
init_guess = [n_A_in/2.0, n_A_in/2.0, T_in_K - 10.0]

# Solve the ATEs
soln = sp.optimize.root(eval_residuals,init_guess)

# Check that the solution is converged
if not(soln.success):
    print(f"A solution was NOT obtained: {soln.message}")

# Extract the solution
n_A_out = soln.x[0]
n_Z_out = soln.x[1]
temp_out_K = soln.x[2]

# Display the solution
print(' ')
print(f'Flow Rate of A: {n_A_out:.3g} mol/h')
print(f'Flow Rate of Z: {n_Z_out:.3g} mol/h')
print(f'Temperature: {temp_out_K - 273.15:.3g} °C')
print(' ')
