""" Calculations for Reaction Engineering Basics Example J.7"""

# Import libraries
import math
from score_utils import solve_ivodes
import numpy as np
import pandas as pd
import scipy as sp

# Set filepath
filepath_to_results = './reb_J_7_1/python/results/'

# Given and known constants
V = 2.0 # m^3
dH_1 = -3470 # cal/mol
Cp_A = 7.5 # cal/mol/K
Cp_B = 8.5 # cal/mol/K
Cp_Y = 12.1 # cal/mol/K
Cp_Z = 5.7 # cal/mol/K
k_0_1 = 83 # m^3 /mol /h
E_1 = 10200 # cal/mol
# Known
Re = 1.987 # cal/mol
Rw = 8.206E-5 # m^3 atm/mol/K

# Derivatives function
def eval_derivs(cur_indep, cur_dep):
    
    # Extract the dependent variables for the current integration step
    n_A_cur = cur_dep[0]
    n_B_cur = cur_dep[1]
    n_Y_cur = cur_dep[2]
    n_Z_cur = cur_dep[3]
    T_cur = cur_dep[4]

    # Create mass matrix, setting all elements to zero
    mass_matrix = np.zeros((6,6))

    # Add 1 on the diagonal for the first 4 rows
    mass_matrix[0,0] = 1.0
    mass_matrix[1,1] = 1.0
    mass_matrix[2,2] = 1.0
    mass_matrix[3,3] = 1.0

    # Add the elements for the energy balance
    mass_matrix[4,4] = n_A_cur*Cp_A + n_B_cur*Cp_B \
        + n_Y_cur*Cp_Y + n_Z_cur*Cp_Z
    mass_matrix[4,5] = -V*Re/Rw

    # Add the elements for the ideal gas law equation
    mass_matrix[5,0] = Rw*T_cur
    mass_matrix[5,1] = Rw*T_cur
    mass_matrix[5,2] = Rw*T_cur
    mass_matrix[5,3] = Rw*T_cur
    mass_matrix[5,4] = Rw*(n_A_cur + n_B_cur + n_Y_cur + n_Z_cur)
    mass_matrix[5,5] = -V

    # Calculate the rate
    C_A = n_A_cur/V
    C_B = n_B_cur/V
    r = k_0_1*math.exp(-E_1/Re/T_cur)*C_A*C_B

    # Create right side vector
    rhs = np.array([-V*r, -V*r, V*r, V*r, -V*r*dH_1, 0])

    # Evaluate the derivatives
    derivs = sp.linalg.solve(mass_matrix, rhs)

    # Return the derivatives
    return derivs

# Initial values
ind_0 = 0.0
dep_0 = np.array([190.0, 190.0, 0.0, 0.0, 450.0, 7.0])

# Stopping criterion
f_var = 0
f_val = 2.0

# Solve the IVODEs
t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , eval_derivs)

# Check that the solution was found
if not(success):
    print(f"A solution was NOT obtained: {message}")

# Extract the results
n_A = dep[0,:]
n_B = dep[1,:]
n_Y = dep[2,:]
n_Z = dep[3,:]
T = dep[4,:]
P = dep[5,:]

# Tabulate the results
results_df = pd.DataFrame({'t':t , 'n_A':n_A, 'n_B':n_B, 'n_Y':n_Y
                   ,'n_Z':n_Z, 'T':T, 'P':P})

# Display the results
print(' ')
print(results_df.to_string(index=False))
print(' ')

# Save the results
results_df.to_csv(filepath_to_results + 'reb_J_7_1_results.csv', index=False)