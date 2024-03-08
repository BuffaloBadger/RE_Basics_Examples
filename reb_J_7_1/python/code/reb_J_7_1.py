"""Calculations for Example J.7.1 of Reaction Engineering Basics"""

# import libraries
import math
from score_utils import solve_ivodes
import numpy as np
import pandas as pd
import scipy as sp

# given and known constants
V = 2.0 # m^3
dH_1 = -3470 # cal/mol
Cp_A = 7.5 # cal/mol/K
Cp_B = 8.5 # cal/mol/K
Cp_Y = 12.1 # cal/mol/K
Cp_Z = 5.7 # cal/mol/K
k_0_1 = 83 # m^3 /mol /h
E_1 = 10200 # cal/mol
Re = 1.987 # cal/mol
Rw = 8.206E-5 # m^3 atm/mol/K

# reactor design equations as derivative expressions
def derivatives(t, dep):
    # Extract the dependent variables for this integration step
    n_A = dep[0]
    n_B = dep[1]
    n_Y = dep[2]
    n_Z = dep[3]
    T = dep[4]

    # Create mass matrix, setting all elements to zero
    mass_matrix = np.zeros((6,6))

    # Add 1 on the diagonal for the first 4 rows
    mass_matrix[0,0] = 1.0
    mass_matrix[1,1] = 1.0
    mass_matrix[2,2] = 1.0
    mass_matrix[3,3] = 1.0

    # Add the elements for the energy balance
    mass_matrix[4,4] = n_A*Cp_A + n_B*Cp_B \
        + n_Y*Cp_Y + n_Z*Cp_Z
    mass_matrix[4,5] = -V*Re/Rw

    # Add the elements for the ideal gas law equation
    mass_matrix[5,0] = Rw*T
    mass_matrix[5,1] = Rw*T
    mass_matrix[5,2] = Rw*T
    mass_matrix[5,3] = Rw*T
    mass_matrix[5,4] = Rw*(n_A + n_B + n_Y + n_Z)
    mass_matrix[5,5] = -V

    # Calculate the rate
    r = other_ivode_variables(T, n_A, n_B)

    # Create right side vector
    rhs = np.array([-V*r, -V*r, V*r, V*r, -V*r*dH_1, 0])

    # Evaluate the derivatives
    derivs = sp.linalg.solve(mass_matrix, rhs)

    # Return the derivatives
    return derivs

# calculate other IVODE variables
def other_ivode_variables(T, n_A, n_B):
    C_A = n_A/V
    C_B = n_B/V
    r = k_0_1*math.exp(-E_1/Re/T)*C_A*C_B
    return r

# calculate IVODE initial and final values
def initial_and_final_values():
    ind_0 = 0.0
    dep_0 = np.array([190.0, 190.0, 0.0, 0.0, 450.0, 7.0])
    f_var = 0
    f_val = 2.0
    return ind_0, dep_0, f_var, f_val

# solve the reactor design equations
def profiles():
    # get the initial values and stopping criterion
    ind_0, dep_0, f_var, f_val = initial_and_final_values()

    # solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0
                            , f_var, f_val, derivatives)
    
    # check that a solution was found
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")

    # extract dependent variable profiles
    n_A = dep[0,:]
    n_B = dep[1,:]
    n_Y = dep[2,:]
    n_Z = dep[3,:]
    T = dep[4,:]
    P = dep[5,:]

    # return all profiles
    return t, n_A, n_B, n_Y, n_Z, T, P

# complete the assignment
def complete_the_assignment():
    # get the solution of the reactor design equations
    t, n_A, n_B, n_Y, n_Z, T, P = profiles()

    # tabulate the results
    results_df = pd.DataFrame({'t':t , 'n_A':n_A, 'n_B':n_B
                   , 'n_Y':n_Y, 'n_Z':n_Z, 'T':T, 'P':P})
    
    # display the results
    print(' ')
    print(results_df.to_string(index=False))
    print(' ')

    # save the results
    file_spec = './reb_J_7_1/python/results/reb_J_7_1_results.csv'
    results_df.to_csv(file_spec, index=False)

if __name__=="__main__":
    complete_the_assignment()
