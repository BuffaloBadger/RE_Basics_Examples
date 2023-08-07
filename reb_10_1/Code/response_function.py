import pandas as pd
import numpy as np
import math
import rebutils as reb
import matplotlib.pyplot as plt

# Define the response function
def response_function(inputs, par):
    # Given and known constants
    P = 1.0 # atm
    D = 1.0 # cm
    L = 10.0 # cm
    T = 1500. # K
    mol_per_cc_at_stp = 1.0E-3/22.4

    # extract the rate coefficient
    k_test = 10.0**par

    # allocate storage for the responses
    shape = np.shape(inputs)
    n_data = shape[0]
    response = np.zeros(n_data)

    # loop through the data points
    for i, input in enumerate(inputs):
        # define the reactor model
        def reactor_model(z,nDot):
            nDot_A = nDot[0]
            nDot_Y = nDot[1]
            nDot_Z = nDot[2]
            PA = nDot_A/(nDot_A + nDot_Y + nDot_Z)*P
            r = k_test*PA
            ddz = np.array([-r, r, r]) * math.pi *D**2/4
            return ddz

        # Solve the reactor model equations
        z0 = 0.0
        n0 = np.array([input[0], input[1], input[2]])*mol_per_cc_at_stp
        f_var = 0
        f_val = L
        soln = reb.solveIVODEs(z0, n0, f_var, f_val, reactor_model)
        if not(soln.success):
            print("A solution was NOT obtained:")
            print(soln.message)
        # Calculate the response
        nA = soln.y[0,-1]
        nAin = input[0]*mol_per_cc_at_stp
        response[i] = 100*(nAin - nA)/nAin
    return response

