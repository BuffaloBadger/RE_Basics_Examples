import pandas as pd
import numpy as np
import rebutils as reb
import matplotlib.pyplot as plt

# Define the response function
def response_function(inputs, par1, par2, par3, par4, par5):
    # Given and known constants
    P = 1.0 # atm
    m = 3.0 # g
    VFR_in = 0.85 # L/min
    T = 400. + 273.15 # K
    K = 12.2

    # extract the rate coefficient
    k_test = 10.0**par1
    alpha_A = par2
    alpha_B = par3
    alpha_Y = par4
    alpha_Z = par5

    # allocate storage for the responses
    shape = np.shape(inputs)
    n_data = shape[0]
    response = np.zeros(n_data)

    # loop through the data points
    for i, input in enumerate(inputs):
        # define the reactor model
        def reactor_model(z,n):
            nA = n[0]
            nB = n[1]
            nY = n[2]
            nZ = n[3]
            ntot = nA + nB + nY + nZ
            PA = nA/ntot*P
            PB = nB/ntot*P
            PY = nY/ntot*P
            PZ = nZ/ntot*P
            r = k_test*PA**alpha_A*PB**alpha_B*PY**alpha_Y*PZ**alpha_Z*\
                (1 - PY*PZ/K/PA/PB)
            ddm = np.array([-r, -r, r, r])
            return ddm

        # Solve the reactor model equations
        m0 = 0.0
        n0 = np.array([input[0], input[1], input[2], input[3]])*VFR_in/22.4
        f_var = 0
        f_val = m
        soln = reb.solveIVODEs(m0, n0, f_var, f_val, reactor_model)
        if not(soln.success):
            print("A solution was NOT obtained:")
            print(soln.message)
        # Calculate the response
        n_A = soln.y[0,-1]
        n_B = soln.y[1,-1]
        n_Y = soln.y[2,-1]
        n_Z = soln.y[3,-1]
        n_tot = n_A + n_B + n_Y + n_Z
        response[i] = n_A/n_tot*P
    return response
