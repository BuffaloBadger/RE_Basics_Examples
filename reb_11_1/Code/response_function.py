import numpy as np
import pandas as pd
import scipy as sp

# define the response function
def response_function(inputs, log_k_guess):
    # Given and known constants
    V = 0.1 # L
    T = 35 + 273.15 # K

    # extract the guess for k
    k_guess = 10**log_k_guess

    # allocate storage for the responses
    shape = np.shape(inputs)
    n_data = shape[0]
    response = np.zeros(n_data)

    # read the experimental responses to use as guesses when solving the reactor
    # design equations  
    df = pd.read_csv('./reb_11_1/Data/reb_11_1_data.csv')
    CYout = df['CYout'].to_numpy()

    # loop through all of the experiments
    for i in range(0,n_data):
        # extract the inputs
        Vdot = inputs[i,0]
        CA0 = inputs[i,1]
        CB0 = inputs[i,2]
        CY0 = inputs[i,3]
        CZ0 = inputs[i,4]
        # calculate inlet molar flow rates
        nAin = CA0*Vdot
        nBin = CB0*Vdot
        nYin = CY0*Vdot
        nZin = CZ0*Vdot

        # define the reactor model residuals
        def reactor_model(x):
            # extract the test values for the unknowns
            nA = x[0]
            nB = x[1]
            nY = x[2]
            nZ = x[3]

            # calculate the rate
            CA = nA/Vdot
            CB = nB/Vdot
            r = k_guess*CA*CB

            # evaluate the reactor model residuals
            r1 = nAin - nA - r*V
            r2 = nBin - nB - r*V
            r3 = nYin - nY + r*V
            r4 = nZin - nZ + r*V
            return [r1, r2, r3, r4]
        
        # guess the solution to the reactor model
        nY_guess = CYout[i]*Vdot
        extent = nY_guess - nYin
        nA_guess = nAin - extent
        nB_guess = nBin - extent
        nZ_guess = nZin + extent
        guess = [nA_guess, nB_guess, nY_guess, nZ_guess]

        # solve the reactor model equations
        soln = sp.optimize.root(reactor_model,guess)
        if not(soln.success):
            print("A solution was NOT obtained:")
            print(soln.message)
        
        # calculate the response, CYout
        response[i] = soln.x[2]/Vdot
    return response