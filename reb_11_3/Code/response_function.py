import pandas as pd
import numpy as np
import scipy as sp
import rebutils as reb
import matplotlib.pyplot as plt

# define the response function
def response_function(inputs, log_kf_guess, log_kr_guess):

    # Given and known constants
    V = 3.0 # gal

    # Read the experimental responses to use as guesses when solving the reactor
    # design equations
    df = pd.read_csv('./reb_10_3/Data/reb_10_3_data.csv')
    CAout = df['CA'].to_numpy()

    # allocate storage for the responses
    shape = np.shape(inputs)
    n_data = shape[0]
    response = np.zeros(n_data)

    # extract the guesses for kf and kr
    kf_guess = 10.0**log_kf_guess
    kr_guess = 10.0**log_kr_guess

    # loop through all of the experiments
    for i in range(0,n_data):
        # extract the inputs
        Vdot = inputs[i,0]
        CA0 = inputs[i,1]
        CY0 = inputs[i,2]
        CZ0 = inputs[i,3]
        # calculate inlet molar flow rates
        nAin = CA0*Vdot
        nYin = CY0*Vdot
        nZin = CZ0*Vdot

        # define the reactor model residuals
        def reactor_model(x):
            # extract the test values for the unknowns
            nA = x[0]
            nY = x[1]
            nZ = x[2]

            # calculate the rate
            CA = nA/Vdot
            CY = nY/Vdot
            CZ = nZ/Vdot
            r = kf_guess*CA**2 - kr_guess*CY*CZ

            # evaluate the reactor model residuals
            r1 = nAin - nA - r*V
            r2 = nYin - nY + r*V
            r3 = nZin - nZ + r*V
            return [r1, r2, r3]
        
        # guess the solution to the reactor model
        nA_guess = CAout[i]*Vdot
        extent = nAin - nA_guess
        nY_guess = nYin + extent
        nZ_guess = nZin + extent
        guess = [nA_guess, nY_guess, nZ_guess]

        # solve the reactor model equations
        soln = sp.optimize.root(reactor_model,guess)
        if not(soln.success):
            print("A solution was NOT obtained:")
            print(soln.message)
        
        # calculate the response, CYout
        response[i] = soln.x[0]/Vdot
    return response
