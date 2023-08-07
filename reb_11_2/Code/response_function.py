import pandas as pd
import numpy as np
import scipy as sp

# define the response function
def response_function(inputs, log_k_guess, alphaA_guess, alphaB_guess):

    # Given and known constants
    P = 3.0 # atm
    T = 450 + 273.15 # K
    V = 1.0 # L (basis)
    R = 0.08206 # L*atm/mol/K

    # extract the guess for k
    k_guess = 10**log_k_guess

    # allocate storage for the responses
    shape = np.shape(inputs)
    n_data = shape[0]
    response = np.zeros(n_data)

    # read the experimental responses to use as guesses when solving the reactor
    # design equations
    df = pd.read_csv('./reb_11_2/Data/reb_11_2_data.csv')
    CZout = df['CZout'].to_numpy()/1000.0

    for i in range(0,n_data):
        # extract the inputs
        space_time = inputs[i,0]
        yA0 = inputs[i,1]
        # calculate inlet molar flow rates
        VFRin = V/space_time
        nAin = yA0*P*VFRin/R/T
        nBin = (1-yA0)*P*VFRin/R/T

        # define the reactor model residuals
        def reactor_model(x):
            # extract the test values for the unknowns
            nA = x[0]
            nB = x[1]
            nZ = x[2]

            # calculate the rate
            PA = nA*P/(nA+nB+nZ)
            PB = nB*P/(nA+nB+nZ)
            r = k_guess*PA**alphaA_guess*PB**alphaB_guess

            # evaluate the reactor model residuals
            r1 = nAin - nA - r*V
            r2 = nBin - nB - r*V
            r3 = -nZ + r*V
            return [r1, r2, r3]
        
        # guess the solution to the reactor model
        nZ_guess = CZout[i]*VFRin
        nA_guess = nAin - nZ_guess
        nB_guess = nBin - nZ_guess
        guess = [nA_guess, nB_guess, nZ_guess]

        # solve the reactor model equations
        soln = sp.optimize.root(reactor_model,guess)
        if not(soln.success):
            print("A solution was NOT obtained:")
            print(soln.message)
        
        # calculate the response, CYout
        nAf = soln.x[0]
        nBf = soln.x[1]
        nZf = soln.x[2]
        response[i] = nZf*P/(nAf + nBf + nZf)/R/T
    return response
