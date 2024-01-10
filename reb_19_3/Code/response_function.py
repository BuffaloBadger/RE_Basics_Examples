import numpy as np
import rebutils as reb

def response(inputs, param):
    # extract the rate coefficient
    k = 10.0**param

    # given and known constants
    V = 100.0 # cm^3
    P0 = 6.0 # atm
    T = 275 + 273.15 # K
    R = 82.06 # cm^3 atm/mol/K

    # allocate storage for the responses
    shape = np.shape(inputs)
    n_data = shape[0]
    response = np.zeros(n_data)

    # define the reactor model
    def reactor_model(t_current, n_current):
        # extract the molar amounts of A and B
        nA = n_current[0]
        nB = n_current[1]

        # calculate the rate
        PA = nA*R*T/V # eqn 6
        PB = nB*R*T/V # eqn 7
        if PA <= 0:
            r = 0.0
        elif PB <=0:
            r = 0.0
        else:
            r = k*PA*np.sqrt(PB) # eqn 2

        # evaluate the derivatives
        dnAdt = -r*V # eqn 3
        dnBdt = -r*V # eqn 4
        dnZdt = r*V # eqn 5

        # collect and return the derivatives
        ddt = np.array([dnAdt, dnBdt, dnZdt])
        return ddt

    # loop through the data points
    for i in range(0,n_data):
        # extract the adjusted inputs
        PA0 = inputs[i,0]
        tRxn = inputs[i,1]

        # solve the reactor model equations
        ind0 = 0.0
        nA0 = PA0*V/R/T # eqn 8
        PB0 = P0 - PA0 # eqn 10
        nB0 = PB0*V/R/T # eqn 9
        dep0 = [nA0, nB0, 0.0]
        f_var = 0
        f_val = tRxn
        soln = reb.solveIVODEs(ind0, dep0, f_var, f_val, reactor_model)
        if not(soln.success):
            print("Error when solving the design equations: " + soln.message)

        # calculate the response
        nAout = soln.y[0,-1]
        nBout = soln.y[1,-1]
        nZout = soln.y[2,-1]
        response[i] = (nAout + nBout + nZout)*R*T/V # eqn 11

    # return the full set of responses
    return response