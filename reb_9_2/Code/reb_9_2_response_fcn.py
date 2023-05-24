import numpy as np
import rebutils as reb

def response(adj_inputs, k_log_10):
    # given or known
    V = 500.0 # cc
    Rpv = 82.06 # cc atm/mol/K

    # allocate storage for the responses
    nExpts = len(adj_inputs)
    resp = np.zeros(nExpts)

    # get the rate coefficient
    k = 10.**k_log_10

    # loop through all of the sets of adjusted inputs
    for i, input in enumerate(adj_inputs):
        # get the adjusted inputs
        T = input[0]
        PA0 = input[1]
        PB0 = input[2]
        tRxn = input[3]

        # define the reactor model
        def reactor_model(t, dep):
            # get the dependent variables
            nA = dep[0]
            nB = dep[1]
            nY = dep[2]
            nZ = dep[3]

            # calculate the partial pressures, eqns 7 and 8
            PA = nA*Rpv*T/V
            PB = nB*Rpv*T/V

            # calculate the rate, eqn 2
            r = k*PA*PB

            # calculate the time derivatives of the dependent variables
            dnAdt = -r*V
            dnBdt = -r*V
            dnYdt = r*V
            dnZdt = r*V

            # return an array containing the derivatives
            ddt = np.array([dnAdt, dnBdt, dnYdt, dnZdt])
            return ddt
        
        # set the initial values, eqns 9 and 10
        t0 = 0
        nA0 = PA0*V/Rpv/T
        nB0 = PB0*V/Rpv/T
        n0 = np.array([nA0, nB0, 0.0, 0.0])

        # set the stopping criterion
        stop_var = 0
        stop_val = tRxn

        # solve the reactor design equations
        soln = reb.solveIVODEs(t0, n0, stop_var, stop_val, reactor_model)

        # print a warning if there was a problem solving the design equations
        if not(soln.success):
            print('Error solving the design equations: ' + soln.message)
        
        # calculate the model-predicted response, eqn 11
        nAf = soln.y[0,-1]
        fA = (nA0 - nAf)/nA0
        resp[i] = fA

    # return the responses
    return resp
