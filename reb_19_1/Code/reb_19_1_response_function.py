import numpy as np
import rebutils as reb

def response(adj_inputs, k_log_10):
    # given or known
    V = 1.0 # L

    # allocate storage for the responses
    nExpts = len(adj_inputs)
    resp = np.zeros(nExpts)

    # get the rate coefficient
    k = 10.**k_log_10

    # loop through all of the sets of adjusted inputs
    for i, input in enumerate(adj_inputs):
        # get the adjusted inputs
        T = input[0]
        CA0 = input[1]
        tRxn = input[2]

        # define the reactor model
        def reactor_model(t, dep):
            # get the dependent variables
            nA = dep[0]
            nZ = dep[1]

            # calculate the concentration, eqn 5
            CA = nA/V

            # calculate the rate, eqn 2
            r = k*CA

            # calculate the time derivatives of the dependent variables
            dnAdt = -r*V
            dnZdt = r*V

            # return an array containing the derivatives
            ddt = np.array([dnAdt, dnZdt])
            return ddt
        
        # set the initial values
        t0 = 0
        nA0 = CA0*V
        n0 = np.array([nA0, 0.0])

        # set the stopping criterion
        stop_var = 0
        stop_val = tRxn

        # solve the reactor design equations
        soln = reb.solveIVODEs(t0, n0, stop_var, stop_val, reactor_model)

        # print a warning if there was a problem solving the design equations
        if not(soln.success):
            print('Error solving the design equations: ' + soln.message)
        
        # calculate the model-predicted response
        nAf = soln.y[0,-1]
        CAf = nAf/V
        resp[i] = CAf

    # return the responses
    return resp
