import numpy as np
import rebutils as reb

def response(inputs, log_Vmax_guess, log_Km_guess):
    # given or known
    V = 50.0E-3 # L

    # extract V_max_guess and K_m_guess
    V_max_guess = 10.0**log_Vmax_guess
    K_m_guess = 10.0**log_Km_guess

    # define the reactor model
    def reactor_model(t,y):
        nS = y[0]
        CS = nS/V # eqn 5
        r = V_max_guess*CS/(K_m_guess + CS)
        dnSdt = -r*V
        dnPdt = r*V
        ddt = np.array([dnSdt, dnPdt])
        return ddt

    # allocate storage for the responses
    shape = np.shape(inputs)
    n_data = shape[0]
    resp = np.zeros(n_data)

    # loop through the data points
    for i in range(0,n_data):
        # get the adjusted inputs
        CS0 = inputs[i,0]
        tRxn = inputs[i,1]
    
        # Solve the reactor model equations
        t0 = 0.0
        nS0 = CS0*V # eqn 6
        n0 = np.array([nS0, 0.0])
        f_var = 0
        f_val = tRxn
        soln = reb.solveIVODEs(t0, n0, f_var, f_val, reactor_model)
        if not(soln.success):
            print("A solution was NOT obtained:")
            print(soln.message)
        # Calculate the response
        nPf = soln.y[1,-1]
        resp[i] = nPf/V # eqn 7
    return resp