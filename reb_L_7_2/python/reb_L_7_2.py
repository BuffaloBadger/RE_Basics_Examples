"""Calculations for Example L.7.2 of Reaction Engineering Basics"""

# import libraries
import numpy as np
import scipy as sp
import pandas as pd
from score_utils import solve_ivodes
from score_utils import fit_to_SR_data
import matplotlib.pyplot as plt

# given and known constants
CA_0 = 1.0 # M
CB_0 = 0.9 # M
V = 1.0 # L (basis)

# make Vdot available to all functions
k_current = float('nan')

# derivatives function
def derivatives(t,dep):
    # extract the dependent variables for this integration step
    nA = dep[0]
    nB = dep[1]

    # calculate the concentrations of A and B
    CA = nA/V
    CB = nB/V
    
    # evaluate the derivatives
    dnAdt = -k_current*V*CA*CB
    dnBdt = -k_current*V*CA*CB
    dnYdt = k_current*V*CA*CB
    dnZdt = k_current*V*CA*CB

    # return the derivatives
    return [dnAdt, dnBdt, dnYdt, dnZdt]

# BSTR model function
def profiles(t_f):
    # set the initial values
    nA_0 = CA_0*V
    nB_0 = CB_0*V
    ind_0 = 0.0
    dep_0 = np.array([nA_0, nB_0, 0.0, 0.0])

    # define the stopping criterion
    f_var = 0
    f_val = t_f
 
    # solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0
                            , f_var, f_val, derivatives)

    # check that a solution was found
    if not(success):
        print(f"The profiles were NOT found: {message}")

    # extract dependent variable profiles
    nA = dep[0,:]
    nB = dep[1,:]
    nY = dep[2,:]
    nZ = dep[3,:]

    # return all profiles
    return t, nA, nB, nY, nZ

# predicted responses function
def predicted_responses(t_f, beta):
    # set the current value of k
    global k_current
    k_current = 10**beta

    # get the number of experiments
    shape = np.shape(t_f)
    n_expt = shape[0]

    # allocate storage for the predicted responses
    CY_model = np.zeros(n_expt)

    # loop through the experiments
    for i in range(0,n_expt):
        # solve the BSTR design equations
        t, nA, nB, nY, nZ = profiles(t_f[i])

        # calculate the response
        CY_model[i] = nY[-1]/V
    
    return CY_model

# function that performs the calculations
def perform_the_calculations():
    # read the experimental data into a dataframe
    df = pd.read_csv("reb_L_7_2/reb_L_7_2_data.csv")

    # extract the data as arrays
    t_f = df["t (min)"].to_numpy()
    CY_expt = df["CY (M)"].to_numpy()

    # guess beta, the base 10 log of k
    guess = 0.0

    # fit the predicted responses model to the data
    beta, beta_ci, r_squared = fit_to_SR_data(guess, t_f, CY_expt, 
        predicted_responses, False)

    # extract the results
    k = 10**beta[0]
    k_lower = 10**beta_ci[0,0]
    k_upper = 10**beta_ci[0,1]

    # calculate the model-predicted responses and the residuals
    CY_model = predicted_responses(t_f,beta[0])
    residuals = CY_expt - CY_model

    # generate, show and save a parity plot
    plt.figure(1) 
    plt.plot(CY_expt, CY_model, color = 'k', marker='o', ls=''
             , fillstyle='none', markeredgewidth=2)
    plt.plot([np.min(CY_expt), np.max(CY_expt)]
             , [np.min(CY_expt), np.max(CY_expt)], color = 'r')
    plt.xlabel("Experimental C$_Y$ (M)")
    plt.ylabel("Predicted C$_Y$ (M)")
    plt.savefig('reb_L_7_2/python/parity_plot.png')
    plt.show()

    # generate, show and save a residuals plot
    plt.figure(2) 
    plt.plot(t_f, residuals, color = 'k', marker='o', ls='', fillstyle='none'
             , markeredgewidth=2)
    plt.axhline(y=0, color = 'r')
    plt.xlabel("Reaction time (min)")
    plt.ylabel("Residual (M)")
    plt.savefig('reb_L_7_2/python/residuals_vs_tf_plot.png')
    plt.show()

    # tabulate, show and save the results
    data = [['k', f'{k}', 'L mol^-1^ min^-1^'],
            ['k_lower_CI_limit',f'{k_lower}', 'L mol^-1^ min^-1^'],
            ['k_upper_CI_limit',f'{k_upper}', 'L mol^-1^ min^-1^'],
            ['R_squared',f'{r_squared}','']]
    result_df = pd.DataFrame(data, columns=['item','value','units'])
    result_df.to_csv("reb_L_7_2/python/results.csv", index=False)
    print(' ')
    print(f'k: {k:.3g} L/mol/min, 95% CI [{k_lower:.3g}, {k_upper:.3g}]')
    print(f'R-squared: {r_squared:.3g}')
    print(' ')

if __name__=="__main__":
    perform_the_calculations()