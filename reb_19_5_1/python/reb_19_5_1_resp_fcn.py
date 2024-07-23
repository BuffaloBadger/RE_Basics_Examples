"""Calculations for Reaction Engineering Basics Example 19.5.1 using a response function"""

#import libraries
import pandas as pd
import numpy as np
import scipy as sp
from score_utils import solve_ivodes
from score_utils import fit_to_SR_data
import matplotlib.pyplot as plt

# given and known constants
V = 1.0 # L
R = 8.314E-3 # kJ/mol/K

# make k_current available to all functions
k_current = float('nan')

# derivatives function
def derivatives(t,dep):
    # get the dependent variables
    nA = dep[0]
    nZ = dep[1]

    # calculate the concentration
    CA = nA/V

    # calculate the rate
    r = k_current*CA

    # calculate the time derivatives of the dependent variables
    dnAdt = -r*V
    dnZdt = r*V

    # return an array containing the derivatives
    ddt = np.array([dnAdt, dnZdt])
    return ddt

# BSTR model function
def profiles(n0,tf):
    # set initial values and stopping criterion
    t0 = 0
    stop_var = 0

    # solve the BSTR reactor design equations
    t, dep, success, message = solve_ivodes(t0, n0, stop_var, tf, derivatives
                                            ,True)

    # print a warning if there was a problem solving the design equations
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")
    
    # extract the dependent variable profiles
    nA = dep[0,:]
    nZ = dep[1,:]

    # return the profiles
    return t, nA, nZ

# predicted responses function
def predicted_responses(adj_inputs, k_log_10):
    # allocate storage for the responses
    nExpts = len(adj_inputs)
    resp = np.zeros(nExpts)

    # make the rate coefficient available to other functions
    global k_current
    k_current = 10.**k_log_10

    # loop through the experiments in the data set
    for i, input in enumerate(adj_inputs):
        # get the adjusted inputs
        CA0 = input[1]
        tf = input[2]
        
        # set the initial values
        nA0 = CA0*V
        n0 = np.array([nA0, 0.0])

        # solve the reactor design equations
        t, nA, nZ = profiles(n0,tf)
        
        # calculate the model-predicted response
        nAf = nA[-1]
        CAf = nAf/V
        resp[i] = CAf

    # return the responses
    return resp

# function that performs the calculations
def perform_the_calculations():
    # read the data filenames
    df = pd.read_csv('reb_19_5_1/python/reb_19_5_1_data_filenames.csv')
    filenames = df['Data Filename'].to_numpy()

    # create a dataframe for the fitting results
    parameter_estimates = pd.DataFrame(columns=['T', 'k', 'k_ll', 'k_ul', 'R_sq'])

    # process the data blocks
    for filename in filenames:
        # read the file
        df = pd.read_csv('reb_19_5_1/python/' + filename)

        # extract the data as arrays
        T = df['T'].to_numpy()
        CA0 = df['CA0'].to_numpy()
        tf = df['tf'].to_numpy()
        CAexpt = df['CAf'].to_numpy()
        T_K = T[0] + 273.15

        # combine the adjusted inputs as a matrix
        adjusted_inputs = np.transpose(np.array([T, CA0, tf]))

        # make a guess for log_10 of k
        par_guess = [0.0]

        # estimate log_10 of k
        beta, beta_ci, r_squared = fit_to_SR_data(par_guess, adjusted_inputs, CAexpt, 
            predicted_responses, False)
        
        # extract the results
        k = 10.**beta[0]
        k_ll = 10.**beta_ci[0,0]
        k_ul = 10.**beta_ci[0,1]
        
        # add the fitting results to the dataframe
        parameter_estimates.loc[len(parameter_estimates)] = [T_K - 273.15, k, k_ll,
                k_ul, r_squared]
        
        # calculate the model-predicted y and the residuals
        y_model = predicted_responses(adjusted_inputs, beta[0])
        residual = CAexpt - y_model
        
        # create, show, and save a parity plot
        plt.figure() 
        plt.plot(CAexpt, y_model, color = 'r', marker='o', ls='')
        plt.plot([min(CAexpt),max(CAexpt)],[min(CAexpt),max(CAexpt)], \
            color = 'k', ls = '-')
        plt.xlabel("$C_{A, expt}$ (M)")
        plt.ylabel("$C_{A, model}$ (M)")
        T_as_text = format(T_K-273.15,'.0f')
        plt.title('T = ' + T_as_text + ' °C')
        plt.tight_layout()
        name = 'reb_19_5_1_resp_fcn_parity_' + T_as_text + '.png'
        plt.savefig('reb_19_5_1/python/' + name)
        plt.show()

        # create, show and save residuals plots
        plt.figure() 
        plt.plot(CA0, residual, color = 'r', marker='o', ls='')
        plt.axhline(y=0, color = 'k')
        plt.xlabel("$C_{A,0}$ (M)")
        plt.ylabel("Residual (M)")
        T_as_text = format(T_K-273.15,'.0f')
        plt.title('T = ' + T_as_text + ' °C')
        plt.tight_layout()
        name = 'reb_19_5_1_resp_fcn_residual_CA0_' + T_as_text + '.png'
        plt.savefig('reb_19_5_1/python/' + name)
        plt.show()

        plt.figure() 
        plt.plot(tf, residual, color = 'r', marker='o', ls='')
        plt.axhline(y=0, color = 'k')
        plt.xlabel("$t_{rxn}$ (min)")
        plt.ylabel("Residual (M)")
        T_as_text = format(T_K-273.15,'.0f')
        plt.title('T = ' + T_as_text + ' °C')
        plt.tight_layout()
        name = 'reb_19_5_1_resp_fcn_residual_tf_' + T_as_text + '.png'
        plt.savefig('reb_19_5_1/python/' + name)
        plt.show()
    
    # show and save the fitting results
    print('\nParameter Estimates:\n')
    print(parameter_estimates)
    parameter_estimates.to_csv(
        'reb_19_5_1/python/reb_19_5_1_resp_fcn_params.csv', index=False)

if __name__=="__main__":
    perform_the_calculations()