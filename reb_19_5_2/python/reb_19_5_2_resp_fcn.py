"""Calculations for Reaction Engineering Basics Example 19.5.2"""

# import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from score_utils import solve_ivodes
from score_utils import fit_to_SR_data

# given and known constants
V = 500.0 # cc
Re = 1.987E-3 # kcal/mol/K
Rpv = 82.06 # cc atm/mol/K

# make k_current and T_current available to all functions
k_current = float('nan')
T_current = float('nan')

# derivatives function
def derivatives(t,dep):
    # get the dependent variables
    nA = dep[0]
    nB = dep[1]
    nY = dep[2]
    nZ = dep[3]

    # calculate the partial pressures
    PA = nA*Rpv*T_current/V
    PB = nB*Rpv*T_current/V

    # calculate the rate
    r = k_current*PA*PB

    # calculate the time derivatives of the dependent variables
    dnAdt = -r*V
    dnBdt = -r*V
    dnYdt = r*V
    dnZdt = r*V

    # return an array containing the derivatives
    ddt = np.array([dnAdt, dnBdt, dnYdt, dnZdt])
    return ddt

# BSTR model function
def profiles(n0,tf):
    # set the initial values and stopping criterion
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
    nB = dep[1,:]
    nY = dep[2,:]
    nZ = dep[3,:]

    # return the profiles
    return t, nA, nB, nY, nZ

# predicted responses function
def predicted_responses(adj_inputs, k_log_10):
    # allocate storage for the responses
    nExpts = len(adj_inputs)
    resp = np.zeros(nExpts)

    # make the rate coefficient available to all functions
    global k_current
    k_current = 10.**k_log_10

    # loop through all of the experiments
    for i, input in enumerate(adj_inputs):
        # get the adjusted inputs and make T available to all functions
        global T_current
        T_current = input[0]
        PA0 = input[1]
        PB0 = input[2]
        tf = input[3]

        # set the initial values
        nA0 = PA0*V/Rpv/T_current
        nB0 = PB0*V/Rpv/T_current
        n0 = np.array([nA0, nB0, 0.0, 0.0])

        # solve the BSTR design equations
        t, nA, nB, nY, nZ = profiles(n0,tf)

        # calculate the response
        nAf = nA[-1]
        fA = (nA0 - nAf)/nA0
        resp[i] = fA
    
    # return the responses
    return resp

# function that performs the calculations
def perform_the_calculations():
    # read the data filenames
    df = pd.read_csv('reb_19_5_2/python/reb_19_5_2_data_filenames.csv')
    filenames = df['Data Filename'].to_numpy()

    # create a dataframe for the fitting results
    parameter_estimates = pd.DataFrame(columns=['T', 'k', 'k_ll', 'k_ul', 'R_sq'])

    # process the data files
    for filename in filenames:
        # read the file
        df = pd.read_csv('reb_19_5_2/python/' + filename)

        # extract the data as arrays
        T = df['T'].to_numpy() + 273.15
        PA0 = df['PA0'].to_numpy()
        PB0 = df['PB0'].to_numpy()
        tf = df['tf'].to_numpy()
        f_A = df['fA'].to_numpy()
        T_K = T[0]

        # combine the adjusted inputs as a matrix
        adjusted_inputs = np.transpose(np.array([T, PA0, PB0, tf]))

        # make a guess for log_10 of k
        par_guess = [-5.0]

        # check the guess
        check = predicted_responses(adjusted_inputs,par_guess[0])
        mean_fA = np.mean(check)
        print('\nMean conversion using guess: %.3g'%(mean_fA))

        # estimate log_10 of k
        beta, beta_ci, r_squared = fit_to_SR_data(par_guess, adjusted_inputs, f_A, 
            predicted_responses, False)
        
        # extract the results
        k = 10.**beta[0]
        k_ll = 10.**beta_ci[0,0]
        k_ul = 10.**beta_ci[0,1]
        
        # add the fitting results to the dataframe
        parameter_estimates.loc[len(parameter_estimates)] = [T_K - 273.15, k, k_ll,
                k_ul, r_squared]
        
        # calculate the model-predicted response and the residuals
        fA_model = predicted_responses(adjusted_inputs, beta[0])
        residual = f_A - fA_model
        
        # create, show, and save a parity plot
        plt.figure() 
        plt.plot(f_A, fA_model, color = 'r', marker='o', ls='')
        plt.plot([min(f_A),max(f_A)],[min(f_A),max(f_A)], \
            color = 'k', ls = '-')
        plt.xlabel("$f_{A, expt}$")
        plt.ylabel("$f_{A, model}$")
        T_as_text = format(T_K-273.15,'.0f')
        plt.title('T = ' + T_as_text + ' °C')
        plt.tight_layout()
        name = 'reb_19_5_2_parity_' + T_as_text + '.png'
        plt.savefig('reb_19_5_2/python/' + name,)
        plt.show()

        # create, show and save residuals plots
        plt.figure() 
        plt.plot(PA0, residual, color = 'r', marker='o', ls='')
        plt.axhline(y=0, color = 'k')
        plt.xlabel("$P_{A,0}$ (atm)")
        plt.ylabel("Residual")
        T_as_text = format(T_K-273.15,'.0f')
        plt.title('T = ' + T_as_text + ' °C')
        plt.tight_layout()
        name = 'reb_19_5_2_residual_PA0_' + T_as_text + '.png'
        plt.savefig('reb_19_5_2/python/' + name)
        plt.show()

        plt.figure() 
        plt.plot(PB0, residual, color = 'r', marker='o', ls='')
        plt.axhline(y=0, color = 'k')
        plt.xlabel("$P_{B,0}$ (atm)")
        plt.ylabel("Residual")
        T_as_text = format(T_K-273.15,'.0f')
        plt.title('T = ' + T_as_text + ' °C')
        plt.tight_layout()
        name = 'reb_19_5_2_residual_PB0_' + T_as_text + '.png'
        plt.savefig('reb_19_5_2/python/' + name)
        plt.show()

        plt.figure() 
        plt.plot(tf, residual, color = 'r', marker='o', ls='')
        plt.axhline(y=0, color = 'k')
        plt.xlabel("$t_{rxn}$ (min)")
        plt.ylabel("Residual")
        T_as_text = format(T_K-273.15,'.0f')
        plt.title('T = ' + T_as_text + ' °C')
        plt.tight_layout()
        name = 'reb_19_5_2_residual_tf_' + T_as_text + '.png'
        plt.savefig('reb_19_5_2/python/' + name)
        plt.show()

    # show and save the fitting results
    print('\nParameter Estimates:\n')
    print(parameter_estimates)
    parameter_estimates.to_csv('reb_19_5_2/python/reb_19_5_2_resp_fcn_params.csv',\
        index=False)

if __name__=="__main__":
    perform_the_calculations()