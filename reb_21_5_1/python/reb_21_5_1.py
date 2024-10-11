import pandas as pd
import numpy as np
from score_utils import solve_ivodes
from score_utils import fit_to_SR_data
import matplotlib.pyplot as plt

# given and known constants
P = 1.0 # atm
D = 1.0 # cm
L = 10.0 # cm
T = 1500. # K
mol_per_cc_at_stp = 1.0E-3/22.4

# globally available variable
global gk
gk = float('nan')

# derivatives function
def derivatives(z,dep):
    # extract the dependent variables
    nDot_A = dep[0]
    nDot_Y = dep[1]
    nDot_Z = dep[2]

    # calculate the rate
    PA = nDot_A/(nDot_A + nDot_Y + nDot_Z)*P
    r = gk*PA

    # evaluate the derivatives
    dnAdz = -np.pi * (D**2/4)*r
    dnBdz = np.pi * (D**2/4)*r
    dnZdz = np.pi * (D**2/4)*r

    # collect and return the derivatives
    ddz = np.array([dnAdz, dnBdz, dnZdz])
    return ddz

# PFR model function
def profiles(Vdot_A0, Vdot_Y0, Vdot_Z0, beta):
    # make k available to the derivatives function
    global gk
    gk = 10**beta

    # initial values and stopping criterion
    ind0 = 0
    nA0 = Vdot_A0*mol_per_cc_at_stp
    nY0 = Vdot_Y0*mol_per_cc_at_stp
    nZ0 = Vdot_Z0*mol_per_cc_at_stp
    dep0 =[nA0, nY0, nZ0]
    stop_var = 0
    stop_val = L

    # solve the PFR design equations
    z, dep, success, message = solve_ivodes(ind0, dep0, stop_var, stop_val
                                            , derivatives, False)

    # print a warning if there was a problem solving the design equations
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")
    
    # extract the dependent variable profiles
    nA = dep[0,:]
    nY = dep[1,:]
    nZ = dep[2,:]

    # return the profiles
    return z, nA, nY, nZ

# predicted responses function
def predicted_responses(adj_inputs, beta):
    # allocate storage for the responses
    nExpts = len(adj_inputs)
    fA_model = np.zeros(nExpts)

    # loop through the data points
    for i, input in enumerate(adj_inputs):
        # get the adjusted inputs
        Vdot_A0 = input[0]
        Vdot_Y0 = input[1]
        Vdot_Z0 = input[2]

        # Solve the PFR design equations
        z, nA, nY, nZ = profiles(Vdot_A0, Vdot_Y0, Vdot_Z0, beta)

        # Calculate the response
        nA0 = Vdot_A0*mol_per_cc_at_stp
        fA_model[i] = 100*(nA0 - nA[-1])/nA0
    
    # return the responses
    return fA_model

# quantities of interest function
def quantities_of_interest(adj_inputs, fA, par_guess):

    # estimate log_10 of k
    beta, beta_ci, r_squared = fit_to_SR_data(par_guess, adj_inputs 
        , fA, predicted_responses, False)
    
    # extract the results
    k = 10.**beta[0]
    k_CI = 10.**beta_ci[0,:]
        
    # calculate the model-predicted response and the residuals
    fA_model = predicted_responses(adj_inputs, beta[0])
    epsilon_expt = fA - fA_model

    return k, k_CI, r_squared, fA_model, epsilon_expt

# master function
def perform_the_calculations():
    # read the data from the .csv file
    df = pd.read_csv('reb_21_5_1/reb_21_5_1_data.csv')
        # columns: Vdot_A0, Vdot_Y0, Vdot_Z0, fA

    # extract the data as arrays
    Vdot_A0 = df['A_in'].to_numpy()
    Vdot_Y0 = df['Y_in'].to_numpy()
    Vdot_Z0 = df['Z_in'].to_numpy()
    fA = df['fA'].to_numpy()

    # combine the adjusted inputs as a matrix
    adj_inputs = np.transpose(np.array([Vdot_A0, Vdot_Y0, Vdot_Z0]))

    # make a guess for log_10 of k
    beta_guess = [-2.0]

    # calculate the quantities of interest
    k, k_CI, r_squared, fA_model, epsilon_expt\
          = quantities_of_interest(adj_inputs, fA, beta_guess)
        
    # create, show, and save a parity plot
    plt.figure() 
    plt.plot(fA, fA_model, color = 'k', marker='o', ls=''
             , label = 'Data')
    plt.plot([min(fA),max(fA)],[min(fA),max(fA)]
             , color = 'r', ls = '-', label = 'Parity Line')
    plt.xlabel("$f_A$ (%)")
    plt.ylabel("$f_{A,model}$ (%)")
    plt.legend
    plt.tight_layout()
    plt.savefig('reb_21_5_1/python/reb_21_5_1_parity.png')
    plt.show()

    # create, show and save residuals plots
    plt.figure() 
    plt.plot(Vdot_A0, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$\dot{V}_{A,0}$ (sccm)")
    plt.ylabel("Residual (%)")
    plt.tight_layout()
    plt.savefig('reb_21_5_1/python/reb_21_5_1_A_residual.png')
    plt.show()

    plt.figure() 
    plt.plot(Vdot_Y0, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$\dot{V}_{Y,0}$ (sccm)")
    plt.ylabel("Residual (%)")
    plt.tight_layout()
    plt.savefig('reb_21_5_1/python/reb_21_5_1_Y_residual.png')
    plt.show()

    plt.figure() 
    plt.plot(Vdot_Z0, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$\dot{V}_{Z,0}$ (sccm)")
    plt.ylabel("Residual (%)")
    plt.tight_layout()
    plt.savefig('reb_21_5_1/python/reb_21_5_1_Z_residual.png')
    plt.show()

    # save the results to a .csv file
    data = [['k', f'{k:.3g}', 'mol cm^-3^ min^-1^ atm^-1^'],
        ['k_lower_limit', f'{k_CI[0]:.3g}', 'mol cm^-3^ min^-1^ atm^-1^'],
        ['k_upper_limit', f'{k_CI[1]:.3g}', 'mol cm^-3^ min^-1^ atm^-1^'],
        ['R_squared', f'{r_squared:.3g}', '']]
    result = pd.DataFrame(data, columns=['item','value','units'])
    print(" ")
    print(result)
    result.to_csv("reb_21_5_1/python/reb_21_5_1_results.csv", index=False)

if __name__=="__main__":
    perform_the_calculations()