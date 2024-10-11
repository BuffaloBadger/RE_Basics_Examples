import pandas as pd
import numpy as np
from score_utils import solve_ivodes
from score_utils import fit_to_SR_data
import matplotlib.pyplot as plt

# given and known constants
P = 1.0 # atm
mTotal = 3.0 # g
VFR_in = 0.85 # L/min
T = 400. + 273.15 # K
K = 12.2

# globally available variables
global gk, ga, gb, gy, gz
gk = float('nan')
ga = float('nan')
gb = float('nan')
gy = float('nan')
gz = float('nan')

# derivatives function
def derivatives (m, dep):
    nA = dep[0]
    nB = dep[1]
    nY = dep[2]
    nZ = dep[3]
    ntot = nA + nB + nY + nZ
    PA = nA/ntot*P
    PB = nB/ntot*P
    PY = nY/ntot*P
    PZ = nZ/ntot*P
    r = gk* PA**ga* PB**gb* PY**gy* PZ**gz*\
        (1 - PY*PZ/K/PA/PB)
    ddm = np.array([-r, -r, r, r])
    return ddm

# PFR model function
def profiles(yA0, yB0, yY0, yZ0, k, alphaA, alphaB, alphaY, alphaZ):
    # make the rate expression parameters available to the derivatives function
    global gk, ga, gb, gy, gz
    gk = k
    ga = alphaA
    gb = alphaB
    gy = alphaY
    gz = alphaZ

    #initial values and stopping criterion
    ind0 = 0
    nA0 = yA0*VFR_in/22.4
    nB0 = yB0*VFR_in/22.4
    nY0 = yY0*VFR_in/22.4
    nZ0 = yZ0*VFR_in/22.4
    dep0 =[nA0, nB0, nY0, nZ0]
    stop_var = 0
    stop_val = mTotal

    # solve the PFR design equations
    m, dep, success, message = solve_ivodes(ind0, dep0, stop_var, stop_val
                                            , derivatives, False)

    # print a warning if there was a problem solving the design equations
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    nA = dep[0,:]
    nB = dep[1,:]
    nY = dep[2,:]
    nZ = dep[3,:]

    # return the profiles
    return m, nA, nB, nY, nZ

# predicted responses function
def predicted_responses(adj_inputs, kLog, alphaA, alphaB, alphaY, alphaZ):
    # allocate storage for the responses
    nExpts = len(adj_inputs)
    PA_1_model = np.zeros(nExpts)

    # extract k
    k = 10**kLog

    # loop through the data points
    for i, input in enumerate(adj_inputs):
        # extract the adjusted inputs
        yA0 = adj_inputs[i,0]
        yB0 = adj_inputs[i,1]
        yY0 = adj_inputs[i,2]
        yZ0 = adj_inputs[i,3]

        # solve the PFR design equations
        m, nA, nB, nY, nZ = profiles(yA0, yB0, yY0, yZ0, k, alphaA
                                     , alphaB, alphaY, alphaZ)
        
        # calculate the response
        nA_1 = nA[-1]
        nB_1 = nB[-1]
        nY_1 = nY[-1]
        nZ_1 = nZ[-1]
        PA_1_model[i] = nA_1/(nA_1 + nB_1 + nY_1 + nZ_1)*P
    
    # return the responses
    return PA_1_model

# quantities of interest function
def quantities_of_interest(adj_inputs, PA_1, par_guess):
    # estimate the parameters
    beta, beta_ci, r_squared = fit_to_SR_data(par_guess, adj_inputs 
        , PA_1, predicted_responses, False)

    # extract the results
    k = 10**beta[0]
    k_CI = 10**beta_ci[0,:]
    alphaA = beta[1]
    alphaA_CI = beta_ci[1,:]
    alphaB = beta[2]
    alphaB_CI = beta_ci[2,:]
    alphaY = beta[3]
    alphaY_CI = beta_ci[3,:]
    alphaZ = beta[4]
    alphaZ_CI = beta_ci[4,:]
        
    # calculate the model-predicted response and the residuals
    PA_1_model = predicted_responses(adj_inputs, beta[0], alphaA, alphaB, alphaY
                                     , alphaZ)
    epsilon_expt = PA_1 - PA_1_model



    return k, k_CI, alphaA, alphaA_CI, alphaB, alphaB_CI, alphaY, alphaY_CI\
         , alphaZ, alphaZ_CI, r_squared, PA_1_model, epsilon_expt

# master function
def perform_the_calculations():
    # read the data from the .csv file
    df = pd.read_csv('reb_21_5_2/reb_21_5_2_data.csv')
        # columns: yA, yB, yY, yZ, PA

    # extract the data as arrays
    yA0 = df['yA'].to_numpy()
    yB0 = df['yB'].to_numpy()
    yY0 = df['yY'].to_numpy()
    yZ0 = df['yZ'].to_numpy()
    PA_1 = df['PA'].to_numpy()

    # combine the adjusted inputs as a matrix
    adj_inputs = np.transpose(np.array([yA0, yB0, yY0, yZ0]))

    # make a guess the parameters
    #par_guess = [-4.0, 1.0, 1.0, 1.0, 1.0]
    par_guess = [np.log10(0.00126), 0.98, 0.32, -0.57, -0.0094]

    # calculate the quantities of interest
    k, k_CI, alphaA, alphaA_CI, alphaB, alphaB_CI, alphaY, alphaY_CI, alphaZ\
         , alphaZ_CI, r_squared, PA_1_model, epsilon_expt\
         = quantities_of_interest(adj_inputs, PA_1, par_guess)
    
    order = format(alphaA + alphaB + alphaY + alphaZ,'.2f')
    kUnits = 'mol g^-1^ min^-1^ atm^-' + order + '^'
    
    # save the results to a .csv file
    data = [['k', f'{k:.3g}', kUnits],
        ['k_lower_limit', f'{k_CI[0]:.3g}', kUnits],
        ['k_upper_limit', f'{k_CI[1]:.3g}', kUnits],
        ['alphaA', f'{alphaA:.2g}',''],
        ['alphaA_lower_limit', f'{alphaA_CI[0]:.2g}',''],
        ['alphaA_upper_limit', f'{alphaA_CI[1]:.2g}',''],
        ['alphaB', f'{alphaB:.2g}',''],
        ['alphaB_lower_limit', f'{alphaB_CI[0]:.2g}',''],
        ['alphaB_upper_limit', f'{alphaB_CI[1]:.2g}',''],
        ['alphaY', f'{alphaY:.2g}',''],
        ['alphaY_lower_limit', f'{alphaY_CI[0]:.2g}',''],
        ['alphaY_upper_limit', f'{alphaY_CI[1]:.2g}',''],
        ['alphaZ', f'{alphaZ:.2g}',''],
        ['alphaZ_lower_limit', f'{alphaZ_CI[0]:.2g}',''],
        ['alphaZ_upper_limit', f'{alphaZ_CI[1]:.2g}',''],
        ['R_squared', f'{r_squared:.3g}', '']]
    result = pd.DataFrame(data, columns=['item','value','units'])
    print(" ")
    print(result)
    result.to_csv("reb_21_5_2/python/reb_21_5_2_results.csv", index=False)

    # create, show, and save a parity plot
    plt.figure() 
    plt.plot(PA_1, PA_1_model, color = 'k', marker='o', ls=''
             , label = 'Data')
    plt.plot([min(PA_1),max(PA_1)],[min(PA_1),max(PA_1)]
             , color = 'r', ls = '-', label = 'Parity Line')
    plt.xlabel("$P_{A,1}$ (atm)")
    plt.ylabel("$P_{A,1,model}$ (atm)")
    plt.legend
    plt.tight_layout()
    plt.savefig('reb_21_5_2/python/reb_21_5_2_parity.png')
    plt.show()

    # create, show and save residuals plots
    plt.figure() 
    plt.plot(yA0, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$y_{A,0}$")
    plt.ylabel("Residual (atm)")
    plt.tight_layout()
    plt.savefig('reb_21_5_2/python/reb_21_5_2_A_residual.png')
    plt.show()

    plt.figure() 
    plt.plot(yB0, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$y_{B,0}$")
    plt.ylabel("Residual (atm)")
    plt.tight_layout()
    plt.savefig('reb_21_5_2/python/reb_21_5_2_B_residual.png')
    plt.show()

    plt.figure() 
    plt.plot(yY0, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$y_{Y,0}$")
    plt.ylabel("Residual (atm)")
    plt.tight_layout()
    plt.savefig('reb_21_5_2/python/reb_21_5_2_Y_residual.png')
    plt.show()

    plt.figure() 
    plt.plot(yZ0, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$y_{Z,0}$")
    plt.ylabel("Residual (atm)")
    plt.tight_layout()
    plt.savefig('reb_21_5_2/python/reb_21_5_2_Z_residual.png')
    plt.show()

if __name__=="__main__":
    perform_the_calculations()