import pandas as pd
import numpy as np
import rebutils as reb
import reb_19_2_response_fcn as rf
import matplotlib.pyplot as plt

# set filepaths
filepath_to_data = './reb_19_2/Data/'
filepath_to_results = './reb_19_2/Results/'
filepath_to_figures = '../RE_Basics/Graphics/'

# Constant inputs
V = 500.0 # cc
R = 1.987E-3 # kcal/mol/K
Rpv = 82.06 # cc atm/mol/K

# read the data filenames
df = pd.read_csv(filepath_to_data + 'reb_19_2_data_filenames.csv')
filenames = df['Data Filename'].to_numpy()

# create a dataframe for the fitting results
parameter_estimates = pd.DataFrame(columns=['T', 'k', 'k_ll', 'k_ul', 'R_sq'])

# process the data files
for filename in filenames:
    # read the file
    df = pd.read_csv(filepath_to_data + filename)
    n_data = len(df.index)

    # extract the data as arrays
    T = df['T'].to_numpy() + 273.15
    PA0 = df['PA0'].to_numpy()
    PB0 = df['PB0'].to_numpy()
    t = df['t'].to_numpy()
    f_A = df['fA'].to_numpy()
    T_K = T[0]

    # combine the adjusted inputs as a matrix
    adjusted_inputs = np.transpose(np.array([T, PA0, PB0, t]))

    # make a guess for log_10 of k
    par_guess = [-5.0]

    # check the guess
    check = rf.response(adjusted_inputs,par_guess[0])
    mean_fA = np.mean(check)
    print('\nMean conversion using guess: %.3g'%(mean_fA))

    # estimate log_10 of k
    beta, beta_ci, r_squared = reb.fitNLSR(par_guess, adjusted_inputs, f_A, 
        rf.response, False)
    
    # extract the results
    k = 10.**beta[0]
    k_ll = 10.**beta_ci[0,0]
    k_ul = 10.**beta_ci[0,1]
    
    # add the fitting results to the dataframe
    parameter_estimates.loc[len(parameter_estimates)] = [T_K - 273.15, k, k_ll,
            k_ul, r_squared]
    
    # calculate the model-predicted response and the residuals
    fA_model = rf.response(adjusted_inputs, beta[0])
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
    name = 'reb_19_2_parity_' + T_as_text + '.png'
    plt.savefig(filepath_to_results  + name)
    plt.savefig(filepath_to_figures + name,)
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
    name = 'reb_19_2_residual_PA0_' + T_as_text + '.png'
    plt.savefig(filepath_to_results  + name)
    plt.savefig(filepath_to_figures + name)
    plt.show()

    plt.figure() 
    plt.plot(PB0, residual, color = 'r', marker='o', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("$P_{B,0}$ (atm)")
    plt.ylabel("Residual")
    T_as_text = format(T_K-273.15,'.0f')
    plt.title('T = ' + T_as_text + ' °C')
    plt.tight_layout()
    name = 'reb_19_2_residual_PB0_' + T_as_text + '.png'
    plt.savefig(filepath_to_results  + name)
    plt.savefig(filepath_to_figures + name)
    plt.show()

    plt.figure() 
    plt.plot(t, residual, color = 'r', marker='o', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("$t_{rxn}$ (min)")
    plt.ylabel("Residual")
    T_as_text = format(T_K-273.15,'.0f')
    plt.title('T = ' + T_as_text + ' °C')
    plt.tight_layout()
    name = 'reb_19_2_residual_t_rxn_' + T_as_text + '.png'
    plt.savefig(filepath_to_results  + name)
    plt.savefig(filepath_to_figures + name)
    plt.show()

# show and save the fitting results
print('\nParameter Estimates:\n')
print(parameter_estimates)
print('\nThese results are being saved as ' + filepath_to_results + \
      'reb_19_2_resp_fcn_params.csv')
parameter_estimates.to_csv(filepath_to_results + 'reb_19_2_resp_fcn_params.csv',\
    index=False)