import pandas as pd
import numpy as np
import rebutils as reb
import reb_9_1_response_function as rf
import matplotlib.pyplot as plt

# set filepaths
filepath_to_data = './reb_9_1/Data/'
filepath_to_results = './reb_9_1/Results/'
filepath_to_figures = '../RE_Basics/Graphics/'

# given or known
V = 1.0 # L
R = 8.314E-3 # kJ/mol/K

# read the data filenames
df = pd.read_csv(filepath_to_data + 'reb_9_1_data_filenames.csv')
filenames = df['Data Filename'].to_numpy()

# create a dataframe for the fitting results
parameter_estimates = pd.DataFrame(columns=['T', 'k', 'k_ll', 'k_ul', 'R_sq'])

# process the data files
for filename in filenames:
    # read the file
    df = pd.read_csv(filepath_to_data + filename)
    n_data = len(df.index)

    # extract the data as arrays
    T = df['T'].to_numpy()
    CA0 = df['CA0'].to_numpy()
    t = df['t'].to_numpy()
    CAexpt = df['CA'].to_numpy()
    T_K = T[0] + 273.15

    # combine the adjusted inputs as a matrix
    adjusted_inputs = np.transpose(np.array([T, CA0, t]))

    # make a guess for log_10 of k
    par_guess = [0.0]

    # estimate log_10 of k
    beta, beta_ci, r_squared = reb.fitNLSR(par_guess, adjusted_inputs, CAexpt, 
        rf.response, False)
    
    # extract the results
    k = 10.**beta[0]
    k_ll = 10.**beta_ci[0,0]
    k_ul = 10.**beta_ci[0,1]
    
    # add the fitting results to the dataframe
    parameter_estimates.loc[len(parameter_estimates)] = [T_K - 273.15, k, k_ll,
            k_ul, r_squared]
    
    # calculate the model-predicted y and the residuals
    y_model = rf.response(adjusted_inputs, beta[0])
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
    name = 'reb_9_1_resp_fcn_parity_' + T_as_text + '.png'
    plt.savefig(filepath_to_results  + name)
    plt.savefig(filepath_to_figures + name,)
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
    name = 'reb_9_1_resp_fcn_residual_CA0_' + T_as_text + '.png'
    plt.savefig(filepath_to_results  + name)
    plt.savefig(filepath_to_figures + name)
    plt.show()

    plt.figure() 
    plt.plot(t, residual, color = 'r', marker='o', ls='')
    plt.axhline(y=0, color = 'k')
    plt.xlabel("$t_{rxn}$ (min)")
    plt.ylabel("Residual (M)")
    T_as_text = format(T_K-273.15,'.0f')
    plt.title('T = ' + T_as_text + ' °C')
    plt.tight_layout()
    name = 'reb_9_1_resp_fcn_residual_t_rxn_' + T_as_text + '.png'
    plt.savefig(filepath_to_results  + name)
    plt.savefig(filepath_to_figures + name)
    plt.show()

# show and save the fitting results
print('\nParameter Estimates:\n')
print(parameter_estimates)
print('\nThese results are being saved as ' + filepath_to_results + \
      'reb_9_1_resp_fcn_params.csv')
parameter_estimates.to_csv(filepath_to_results + 'reb_9_1_resp_fcn_params.csv',\
    index=False)
