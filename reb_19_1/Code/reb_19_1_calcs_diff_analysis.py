import pandas as pd
import numpy as np
import rebutils as reb
import matplotlib.pyplot as plt

# filepaths
filepath_to_data = './reb_19_1/Data/'
filepath_to_results = './reb_19_1/Results/'
filepath_to_figures = '../RE_Basics/Graphics/'

# given or known
V = 1.0 # L
R = 8.314E-3 # kJ/mol/K

# read the data filenames
df = pd.read_csv(filepath_to_data + 'reb_19_1_data_filenames.csv')
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
    T_K = T[0] + 273.15
    CA0 = df['CA0'].to_numpy()
    t = df['t'].to_numpy()
    x = df['CA'].to_numpy()*V

    # allocate storage for y
    y = np.zeros(len(df.index))
    
    # calculate y noting that x_i = V*CA_i
    for i in range(0,len(x)):
        if t[i] == 5.0:
            x_minus = CA0[i]*V
            t_minus = 0
        else:
            x_minus = x[i-1]
            t_minus = t[i-1]
        y[i] = (x_minus - x[i])/(t[i] - t_minus)

    # fit y = mx to the data
    beta, beta_ci, r_squared = reb.fitLinSR(x,y,False)

    # add the fitting results to the dataframe
    parameter_estimates.loc[len(parameter_estimates)] = [T_K - 273.15, beta[0],
            beta_ci[0,0], beta_ci[0,1], r_squared]
    
    # calculate the model-predicted y
    y_model = beta[0]*x
    
    # create, show, and save a model plot
    plt.figure() 
    plt.plot(x, y_model, color = 'k', ls='-' , label = 'model')
    plt.plot(x, y, color = 'r', marker='o', ls='', label = 'experiment')
    plt.xlabel("x")
    plt.ylabel("y")
    T_as_text = format(T_K-273.15,'.0f')
    plt.legend(title='T = ' + T_as_text + '°C')
    plt.tight_layout()
    f_name = 'reb_19_1_model_' + T_as_text + '.png'
    plt.savefig(filepath_to_results  + f_name)
    plt.savefig(filepath_to_figures + f_name)
    plt.show()

# show and save the fitting results
print('\nParameter Estimates:\n')
print(parameter_estimates)
print('\nThese results are being saved as ' + filepath_to_results + \
      'reb_19_1_diff_params.csv')
parameter_estimates.to_csv(filepath_to_results + 'reb_19_1_diff_params.csv',\
    index=False)

