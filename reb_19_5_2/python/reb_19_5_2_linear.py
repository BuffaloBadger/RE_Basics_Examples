"""Calculations for Reaction Engineering Basics Example 19.5.2"""

# import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

# given or known
V = 500.0 # cc
R = 1.987E-3 # kcal/mol/K
Rpv = 82.06 # cc atm/mol/K

# read the data filenames
df = pd.read_csv('reb_19_5_2/python/reb_19_5_2_data_filenames.csv')
filenames = df['Data Filename'].to_numpy()

# create a dataframe for the fitting results
parameter_estimates = pd.DataFrame(columns=['T', 'k', 'k_ll', 'k_ul', 'R_sq'])

# process the data files
for filename in filenames:
    # read the file
    df = pd.read_csv('reb_19_5_2/python/' + filename)
    n_data = len(df.index)

    # extract the data as arrays
    T = df['T'].to_numpy() + 273.15
    T_K = T[0]
    PA0 = df['PA0'].to_numpy()
    PB0 = df['PB0'].to_numpy()
    t = df['tf'].to_numpy()
    f_A = df['fA'].to_numpy()

    # allocate storage for x and y
    x = np.zeros(n_data)
    y = np.zeros(n_data)
    
    # calculate x and y 
    for i in range(0,n_data):
        nA0 = PA0[i]*V/Rpv/T_K
        nB0 = PB0[i]*V/Rpv/T_K
        nA = nA0*(1-f_A[i])
        nB = nB0 - nA0 + nA
        x[i] = -t[i]*(Rpv*T_K)**2/V
        if nA0 == nB0:
            y[i] = 1/nA0 - 1/nA
        else:
            y[i] = 1/(nA0-nB0)*np.log(nA0*nB/nB0/nA)

    # fit y = mx to the data
        model = sm.OLS(y,x)
        res = model.fit()
        m = res.params[0]
        m_ci = res.conf_int(alpha=0.05)

        # add the fitting results to the dataframe
        parameter_estimates.loc[len(parameter_estimates)] = [T_K - 273.15
            , m, m_ci[0,0], m_ci[0,1], res.rsquared]
    
        # calculate the model-predicted y
        y_model = m*x
    
    # create, show, and save a model plot
    plt.figure() 
    plt.plot(x, y_model, color = 'k', ls='-' , label = 'model')
    plt.plot(x, y, color = 'r', marker='o', ls='', label = 'experiment')
    plt.xlabel("x")
    plt.ylabel("y")
    T_as_text = format(T_K-273.15,'.0f')
    plt.legend(title='T = ' + T_as_text + '°C')
    plt.ticklabel_format(axis='x',style='plain')
    plt.xticks(rotation=25)
    plt.tight_layout()
    filename = 'reb_19_5_2_linear_' + T_as_text + '.png'
    plt.savefig('reb_19_5_2/python/' + filename)
    plt.show()

# show and save the results
print('\nParameter Estimates:\n')
print(parameter_estimates)
parameter_estimates.to_csv('reb_19_5_2/python/reb_19_5_2_lin_params.csv',\
    index=False)