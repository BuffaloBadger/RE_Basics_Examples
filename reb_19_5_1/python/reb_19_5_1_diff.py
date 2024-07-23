"""Calculations for Reaction Engineering Basics Example 19.5.1 using differential analysis"""

# import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm


# given and known constants
V = 1.0 # L
R = 8.314E-3 # kJ/mol/K

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

        # calculate x
        x = df['CAf'].to_numpy()*V

        # allocate storage for y
        y = np.zeros(len(df.index))

        # calculate y
        for i in range(0,len(x)):
            if tf[i] == 5.0:
                x_minus = CA0[i]*V
                t_minus = 0
            else:
                x_minus = x[i-1]
                t_minus = tf[i-1]
            y[i] = (x_minus - x[i])/(tf[i] - t_minus)
        
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
        plt.tight_layout()
        f_name = 'reb_19_5_1_model_' + T_as_text + '.png'
        plt.savefig('reb_19_5_1/python/' + f_name)
        plt.show()

    # show and save the fitting results
    print('\nParameter Estimates:\n')
    print(parameter_estimates)
    parameter_estimates.to_csv(
        'reb_19_5_1/python/reb_19_5_1_diff_params.csv', index=False)

if __name__=="__main__":
    perform_the_calculations()