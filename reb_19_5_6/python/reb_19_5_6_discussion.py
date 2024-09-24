"""Calculations for Reaction Engineering Basics Example 19.5.6 Discussion"""

#import libraries
import pandas as pd
import numpy as np
from score_utils import solve_ivodes
from score_utils import fit_to_SR_data
from score_utils import Arrhenius_parameters
import matplotlib.pyplot as plt

# given and known constants
V = 1.0 # L
R = 8.314E-3 # kJ/mol/K

# make k_current available to all functions
k0_current = float('nan')
E_current = float('nan')
T_current = float('nan')

# derivatives function
def derivatives(t,dep):
    # get the dependent variables
    nA = dep[0]
    nZ = dep[1]

    # calculate the concentration
    CA = nA/V

    # calculate the rate
    k = k0_current*np.exp(-E_current/R/T_current)
    r = k*CA

    # calculate the time derivatives of the dependent variables
    dnAdt = -r*V
    dnZdt = r*V

    # return an array containing the derivatives
    ddt = np.array([dnAdt, dnZdt])
    return ddt

# BSTR model function
def profiles(CA0,tf):
    # set initial values and stopping criterion
    t0 = 0
    nA0 = CA0*V
    n0 = np.array([nA0, 0.0])
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
def predicted_responses(adj_inputs, k0, E):
    # allocate storage for the responses
    nExpts = len(adj_inputs)
    resp = np.zeros(nExpts)

    # make the rate coefficient available to other functions
    global k0_current, E_current
    k0_current = k0
    E_current = E

    # loop through the experiments in the data set
    for i, input in enumerate(adj_inputs):
        # get the adjusted inputs
        global T_current
        T_current = input[0]
        CA0 = input[1]
        tf = input[2]        

        # solve the reactor design equations
        t, nA, nZ = profiles(CA0,tf)
        
        # calculate the model-predicted response
        nAf = nA[-1]
        CAf = nAf/V
        resp[i] = CAf

    # return the responses
    return resp

# function that performs the calculations
def perform_the_calculations():
    # read the data from the .csv file
    df = pd.read_csv('reb_19_5_1/reb_19_5_1_data.csv')

    # extract the data as arrays
    T = df['T'].to_numpy()
    CA0 = df['CA0'].to_numpy()
    tf = df['tf'].to_numpy()
    CAf = df['CAf'].to_numpy()
    T_K = T + 273.15

    # combine the adjusted inputs as a matrix
    adjusted_inputs = np.transpose(np.array([T_K, CA0, tf]))

    # estimated parameters
    k0_resp_fcn = 3.67E8
    E_resp_fcn = 67.6
    k0_excel = np.exp(21.685)
    E_excel = 72.936
        
    # calculate the model-predicted y
    CAf_resp_fcn = predicted_responses(adjusted_inputs, k0_resp_fcn, E_resp_fcn)
    CAf_excel = predicted_responses(adjusted_inputs, k0_excel, E_excel)
        
    # create, show, and save a parity plot
    plt.figure() 
    plt.scatter(CAf, CAf_resp_fcn, facecolor='none', edgecolor = 'b', marker='o'
                , label='Full Model')
    plt.scatter(CAf,CAf_excel, facecolors='none', edgecolors =  'r', marker='o'
                , label = 'Approximate Model')
    plt.plot([min(CAf),max(CAf)],[min(CAf),max(CAf)], color = 'k', ls = '-'
             , label='Parity Line')
    plt.xlabel("$C_{A, expt}$ (M)")
    plt.ylabel("$C_{A, model}$ (M)")
    plt.legend()
    plt.tight_layout()
    plt.savefig('reb_19_5_6/python/reb_19_5_6_parity.png')
    plt.show()

if __name__=="__main__":
    perform_the_calculations()