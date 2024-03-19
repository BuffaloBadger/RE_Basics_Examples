"""Calculations for Reaction Engineering Basics Example 9.6.4"""

# import libraries
import math
import numpy as np
from score_utils import solve_ivodes
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt

# given and known constants
k0_1 = 2.59E9 # /min
E_1 = 16500. # cal /mol
dH_1 = -22200. # cal /mol
CA_0 = 2. # M
T_0 = 23  + 273.15 # K
Cp = 440. # cal /L /K
V = 4.0 # L
Ve = 0.5 # L
A = 0.6 # ft^2
U = 1.13E4/60 # cal /ft^2 /min /K
Te_in = 20 + 273.15 # K
rho_e = 1.0 # g /cm^3
Cp_e = 1.0 # cal /g /K
U_coil = 3.8E4/60 # cal /ft^2 /min /K
A_coil = 0.23 # ft^2
T_coil = 120 + 273.15 # K
Te_0 = 23 + 273.15 # K
T_1 = 50 + 273.15 # K
T_f = 25 + 273.15 # K
t_turn = 25 # min
Re = 1.987 # cal /mol /K

# define a variable for the current coolant flow rate
mDot_e = float('nan')

# derivatives function for first stage of operation
def first_stage_derivatives(ind, dep):
	# extract the dependent variables for this integration step
    nA = dep[0]
    nZ = dep[1]
    T = dep[2]
    Te = dep[3]

	# calculate the rate
    CA = nA/V
    k1 = k0_1*math.exp(-E_1/Re/T)
    r1 = k1*CA

    # calculate the rates of heat exchange
    Qdot_e = U*A*(Te - T)
    Qdot_coil = U_coil*A_coil*(T_coil - T)

	# evaluate the derivatives
    dnAdt = -V*r1
    dnZdt = V*r1
    dTdt = (Qdot_e + Qdot_coil - V*r1*dH_1)/V/Cp
    dTedt = -Qdot_e/rho_e/Ve/Cp_e
    
	# return the derivatives
    return [dnAdt, dnZdt, dTdt, dTedt]

# derivatives function for second stage of operation
def second_stage_derivatives(ind, dep):
	# extract the dependent variables for this integration step
    nA = dep[0]
    nZ = dep[1]
    T = dep[2]
    Te = dep[3]

	# calculate the rate
    CA = nA/V
    k1 = k0_1*math.exp(-E_1/Re/T)
    r1 = k1*CA

    # calculate the rates of heat exchange
    Qdot_e = U*A*(Te - T)

	# evaluate the derivatives
    dnAdt = -V*r1
    dnZdt = V*r1
    dTdt = (Qdot_e - V*r1*dH_1)/V/Cp
    dTedt = -(Qdot_e + mDot_e*Cp_e*(Te-Te_in))/rho_e/Ve/Cp_e
    
	# return the derivatives
    return [dnAdt, dnZdt, dTdt, dTedt]

# reactor model for first stage of operation
def first_stage_profiles():
	# set the initial values
    nA_0 = CA_0*V
    ind_0 = 0.0
    dep_0 = np.array([nA_0, 0.0, T_0, Te_0])

	# define the stopping criterion
    f_var = 3
    f_val = T_1
     
	# solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , first_stage_derivatives, True)

    # check that a solution was found
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    nA = dep[0,:]
    nZ = dep[1,:]
    T = dep[2,:]
    Te = dep[3,:]

    # return all profiles
    return t, nA, nZ, T, Te

# reactor model for second stage of operation
def second_stage_profiles(t0, nA0, nZ0, T0, Te0):
	# set the initial values
    ind_0 = t0
    dep_0 = np.array([nA0, nZ0, T0, Te0])

	# define the stopping criterion
    f_var = 3
    f_val = T_f
     
	# solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , second_stage_derivatives, True)

    # check that a solution was found
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    nA = dep[0,:]
    nZ = dep[1,:]
    T = dep[2,:]
    Te = dep[3,:]

    # return all profiles
    return t, nA, nZ, T, Te

# perform the analysis
def perform_the_analysis():
    # allow this function to set mDot_e
    global mDot_e

    # choose a range of coolant flow rates
    coolant_flows = np.linspace(100.0, 250.0, 100)

    # calculate the net rate for each coolant flow rate
    net_rates = np.zeros(100)
    for i in range(0,100):
        # make the coolant flow rate available
        mDot_e = coolant_flows[i]

        # solve the reactor design equations
        t1, nA1, nZ1, T1, Te1 = first_stage_profiles()
        t, nA, nZ, T, Te = second_stage_profiles(t1[-1], nA1[-1], nZ1[-1]
                                                , T1[-1], Te1[-1])
        
        # calculate the net rate
        net_rates[i] = nZ[-1]/(t[-1] + t_turn)

    # find the coolant flow where the net rate is maximized
    i_max = np.argmax(net_rates)
    mDot_max = coolant_flows[i_max]

    # solve the reactor design equations using the optimum coolant flow
    mDot_e = mDot_max
    t, nA, nZ, T, Te = first_stage_profiles()
    t2, nA2, nZ2, T2, Te2 = second_stage_profiles(t[-1], nA[-1], nZ[-1]
                                                , T[-1], Te[-1])
    
    # combine the profiles
    t = np.concatenate((t, t2))
    nA = np.concatenate((nA, nA2))
    nZ = np.concatenate((nZ, nZ2))
    T = np.concatenate((T, T2))
    Te = np.concatenate((Te, Te2))
    
    # calculate the quantities of interest
    pct_conversion = 100*(CA_0*V - nA)/(CA_0*V)
    
    # tabulate the results
    max_net_rate = nZ[-1]/(t[-1] + t_turn)
    data = [["Optimum Coolant Flow", f"{mDot_max}", "g min^-1^"],
            ["Maximum Net Rate", f"{max_net_rate}", "mol min^-1^"]]
    results_df = pd.DataFrame(data, columns=['item','value','units'])

    # display the results
    print(" ")
    print(f"Optimum Coolant Flow: {mDot_max} g /min")
    print(f"Maximum Net Rate: {max_net_rate} mol /min")

    # save the results
    results_df.to_csv('reb_9_6_4/python/results/reb_9_6_4_results.csv'
                      , index=False)

    # display and save the graphs
    plt.figure() # net rate vs coolant flow
    plt.plot(coolant_flows, net_rates)
    plt.xlabel("Coolant Flow (g/min)")
    plt.ylabel("Net Rate (mol/min)")
    plt.savefig('reb_9_6_4/python/results/reb_9_6_4_net_rate_vs_coolant_flow.png')
    plt.show()

    plt.figure() # conversion vs reaction time
    plt.plot(t,pct_conversion)
    plt.xlabel("Reaction Time (min)")
    plt.ylabel("Conversion (%)")
    plt.savefig('reb_9_6_4/python/results/reb_9_6_4_conversion_vs_time.png')
    plt.show()

    plt.figure() # temperature vs reaction time
    plt.plot(t,T - 273.15)
    plt.xlabel("Reaction Time (min)")
    plt.ylabel("Temperature (°C)")
    plt.savefig('reb_9_6_4/python/results/reb_9_6_4_temperature_vs_time.png')
    plt.show()

    # calculate and plot instantaneous rate profile for discussion
    r = k0_1*np.exp(-E_1/Re/T)*nA/V*1000.
    plt.figure() # instantantaneoud vs reaction time
    plt.plot(t,r)
    plt.xlabel("Reaction Time (min)")
    plt.ylabel("Instantaneous Rate (mol/L/min)")
    plt.savefig('reb_9_6_4/python/results/reb_9_6_4_inst_rate_vs_time.png')
    plt.show()

    return

if __name__=="__main__":
    perform_the_analysis()
