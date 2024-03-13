"""Calculations for Reaction Engineering Basics Example 9.6.4"""

# import libraries

# given and known constants

# make missing initial value or IVODE constant available to all functions

# derivatives function
def derivatives(ind, dep):
	# extract the dependent variables for this integration step
	# calculate the rate
	# evaluate the derivatives
	# return the derivatives
	return list

# reactor model
def profiles():
	# set the initial values
    ind_0 = 0.0
    dep_0 = np.array([comma separated list])

	# define the stopping criterion
    f_var = 0
    f_val = t_f
     
	# solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , derivatives, True)

    # check that a solution was found
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    # return all profiles
    return list

# residuals function
def constant_residuals(list):
	return list

# perform the analysis
def perform_the_analysis():
    # solve the reactor design equations
    # calculate the other quantities of interest
    # tabulate the results
    # display the results
    # save the results
    # display and save the graphs
    return

if __name__=="__main__":
    perform_the_analysis()
