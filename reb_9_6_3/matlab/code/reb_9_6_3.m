function reb_9_6_3
%REB_9_6_3 Calculations for Example 9.6.3 of Reaction Engineering Basics
    % given and known constants

    % make missing initial value or IVODE constant available to all functions
    ? = nan;

    % derivatives function
    function derivs = derivatives(~, dep)
        % extract the dependent variables for this integration step
        % calculate the rate
        % evaluate the derivatives
        % return the derivatives
        derivs = [column vector];
    end

    % reactor model
    function [row vector] = profiles()
        % set the initial values
        ind_0 = 0.0;
        dep_0 = [column vector];

        % define the stopping criterion
        f_var = ?;
        f_val = ?;
        
        % solve the IVODEs
        odes_are_stiff = false;
        [t, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
            , @derivatives, odes_are_stiff);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: The ODE solution may not be accurate!')
        end

        % extract and return the dependent variable profiles
    end

    % perform the analysis

    % solve the reactor design equations
    [row vector] = profiles();

    % calculate the other quantities of interest
    % tabulate the results
    results_table = table(column_vector_1,column_vector_2);
    %or
    item = ["item1";"item2"];
    value = [var1; var2];
    units = ["units1";"units2"];
    results_table = table(item,value,units);

    % display the results
    % save the results
    writetable(results_table,?file_spec?);

    % display and save the graphs
end