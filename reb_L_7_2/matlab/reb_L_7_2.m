function reb_L_7_2()
%reb_L_7_2 Reaction Engineering Basics Example L.7.2

    % given and known constants
    CA_0 = 1.0; % M
    CB_0 = 0.9; % M
    V = 1.0; % L (basis)

    % make k available to all functions
    k_current = nan;

    % derivatives function
    function derivs = derivatives(~,dep)
        CA = dep(1)/V;
        CB = dep(2)/V;
        dnAdt = -k_current*V*CA*CB;
        dnBdt = -k_current*V*CA*CB;
        dnYdt = k_current*V*CA*CB;
        dnZdt = k_current*V*CA*CB;
        derivs = [dnAdt; dnBdt; dnYdt; dnZdt];
    end

    % BSTR model function
    function [t, nA, nB, nY, nZ] = profiles(t_f)
        % initial values
        ind_0 = 0.0;
        dep_0 = [CA_0*V; CB_0*V; 0.0; 0.0];

        % stopping criterion
        f_var = 0;
        f_val = t_f;
        
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
        nA = dep(:,1);
        nB = dep(:,2);
        nY = dep(:,3);
        nZ = dep(:,4);
    end

    % predicted responses function
    function CY_model = predicted_responses(log_k_guess, expt_inputs)
        % set the current value of k
        k_current = 10^log_k_guess;

        % get the number of experiments
        n_expt = length(expt_inputs);

        % allocate storage for the predicted responses
        CY_model = nan(n_expt,1);

        % loop through the experiments
        for i = 1:n_expt
            % solve the BSTR design equations
            [~, ~, ~, nY, ~] = profiles(expt_inputs(i));

            % calculate the response
            CY_model(i) = nY(end)/V;
        end
    end

    % function that performs the calculations
    function perform_the_calculations()
        % read the experimental data
        data_file = '../reb_L_7_2_data.csv';
        data_table = readtable(data_file, 'VariableNamingRule'...
            , 'preserve');
        t_f = table2array(data_table(:,1)); % min
        CY_expt = table2array(data_table(:,2)); % M

        % guess beta, the base 10 log of k
        guess = 0.0;

        % fit the predicted responses model to the data
        useRelErr = false;
        [beta, betaCI, r_squared] = fit_to_SR_data(guess, t_f...
            , CY_expt, @predicted_responses, useRelErr);
        
        % extract the results
        k_current = 10^beta;
        k_CI_lower = 10^betaCI(1);
        k_CI_upper = 10^betaCI(2);

        % generate, show, and save a parity plot
        CY_model = predicted_responses(beta, t_f);
        figure
        hold on
        plot([min(CY_expt),max(CY_expt)], [min(CY_expt),max(CY_expt)]...
            ,'r','LineWidth',2)
        plot(CY_expt, CY_model,'ok','MarkerSize',10,'LineWidth',2)
        hold off
        set(gca, 'FontSize', 14);
        xlabel('Experimental C_Y (M)','FontSize', 14)
        ylabel('Predicted C_Y (M)','FontSize', 14)
        saveas(gcf,"parity_plot.png")

        % generate, show, and save a resituals plot
        resid = CY_expt - CY_model;
        figure
        plot(t_f, resid,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Elapsed Time (min)','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,"residuals_vs_tf_plot.png")

        % tabulate, show, and save the results
        results_file ="results.csv";
        item = ["k";"k_CI_lower"; "k_CI_upper";"R squared"];
        value = [k_current; k_CI_lower; k_CI_upper; r_squared];
        units = ["L mol^-1^ min^-1^"; "L mol^-1^ min^-1^"; 
            "L mol^-1^ min^-1^"; ""];
        results_table = table(item, value, units);
        disp(results_table)
        writetable(results_table,results_file);
    end

    % perform the calculations
    perform_the_calculations()
end