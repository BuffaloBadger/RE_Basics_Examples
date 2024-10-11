function reb_21_5_1()
%reb_21_5_1 Reaction Engineering Basics Example 21.5.1

    % given and known constants
    P = 1.0; % atm
    D = 1.0; % cm
    L = 10.0; % cm
    scc2mol = 1.0E-3/22.4; % mol/std cc

    % globally available variable
    gk = nan; % rate coefficient

    % derivatives function
    function ddz = derivatives(~,dep)
        % get dependent variables that are needed
        nA = dep(1);
        nY = dep(2);
        nZ = dep(3);

        % calculate the rate
        PA = nA/(nA + nY + nZ)*P;
        r = gk*PA;
    
        % evaluate the residuals
        dnAdz = -pi()*D^2/4*r;
        dnYdz = pi()*D^2/4*r;
        dnZdz = pi()*D^2/4*r;

        % return the derivatives
        ddz = [dnAdz; dnYdz; dnZdz];
    end

    % PFR model function
    function [z, nA, nY, nZ] = profiles(VA0, VY0, VZ0, beta)
        % make rate coefficient available to the derivatives function
        gk = 10^beta;

        % initial values and stopping criterion
        ind_0 = 0.0;
        nA0 = VA0*scc2mol;
        nY0 = VY0*scc2mol;
        nZ0 = VZ0*scc2mol;
        dep_0 = [nA0; nY0; nZ0];

        % stopping criterion
        f_var = 0;
        f_val = L;
        
        % solve the IVODEs
        odes_are_stiff = false;
        [z, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
            , @derivatives, odes_are_stiff);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: The ODE solution may not be accurate!')
        end

        % extract the dependent variable profiles
        nA = dep(:,1);
        nY = dep(:,2);
        nZ = dep(:,3);
    end

    % predicted responses function
    function fA_model = predicted_responses(params, expt_inputs)
        % get the number of experiments
        n_expt = size(expt_inputs,1);

        % allocate storage for the predicted responses
        fA_model = nan(n_expt,1);

        % loop through the experiments
        for i = 1:n_expt
            % extract the inputs for this experiment
            VA0 = expt_inputs(i,1);
            VY0 = expt_inputs(i,2);
            VZ0 = expt_inputs(i,3);
            % solve the PFR design equations
            [~, nA, ~, ~] = profiles(VA0, VY0, VZ0, params);

            % calculate the response
            fA_model(i) = 100*(VA0*scc2mol - nA(end))/(VA0*scc2mol);
        end
    end

    % quantities of interest function
    function [k, k_CI, r_squared, fA_model, epsilon_expt]...
            = quantities_of_interest(adj_inputs, fA)
        % guess the parameters
        par_guess = -2.0;
        %par_guess = [log10(5.34E5); 11800; log10(6.44E7); 18500];

        % estimate the parameters
        useRelErr = false;
        [beta, beta_ci, r_squared] = fit_to_SR_data(par_guess...
                , adj_inputs, fA, @predicted_responses...
                , useRelErr);

        % extract the results
        k = 10.^beta(1);
        k_CI = 10.^beta_ci(1,:);

        % calculate the model-predicted response and the residuals
        fA_model = predicted_responses(beta, adj_inputs);
        epsilon_expt = fA - fA_model;
    end

    % master function
    function perform_the_calculations()
        % read the experimental data
        data_table = readtable('../reb_21_5_1_data.csv'...
            , 'VariableNamingRule', 'preserve');
        data = table2array(data_table(:,:));

        % extract the adjusted inputs and response
        VA0 = data(:,1);
        VY0 = data(:,2);
        VZ0 = data(:,3);
        fA = data(:,4);
        adj_inputs = data(:,1:3);

        % calculate the quantities of interest
        [k, k_CI, r_squared, fA_model, epsilon_expt]...
            = quantities_of_interest(adj_inputs, fA);

        % tabulate, show, and save the results
        item = ["k"; "k_lower_limit"; "k_upper_limit"; "R_squared"];
        value = [k; k_CI(1); k_CI(2); r_squared];
        units = ["mol cm^-3^ min^-1^ atm^-1^"...
            ; "mol cm^-3^ min^-1^ atm^-1^"...
            ; "mol cm^-3^ min^-1^ atm^-1^"; ""];
        results_table = table(item, value, units);
        disp(results_table)
        writetable(results_table,'reb_21_5_1_results.csv')
        
        % create, show, and save a parity plot
        figure
        hold on
        plot(fA, fA_model,'ok','MarkerSize',10,'LineWidth',2)
        plot([min(fA),max(fA)], [min(fA), max(fA)],'r','LineWidth',2)
        hold off
        set(gca, 'FontSize', 14);
        xlabel('Experimental Conversion (%)','FontSize', 14)
        ylabel('Predicted Conversion (%)','FontSize', 14)
        legend({'Data','Parity Line'},'Location', 'northwest' ...
            ,'FontSize',14)
        saveas(gcf,'reb_21_5_1_parity.png')

        % create show, and save residuals plots
        figure
        plot(VA0, epsilon_expt,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('A Feed Rate (sccm)','FontSize', 14)
        ylabel('Residual (%)','FontSize', 14)
        saveas(gcf,'reb_21_5_1_A_residual.png')

        figure
        plot(VY0, epsilon_expt,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Y Feed Rate (sccm)','FontSize', 14)
        ylabel('Residual (%)','FontSize', 14)
        saveas(gcf,'reb_21_5_1_Y_residual.png')

        figure
        plot(VZ0, epsilon_expt,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Z Feed Rate (sccm)','FontSize', 14)
        ylabel('Residual (%)','FontSize', 14)
        saveas(gcf,'reb_21_5_1_Z_residual.png')
    end

    % perform the calculations
    perform_the_calculations()
end