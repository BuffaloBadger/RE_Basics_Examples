function reb_19_5_1_resp_fcn()
%reb_L_7_2 Reaction Engineering Basics Example 19.5.1

    % given and known constants
    V = 1.0; % 

    % make k available to all functions
    k_current = nan;

    % derivatives function
    function derivs = derivatives(~,dep)
        % extract nA
        nA = dep(1);

        % calculate the concentration of A
        CA = nA/V;

        % calculate the rate
        r = k_current*CA;
   
        % calculate the derivatives
        dnAdt = -r*V;
        dnZdt = r*V;

        % return the derivatives
        derivs = [dnAdt; dnZdt];
    end

    % BSTR model function
    function [t, nA, nZ] = profiles(CA0,t_f)
        % initial values
        ind_0 = 0.0;
        nA0 = CA0*V;
        dep_0 = [nA0; 0.0];

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
        nZ = dep(:,2);
    end

    % predicted responses function
    function CAf_model = predicted_responses(log_k_guess, expt_inputs)
        % set the current value of k
        k_current = 10^log_k_guess;

        % get the number of experiments
        n_expt = length(expt_inputs);

        % allocate storage for the predicted responses
        CAf_model = nan(n_expt,1);

        % loop through the experiments
        for i = 1:n_expt
            % solve the BSTR design equations
            [~, nA, ~] = profiles(expt_inputs(i,2),expt_inputs(i,3));

            % calculate the response
            nAf = nA(end);
            CAf_model(i) = nAf/V;
        end
    end

    % function that performs the calculations
    function perform_the_calculations()
        % read the data
        data_table = readtable('../python/reb_19_5_1_data.csv'...
            , 'VariableNamingRule', 'preserve');
        data = table2array(data_table(:,2:end));

        % get the block temperatures
        block_temps = unique(data(:,1));
        n_blocks = length(block_temps);

        % create arrays to store the fitting results
        T_fit = nan(n_blocks,1);
        k_fit = nan(n_blocks,1);
        k_ll_fit = nan(n_blocks,1);
        k_ul_fit = nan(n_blocks,1);
        r_sq_fit = nan(n_blocks,1);

        % process the data blocks
        for iBlock = 1:n_blocks
            % filter the data
            idx = data(:,1) == block_temps(iBlock);
            T = data(idx,1); % °C
            CA0 = data(idx,2); % M
            tf = data(idx,3); % min
            CAf = data(idx,4); % M

            % combine the adjusted inputs into a matrix
            adj_inputs = [T, CA0, tf];

            % guess beta, the base 10 log of k
            if iBlock == 1
                guess = 0.0;
            else
                guess = log10(k_fit(iBlock-1));
            end

            % estimate log_10 of k
            useRelErr = false;
            [beta, betaCI, r_squared] = fit_to_SR_data(guess...
                , adj_inputs, CAf, @predicted_responses, useRelErr);

            % extract the results
            T_fit(iBlock) = T_k - 273.15;
            k_fit(iBlock) = 10^beta;
            k_ll_fit(iBlock) = 10^betaCI(1);
            k_ul_fit(iBlock) = 10^betaCI(2);
            r_sq_fit(iBlock) = r_squared;

            % calculate the model-predicted responses and residuals
            CAf_model = predicted_responses(beta, adj_inputs);
            residual = CAf - CAf_model;

            % create, show, and save a parity plot
            figure
            hold on
            plot([min(CAf),max(CAf)], [min(CAf),max(CAf)]...
                ,'r','LineWidth',2)
            plot(CAf, CAf_model,'ok','MarkerSize',10,'LineWidth',2)
            hold off
            set(gca, 'FontSize', 14);
            xlabel('Experimental C_A_f (M)','FontSize', 14)
            ylabel('Predicted C_A_f (M)','FontSize', 14)
            filename = strcat('reb_19_5_1_resp_fcn_parity_'...
                ,int2str(T(1)),'.png');
            saveas(gcf,filename)

            % generate, show, and save residuals plots
            figure
            plot(tf, residual,'ok','MarkerSize',10,'LineWidth',2)
            yline(0.0,'r','LineWidth',2)
            set(gca, 'FontSize', 14);
            xlabel('Elapsed Time (min)','FontSize', 14)
            ylabel('Residual (M)','FontSize', 14)
            filename = strcat('reb_19_5_1_resp_fcn_residual_tf_'...
                ,int2str(T(1)),'.png');
            saveas(gcf,filename)

            figure
            plot(CA0, residual,'ok','MarkerSize',10,'LineWidth',2)
            yline(0.0,'r','LineWidth',2)
            set(gca, 'FontSize', 14);
            xlabel('Elapsed Time (min)','FontSize', 14)
            ylabel('Residual (M)','FontSize', 14)
            filename = strcat('reb_19_5_1_resp_fcn_residual_CA0_'...
                ,int2str(T(1)),'.png');
            saveas(gcf,filename)
        end

        % tabulate, show, and save the results
        results_table = table(T_fit, k_fit, k_ll_fit, k_ul_fit...
            , r_sq_fit);
        results_table.Properties.VariableNames = ["T", "k", "k_ll"...
            ,"k_ul", "R_sq"];
        results_file ="reb_19_5_1_resp_fcn_params.csv";        
        disp(results_table)
        writetable(results_table,results_file);
    end

    % perform the calculations
    perform_the_calculations()
end