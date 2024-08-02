function reb_20_5_1()
%reb_20_5_1 Reaction Engineering Basics Example 20.5.1

    % given and known constants
    V = 0.1; % L

    % globally available variables
    Vdot = nan;
    CA_0 = nan;
    CB_0 = nan;
    CY_0 = nan;
    CZ_0 = nan;
    CY_1 = nan;
    i_expt_current = -1;
    k_current = nan;

    % residuals function
    function resids = residuals(guess)
        % get dependent variables that are needed
        nA_1 = guess(1);
        nB_1 = guess(2);
        nY_1 = guess(3);
        nZ_1 = guess(4);

        % calculate the other unknown quantities
        nA_0 = CA_0(i_expt_current)*Vdot(i_expt_current);
        nB_0 = CB_0(i_expt_current)*Vdot(i_expt_current);
        nY_0 = CY_0(i_expt_current)*Vdot(i_expt_current);
        nZ_0 = CZ_0(i_expt_current)*Vdot(i_expt_current);
        CA_1 = nA_1/Vdot(i_expt_current);
        CB_1 = nB_1/Vdot(i_expt_current);
        r = k_current*CA_1*CB_1;
    
        % evaluate the residuals
        epsilon_1 = nA_0 - nA_1 - V*r;
        epsilon_2 = nB_0 - nB_1 - V*r;
        epsilon_3 = nY_0 - nY_1 + V*r;
        epsilon_4 = nZ_0 - nZ_1 + V*r;

        % return the derivatives
        resids = [epsilon_1; epsilon_2; epsilon_3; epsilon_4];
    end

    % CSTR model function
    function [nA_1, nB_1, nY_1, nZ_1] = unknowns()
        % guess the solution
        nA_0 = CA_0(i_expt_current)*Vdot(i_expt_current);
        nB_0 = CB_0(i_expt_current)*Vdot(i_expt_current);
        nY_0 = CY_0(i_expt_current)*Vdot(i_expt_current);
        nZ_0 = CZ_0(i_expt_current)*Vdot(i_expt_current);
        nY_1 = CY_1(i_expt_current)*Vdot(i_expt_current);
        xi = nY_1 - nY_0;
        nA_1 = nA_0 - xi;
        nB_1 = nB_0 - xi;
        nZ_1 = nZ_0 + xi;
        initial_guess = [nA_1; nB_1; nY_1; nZ_1];
         
	    % solve the ATEs
        [soln, flag, message] = solve_ates(@residuals, initial_guess);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp(['WARNING: The ATE solver did not converge: ',message])
        end
    
        % extract the unknowns
        nA_1 = soln(1);
        nB_1 = soln(2);
        nY_1 = soln(3);
        nZ_1 = soln(4);
    end

    % predicted responses function
    function CY_1_model = predicted_responses(params, expt_inputs)
        % set the current value of k
        k_current = 10^params(1);

        % get the number of experiments
        n_expt = size(expt_inputs,1);

        % allocate storage for the predicted responses
        CY_1_model = nan(n_expt,1);

        % loop through the experiments
        for i = 1:n_expt
            % make the experiment index globally available
            i_expt_current = i;

            % solve the BSTR design equations
            [~, ~, nY_1,~] = unknowns();

            % calculate the response
            CY_1_model(i) = nY_1/Vdot(i);
        end
    end

    % function that performs the calculations
    function perform_the_calculations()
        % read the experimental data
        data_table = readtable('../reb_20_5_1_data.csv'...
            , 'VariableNamingRule', 'preserve');
        data = table2array(data_table(:,:));

        % extract the adjusted inputs and response
        Vdot = data(:,1)*1.0E-3;
        CA_0 = data(:,2);
        CB_0 = data(:,3);
        CY_0 = data(:,4);
        CZ_0 = data(:,5);
        CY_1 = data(:,6);
        adj_inputs = data(:,1:5);

        % guess the parameters
        par_guess = 0.0;

        % estimate the parameters
        useRelErr = false;
        [beta, betaCI, r_squared] = fit_to_SR_data(par_guess...
                , adj_inputs, CY_1, @predicted_responses, useRelErr);

        % extract the results
        k = 10^beta(1);
        k_ll = 10^betaCI(1,1);
        k_ul = 10^betaCI(1,2);
        
        % tabulate, show, and save the results
        item = ["k"; "k_lower_limit"; "k_upper_limit"; "R_squared"];
        value = [k; k_ll; k_ul; r_squared];
        units = ["L mol^-1^ min^-1^"; "L mol^-1^ min^-1^";
            "L mol^-1^ min^-1^"; ""];
        results_table = table(item, value, units);
        disp(results_table)
        writetable(results_table,'reb_20_5_1_results.csv')
        
        % calculate the model-predicted response and the residuals
        CY_1_model = predicted_responses(beta, adj_inputs);
        residual = CY_1 - CY_1_model;

        % create, show, and save a parity plot
        figure
        hold on
        plot([min(CY_1),max(CY_1)], [min(CY_1),max(CY_1)],'r'...
            ,'LineWidth',2)
        plot(CY_1, CY_1_model,'ok','MarkerSize',10,'LineWidth',2)
        hold off
        set(gca, 'FontSize', 14);
        xlabel('Experimental Outlet Y Concentration (M)','FontSize', 14)
        ylabel('Predicted Outlet Y Concentration (M)','FontSize', 14)
        saveas(gcf,'reb_20_5_1_parity.png')

        % create show, and save residuals plots
        figure
        plot(Vdot, residual,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Volumetric Flow Rate (cm^3 min^-^1)','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,'reb_20_5_1_Vdot_residual.png')

        figure
        plot(CA_0, residual,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet Concentration of A (M)','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,'reb_20_5_1_CA0_residual.png')

        figure
        plot(CB_0, residual,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet Concentration of B (M)','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,'reb_20_5_1_CB0_residual.png')

        figure
        plot(CY_0, residual,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet Concentration of Y (M)','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,'reb_20_5_1_CY0_residual.png')

        figure
        plot(CZ_0, residual,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet Concentration of Z (M)','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,'reb_20_5_1_CZ0_residual.png')

    end

    % perform the calculations
    perform_the_calculations()
end