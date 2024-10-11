function reb_21_5_3()
%reb_21_5_3 Reaction Engineering Basics Example 21.5.3

    % given and known constants
    P = 1.0; % atm
    mTotal = 3.0; % g
    VFR_in = 0.85; % L/min
    K = 12.2;
    sL2mol = 22.4; % mol/std cc

    % globally available variables
    gk = nan; % rate coefficient
    gKa = nan; % alpha A
    gKb = nan; % alpha B
    gKy = nan; % alpha Y
    gKz = nan; % alpha Z

    % derivatives function
    function ddm = derivatives(~,dep)
        % get dependent variables that are needed
        nA = dep(1);
        nB = dep(2);
        nY = dep(3);
        nZ = dep(4);

        % calculate the rate
        ntot = nA + nB + nY + nZ;
        PA = nA/ntot*P;
        PB = nB/ntot*P;
        PY = nY/ntot*P;
        PZ = nZ/ntot*P;
        r = gk*PA*PB/(1 + gKa*PA + gKb*PB + gKy*PY + gKz*PZ)...
            *(1 - PY*PZ/K/PA/PB);

        dnAdm = -r;
        dnBdm = -r;
        dnYdm = r;
        dnZdm = r;

        % return the derivatives
        ddm = [dnAdm; dnBdm; dnYdm; dnZdm];
    end

    % PFR model function
    function [m, nA, nB, nY, nZ] = profiles(yA0, yB0, yY0, yZ0, k...
            , KA, KB, KY, KZ)
        % make rate coefficient available to the derivatives function
        gk = k;
        gKa = KA;
        gKb = KB;
        gKy = KY;
        gKz = KZ;

        % initial values and stopping criterion
        ind_0 = 0.0;
        nA0 = yA0*VFR_in/sL2mol;
        nB0 = yB0*VFR_in/sL2mol;
        nY0 = yY0*VFR_in/sL2mol;
        nZ0 = yZ0*VFR_in/sL2mol;
        dep_0 =[nA0, nB0, nY0, nZ0];

        % stopping criterion
        f_var = 0;
        f_val = mTotal;
        
        % solve the IVODEs
        odes_are_stiff = false;
        [m, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
            , @derivatives, odes_are_stiff);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: The ODE solution may not be accurate!')
        end

        % extract the dependent variable profiles
        nA = dep(:,1);
        nB = dep(:,2);
        nY = dep(:,3);
        nZ = dep(:,4);
    end

    % predicted responses function
    function PA_1_model = predicted_responses(params, expt_inputs)
        % get the number of experiments
        n_expt = size(expt_inputs,1);

        % allocate storage for the predicted responses
        PA_1_model = nan(n_expt,1);

        % extract the rate expression parameters
        k = 10^params(1);
        KA = 10.^params(2);
        KB = 10.^params(3);
        KY = 10.^params(4);
        KZ = 10.^params(5);

        % loop through the experiments
        for i = 1:n_expt
            % extract the inputs for this experiment
            yA0 = expt_inputs(i,1);
            yB0 = expt_inputs(i,2);
            yY0 = expt_inputs(i,3);
            yZ0 = expt_inputs(i,4);
            % solve the PFR design equations
            [~, nA, nB, nY, nZ] = profiles(yA0, yB0, yY0, yZ0, k...
                , KA, KB, KY, KZ);

            % calculate the response
            nA_1 = nA(end);
            nB_1 = nB(end);
            nY_1 = nY(end);
            nZ_1 = nZ(end);
            PA_1_model(i) = nA_1/(nA_1 + nB_1 + nY_1 + nZ_1)*P;
        end
    end

    % quantities of interest function
    function [k, k_CI, KA, KA_CI, KB, KB_CI, KY...
            , KY_CI, KZ, KZ_CI, r_squared, PA_1_model...
            , epsilon_expt] = quantities_of_interest(adj_inputs...
            , PA_1, par_guess)

        % estimate the parameters
        useRelErr = false;
        [beta, beta_ci, r_squared] = fit_to_SR_data(par_guess...
                , adj_inputs, PA_1, @predicted_responses...
                , useRelErr);

        % extract the results
        k = 10.^beta(1);
        k_CI = 10.^beta_ci(1,:);
        KA = 10.^beta(2);
        KA_CI = 10.^beta_ci(2,:);
        KB = 10.^beta(3);
        KB_CI = 10.^beta_ci(3,:);
        KY = 10.^beta(4);
        KY_CI = 10.^beta_ci(4,:);
        KZ = 10.^beta(5);
        KZ_CI = 10.^beta_ci(5,:);

        % calculate the model-predicted response and the residuals
        PA_1_model = predicted_responses(beta, adj_inputs);
        epsilon_expt = PA_1 - PA_1_model;
    end

    % master function
    function perform_the_calculations()
        % read the experimental data
        data_table = readtable('../../reb_21_5_2/reb_21_5_2_data.csv'...
            , 'VariableNamingRule', 'preserve');
        data = table2array(data_table(:,:));

        % extract the adjusted inputs and response
        yA0 = data(:,1);
        yB0 = data(:,2);
        yY0 = data(:,3);
        yZ0 = data(:,4);
        PA_1 = data(:,5);
        adj_inputs = data(:,1:4);

        % guess the parameters
        %par_guess = [1.0; 1.0; 1.0; 1.0; 1.0];
        par_guess = [-3.0; -3.0; 0.1; 0.4; -5.0];

        % calculate the quantities of interest
        [k, k_CI,KA, KA_CI, KB, KB_CI, KY...
            , KY_CI, KZ, KZ_CI, r_squared, PA_1_model...
            , epsilon_expt] = quantities_of_interest(adj_inputs...
            , PA_1, par_guess);

        % tabulate, show, and save the results
        kUnits = "mol g^-1^ min^-1 atm^-2";
        Kunits = "atm^-1^";
        item = ["k"; "k_lower_limit"; "k_upper_limit"...
            ;"KA"; "KA_lower_limit"; "KA_upper_limit"...
            ;"KB"; "KB_lower_limit"; "KB_upper_limit"...
            ;"KY"; "KY_lower_limit"; "KY_upper_limit"...
            ;"KZ"; "KZ_lower_limit"; "KZ_upper_limit"...
            ; "R_squared"];
        value = [k; k_CI(1); k_CI(2); KA; KA_CI(1)...
            ; KA_CI(2); KB; KB_CI(1); KB_CI(2)...
            ; KY; KY_CI(1); KY_CI(2); KZ...
            ; KZ_CI(1); KZ_CI(2); r_squared];
        units = [kUnits; kUnits; kUnits;Kunits;Kunits;Kunits;Kunits;Kunits...
            ;Kunits;Kunits;Kunits;Kunits;Kunits;Kunits;Kunits;""];
        results_table = table(item, value, units);
        disp(results_table)
        writetable(results_table,'reb_21_5_3_results.csv')
        
        % create, show, and save a parity plot
        figure
        hold on
        plot(PA_1, PA_1_model,'ok','MarkerSize',10,'LineWidth',2)
        plot([min(PA_1),max(PA_1)], [min(PA_1), max(PA_1)],'r'...
            ,'LineWidth',2)
        hold off
        set(gca, 'FontSize', 14);
        xlabel('Experimental P_A (atm)','FontSize', 14)
        ylabel('Predicted P_A (atm)','FontSize', 14)
        legend({'Data','Parity Line'},'Location', 'northwest' ...
            ,'FontSize',14)
        saveas(gcf,'reb_21_5_3_parity.png')

        % create show, and save residuals plots
        figure
        plot(yA0, epsilon_expt,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet A Mole Fraction','FontSize', 14)
        ylabel('Residual (atm)','FontSize', 14)
        saveas(gcf,'reb_21_5_3_A_residual.png')

        figure
        plot(yB0, epsilon_expt,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet B Mole Fraction','FontSize', 14)
        ylabel('Residual (atm)','FontSize', 14)
        saveas(gcf,'reb_21_5_3_B_residual.png')

        figure
        plot(yY0, epsilon_expt,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet Y Mole Fraction','FontSize', 14)
        ylabel('Residual (atm)','FontSize', 14)
        saveas(gcf,'reb_21_5_3_Y_residual.png')

        figure
        plot(yZ0, epsilon_expt,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet Z Mole Fraction','FontSize', 14)
        ylabel('Residual (atm)','FontSize', 14)
        saveas(gcf,'reb_21_5_2_3_residual.png')
    end

    % perform the calculations
    perform_the_calculations()
end