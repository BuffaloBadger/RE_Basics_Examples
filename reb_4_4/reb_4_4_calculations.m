function reb_4_4_calculations
    % given
    R = 8.314E-3; % kJ/mol/K

    % read the data file
    data_file = "./reb_4_4_data.csv";
    data_table = readtable(data_file,'VariableNamingRule','preserve');
    T = table2array(data_table(:,2)); % °C
    k = table2array(data_table(:,3)); % L/mol/min

    % create adjusted input vector, eqn 4
    x = T + 273.15;

    % create experimental response vector, eqn 5
    y_expt = k;

    % response model function, eqn 3
    function yPred = calcY(kin_par, T_K)
        k0 = 10^kin_par(1);
        E = kin_par(2);
        yPred = k0*exp(-E/R./T_K);
    end

    % provide a guess for the parameters and fit the model to the data
    kin_par_guess = [
        0 % log(k0)
        20 % E in kJ/mol
    ];
    use_relative_errors = false;
    [beta, betaCI, R_squared] = fitNLSR(kin_par_guess, x, y_expt, ...
        @calcY, use_relative_errors);
    
    % extract the results
    k0 = 10^beta(1); % eqn 6
    k0_lower_limit = 10^betaCI(1,1); % eqn 6
    k0_upper_limit = 10^betaCI(1,2); % eqn 6
    E = beta(2);
    E_lower_limit = betaCI(2,1);
    E_upper_limit = betaCI(2,2);

    % report the results
    disp(' ')
    disp(['k0 = ',num2str(k0,3),' [',num2str(k0_lower_limit,3),...
        ', ',num2str(k0_upper_limit,3),'] (95%CI) L/mol/min'])
    disp(['E = ',num2str(E,3),' [',num2str(E_lower_limit,3),...
        ', ',num2str(E_upper_limit,3),'] (95%CI) kJ/mol'])
    disp(['R_squared = ',num2str(R_squared,3)])

    % save the results to a .csv file
    if(use_relative_errors)
        results_file ="reb_4_4_Matlab_results_rel.csv";
    else
        results_file ="reb_4_4_Matlab_results_abs.csv";
    end
    item = ["k0";"k0_lower_limit";"k0_upper_limit"
        "E";"E_lower_limit";"E_upper_limit";"R_squared"];
    value = round([k0;k0_lower_limit;k0_upper_limit;E;E_lower_limit
        E_upper_limit;R_squared],3,'significant');
    units = ["L/mol/min";"L/mol/min";"L/mol/min";"kJ/mol";"kJ/mol"
        "kJ/mol";""];
    results_table = table(item,value,units);
    writetable(results_table,results_file);


    % generate the parity plot
    y_pred = calcY(beta, x); % eqn 7
    parity_range = [0.9*min(y_expt), 1.1*max(y_expt)];
    figure
    loglog(y_expt,y_pred,'or',parity_range,parity_range,'k',...
        'MarkerFaceColor','r','MarkerSize',10,'LineWidth',2)
    set(gca, 'FontSize', 14);
    xlabel('Measured k (L mol^-^1 min^-^1)','FontSize', 14)
    ylabel('Predicted k (L mol^-^1 min^-^1)','FontSize', 14)
    legend({'Actual Model','Perfect Model'},'Location',...
        'northwest','FontSize',14)
    % save the graph
    if(use_relative_errors)
        saveas(gcf,"reb_4_4_Matlab_parity_plot_rel.png")
    else
        saveas(gcf,"reb_4_4_Matlab_parity_plot_abs.png")
    end

    % generate the residuals plot
    residual = y_expt - y_pred; % eqn 8
    figure
    plot(x,residual,'or','MarkerFaceColor','r','MarkerSize',10)
    set(gca, 'FontSize', 14);
    yline(0,'LineWidth',2)
    xlabel('T (K)','FontSize', 14)
    ylabel('Residual (L mol^-^1 min^-^1)','FontSize', 14)
    if(use_relative_errors)
        saveas(gcf,"reb_4_4_Matlab_residual_plot_rel.png")
    else
        saveas(gcf,"reb_4_4_Matlab_residual_plot_abs.png")
    end
end