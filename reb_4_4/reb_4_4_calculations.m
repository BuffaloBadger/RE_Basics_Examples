function reb_4_4_calculations
    % given
    R = 8.314E-3; % kJ/mol/K

    % read the data file
    data_file = "./reb_4_4_data.csv";
    data_table = readtable(data_file,'VariableNamingRule','preserve');
    T_C = table2array(data_table(:,1)); % °C
    k = table2array(data_table(:,2)); % L/mol/min

    % convert temperatures to K
    T_K = T_C + 273.15;

    % function to fit to the data
    function yPred = calcY(kin_par, x)
        k0 = 10^kin_par(1);
        E = kin_par(2);
        yPred = k0*exp(-E/R./x);
    end

    % provide a guess for the parameters and fit the model to the data
    kin_par_guess = [
        0 % log(k0)
        20 % E in kJ/mol
    ];
    [beta, betaCI, R_squared] = fitNLSR(kin_par_guess, T_K, k, @calcY,...
        true);
    
    % extract the results
    k0 = 10^beta(1);
    k0_lower_limit = 10^betaCI(1,1);
    k0_upper_limit = 10^betaCI(1,2);
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
    results_file ="reb_4_4_Matlab_results.csv";
    item = ["k0";"k0_lower_limit";"k0_upper_limit"
        "E";"E_lower_limit";"E_upper_limit";"R_squared"];
    value = round([k0;k0_lower_limit;k0_upper_limit;E;E_lower_limit
        E_upper_limit;R_squared],3,'significant');
    units = ["L/mol/min";"L/mol/min";"L/mol/min";"kJ/mol";"kJ/mol"
        "kJ/mol";""];
    results_table = table(item,value,units);
    writetable(results_table,results_file);


    % plot the results
    T = linspace(min(T_K), max(T_K));
    yPred = calcY(beta, T);
    figure
    semilogy(T_K,k,'or',T,yPred,'k','LineWidth',2)
    set(gca, 'FontSize', 14);
    xlabel('T (K)','FontSize', 14)
    ylabel('Rate Coefficient (L mol^-^1 min^-^1)','FontSize', 14)
    legend({'Measured','Predicted'},'Location',...
        'northwest','FontSize',14)
    % save the graph
    saveas(gcf,"reb_4_4_Matlab_fig_1.png")
end