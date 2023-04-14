function reb_10_3_calculations
    % given
    P = 1.0; % atm
    m = 3.0; % g
    VFR_in = 0.85; % L/min
    K = 12.2;

    % read the data file
    data_file = "../reb_10_2/reb_10_2_data.csv";
    data_table = readtable(data_file,'VariableNamingRule','preserve');

    % create adjusted input matrix
    adj_inp = table2array(data_table(:,1:4));
    nExpt = size(adj_inp,1);

    % create experimental response vector, eqn 5
    y_expt = table2array(data_table(:,5)); % atm

    % provide a guess for the parameters and fit the model to the data
    kin_par_guess = [
        5.09352415
        -4.30658933
        7.56602674
        7.3611074
        -2.85303979
    ];


    use_relative_errors = false;
    [beta, betaCI, R_squared] = fitNLSR(kin_par_guess, adj_inp, y_expt, ...
        @calcY, use_relative_errors);

    beta

    % extract the results
    k = 10^beta(1);
    k_lower_limit = 10^betaCI(1,1); % eqn 6
    k_upper_limit = 10^betaCI(1,2); % eqn 6
    KA = 10^beta(2);
    KA_lower_limit = 10^betaCI(2,1);
    KA_upper_limit = 10^betaCI(2,2);
    KB = 10^beta(3);
    KB_lower_limit = 10^betaCI(3,1);
    KB_upper_limit = 10^betaCI(3,2);
    KY = 10^beta(4);
    KY_lower_limit = 10^betaCI(4,1);
    KY_upper_limit = 10^betaCI(4,2);
    KZ = 10^beta(5);
    KZ_lower_limit = 10^betaCI(5,1);
    KZ_upper_limit = 10^betaCI(5,2);

    % report the results
    disp(' ')
    disp(['k = ',num2str(k,3),' [',num2str(k_lower_limit,3),...
        ', ',num2str(k_upper_limit,3),'] (95%CI) L/mol/min'])
    disp(['KA = ',num2str(KA,3),' [',num2str(KA_lower_limit,3),...
        ', ',num2str(KA_upper_limit,3),'] (95%CI) kJ/mol'])
    disp(['KB = ',num2str(KB,3),' [',num2str(KB_lower_limit,3),...
        ', ',num2str(KB_upper_limit,3),'] (95%CI) kJ/mol'])
    disp(['KY = ',num2str(KY,3),' [',num2str(KY_lower_limit,3),...
        ', ',num2str(KY_upper_limit,3),'] (95%CI) kJ/mol'])
    disp(['KZ = ',num2str(KZ,3),' [',num2str(KZ_lower_limit,3),...
        ', ',num2str(KZ_upper_limit,3),'] (95%CI) kJ/mol'])
    disp(['R_squared = ',num2str(R_squared,3)])

    % generate the parity plot
    y_pred = calcY(beta, adj_inp);
    parity_range = [0.9*min(y_expt), 1.1*max(y_expt)];
    figure
    loglog(y_expt,y_pred,'or',parity_range,parity_range,'k',...
        'MarkerFaceColor','r','MarkerSize',10,'LineWidth',2)
    set(gca, 'FontSize', 14);
    xlabel('Measured k (L mol^-^1 min^-^1)','FontSize', 14)
    ylabel('Predicted k (L mol^-^1 min^-^1)','FontSize', 14)
    legend({'Actual Model','Perfect Model'},'Location',...
        'northwest','FontSize',14)


    % response model function
    function ypred = calcY(kin_par, inputs)
        k_test = 10^kin_par(1);
        KA = 10^kin_par(2);
        KB = 10^kin_par(3);
        KY = 10^kin_par(4);
        KZ = 10^kin_par(5);

        ypred = nan(nExpt,1);


        for i = 1:nExpt
            m0 = 0.0;
            n0 = [inputs(i,1); inputs(i,2); inputs(i,3); inputs(i,4)]*...
                VFR_in/22.4;
            fVar = 0;
            fVal = m;
            [~, dep, flag] = solveIVODEs(m0, n0, fVar, fVal, @reactor_model, ...
                false, false, 1);
            if flag <= 0
                disp('The solution to the ODEs may not be accurate.')
            end
            n_A = dep(end,1);
            n_B = dep(end,2);
            n_Y = dep(end,3);
            n_Z = dep(end,4);
            n_tot = n_A + n_B + n_Y + n_Z;
            ypred(i) = n_A/n_tot*P;
        end
        function ddm = reactor_model(~,n)
            nA = n(1);
            nB = n(2);
            nY = n(3);
            nZ = n(4);
            ntot = nA + nB + nY + nZ;
            PA = nA/ntot*P;
            PB = nB/ntot*P;
            PY = nY/ntot*P;
            PZ = nZ/ntot*P;
            r = k_test*PA*PB/(1+KA*PA+KB*PB+KY*PY+KZ*PZ)*(1 - PY*PZ/K/PA/PB);
            ddm = [-r; -r; r; r];
        end
    end

%{
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
%}


%{
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
%}


%}
end