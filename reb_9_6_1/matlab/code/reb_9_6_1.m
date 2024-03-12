function reb_9_6_1
%REB_9_6_1 Calculations for Example 9.6.1 of Reaction Engineering Basics
    % given and known constants
    dH = -101.2E3; % J/mol
    k_0 = 5.11e4 * 3600; % L/mol/h
    E = 74.8e3; % J/mol
    T_0 = 180 + 273.15; % K
    V = 1900.0; % L
    CA_0 = 2.9; % mol/L
    CB_0 = 3.2; % mol/L
    Cp = 1.23 * 4.184; % J/g/K
    rho = 1.02 * 1000.0; % g/L
    t_f = 2.0; % h
    Re = 8.314; % J/mol/K

    % derivatives function
    function derivs = derivatives(~, dep)
        % Extract the dependent variables for this integration step
        nA = dep(1);
        nB = dep(2);
        T = dep(5);
        
        % Calculate the rate
        CA = nA/V;
        CB = nB/V;
        k = k_0*exp(-E/Re/T);
        r = k*CA*CB;

        % Evaluate the derivatives
        dAdt = -V*r;
        dBdt = -V*r;
        dYdt = V*r;
        dZdt = V*r;
        dTdt = -r*dH/rho/Cp;

        % Return the derivatives
        derivs = [dAdt; dBdt; dYdt; dZdt; dTdt];
    end

    % reactor model
    function [t, nA, nB, nY, nZ, T] = profiles()
        % set the initial values
        nA_0 = CA_0*V;
        nB_0 = CB_0*V;
        ind_0 = 0.0;
        dep_0 = [nA_0; nB_0; 0.0; 0.0; T_0];

        % define the stopping criterion
        f_var = 0;
        f_val = t_f;
        
        % solve the IVODEs
        odes_are_stiff = false;
        [t, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
            , @derivatives, odes_are_stiff);
    
        % Check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: The ODE solution may not be accurate!')
        end

        % extract and return the dependent variable profiles
        nA = dep(:,1);
        nB = dep(:,2);
        nY = dep(:,3);
        nZ = dep(:,4);
        T = dep(:,5);
    end

    % perform the analysis

    % solve the reactor design equations
    [t, nA, nB, nY, nZ, T] = profiles();

    % Calculate the other quantities of interest
    CA = nA/V;
    CB = nB/V;
    CY = nY/V;
    CZ = nZ/V;
    k = k_0*exp(-E/Re./T);
    r = k.*CA.*CB;
    T_C = T - 273.15;

    % Tabulate the results
    results_table = table(t,nA,nB,nY,nZ,T_C,r);

    % Display the results
    conc_plot = figure;
    hold("on")
    plot(t,CA,t,CB,t,CY,'LineWidth',2)
    plot(t,CZ,':','LineWidth',4)
    set(gca, 'FontSize', 14);
    xlabel('Time (h)','FontSize', 14)
    ylabel('Concentration (mol/L)','FontSize', 14)
    legend({'A','B','Y','Z'}, 'Location', 'east', 'FontSize', 14)

    temp_plot = figure;
    plot(t,T_C,'k','LineWidth',2)
    set(gca, 'FontSize', 14);
    xlabel('Time (h)','FontSize', 14)
    ylabel('Temperature (°C)','FontSize', 14)

    rate_plot = figure;
    plot(t,r,'k','LineWidth',2)
    set(gca, 'FontSize', 14);
    xlabel('Time (h)','FontSize', 14)
    ylabel('Rate (mol/L/h)','FontSize', 14)

    % Save the results
    saveas(conc_plot, '../results/reb_9_6_1_concentrations.png')
    saveas(temp_plot, '../results/reb_9_6_1_temperature.png')
    saveas(rate_plot, '../results/reb_9_6_1_rate.png')
    writetable(results_table, '../results/reb_9_6_1_results.csv');
end