function reb_9_6_1_matlab
%REB_9_6_1_MATLAB calculations for Example 9.6.1 in Reaction Engineering
%   Basics using Matlab

    % Given and known constants
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

    % Derivatives function
    function derivs = eval_derivs(~, cur_dep)
        % Extract the dependent variables for the current 
        %   integration step
        nA_cur = cur_dep(1);
        nB_cur = cur_dep(2);
        T_cur = cur_dep(5);
        
        % Calculate the rate
        CA = nA_cur/V;
        CB = nB_cur/V;
        k = k_0*exp(-E/Re/T_cur);
        r_cur = k*CA*CB;

        % Evaluate the derivatives
        dAdt = -V*r_cur;
        dBdt = -V*r_cur;
        dYdt = V*r_cur;
        dZdt = V*r_cur;
        dTdt = -r_cur*dH/rho/Cp;

        % Return the derivatives
        derivs = [dAdt; dBdt; dYdt; dZdt; dTdt];
    end
    
    % Inlet molar flow rates
    nA_0 = CA_0*V;
    nB_0 = CB_0*V;

    % Initial values
    ind_0 = 0.0;
    dep_0 = [nA_0; nB_0; 0.0; 0.0; T_0];

    % Stopping criterion
    f_var = 0;
    f_val = t_f;

    % Solve the IVODEs
    odes_are_stiff = false;
    [t, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
        , @eval_derivs, odes_are_stiff);

    % Check that the solution was found
    if flag <= 0
        disp(' ')
        disp('WARNING: The ODE solution may not be accurate!')
    end

    % Extract the results
    nA = dep(:,1);
    nB = dep(:,2);
    nY = dep(:,3);
    nZ = dep(:,4);
    T = dep(:,5);

    % Calculate the other quantities of interest
    CA = nA/V;
    CB = nB/V;
    CY = nY/V;
    CZ = nZ/V;
    k = k_0*exp(-E/Re./T);
    r = k.*CA.*CB;

    % Plot the results as requested
    figure;
    hold("on")
    plot(t,CA,t,CB,t,CY,'LineWidth',2)
    plot(t,CZ,':','LineWidth',4)
    set(gca, 'FontSize', 14);
    xlabel('Time (h)','FontSize', 14)
    ylabel('Concentration (mol/L)','FontSize', 14)
    legend({'A','B','Y','Z'}, 'Location', 'east', 'FontSize', 14)

    figure;
    plot(t,T-273.15,'k','LineWidth',2)
    set(gca, 'FontSize', 14);
    xlabel('Time (h)','FontSize', 14)
    ylabel('Temperature (°C)','FontSize', 14)

    figure;
    plot(t,r,'k','LineWidth',2)
    set(gca, 'FontSize', 14);
    xlabel('Time (h)','FontSize', 14)
    ylabel('Rate (mol/L/h)','FontSize', 14)
end