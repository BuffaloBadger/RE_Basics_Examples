function reb_9_6_2
%REB_9_6_2 Calculations for Example 9.6.2 of Reaction Engineering Basics
    % given and known constants
    V = 10.0E3; % cm^3
    Ve = 1.4E3; % cm^3
    U = 138.; % cal /ft^2 /min /K
    A = 1200./929.; % ft^2
    Te_in = 40. + 273.15; % K
    mDot_e = 100.; % g /min
    rho = 1.0; % g /cm^3
    rho_e = 1.0; % g /cm^3
    Cp = 1.0; % cal /g /K
    Cp_e = 1.0; % cal /g /K
    CA_0 = 5.0E-3; % mol /cm^3
    CB_0 = 7.0E-3; % mol /cm^3
    dH_1 = -16.7E3; % cal /mol
    dH_2 = -14.3E3; % cal /mol
    k0_1 = 9.74E12; % cm^3 /mol /min
    E_1 = 20.1E3; % cal /mol
    k0_2 = 2.38E13; % /min
    E_2 = 25.3E3; % cal /mol
    t_rxn = 30.; % min
    fA_f = 0.45;
    Re = 1.987; % cal /mol /K

    % make missing initial value or IVODE constant available to all functions
    T0 = nan;

    % derivatives function
    function derivs = derivatives(~, dep)
        % extract the dependent variables for this integration step
        nA = dep(1);
        nB = dep(2);
        nX = dep(3);
        nY = dep(4);
        nZ = dep(5);
        T = dep(6);
        Te = dep(7);

        % calculate the rate
        CA = nA/V;
        CB = nB/V;
        k_1 = k0_1*exp(-E_1/Re/T);
        r1 = k_1*CA*CB;
        k_2 = k0_2*exp(-E_2/Re/T);
        r2 = k_2*CA;

        % calculate the rate of heat exchange
        Qdot = U*A*(Te-T);

        % evaluate the derivatives
        dnAdt = -V*(r1 + r2);
        dnBdt = -V*r1;
        dnXdt = V*r1;
        dnYdt = V*r1;
        dnZdt = V*r2;
        dTdt = (Qdot - V*(r1*dH_1 + r2*dH_2))/rho/V/Cp;
        dTedt = (-Qdot - mDot_e*Cp_e*(Te-Te_in))/rho_e/Ve/Cp_e;

        % return the derivatives
        derivs = [dnAdt; dnBdt; dnXdt; dnYdt; dnZdt; dTdt; dTedt];
    end

    % reactor model
    function [t, nA, nB, nX, nY, nZ, T, Te] = profiles()
        % set the initial values
        nA_0 = CA_0*V;
        nB_0 = CB_0*V;
        ind_0 = 0.0;
        dep_0 = [nA_0; nB_0; 0.0; 0.0; 0.0; T0; Te_in];

        % define the stopping criterion
        f_var = 0;
        f_val = t_rxn;
        
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
        nX = dep(:,3);
        nY = dep(:,4);
        nZ = dep(:,5);
        T = dep(:,6);
        Te = dep(:,7);
    end

    % implicit equation for IVODE initial value as residual
    function resid = residual(guess)
        % make the guess available to all functions
        T0 = guess;

        % solve the reactor design equations
        [t, nA, nB, nX, nY, nZ, T, Te] = profiles();

        % extract the calculated final molar amount of A
        nA_f = nA(end);

        % evaluate and return the residual
        resid = nA_f - CA_0*V*(1-fA_f);
    end

    % perform the analysis

    % initial guess for T0
    initial_guess = Te_in + 20;

    % calculate T0
    [T0, flag, message] = solve_ates(@residual, initial_guess);

    % check that the solution converged
    if flag <= 0
        disp(' ')
        disp(['The ATE solver did not converge: ',message])
    end

    % solve the reactor design equations
    [t, nA, nB, nX, nY, nZ, T, Te] = profiles();

    % calculate the other quantities of interest
    T0 = T0 - 273.15;
    T_f = T(end) - 273.15;
    Te_out = Te(end) - 273.15;
    sel_X_Z = nX(end)/nZ(end);

    % tabulate the results
    item = ["T0";"T_f";"Te_f";"sel_X_Z"];
    value = [T0; T_f; Te_out; sel_X_Z];
    units = ["°C"; "°C"; "°C"; " "];
    results_table = table(item,value,units);

    % display the results
    disp(' ')
    disp(['Initial Temperature: ',num2str(T0,3), ' °C'])
    disp(['Final Temperature: ',num2str(T_f,3), ' °C'])
    disp(['Outlet Coolant Temperature: ',num2str(Te_out,3), ' °C'])
    disp(['Selectivity: ',num2str(sel_X_Z,3), ' mol X per mol Z'])
    disp(' ')

    % save the results
    writetable(results_table,'../results/reb_9_6_2_results.csv');
end