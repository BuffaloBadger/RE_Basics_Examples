function reb_23_3_1
% Calculations for Reaction Engineering Basics Example 23.3.1
    % make given and known constants available to all functions
    %V_nonIdeal = 10.0; % L
    P = 1.0; % atm
    VFR_feed = 3.0; % L/min
    T_feed = 300.0; % K
    CA_feed = 0.5; % mol/L
    k0 = 5.29E9; % L/mol/min
    E = 12100; % cal/mol
    dH = -24400; % cal/mol
    Cp = 1000.0; % cal/L/K
    rho = 1000.0; % g/L
    % basis
    delta_tB = 1.0; % min
    % universal constants
    Ren = 1.987; % cal/mol/K
    Rpv = 0.0821; % L atm /mol /K
    
    % read F vs. lambda and make it available to all functions
    data_table = readtable(...
        '../../reb_22_5_1/matlab/cum_age_dist_fcn.csv'...
        ,'VariableNamingRule','preserve');
    lambda = table2array(data_table(:,1)); % min
    F = table2array(data_table(:,2));

    % derivatives function
    function ddt = derivatives(t, dep)
        nA0 = dep(1);
        nZ0 = dep(2);
        T0 = dep(3);
        T1 = dep(6);
        V00 = VFR_feed*delta_tB;
        m00 = rho*V00;
        lambdale = lambda <= t;
        [~, i] = max(lambdale);
        if i == 1
            F_lambda = 0.0;
            dFdlambda = (F(2) - F(1))/(lambda(2) - lambda(1));
        else
            dFdlambda = (F(i) - F(i-1))/(lambda(i) - lambda(i-1));
            F_lambda = F(i-1) + dFdlambda*(t-lambda(i-1));
        end
        m0 = m00*(1-F_lambda);
        mDot = m00*dFdlambda;
        nAdot = nA0/m0*mDot;
        nZdot = nZ0/m0*mDot;
        V0 = m0/rho;
        V1 = V00 - V0;
        k = k0*exp(-E/Ren/T0);
        CA = nA0/V0;
        r = k*CA^2;
        dnA0dt = -nAdot - r*V0;
        dnZ0dt = -nZdot + r*V0;
        dT0dt = -(mDot*Ren/Rpv*P + rho*V0*r*dH)/rho/V0/Cp;
        dnA1dt = nAdot;
        dnZ1dt = nZdot;
        dT1dt = (mDot*Ren/Rpv*P + mDot*Cp*(T1-T0))/rho/V1/Cp;
        ddt = [dnA0dt; dnZ0dt; dT0dt; dnA1dt; dnZ1dt; dT1dt];
    end

    % early mixing segregated flow reactor function
    function emsfm()
        % initial values
        ind_0 = 0;
        dep_0 = [VFR_feed*CA_feed*delta_tB; 0; T_feed; 0.0; 0.0 ...
            ; T_feed];

        % stopping criterion
        f_var = 0;
        f_val = 25.0;

        % solve the design 
        odes_are_stiff = false;
        [~, dep, flag, message] = solve_ivodes(ind_0, dep_0, f_var...
            , f_val, @derivatives, odes_are_stiff);

        % check that the solver was successful
        if (flag <= 0.0)
            disp(" ")
            disp(["error solving the ivodes: ",message])
        end
        nADot_product = dep(end,4)/delta_tB;
        T_product = dep(end,6);
        fA_product = 100*(VFR_feed*CA_feed*delta_tB - nADot_product)...
            /(VFR_feed*CA_feed*delta_tB);

        disp(' ')
        disp(['T = ',num2str(T_product,3),' K'])
        disp(['fA = ',num2str(fA_product,3),'%'])
    end

    % quantities of interest function
    function quantities_of_interest()
        emsfm()
    end

    % calculate the quantities of interest
    quantities_of_interest();
end