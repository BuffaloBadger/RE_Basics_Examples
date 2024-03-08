function reb_J_7_3
%REB_J_7_3 Calculations for Example J.7.3 of Reaction Engineering Basics
    % given and known constants
    T_out = 400.; % K
    D = 1. ; % in
    k_0 = 7.49E9*61.02; % in^3 /mol /min
    E = 15300.; % cal /mol
    P = 4.; % atm
    dH = -14500.; % cal /mol
    Cp_A = 10.9; % cal /mol /K
    Cp_Z = 21.8; % cal /mol /K
    L = 100.; % in
    nDot_A_in = 1.5; % mol /min
    Re = 1.987; % cal /mol /K
    Rw = 0.08206*61.02; % in^3 atm /mol /K

    % reactor design equations as derivative expressions
    function derivs = derivatives(~, dep)
        % extract ind and dep vars for the current integration step
        nDot_A = dep(1);
        nDot_Z = dep(2);
        T = dep(3);

        % calculate rate
        [r] = other_ivode_variables(T, nDot_A, nDot_Z);

        % evaluate the derivatives
        dnDotAdz = -2*pi()*D^2/4*r;
        dnDotZdz = pi()*D^2/4*r;
        dTdz = -pi()*D^2/4*r*dH/(nDot_A*Cp_A + nDot_Z*Cp_Z);

        % return the derivatives
        derivs = [dnDotAdz; dnDotZdz; dTdz];

    end

    % calculate other IVODE variables
    function [r] = other_ivode_variables(T, nDot_A, nDot_Z)
        r = k_0*exp(-E/Re/T)*(nDot_A*P/Rw/T/(nDot_A + nDot_Z))^2;
    end

    % calculate IVODE initial and final values
    function [ind_0, dep_0, f_var, f_val] ...
            = initial_and_final_values(T_in_guess)
        % initial values
        ind_0 = 0.0;
        dep_0 = [nDot_A_in; 0.0; T_in_guess];

        % stopping criterion
        f_var = 0;
        f_val = L;
    end

    % implicit equation for IVODE initial value as residual
    function resid = residual(guess)
        % solve the reactor design equations using the guess
        T_in_guess = guess;
        [~, ~, ~, T] = profiles(T_in_guess);

        % extract the calculated final temperature
        T_f = T(end);

        % evaluate and return the residual
        resid = T_f - T_out;
    end

    % solve the reactor design equations
    function [z, nDot_A, nDot_Z, T] = profiles(T_in)
        % get the initial values and stopping criterion
        [ind_0, dep_0, f_var, f_val] ...
            = initial_and_final_values(T_in);

        % solve the IVODEs
        odes_are_stiff = false;
        [z, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
            , @derivatives, odes_are_stiff);
    
        % Check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: The ODE solution may not be accurate!')
        end

        % extract the dependent variable profiles
        nDot_A = dep(:,1);
        nDot_Z = dep(:,2);
        T = dep(:,3);
        
    end

    % complete the assignment

    % calculate the inlet temperature
    % initial guess
    initial_guess = T_out - 100.0;

    % solve the implicit equation for Tin
    [T_in, flag, message] = solve_ates(@residual, initial_guess);

    % check that the solution converged
    if flag <= 0
        disp(' ')
        disp(['The ATE solver did not converge: ',message])
    end

    % get the solution of the reactor design equations
    [z, nDot_A, nDot_Z, T] = profiles(T_in);

    % tabulate the results
    item = "T_in";
    value = T_in;
    units = "K";
    Tin_results_table = table(item,value,units);
    profile_results_table = table(z, nDot_A, nDot_Z, T);

    % display the results
    disp(' ')
    disp(['Inlet Temperature: ', num2str(T_in,3), ' K'])
    disp(' ')
    disp('Molar flow and Temperature Profiles')
    disp(profile_results_table)

    % Save the results
    Tin_results_file = "../results/reb_J_7_3_Tin_results.csv";
    writetable(Tin_results_table,Tin_results_file);
    profile_results_file ...
        = "../results/reb_J_7_3_profile_results.csv";
    writetable(profile_results_table,profile_results_file);
end