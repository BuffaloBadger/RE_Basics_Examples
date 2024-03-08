function reb_J_7_3
%REB_J_7_3 Calculations for Example J.7.3 of Reaction Engineering Basics
    % given and known constants
    T_f = 325.; % K
    D = 5.; % cm
    dH = -14000.; % cal /mol
    Cp = 1.3; % cal /cm^3 /K
    nDot_A_in = 1.0; % mol /min
    nDot_Z_in = 0.0; % mol /min
    T_in = 300.; % K
    L = 50.0; % cm
    k_0 = 4.2E15; % cm^3 /mol /min
    E = 18000.; % cal /mol
    Re = 1.987; % cal /mol /K
    
    % make IVODE constant available to all functions
    Vdot = nan;

    % reactor design equations as derivative expressions
    function derivs = derivatives(~, dep)
        % extract the dependent variables for this integration step
        nDot_A = dep(1);
        nDot_Z = dep(2);
        T = dep(3);
    
        % calculate rate
        r = other_ivode_variables(T, nDot_A);
    
        % evaluate the derivatives
        dnDotAdz = -pi()*D^2/4*r;
        dnDotZdz = pi()*D^2/4*r;
        dTdz = -pi()*D^2/4*r*dH/Vdot/Cp;
    
        % return the derivatives
        derivs = [dnDotAdz; dnDotZdz; dTdz];
    end

    % calculate other IVODE variables
    function r = other_ivode_variables(T, nDot_A)
        r = k_0*exp(-E/Re/T)*(nDot_A^2)/Vdot^2;
    end

    % calculate unknown IVODE constant
    function Vdot = ivode_constant()
        % initial guess
        initial_guess = 1150.;
    
        % solve the ATE
        [Vdot, flag, message] = solve_ates(@residual, initial_guess);
    
        % check that the solution converged
        if flag <= 0
            disp(' ')
            disp(['The ATE solver did not converge: ',message])
        end
    end

    % implicit equation for the IVODE constant as residual
    function resid = residual(guess)
        % set Vdot
        Vdot = guess;

        % solve the reactor design equations
        [~, ~, ~, T] = profiles();

        % extract the calculated final T
        T_f_calc = T(end);

        % evaluate and return the residual
        resid = T_f - T_f_calc;
    end

    % calculate IVODE initial and final values
    function [ind_0, dep_0, f_var, f_val] = initial_and_final_values()
        % initial values
        ind_0 = 0.0;
        dep_0 = [nDot_A_in; nDot_Z_in; T_in];

        % stopping criterion
        f_var = 0;
        f_val = L;
    end

    % solve the reactor design equations
    function [z, nDot_A, nDot_Z, T] = profiles()
        % get the initial values and stopping criterion
        [ind_0, dep_0, f_var, f_val] = initial_and_final_values();

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

    % calculate other quantities of interest
    function [nDot_A_f, nDot_Z_f] ...
            = other_quantities_of_interest(nDot_A, nDot_Z)
        % calculate the outlet molar flow rates of A and Z
        nDot_A_f = nDot_A(end);
        nDot_Z_f = nDot_Z(end);
    end

    % complete the assignment
    % calculate Vdot
    Vdot = ivode_constant();

    % get the solution of the reactor design equations
    [~, nDot_A, nDot_Z, ~] = profiles();

    % get the quantities of interest
    [nDot_A_f, nDot_Z_f] = other_quantities_of_interest(nDot_A, nDot_Z);

    % Tabulate the results
    item = ["$\dot{V}$";"$\dot{n}_{A,f}$";"$\dot{n}_{Z,f}$"];
    value = [Vdot; nDot_A_f; nDot_Z_f];
    units = ["cm^3^ min^-1^";"mol min^-1^";"mol min^-1^"];
    results_table = table(item,value,units);

    % Display the results
    disp(' ')
    disp(['Volumetric Flow Rate: ', num2str(Vdot,3), ' cm^3^ min^-1^'])
    disp(['Final Molar Flow of A: ', num2str(nDot_A_f,3), ' mol min^-1^'])
    disp(['Final Molar Flow of Z: ', num2str(nDot_Z_f,3), ' mol min^-1^'])

    % Save the results
    results_file = "../results/reb_J_7_3_results.csv";
    writetable(results_table,results_file);
end