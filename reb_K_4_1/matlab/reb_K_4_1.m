function reb_K_4_1
    %REB_K_4_1 Calculations for Example K.4.1 of Reaction Engineering Basics

    % constants available to all functions
	% given
    dH_1 = -14000.; % cal /mol
    k0_1 = 4.2E15; % cm^3 /mol /min
    E_1 = 18000.; % cal /mol
    Cp = 1.3; % cal /cm^3 /K
    CA_0 = 2.0E-3; % mol /cm^3
    CZ_0 = 0; % mol /cm^3
    Vdot_0 = 500.; % cm^3 /min
    T_0 = 300.; % K
    R_R = 1.3; %
    D = 5.; % cm
    L = 50.; % cm
    % known
    R = 1.987; % cal /mol /K

    % derivatives function
    function derivs = derivatives(~, dep)
        % extract ind and dep vars for the current integration step
        nDot_A = dep(1);
        nDot_Z = dep(2);
        T = dep(3);

        % calculate other unknown quantities
        Vdot_3 = Vdot_0;
        Vdot_4 = R_R*Vdot_3;
        Vdot = Vdot_3 + Vdot_4;
        k_1 = k0_1*exp(-E_1/R/T);
        CA = nDot_A/Vdot;
        CZ = nDot_Z/Vdot;
        r_1 = k_1*CA*CZ;

        % evaluate the derivatives
        dnDotAdz = -pi()*D^2/4*r_1;
        dnDotZdz = pi()*D^2/4*r_1;
        dTdz = -pi()*D^2/4*r_1*dH_1/Vdot/Cp;

        % return the derivatives
        derivs = [dnDotAdz; dnDotZdz; dTdz];
    end

    % reactor model
    function [z, nDotA, nDotZ, T] = profiles(dep_0)
        % set the initial values
        ind_0 = 0.0;

        % define the stopping criterion
        f_var = 0;
        f_val = L;
        
        % solve the IVODEs
        odes_are_stiff = false;
        [z, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
            , @derivatives, odes_are_stiff);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: The ODE solution may not be accurate!')
        end

        % extract and return the dependent variable profiles
        nDotA = dep(:,1);
        nDotZ = dep(:,2);
        T = dep(:,3);
    end

    % splitter and mixer residuals function
    function resid = residuals(guess)
        % extract the individual guesses
        nDotA_1 = guess(1);
        nDotZ_1 = guess(2);
        T_1 = guess(3);
        T_4 = guess(4);
        nDotA_3 = guess(5);
        nDotZ_3 = guess(6);
        T_3 = guess(7);

        % solve the reactor design equations
        initial_values = [nDotA_1; nDotZ_1; T_1];
        [~, nDotA, nDotZ, T] = profiles(initial_values);

        % calculate the other unknown quantities
        nDotA_2 = nDotA(end);
        nDotZ_2 = nDotZ(end);
        T_2 = T(end);
        nDotA_0 = Vdot_0*CA_0;
        nDotZ_0 = Vdot_0*CZ_0;
        Vdot_3 = Vdot_0;
        Vdot_4 = R_R*Vdot_3;
        nDotA_4 = R_R*nDotA_3;
        nDotZ_4 = R_R*nDotZ_3;

        % evaluate the residuals
        eps_1 = nDotA_0 + nDotA_4 - nDotA_1;
        eps_2 = nDotZ_0 + nDotZ_4 - nDotZ_1;
        eps_3 = Vdot_0*Cp*(T_1 - T_0) + Vdot_4*Cp*(T_1 - T_4);
        eps_4 = nDotA_2 - nDotA_4 - nDotA_3;
        eps_5 = nDotZ_2 - nDotZ_4 - nDotZ_3;
        eps_6 = T_2 - T_4;
        eps_7 = T_2 - T_3;

        % return the residuals
        resid = [eps_1; eps_2; eps_3; eps_4; eps_5; eps_6; eps_7];
    end

    % splitter and mixer model
    function [nDotA_1, nDotZ_1, T_1, T_4, nDotA_3, nDotZ_3, T_3] ...
            = unknowns(initial_guess)

        % solve the other equipment mole and energy balances
        [soln, flag, message] = solve_ates(@residuals, initial_guess);

        % check that the solution converged
        if flag <= 0
            disp(' ')
            disp(['The ATE solver did not converge: ',message])
        end

        % extract and return the results
        nDotA_1 = soln(1);
        nDotZ_1 = soln(2);
        T_1 = soln(3);
        T_4 = soln(4);
        nDotA_3 = soln(5);
        nDotZ_3 = soln(6);
        T_3 = soln(7);
    end

    % function that performs the analysis
	function perform_the_analysis()

        % initial guess for the unknowns
        nDotA_0 = Vdot_0*CA_0;
        initial_guess = [0.9*nDotA_0; 0.1*nDotA_0; T_0 + 5
            T_0 + 10; 0.1*nDotA_0; 0.1*nDotA_0; T_0 + 10];

        % calculate the unknowns
        [nDotA_1, nDotZ_1, T_1, T_4, nDotA_3, nDotZ_3, T_3] ...
            = unknowns(initial_guess);

        % tabulate the results
        item = ["A 1"; "Z 1"; "T 1"; "T 4"; "A 3"; "Z 3"; "T 3"];
        value = [nDotA_1; nDotZ_1; T_1; T_4; nDotA_3; nDotZ_3; T_3];
        units = ["mol min^-1^"; "mol min^-1^"; "K"; "K"; "mol min^-1^"
            "mol min^-1^"; "K"];
        results_table = table(item,value,units);

        % display the results
        disp(' ')
        disp(results_table)

        % save the results
        writetable(results_table,"results.csv");
    end

    % perform the analysis
    perform_the_analysis()

end