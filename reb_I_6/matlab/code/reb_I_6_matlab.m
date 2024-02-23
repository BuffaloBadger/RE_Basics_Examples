function reb_I_6_matlab
%REB_I_6_MATLAB solve 3 ATEs in Reaction Engineering Basics 
%   Example I.6

    % Given and known constants
    n_A_in = 500.; % mol/h
    V_fluid = 500.; % L
    n_Z_in = 0.; % mol/h
    Cp_vol = 1170.; % cal/L/K
    V_flow = 250.; % L/h
    T_in_K = 423.; % K
    dH_rxn = 18200.; % cal/mol
    k_0 = 1.14E9; % L/mol/h
    E = 16200.; % cal/mol
    Re = 1.987; % cal/mol/K

    % Residuals function
    function residuals = eval_resids(guess)
        % Extract the guess values
        n_A_guess = guess(1);
        n_Z_guess = guess(2);
        T_guess_K = guess(3);

        % Calculate the concentration of A
        C_A = n_A_guess/V_flow;

        % Calculate the rate
        r = k_0*exp(-E/Re/T_guess_K)*C_A^2;
        
        % Evaluate the residuals
        residual_1 = n_A_in - n_A_guess - V_fluid*r;
        residual_2 = n_Z_in - n_Z_guess + V_fluid*r;
        residual_3 = Cp_vol*V_flow*(T_guess_K - T_in_K)...
            + V_fluid*r*dH_rxn;

        % Return the residuals
        residuals = [residual_1; residual_2; residual_3];
    end

    % Initial guess
    init_guess = [n_A_in/2; n_A_in/2; T_in_K - 1.];

    % Solve the ATEs
    [soln, flag, message] = solve_ates(@eval_resids, init_guess);

    % Check that the solution is converged
    if flag <= 0
        disp(' ')
        disp(['The ATE solver did not converge: ',message])
    end

    % Extract the solution
    n_A_out = soln(1);
    n_Z_out = soln(2);
    T_out_K = soln(3);

    % Display the solution
    disp(' ')
    disp(['Flow Rate of A: ',num2str(n_A_out,3),' mol/h'])
    disp(['Flow Rate of Z: ',num2str(n_Z_out,3),' mol/h'])
    disp(['Temperature: ',num2str(T_out_K-273.15,3),' °C'])
end