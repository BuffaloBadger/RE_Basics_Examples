function reb_I_6_calculations
%REB_I_6_CALCULATIONS solve 3 ATEs in Reaction Engineering Basics 
%   Example I.6

    % Set filepath
    filepath_to_results = '../results/';

    % Given and known constants
    n_A_in = 500.; % mol/h
    v_fluid = 500.; % L
    n_Z_in = 0.; % mol/h
    cp_vol = 1170.; % cal/L/K
    v_flow = 250.; % L/h
    temp_in_K = 423.; % K
    dh_rxn = 18200.; % cal/mol
    k0 = 1.14E9; % L/mol/h
    e_act = 16200.; % cal/mol
    gas_const_energy = 1.987; % cal/mol/K

    % Initial guess for the solution
    init_guess = [n_A_in/2; n_A_in/2; temp_in_K - 1.];

    % Solve the ATEs
    [solution, flag, message] = solve_ates(@eval_resids...
        , init_guess);

    % Check that the solution is converged
    if flag <= 0
        disp(' ')
        disp(['The ATE solver did not converge: ',message])
    end

    % calculate the residuals
    residuals = eval_resids(solution);

    % Display and save the results
    n_A_out = solution(1);
    n_Z_out = solution(2);
    temp_out_K = solution(3);
    disp(' ')
    disp(['Flow Rate of A: ',num2str(n_A_out,3),' mol/h'])
    disp(['Flow Rate of Z: ',num2str(n_Z_out,3),' mol/h'])
    disp(['Temperature: ',num2str(temp_out_K-273.15,3),' °C'])

    results_file = strcat(filepath_to_results,"reb_I_6_results.csv");
    item = ["Flow Rate of A";"Flow Rate of Z";"Temperature"...
        ;"Residual 1"; "Residual 2"; "Residual 3"];
    value = [n_A_out;n_Z_out;temp_out_K-273.15; residuals];
    units = ["mol h^-1^";"mol h^-1^";"°C";"";"";""];
    results_table = table(item,value,units);
    writetable(results_table,results_file);

    % Function that evaluates the residuals functions.
    function residuals = eval_resids(guess)
        % Extract the guess values
        n_A_guess = guess(1);
        n_Z_guess = guess(2);
        temp_guess_K = guess(3);

        % Calculate the concentration of A
        conc_A = n_A_guess/v_flow;

        % Calculate the rate
        r = k0*exp(-e_act/gas_const_energy/temp_guess_K)*conc_A^2;
        
        % Evaluate the residuals
        residual_1 = n_A_in - n_A_guess - v_fluid*r;
        residual_2 = n_Z_in - n_Z_guess + v_fluid*r;
        residual_3 = cp_vol*v_flow*(temp_guess_K - temp_in_K)...
            + v_fluid*r*dh_rxn;
        residuals = [residual_1; residual_2; residual_3];
    end
end