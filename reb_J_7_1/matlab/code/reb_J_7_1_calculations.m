function reb_J_7_1_calculations
%REB_J_7_CALCULATIONS solve IVODEs in Reaction Engineering Basics 
%   Example J.7
    
    % Set filepath
    filepath_to_results = '../results/';

    % Given and known constants
    V = 2.0; % m^3
    dH_1 = -3470; % cal/mol
    Cp_A = 7.5; % cal/mol/K
    Cp_B = 8.5; % cal/mol/K
    Cp_Y = 12.1; % cal/mol/K
    Cp_Z = 5.7; % cal/mol/K
    k_0_1 = 83; % m^3 /mol /h
    E_1 = 10200; % cal/mol
    % Known
    Re = 1.987; % cal/mol
    Rw = 8.206E-5; % m^3 atm/mol/K

    % Derivatives function
    function derivs = eval_derivs(~, cur_dep)
        % Extract the dependent variables for the current 
        %   integration step
        n_A_cur = cur_dep(1);
        n_B_cur = cur_dep(2);
        n_Y_cur = cur_dep(3);
        n_Z_cur = cur_dep(4);
        T_cur = cur_dep(5);
        
        % Create mass matrix, settinng all elements to zero
        mass_matrix = zeros(6,6);

        % Add 1 on the diagonal for the first 4 rows
        mass_matrix(1,1) = 1;
        mass_matrix(2,2) = 1;
        mass_matrix(3,3) = 1;
        mass_matrix(4,4) = 1;

        % Add the elements for the energy balance
        mass_matrix(5,5) = n_A_cur*Cp_A + n_B_cur*Cp_B + n_Y_cur*Cp_Y ...
            + n_Z_cur*Cp_Z;
        mass_matrix(5,6) = -V*Re/Rw;

        % Add the elements for the ideal gas law equation
        mass_matrix(6,1) = Rw*T_cur;
        mass_matrix(6,2) = Rw*T_cur;
        mass_matrix(6,3) = Rw*T_cur;
        mass_matrix(6,4) = Rw*T_cur;
        mass_matrix(6,5) = Rw*(n_A_cur + n_B_cur + n_Y_cur + n_Z_cur);
        mass_matrix(6,6) = -V;

        % Calculate the rate
        C_A = n_A_cur/V;
        C_B = n_B_cur/V;
        r = k_0_1*exp(-E_1/Re/T_cur)*C_A*C_B;

        % Create right side vector
        rhs = [-V*r; -V*r; V*r; V*r; -V*r*dH_1; 0];

        % Evaluate the derivatives
        derivs = mass_matrix\rhs;
    end

    % Initial values
    ind_0 = 0.0;
    dep_0 = [190.0; 190.0; 0.0; 0.0; 450.0; 7.0];

    % Stopping criterion
    f_var = 0;
    f_val = 2.0;

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
    n_A = dep(:,1);
    n_B = dep(:,2);
    n_Y = dep(:,3);
    n_Z = dep(:,4);
    T = dep(:,5);
    P = dep(:,6);

    % Tabulate the results
    results_table = table(t,n_A,n_B,n_Y,n_Z,T,P);

    % Display the results
    disp(' ')
    disp(results_table)

    % Save the results
    results_file =strcat(filepath_to_results,"reb_J_7_1_results.csv");
    writetable(results_table,results_file);
end